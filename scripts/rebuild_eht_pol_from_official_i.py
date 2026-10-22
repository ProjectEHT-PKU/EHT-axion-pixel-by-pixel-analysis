"""Alternate P/D using the actual official-I pipeline checkpoint, modern ehtim1.3.2.

This is a documented implementation of PaperVII parameter choices, not recovered
original PaperVII executable code or proof of identical results.
"""
import argparse
from importlib.metadata import version
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.pipeline.common import (
    stage_directory, verify_stage, require_declared_outputs, manifest_path, verification_session,
)
from rebuild_eht_fiducial import optimizer_records, pol_noise, pol_chi2


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--i-root', type=Path, required=True)
    p.add_argument('--feed-file', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--rounds', type=int, default=50)
    args = p.parse_args(argv)
    if args.output_root.exists() or not 1 <= args.rounds <= 50:
        p.error('New output root / 1..50 rounds required')
    if version('ehtim') != '1.3.2':
        raise ValueError('Pinned ehtim1.3.2 required')
    import ehtim as eh
    from ehtim.calibrating import pol_cal
    from astropy.utils import iers
    iers.conf.auto_download = False
    with verification_session():
        verify_stage(args.i_root, recursive=True, allow_imported=False)
        raw = args.i_root / 'official-I-pre-restoration.fits'
        uv = args.i_root / 'official-I-final-selfcal.uvfits'
        require_declared_outputs(args.i_root, [raw, uv])
        inputs = [Path(__file__).resolve(), Path(__file__).with_name('rebuild_eht_fiducial.py').resolve(),
                  manifest_path(args.i_root), raw, uv, args.feed_file.resolve()]
        inputs += sorted(Path(eh.__file__).resolve().parent.rglob('*.py'))
        config = dict(ehtim='1.3.2', rounds=args.rounds, maxit=100, fov_uas=120, npix=64,
                      pvis=1, m=1, hw=100, ptv=1, fsys=.005, leakage_tol=.1,
                      station_solve=['AZ', 'LM', 'PV'], restart_blur_uas=20,
                      adaptation='official released I checkpoint128uas to120uas; compact flux retained',
                      claim='official-inspired P/D reconstruction, NOT certified official fiducial')
        with stage_directory(args.output_root, 'eht-fiducial-rebuild', config, inputs) as out:
            fixed = eh.image.load_fits(str(raw)).regrid_image(120*eh.RADPERUAS, 64)
            # Regridding conserves surface brightness; preserve declared compact flux0.6 after crop.
            fixed.imvec *= .6/fixed.imvec.sum()
            if np.any(fixed.imvec <= 0) or not np.isfinite(fixed.imvec).all():
                raise ValueError('Positive finite I checkpoint required')
            base = eh.obsdata.load_uvfits(str(uv), polrep='stokes').flag_uvdist(uv_min=1e8)
            table = {r['site']: r for r in json.loads(args.feed_file.read_text())['table']}
            for row in base.tarr:
                for key in ['fr_par', 'fr_elev', 'fr_off']:
                    row[key] = table[str(row['site'])][key]
            base.tarr['dr'] = base.tarr['dl'] = 0j
            base.frcal, base.dcal = True, False
            image = fixed.copy()
            image.add_qu((.05*image.imvec).reshape(64, 64), np.zeros((64, 64)))
            image.add_pol_image(np.zeros((64, 64)), 'V')
            current = pol_noise(base, .005)
            previous_d = {str(s): np.zeros(2, complex) for s in base.tarr['site']}
            steps = []
            for k in range(1, args.rounds+1):
                previous = image
                init = image.copy() if k == 1 else image.blur_circ(20*eh.RADPERUAS, 20*eh.RADPERUAS)
                init.imvec = fixed.imvec.copy()
                power = np.hypot(init.qvec, init.uvec)
                scaling = np.minimum(1, (1-1e-8)*init.imvec/np.maximum(power, np.finfo(float).tiny))
                init.qvec *= scaling
                init.uvec *= scaling
                imgr = eh.imager.Imager(current, init, init, .6, pol='P', transform=['mcv'],
                    data_term={'pvis': 1., 'm': 1.}, reg_term={'hw': 100., 'ptv': 1.},
                    maxit=100, ttype='direct', norm_reg=True, norm_init=False, stop=1e-4)
                with optimizer_records() as actual:
                    image = imgr.make_image_P(show_updates=False)
                    current = pol_cal.leakage_cal(pol_noise(base, .005), image,
                        sites=['AZ', 'LM', 'PV'], leakage_tol=.1, ttype='direct', show_solution=False)
                np.testing.assert_allclose(image.imvec, fixed.imvec, rtol=1e-12, atol=1e-15)
                d = {str(row['site']): np.array([row['dr'], row['dl']]) for row in current.tarr}
                dd = max(float(np.max(abs(v-previous_d[s]))) for s, v in d.items())
                rel = float(np.linalg.norm(image.qvec+1j*image.uvec-previous.qvec-1j*previous.uvec)
                            /max(np.linalg.norm(previous.qvec+1j*previous.uvec), 1e-30))
                row = dict(stage='P+D', round=k, relative_qu_change=rel, max_delta_d=dd,
                           restart_projected_pixels=int(np.count_nonzero(scaling < 1)),
                           chi2_after_D=pol_chi2(current, image, 'direct'), optimizers=actual)
                steps.append(row)
                image.save_fits(str(out/f'round-{k:02d}-raw-IQU.fits'))
                (out/'steps.json').write_text(json.dumps(steps, indent=2)+'\n')
                (out/f'round-{k:02d}-dterms.json').write_text(json.dumps({s:[[float(z.real),float(z.imag)]
                    for z in v] for s,v in d.items()},indent=2)+'\n')
                print('ROUND', k, 'chi2', row['chi2_after_D'], 'deltaD', dd, 'QUchange', rel, flush=True)
                previous_d = d
            image.save_fits(str(out/'final-raw-IQU.fits'))
            image.blur_circ(19*eh.RADPERUAS, 19*eh.RADPERUAS).save_fits(
                str(out/'final-display-IQU-additional19uas.fits'))
            (out/'termination.json').write_text(json.dumps(dict(rounds=args.rounds,
                engineering_stable=rel<.001 and dd<.0001, scientifically_converged=False,
                reason='round_budget_exhausted'),indent=2)+'\n')
            (out/'observation.json').write_text(json.dumps(dict(rows=len(base.data), rf=float(base.rf),
                source=base.source, fixed_I_flux=float(fixed.imvec.sum())),indent=2)+'\n')
            (out/'plan.json').write_text(json.dumps(config,indent=2)+'\n')


if __name__ == '__main__':
    main()
