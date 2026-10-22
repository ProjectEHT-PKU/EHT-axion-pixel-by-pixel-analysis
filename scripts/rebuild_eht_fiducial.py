#!/usr/bin/env python3
"""Isolated official-inspired UVFITS rebuild; NOT bitwise official or a Top Set.

--plan hashes/validates inputs only (no ehtim import or reconstruction). Run in
pixi's optional eht-rebuild environment. --max-rounds 2 --maxit 30 is an
engineering-partial pilot, never a converged scientific fiducial claim.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "b71399714fdfe13579c741075d757b4ce65fd3d1"
NAME = "SR2_M87_2017_101_lo_hops_zbl-dtcal+selfcal.uvfits"
NUMBERS = dict(npix=64, fov_uas=120., flux_jy=.6, prior_fwhm_uas=40., i_maxit=100,
               i_blur_cycles=3, i_systematic_noise=.02, max_rounds=50, maxit=100,
               fsys=.005, leakage_tol=.1, restart_blur_uas=20., display_blur_uas=19.,
               map_rtol=.001, dterm_atol=.0001, seed=0)
I_DATA = {'amp': .2, 'cphase': 1., 'logcamp': 1.}
I_REG = {'simple': 100., 'tv': 1., 'tv2': 1., 'l1': 0., 'flux': 1e4}
P_DATA, P_REG = {'pvis': 1., 'm': 1.}, {'hw': 100., 'ptv': 1.}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=ROOT / 'configs/eht-fiducial.yaml')
    p.add_argument('--output-root', required=True, type=Path, help='new stage directory; never reused')
    p.add_argument('--plan', action='store_true')
    p.add_argument('--max-rounds', type=int)
    p.add_argument('--maxit', type=int, help='override BOTH I and polarization optimizer limits')
    p.add_argument('--ttype', choices=['direct', 'nfft'])
    return p


def build_plan(args):
    import yaml

    c = yaml.safe_load(args.config.read_text())
    extra = {'input_uvfits', 'input_receipt', 'input_sha256', 'commit', 'solve_sites', 'fixed_sites', 'ttype',
             'feed_rotation_file'}
    if not isinstance(c, dict) or set(c) != set(NUMBERS) | extra:
        raise ValueError('config must have exactly the documented fiducial keys')
    for key, default in NUMBERS.items():
        v = c[key]
        if type(v) not in (int, float) or not math.isfinite(v) or v < 0:
            raise ValueError(f'{key}: finite nonnegative numeric value required (no strings/bools)')
        if type(default) is int and (type(v) is not int or (key != 'seed' and v < 1)):
            raise ValueError(f'{key}: positive integer required')
        if v != default:
            raise ValueError(f'{key}: focused fiducial config requires {default}; use explicit CLI limits')
    if (c['commit'] != COMMIT or c['solve_sites'] != ['AZ', 'LM', 'PV']
            or c['fixed_sites'] != ['AA', 'AP', 'JC', 'SM'] or c['ttype'] not in ('direct', 'nfft')):
        raise ValueError('unsupported provenance, sites, or Fourier backend')
    for key in ('max_rounds', 'maxit'):
        v = getattr(args, key)
        if v is not None:
            if not 1 <= v <= NUMBERS[key]:
                raise ValueError(f'{key} override must be in [1,{NUMBERS[key]}]')
            c[key] = v
    if args.maxit is not None:
        c['i_maxit'] = args.maxit
    c['ttype'] = args.ttype or c['ttype']
    inputs = [args.config.resolve(), Path(__file__).resolve()]
    for key in ('feed_rotation_file', 'input_uvfits', 'input_receipt'):
        if not isinstance(c[key], str) or not Path(c[key]).is_absolute():
            raise ValueError(f'{key} must be absolute')
        inputs.append(Path(c[key]).resolve(strict=True))
    receipt = json.loads(inputs[-1].read_text())
    if (inputs[-2].name != NAME or receipt['commit'] != COMMIT
            or receipt['relative_path'] != 'hops_data/April11/' + NAME
            or receipt['source'] != 'https://github.com/eventhorizontelescope/2023-D01-01'
            or receipt['bytes'] != inputs[-2].stat().st_size
            or receipt['sha256'] != c['input_sha256'] or digest(inputs[-2]) != c['input_sha256']):
        raise ValueError('input does not match pinned HOPS April11 receipt (CASA not accepted)')
    out = args.output_root.absolute()
    if any(p.is_symlink() for p in (out, *out.parents)):
        raise ValueError('output symlinks forbidden')
    out = out.resolve()
    if any(out == p or out in p.parents or p in out.parents for p in inputs):
        raise ValueError('output/input overlap')
    if out.exists():
        raise FileExistsError(f'new output root required: {out}')
    partial = c['max_rounds'] < 50 or c['maxit'] < 100
    return dict(label='engineering-partial' if partial else 'official-inspired-full-budget',
                parameters=c, output_root=str(out), inputs=[str(p) for p in inputs],
                hashes={str(p): digest(p) for p in inputs},
                i_data=I_DATA, i_reg=I_REG, p_data=P_DATA, p_reg=P_REG,
                i_initialization='official 1e-3 shifted Gaussian at (40,40) uas; renormalize to 0.6 Jy',
                deviations=['limited I amplitude/closure rebuild; not full Paper IV gain/taper pipeline',
                            'supplied SMILI I selfcal retained; additional phase-only selfcal to rebuilt I',
                            'engineering stability tolerances are not scientific convergence proof'],
                display='additional 19 uas circular kernel on I,Q,U; NOT effective PSF',
                dataset_validation='pinned file hash/receipt; observation metadata checked only on run')


@contextmanager
def optimizer_records():
    """Upstream drops OptimizeResult; capture real nit/status for I, P and leakage."""
    import scipy.optimize as opt

    original, records = opt.minimize, []

    def capture(*a, **kw):
        import numpy as np

        fun, x0 = a[:2]
        initial = float(fun(x0))
        jac = kw.get('jac')
        if not math.isfinite(initial) or (callable(jac) and not np.isfinite(jac(x0)).all()):
            raise RuntimeError('nonfinite initial objective/gradient; optimizer NOT started')
        r = original(*a, **kw)
        if not math.isfinite(float(r.fun)):
            raise RuntimeError('optimizer returned nonfinite objective')
        records.append(dict(initial_objective=initial, nit=int(r.nit), success=bool(r.success),
                            status=int(r.status), message=str(r.message), objective=float(r.fun)))
        return r

    with patch.object(opt, 'minimize', capture):
        yield records


def pol_noise(obs, fsys):
    """Upstream ignores systematic_noise for P: add fsys*|VI| to psigma.

    In Stokes representation psigma²=qsigma²+usigma², so split variance
    equally. I errors stay thermal; upstream propagates these into msigma.
    Apply once to a fresh calibrated copy, never cumulatively.
    """
    import numpy as np

    noisy = obs.switch_polrep('stokes').copy()
    variance = (fsys * np.abs(noisy.data['vis']))**2 / 2
    for field in ('qsigma', 'usigma'):
        noisy.data[field] = np.sqrt(noisy.data[field]**2 + variance)
    return noisy


def pol_chi2(obs, image, ttype):
    # Obsdata.polchisq in 1.3.2 has an undefined vvec; use its lower-level APIs.
    import numpy as np
    from ehtim.imaging import pol_imager_utils as pu

    p = image.qvec + 1j * image.uvec
    arr = np.array([image.imvec, np.abs(p) / image.imvec, np.angle(p), np.zeros_like(image.imvec)])
    return {k: float(pu.polchisq(arr, a, data, sigma, k, ttype=ttype))
            for k in P_DATA
            for data, sigma, a in [pu.polchisqdata(obs, image, [], k, ttype=ttype)]}


def positive_i(image):
    import numpy as np

    if not np.isfinite(image.imvec).all() or not np.all(image.imvec > 0):
        raise ValueError('initial I must be strictly positive and finite')
    return image


def i_prior(obs, c, eh):
    prior = eh.image.make_square(obs, c['npix'], c['fov_uas'] * eh.RADPERUAS)
    w, flux = c['prior_fwhm_uas'] * eh.RADPERUAS, c['flux_jy']
    prior = prior.add_gauss(flux, (w, w, 0, 0, 0))
    # Exact official-I.py symmetry breaker; renormalization preserves requested flux.
    prior = prior.add_gauss(flux * 1e-3, (w, w, 0, w, w))
    prior.imvec *= flux / prior.imvec.sum()
    return positive_i(prior)


def run_rebuild(plan, eh, pol_cal, out):
    import numpy as np
    from astropy.utils import iers

    iers.conf.auto_download = False
    c, steps = plan['parameters'], []
    np.random.seed(c['seed'])
    obs = eh.obsdata.load_uvfits(c['input_uvfits'], polrep='stokes')
    if obs.mjd != 57854 or 'M87' not in obs.source.upper() or not 220e9 < obs.rf < 240e9:
        raise ValueError('expected M87 April11 2017 230GHz observation')
    obs.add_scans()
    obs = obs.avg_coherent(0., scan_avg=True)
    sites = set(obs.data['t1']) | set(obs.data['t2'])
    if not set(c['solve_sites']) <= sites or not sites <= set(c['solve_sites'] + c['fixed_sites']):
        raise ValueError(f'unexpected stations or missing solve sites: {sites}')
    feed = json.loads(Path(c['feed_rotation_file']).read_text())
    coefficients = {row['site']: row for row in feed['table']}
    for row in obs.tarr:
        for key in ('fr_par', 'fr_elev', 'fr_off'):
            row[key] = coefficients[str(row['site'])][key]
    dump(out / 'feed-rotation.json', feed)
    obs.frcal = True  # release explicitly already field-rotation and absolute-EVPA corrected
    obs.tarr['dr'] = obs.tarr['dl'] = 0j  # pre-corrected sites must remain fixed zero
    obs.reorder_tarr_snr()
    # Exclude extended flux / zero baselines in this limited I rebuild and in P.
    obs = obs.flag_uvdist(uv_min=1e8)
    dump(out / 'observation.json', dict(mjd=int(obs.mjd), source=obs.source, rf=float(obs.rf),
                                      rows=len(obs.data), sites=sorted(sites), zero_baselines_removed=True,
                                      fixed_sites_present=sorted(sites & set(c['fixed_sites'])),
                                      dcal_flag_not_proof_of_full_correction=True))
    prior = i_prior(obs, c, eh)
    init = prior
    for cycle in range(c['i_blur_cycles'] + 1):
        positive_i(init)
        imgr = eh.imager.Imager(obs, init, prior, c['flux_jy'], data_term=I_DATA, reg_term=I_REG,
                               maxit=c['i_maxit'], ttype=c['ttype'], norm_reg=True,
                               systematic_noise=c['i_systematic_noise'], transform=['log'], stop=1e-4)
        with optimizer_records() as actual:
            fixed = imgr.make_image_I(show_updates=False)
        steps.append(dict(stage='I', cycle=cycle, optimizers=actual))
        fixed.save_fits(str(out / f'I-cycle-{cycle:02d}.fits'))
        dump(out / 'steps.json', steps)
        init = fixed.blur_circ(obs.res())
    # Always solve total leakage against this SAME selfcal original, never corrected previous-round data.
    with optimizer_records() as actual:
        base = eh.selfcal(obs, fixed, method='phase', ttype=c['ttype'], solution_interval=0., processes=1)
    steps.append(dict(stage='phase-selfcal', optimizers=actual,
                      upstream_budget=10000, iers_auto_download=False))
    base.tarr['dr'] = base.tarr['dl'] = 0j
    base.dcal = False
    current = pol_noise(base, c['fsys'])
    image = fixed.copy()
    image.add_qu(.05 * fixed.imvec.reshape(c['npix'], c['npix']),
                 np.zeros((c['npix'], c['npix'])))
    image.add_pol_image(np.zeros((c['npix'], c['npix'])), 'V')
    previous_d = np.zeros((len(base.tarr), 2), complex)
    stable = False
    for round_id in range(1, c['max_rounds'] + 1):
        previous = image
        init = image.blur_circ(c['restart_blur_uas'] * eh.RADPERUAS,
                               c['restart_blur_uas'] * eh.RADPERUAS) if round_id > 1 else image.copy()
        init.imvec = fixed.imvec.copy()  # blur Q/U only effectively; restore SAME fixed I
        # Preserve physical mcv domain: smoothing P then replacing I can violate P<I.
        polarization = np.hypot(init.qvec, init.uvec)
        ceiling = (1 - 1e-8) * init.imvec
        projected = polarization >= ceiling
        scale = np.minimum(1, ceiling / np.maximum(polarization, np.finfo(float).tiny))
        init.qvec *= scale
        init.uvec *= scale
        fitobs = current.copy()
        imgr = eh.imager.Imager(fitobs, init, init, c['flux_jy'], pol='P', transform=['mcv'],
                               data_term=P_DATA, reg_term=P_REG, maxit=c['maxit'],
                               ttype=c['ttype'], norm_reg=True, norm_init=False, stop=1e-4)
        with optimizer_records() as actual:
            image = imgr.make_image_P(show_updates=False)
            chi_before = pol_chi2(fitobs, image, c['ttype'])
            current = pol_cal.leakage_cal(pol_noise(base, c['fsys']), image, sites=c['solve_sites'],
                                         leakage_tol=c['leakage_tol'], ttype=c['ttype'], show_solution=False)
        np.testing.assert_allclose(image.imvec, fixed.imvec, rtol=1e-12, atol=1e-15)
        d = np.column_stack((current.tarr['dr'], current.tarr['dl']))
        if np.any(d[[s in c['fixed_sites'] for s in current.tarr['site']]] != 0):
            raise ValueError('pre-corrected sites acquired nonzero leakage')
        relmap = float(np.linalg.norm((image.qvec + 1j*image.uvec) - (previous.qvec + 1j*previous.uvec))
                       / max(np.linalg.norm(previous.qvec + 1j*previous.uvec), 1e-30))
        delta = float(np.max(np.abs(d - previous_d)))
        # current already carries the once-inflated errors from the fresh base solve.
        diagnostic_obs = current
        chi = pol_chi2(diagnostic_obs, image, c['ttype'])
        stable = round_id >= 3 and relmap < c['map_rtol'] and delta < c['dterm_atol']
        row = dict(stage='P+D', round=round_id, relative_qu_change=relmap, max_delta_d=delta,
                   restart_projected_pixels=int(projected.sum()), leakage_upstream_maxit=1000,
                   chi2_before_D=chi_before, chi2_after_D=chi, optimizers=actual, engineering_stable=stable)
        steps.append(row)
        image.save_fits(str(out / f'round-{round_id:02d}-raw-IQU.fits'))
        dump(out / f'round-{round_id:02d}.json', row)
        dump(out / f'round-{round_id:02d}-dterms.json',
             {str(s): [[float(z.real), float(z.imag)] for z in pair]
              for s, pair in zip(current.tarr['site'], d)})
        dump(out / 'steps.json', steps)
        previous_d = d.copy()
        if stable:
            break
    image.save_fits(str(out / 'final-raw-IQU.fits'))
    beam = c['display_blur_uas'] * eh.RADPERUAS
    image.blur_circ(beam, beam).save_fits(str(out / 'final-display-IQU-additional19uas.fits'))
    dump(out / 'termination.json', dict(rounds=round_id, engineering_stable=stable,
         reason='engineering_stability' if stable else 'round_budget_exhausted',
         scientifically_converged=False))


def main(argv=None):
    args = parser().parse_args(argv)
    plan = build_plan(args)
    if args.plan:
        print(json.dumps(plan, indent=2))
        return plan
    import _bootstrap  # noqa: F401
    from eht_axion.pipeline.common import reject_directory_symlinks, stage_directory

    if metadata.version('ehtim') != '1.3.2':
        raise RuntimeError('requires official PyPI ehtim==1.3.2 in optional pixi eht-rebuild')
    import ehtim as eh
    from ehtim.calibrating import pol_cal

    package = Path(eh.__file__).resolve().parent
    dist = metadata.distribution('ehtim')
    code = sorted(package.rglob('*.py'))
    receipts = [Path(dist.locate_file(f)).resolve() for f in dist.files or []
                if Path(str(f)).name in ('METADATA', 'RECORD', 'direct_url.json')]
    plan['runtime'] = dict(version=dist.version, module_version=str(getattr(eh, '__version__', 'undeclared')),
                           loaded_path=str(package),
                           code_hashes={str(p): digest(p) for p in code})
    inputs = [Path(p) for p in plan['inputs']] + code + receipts
    destination = Path(plan['output_root'])
    reject_directory_symlinks(destination)
    if destination.exists():
        raise FileExistsError(destination)
    with stage_directory(destination, 'eht-fiducial-rebuild', plan, inputs) as out:
        if out is None:
            raise FileExistsError('rebuild stages may not reuse existing outputs')
        dump(out / 'plan.json', plan)
        run_rebuild(plan, eh, pol_cal, out)
    return plan


if __name__ == '__main__':
    main()
