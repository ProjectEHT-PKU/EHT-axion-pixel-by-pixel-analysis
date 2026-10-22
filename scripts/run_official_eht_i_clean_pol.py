"""Execute a frozen official PaperIV I pipeline, with explicit backend-only adaptation.

No rewritten I optimizer. Preserves published stages/self-calibration/parameters.
Uses newer ehtim1.3.2, so this is not a bitwise historical reproduction.
"""
from pathlib import Path
import argparse
import hashlib
import json
from importlib.metadata import version
import runpy
import sys

import _bootstrap  # noqa: F401
from eht_axion.pipeline.common import stage_directory


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--official-script', type=Path, required=True)
    p.add_argument('--input-uvfits', type=Path, required=True)
    p.add_argument('--input-receipt', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    args = p.parse_args(argv)
    if args.output_root.exists():
        raise FileExistsError('New official-I output directory required')
    if version('ehtim') != '1.3.2':
        raise ValueError('Requires pinned ehtim1.3.2')
    text = args.official_script.read_text()
    anchor = "ttype     = 'nfft'"
    if text.count(anchor) != 1:
        raise ValueError('Official backend adaptation anchor differs')
    adapted = text.replace(anchor, "ttype     = 'direct'", 1)
    import ehtim as eh
    from astropy.utils import iers
    iers.conf.auto_download = False
    code = sorted(Path(eh.__file__).resolve().parent.rglob('*.py'))
    identity = dict(ehtim_version=version('ehtim'),
                    official_script_sha256=hashlib.sha256(text.encode()).hexdigest(),
                    adaptations=['nfft to direct DTFT; IERS offline'],
                    scope='official released I stages with modern library; not bitwise historical')
    with stage_directory(args.output_root, 'research:official-eht-I', identity,
                         [Path(__file__).resolve(), args.official_script.resolve(),
                          args.input_uvfits.resolve(), args.input_receipt.resolve()] + code) as out:
        snapshot = out / 'official-I-original.py'
        snapshot.write_text(text)
        executable = out / 'official-I-direct.py'
        executable.write_text(adapted)
        sys.argv = [str(executable), '-i', str(args.input_uvfits.resolve()),
                    '-o', str(out / 'official-I-restored.fits')]
        ns = runpy.run_path(str(executable), run_name='__main__')
        # Preserve raw I and the pipeline's final self-calibrated data for P stage.
        ns['imgr'].out_last().save_fits(str(out / 'official-I-pre-restoration.fits'))
        ns['obs_sc'].save_uvfits(str(out / 'official-I-final-selfcal.uvfits'))
        # Apply final gain table to original scan-averaged observations, which
        # have neither reverse taper nor the I-only 2% noise inflation.
        clean_source = ns['obs_orig'].copy()
        # obs_orig predates SNR-based station sorting; align by station identity,
        # never by array position, before applying the same final gains.
        order = [list(clean_source.tarr['site']).index(site) for site in ns['caltab'].tarr['site']]
        clean_source.tarr = clean_source.tarr[order].copy()
        clean_pol = ns['caltab'].applycal(clean_source, interp='nearest', extrapolate=True)
        clean_pol.save_uvfits(str(out / 'official-I-gaincal-clean-pol.uvfits'))
        (out / 'polarization-input.json').write_text(json.dumps({
            'source': 'obs_orig: before reverse taper / I systematic noise',
            'gains': 'final caltab solved by released I pipeline',
            'noise': 'original errors propagated by gain application; no I2% additive inflation',
            'remaining_D_sites': ['AZ', 'LM', 'PV'], 'frcal_from_release': True}, indent=2)+'\n')
        metric = {key: float(ns['obs_sc'].chisq(ns['imgr'].out_last(), dtype=key, ttype='direct'))
                  for key in ['vis', 'amp', 'cphase', 'logcamp']}
        (out / 'report.json').write_text(json.dumps(dict(identity, chi2_final_selfcal=metric,
            original_fov_uas=float(ns['fov']/eh.RADPERUAS), pixels=int(ns['npix']),
            note='No polarimetric D-term solve in this I-only official script.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
