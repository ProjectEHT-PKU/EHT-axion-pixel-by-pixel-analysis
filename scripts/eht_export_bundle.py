#!/usr/bin/env python3
"""EHT imaging: export the verified fiducial inputs of the D-term Monte Carlo into one bundle.

For each of the four days, copy the fixed I image and gain-calibrated data of the
I stage (eht_reconstruct.py --stage I) and the refined fiducial D terms of the
tightD stage (eht_refine.py), plus their manifests/plans as lineage, and the
feed-rotation table recorded by the I stages. bundle.json lists every file with
its SHA-256 and the science contract checked by eht_sample.py (fsys, D-term
statistics convention, sky grid, imaging weights) taken from experiment.yaml.

    pixi run -e eht-rebuild eht-export-bundle --i-root I95 ... --d-root D95 ... --output-root BUNDLE

Both stage kinds must be complete, lineage-verified, of the same day set, made
with the eht_imaging parameters of experiment.yaml, and every P solve of the
tightD stage must be numerically certified. Source stages are only read.
"""
from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path
import shutil

from _bootstrap import ROOT  # noqa: F401
from eht_axion import experiment
from eht_axion import provenance

DAYS = experiment.EHT_DAYS
I_FILES = {'final-fixedI.fits': 'fixedI.fits', 'clean-gaincal-pol.uvfits': 'gaincal-pol.uvfits'}
D_FILES = ('final-dterms.json', 'termination.json', 'steps.json')
D_LINEAGE = ('manifest.json', 'plan.json', 'termination.json', 'steps.json')
SCHEMA = 'eht-bundle-v2'


def _by_day(roots, stage, config):
    """Verify each root (full lineage) and index it by its recorded day."""
    result = {}
    for root in roots:
        root = Path(root).resolve(strict=True)
        manifest = provenance.verify_stage(root, recursive=True)
        if manifest['stage'] != stage:
            raise ValueError(f'{root}: expected stage {stage}, found {manifest["stage"]}')
        parameters = manifest['config']['parameters']
        locations = {key: parameters[key] for key in experiment.EHT_LOCATIONS}
        expected = json.loads(json.dumps(experiment.eht_imaging_parameters(config, locations)))
        if parameters != expected:
            raise ValueError(f'{root}: imaging parameters differ from experiment.yaml eht_imaging')
        day = manifest['config']['day']
        if day in result:
            raise ValueError(f'two {stage} roots for day {day}')
        result[day] = (root, manifest)
    if set(result) != set(DAYS):
        raise ValueError(f'{stage}: need exactly the days {DAYS}, got {sorted(result)}')
    return result


def verify_sources(i_roots, d_roots, config):
    """Return {day: (I root, D root)} and the feed-rotation table recorded by the I stages."""
    images = _by_day(i_roots, 'eht:i', config)
    dterms = _by_day(d_roots, 'eht:refine', config)
    feeds = {manifest['config']['parameters']['feed_rotation_file'] for _, manifest in images.values()}
    if len(feeds) != 1:
        raise ValueError('the four I stages used different feed-rotation tables')
    sources = {}
    for day in DAYS:
        i_root, i_manifest = images[day]
        d_root, d_manifest = dterms[day]
        if d_manifest['config']['parameters'] != i_manifest['config']['parameters']:
            raise ValueError(f'day {day}: I and D stage parameters differ')
        provenance.require_declared_outputs(i_root, [i_root / name for name in I_FILES])
        provenance.require_declared_outputs(d_root, [d_root / name for name in D_FILES])
        closure = set(provenance.source_closure([provenance.manifest_path(d_root)]))
        if any((i_root / name).resolve() not in closure for name in I_FILES):
            raise ValueError(f'day {day}: the D stage does not descend from this I stage')
        termination = json.loads((d_root / 'termination.json').read_text())
        if termination.get('all_P_numerically_certified') is not True:
            raise ValueError(f'day {day}: not every P solve of the D stage is certified')
        sources[day] = (i_root, d_root)
    feed = Path(feeds.pop()).resolve(strict=True)
    recorded = {Path(e['path']).resolve() for _, m in images.values() for e in m['inputs']}
    if feed not in recorded:
        raise ValueError('feed-rotation table is not a declared input of the I stages')
    return sources, feed


def bundle_record(config):
    """The science contract of the bundle (checked by eht_sample.verify_bundle)."""
    c, image = config['eht_imaging'], config['shared']['image']
    return dict(schema=SCHEMA, created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                days={}, files=[], feed_file='feed-rotation.json',
                covariance_ddof=c['covariance_ddof'], residual_component_sigma=c['residual_d_sigma'],
                fsys=c['fsys'],
                image=dict(npix=image['shape'], fov_uas=image['fov_uas'],
                           blur_fwhm_uas=image['blur_fwhm_uas']),
                objective=dict(pvis=c['p_data']['pvis'], m=c['p_data']['m'], hw=c['p_reg']['hw'],
                               ptv=c['p_reg']['ptv'], norm_reg=c['norm_reg']),
                disclaimer='D calibration uncertainty under declared four-day scatter model, '
                           'not all systematics')


def export(config, i_roots, d_roots, output_root):
    output_root = Path(output_root).absolute()
    provenance.reject_directory_symlinks(output_root)
    if output_root.exists():
        raise FileExistsError(f'new output root required: {output_root}')
    with provenance.verification_session():
        sources, feed = verify_sources(i_roots, d_roots, config)
        record = bundle_record(config)
        output_root.mkdir(parents=True)

        def copy(source, relative):
            destination = output_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            destination.chmod(0o444)
            digest = provenance.sha256(destination)
            if digest != provenance.sha256(source):
                raise ValueError(f'copy differs from its source: {source}')
            record['files'].append(dict(path=relative, sha256=digest, size=destination.stat().st_size,
                                        original_path=str(source)))

        for day in DAYS:
            i_root, d_root = sources[day]
            for name, target in I_FILES.items():
                copy(i_root / name, f'{day}/{target}')
            copy(d_root / 'final-dterms.json', f'{day}/fiducial-dterms.json')
            for name in D_LINEAGE:
                copy(d_root / name, f'lineage/{day}/D-{name}')
            copy(i_root / 'manifest.json', f'lineage/{day}/I-manifest.json')
            record['days'][str(day)] = {'I': f'{day}/fixedI.fits', 'obs': f'{day}/gaincal-pol.uvfits',
                                        'dterms': f'{day}/fiducial-dterms.json', 'I_stage': str(i_root),
                                        'D_stage': str(d_root), 'all_P_numerically_certified': True}
        copy(feed, 'feed-rotation.json')
        copy(Path(__file__).resolve(), Path(__file__).name)
    (output_root / 'bundle.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=experiment.DEFAULT_PATH, help='default: experiment.yaml')
    p.add_argument('--i-root', type=Path, action='append', required=True,
                   help='completed I stage; give all four days')
    p.add_argument('--d-root', type=Path, action='append', required=True,
                   help='completed eht_refine (tightD) stage; give all four days')
    p.add_argument('--output-root', type=Path, required=True, help='new bundle directory')
    args = p.parse_args(argv)
    record = export(experiment.read_experiment(args.config), args.i_root, args.d_root, args.output_root)
    print(f"EXPORTED {args.output_root}: {len(record['files'])} files, "
          f"{sum(f['size'] for f in record['files'])} bytes")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
