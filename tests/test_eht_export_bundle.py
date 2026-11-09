"""eht_export_bundle: synthetic I and refine stages -> bundle accepted by eht_sample (no imaging)."""
import json
from pathlib import Path

import pytest

import eht_export_bundle as w
import eht_sample
from eht_axion import experiment
from eht_axion import provenance

ROOT = Path(__file__).resolve().parents[1]


def stages(tmp_path, config, *, certified=True, days=w.DAYS, feed=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    if feed is None:
        feed = tmp_path / 'feed-rotation.json'
        feed.write_text('{"table": []}')
    locations = dict(input_root='/data/uvfits', input_receipt='/data/receipt.json',
                     feed_rotation_file=str(feed), official_script='/data/official.py')
    parameters = experiment.eht_imaging_parameters(config, locations)
    i_roots, d_roots = [], []
    for day in days:
        i_root = tmp_path / f'I-{day}'
        with provenance.stage_directory(i_root, 'eht:i',
                                        dict(day=day, parameters=parameters), [feed]) as out:
            (out / 'final-fixedI.fits').write_bytes(f'I {day}'.encode())
            (out / 'clean-gaincal-pol.uvfits').write_bytes(f'obs {day}'.encode())
        d_root = tmp_path / f'D-{day}'
        inputs = [i_root / 'manifest.json', i_root / 'final-fixedI.fits', i_root / 'clean-gaincal-pol.uvfits']
        with provenance.stage_directory(d_root, 'eht:refine',
                                    dict(day=day, parameters=parameters), inputs) as out:
            (out / 'final-dterms.json').write_text(json.dumps({'AZ': [[.01, 0], [0, .01]]}))
            (out / 'termination.json').write_text(json.dumps(dict(all_P_numerically_certified=certified)))
            (out / 'steps.json').write_text('[]')
            (out / 'plan.json').write_text('{}')
        i_roots.append(i_root)
        d_roots.append(d_root)
    return i_roots, d_roots, feed


def test_bundle_contract_files_and_eht_sample_acceptance(tmp_path):
    config = experiment.read_experiment()
    i_roots, d_roots, feed = stages(tmp_path, config)
    before = {p: p.read_bytes() for root in i_roots + d_roots for p in root.iterdir()}
    record = w.export(config, list(reversed(i_roots)), d_roots, tmp_path / 'bundle')
    bundle = tmp_path / 'bundle'
    assert json.loads((bundle / 'bundle.json').read_text()) == record
    c = config['eht_imaging']
    assert (record['fsys'], record['covariance_ddof'], record['residual_component_sigma']) == (
        c['fsys'], c['covariance_ddof'], c['residual_d_sigma'])
    assert record['image'] == dict(npix=64, fov_uas=120.0, blur_fwhm_uas=19.0)
    assert record['objective'] == dict(pvis=1.0, m=1.0, hw=100.0, ptv=1.0, norm_reg=True)
    assert list(record['days']) == ['95', '96', '100', '101']
    assert record['days']['96']['D_stage'] == str(d_roots[1].resolve())
    paths = [f['path'] for f in record['files']]
    assert paths[:8] == ['95/fixedI.fits', '95/gaincal-pol.uvfits', '95/fiducial-dterms.json',
                         'lineage/95/D-manifest.json', 'lineage/95/D-plan.json',
                         'lineage/95/D-termination.json', 'lineage/95/D-steps.json',
                         'lineage/95/I-manifest.json']
    assert paths[-2:] == ['feed-rotation.json', 'eht_export_bundle.py']
    assert (bundle / '100/fixedI.fits').read_bytes() == b'I 100'
    assert not (bundle / '100/fixedI.fits').stat().st_mode & 0o222
    assert {p: p.read_bytes() for p in before} == before  # sources untouched
    science = eht_sample.science_settings(ROOT / 'experiment.yaml')
    verified, lineage = eht_sample.verify_bundle(bundle, science)
    assert verified == record and all(f['readonly_mode'] for f in lineage['files'].values())
    with pytest.raises(FileExistsError):
        w.export(config, i_roots, d_roots, bundle)


def test_rejects_uncertified_missing_day_and_foreign_parameters(tmp_path):
    config = experiment.read_experiment()
    i_roots, d_roots, _ = stages(tmp_path / 'a', config, certified=False)
    with pytest.raises(ValueError, match='certified'):
        w.export(config, i_roots, d_roots, tmp_path / 'out1')
    with pytest.raises(ValueError, match='exactly the days'):
        w.export(config, i_roots[:3], d_roots, tmp_path / 'out2')
    changed = json.loads(json.dumps(config))
    changed['eht_imaging']['fsys'] = 0.01
    i_roots, d_roots, feed = stages(tmp_path / 'b', config)
    with pytest.raises(ValueError, match='differ from experiment.yaml'):
        w.export(changed, i_roots, d_roots, tmp_path / 'out3')
    _, other_d_roots, _ = stages(tmp_path / 'c', config, feed=feed)
    with pytest.raises(ValueError, match='does not descend'):
        w.export(config, i_roots, other_d_roots, tmp_path / 'out4')
    assert not any((tmp_path / f'out{k}').exists() for k in range(1, 5))
