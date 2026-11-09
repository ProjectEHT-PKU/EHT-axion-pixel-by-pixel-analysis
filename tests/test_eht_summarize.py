"""eht_summarize: experiment-derived settings and verification of an eht_sample output root (synthetic)."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('eht_summarize_test', ROOT / 'scripts/eht_summarize.py')
s = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_settings_come_from_experiment():
    from eht_axion.experiment import read_experiment

    settings = s.summary_settings(read_experiment())
    assert settings == dict(draws=1000, shape=64, fov_uas=120.0, rmax_uas=60.0, window_degrees=10)


SCIENCE = {'draws_per_day': 2, 'seed': 1}


def sample_root(tmp_path, count=2):
    results = tmp_path / 'mc'
    results.mkdir(parents=True)
    draws = {d: {'draws': [{'AA': {'R': [0., 0.], 'L': [0., 0.]}} for _ in range(count)]} for d in s.en.DAYS}
    (results / 'all-draws.json').write_text(json.dumps(draws))
    tasks, receipts = [], {}
    identity_plan = {'schema': 'eht-sample-batch-v2', 'science': SCIENCE}
    encoded = json.dumps(identity_plan, sort_keys=True, separators=(',', ':')).encode()
    identity = hashlib.sha256(encoded).hexdigest()
    for d in s.en.DAYS:
        for k in range(count):
            task = dict(id=f'{d}-{k:04d}', day=d, draw=k, dterms=draws[d]['draws'][k])
            tasks.append(task)
            directory = results / 'attempts' / task['id'] / 'attempt-0001'
            directory.mkdir(parents=True)
            (directory / 'best-raw-IQU.fits').write_bytes(b'raw')
            (directory / 'best-display-IQU-additional19uas.fits').write_bytes(b'display')
            (directory / 'certificate.json').write_text(json.dumps(
                dict(solver=dict(numerically_certified=True, best_gap=0.), gap_atol=1e-6)))
            names = ('best-raw-IQU.fits', 'best-display-IQU-additional19uas.fits', 'certificate.json')
            outputs = {n: sha(directory / n) for n in names}
            (directory / 'taskreceipt.json').write_text(json.dumps(dict(
                task=task, identity=identity, status='certified', schema='eht-sample-task-v2',
                outputs=outputs)))
            receipts[task['id']] = dict(directory=str(directory.relative_to(results)),
                                        sha256=sha(directory / 'taskreceipt.json'))
    (results / 'tasks.jsonl').write_text(''.join(json.dumps(t) + '\n' for t in tasks))
    (results / 'batch.json').write_text(json.dumps(dict(identity_plan, identity=identity)))
    finish = dict(identity=identity, schema='eht-sample-finish-v2', numerically_certified=4 * count,
                  receipts=receipts)
    (results / 'finish.json').write_text(json.dumps(finish))
    return results


def test_verify_batch_reads_the_sample_root_directly(tmp_path):
    root = sample_root(tmp_path, count=2)
    with s.provenance.verification_session():
        selected, paths, identity = s.verify_batch(root, dict(draws=2), SCIENCE)
    assert {d: sorted(v) for d, v in selected.items()} == {d: [0, 1] for d in s.en.DAYS}
    assert identity == json.loads((root / 'batch.json').read_text())['identity']
    assert root / 'batch.json' in paths and root / 'finish.json' in paths
    with pytest.raises(ValueError, match='inconsistent eht_sample batch'):
        with s.provenance.verification_session():
            s.verify_batch(root, dict(draws=3), SCIENCE)
    with pytest.raises(ValueError, match='different experiment values'):
        with s.provenance.verification_session():
            s.verify_batch(root, dict(draws=2), dict(SCIENCE, seed=2))
    batch = json.loads((root / 'batch.json').read_text())
    (root / 'batch.json').write_text(json.dumps(dict(batch, science=dict(SCIENCE, seed=2))))
    with pytest.raises(ValueError, match='inconsistent eht_sample batch'):  # identity no longer matches
        with s.provenance.verification_session():
            s.verify_batch(root, dict(draws=2), dict(SCIENCE, seed=2))
