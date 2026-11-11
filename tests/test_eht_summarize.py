"""scripts/eht/summarize.py: experiment-derived settings, checks of a Monte Carlo root (synthetic)."""
import json

import numpy as np
from pathlib import Path

import pytest

from scripts.eht import summarize as s

ROOT = Path(__file__).resolve().parents[1]


def test_settings_come_from_experiment():
    from eht_axion.experiment import read_experiment

    settings = s.summary_settings(read_experiment())
    assert settings == dict(draws=1000, shape=64, fov_uas=120.0)


SCIENCE = {'draws_per_day': 2, 'seed': 1}


def sample_root(tmp_path, count=2):
    results = tmp_path / 'mc'
    results.mkdir(parents=True)
    draws = {d: {'draws': [{'AA': {'R': [0., 0.], 'L': [0., 0.]}} for _ in range(count)]} for d in s.en.DAYS}
    (results / 'all-draws.json').write_text(json.dumps(draws))
    tasks, receipts = [], {}
    for d in s.en.DAYS:
        for k in range(count):
            task = dict(id=f'{d}-{k:04d}', day=d, draw=k, dterms=draws[d]['draws'][k])
            tasks.append(task)
            directory = results / 'attempts' / task['id'] / 'attempt-0001'
            directory.mkdir(parents=True)
            (directory / 'best-blurred-IQU.fits').write_bytes(b'blurred')
            (directory / 'certificate.json').write_text(json.dumps(
                dict(solver=dict(numerically_certified=True, best_gap=0.), gap_atol=1e-6)))
            (directory / 'taskreceipt.json').write_text(json.dumps(dict(
                task=task, status='certified', schema='eht-sample-task',
                outputs=['best-blurred-IQU.fits', 'certificate.json'])))
            receipts[task['id']] = str(directory.relative_to(results))
    (results / 'batch.json').write_text(json.dumps({'schema': 'eht-sample-batch', 'science': SCIENCE}))
    finish = dict(schema='eht-sample-finish', numerically_certified=4 * count, receipts=receipts)
    (results / 'finish.json').write_text(json.dumps(finish))
    return results


def test_verify_batch_reads_the_sample_root_directly(tmp_path):
    root = sample_root(tmp_path, count=2)
    selected, paths = s.verify_batch(root, dict(draws=2), SCIENCE)
    assert {d: sorted(v) for d, v in selected.items()} == {d: [0, 1] for d in s.en.DAYS}
    assert root / 'batch.json' in paths and root / 'finish.json' in paths
    assert not any(p.name == 'tasks.jsonl' for p in paths)
    with pytest.raises(ValueError, match='inconsistent Monte Carlo batch'):
        s.verify_batch(root, dict(draws=3), SCIENCE)
    with pytest.raises(ValueError, match='different experiment values'):
        s.verify_batch(root, dict(draws=2), dict(SCIENCE, seed=2))


def test_run_writes_one_ensemble_of_all_days(tmp_path, monkeypatch):
    n, draws = 4, 2
    east, north = -(np.arange(n) - 1.5) * 30., (np.arange(n) - 1.5) * 30.

    def load(path):
        stokes = np.zeros((n, n, 4))
        stokes[..., 0], stokes[..., 1] = 1., float(path)
        return dict(stokes=stokes, east_uas=east, north_uas=north, frequency_hz=227e9, bunit='JY/PIXEL')

    monkeypatch.setattr(s, 'load_stokes_fits', load)
    selected = {d: {k: dict(blurred=f'0.{i}{k}') for k in range(draws)} for i, d in enumerate(s.en.DAYS)}
    s.run(selected, tmp_path, dict(draws=draws, shape=n, fov_uas=120.))
    assert {p.name for p in tmp_path.iterdir()} == {'stokes.npy', 'global.npz', 'axes.npz', 'report.json'}
    stokes = np.load(tmp_path / 'stokes.npy')
    assert stokes.shape == (4, draws, n, n, 4) and stokes[2, 1, 0, 0, 1] == .21
    with np.load(tmp_path / 'global.npz') as z:
        assert z['values'].shape == (4, draws, len(s.en.GLOBAL_COLUMNS))
        assert list(z['days']) == list(s.en.DAYS)
    assert json.loads((tmp_path / 'report.json').read_text())['map_order'].startswith('day,draw,')
