"""Synthetic snapshots and tiny exact PDHG only; no real observations/MC launch."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'eht_sample.py'
spec = importlib.util.spec_from_file_location('uncertainty_batch', SCRIPT)
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)


def snapshot(tmp_path):
    root = tmp_path / 'inputs'
    root.mkdir()
    days, files = {}, []
    for day in batch.DAYS:
        days[day] = dict(I=f'{day}/fixedI.fits', obs=f'{day}/gaincal-pol.uvfits',
                         dterms=f'{day}/fiducial-dterms.json')
        for kind, rel in days[day].items():
            p = root / rel
            p.parent.mkdir(exist_ok=True)
            if kind == 'dterms':
                value = {s: [[.01 * int(day), .002], [.003, -.004 * int(day)]]
                         for s in ('AZ', 'LM', 'PV')}
                p.write_text(json.dumps(value))
            else:
                p.write_bytes((day + kind).encode())
            p.chmod(0o444)
            files.append(dict(path=rel, sha256=batch.sha(p)))
    feed = root / 'feed-rotation.json'
    feed.write_text('{}')
    feed.chmod(0o444)
    files.append(dict(path=feed.name, sha256=batch.sha(feed)))
    b = dict(schema='eht-bundle-v2', days=days, files=files, feed_file=feed.name,
             covariance_ddof=1, residual_component_sigma=.01, fsys=.005,
             objective=dict(pvis=1.0, m=1.0, hw=100.0, ptv=1.0, norm_reg=True))
    (root / 'bundle.json').write_text(json.dumps(b))
    return root, b


def test_bundle_bytes_seeds_and_4000_unique_tasks(tmp_path):
    root, b = snapshot(tmp_path)
    science = batch.science_settings(SCRIPT.parents[1] / 'experiment.yaml')
    assert science['sites'] == ['AA', 'AP', 'AZ', 'JC', 'LM', 'PV', 'SM']
    assert (science['draws_per_day'], science['seed'], science['fsys']) == (1000, 20210405, .005)
    # The settings recorded in the batch config of eht-uncertainty-remote-v1.
    assert science['monte_carlo'] == dict(initial_p_fraction=.05, gap_atol=1e-6, maxiter=50000,
                                          check_every=500, time_limit=600., retry_maxiter=200000,
                                          retry_time_limit=1800.)
    verified, lineage = batch.verify_bundle(root, science)
    assert verified == b
    assert all(f['readonly_mode'] for f in lineage['files'].values())
    assert lineage['bundle_sha256'] == batch.sha(root / 'bundle.json')
    draws = batch.make_draws(root, b, 1000, 20210405, science)
    tasks = batch.enumerate_tasks(draws)
    assert len(tasks) == len({t['id'] for t in tasks}) == 4000
    assert all(sum(t['day'] == d for t in tasks[:64]) == 16 for d in batch.DAYS)
    assert {t['day'] for t in tasks} == set(batch.DAYS)
    assert all(len(draws[d]['draws']) == 1000 for d in batch.DAYS)
    assert all(set(t['dterms']) == set(science['sites']) for t in tasks)
    repeated = batch.make_draws(root, dict(b, days=dict(reversed(list(b['days'].items())))), 1000, 20210405,
                                science)
    assert repeated == draws
    assert len({draws[d]['seed'] for d in batch.DAYS}) == 4
    assert draws['95']['draws'][0] != draws['96']['draws'][0]
    p = root / b['days']['95']['I']
    p.chmod(0o644)
    p.write_bytes(b'tampered')
    with pytest.raises(ValueError, match='hash mismatch'):
        batch.verify_bundle(root, science)
    with pytest.raises(ValueError, match='science contract'):
        batch.verify_bundle(root, dict(science, fsys=.01))
    with pytest.raises(ValueError, match='science contract'):
        batch.verify_bundle(root, dict(science, p_reg=dict(hw=50.0, ptv=1.0)))
    with pytest.raises(ValueError, match='Unsafe'):
        batch.contained(root, '../escape')


def test_receipts_retry_and_tiny_exact_solver(tmp_path):
    from eht_axion.eht import convex as opt
    intensity = np.ones(4)
    A = np.eye(4, dtype=complex)
    problem = opt.CartesianProblem(intensity, (2, 2),
                {name: (np.zeros(4, complex), np.ones(4), A) for name in ('pvis', 'm')}, 1., 4.,
                weights=dict(pvis=1., m=1., hw=100., ptv=1.))
    task = dict(id='95-0000', day='95', draw=0, dterms={})
    config = dict(maxiter=1, time_limit=10, gap_atol=1e-8, retry_maxiter=7, retry_time_limit=20.)
    assert batch.attempt_budget(config, 1) == (1, 10)
    assert batch.attempt_budget(config, 2) == (7, 20.)
    parent = tmp_path / 'attempts' / task['id']
    interrupted = parent / 'attempt-0001'
    interrupted.mkdir(parents=True)
    (interrupted / 'job.log').write_text('interrupted; must remain unchanged')
    for number, iterations in ((2, 1), (3, 6000)):
        directory = parent / f'attempt-{number:04d}'
        directory.mkdir()
        P, report, _ = opt.solve_problem(problem, .05 * intensity.astype(complex),
                    gap_atol=1e-8, maxiter=iterations, time_limit=10, check_every=10)
        assert report['numerically_certified'] == (number == 3)
        batch.atomic(directory / 'certificate.json', dict(solver=report, gap_atol=1e-8))
        batch.atomic(directory / 'taskreceipt.json',
            dict(task=task, identity='identity', status='certified' if number == 3 else 'uncertified',
                 scientifically_certified=False,
                 outputs={'certificate.json': batch.sha(directory / 'certificate.json')}))
        completed, numerical, errors = batch.scan_task(tmp_path, task, 'identity')
        assert numerical == 1 and errors == 0
        assert len(completed) == (number == 3)
        assert np.all(abs(P) <= intensity)
    assert (interrupted / 'job.log').read_text() == 'interrupted; must remain unchanged'
    with pytest.raises(ValueError, match='identity mismatch'):
        batch.scan_task(tmp_path, task, 'different')
    (parent / 'attempt-0003' / 'certificate.json').write_text('{}')
    with pytest.raises(ValueError, match='hash mismatch'):
        batch.scan_task(tmp_path, task, 'identity')
    assert 'import numpy' not in SCRIPT.read_text().split('def make_draws')[0]


def test_pool_memory_reservation_and_idle_timeout():
    gib = 2**30
    before = batch.memory_projection(25*gib, 0, 64)
    assert before['unmaterialized_reserve_bytes'] == 64*gib
    assert before['projected_memory_bytes'] == 89*gib
    # RSS already in account usage is subtracted from the reservation, not doubled.
    halfway = batch.memory_projection(57*gib, 32*gib, 64)
    full = batch.memory_projection(89*gib, 64*gib, 64)
    assert halfway['projected_memory_bytes'] == full['projected_memory_bytes'] == 89*gib
    over = batch.memory_projection(105*gib, 80*gib, 64)
    assert over['unmaterialized_reserve_bytes'] == 0
    assert over['projected_memory_bytes'] > 100*gib
    assert batch.memory_projection(40*gib, 0, 64)['projected_memory_bytes'] > 100*gib
    # Physical admission also needs all unmaterialized memory plus 2 GiB.
    assert 65*gib < before['unmaterialized_reserve_bytes'] + 2*gib
    since, expired = batch.resource_wait_expired(None, 10., 0, 12, False, 300.)
    assert (since, expired) == (10., False)
    assert batch.resource_wait_expired(since, 309., 0, 12, False, 300.) == (10., False)
    assert batch.resource_wait_expired(since, 310., 0, 12, False, 300.) == (10., True)
    for active, pending, admitted in ((1, 12, False), (0, 0, False), (0, 12, True)):
        assert batch.resource_wait_expired(since, 1000., active, pending, admitted, 300.) == (None, False)
