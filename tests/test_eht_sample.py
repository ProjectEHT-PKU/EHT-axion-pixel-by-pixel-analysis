"""Synthetic snapshots and tiny exact PDHG only; no real observations/MC launch."""
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.eht import sample as batch

SCRIPT = Path(batch.__file__)


def sources(tmp_path):
    result = {}
    for day in batch.DAYS:
        p = tmp_path / f'{day}-final-dterms.json'
        terms = {s: [[.01 * int(day), .002], [.003, -.004 * int(day)]] for s in ('AZ', 'LM', 'PV')}
        p.write_text(json.dumps(terms))
        result[day] = dict(I=str(tmp_path / f'{day}.fits'), obs=str(tmp_path / f'{day}.uvfits'),
                           dterms=str(p))
    return result


def test_settings_seeds_and_4000_unique_tasks(tmp_path):
    inputs = sources(tmp_path)
    science = batch.science_settings(SCRIPT.parents[2] / 'configs' / 'experiment.yaml')
    assert science['sites'] == ['AA', 'AP', 'AZ', 'JC', 'LM', 'PV', 'SM']
    assert (science['draws_per_day'], science['seed'], science['fsys']) == (1000, 20210405, .005)
    assert science['feed_rotation']['SM'] == [1, -1, 45]
    assert science['monte_carlo'] == dict(initial_p_fraction=.05, gap_atol=1e-6, maxiter=50000,
                                          check_every=500, time_limit=600., retry_maxiter=200000,
                                          retry_time_limit=1800.)
    draws = batch.make_draws(inputs, 1000, 20210405, science)
    tasks = batch.enumerate_tasks(draws)
    assert len(tasks) == len({t['id'] for t in tasks}) == 4000
    assert all(sum(t['day'] == d for t in tasks[:64]) == 16 for d in batch.DAYS)
    assert {t['day'] for t in tasks} == set(batch.DAYS)
    assert all(len(draws[d]['draws']) == 1000 for d in batch.DAYS)
    assert all(set(t['dterms']) == set(science['sites']) for t in tasks)
    repeated = batch.make_draws(dict(reversed(list(inputs.items()))), 1000, 20210405, science)
    assert repeated == draws
    assert len({draws[d]['seed'] for d in batch.DAYS}) == 4
    assert draws['95']['draws'][0] != draws['96']['draws'][0]
    with pytest.raises(ValueError, match='Unsafe'):
        batch.contained(tmp_path, '../escape')


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
            dict(task=task, status='certified' if number == 3 else 'uncertified',
                 scientifically_certified=False, outputs=['certificate.json']))
        completed, numerical, errors = batch.scan_task(tmp_path, task)
        assert numerical == 1 and errors == 0
        assert len(completed) == (number == 3)
        assert np.all(abs(P) <= intensity)
    assert (interrupted / 'job.log').read_text() == 'interrupted; must remain unchanged'
    with pytest.raises(ValueError, match='task mismatch'):
        batch.scan_task(tmp_path, dict(task, draw=1))
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
