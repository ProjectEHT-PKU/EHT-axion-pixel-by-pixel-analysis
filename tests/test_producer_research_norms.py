"""Positive-norm producer generalization; synthetic inputs, never real RT."""
from copy import deepcopy
from pathlib import Path

import pytest

from eht_axion.pipeline import common
from eht_axion.producer.config import read_production_config, validate_production_config
from eht_axion.producer.jobs import _prepare_task, _validate_plan_binding, serialize_parameters
from eht_axion.producer.plan import plan_generation
from test_producer_plan import write_snapshot, write_table
from test_producer_jobs import production, setup  # noqa: F401 -- shared synthetic fixtures

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('norm', [.05, .1, .2, .4, .8, 1.6, 2.0])
def test_positive_finite_norm_and_explicit_model_binding(tmp_path, norm):
    c = read_production_config(ROOT / 'configs/production.yaml')
    c['producer']['snapshots'] = {'start': 2300, 'stop': 2300, 'step': 10}
    c['producer']['response']['reference_norm'] = norm
    c['data']['alphas'] = [.3]
    c['producer']['alpha_table'] = {'0.30': 'low'}
    write_snapshot(tmp_path / 'data002300.dat', time=15080)
    table = tmp_path / 'radial.csv'
    write_table(table)
    p = plan_generation(c, tmp_path, {'low': table, 'high': table})
    assert p['config']['producer']['response']['reference_norm'] == norm
    for t in p['tasks']:
        expected = norm if t['kind'] == 'response' else 0.
        assert t['norm'] == t['parameters']['AXION_NORM']['value'] == expected
        assert float(serialize_parameters(t['parameters']).split('AXION_NORM')[1].split()[1]) == expected
    default = read_production_config(ROOT / 'configs/production.yaml')
    assert default['producer']['response']['reference_norm'] == .1


@pytest.mark.parametrize('norm', [0, -1, True, False, '.4', None, float('nan'), float('inf')])
def test_invalid_norm_still_rejected(norm):
    c = read_production_config(ROOT / 'configs/production.yaml')
    c['producer']['response']['reference_norm'] = norm
    with pytest.raises(ValueError):
        validate_production_config(c)


def test_norm_changes_response_not_baseline_identity(tmp_path):
    c = read_production_config(ROOT / 'configs/production.yaml')
    c['producer']['snapshots'] = {'start': 1600, 'stop': 1600, 'step': 10}
    write_snapshot(tmp_path / 'data001600.dat')
    table = tmp_path / 'radial.csv'
    write_table(table)
    p = plan_generation(c, tmp_path, {'low': table, 'high': table})
    c['producer']['response']['reference_norm'] = .8
    q = plan_generation(c, tmp_path, {'low': table, 'high': table})
    assert [t['id'] for t in p['tasks'] if t['kind'] == 'baseline'] == [
        t['id'] for t in q['tasks'] if t['kind'] == 'baseline']
    assert set(t['id'] for t in p['tasks'] if t['kind'] == 'response').isdisjoint(
        t['id'] for t in q['tasks'] if t['kind'] == 'response')


@pytest.mark.parametrize('drift', ['config', 'model', 'task-and-model'])
def test_norm_semantic_drift_still_fails_production_binding(production, drift):  # noqa: F811
    fixture, p = production
    t = deepcopy(next(t for t in p['tasks'] if t['kind'] == 'response'))
    p['tasks'] = [t]
    if drift == 'config':
        p['config']['producer']['response']['reference_norm'] = .8
    elif drift == 'model':
        t['parameters']['AXION_NORM']['value'] = .8
    else:
        t['norm'] = .8
        t['parameters']['AXION_NORM']['value'] = .8
    with common.verification_session(), pytest.raises(ValueError):
        prepared = [_prepare_task(t)]
        m = common.verify_stage(fixture['build'], recursive=True)
        _validate_plan_binding(p, prepared, m, fixture['build'], fixture['build'] / 'grid_mks.in')
    assert not fixture['run'].exists()
