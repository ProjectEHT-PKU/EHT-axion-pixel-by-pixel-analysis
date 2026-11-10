"""RAPTOR tasks with a positive axion norm: identity, binding and drift checks; synthetic inputs only."""
from copy import deepcopy

import pytest

from eht_axion import provenance
from eht_axion.experiment import read_experiment, validate
from eht_axion.raptor.jobs import _prepare_task, _validate_plan_binding, serialize_parameters
from eht_axion.raptor.plan import plan_generation
from test_raptor_plan import write_snapshot, write_table
from test_raptor_jobs import bound_plan, setup  # noqa: F401 -- shared synthetic fixtures

@pytest.mark.parametrize('norm', [.05, .1, .2, .4, .8, 1.6, 2.0])
def test_positive_finite_norm_and_explicit_model_binding(tmp_path, norm):
    c = read_experiment()
    c['raptor']['snapshots'] = {'start': 2300, 'stop': 2300, 'step': 10}
    c['raptor']['response']['reference_norm'] = norm
    c['shared']['alphas'] = [.3]
    c['raptor']['alpha_table'] = {'0.30': 'low'}
    write_snapshot(tmp_path / 'data002300.dat', time=15080)
    table = tmp_path / 'radial.csv'
    write_table(table)
    p = plan_generation(c, tmp_path, {'low': table, 'high': table})
    assert p['config']['raptor']['response']['reference_norm'] == norm
    for t in p['tasks']:
        expected = norm if t['kind'] == 'response' else 0.
        assert t['norm'] == t['parameters']['AXION_NORM']['value'] == expected
        assert float(serialize_parameters(t['parameters']).split('AXION_NORM')[1].split()[1]) == expected
    default = read_experiment()
    assert default['raptor']['response']['reference_norm'] == .1


@pytest.mark.parametrize('norm', [0, -1, True, False, '.4', None, float('nan'), float('inf')])
def test_invalid_norm_still_rejected(norm):
    c = read_experiment()
    c['raptor']['response']['reference_norm'] = norm
    with pytest.raises(ValueError):
        validate(c)


def test_norm_changes_response_not_baseline_identity(tmp_path):
    c = read_experiment()
    c['raptor']['snapshots'] = {'start': 1600, 'stop': 1600, 'step': 10}
    write_snapshot(tmp_path / 'data001600.dat')
    table = tmp_path / 'radial.csv'
    write_table(table)
    p = plan_generation(c, tmp_path, {'low': table, 'high': table})
    c['raptor']['response']['reference_norm'] = .8
    q = plan_generation(c, tmp_path, {'low': table, 'high': table})
    assert [t['id'] for t in p['tasks'] if t['kind'] == 'baseline'] == [
        t['id'] for t in q['tasks'] if t['kind'] == 'baseline']
    assert set(t['id'] for t in p['tasks'] if t['kind'] == 'response').isdisjoint(
        t['id'] for t in q['tasks'] if t['kind'] == 'response')


@pytest.mark.parametrize('drift', ['config', 'model', 'task-and-model'])
def test_norm_semantic_drift_still_fails_plan_binding(bound_plan, drift):  # noqa: F811
    fixture, p = bound_plan
    t = deepcopy(next(t for t in p['tasks'] if t['kind'] == 'response'))
    p['tasks'] = [t]
    if drift == 'config':
        p['config']['raptor']['response']['reference_norm'] = .8
    elif drift == 'model':
        t['parameters']['AXION_NORM']['value'] = .8
    else:
        t['norm'] = .8
        t['parameters']['AXION_NORM']['value'] = .8
    with provenance.verification_session(), pytest.raises(ValueError):
        prepared = [_prepare_task(t)]
        m = provenance.verify_stage(fixture['build'], recursive=True)
        _validate_plan_binding(p, prepared, m, fixture['build'], fixture['build'] / 'grid_mks.in')
    assert not fixture['run'].exists()
