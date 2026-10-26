"""Synthetic background-only oracles; never load EHT or run production science."""
import numpy as np
import pytest

from eht_axion.research.background_validation import (
    MAIN_MODELS, DIAGNOSTIC_MODELS, evaluate_folds, model_covariance,
    purged_folds, support_indices,
)


def test_interpolation_brackets_and_exact_knots():
    np.testing.assert_array_equal(
        support_indices([2., 4.], np.arange(10.), [0.], width=2.),
        [[1, 3], [3, 5]])
    np.testing.assert_array_equal(
        support_indices([2.2, 4.2], np.arange(10.), [0.], width=2.),
        [[1, 4], [3, 6]])
    with pytest.raises(ValueError, match='outside'):
        support_indices([0.], np.arange(10.), [0.], width=2.)


def test_shared_single_raw_knot_is_purged_even_without_exposure_overlap():
    starts = np.array([1.2, 2.8, 4.2, 5.8, 7.2, 8.8])
    folds = purged_folds(starts, np.arange(11.), [0.], .2, 3)
    for fold in folds:
        s = fold['support_indices']
        test_raw = set().union(*(set(range(a, b+1)) for a, b in s[fold['test']]))
        for i in fold['train']:
            assert test_raw.isdisjoint(range(s[i, 0], s[i, 1]+1))
    # Nonadjacent, genuinely disjoint support is retained.
    assert 0 in folds[1]['train']  # nonadjacent support retained
    toy = purged_folds([1.2, 1.8, 4.2, 5.8], np.arange(8.), [0.], .1, 4)
    assert 1 in toy[0]['purged']  # both use knots 1 and 2


def test_actual_geometry_and_prepared_centers_equivalence():
    cadence = 8.89324754
    times = np.arange(141)*cadence
    starts = 3 + np.arange(124)*cadence
    days = np.array([0., 24., 120., 144.])
    a = support_indices(starts, times)
    b = support_indices(starts, times, window_centers=starts[:, None]+days)
    np.testing.assert_array_equal(a, b)
    assert np.all(a[:, 1]-a[:, 0] == 17)
    folds = purged_folds(starts, times)
    assert [len(f['test']) for f in folds] == [25, 25, 25, 25, 24]
    assert [len(f['train']) for f in folds] == [82, 65, 65, 65, 83]
    assert 123 in folds[0]['train']  # explicitly two-sided, not forecasting
    sensitivity = purged_folds(starts, times, n_splits=3)
    assert [len(f['test']) for f in sensitivity] == [42, 41, 41]
    assert [len(f['train']) for f in sensitivity] == [65, 49, 66]


def test_covariance_formulas_and_centered_population_convention():
    b = np.random.default_rng(2).normal(size=(40, 3, 4)) + 7
    x = (b-b.mean(axis=0)).reshape(40, 12)
    full = x.T@x/40
    t = np.einsum('mpn,mqn->pq', b-b.mean(axis=0), b-b.mean(axis=0))/160
    expected = {
        'iso': np.eye(12)*np.trace(full)/12,
        'temporal_iso': np.kron(t, np.eye(4)),
        'partial_diag': np.diag(.5*np.diag(full)+.5*np.trace(full)/12),
        'diagonal': np.diag(np.diag(full)), 'full': full,
        'shrink_full_0.25': .75*np.kron(t, np.eye(4))+.25*full,
    }
    for name, want in expected.items():
        cov, meta = model_covariance(b, name)
        np.testing.assert_allclose(cov, want)
        assert meta['scoreable'] and meta['predictive_mean'] == 'zero'
        np.testing.assert_allclose(meta['train_mean'], b.mean(axis=0))


def test_singular_full_is_unscoreable_without_rescue():
    b = np.random.default_rng(3).normal(size=(5, 3, 4))
    full, meta = model_covariance(b, 'full')
    assert meta['rank'] <= 4 and not meta['scoreable']
    assert meta['empirical_rank_bound'] == 4
    assert model_covariance(b, 'shrink_full_0.25')[1]['scoreable']
    for name in MAIN_MODELS + DIAGNOSTIC_MODELS:
        cov, meta = model_covariance(np.ones((10, 3, 4)), name)
        assert not meta['scoreable']
        np.testing.assert_array_equal(cov, 0)


def toy_evaluation(b):
    return evaluate_folds(b, np.arange(len(b))+1., np.arange(len(b)+3.),
                          [0, 2], day_hours=[0.], width=.2)


def test_zero_mean_gaussian_score_oracle_and_weighting():
    b = np.random.default_rng(7).normal(size=(31, 3, 2)) + 2
    result = toy_evaluation(b)
    assert result['selected_model'] is None
    scores, sizes = [], []
    for r in result['records']:
        if r['model'] != 'iso':
            continue
        fold = result['folds'][r['fold']]
        cov, _ = model_covariance(b[fold['train']][:, [0, 2]], 'iso')
        test = b[fold['test']][:, [0, 2]].reshape(r['n_test'], -1)
        q = np.einsum('mi,ij,mj->m', test, np.linalg.inv(cov), test)
        expected = .5*(4*np.log(2*np.pi)+np.linalg.slogdet(cov)[1]+q)/4
        np.testing.assert_allclose(r['nll_per_dimension'], expected)
        scores.append(expected.mean())
        sizes.append(len(expected))
    assert result['summary']['iso']['weighted_mean_nll_per_dimension'] == pytest.approx(
        np.average(scores, weights=sizes))


def test_heldout_values_never_change_training_moments():
    b = np.random.default_rng(12).normal(size=(30, 3, 2))
    first = toy_evaluation(b)
    changed = b.copy()
    changed[first['folds'][0]['test']] += 1e4
    second = toy_evaluation(changed)
    for a, c in zip(first['records'][:6], second['records'][:6]):
        np.testing.assert_array_equal(a['train_mean'], c['train_mean'])
        assert a['condition_number'] == c['condition_number']
        assert c['mean_nll_per_dimension'] > a['mean_nll_per_dimension']


def test_no_silent_aggregate_of_only_scoreable_folds():
    b = np.random.default_rng(4).normal(size=(20, 3, 10))
    result = toy_evaluation(b)
    assert not result['summary']['full']['all_folds_scoreable']
    assert 'weighted_mean_nll_per_dimension' not in result['summary']['full']
    assert result['summary']['iso']['all_folds_scoreable']


@pytest.mark.parametrize('bad', [np.zeros((1, 2, 3)), np.zeros((3, 2)),
                                  np.full((3, 2, 1), np.nan),
                                  np.ones((3, 2, 1), dtype=complex)])
def test_bad_samples_rejected(bad):
    with pytest.raises(ValueError):
        model_covariance(bad, 'iso')
