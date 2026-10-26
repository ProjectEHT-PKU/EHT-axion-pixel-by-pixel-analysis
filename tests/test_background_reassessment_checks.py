"""Regression checks for embargo and leakage-sensitive predictive scoring."""
import numpy as np

from check_background_embargo import embargo_training
from eht_axion.research.background_validation import evaluate_folds, purged_folds


def test_extra_embargo_never_restores_a_purged_raw_frame():
    times = np.arange(141)*8.893247544352288
    starts = times[:124]+3
    for fold in purged_folds(starts, times, n_splits=5):
        train = embargo_training(fold['support_indices'], fold['train'], fold['test'], times, 24.)
        assert set(train) <= set(fold['train'])
        lo, hi = times[fold['test_support']]
        supports = fold['support_indices'][train]
        assert np.all((times[supports[:, 1]] < lo-24) | (times[supports[:, 0]] > hi+24))
        assert len(train) < len(fold['train'])


def test_test_only_shift_is_not_silently_centered_out():
    rng = np.random.default_rng(17)
    times = np.arange(141)*8.893247544352288
    starts = times[:124]+3
    data = rng.normal(0, .02, (124, 3, 2))
    folds = purged_folds(starts, times, n_splits=5)
    original = evaluate_folds(data, starts, times, [0, 1, 2])
    changed = data.copy()
    changed[folds[0]['test']] += .4
    altered = evaluate_folds(changed, starts, times, [0, 1, 2])
    before = next(r for r in original['records'] if r['fold'] == 0 and r['model'] == 'iso')
    after = next(r for r in altered['records'] if r['fold'] == 0 and r['model'] == 'iso')
    assert after['mean_nll_per_dimension'] > before['mean_nll_per_dimension']+50
    np.testing.assert_array_equal(after['train_mean'], before['train_mean'])
