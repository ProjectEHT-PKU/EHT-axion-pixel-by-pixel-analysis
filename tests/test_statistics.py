import numpy as np
import pytest

from eht_axion.statistics import (
    combine_independent_log_likelihoods,
    gaussian_log_likelihood,
    normalized_posterior_cdf,
    posterior_upper_limit,
)


def test_gaussian_log_likelihood_includes_parameter_dependent_logdet():
    residual = np.array([1.0, 2.0])
    covariance = np.diag([4.0, 9.0])
    expected = -0.5 * (1 / 4 + 4 / 9 + np.log(4 * 9))
    assert gaussian_log_likelihood(residual, covariance) == pytest.approx(expected)


def test_independent_likelihoods_multiply_in_log_space():
    assert combine_independent_log_likelihoods([-2, -3]) == -5


def test_log_uniform_constant_likelihood_is_uniform_in_log():
    grid = np.logspace(-2, 2, 10001)
    limit = posterior_upper_limit(grid, np.zeros_like(grid), prior="log-uniform")
    assert np.log10(limit) == pytest.approx(1.8, abs=2e-3)


def test_uniform_prior_constant_likelihood():
    grid = np.linspace(1, 101, 10001)
    limit = posterior_upper_limit(grid, np.zeros_like(grid), prior="uniform")
    assert limit == pytest.approx(96.0, abs=1e-3)


def test_cdf_is_normalized_and_monotonic():
    grid = np.logspace(-2, 2, 101)
    _, cdf = normalized_posterior_cdf(grid, -grid, prior="log-uniform")
    assert cdf[0] == 0
    assert cdf[-1] == 1
    assert np.all(np.diff(cdf) >= 0)
