import numpy as np
import pytest

from eht_axion.statistics import normalized_posterior_cdf, posterior_upper_limit


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


@pytest.mark.parametrize("slope", [0.0, 0.7, -0.3])
@pytest.mark.parametrize("credibility", [0.05, 0.5, 0.95])
def test_coarse_grid_quantile_inverts_integrated_linear_density(slope, credibility):
    grid = np.array([1.0, 1.4, 3.0])
    density = 1 + slope * (grid - 1)
    limit = posterior_upper_limit(grid, np.log(density), prior="uniform",
                                  credibility=credibility)
    width = limit - 1
    mass = width + slope * width**2 / 2
    normalization = 2 + slope * 2
    assert mass / normalization == pytest.approx(credibility, abs=2e-15)


def test_quantile_handles_zero_density_plateaus():
    grid = np.arange(1.0, 8.0)
    loglike = np.array([-np.inf, -np.inf, 0, -np.inf, -np.inf, 0, -np.inf])
    limit = posterior_upper_limit(grid, loglike, prior="uniform", credibility=0.5)
    assert limit == pytest.approx(4.0)
    density, _ = normalized_posterior_cdf(grid, loglike, prior="uniform")
    assert density[0] == 0


@pytest.mark.parametrize("grid, logs", [
    ([1.0], [0.0]), ([1.0, np.inf], [0.0, 0.0]),
    ([1.0, np.nan], [0.0, 0.0]), ([1.0, 1.0], [0.0, 0.0]),
    ([0.0, 1.0], [0.0, 0.0]), ([2.0, 1.0], [0.0, 0.0]),
    ([1.0, 2.0], [0.0, np.nan]), ([1.0, 2.0], [0.0, np.inf]),
    ([1.0, 2.0], [-np.inf, -np.inf]), ([1.0, 2.0], [[0.0, 0.0]]),
])
def test_posterior_shape_and_finite_guards(grid, logs):
    with pytest.raises(ValueError):
        normalized_posterior_cdf(grid, logs)
