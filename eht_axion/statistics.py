"""Numerically stable likelihood and posterior utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import cho_factor, cho_solve
from scipy.special import logsumexp


def gaussian_log_likelihood(
    residual: ArrayLike,
    covariance: ArrayLike,
    *,
    include_logdet: bool = True,
) -> float:
    """Return a multivariate Gaussian log-likelihood up to ``n*log(2*pi)``.

    The log determinant must be retained whenever covariance depends on a model
    parameter.  Cholesky factorization both validates positive definiteness and
    avoids explicitly constructing an inverse.
    """
    residual_array = np.asarray(residual, dtype=float)
    covariance_array = np.asarray(covariance, dtype=float)
    if covariance_array.shape != (residual_array.size, residual_array.size):
        raise ValueError("covariance shape does not match residual")
    factor, lower = cho_factor(covariance_array, lower=True, check_finite=False)
    solved = cho_solve((factor, lower), residual_array, check_finite=False)
    value = float(residual_array @ solved)
    if include_logdet:
        value += float(2 * np.log(np.diag(factor)).sum())
    return -0.5 * value


def combine_independent_log_likelihoods(log_likelihoods: ArrayLike) -> float:
    """Combine independent observations by multiplying their likelihoods."""
    return float(np.asarray(log_likelihoods, dtype=float).sum())


def marginalize_log_likelihood(log_likelihoods: ArrayLike, *, axis: int = -1):
    """Average likelihood over a uniformly sampled nuisance parameter."""
    values = np.asarray(log_likelihoods, dtype=float)
    return logsumexp(values, axis=axis) - np.log(values.shape[axis])


def normalized_posterior_cdf(
    grid: ArrayLike,
    log_likelihood: ArrayLike,
    *,
    prior: str = "log-uniform",
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    """Return normalized posterior density and CDF on an arbitrary grid.

    ``prior='log-uniform'`` means p(g) proportional to 1/g.  Integration is
    performed with trapezoids in g, so prior choice is explicit rather than an
    accidental consequence of a logarithmic grid.
    """
    grid_array = np.asarray(grid, dtype=float)
    log_like = np.asarray(log_likelihood, dtype=float)
    if grid_array.ndim != 1 or log_like.shape != grid_array.shape:
        raise ValueError("grid and log_likelihood must be matching 1-D arrays")
    if np.any(grid_array <= 0) or np.any(np.diff(grid_array) <= 0):
        raise ValueError("grid must be positive and strictly increasing")
    if prior == "log-uniform":
        log_prior = -np.log(grid_array)
    elif prior == "uniform":
        log_prior = np.zeros_like(grid_array)
    else:
        raise ValueError("prior must be 'uniform' or 'log-uniform'")
    log_density = log_like + log_prior
    density = np.exp(log_density - np.max(log_density))
    increments = np.diff(grid_array) * (density[:-1] + density[1:]) / 2
    normalization = increments.sum()
    if not np.isfinite(normalization) or normalization <= 0:
        raise ValueError("posterior cannot be normalized on this grid")
    density /= normalization
    cdf = np.concatenate(([0.0], np.cumsum(increments) / normalization))
    cdf[-1] = 1.0
    return density, cdf


def posterior_upper_limit(
    grid: ArrayLike,
    log_likelihood: ArrayLike,
    *,
    credibility: float = 0.95,
    prior: str = "log-uniform",
) -> float:
    """Interpolate a posterior upper limit from its normalized CDF."""
    if not 0 < credibility < 1:
        raise ValueError("credibility must lie strictly between zero and one")
    grid_array = np.asarray(grid, dtype=float)
    _, cdf = normalized_posterior_cdf(grid_array, log_likelihood, prior=prior)
    return float(np.interp(credibility, cdf, grid_array))
