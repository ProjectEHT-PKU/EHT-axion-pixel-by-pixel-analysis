"""Numerically stable likelihood and posterior utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _validated_log_likelihood(values: ArrayLike) -> NDArray[np.floating]:
    result = np.asarray(values, dtype=float)
    if result.size == 0 or np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError("log likelihood must be nonempty, with no NaN or +inf")
    return result


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
    log_like = _validated_log_likelihood(log_likelihood)
    if grid_array.ndim != 1 or grid_array.size < 2 or log_like.shape != grid_array.shape:
        raise ValueError("grid and log_likelihood must be matching 1-D arrays")
    if (not np.all(np.isfinite(grid_array)) or np.any(grid_array <= 0)
            or np.any(np.diff(grid_array) <= 0)):
        raise ValueError("grid must be finite, positive and strictly increasing")
    if not np.any(np.isfinite(log_like)):
        raise ValueError("posterior cannot be normalized from all -inf likelihoods")
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
    cdf = np.minimum(np.concatenate(([0.0], np.cumsum(increments) / normalization)), 1.0)
    cdf[-1] = 1.0
    return density, cdf


def posterior_upper_limit(
    grid: ArrayLike,
    log_likelihood: ArrayLike,
    *,
    credibility: float = 0.95,
    prior: str = "log-uniform",
) -> float:
    """Invert the same piecewise-linear density used by trapezoidal integration.

    Linear interpolation of CDF nodes instead assumes piecewise-constant
    density and is inconsistent, particularly on a coarse logarithmic grid.
    """
    if not 0 < credibility < 1:
        raise ValueError("credibility must lie strictly between zero and one")
    grid_array = np.asarray(grid, dtype=float)
    density, cdf = normalized_posterior_cdf(grid_array, log_likelihood, prior=prior)
    upper = int(np.searchsorted(cdf, credibility, side="left"))
    lower = upper - 1
    width = grid_array[upper] - grid_array[lower]
    # Solve d0*u + (d1-d0)*u**2/2 = local_mass/width. This
    # rationalized root avoids cancellation for nearly constant densities.
    d0, d1 = density[lower], density[upper]
    local = (credibility - cdf[lower]) / width
    scale = max(d0, d1)
    d0, d1, local = d0 / scale, d1 / scale, local / scale
    discriminant = max(0.0, d0**2 + 2 * (d1 - d0) * local)
    denominator = d0 + np.sqrt(discriminant)
    fraction = 0.0 if local == 0 else 2 * local / denominator
    return float(grid_array[lower] + width * np.clip(fraction, 0, 1))
