"""Vectorized circular statistics used by the EVPA analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def circular_difference(
    first: ArrayLike,
    second: ArrayLike,
    *,
    period: float = np.pi,
) -> NDArray[np.floating] | np.floating:
    """Return ``first - second`` wrapped to ``[-period/2, period/2)``.

    EVPA has period pi rather than 2*pi.  The implementation broadcasts and
    supports both scalar and array inputs, unlike the historical ``CircleCalc``
    helper.
    """
    if not np.isfinite(period) or period <= 0:
        raise ValueError("period must be finite and positive")
    result = (np.asarray(first) - np.asarray(second) + period / 2) % period - period / 2
    return result[()] if result.ndim == 0 else result


def circular_mean(
    angles: ArrayLike,
    *,
    axis: int | tuple[int, ...] | None = None,
    period: float = np.pi,
) -> NDArray[np.floating] | np.floating:
    """Return the circular mean while ignoring NaNs."""
    angles_array = np.asarray(angles)
    phase = angles_array * (2 * np.pi / period)
    vector = np.nanmean(np.exp(1j * phase), axis=axis)
    result = np.angle(vector) * (period / (2 * np.pi))
    return result[()] if np.ndim(result) == 0 else result


def unwrap_around(
    angles: ArrayLike,
    center: ArrayLike,
    *,
    period: float = np.pi,
) -> NDArray[np.floating] | np.floating:
    """Map circular values to the continuous interval centered on ``center``."""
    result = np.asarray(center) + circular_difference(angles, center, period=period)
    return result[()] if np.ndim(result) == 0 else result
