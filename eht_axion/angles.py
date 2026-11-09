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
