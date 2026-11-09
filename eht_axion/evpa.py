"""EVPA of linear polarization and pi-periodic angle differences (shared by eht/ and analysis/)."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def finite_evpa(stokes_q: ArrayLike, stokes_u: ArrayLike) -> NDArray[np.float64]:
    """Return EVPA in radians (period pi), or NaN for invalid/zero Q+iU.

    Q and U broadcast under NumPy rules and are converted to float64. Neither
    intensity nor V participates in this observable. Test Q/U directly rather
    than forming P squared, which can overflow or underflow for finite inputs.
    """
    q, u = np.broadcast_arrays(
        np.asarray(stokes_q, dtype=np.float64), np.asarray(stokes_u, dtype=np.float64)
    )
    valid = np.isfinite(q) & np.isfinite(u) & ((q != 0) | (u != 0))
    return np.where(valid, 0.5 * np.arctan2(u, q), np.nan)


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
