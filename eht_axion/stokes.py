"""Shared validity rules for linear-polarization observables."""

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
