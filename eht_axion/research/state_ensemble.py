"""Small, I/O-free adapters for the four-day research input study.

Local chronological lifting chooses minimal consecutive EVPA differences; it
is a branch convention, not recovery of physical winding across 96 hours.
"""

from __future__ import annotations

import numpy as np

from ..pipeline.raptor import _integrate
from .time_sampling import DAY_HOURS, average_modulated_response, evpa_contrasts


def selected_background(days, indices, dimension, operator=None, weighting="intensity"):
    """Project already exposure-averaged Stokes, retaining existing RAPTOR axes."""
    days = np.asarray(days)
    if days.ndim == 5 and days.shape[1] == 1:
        days = days[:, 0]
    if days.ndim != 4 or days.shape[0] != 4 or days.shape[-1] != 4:
        raise ValueError("need (4,[1,]X,Y,4) averaged Stokes")
    if dimension == "1d":
        values = _integrate(days, operator, weighting)
    elif dimension == "2d":
        q, u = days[..., 1], days[..., 2]
        values = np.where(np.hypot(q, u) > 0, 0.5 * np.arctan2(u, q), np.nan)
    else:
        raise ValueError("dimension must be 1d or 2d")
    result = values.reshape(4, -1)[:, indices]
    if not np.isfinite(result).all():
        raise ValueError(
            "nonfinite selected background: no exclusion or imputation allowed"
        )
    return result


def local_branches(daily):
    """Retain raw and principal edges plus a consistently lifted local branch."""
    daily = np.asarray(daily, dtype=float)
    diagnostics = evpa_contrasts(daily)
    lifted = np.concatenate(
        (
            daily[..., :1, :],
            daily[..., :1, :] + np.cumsum(diagnostics["principal"], axis=-2),
        ),
        axis=-2,
    )
    return lifted, diagnostics


def response_days(times, coefficients, starts, omega, exposure=6.0):
    """Frozen and evolving envelopes; phase zero at first day center.

    coefficients already include reference normalization and cloud scaling.
    Frozen states use every native snapshot, not just valid evolving starts.
    """
    times, coefficients, starts = map(np.asarray, (times, coefficients, starts))
    if coefficients.ndim != 2 or len(coefficients) != len(times):
        raise ValueError("need (snapshot,selected feature) complex coefficients")
    if not np.isfinite(coefficients).all():
        raise ValueError("nonfinite selected response: no imputation allowed")
    frozen = (
        coefficients[:, None, :]
        * (np.exp(1j * omega * DAY_HOURS) * np.sinc(omega * exposure / (2 * np.pi)))[
            None, :, None
        ]
    )
    evolving = average_modulated_response(
        starts[:, None] + DAY_HOURS, times, coefficients, omega, exposure
    )
    evolving *= np.exp(-1j * omega * starts)[:, None, None]
    return frozen, evolving
