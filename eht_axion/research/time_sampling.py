"""Physical-hour, piecewise-linear research sampling; no production I/O.

Exposure integration is exact for the *interpolant*, not unresolved GRMHD.
Stokes averaging must precede aperture projection/EVPA extraction.
"""
from __future__ import annotations

import numpy as np

DAY_HOURS = np.array([0., 24., 120., 144.])
H3 = np.array([[-1., 1., 0., 0.], [0., -1., 1., 0.], [0., 0., -1., 1.]])
H2 = H3[[0, 2]].copy()


def _times_values(times, values):
    t = np.asarray(times, dtype=float)
    v = np.asarray(values)  # preserves memmaps; never copies the full cube
    if (t.ndim != 1 or len(t) < 2 or not np.isfinite(t).all()
            or np.any(np.diff(t) <= 0) or v.ndim < 1 or len(v) != len(t)):
        raise ValueError('need increasing finite times and matching leading value axis')
    return t, v


def _bounds(t, centers, width):
    c = np.asarray(centers, dtype=float)
    if (not np.isfinite(width) or width < 0 or not np.isfinite(c).all()
            or np.any(c - width / 2 < t[0]) or np.any(c + width / 2 > t[-1])):
        raise ValueError('finite nonnegative exposure must lie wholly inside time range')
    return c


def _at(t, v, x):
    j = min(np.searchsorted(t, x, side='right') - 1, len(t) - 2)
    f = (x - t[j]) / (t[j + 1] - t[j])
    return (1 - f) * v[j] + f * v[j + 1]


def _average(times, values, centers, width, omega):
    t, v = _times_values(times, values)
    c = _bounds(t, centers, width)
    if not np.isfinite(omega):
        raise ValueError('omega must be finite radians/hour')
    dtype = np.result_type(v.dtype, np.float64, complex if omega else float)
    out = np.empty((c.size,) + v.shape[1:], dtype=dtype)
    for k, center in enumerate(c.ravel()):
        if width == 0:
            out[k] = _at(t, v, center) * np.exp(1j * omega * center) if omega else _at(t, v, center)
            continue
        lo, hi = center - width / 2, center + width / 2
        knots = np.r_[lo, t[(t > lo) & (t < hi)], hi]
        total = np.zeros(v.shape[1:], dtype=dtype)
        for a, b in zip(knots[:-1], knots[1:]):
            va, vb = _at(t, v, a), _at(t, v, b)
            h, mid = (b - a) / 2, (a + b) / 2
            x = omega * h
            # sinc(x) and its derivative with respect to x; stable at omega=0.
            sinc = np.sinc(x / np.pi)
            derivative = (-x / 3 + x**3 / 30 - x**5 / 840 if abs(x) < 1e-3
                          else (x * np.cos(x) - np.sin(x)) / x**2)
            segment = (b - a) * (va + vb) / 2 * sinc
            if omega:
                segment = np.exp(1j * omega * mid) * (
                    segment - 1j * h * (vb - va) * derivative)
            total += segment
        out[k] = total / width
    return out.reshape(c.shape + v.shape[1:])


def integrate_stokes_at(times, stokes, centers, width=6.):
    """Exact rectangular average of linear Stokes interpolant, no extrapolation.

    All times/widths in hours; shape centers.shape + stokes.shape[1:]. Spatial
    and channel axes are untouched (including singleton phase/state0 axis).
    Only the current exposure's endpoint frames are materialized.
    """
    return _average(times, stokes, centers, width, 0.)


def average_modulated_response(centers, coefftimes, complexcoeff, omega, width=6.):
    """Average z(t)*exp(i*omega*t), z=A*exp(i*phase), exactly per linear segment.

    omega is radians/hour, NOT dimensionless omega(alpha). A phase realization
    is real(output*exp(i*psi)). Caller supplies fit normalization/cloud scaling;
    this function deliberately does not normalize, select columns or add noise.
    """
    return _average(coefftimes, complexcoeff, centers, width, omega)


def four_day_stokes(times, stokes, width=6., starts=None):
    """Yield (start_hour, averaged_days) one window at a time.

    Default starts are times[0]+width/2+k*native cadence; irregular inputs use
    their own times+width/2. Overlapping windows are model states, not iid data.
    Explicit starts are strictly checked; invalid windows are never truncated.
    """
    t, v = _times_values(times, stokes)
    if not np.isfinite(width) or width < 0:
        raise ValueError('width must be finite and nonnegative')
    if starts is None:
        starts = t + width / 2
        starts = starts[starts + DAY_HOURS[-1] + width / 2 <= t[-1]]
    starts = np.atleast_1d(np.asarray(starts, dtype=float))
    if starts.ndim != 1:
        raise ValueError('starts must be one dimensional')
    _bounds(t, starts[:, None] + DAY_HOURS, width)
    for start in starts:
        yield float(start), integrate_stokes_at(t, v, start + DAY_HOURS, width)


def contrast_covariance(day_covariances, h=H3):
    """Conditional-independent day errors; flattened order (contrast, spatial).

    Return (K*N,K*N), including negative shared-day off-diagonal blocks.
    This is observational covariance only; never add model-state scatter here.
    """
    c, h = np.asarray(day_covariances), np.asarray(h, dtype=float)
    if (c.ndim != 3 or c.shape[1] != c.shape[2] or h.ndim != 2
            or h.shape[1] != c.shape[0] or not np.isfinite(c).all()
            or not np.isfinite(h).all()):
        raise ValueError('need finite (day,N,N) covariances and (contrast,day) H')
    return np.einsum('ad,dij,bd->aibj', h, c, h).reshape(
        h.shape[0] * c.shape[1], h.shape[0] * c.shape[1])


def evpa_contrasts(chi_days):
    """Last axes (4,N); return raw/principal H3, winding and six-edge closure.

    No chronological unwrap across long gaps and no automatic sample exclusion.
    Raw angles must also be retained by callers. Linear Gaussian equivalence
    needs a consistent local branch; principal contrasts alone do not prove it.
    """
    chi = np.asarray(chi_days, dtype=float)
    if chi.ndim < 2 or chi.shape[-2] != 4:
        raise ValueError('need (...,4,N) EVPA radians')
    raw = np.einsum('ad,...dn->...an', H3, chi)
    principal = (raw + np.pi / 2) % np.pi - np.pi / 2
    winding = np.rint((raw - principal) / np.pi)
    closure = []
    for i in range(4):
        for j in range(i + 1, 4):
            direct = (chi[..., j, :] - chi[..., i, :] + np.pi / 2) % np.pi - np.pi / 2
            closure.append(np.rint((principal[..., i:j, :].sum(axis=-2) - direct) / np.pi))
    return dict(raw=raw, principal=principal, winding=winding,
                closure_winding=np.stack(closure, axis=-2))
