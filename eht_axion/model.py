"""Axion EVPA templates reconstructed from saved GRMHD fit products."""

from __future__ import annotations

from typing import Mapping

import numpy as np
from numpy.typing import NDArray

def omega(alpha: float, n_principal: int) -> float:
    """Leading small-alpha M*omega_R; principal N=2 for the scalar 211 mode.

    N = radial_overtone + l + 1, not N+l+1 or n+l+m. This is not
    an exact Kerr eigenfrequency or the interpolated radial-generator table.
    """
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive and finite")
    if isinstance(n_principal, bool) or not isinstance(n_principal, (int, np.integer)) or n_principal < 1:
        raise ValueError("n_principal must be a positive integer")
    return alpha * (1 - alpha**2 / (2 * n_principal**2))


def cloud_amplitude_scaling(
    alpha: float,
    *,
    coefficient: float,
    pivot_alpha: float,
    exponent: float,
) -> float:
    """Cloud amplitude b(alpha) = coefficient * (alpha / pivot_alpha) ** exponent (analysis.cloud)."""
    for name, value in (("alpha", alpha), ("coefficient", coefficient),
                        ("pivot_alpha", pivot_alpha)):
        if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if isinstance(exponent, (bool, np.bool_)) or not np.isfinite(exponent):
        raise ValueError("exponent must be finite")
    return coefficient * (alpha / pivot_alpha) ** exponent


def interpolate_fit_parameters(
    alpha: float,
    values: Mapping[float, NDArray[np.floating]],
    *,
    circular: bool,
) -> NDArray[np.floating]:
    """Interpolate fit products in alpha, treating phase means circularly."""
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive and finite")
    if not values:
        raise ValueError("fit parameter grid cannot be empty")
    grid = np.array(sorted(values), dtype=float)
    if not np.all(np.isfinite(grid)) or np.any(grid <= 0):
        raise ValueError("fit alpha nodes must be positive and finite")
    arrays = [np.asarray(values[key], dtype=float) for key in sorted(values)]
    shape = arrays[0].shape
    if len(shape) < 2 or shape[0] != 2 or any(size == 0 for size in shape):
        raise ValueError("fit parameters must have shape (2, spatial...)")
    for array in arrays:
        if array.shape != shape or not np.all(np.isfinite(array)):
            raise ValueError("fit parameters must have matching shapes and finite values")
        if np.any(array[1] < 0):
            raise ValueError("fit standard deviations must be nonnegative")
    # Exact products must survive unchanged, including phases outside [-pi, pi].
    if alpha in values:
        return np.array(values[alpha], dtype=float, copy=True)
    if alpha <= grid[0]:
        return np.array(values[float(grid[0])], copy=True)
    if alpha >= grid[-1]:
        return np.array(values[float(grid[-1])], copy=True)
    lower_index = np.searchsorted(grid, alpha) - 1
    lower, upper = grid[lower_index : lower_index + 2]
    weight = (alpha - lower) / (upper - lower)
    first, second = arrays[lower_index], arrays[lower_index + 1]
    if not circular:
        return (1 - weight) * first + weight * second
    # The first row is a phase mean; interpolate it on the unit circle.  The
    # second row is an ordinary standard deviation and remains linear.
    result = np.empty_like(first)
    vector = (1 - weight) * np.exp(1j * first[0]) + weight * np.exp(1j * second[0])
    result[0] = np.angle(vector)
    result[1] = (1 - weight) * first[1] + weight * second[1]
    return result


def evpa_templates(
    alpha: float,
    amplitude_parameters: Mapping[float, NDArray[np.floating]],
    phase_parameters: Mapping[float, NDArray[np.floating]],
    *,
    phase_count: int,
    time_conversion: float,
    exposure_hours: float,
    day_gap_hours: float,
    pair_separation_hours: float,
    n_principal: int,
    reference_norm: float,
) -> NDArray[np.floating]:
    """Return the two-epoch signal templates.

    Signals have unit ``gagg * cloud_amplitude_scaling(alpha)`` and shape
    ``(phase, epoch_pair, spatial...)``. ``day_gap_hours`` is the within-pair
    observation separation (one day); ``pair_separation_hours`` separates the
    two pair midpoints (five days). Phases advance by omega*time on a uniform
    phase grid without a duplicated endpoint; fit phases are interpolated in
    alpha on the unit circle. Exposure is a rectangular time average. The first
    pair midpoint is the phase origin, absorbing the unobservable common
    cosine-to-sine shift.

    ``reference_norm`` is the RAPTOR axion normalisation of the fit products;
    the signal scales inversely with it.
    """
    if (isinstance(reference_norm, (bool, np.bool_))
            or not np.isfinite(reference_norm) or reference_norm <= 0):
        raise ValueError("reference_norm must be positive and finite")
    if isinstance(phase_count, bool) or not isinstance(phase_count, (int, np.integer)) or phase_count < 2:
        raise ValueError("phase_count must be an integer >= 2")
    for name, value in (("time_conversion", time_conversion),
                        ("day_gap_hours", day_gap_hours)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    for name, value in (("exposure_hours", exposure_hours),
                        ("pair_separation_hours", pair_separation_hours)):
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be nonnegative and finite")
    if set(amplitude_parameters) != set(phase_parameters):
        raise ValueError("amplitude and phase grids differ")
    amplitude = interpolate_fit_parameters(alpha, amplitude_parameters, circular=False)
    phase = interpolate_fit_parameters(alpha, phase_parameters, circular=True)
    if amplitude.shape != phase.shape:
        raise ValueError("amplitude and phase shapes differ")
    initial_phase = np.linspace(0, 2 * np.pi, phase_count, endpoint=False).reshape(
        (phase_count,) + (1,) * (amplitude.ndim - 1)
    )
    angular_frequency = omega(alpha, n_principal)
    exposure = time_conversion * angular_frequency / 2 * exposure_hours
    constant = (
        2
        * np.sin(time_conversion * angular_frequency / 2 * day_gap_hours)
        * np.sinc(exposure / np.pi)
        / reference_norm
    )
    phase_first = initial_phase + phase[0]
    phase_second = phase_first + time_conversion * angular_frequency * pair_separation_hours
    signal = constant * amplitude[0] * np.stack((np.sin(phase_first), np.sin(phase_second)), axis=1)
    if not np.all(np.isfinite(signal)):
        raise ValueError("template calculation produced nonfinite values")
    return signal
