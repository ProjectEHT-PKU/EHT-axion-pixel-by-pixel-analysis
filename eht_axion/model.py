"""Axion EVPA templates reconstructed from saved GRMHD fit products."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
from numpy.typing import NDArray

ALPHAS = np.array((0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44))
# Historical fit-product normalization; active pipelines pass producer metadata.
AXION_NORM = 0.1
TIME_CONVERSION = 0.112411


def omega(alpha: float, n_principal: int = 2) -> float:
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
    coefficient: float = 0.03,
    pivot_alpha: float = 0.05,
    exponent: float = 7 / 4,
) -> float:
    """Cloud amplitude model with historical defaults for API compatibility.

    Active inference supplies these model choices explicitly from configuration;
    they are not derived from the black-hole mass or distance.
    """
    for name, value in (("alpha", alpha), ("coefficient", coefficient),
                        ("pivot_alpha", pivot_alpha)):
        if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if isinstance(exponent, (bool, np.bool_)) or not np.isfinite(exponent):
        raise ValueError("exponent must be finite")
    return coefficient * (alpha / pivot_alpha) ** exponent


def _load_dictionary(path: Path) -> Mapping[float, NDArray[np.floating]]:
    result = np.load(path, allow_pickle=True).item()
    if not isinstance(result, dict):
        raise ValueError(f"{path} does not contain a dictionary")
    return result


def load_fit_parameters(root: Path, dimension: str, data_dir: Path | None = None):
    if dimension not in {"1d", "2d"}:
        raise ValueError("dimension must be '1d' or '2d'")
    suffix = "_1d" if dimension == "1d" else ""
    data_dir = data_dir or root / "GRMHD" / "data_out"
    amplitude = _load_dictionary(data_dir / f"A_params{suffix}_dict.npy")
    phase = _load_dictionary(data_dir / f"phi0_params{suffix}_dict.npy")
    if set(amplitude) != set(phase):
        raise ValueError("amplitude and phase grids differ")
    return amplitude, phase


def interpolate_fit_parameters(
    alpha: float,
    values: Mapping[float, NDArray[np.floating]],
    *,
    circular: bool = False,
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
    phase_count: int | None = None,
    legacy: bool = False,
    duplicate_endpoint: bool | None = None,
    circular_phase_interpolation: bool | None = None,
    legacy_shared_uncertainty: bool | None = None,
    time_conversion: float = TIME_CONVERSION,
    exposure_hours: float = 6.0,
    day_gap_hours: float = 24.0,
    pair_separation_hours: float = 120.0,
    n_principal: int = 2,
    reference_norm: float = AXION_NORM,
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    """Return two-epoch signal templates and their propagated uncertainties.

    Signals and propagated standard deviations have unit
    ``gagg * cloud_amplitude_scaling(alpha)`` and shape
    ``(phase, epoch_pair, spatial...)``. ``day_gap_hours`` is the within-pair
    observation separation (one day); ``pair_separation_hours`` separates the
    two pair midpoints (five days).
    Corrected phases advance by omega*time; legacy retains the extra /2.
    Exposure is a rectangular time average. The first pair midpoint is the
    phase origin, absorbing the unobservable common cosine-to-sine shift.

    ``reference_norm`` is the positive fit-product amplitude normalization;
    both signal and uncertainty scale inversely with it. Its historical default
    is retained for API compatibility, not fixed for active inference, which
    passes the producer metadata explicitly.
    """
    if (isinstance(reference_norm, (bool, np.bool_))
            or not np.isfinite(reference_norm) or reference_norm <= 0):
        raise ValueError("reference_norm must be positive and finite")
    phase_count = (16 if legacy else 64) if phase_count is None else phase_count
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
    duplicate_endpoint = legacy if duplicate_endpoint is None else duplicate_endpoint
    circular_phase_interpolation = (
        not legacy if circular_phase_interpolation is None else circular_phase_interpolation
    )
    legacy_shared_uncertainty = (
        legacy if legacy_shared_uncertainty is None else legacy_shared_uncertainty
    )
    if set(amplitude_parameters) != set(phase_parameters):
        raise ValueError("amplitude and phase grids differ")
    amplitude = interpolate_fit_parameters(alpha, amplitude_parameters)
    phase = interpolate_fit_parameters(
        alpha, phase_parameters, circular=circular_phase_interpolation
    )
    if amplitude.shape != phase.shape:
        raise ValueError("amplitude and phase shapes differ")
    initial_phase = np.linspace(0, 2 * np.pi, phase_count, endpoint=duplicate_endpoint).reshape(
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
    phase_second = (
        phase_first + time_conversion * angular_frequency / (2 if legacy else 1)
        * pair_separation_hours
    )
    signal = constant * amplitude[0] * np.stack((np.sin(phase_first), np.sin(phase_second)), axis=1)
    uncertainty = np.stack(
        tuple(
            np.sqrt(
                (constant * np.sin(epoch_phase) * amplitude[1]) ** 2
                + (constant * amplitude[0] * np.cos(epoch_phase) * phase[1]) ** 2
            )
            for epoch_phase in (phase_first, phase_second)
        ),
        axis=1,
    )
    if legacy_shared_uncertainty:
        uncertainty[:, 1] = uncertainty[:, 0]
    if not np.all(np.isfinite(signal)) or not np.all(np.isfinite(uncertainty)):
        raise ValueError("template calculation produced nonfinite values")
    return signal, uncertainty
