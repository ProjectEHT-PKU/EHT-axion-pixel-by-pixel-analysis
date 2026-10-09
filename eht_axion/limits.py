"""Final limits, with a pure product API and a compatible historical loader."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import cho_factor, cho_solve, eigh, solve_triangular
from scipy.special import logsumexp

from .model import ALPHAS, TIME_CONVERSION, cloud_amplitude_scaling, evpa_templates, load_fit_parameters
from .observations import (
    ObservationData,
    diagonal_pair_covariances,
    load_1d_observations,
    load_2d_observations,
)
from .statistics import normalized_posterior_cdf, posterior_upper_limit


@dataclass(frozen=True)
class LimitOptions:
    legacy: bool = False
    prior: str = "log-uniform"
    coupling_points: int = 60
    coupling_min: float = 1 / 137
    coupling_max: float = 1_000
    grmhd_sigma_degrees: float = 1.9752611467515768
    phase_count: int | None = None  # corrected: 64, legacy: 16
    time_conversion: float = TIME_CONVERSION
    exposure_hours: float = 6.0
    day_gap_hours: float = 24.0
    pair_separation_hours: float = 120.0
    n_principal: int = 4  # Historical omega4; pending theoretical confirmation.
    # New products use the diagonal of the same C as cor. The old loader
    # explicitly selects historical marginal sigma instead for compatibility.
    std_covariance: str = "covariance"
    duplicate_phase_endpoint: bool | None = None
    circular_phase_interpolation: bool | None = None
    shared_epoch_uncertainty: bool | None = None
    combine_epochs_by_product: bool | None = None
    include_parameter_logdet: bool | None = None
    posterior_method: str | None = None


def _legacy_quadratic(residual: NDArray, covariance: NDArray) -> float:
    factor = cho_factor(covariance, lower=True, check_finite=False)
    return float(residual @ cho_solve(factor, residual, check_finite=False))


def _legacy_upper_limit(couplings: NDArray, likelihood: NDArray) -> float:
    """Reproduce the notebook's equal-log-grid discrete cumulative sum."""
    cumulative = np.cumsum(likelihood) / likelihood.size
    target = cumulative[-1] * 0.95
    upper = int(np.searchsorted(cumulative, target, side="left"))
    if upper == 0:
        return float(np.log10(couplings[0]))
    lower = upper - 1
    fraction = (target - cumulative[lower]) / (cumulative[upper] - cumulative[lower])
    return float(
        np.log10(couplings[lower])
        + fraction * (np.log10(couplings[upper]) - np.log10(couplings[lower]))
    )


def _validate_request(dimension: str, mode: str, options: LimitOptions) -> None:
    if dimension not in {"1d", "2d"}:
        raise ValueError("dimension must be '1d' or '2d'")
    if mode not in {"std", "cor", "cor_grmhd"}:
        raise ValueError("mode must be std, cor, or cor_grmhd")
    if (isinstance(options.coupling_points, bool)
            or not isinstance(options.coupling_points, (int, np.integer))
            or options.coupling_points < 2):
        raise ValueError("coupling_points must be an integer >= 2")
    if (not np.isfinite(options.coupling_min) or not np.isfinite(options.coupling_max)
            or not 0 < options.coupling_min < options.coupling_max):
        raise ValueError("coupling bounds must be finite and satisfy 0 < min < max")
    if not np.isfinite(options.grmhd_sigma_degrees) or options.grmhd_sigma_degrees < 0:
        raise ValueError("grmhd_sigma_degrees must be nonnegative and finite")
    if options.std_covariance not in {"covariance", "legacy-sigma"}:
        raise ValueError("std_covariance must be 'covariance' or 'legacy-sigma'")
    if options.prior not in {"uniform", "log-uniform"}:
        raise ValueError("prior must be 'uniform' or 'log-uniform'")
    if options.posterior_method not in {None, "legacy", "trapezoid"}:
        raise ValueError("posterior_method must be 'legacy' or 'trapezoid'")
    posterior_method = options.posterior_method or ("legacy" if options.legacy else "trapezoid")
    if posterior_method == "legacy" and options.prior != "log-uniform":
        raise ValueError("legacy posterior_method requires prior='log-uniform' (equal-log-grid CDF)")


def calculate_from_products(
    observations: ObservationData,
    amplitude: Mapping[float, NDArray],
    phase: Mapping[float, NDArray],
    alphas=ALPHAS,
    dimension: str = "1d",
    mode: str = "std",
    options: LimitOptions | None = None,
) -> tuple[NDArray, NDArray, dict]:
    """Compute log10(c95) without loading, printing, or saving anything.

    Inference always assumes zero observed EVPA day differences (null).
    Selected pair differences, in radians, are diagnostic-only: their schema
    and shape are validated for compatibility but they never enter the likelihood.
    Selected observation covariance blocks, in radians squared, are retained;
    selected_indices index the *full* fit-product grid.

    cor_grmhd uses C_epoch + diag(sigma_GRMHD**2 + (g*b*sigma_model)**2).
    Model errors are assumed spatially diagonal and independent between the
    two epochs, while their phase nuisance parameter is shared. Corrected mode
    retains the full Gaussian log determinant. For dense C, B=C+sigma_GRMHD² I
    is factored once per epoch, then L^-1 D L^-T is diagonalized once per phase.
    All couplings use the exact eigenbasis quadratic and log1p determinant;
    zero eigenvalues are allowed. Diagonal cases use only variance vectors.

    Diagnostics contain per-alpha coupling log likelihoods, density/CDF, and
    edge-bin masses (truncation indicators, not proof of convergence). Legacy
    CDFs are explicitly the historical discrete equal-log-grid sums.
    """
    options = options or LimitOptions()
    _validate_request(dimension, mode, options)
    alpha_grid = np.asarray(alphas, dtype=float)
    if (alpha_grid.ndim != 1 or alpha_grid.size == 0 or not np.all(np.isfinite(alpha_grid))
            or np.any(alpha_grid <= 0) or np.any(np.diff(alpha_grid) <= 0)):
        raise ValueError("alphas must be finite, positive and strictly increasing")
    indices = np.asarray(observations.selected_indices)
    if (indices.ndim != 1 or indices.size == 0 or indices.dtype.kind not in "iu"
            or np.any(indices < 0) or np.unique(indices).size != indices.size):
        raise ValueError("selected_indices must be nonempty unique nonnegative integers")
    size = indices.size
    if len(observations.pair_differences) != 2 or len(observations.pair_covariances) != 2:
        raise ValueError("exactly two epoch pairs are required")
    # Compatibility validation only; measured differences do not set the center.
    pair_differences = np.asarray(observations.pair_differences, dtype=float)
    if pair_differences.shape != (2, size) or not np.all(np.isfinite(pair_differences)):
        raise ValueError("pair differences must be finite with shape (2, selected size)")
    covariances = tuple(np.asarray(c, dtype=float) for c in observations.pair_covariances)
    for covariance in covariances:
        if covariance.shape != (size, size) or not np.all(np.isfinite(covariance)):
            raise ValueError("pair covariance must be finite with shape (selected size, selected size)")
        if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-14):
            raise ValueError("pair covariance must be symmetric")
    if mode == "std":
        if options.std_covariance == "legacy-sigma":
            if len(observations.standard_deviations) != 4:
                raise ValueError("legacy-sigma requires four daily standard deviations")
            for std in observations.standard_deviations:
                if (np.shape(std) != (size,) or not np.all(np.isfinite(std))
                        or np.any(np.asarray(std) < 0)):
                    raise ValueError("daily standard deviations must be finite and nonnegative")
            covariances = diagonal_pair_covariances(observations)
        variances = tuple(np.diag(c).copy() for c in covariances)
    else:
        variances = None
    diagonal = mode == "std" or all(
        np.array_equal(c, np.diag(np.diag(c))) for c in covariances
    )
    if diagonal and variances is None:
        variances = tuple(np.diag(c).copy() for c in covariances)
    grmhd_variance = np.deg2rad(options.grmhd_sigma_degrees) ** 2
    if diagonal:
        offset = grmhd_variance if mode == "cor_grmhd" else 0.0
        if any(np.any(v + offset <= 0) for v in variances):
            raise ValueError("covariance diagonal must be positive")
        factors = None
    else:
        # Validate and factor fixed C once. In the varying branch only C+sigma
        # needs to be positive definite; no inverse of the base C is required.
        factors = (
            tuple(cho_factor(c, lower=True, check_finite=False) for c in covariances)
            if mode != "cor_grmhd" else None
        )
    whiteners = []
    if mode == "cor_grmhd" and not diagonal:
        for covariance in covariances:
            background = covariance.copy()
            background[np.diag_indices(size)] += grmhd_variance
            factor, _ = cho_factor(background, lower=True, check_finite=False)
            inverse_lower = solve_triangular(factor, np.eye(size), lower=True, check_finite=False)
            whiteners.append((inverse_lower, 2 * np.log(np.diag(factor)).sum()))
    couplings = np.logspace(
        np.log10(options.coupling_min), np.log10(options.coupling_max), options.coupling_points
    )
    combine_product = (not options.legacy if options.combine_epochs_by_product is None
                       else options.combine_epochs_by_product)
    include_logdet = (not options.legacy if options.include_parameter_logdet is None
                      else options.include_parameter_logdet)
    posterior_method = options.posterior_method or ("legacy" if options.legacy else "trapezoid")
    cuts, loglikes, densities, cdfs, boundaries = [], [], [], [], []
    resolved_phase_count = None
    # Select before model validation: excluded zero-response/NaN pixels are not data.
    def selected_parameters(parameters):
        result = {}
        expected_ndim = 2 if dimension == "1d" else 3
        for key, value in parameters.items():
            value = np.asarray(value)
            if value.ndim != expected_ndim or value.shape[0] != 2:
                raise ValueError("fit spatial shape does not match dimension")
            flat = value.reshape(2, -1)
            if np.any(indices >= flat.shape[-1]):
                raise ValueError("selected_indices exceed fit spatial grid")
            selected = flat[:, indices]
            result[key] = selected if dimension == "1d" else selected[:, None, :]
        return result

    amplitude = selected_parameters(amplitude)
    phase = selected_parameters(phase)
    for alpha in alpha_grid:
        signal, signal_std = evpa_templates(
            float(alpha), amplitude, phase,
            phase_count=options.phase_count, legacy=options.legacy,
            duplicate_endpoint=options.duplicate_phase_endpoint,
            circular_phase_interpolation=options.circular_phase_interpolation,
            legacy_shared_uncertainty=options.shared_epoch_uncertainty,
            time_conversion=options.time_conversion, exposure_hours=options.exposure_hours,
            day_gap_hours=options.day_gap_hours, pair_separation_hours=options.pair_separation_hours,
            n_principal=options.n_principal,
        )
        if signal.ndim != (3 if dimension == "1d" else 4):
            raise ValueError("fit spatial shape does not match dimension")
        signal = signal.reshape(signal.shape[0], 2, -1)
        signal_std = signal_std.reshape(signal_std.shape[0], 2, -1)
        scale = cloud_amplitude_scaling(float(alpha))
        signal = signal * scale
        signal_std = signal_std * scale
        resolved_phase_count = signal.shape[0]
        if mode != "cor_grmhd":
            # Null residual r=g*t: only the g² coefficient is needed.
            # Solve all phase RHS together using the observation covariance.
            coefficients = []
            for epoch in range(2):
                template = signal[:, epoch]
                if diagonal:
                    solved_template = template / variances[epoch]
                else:
                    solved_template = cho_solve(factors[epoch], template.T, check_finite=False).T
                coefficients.append(np.sum(template * solved_template, axis=1))
            epoch_loglikes = np.stack([
                -0.5 * couplings[:, None] ** 2 * coefficient
                for coefficient in coefficients
            ], axis=-1)
            # Fixed-C log determinants are constant in g and phase, so omitted
            # in both paths, matching historical likelihood normalization.
        else:
            epoch_loglikes = np.empty((couplings.size, resolved_phase_count, 2))
            for epoch in range(2):
                if not diagonal:
                    inverse_lower, background_logdet = whiteners[epoch]
                    white_signal = signal[:, epoch] @ inverse_lower.T
                for phase_index in range(resolved_phase_count):
                    model_variance = signal_std[phase_index, epoch] ** 2
                    if diagonal:
                        variance = (variances[epoch] + grmhd_variance
                                    + couplings[:, None] ** 2 * model_variance)
                        residual = couplings[:, None] * signal[phase_index, epoch]
                        quadratic = np.sum(residual**2 / variance, axis=1)
                        logdet = np.log(variance).sum(axis=1)
                    else:
                        white_variance = (inverse_lower * model_variance) @ inverse_lower.T
                        eigenvalues, eigenvectors = eigh(white_variance, check_finite=False)
                        # D is PSD. Only roundoff-sized negative eigenvalues may
                        # be clipped, never a genuinely indefinite model error.
                        tolerance = 1e-12 * max(1.0, float(np.max(np.abs(eigenvalues))))
                        if eigenvalues[0] < -tolerance:
                            raise ValueError("whitened model covariance is not positive semidefinite")
                        eigenvalues = np.maximum(eigenvalues, 0)
                        # An invertible whitening preserves rank(D). Set the
                        # known nullspace exactly to zero, avoiding g² times
                        # roundoff-positive eigenvalues at large couplings.
                        nullity = size - np.count_nonzero(model_variance)
                        eigenvalues[:nullity] = 0
                        template = white_signal[phase_index] @ eigenvectors
                        addition = couplings[:, None] ** 2 * eigenvalues
                        residual = couplings[:, None] * template
                        quadratic = np.sum(residual**2 / (1 + addition), axis=1)
                        logdet = background_logdet + np.log1p(addition).sum(axis=1)
                    epoch_loglikes[:, phase_index, epoch] = -0.5 * (
                        quadratic + (logdet if include_logdet else 0)
                    )
        phase_loglikes = (epoch_loglikes.sum(axis=-1) if combine_product
                          else logsumexp(epoch_loglikes, axis=-1))
        log_likelihood = logsumexp(phase_loglikes, axis=1) - np.log(resolved_phase_count)
        if not np.all(np.isfinite(log_likelihood)):
            raise ValueError("likelihood calculation produced nonfinite values")
        if posterior_method == "legacy":
            likelihood = np.exp(log_likelihood - np.max(log_likelihood))
            cut = _legacy_upper_limit(couplings, likelihood)
            density = likelihood / likelihood.sum()  # probability mass, not density in g
            cdf = np.minimum(np.cumsum(density), 1.0)
            cdf[-1] = 1.0
            lower_mass, upper_mass = density[0], density[-1]
        else:
            density, cdf = normalized_posterior_cdf(couplings, log_likelihood, prior=options.prior)
            cut = np.log10(posterior_upper_limit(couplings, log_likelihood, prior=options.prior))
            lower_mass, upper_mass = cdf[1], 1 - cdf[-2]
        cuts.append(cut)
        loglikes.append(log_likelihood)
        densities.append(density)
        cdfs.append(cdf)
        boundaries.append({
            "lower_bin_mass": float(lower_mass), "upper_bin_mass": float(upper_mass),
            "peak_at_lower_bound": bool(np.argmax(log_likelihood) == 0),
            "peak_at_upper_bound": bool(np.argmax(log_likelihood) == couplings.size - 1),
            "limit_in_last_bin": bool(cut >= np.log10(couplings[-2])),
        })
    diagnostics = {
        "alphas": alpha_grid.tolist(), "couplings": couplings.tolist(),
        "log_likelihood": np.asarray(loglikes).tolist(),
        "posterior_density": np.asarray(densities).tolist(),
        "cdf": np.asarray(cdfs).tolist(), "boundaries": boundaries,
        "posterior_method": posterior_method, "prior": options.prior,
        "density_measure": "discrete-log-grid-mass" if posterior_method == "legacy" else "dg",
        "phase_count": resolved_phase_count, "legacy": options.legacy,
        "model_error_assumption": "spatially diagonal, independent between epochs",
        "fixed_covariance_logdet_omitted": mode != "cor_grmhd",
        "std_covariance": options.std_covariance,
    }
    return alpha_grid.copy(), np.asarray(cuts), diagnostics


def calculate_limits(
    root: Path,
    *,
    dimension: str,
    mode: str,
    options: LimitOptions,
    observation_dir: Path | None = None,
    fit_parameter_dir: Path | None = None,
) -> tuple[NDArray, NDArray]:
    """Compatibility loader; preserves historical std sigma, grid and return shape."""
    _validate_request(dimension, mode, options)
    observations = (
        load_1d_observations(root, observation_dir)
        if dimension == "1d" else load_2d_observations(root, observation_dir)
    )
    amplitude, phase = load_fit_parameters(root, dimension, fit_parameter_dir)
    alphas, cuts, _ = calculate_from_products(
        observations, amplitude, phase, ALPHAS, dimension, mode,
        replace(options, std_covariance="legacy-sigma"),
    )
    for alpha, cut in zip(alphas, cuts):
        print(f"{dimension} {mode} alpha={alpha:.2f}: log10(c95)={cut:.8f}", flush=True)
    return alphas, cuts
