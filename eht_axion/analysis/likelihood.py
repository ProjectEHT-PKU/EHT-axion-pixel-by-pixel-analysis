"""Gaussian likelihood of the observed EHT day-pair EVPA differences and the limit on c.

LimitOptions holds every inference value (from experiment.yaml via
options_from_config of scripts/analysis/limits.py). calculate_observed multiplies the two
day-pair Gaussian likelihoods at a shared axion phase, averages over a uniform
phase grid and returns the posterior upper limit on the coupling c per alpha.

Scope: the covariances are the scatter of the D-term Monte Carlo images, not
calibrated measurement errors; response-model variance, cross-pair covariance and
covariance-estimation uncertainty are not included. The data are principal-branch
EVPA differences in radians and the residuals d - g*t are not wrapped again (a
Gaussian on that branch, not a wrapped-normal likelihood). Nothing is written.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from numbers import Real
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import solve_triangular
from scipy.special import logsumexp

from .model import cloud_amplitude_scaling, evpa_templates
from .observations import ObservationData
from .posterior import normalized_posterior_cdf, posterior_upper_limit

# Covariance modes: diagonal (std), full D-scatter (cor), full + fixed GRMHD background.
MODES = ("std", "cor", "cor_fixed_background")


@dataclass(frozen=True, kw_only=True)
class LimitOptions:
    """Inference options. Every value is supplied by the caller
    (options_from_config of scripts/analysis/limits.py reads the experiment file); there
    are no defaults and no alternative methods. Field order is the
    recorded JSON order.
    """

    prior: str
    credibility: float
    coupling_points: int
    coupling_min: float
    coupling_max: float
    grmhd_sigma_degrees: float
    phase_count: int
    time_conversion: float
    exposure_hours: float
    day_gap_hours: float
    pair_separation_hours: float
    n_principal: int
    reference_norm: float
    cloud_coefficient: float
    cloud_pivot_alpha: float
    cloud_exponent: float


def validate_options(dimension: str, options: LimitOptions) -> None:
    if dimension not in {"1d", "2d"}:
        raise ValueError("dimension must be '1d' or '2d'")
    if (isinstance(options.coupling_points, bool)
            or not isinstance(options.coupling_points, (int, np.integer))
            or options.coupling_points < 2):
        raise ValueError("coupling_points must be an integer >= 2")
    if (not np.isfinite(options.coupling_min) or not np.isfinite(options.coupling_max)
            or not 0 < options.coupling_min < options.coupling_max):
        raise ValueError("coupling bounds must be finite and satisfy 0 < min < max")
    if not np.isfinite(options.grmhd_sigma_degrees) or options.grmhd_sigma_degrees < 0:
        raise ValueError("grmhd_sigma_degrees must be nonnegative and finite")
    if options.prior not in {"uniform", "log-uniform"}:
        raise ValueError("prior must be 'uniform' or 'log-uniform'")
    if not 0 < options.credibility < 1:
        raise ValueError("credibility must lie strictly between zero and one")


def _validate_options(dimension: str, mode: str, options: LimitOptions) -> None:
    if not isinstance(options, LimitOptions):
        raise TypeError("options must be a LimitOptions instance")
    if mode not in {"std", "cor", "cor_fixed_background"}:
        raise ValueError("mode must be std, cor, or cor_fixed_background")
    for name in ("coupling_min", "coupling_max", "grmhd_sigma_degrees", "time_conversion",
                 "exposure_hours", "day_gap_hours", "pair_separation_hours", "reference_norm",
                 "cloud_coefficient", "cloud_pivot_alpha", "cloud_exponent", "credibility"):
        value = getattr(options, name)
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
            raise ValueError(f"{name} must be a finite real number")
    validate_options(dimension, options)


def _selected_parameters(parameters: Mapping, indices: NDArray, dimension: str) -> dict:
    if not isinstance(parameters, Mapping) or not parameters:
        raise ValueError("fit parameter grid must be a nonempty mapping")
    result = {}
    for alpha, value in parameters.items():
        if (isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, Real)
                or not np.isfinite(alpha) or alpha <= 0):
            raise ValueError("fit alpha nodes must be positive and finite")
        value = np.asarray(value, dtype=float)
        if (value.ndim != (2 if dimension == "1d" else 3) or value.shape[0] != 2
                or any(size == 0 for size in value.shape)):
            raise ValueError("fit spatial shape does not match dimension")
        flat = value.reshape(2, -1)
        if np.any(indices >= flat.shape[-1]):
            raise ValueError("selected_indices exceed fit spatial grid")
        # Indices address the full flattened template grid.
        # Unselected NaN/zero-response pixels are not part of this likelihood.
        selected = flat[:, indices]
        result[alpha] = selected if dimension == "1d" else selected[:, None, :]
    return result


def _factor_covariance(covariance: NDArray, size: int) -> NDArray:
    if covariance.shape != (size, size) or not np.all(np.isfinite(covariance)):
        raise ValueError("pair covariance must be finite with shape (selected size, selected size)")
    if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-14):
        raise ValueError("pair covariance must be symmetric")
    try:
        return np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("pair covariance must be strictly positive definite") from exc


def calculate_observed(
    observations: ObservationData,
    amplitude: Mapping[float, NDArray],
    phase: Mapping[float, NDArray],
    alphas,
    dimension: str = "1d",
    mode: str = "std",
    *,
    options: LimitOptions,
) -> tuple[NDArray, NDArray, dict]:
    """Return (alpha, log10(g95), diagnostics) using actual two-pair data.

    ``observations.pair_differences`` must have shape (2, N), in radians;
    covariance blocks are already selected, in radians squared. Fit products
    retain their full spatial grids; selected_indices index the flattened grid.

    ``std`` uses diag(C_pair); ``cor`` uses C_pair; ``cor_fixed_background`` uses C_pair + sigma_bg**2 I,
    with sigma_bg = radians(options.grmhd_sigma_degrees), fixed for all g/phi.
    All supplied observation blocks must be finite SPD, even in std mode.
    There is NO g**2 response variance in any mode.

    For each g, multiply the two independent pair likelihoods at the SAME phi,
    then average over a uniform, nonduplicated phase grid. Fixed-C Gaussian
    logdet and normalization constants are omitted because they depend on
    neither g nor phi. The posterior uses an explicit uniform or log-uniform
    prior truncated to [coupling_min, coupling_max], integrated by trapezoids
    in g. The returned 95% bound is conditional on this finite grid/prior, not
    a frequentist coverage claim. Edge diagnostics do not certify convergence.
    See the module docstring for the scope of the covariance, the missing model
    variance and the EVPA branch.
    """
    _validate_options(dimension, mode, options)
    alpha_grid = np.asarray(alphas, dtype=float)
    if (alpha_grid.ndim != 1 or alpha_grid.size == 0 or not np.all(np.isfinite(alpha_grid))
            or np.any(alpha_grid <= 0) or np.any(np.diff(alpha_grid) <= 0)):
        raise ValueError("alphas must be finite, positive and strictly increasing")
    if not isinstance(observations, ObservationData):
        raise TypeError("observations must be ObservationData")
    indices = np.asarray(observations.selected_indices)
    if (indices.ndim != 1 or indices.size == 0 or indices.dtype.kind not in "iu"
            or np.any(indices < 0) or np.unique(indices).size != indices.size):
        raise ValueError("selected_indices must be nonempty unique nonnegative integers")
    size = indices.size
    data = np.asarray(observations.pair_differences, dtype=float)
    if data.shape != (2, size) or not np.all(np.isfinite(data)):
        raise ValueError("pair differences must be finite with shape (2, selected size)")
    if len(observations.pair_covariances) != 2:
        raise ValueError("exactly two epoch pair covariances are required")
    covariances = tuple(np.asarray(c, dtype=float).copy() for c in observations.pair_covariances)
    factors = tuple(_factor_covariance(c, size) for c in covariances)
    if mode == "std":
        covariances = tuple(np.diag(np.diag(c)) for c in covariances)
    background_variance = (float(np.deg2rad(options.grmhd_sigma_degrees) ** 2)
                           if mode == "cor_fixed_background" else 0.0)
    if background_variance:
        for covariance in covariances:
            covariance[np.diag_indices(size)] += background_variance
    if mode == "std" or background_variance:
        factors = tuple(_factor_covariance(c, size) for c in covariances)
    white_data = tuple(solve_triangular(f, d, lower=True, check_finite=False)
                       for f, d in zip(factors, data))
    data_quadratic = sum(float(d @ d) for d in white_data)
    amplitude = _selected_parameters(amplitude, indices, dimension)
    phase = _selected_parameters(phase, indices, dimension)
    couplings = np.logspace(np.log10(options.coupling_min), np.log10(options.coupling_max),
                            options.coupling_points)
    if not np.all(np.isfinite(couplings)) or np.any(np.diff(couplings) <= 0):
        raise ValueError("coupling grid must be finite and strictly increasing")
    cuts, loglikes, phase_loglikes, densities, cdfs, boundaries, fits = [], [], [], [], [], [], []
    for alpha in alpha_grid:
        signal = evpa_templates(
            float(alpha), amplitude, phase, phase_count=options.phase_count,
            time_conversion=options.time_conversion, exposure_hours=options.exposure_hours,
            day_gap_hours=options.day_gap_hours, pair_separation_hours=options.pair_separation_hours,
            n_principal=options.n_principal, reference_norm=options.reference_norm,
        )
        signal = signal.reshape(signal.shape[0], 2, size) * cloud_amplitude_scaling(
            float(alpha), coefficient=options.cloud_coefficient,
            pivot_alpha=options.cloud_pivot_alpha, exponent=options.cloud_exponent,
        )
        phase_count = signal.shape[0]
        cross, power = np.zeros(phase_count), np.zeros(phase_count)
        for epoch in range(2):
            white_template = solve_triangular(factors[epoch], signal[:, epoch].T,
                                             lower=True, check_finite=False).T
            cross += white_template @ white_data[epoch]
            power += np.einsum("pi,pi->p", white_template, white_template)
        # Only (g, phi) arrays, never a (g, phi, pair, pixel) residual cube.
        with np.errstate(over="ignore", invalid="ignore"):
            quadratic = (data_quadratic - 2 * couplings[:, None] * cross
                         + couplings[:, None] ** 2 * power)
        if not np.all(np.isfinite(quadratic)):
            raise ValueError("likelihood calculation produced nonfinite values")
        joint = -0.5 * np.maximum(quadratic, 0.0)  # protect exact-fit roundoff
        log_likelihood = logsumexp(joint, axis=1) - np.log(phase_count)
        density, cdf = normalized_posterior_cdf(couplings, log_likelihood, prior=options.prior)
        cut = np.log10(posterior_upper_limit(couplings, log_likelihood, credibility=options.credibility,
                                             prior=options.prior))
        g_index, p_index = np.unravel_index(np.argmax(joint), joint.shape)
        prediction = couplings[g_index] * signal[p_index]
        fits.append({
            "joint_grid_coupling": float(couplings[g_index]), "joint_grid_phase_index": int(p_index),
            "joint_grid_phase_radians": float(2 * np.pi * p_index / phase_count),
            "joint_grid_chi2": float(-2 * joint[g_index, p_index]),
            "prediction_radians": prediction.tolist(), "residual_radians": (data - prediction).tolist(),
            "marginal_mle_coupling": float(couplings[np.argmax(log_likelihood)]),
            "data_quadratic": data_quadratic, "template_power_by_phase": power.tolist(),
            "data_template_cross_by_phase": cross.tolist(),
        })
        cuts.append(cut)
        loglikes.append(log_likelihood)
        phase_loglikes.append(joint)
        densities.append(density)
        cdfs.append(cdf)
        boundaries.append({
            "lower_bin_mass": float(cdf[1]), "upper_bin_mass": float(1 - cdf[-2]),
            "peak_at_lower_bound": bool(np.argmax(log_likelihood) == 0),
            "peak_at_upper_bound": bool(np.argmax(log_likelihood) == couplings.size - 1),
            "limit_in_last_bin": bool(cut >= np.log10(couplings[-2])),
        })
    effective = {key: value.item() if isinstance(value, np.generic) else value
                 for key, value in asdict(options).items()}
    diagnostics = {
        "mode": mode, "dimension": dimension,
        "data_nonzero_truth": bool(np.any(data != 0)), "data_radians": data.tolist(),
        "selected_indices": indices.tolist(), "effective_options": effective,
        "effective_pair_covariances": [c.tolist() for c in covariances],
        "background_variance_radians2": background_variance,
        "alphas": alpha_grid.tolist(), "couplings": couplings.tolist(),
        "log_likelihood": np.asarray(loglikes).tolist(),
        "phase_log_likelihood": np.asarray(phase_loglikes).tolist(),
        "phase_count": phase_count, "phase_prior": "uniform discrete shared phase; endpoint excluded",
        "posterior_density": np.asarray(densities).tolist(), "cdf": np.asarray(cdfs).tolist(),
        "boundaries": boundaries, "fits": fits, "prior": options.prior,
        "coupling_prior_bounds": [float(couplings[0]), float(couplings[-1])],
        "credibility": options.credibility, "posterior_method": "trapezoid", "density_measure": "dg",
        "fixed_covariance_logdet_omitted": True, "response_model_variance_included": False,
        "warnings": [
            "Working image-ensemble covariance is not calibrated observational sampling covariance.",
            "Excludes response model variance, cross-pair covariance and covariance-estimation uncertainty.",
            "Gaussian approximation on supplied principal EVPA branch; d-g*t residual is not rewrapped.",
            "Credible limit conditional on the finite grid and prior; "
            "not a median expected limit or a coverage statement.",
        ],
    }
    return alpha_grid.copy(), np.asarray(cuts), diagnostics
