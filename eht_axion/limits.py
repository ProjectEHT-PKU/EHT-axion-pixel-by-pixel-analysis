"""Final 1-D and 2-D EHT axion upper-limit calculations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import cho_factor, cho_solve
from scipy.special import logsumexp

from .model import cloud_amplitude_scaling, evpa_templates, load_fit_parameters
from .observations import (
    diagonal_pair_covariances,
    load_1d_observations,
    load_2d_observations,
)
from .statistics import gaussian_log_likelihood, posterior_upper_limit

ALPHAS = np.array((0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44))


@dataclass(frozen=True)
class LimitOptions:
    legacy: bool = False
    observed: bool = False
    prior: str = "log-uniform"
    coupling_points: int = 60
    coupling_min: float = 1 / 137
    coupling_max: float = 1_000
    grmhd_sigma_degrees: float = 1.9752611467515768
    # Optional switches support controlled ablation studies.  None selects the
    # corrected behavior unless ``legacy`` is enabled.
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


def calculate_limits(
    root: Path,
    *,
    dimension: str,
    mode: str,
    options: LimitOptions,
    observation_dir: Path | None = None,
    fit_parameter_dir: Path | None = None,
) -> tuple[NDArray, NDArray]:
    if dimension not in {"1d", "2d"}:
        raise ValueError("dimension must be '1d' or '2d'")
    if mode not in {"std", "cor", "cor_grmhd"}:
        raise ValueError("mode must be std, cor, or cor_grmhd")
    observations = (
        load_1d_observations(root, observation_dir)
        if dimension == "1d"
        else load_2d_observations(root, observation_dir)
    )
    amplitude, phase = load_fit_parameters(root, dimension, fit_parameter_dir)
    couplings = np.logspace(
        np.log10(options.coupling_min),
        np.log10(options.coupling_max),
        options.coupling_points,
    )
    base_covariances = (
        diagonal_pair_covariances(observations) if mode == "std" else observations.pair_covariances
    )
    grmhd_variance = np.deg2rad(options.grmhd_sigma_degrees) ** 2
    # std/cor covariance is independent of alpha, coupling, and axion phase.
    # Factor it once instead of repeating thousands of O(N^3) decompositions.
    constant_factors = (
        tuple(
            cho_factor(covariance, lower=True, check_finite=False)
            for covariance in base_covariances
        )
        if mode != "cor_grmhd"
        else None
    )
    cuts = []

    for alpha in ALPHAS:
        duplicate_endpoint = (
            options.legacy
            if options.duplicate_phase_endpoint is None
            else options.duplicate_phase_endpoint
        )
        circular_interpolation = (
            not options.legacy
            if options.circular_phase_interpolation is None
            else options.circular_phase_interpolation
        )
        shared_uncertainty = (
            options.legacy
            if options.shared_epoch_uncertainty is None
            else options.shared_epoch_uncertainty
        )
        combine_by_product = (
            not options.legacy
            if options.combine_epochs_by_product is None
            else options.combine_epochs_by_product
        )
        include_logdet = (
            not options.legacy
            if options.include_parameter_logdet is None
            else options.include_parameter_logdet
        )
        posterior_method = options.posterior_method or ("legacy" if options.legacy else "trapezoid")
        signal, signal_std = evpa_templates(
            float(alpha),
            amplitude,
            phase,
            duplicate_endpoint=duplicate_endpoint,
            circular_phase_interpolation=circular_interpolation,
            legacy_shared_uncertainty=shared_uncertainty,
        )
        if dimension == "2d":
            signal = signal.reshape(signal.shape[0], 2, -1)[:, :, observations.selected_indices]
            signal_std = signal_std.reshape(signal_std.shape[0], 2, -1)[
                :, :, observations.selected_indices
            ]
        else:
            signal = signal[:, :, observations.selected_indices]
            signal_std = signal_std[:, :, observations.selected_indices]

        scale = cloud_amplitude_scaling(float(alpha))
        log_likelihood = np.empty(couplings.size)
        for coupling_index, coupling in enumerate(couplings):
            phase_log_likelihood = np.empty(signal.shape[0])
            for phase_index in range(signal.shape[0]):
                epoch_log_likelihoods = []
                for pair_index in range(2):
                    covariance = np.array(base_covariances[pair_index], copy=True)
                    if mode == "cor_grmhd":
                        diagonal = np.diag_indices_from(covariance)
                        covariance[diagonal] += grmhd_variance
                        covariance[diagonal] += (
                            signal_std[phase_index, pair_index] * coupling * scale
                        ) ** 2
                    prediction = signal[phase_index, pair_index] * coupling * scale
                    residual = prediction - (
                        observations.pair_differences[pair_index] if options.observed else 0
                    )
                    if constant_factors is not None:
                        solved = cho_solve(
                            constant_factors[pair_index], residual, check_finite=False
                        )
                        epoch_log_likelihoods.append(-0.5 * float(residual @ solved))
                    elif not include_logdet:
                        epoch_log_likelihoods.append(-0.5 * _legacy_quadratic(residual, covariance))
                    else:
                        epoch_log_likelihoods.append(
                            gaussian_log_likelihood(
                                residual,
                                covariance,
                                include_logdet=True,
                            )
                        )
                if combine_by_product:
                    # Independent observations have a product likelihood.
                    phase_log_likelihood[phase_index] = sum(epoch_log_likelihoods)
                else:
                    # Historical code adds the two likelihood values.
                    phase_log_likelihood[phase_index] = logsumexp(epoch_log_likelihoods)
            log_likelihood[coupling_index] = logsumexp(phase_log_likelihood) - np.log(
                signal.shape[0]
            )

        if posterior_method == "legacy":
            cut = _legacy_upper_limit(couplings, np.exp(log_likelihood - np.max(log_likelihood)))
        elif posterior_method == "trapezoid":
            cut = np.log10(
                posterior_upper_limit(
                    couplings,
                    log_likelihood,
                    prior=options.prior,
                )
            )
        else:
            raise ValueError("posterior_method must be 'legacy' or 'trapezoid'")
        cuts.append(cut)
        print(
            f"{dimension} {mode} alpha={alpha:.2f}: log10(c95)={cut:.8f}",
            flush=True,
        )
    return ALPHAS.copy(), np.asarray(cuts)
