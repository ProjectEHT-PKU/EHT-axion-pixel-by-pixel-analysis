"""Likelihood/prior options shared by the limit calculation (research.observed)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class LimitOptions:
    """Inference options. Every scientific value is supplied by the caller
    (scripts/analysis_limits.options_from_config reads experiment.yaml); only
    method switches have defaults. Field order is the recorded JSON order.
    """

    legacy: bool = False
    prior: str
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
    std_covariance: str = "covariance"
    duplicate_phase_endpoint: bool | None = None
    circular_phase_interpolation: bool | None = None
    shared_epoch_uncertainty: bool | None = None
    combine_epochs_by_product: bool | None = None
    include_parameter_logdet: bool | None = None
    posterior_method: str | None = None
    reference_norm: float
    cloud_coefficient: float
    cloud_pivot_alpha: float
    cloud_exponent: float


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
