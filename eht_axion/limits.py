"""Likelihood/prior options of the limit calculation (research.observed)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, kw_only=True)
class LimitOptions:
    """Inference options. Every value is supplied by the caller
    (scripts/analysis_limits.options_from_config reads experiment.yaml); there
    are no defaults and no alternative (historical) methods. Field order is the
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
