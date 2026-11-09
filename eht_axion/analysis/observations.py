"""Selected EHT features -> observations: daily circular moments and the two day-pair differences.

The D-term Monte Carlo draws give, per day and feature, a circular mean EVPA and
the population scatter of the draws (not the error of the mean). Days are
treated as independent: the pair covariance is the sum of the two daily
covariances and equal draw indices have no paired meaning. The zero-data control
has the same features and covariance and day differences set to zero.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from ..evpa import circular_difference

DAYS = (95, 96, 100, 101)


@dataclass(frozen=True)
class ObservationData:
    means: tuple[NDArray[np.floating], ...]
    standard_deviations: tuple[NDArray[np.floating], ...]
    covariances: tuple[NDArray[np.floating], ...]
    pair_differences: tuple[NDArray[np.floating], ...]
    pair_covariances: tuple[NDArray[np.floating], ...]
    selected_indices: NDArray[np.integer]


def moments(values):
    """EVPA circular mean, centered local population covariance and resultant.

    Covariance is scatter of individual draws (ddof=0), NOT covariance of the
    mean: there is deliberately no extra division by the Monte Carlo count.
    """
    values = np.asarray(values)
    if (values.dtype.kind not in "fiu" or values.ndim != 2
            or values.shape[0] < 2 or values.shape[1] < 1
            or not np.isfinite(values).all()):
        raise ValueError("EVPA samples must be finite real (draw>=2, feature>=1)")
    values = values.astype(float)
    phasor = np.mean(np.exp(2j * values), axis=0)
    center = .5 * np.angle(phasor)
    residual = circular_difference(values, center)
    shift = residual.mean(axis=0)
    residual -= shift
    return center, residual.T @ residual / len(values), np.abs(phasor)


def covariance_diagnostics(covariances):
    """Unregularized numerical diagnostics for the effective likelihood blocks."""
    return [dict(eigenvalue_min=float(np.linalg.eigvalsh(c)[0]),
                 eigenvalue_max=float(np.linalg.eigvalsh(c)[-1]),
                 condition_number=float(np.linalg.cond(c))) for c in covariances]


def branch_risk(pairs, covariances):
    """Report branch proximity without data-dependent exclusion or tuning."""
    report = []
    for difference, covariance in zip(pairs, covariances):
        sigma = np.sqrt(np.diag(covariance))
        distance = np.pi / 2 - np.abs(difference)
        relative = distance / sigma
        report.append(dict(
            max_abs_difference_radians=float(np.max(np.abs(difference))),
            max_abs_difference_degrees=float(np.rad2deg(np.max(np.abs(difference)))),
            marginal_sigma_radians=dict(min=float(sigma.min()), median=float(np.median(sigma)),
                                        max=float(sigma.max())),
            marginal_sigma_degrees=dict(min=float(np.rad2deg(sigma.min())),
                median=float(np.rad2deg(np.median(sigma))), max=float(np.rad2deg(sigma.max()))),
            distance_to_seam_radians=distance.tolist(), distance_to_seam_in_sigma=relative.tolist(),
            count_within_one_sigma=int(np.count_nonzero(relative <= 1)),
            count_within_three_sigma=int(np.count_nonzero(relative <= 3)),
            warning='Conditional local-branch Gaussian, not wrapped normal; seam proximity is diagnostic '
                    'only and does not trigger feature exclusions.'))
    return report


def build_observations(samples, selected_indices, *, resultant_min):
    """Build the observed later-minus-earlier pair differences and the zero-data control.

    Selection depends only on concentration, not the size of day differences.
    No covariance ridge, eigenvalue clipping, or data-dependent tuning is used.
    Pair C=C_earlier+C_later; the two pairs have zero cross-block covariance by
    the declared independent-day working model. No random-index subtraction.
    """
    if (isinstance(resultant_min, (bool, np.bool_))
            or not np.isscalar(resultant_min) or not np.isreal(resultant_min)
            or not np.isfinite(resultant_min) or not 0 < resultant_min <= 1):
        raise ValueError("resultant_min must be finite in (0,1]")
    indices = np.asarray(selected_indices)
    if (indices.ndim != 1 or not indices.size or indices.dtype.kind not in "iu"
            or np.any(indices < 0) or len(np.unique(indices)) != len(indices)):
        raise ValueError("selected_indices must be unique nonnegative integers")
    if len(samples) != 4:
        raise ValueError("Exactly four day ensembles in order 95,96,100,101 required")
    daily = [moments(v) for v in samples]
    if any(len(v[0]) != len(indices) for v in daily):
        raise ValueError("Daily feature count differs from selected_indices")
    good = np.logical_and.reduce([v[2] >= resultant_min for v in daily])
    if not good.any():
        raise ValueError("No concentrated features remain")
    means = tuple(v[0][good] for v in daily)
    covs = tuple(v[1][np.ix_(good, good)] for v in daily)
    pairs = tuple(circular_difference(means[b], means[a]) for a, b in ((0, 1), (2, 3)))
    pair_covs = (covs[0] + covs[1], covs[2] + covs[3])
    # Fail rather than silently repair singular/indefinite working covariances.
    for covariance in pair_covs:
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as exc:
            raise ValueError("Pair working covariance is not strictly positive definite") from exc
    observed = ObservationData(means, tuple(np.sqrt(np.diag(c)) for c in covs), covs,
                               pairs, pair_covs, indices[good])
    zero_data = replace(observed, pair_differences=tuple(np.zeros_like(v) for v in pairs))
    metadata = {
        "days": list(DAYS), "sample_counts": [len(v) for v in samples],
        "resultant": np.stack([v[2][good] for v in daily]).tolist(),
        "retained_initial_positions": np.flatnonzero(good).tolist(),
        "mean_definition": "pi-periodic circular mean; not arithmetic local-branch mean",
        "local_branch_mean_minus_circular_radians": [
            circular_difference(np.asarray(v)[:, good], mean).mean(axis=0).tolist()
            for v, mean in zip(samples, means)],
        "covariance_definition": "centered local-branch population D-draw scatter (ddof=0)",
        "covariance_of_mean": False, "cross_day_covariance": "zero, assumed conditional independence",
        "pair_covariance": "C95+C96 and C100+C101; block diagonal across pairs",
        "pair_eigenvalue_min": [float(np.linalg.eigvalsh(c)[0]) for c in pair_covs],
        "pair_condition_numbers": [float(np.linalg.cond(c)) for c in pair_covs],
        "branch_risk": branch_risk(pairs, pair_covs),
        "regularization": "none", "resultant_min": float(resultant_min),
    }
    return observed, zero_data, metadata
