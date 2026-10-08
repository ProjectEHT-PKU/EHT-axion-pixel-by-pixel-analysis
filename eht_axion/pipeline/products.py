"""Named, non-pickled inference products with a consistent ensemble covariance."""

from pathlib import Path
import json
import numpy as np
from ..angles import circular_difference
from ..observations import ObservationData, DAYS, DAY_PAIRS


def load_templates(path, dimension):
    with np.load(path, allow_pickle=False) as z:
        alphas = z["alphas"].copy()
        keys = ["amplitude_mean", "amplitude_std", "phase_mean", "phase_std"]
        values = [z[k].copy() for k in keys]
    if dimension not in {"1d", "2d"} or alphas.ndim != 1 or not len(alphas):
        raise ValueError("Invalid template dimension/alpha grid")
    if not np.isfinite(alphas).all() or np.any(np.diff(alphas) <= 0):
        raise ValueError("Template alpha grid must be finite and increasing")
    expected_ndim = 2 if dimension == "1d" else 3
    if any(v.shape != values[0].shape or v.ndim != expected_ndim for v in values):
        raise ValueError("Template summary shapes differ")
    if values[0].shape[0] != len(alphas):
        raise ValueError("Template alpha axis mismatch")
    amplitude = {float(a): np.stack((values[0][k], values[1][k])) for k, a in enumerate(alphas)}
    phase = {float(a): np.stack((values[2][k], values[3][k])) for k, a in enumerate(alphas)}
    valid = np.ones(values[0].shape[1:], dtype=bool)
    for v in values:
        valid &= np.isfinite(v).all(axis=0)
    valid &= (values[1] >= 0).all(axis=0) & (values[3] >= 0).all(axis=0)
    return alphas, amplitude, phase, valid


def build_observations(
    eht_dir,
    template_path,
    dimension,
    output_path,
    *,
    angular_offset=4,
    spacing=10,
    covariance_ridge=0.0,
):
    """Center and scatter from the same local unwrapped ensemble.

    Std is exactly diag(cor), so mode comparisons isolate spatial association.
    Observed pair difference is later minus earlier. Cross-day covariance is
    assumed zero, not inferred from common hyperparameter labels.
    """
    directory = Path(eht_dir)
    _, _, _, template_valid = load_templates(template_path, dimension)
    with np.load(directory / f"summary_{dimension}.npz", allow_pickle=False) as summary:
        means = summary["mean"].copy()
        stds = summary["std"].copy()
        centers = summary["center"].copy()
        unwrapped_means = summary["unwrapped_mean"].copy()
    if centers.shape != means.shape or unwrapped_means.shape != means.shape:
        raise ValueError("Missing or incompatible circular branch metadata")
    if means.shape != stds.shape or means.shape[0] != 4 or means.shape[1:] != template_valid.shape:
        raise ValueError("Observation/template spatial shapes differ")
    valid = template_valid & np.isfinite(means).all(axis=0) & np.isfinite(stds).all(axis=0)
    if dimension == "2d":
        mask = np.load(directory / "mask.npy", allow_pickle=False)
        if mask.dtype != bool or mask.shape != valid.shape:
            raise ValueError("Invalid 2D mask")
        valid &= mask
    else:
        if spacing < 1 or not 0 <= angular_offset < spacing:
            raise ValueError("Invalid angular spacing/offset")
        valid &= np.arange(valid.size) % spacing == angular_offset
    selected = np.flatnonzero(valid.ravel())
    if not len(selected):
        raise ValueError("No valid observations after template/mask selection")
    if not np.isfinite(covariance_ridge) or not 0 <= covariance_ridge < 1:
        raise ValueError("Ridge must be a finite fraction in [0,1)")
    covariances, fitted_means, biases, n = [], [], [], []
    for k, day in enumerate(DAYS):
        if dimension == "2d":
            data = np.load(directory / f"stokes_{day}.npy", mmap_mode="r", allow_pickle=False)
            if data.ndim != 4 or data.shape[1:3] != means.shape[1:] or data.shape[-1] != 4:
                raise ValueError("Stokes shape does not match summary")
            q = data[..., 1].reshape(len(data), -1)[:, selected]
            u = data[..., 2].reshape(len(data), -1)[:, selected]
            theta = 0.5 * np.arctan2(u, q)
            if np.any((q == 0) & (u == 0)):
                raise ValueError("Selected pixel has undefined zero polarization")
        else:
            data = np.load(directory / f"evpa_1d_{day}.npy", mmap_mode="r", allow_pickle=False)
            if data.ndim != 2 or data.shape[1:] != means.shape[1:]:
                raise ValueError("1D ensemble shape mismatch")
            theta = data[:, selected]
        if len(theta) < 2 or not np.isfinite(theta).all():
            raise ValueError(f"Day {day}: finite ensemble with at least two samples required")
        center = centers[k].ravel()[selected]
        unwrapped_mean = unwrapped_means[k].ravel()[selected]
        residual = center + circular_difference(theta, center) - unwrapped_mean
        bias = residual.mean(axis=0)
        if not np.allclose(bias, 0, atol=1e-12, rtol=0):
            raise ValueError("Ensemble does not match saved circular branch/mean")
        covariance = residual.T @ residual / len(residual)
        diagonal = np.diag_indices_from(covariance)
        if (covariance[diagonal] <= 0).any():
            raise ValueError("Zero ensemble variance; cannot infer with singular noise")
        covariance[diagonal] += covariance_ridge * np.median(covariance[diagonal])
        covariances.append(covariance)
        fitted_means.append(means[k].ravel()[selected])
        biases.append(bias)
        n.append(len(theta))
    pairs = [circular_difference(fitted_means[j], fitted_means[i]) for i, j in DAY_PAIRS]
    pair_cov = [covariances[i] + covariances[j] for i, j in DAY_PAIRS]
    # Fail here with a readable message rather than deep in Cholesky inference.
    for covariance in pair_cov:
        np.linalg.cholesky(covariance)
    np.savez(
        output_path,
        means=np.array(fitted_means),
        standard_deviations=np.sqrt(np.diagonal(covariances, axis1=1, axis2=2)),
        covariances=np.array(covariances),
        pair_differences=np.array(pairs),
        pair_covariances=np.array(pair_cov),
        selected_indices=selected,
        dimension=dimension,
        spatial_shape=np.array(valid.shape),
    )
    report = {
        "dimension": dimension,
        "selected_points": len(selected),
        "sample_counts": n,
        "covariance_definition": "centered local ensemble population moment (ddof=0)",
        "std_definition": "diagonal of the same covariance used by cor",
        "cross_day_observation_covariance": "zero (assumed)",
        "ridge_fraction_of_median_variance": covariance_ridge,
        "max_recenter_degrees": float(np.rad2deg(np.max(np.abs(biases)))),
        "pair_covariance_condition_numbers": [float(np.linalg.cond(c)) for c in pair_cov],
    }
    Path(output_path).with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def load_observations(path):
    with np.load(path, allow_pickle=False) as z:
        arrays = {
            k: z[k].copy()
            for k in (
                "means",
                "standard_deviations",
                "covariances",
                "pair_differences",
                "pair_covariances",
                "selected_indices",
            )
        }
    if not len(arrays["selected_indices"]) or any(not np.isfinite(v).all() for v in arrays.values()):
        raise ValueError("Empty/non-finite observation product")
    return ObservationData(**{k: tuple(v) if k != "selected_indices" else v for k, v in arrays.items()})
