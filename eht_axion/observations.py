"""Load EHT topset products and construct EVPA covariance matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .angles import circular_difference

DAYS = (95, 96, 100, 101)
DAY_PAIRS = ((0, 1), (2, 3))


@dataclass(frozen=True)
class ObservationData:
    means: tuple[NDArray[np.floating], ...]
    standard_deviations: tuple[NDArray[np.floating], ...]
    covariances: tuple[NDArray[np.floating], ...]
    pair_differences: tuple[NDArray[np.floating], ...]
    pair_covariances: tuple[NDArray[np.floating], ...]
    selected_indices: NDArray[np.integer]


def _load_summaries(data_dir: Path, suffix: str = ""):
    means, standard_deviations = [], []
    for day in DAYS:
        with np.load(data_dir / f"{day}{suffix}.npz") as data:
            means.append(np.asarray(data["arr_0"]))
            standard_deviations.append(np.asarray(data["arr_1"]))
    return tuple(means), tuple(standard_deviations)


def load_2d_observations(root: Path, data_dir: Path | None = None) -> ObservationData:
    """Construct only the selected covariance blocks.

    The historical notebook materialized four 4096x4096 matrices (~512 MiB),
    only to subselect 770 pixels.  Selecting first cuts memory and matrix-work
    by roughly 28x while producing the same covariance entries.
    """
    data_dir = data_dir or root / "EHT-Imaging" / "data_out"
    means, standard_deviations = _load_summaries(data_dir)
    condition = np.load(data_dir / "condition.npy").astype(bool)
    valid = condition.copy()
    for mean, std in zip(means, standard_deviations):
        valid &= np.isfinite(mean) & np.isfinite(std)
    selected = np.flatnonzero(valid.ravel())

    covariances = []
    selected_means = []
    selected_stds = []
    for day, mean, std in zip(DAYS, means, standard_deviations):
        stack = np.load(data_dir / f"{day}_all_data.npy", mmap_mode="r")
        flat_q = stack[:, :, :, 1].reshape(stack.shape[0], -1)[:, selected]
        flat_u = stack[:, :, :, 2].reshape(stack.shape[0], -1)[:, selected]
        theta = 0.5 * np.arctan2(flat_u, flat_q)
        mean_selected = mean.ravel()[selected]
        residual = circular_difference(theta, mean_selected)
        covariances.append(residual.T @ residual / residual.shape[0])
        selected_means.append(mean_selected)
        selected_stds.append(std.ravel()[selected])

    pair_differences = tuple(
        circular_difference(selected_means[second], selected_means[first])
        for first, second in DAY_PAIRS
    )
    pair_covariances = tuple(
        covariances[first] + covariances[second] for first, second in DAY_PAIRS
    )
    return ObservationData(
        means=tuple(selected_means),
        standard_deviations=tuple(selected_stds),
        covariances=tuple(covariances),
        pair_differences=pair_differences,
        pair_covariances=pair_covariances,
        selected_indices=selected,
    )


def select_independent_1d_indices(
    valid: NDArray[np.bool_], *, spacing: int = 10, offset: int | None = None
) -> NDArray[np.integer]:
    """Select independent angular bins for a documented residue class.

    Offset 4 reproduces the pre-existing analysis choice. Passing ``None``
    retains the historical fallback of choosing the class with most valid bins.
    """
    if offset is not None:
        if not 0 <= offset < spacing:
            raise ValueError("offset must satisfy 0 <= offset < spacing")
        return np.flatnonzero(valid & (np.arange(valid.size) % spacing == offset))
    candidates = [
        np.flatnonzero(valid & (np.arange(valid.size) % spacing == candidate))
        for candidate in range(spacing)
    ]
    return max(candidates, key=len)


def load_1d_observations(
    root: Path, data_dir: Path | None = None, *, angular_offset: int = 4
) -> ObservationData:
    data_dir = data_dir or root / "EHT-Imaging" / "data_out"
    means, standard_deviations = _load_summaries(data_dir, "_1d")
    valid = np.ones_like(means[0], dtype=bool)
    for mean, std in zip(means, standard_deviations):
        valid &= np.isfinite(mean) & np.isfinite(std)
    selected = select_independent_1d_indices(valid, offset=angular_offset)

    covariances = []
    selected_means = []
    selected_stds = []
    for day, mean, std in zip(DAYS, means, standard_deviations):
        theta = np.load(data_dir / f"{day}_EVPA_1d.npy", mmap_mode="r")[:, selected]
        mean_selected = mean[selected]
        residual = circular_difference(theta, mean_selected)
        covariances.append(residual.T @ residual / residual.shape[0])
        selected_means.append(mean_selected)
        selected_stds.append(std[selected])

    pair_differences = tuple(
        circular_difference(selected_means[second], selected_means[first])
        for first, second in DAY_PAIRS
    )
    pair_covariances = tuple(
        covariances[first] + covariances[second] for first, second in DAY_PAIRS
    )
    return ObservationData(
        means=tuple(selected_means),
        standard_deviations=tuple(selected_stds),
        covariances=tuple(covariances),
        pair_differences=pair_differences,
        pair_covariances=pair_covariances,
        selected_indices=selected,
    )


def diagonal_pair_covariances(data: ObservationData):
    return tuple(
        np.diag(data.standard_deviations[first] ** 2 + data.standard_deviations[second] ** 2)
        for first, second in DAY_PAIRS
    )
