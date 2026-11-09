"""Load EHT topset products and construct EVPA covariance matrices."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


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


def diagonal_pair_covariances(data: ObservationData):
    return tuple(
        np.diag(data.standard_deviations[first] ** 2 + data.standard_deviations[second] ** 2)
        for first, second in DAY_PAIRS
    )
