"""Selected EHT observations: daily circular moments and the two day-pair differences."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ObservationData:
    means: tuple[NDArray[np.floating], ...]
    standard_deviations: tuple[NDArray[np.floating], ...]
    covariances: tuple[NDArray[np.floating], ...]
    pair_differences: tuple[NDArray[np.floating], ...]
    pair_covariances: tuple[NDArray[np.floating], ...]
    selected_indices: NDArray[np.integer]
