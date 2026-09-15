import numpy as np

from eht_axion.limits import _legacy_upper_limit
from eht_axion.model import evpa_templates


def test_legacy_upper_limit_is_bounded_by_grid():
    grid = np.logspace(-2, 2, 60)
    result = _legacy_upper_limit(grid, np.exp(-grid))
    assert np.log10(grid[0]) <= result <= np.log10(grid[-1])


def test_templates_use_unique_default_phases_and_epoch_specific_errors():
    amplitude = {0.1: np.array([[1.0], [0.2]])}
    phase = {0.1: np.array([[0.3], [0.1]])}
    signal, uncertainty = evpa_templates(0.1, amplitude, phase)
    assert signal.shape == (16, 2, 1)
    assert uncertainty.shape == (16, 2, 1)
    # 16 unique phases have no duplicated first/last signal.
    assert not np.allclose(signal[0], signal[-1])
    # The five-day phase shift must also affect propagated uncertainty.
    assert not np.allclose(uncertainty[:, 0], uncertainty[:, 1])
