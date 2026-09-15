import numpy as np
import pytest

from eht_axion.observations import select_independent_1d_indices


def test_explicit_angular_offset_is_reproducible():
    valid = np.ones(40, dtype=bool)
    np.testing.assert_array_equal(
        select_independent_1d_indices(valid, spacing=10, offset=4), [4, 14, 24, 34]
    )


def test_angular_offset_validation():
    with pytest.raises(ValueError):
        select_independent_1d_indices(np.ones(20, dtype=bool), spacing=10, offset=10)
