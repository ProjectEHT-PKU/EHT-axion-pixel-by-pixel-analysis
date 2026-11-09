"""evpa: EVPA validity and pi-periodic differences."""
import numpy as np
import pytest

from eht_axion.evpa import circular_difference, finite_evpa


def test_evpa_difference_wraps_with_period_pi():
    actual = circular_difference(np.deg2rad(89), np.deg2rad(-89))
    assert np.rad2deg(actual) == pytest.approx(-2)


def test_difference_supports_scalars_and_broadcasting():
    assert circular_difference(1.5, -1.5) == pytest.approx(3 - np.pi)
    result = circular_difference(np.array([0.1, 0.2]), 0.1)
    np.testing.assert_allclose(result, [0, 0.1])


def test_rejects_invalid_period():
    with pytest.raises(ValueError):
        circular_difference(1, 2, period=0)


def test_finite_evpa_rejects_nonfinite_and_zero_but_not_finite_extremes():
    tiny = np.nextafter(0.0, 1.0)
    big = np.finfo(np.float64).max
    q = np.array([0, -0.0, np.nan, np.inf, 1, tiny, big, -big])
    u = np.array([0, -0.0, 1, 1, -np.inf, tiny, big, big])
    actual = finite_evpa(q, u)
    assert np.isnan(actual[:5]).all()
    np.testing.assert_array_equal(actual[5:], [np.pi / 8, np.pi / 8, 3 * np.pi / 8])
    np.testing.assert_array_equal(finite_evpa([[1], [-1]], [0, 1]),
                                  [[0, np.pi / 8], [np.pi / 2, 3 * np.pi / 8]])
