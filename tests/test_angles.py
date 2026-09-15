import numpy as np
import pytest

from eht_axion.angles import circular_difference, circular_mean, unwrap_around


def test_evpa_difference_wraps_with_period_pi():
    actual = circular_difference(np.deg2rad(89), np.deg2rad(-89))
    assert np.rad2deg(actual) == pytest.approx(-2)


def test_difference_supports_scalars_and_broadcasting():
    assert circular_difference(1.5, -1.5) == pytest.approx(3 - np.pi)
    result = circular_difference(np.array([0.1, 0.2]), 0.1)
    np.testing.assert_allclose(result, [0, 0.1])


def test_circular_mean_straddles_evpa_boundary():
    result = circular_mean(np.deg2rad([89, -89]))
    assert abs(abs(np.rad2deg(result)) - 90) < 1e-12


def test_unwrap_around_center():
    values = unwrap_around(np.deg2rad([89, -89]), np.deg2rad(89))
    np.testing.assert_allclose(np.rad2deg(values), [89, 91], atol=1e-12)


def test_rejects_invalid_period():
    with pytest.raises(ValueError):
        circular_difference(1, 2, period=0)
