import numpy as np
import pytest

from eht_axion.model import interpolate_fit_parameters


def test_phase_interpolation_crosses_branch_cut_not_zero():
    first = np.array([[np.deg2rad(179.0)], [0.1]])
    second = np.array([[np.deg2rad(-179.0)], [0.3]])
    result = interpolate_fit_parameters(0.15, {0.1: first, 0.2: second}, circular=True)
    assert abs(abs(np.rad2deg(result[0, 0])) - 180) < 1e-10
    assert result[1, 0] == pytest.approx(0.2)


def test_linear_interpolation_for_amplitudes():
    first = np.array([[1.0], [2.0]])
    second = np.array([[3.0], [4.0]])
    np.testing.assert_allclose(
        interpolate_fit_parameters(0.15, {0.1: first, 0.2: second}),
        [[2.0], [3.0]],
    )
