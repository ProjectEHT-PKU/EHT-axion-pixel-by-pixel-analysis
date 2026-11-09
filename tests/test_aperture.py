"""analysis.aperture: quadrature operator and intensity-weighted aperture EVPA validity rules."""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from eht_axion.analysis.aperture import aperture_operator, integrate_stokes


@pytest.mark.parametrize("sparse", [False, True])
def test_padding_is_zero_flux_and_unused_v_does_not_change_evpa(sparse):
    stokes = np.ones((2, 2, 2, 4))
    stokes[..., 0] = 2
    stokes[:, 0, 0, :] = np.nan
    operator = np.ones((1, 4))
    if sparse:
        operator = csr_matrix(operator)
    reference = integrate_stokes(stokes, operator)
    np.testing.assert_array_equal(reference, np.full((2, 1), np.pi / 8))
    stokes[:, 1, 1, 3] = np.nan
    stokes[:, 1, 0, 3] = np.inf
    np.testing.assert_array_equal(integrate_stokes(stokes, operator), reference)
    assert np.isnan(integrate_stokes(np.full_like(stokes, np.nan), operator)).all()


@pytest.mark.parametrize("channel", [0, 1, 2])
@pytest.mark.parametrize("nonfinite", [np.nan, np.inf, -np.inf])
def test_unknown_individual_iqu_is_rejected_even_with_zero_operator(channel, nonfinite):
    stokes = np.ones((2, 2, 4))
    stokes[0, 0, channel] = nonfinite
    with pytest.raises(ValueError, match="Nonfinite I/Q/U.*padding"):
        integrate_stokes(stokes, np.zeros((1, 4)))


@pytest.mark.parametrize("sparse", [False, True])
def test_finite_float64_observable_matches_original_formula(sparse):
    stokes = np.random.default_rng(827).normal(size=(3, 2, 2, 4))
    operator = np.array([[1.0, 0.25, 0.75, 0.5], [0.0, 1.0, 0.5, 2.0]])
    weights = stokes[..., 0]
    q = (stokes[..., 1] * weights).reshape(3, 4) @ operator.T
    u = (stokes[..., 2] * weights).reshape(3, 4) @ operator.T
    expected = 0.5 * np.arctan2(u, q)
    actual = integrate_stokes(stokes, csr_matrix(operator) if sparse else operator)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-16)


@pytest.mark.parametrize("sparse", [False, True])
def test_integer_weighting_uses_float64_products_not_wrapping(sparse):
    stokes = np.zeros((2, 2, 4), dtype=np.int64)
    stokes[..., 0] = 2**32
    stokes[..., 1] = 2**32
    stokes[..., 2] = 2**31
    operator = np.ones((1, 4))
    if sparse:
        operator = csr_matrix(operator)
    expected = 0.5 * np.arctan2(0.5, 1.0)
    np.testing.assert_array_equal(integrate_stokes(stokes, operator), [expected])
    np.testing.assert_array_equal(integrate_stokes(stokes, operator),
                                  integrate_stokes(stokes.astype(np.float64), operator))


@pytest.mark.parametrize("channel", [1, 2])
def test_finite_inputs_with_overflowing_weighted_product_are_rejected(channel):
    stokes = np.ones((2, 2, 4))
    stokes[..., 0] = 1e200
    stokes[..., channel] = 1e200
    with pytest.raises(ValueError, match="Nonfinite float64 intensity-weighted Stokes product"):
        integrate_stokes(stokes, np.ones((1, 4)))


@pytest.mark.parametrize("sparse", [False, True])
def test_zero_integrated_polarization_is_undefined(sparse):
    stokes = np.ones((2, 2, 4))
    stokes[:, 1, 1:3] = -1
    operator = np.ones((1, 4))
    if sparse:
        operator = csr_matrix(operator)
    assert np.isnan(integrate_stokes(stokes, operator)).all()


def test_nonfinite_integrated_qu_cannot_become_finite_angle():
    stokes = np.ones((2, 2, 4))
    stokes[..., 1] = 1e308
    with np.errstate(over="ignore", invalid="ignore"):
        assert np.isnan(integrate_stokes(stokes, np.ones((1, 4)))).all()


def test_aperture_area_exact_and_intensity_weighting():
    operator = aperture_operator(
        shape=4, fov=4, radial_max=1.9, angular_bins=4, window_degrees=10
    )
    area = 0.5 * 1.9**2 * np.deg2rad(10)
    np.testing.assert_allclose(np.asarray(operator.sum(axis=1)).ravel(), area, rtol=1e-12)
    image = np.zeros((4, 4, 4))
    image[..., 0] = 3
    image[..., 1] = 1
    image[..., 2] = 1
    np.testing.assert_allclose(integrate_stokes(image, operator), np.pi / 8)
    image[...] = 0
    assert np.isnan(integrate_stokes(image, operator)).all()


def test_window_center_symmetry():
    op = aperture_operator(
        shape=16, fov=4, radial_max=1.5, angular_bins=360, window_degrees=10
    )
    # x-positive aperture is exactly symmetric under reflection y -> -y.
    row = op[0].toarray().reshape(16, 16)
    np.testing.assert_allclose(row, row[:, ::-1], atol=1e-14)
