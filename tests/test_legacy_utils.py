import numpy as np
import pytest

pytest.importorskip("statsmodels", reason="optional legacy utilities: install .[legacy]")
from utils.libanaly import CircleCalc, ImageS, fit_gaussian2


def test_circle_sub_handles_scalar_evpa():
    assert CircleCalc.circle_sub_evpa(np.deg2rad(89), np.deg2rad(-89)) == pytest.approx(
        np.deg2rad(-2)
    )


def test_circle_sub_non_symmetric_degree_range():
    result = CircleCalc.circle_sub(np.array([10.0]), np.array([350.0]), 0, 360)
    np.testing.assert_allclose(result, [20.0])


def test_image_data_is_not_corrupted_after_filtering():
    data = np.zeros((8, 8, 4))
    data[..., 0] = 1
    image = ImageS(
        data,
        gaussian_kernel=0,
        rotation_angle=0,
        init_fov=8,
        final_fov=8,
        final_shape=8,
        calc_1D=False,
    )
    assert isinstance(image.data, np.ndarray)
    assert image.data.shape == data.shape


def test_fit_gaussian2_honors_requested_normality_test():
    clearly_non_normal = np.repeat([0.0, 1.0], 50)
    result = fit_gaussian2(clearly_non_normal, plot_failed=False, p_value_threshold=0.05)
    assert np.isnan(result).all()


def test_no_blur_honors_rotation_angle():
    data = np.zeros((8, 8, 4))
    data[..., 0] = 1
    data[..., 1] = 1
    image = ImageS(
        data,
        blur=False,
        rotation_angle=0,
        init_fov=8,
        final_fov=8,
        final_shape=8,
        calc_1D=False,
    )
    np.testing.assert_allclose(image.U, 0, atol=1e-12)
    assert np.all(image.Q > 0.99)
