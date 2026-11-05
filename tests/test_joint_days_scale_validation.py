"""Scale-relative validation and unresolved exposure-width regressions."""
import warnings

import numpy as np
import pytest

from eht_axion.research.joint_days import calculate_joint_days, contrast_observation
from eht_axion.research.time_sampling import (
    average_modulated_response, four_day_stokes, integrate_stokes_at,
)


@pytest.mark.parametrize('scale', [1e-20, 1., 1e20])
def test_asymmetric_daily_covariance_rejected_independently_of_units(scale):
    # At scale=1e-20 even large relative asymmetry passed the old absolute atol.
    cov = np.broadcast_to(scale * np.array([[1., .4], [0., 1.]]), (4, 2, 2))
    with pytest.raises(ValueError, match='symmetric relative'):
        contrast_observation(np.zeros((4, 2)), cov, [[-1., 1., 0., 0.]])


def test_original_tiny_asymmetric_covariance_reproduction_is_rejected():
    cov = np.broadcast_to(np.array([[1e-16, 4e-15], [0., 1e-16]]), (4, 2, 2))
    with pytest.raises(ValueError, match='daily_cov block must be symmetric'):
        contrast_observation(np.zeros((4, 2)), cov, [[-1., 1., 0., 0.]])


@pytest.mark.parametrize('scale', [1e-20, 1., 1e20])
def test_legitimate_scaled_spd_covariance_accepted_without_modification(scale):
    cov = np.broadcast_to(scale * np.array([[2., .3], [.3, 1.]]), (4, 2, 2)).copy()
    saved = cov.copy()
    _, joint = contrast_observation(np.zeros((4, 2)), cov, [[-1., 1., 0., 0.]])
    np.testing.assert_allclose(joint, 2 * cov[0], rtol=1e-14, atol=0)
    np.testing.assert_array_equal(cov, saved)


@pytest.mark.parametrize('scale', [1e-20, 1., 1e20])
def test_roundoff_relative_asymmetry_accepted_without_symmetrization(scale):
    cov = np.broadcast_to(scale * np.array([[2., .3 + 1e-13], [.3, 1.]]), (4, 2, 2))
    _, joint = contrast_observation(np.zeros((4, 2)), cov, [[-1., 1., 0., 0.]])
    np.testing.assert_array_equal(joint, 2 * cov[0])
    assert joint[0, 1] != joint[1, 0]  # tolerance is validation, never a repair


def test_asymmetric_fixed_covariance_also_uses_relative_spd_guard():
    # Leave the separately documented background_cov PSD tolerance unchanged;
    # fixed total covariance must still pass the stricter SPD symmetry guard.
    with pytest.raises(ValueError, match='fixed covariance must be symmetric'):
        calculate_joint_days(
            np.zeros((4, 2)), np.broadcast_to(np.eye(2) * 1e-20, (4, 2, 2)),
            [[-1., 1., 0., 0.]], signal_coeff=np.zeros((1, 1, 2)),
            background_cov=np.array([[0., 1e-15], [0., 0.]]),
            couplings=[.1, 1.], phases=[0.])


@pytest.mark.parametrize('scale', [1e-20, 1., 1e20])
def test_noncontrast_row_cannot_pass_zero_sum_test_by_rescaling(scale):
    with pytest.raises(ValueError, match='independent zero-sum rows'):
        contrast_observation(np.ones((4, 1)), np.ones((4, 1, 1)), [[scale, 0., 0., 0.]])


@pytest.mark.parametrize('scale', [1e-20, 1., 1e20])
def test_legitimate_scaled_contrast_row_is_accepted(scale):
    data, covariance = contrast_observation(
        np.arange(4.)[:, None], np.ones((4, 1, 1)), [[-scale, scale, 0., 0.]])
    np.testing.assert_allclose(data, [[scale]], rtol=1e-14, atol=0)
    np.testing.assert_allclose(covariance, [[2 * scale**2]], rtol=1e-14, atol=0)


def test_zero_contrast_row_rejected_without_divide_by_zero_warning():
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(ValueError, match='independent zero-sum rows'):
            contrast_observation(np.ones((4, 1)), np.ones((4, 1, 1)), [[0., 0., 0., 0.]])


def test_zero_covariance_still_fails_cholesky_naturally():
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(np.linalg.LinAlgError):
            contrast_observation(np.ones((4, 1)), np.zeros((4, 1, 1)), [[-1., 1., 0., 0.]])


@pytest.mark.parametrize('modulated', [False, True])
def test_positive_exposure_width_below_time_resolution_rejected(modulated):
    if modulated:
        with pytest.raises(ValueError, match='exposure width below time resolution'):
            average_modulated_response([1000.], [0., 2000.], np.ones(2, complex), .1, width=1e-14)
    else:
        with pytest.raises(ValueError, match='exposure width below time resolution'):
            integrate_stokes_at([0., 2000.], np.ones(2), [1000.], width=1e-14)


def test_any_collapsed_center_in_multidimensional_request_is_rejected():
    with pytest.raises(ValueError, match='exposure width below time resolution'):
        integrate_stokes_at([0., 2000.], np.ones(2), [[1., 1000.]], width=1e-14)


def test_four_day_stokes_uses_same_collapsed_width_guard():
    with pytest.raises(ValueError, match='exposure width below time resolution'):
        list(four_day_stokes([0., 2000.], np.ones(2), width=1e-14, starts=[1000.]))


def test_zero_width_retains_point_interpolation_and_modulation():
    np.testing.assert_array_equal(integrate_stokes_at([0., 2000.], [1., 1.], [1000.], width=0), [1.])
    actual = average_modulated_response([1000.], [0., 2000.], np.ones(2, complex), .1, width=0)
    np.testing.assert_allclose(actual, [np.exp(100j)], rtol=1e-14, atol=0)


def test_six_hour_constant_exposure_is_unchanged():
    np.testing.assert_array_equal(integrate_stokes_at([0., 2000.], [1., 1.], [1000.], width=6.), [1.])
    actual = average_modulated_response([1000.], [0., 2000.], np.ones(2, complex), .1, width=6.)
    np.testing.assert_allclose(actual, [np.exp(100j) * np.sinc(.3 / np.pi)], rtol=1e-13, atol=0)
