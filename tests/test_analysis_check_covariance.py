"""Finite-MC sensitivity estimates scatter, not uncertainty of a sample mean."""
import numpy as np

from analysis_check_covariance import covariance


def test_centering_wrap_and_no_extra_sample_count_factor():
    samples = np.array([[89., 5.], [-89., 7.], [88., 9.], [-88., 11.]])
    center = np.array([90., 8.])
    actual = covariance(np.deg2rad(samples), np.deg2rad(center))
    residual = np.deg2rad(np.array([[-1., -3.], [1., -1.], [-2., 1.], [2., 3.]]))
    expected = residual.T@residual/len(residual)
    np.testing.assert_allclose(actual, expected, atol=1e-15)
    repeated = covariance(np.repeat(np.deg2rad(samples), 2, axis=0), np.deg2rad(center))
    np.testing.assert_allclose(repeated, actual, atol=1e-15)
