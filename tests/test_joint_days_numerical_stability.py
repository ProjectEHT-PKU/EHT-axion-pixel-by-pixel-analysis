"""Direct-residual numerical regressions; synthetic inputs, no RT/EHT I/O."""
import warnings

import numpy as np
import pytest
from scipy.stats import multivariate_normal

from eht_axion.research.joint_days import calculate_joint_days


def one_feature(*, observed=1., variance=1e-16, template=1. + 1e-8):
    return dict(
        daily_means=np.array([[0.], [observed], [0.], [0.]]),
        daily_cov=np.full((4, 1, 1), variance / 2),
        H=np.array([[-1., 1., 0., 0.]]),
        signal_coeff=np.array([[[template]]], complex),
        couplings=np.array([1., 1. + 1e-8]), phases=np.array([0.]),
    )


@pytest.mark.parametrize('variance,template', [
    (1e-16, 1. + 1e-8),
    # Baseline polynomial is -0.5 at the first node and was clipped to zero.
    (3e-16, 1. + 5.878016072274924e-9),
])
def test_empirical_cancellation_matches_direct_dense_gaussian(variance, template):
    args = one_feature(variance=variance, template=template)
    result = calculate_joint_days(**args)
    expected = [multivariate_normal.logpdf(
        result['data'].ravel(), mean=(coupling*args['signal_coeff'][0].real).ravel(),
        cov=result['covariance']) for coupling in args['couplings']]
    # Whitened subtraction has O(eps*sqrt(q0)) error in the tiny residual,
    # but must not lose its O(1) likelihood contribution as the polynomial did.
    np.testing.assert_allclose(result['phase_loglike'][:, 0], expected, rtol=0, atol=3e-8)
    assert np.all(result['component_min_chi2'] > 0)


def test_null_template_large_finite_couplings_has_constant_finite_likelihood():
    args = one_feature(observed=0., variance=2., template=0.)
    args['couplings'] = np.array([1e155, 2e155])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        result = calculate_joint_days(**args)
    expected = multivariate_normal.logpdf([0.], mean=[0.], cov=[[2.]])
    np.testing.assert_allclose(result['loglike'], [expected, expected], rtol=0, atol=1e-14)
    np.testing.assert_array_equal(result['component_min_chi2'], 0)
    assert np.isfinite(result['c95'])


def test_large_coupling_small_signal_avoids_square_before_product():
    args = one_feature(observed=0., variance=2., template=1e-155)
    args['couplings'] = np.array([1e155, 2e155])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        result = calculate_joint_days(**args)
    expected = [multivariate_normal.logpdf([0.], mean=[mean], cov=[[2.]]) for mean in (1., 2.)]
    np.testing.assert_allclose(result['loglike'], expected, rtol=0, atol=1e-14)


def test_nonfinite_scaled_signal_is_explicit_error_not_clipped():
    args = one_feature(observed=0., variance=2., template=1e155)
    args['couplings'] = np.array([1e155, 2e155])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(ValueError, match='coupling-scaled whitened signal is nonfinite'):
            calculate_joint_days(**args)


def test_nonfinite_squared_residual_is_explicit_error_not_clipped():
    args = one_feature(observed=1e155, variance=2., template=0.)
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(ValueError, match='residual quadratic is nonfinite'):
            calculate_joint_days(**args)
