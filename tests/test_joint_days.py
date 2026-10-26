"""Synthetic tests only: no production EHT/template inputs or driver execution."""
import numpy as np
import pytest
from scipy.special import logsumexp

from eht_axion.research.joint_days import (
    background_moments, calculate_joint_days, contrast_observation,
)


H2 = np.array([[-1., 1, 0, 0], [0, 0, -1, 1]])
H3 = np.array([[-1., 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]])


def inputs(h=H2, m=3, n=2):
    rng = np.random.default_rng(42)
    a = rng.normal(size=(4, n, n))
    return dict(
        daily_means=rng.normal(size=(4, n)),
        daily_cov=a @ a.transpose(0, 2, 1) + np.eye(n), H=h,
        signal_coeff=rng.normal(size=(m, len(h), n)) + 1j*rng.normal(size=(m, len(h), n)),
        couplings=np.array([.01, .3, 1., 2.]),
        phases=np.arange(8)*2*np.pi/8,
    )


def logpdf(y, mean, cov):
    r = y - mean
    return -.5 * (r @ np.linalg.solve(cov, r) + np.linalg.slogdet(cov)[1]
                  + len(y)*np.log(2*np.pi))


@pytest.mark.parametrize("h", [H2, H3])
def test_contrast_covariance_is_linear_transform_of_four_day_blocks(h):
    args = inputs(h)
    y, cov = contrast_observation(args['daily_means'], args['daily_cov'], h)
    block = np.zeros((8, 8))
    for day in range(4):
        block[2*day:2*day+2, 2*day:2*day+2] = args['daily_cov'][day]
    transform = np.kron(h, np.eye(2))
    np.testing.assert_allclose(y.ravel(), transform @ args['daily_means'].ravel())
    np.testing.assert_allclose(cov, transform @ block @ transform.T)
    if len(h) == 3:
        np.testing.assert_allclose(cov[:2, 2:4], args['daily_cov'][0])
    else:
        np.testing.assert_array_equal(cov[:2, 2:], 0)


@pytest.mark.parametrize("model", ['empirical', 'gaussian'])
def test_full_n2_direct_oracle_shared_latents_phase_and_gaussian_logdet(model):
    args = inputs(H3)
    rng = np.random.default_rng(3)
    b = rng.normal(size=(3, 3, 2))
    w = np.array([.2, .3, .5])
    result = calculate_joint_days(**args, backgrounds=b, weights=w, model=model)
    y = result['data'].ravel()
    cov = result['covariance']
    expected = []
    for c in args['couplings']:
        byphase = []
        for phase in args['phases']:
            states = (b + c*np.real(np.exp(1j*phase)*args['signal_coeff'])).reshape(3, -1)
            if model == 'empirical':
                value = logsumexp([logpdf(y, state, cov) + np.log(weight)
                                   for state, weight in zip(states, w)])
            else:
                mean = w @ states
                delta = states - mean
                value = logpdf(y, mean, cov + (delta.T*w) @ delta)
            byphase.append(value)
        expected.append(byphase)
    np.testing.assert_allclose(result['phase_loglike'], expected, atol=1e-12)
    np.testing.assert_allclose(result['loglike'], logsumexp(expected, axis=1)-np.log(8))
    assert result['cdf'][0] == 0 and result['cdf'][-1] == 1
    assert np.all(np.diff(result['cdf']) >= 0)
    assert .01 < result['c95'] < 2
    assert result['diagnostics']['effective_state_count'] == pytest.approx(1/(w@w))


def test_nonangular_invertible_contrast_basis_invariance():
    args = inputs(H3)
    first = calculate_joint_days(**args)
    a = np.array([[2., 1, 0], [0, 1, -.2], [.3, 0, 1]])
    args['H'] = a @ H3
    args['signal_coeff'] = np.einsum('lk,mkn->mln', a, args['signal_coeff'])
    second = calculate_joint_days(**args)
    jacobian = 2*np.linalg.slogdet(a)[1]
    np.testing.assert_allclose(second['loglike'], first['loglike'] - jacobian)
    np.testing.assert_allclose(second['cdf'], first['cdf'], atol=1e-14)
    np.testing.assert_allclose(second['component_min_chi2'], first['component_min_chi2'])


def test_legacy_mean_template_common_phase_and_exact_dft_rotation():
    args = inputs(m=1)
    args['phases'] = np.arange(64)*2*np.pi/64
    # Legacy product of pair densities, with ONE common nuisance phase.
    sigma2 = .7
    result = calculate_joint_days(**args, background_cov=sigma2*np.eye(4))
    expected = []
    for c in args['couplings']:
        terms = []
        for phase in args['phases']:
            mean = c*np.real(np.exp(1j*phase)*args['signal_coeff'][0])
            terms.append(sum(logpdf(result['data'][k], mean[k],
                                    result['covariance'][2*k:2*k+2, 2*k:2*k+2])
                             for k in range(2)))
        expected.append(logsumexp(terms)-np.log(64))
    np.testing.assert_allclose(result['loglike'], expected)
    args['signal_coeff'] *= np.exp(1j*7*2*np.pi/64)
    rotated = calculate_joint_days(**args, background_cov=sigma2*np.eye(4))
    np.testing.assert_allclose(rotated['loglike'], result['loglike'], atol=1e-13)
    # Arbitrary continuous offsets are not asserted exact for finite quadrature.


def test_shared_latent_mixture_is_not_product_of_independent_feature_mixtures():
    args = inputs(m=2, n=1)
    args['daily_means'][:] = 0
    args['signal_coeff'][:] = 0
    b = np.array([[[3.], [3.]], [[-3.], [-3.]]])
    args['daily_means'][1] = 3
    args['daily_means'][3] = -3
    result = calculate_joint_days(**args, backgrounds=b)
    y, cov = result['data'].ravel(), result['covariance']
    independent = sum(logsumexp([logpdf(y[k:k+1], state[k], cov[k:k+1, k:k+1])
                                 for state in b])-np.log(2) for k in range(2))
    assert not np.isclose(independent, result['loglike'][0])
    shuffled = calculate_joint_days(**args, backgrounds=b[:, ::-1])
    # Explicitly breaking B/T pairing in a nonzero template changes the result.
    args['signal_coeff'][:, :, :] = np.array([1, 3])[:, None, None]
    paired = calculate_joint_days(**args, backgrounds=b)
    broken = calculate_joint_days(**args, backgrounds=b[::-1])
    assert not np.allclose(paired['loglike'], broken['loglike'])
    np.testing.assert_allclose(shuffled['loglike'], result['loglike'])


def test_background_population_covariance_and_independent_gaussian_control():
    args = inputs()
    b = np.arange(12.).reshape(3, 2, 2)
    mean, cov = background_moments(b)
    np.testing.assert_allclose(cov, np.cov(b.reshape(3, 4), rowvar=False, ddof=0))
    assert cov[0, 3] > 0  # not block-diagonalized across contrasts
    for mode in ['scalar', 'diagonal', 'full']:
        mean, cov = background_moments(b, covariance=mode)
        independent = calculate_joint_days(**args, backgrounds=b, background_mode=mode)
        explicit = calculate_joint_days(**args, background_mean=mean, background_cov=cov)
        np.testing.assert_allclose(independent['loglike'], explicit['loglike'])
    centered = calculate_joint_days(**args, backgrounds=b, center_background=True)
    manual = calculate_joint_days(**args, backgrounds=b-b.mean(axis=0))
    np.testing.assert_allclose(centered['loglike'], manual['loglike'])


def test_common_phase_is_not_independent_pair_phase_marginalization():
    args = inputs(m=1, n=1)
    args['daily_cov'][:] = .05
    args['daily_means'][:] = 0
    args['daily_means'][1] = 1
    args['daily_means'][3] = -1
    args['signal_coeff'][:] = 1
    result = calculate_joint_days(**args)
    independent = []
    for c in args['couplings']:
        independent.append(sum(
            logsumexp([logpdf(result['data'][k], np.array([c*np.cos(phase)]),
                              result['covariance'][k:k+1, k:k+1])
                       for phase in args['phases']]) - np.log(len(args['phases']))
            for k in range(2)))
    assert not np.allclose(independent, result['loglike'])
    assert np.all(result['loglike'] <= result['phase_loglike'].max(axis=1) + 1e-14)
    assert np.all(result['loglike'] >= result['phase_loglike'].min(axis=1) - 1e-14)


def test_state_permutation_and_input_arrays_unchanged():
    args = inputs()
    b = np.arange(12.).reshape(3, 2, 2)
    saved_b = b.copy()
    saved_t = args['signal_coeff'].copy()
    weights = np.array([1., 2., 4.])
    result = calculate_joint_days(**args, backgrounds=b, weights=weights,
                                  background_mean=np.ones((2, 2)))
    np.testing.assert_array_equal(b, saved_b)
    np.testing.assert_array_equal(args['signal_coeff'], saved_t)
    args['signal_coeff'] = args['signal_coeff'][::-1]
    permuted = calculate_joint_days(**args, backgrounds=b[::-1], weights=weights[::-1],
                                    background_mean=np.ones((2, 2)))
    np.testing.assert_allclose(permuted['loglike'], result['loglike'])


def test_daily_response_projection_matches_direct_coefficients():
    args = inputs(H3)
    coeff = args.pop('signal_coeff')
    response = np.einsum('dk,mkn->mdn', np.linalg.pinv(H3), coeff)
    direct = calculate_joint_days(**args, signal_coeff=coeff)
    projected = calculate_joint_days(**args, responses=response)
    np.testing.assert_allclose(direct['loglike'], projected['loglike'])


@pytest.mark.parametrize('change', [
    {'couplings': [0, 1]}, {'couplings': [1, 1]}, {'phases': [np.nan]},
    {'phases': [True, False]}, {'daily_means': np.zeros((3, 2))},
    {'daily_cov': np.zeros((4, 2, 2))}, {'H': np.ones((2, 4))},
    {'H': np.array([[-1., 1, 0, 0], [-1, 1, 0, 0]])},
    {'signal_coeff': np.zeros((3, 1, 2))}, {'weights': [1, 0, 1]},
    {'background_cov': -np.eye(4)}, {'background_cov': np.ones((3, 3))},
    {'background_cov': np.array([[1, .1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])},
    {'backgrounds': np.zeros((2, 2, 2))}, {'background_mode': 'full'},
    {'backgrounds': np.zeros((3, 2, 2)), 'background_cov': np.eye(4)},
    {'model': 'invalid'}, {'prior': 'invalid'}, {'center_background': 1},
])
def test_validation_fails_without_repairs(change):
    args = inputs()
    args.update(change)
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        calculate_joint_days(**args)
