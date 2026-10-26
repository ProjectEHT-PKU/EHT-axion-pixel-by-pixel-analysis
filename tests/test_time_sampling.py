import numpy as np
import pytest

from eht_axion.research.time_sampling import (
    DAY_HOURS, H2, H3, average_modulated_response, contrast_covariance,
    evpa_contrasts, four_day_stokes, integrate_stokes_at,
)


def test_linear_stokes_and_piecewise_knot():
    t = np.array([0., 4., 8., 12.])
    s = (2 * t + 3)[:, None, None] * np.ones((1, 2, 4))
    np.testing.assert_allclose(integrate_stokes_at(t, s, [3., 9.]), s[[0, 0]] + np.array([6., 18.])[:, None, None])
    np.testing.assert_allclose(integrate_stokes_at([0, 4, 8], np.array([0., 4., 0.]), 4., 8.), 2.)


@pytest.mark.parametrize('omega', [0., 1e-10, .05, -1.7])
def test_constant_response_sinc(omega):
    t = np.arange(0., 21., 2.)
    z = np.full((len(t), 2), 2 + 3j)
    centers = np.array([4., 13.])
    expected = np.exp(1j * omega * centers)[:, None] * z[0] * np.sinc(omega * 3 / np.pi)
    np.testing.assert_allclose(average_modulated_response(centers, t, z, omega), expected, atol=2e-14)


def test_linear_modulation_against_quadrature():
    from scipy.integrate import quad_vec
    t = np.array([0., 4., 9., 15.])
    z = (2 + 1j) + t[:, None] * np.array([.3 - .2j, -.1 + .5j])
    omega = .7
    expected = quad_vec(lambda x: ((2 + 1j) + x * np.array([.3 - .2j, -.1 + .5j])) * np.exp(1j * omega * x), 2, 8)[0] / 6
    np.testing.assert_allclose(average_modulated_response(5., t, z, omega), expected, atol=1e-13)


@pytest.mark.parametrize('times,center,width', [([0, 8, 4], 4, 2), ([0, 4, 8], 1, 6), ([0, 4, 8], 4, -1), ([0, 4, 8], np.nan, 0)])
def test_no_extrapolation(times, center, width):
    with pytest.raises(ValueError):
        integrate_stokes_at(times, np.ones((3, 4)), center, width)


def test_windows_and_zero_width():
    t = np.arange(141) * 8.89324754
    states = list(four_day_stokes(t, t[:, None]))
    assert len(states) == 124
    np.testing.assert_allclose(states[0][1][:, 0], 3 + DAY_HOURS)
    np.testing.assert_allclose(integrate_stokes_at(t, t, [0, t[-1]], 0), [0, t[-1]])


def test_covariance_and_basis_invariance():
    c = np.array([np.diag([j, j + 1]) for j in range(1, 5)])
    full = np.zeros((8, 8))
    for j in range(4):
        full[2*j:2*j+2, 2*j:2*j+2] = c[j]
    k = np.kron(H3, np.eye(2))
    cov = contrast_covariance(c)
    np.testing.assert_allclose(cov, k @ full @ k.T)
    np.testing.assert_allclose(cov[:2, 2:4], -c[1])
    np.testing.assert_allclose(cov[2:4, 4:], -c[2])
    np.testing.assert_allclose(contrast_covariance(c, H2), cov[np.ix_([0, 1, 4, 5], [0, 1, 4, 5])])
    a = np.kron(np.array([[1., 1, 0], [0, 1, 1], [1, 0, 0]]), np.eye(2))
    x = np.arange(6.)
    np.testing.assert_allclose(x @ np.linalg.solve(cov, x), (a @ x) @ np.linalg.solve(a @ cov @ a.T, a @ x))


def test_spin_two_winding_is_reported_not_excluded():
    chi = np.array([[0.], [.4*np.pi], [.8*np.pi], [.1*np.pi]])
    d = evpa_contrasts(chi)
    np.testing.assert_allclose(np.exp(2j*d['raw']), np.exp(2j*d['principal']))
    assert np.any(d['winding'] != 0)
    assert np.any(d['closure_winding'] != 0)
