"""State concentration and contrast-basis checks on finite Gaussian mixtures."""
import numpy as np

from diagnose_joint_day_states import state_weights
from check_joint_day_state_support import subset_states


def test_state_weights_recover_identical_mixture_and_basis():
    data = np.array([[.2, -.1], [.1, .3], [-.2, .1]])
    cov = np.diag(np.arange(1, 7)/10)
    templates = np.ones((3, 3, 2), complex)*(.2+.1j)
    backgrounds = np.zeros((3, 3, 2))
    weights = np.array([.2, .3, .5])
    r, q, ll = state_weights(data, cov, templates, backgrounds, .7, .3, weights)
    np.testing.assert_allclose(r, weights)
    transform = np.tril(np.ones((3, 3)))
    big = np.kron(transform, np.eye(2))
    rr, qq, lll = state_weights(transform@data, big@cov@big.T,
        np.einsum('ij,mjn->min', transform, templates),
        np.einsum('ij,mjn->min', transform, backgrounds), .7, .3, weights)
    np.testing.assert_allclose(rr, r, atol=1e-14)
    np.testing.assert_allclose(qq, q, atol=1e-14)
    assert abs(ll-lll) < 1e-13


def test_state_subset_keeps_background_response_pairing():
    s = {'evolving_daily': np.arange(8*124*4*2).reshape(8, 124, 4, 2),
         'background_h3_principal': np.arange(124*3*2).reshape(124, 3, 2)}
    t, b, w = subset_states(s, 62, 124)
    np.testing.assert_array_equal(t[:, 0], s['evolving_daily'][:, 62])
    np.testing.assert_array_equal(b[0], s['background_h3_principal'][62])
    assert len(w) == 62 and abs(w.sum()-1) < 1e-14
