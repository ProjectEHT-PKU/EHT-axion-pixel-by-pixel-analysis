import numpy as np
import pytest

from eht_axion.polarization import circular_window, integrated_evpa


def test_circular_window_has_requested_width_and_wraps():
    identity = np.eye(12)
    result = circular_window(identity, 4)
    np.testing.assert_allclose(result.sum(axis=1), 4)
    assert np.count_nonzero(result[0]) == 4
    assert result[0, -2] == 1
    assert result[0, -1] == 1
    assert result[0, 0] == 1
    assert result[0, 1] == 1


def test_integrated_evpa_does_not_depend_on_total_intensity():
    # A tiny explicit linear operator is enough to verify Stokes semantics.
    operator = np.array([[1, 1, 1, 1]], dtype=float)
    q = np.ones((2, 2))
    u = np.ones((2, 2))
    result = integrated_evpa(q, u, window_degrees=1, operator=operator)
    assert result[0] == pytest.approx(np.pi / 8)


def test_integrated_evpa_zero_and_nonfinite_integrals_are_undefined():
    operator = np.array([[1, 1, 1, 1]], dtype=float)
    for q, u in ((np.zeros((2, 2)), np.zeros((2, 2))),
                 (np.full((2, 2), np.inf), np.ones((2, 2))),
                 (np.ones((2, 2)), np.full((2, 2), np.nan))):
        assert np.isnan(integrated_evpa(q, u, window_degrees=1, operator=operator)).all()


def test_integrated_evpa_accepts_batches():
    operator = np.array([[1, 0, 0, 0], [0, 0, 0, 1]], dtype=float)
    q = np.ones((3, 2, 2))
    u = np.zeros_like(q)
    result = integrated_evpa(q, u, window_degrees=1, operator=operator)
    assert result.shape == (3, 2)
    np.testing.assert_allclose(result, 0)
