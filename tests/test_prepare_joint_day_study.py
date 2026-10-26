"""New-input preparation adapters; no source data/RT required."""

import numpy as np
import pytest

from eht_axion.pipeline.raptor import _integrate, _operator, fit_phase_dft
from eht_axion.research.state_ensemble import (
    local_branches,
    response_days,
    selected_background,
)
from eht_axion.research.time_sampling import DAY_HOURS, H3, integrate_stokes_at


def test_constant_response_frozen_evolving_and_phase_origin():
    times = np.arange(0.0, 201.0, 10.0)
    z = np.broadcast_to(np.array([0.4 + 0.2j, -0.1 + 0.8j]), (len(times), 2)).copy()
    starts = np.array([3.0, 13.0, 23.0])
    w = 0.12
    frozen, evolving = response_days(times, z, starts, w)
    expected = (
        z[0, None, :] * np.exp(1j * w * DAY_HOURS)[:, None] * np.sinc(w * 3 / np.pi)
    )
    np.testing.assert_allclose(
        frozen, np.broadcast_to(expected, frozen.shape), atol=1e-14
    )
    np.testing.assert_allclose(
        evolving, np.broadcast_to(expected, evolving.shape), atol=1e-14
    )
    assert frozen.shape == (21, 4, 2)
    assert evolving.shape == (3, 4, 2)


def test_dft_sign_is_positive_complex_phase():
    psi = 2 * np.pi * np.arange(16) / 16
    phase = 0.7
    fit = fit_phase_dft(0.05 + 0.02 * np.cos(psi + phase))
    np.testing.assert_allclose(fit[:2], [0.02, phase], atol=1e-14)
    predicted = np.real(fit[0] * np.exp(1j * fit[1]) * np.exp(1j * psi))
    np.testing.assert_allclose(predicted, 0.02 * np.cos(psi + phase), atol=1e-14)


def test_stokes_average_before_evpa_and_operator():
    t = np.arange(0.0, 201.0, 10.0)
    cube = np.zeros((len(t), 1, 4, 4, 4))
    cube[..., 0] = 2
    cube[..., 1] = (1 + 0.01 * t)[:, None, None, None]
    cube[..., 2] = (0.2 + 0.03 * t)[:, None, None, None]
    centers = np.array([3.0, 27.0, 123.0, 147.0])
    days = integrate_stokes_at(t, cube, centers, 6)
    indices = np.array([0, 5, 15])
    selected = selected_background(days, indices, "2d")
    expected = 0.5 * np.arctan2(0.2 + 0.03 * centers, 1 + 0.01 * centers)
    np.testing.assert_allclose(selected, np.broadcast_to(expected[:, None], (4, 3)))
    op = _operator(4, 4, None, 2, 10, 360)
    selected1 = selected_background(days, np.array([4, 14]), "1d", op)
    np.testing.assert_allclose(
        selected1, _integrate(days[:, 0], op, "intensity")[:, [4, 14]]
    )


def test_local_lift_consistent_but_reports_long_gap_ambiguity():
    chi = np.array([[0.0], [0.4], [0.8], [1.2]]) * np.pi
    lifted, diagnostic = local_branches(chi)
    np.testing.assert_allclose(H3 @ lifted, diagnostic["principal"])
    assert np.count_nonzero(diagnostic["closure_winding"]) > 0
    assert (abs(lifted - lifted[:1]) > np.pi / 2).any()


def test_selected_nonfinite_never_imputed():
    days = np.zeros((4, 2, 2, 4))
    with pytest.raises(ValueError, match="nonfinite selected"):
        selected_background(days, np.array([0]), "2d")
    with pytest.raises(ValueError, match="nonfinite selected"):
        response_days(
            np.array([0.0, 200.0]), np.array([[np.nan], [1j]]), np.array([3.0]), 0.1
        )


def test_shared_state_preserves_correlated_columns():
    t = np.arange(0.0, 201.0, 10.0)
    first = (1 + 0.001 * t) * np.exp(0.002j * t)
    coeff = np.stack([first, 3j * first], axis=-1)
    f, e = response_days(t, coeff, np.array([3.0, 13.0, 23.0]), 0.08)
    np.testing.assert_allclose(f[..., 1], 3j * f[..., 0])
    np.testing.assert_allclose(e[..., 1], 3j * e[..., 0])
