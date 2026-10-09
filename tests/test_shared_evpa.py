"""RAPTOR/transport share validity rules without changing zero-flux synthesis."""

import numpy as np
from eht_axion.pipeline.raptor import _evpa
from eht_axion.pipeline.transport import analytical_stokes


def test_raptor_nonfinite_linear_stokes_is_invalid():
    image = np.ones((2, 2, 4))
    for q, u in ((np.inf, 1), (1, -np.inf), (np.nan, 1), (0, 0)):
        image[..., 1], image[..., 2] = q, u
        assert np.isnan(_evpa(image)).all()


def test_screen_zero_polarization_remains_zero_stokes():
    baseline = np.zeros((2, 2, 4))
    baseline[..., 0] = 1
    field = np.full((16, 2, 2), .1)
    result = analytical_stokes(baseline, field)
    np.testing.assert_array_equal(result[..., 1:3], 0)
    assert np.isfinite(result).all()


def test_screen_nonfinite_polarization_never_becomes_finite():
    baseline = np.ones((2, 2, 4))
    baseline[..., 1] = np.inf
    with np.errstate(invalid="ignore"):
        result = analytical_stokes(baseline, np.zeros((16, 2, 2)))
    assert np.isnan(_evpa(result)).all()
