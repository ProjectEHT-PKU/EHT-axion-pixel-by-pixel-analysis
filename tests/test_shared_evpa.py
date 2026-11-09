"""RAPTOR EVPA validity rules."""

import numpy as np
from eht_axion.pipeline.raptor import _evpa


def test_raptor_nonfinite_linear_stokes_is_invalid():
    image = np.ones((2, 2, 4))
    for q, u in ((np.inf, 1), (1, -np.inf), (np.nan, 1), (0, 0)):
        image[..., 1], image[..., 2] = q, u
        assert np.isnan(_evpa(image)).all()
