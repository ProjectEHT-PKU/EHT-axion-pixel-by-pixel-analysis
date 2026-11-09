"""Article 211 means principal N=2, radial overtone=0; m is not in N."""
import dataclasses
import inspect

import numpy as np
import pytest

from eht_axion.limits import LimitOptions
from eht_axion.model import evpa_templates, omega
from eht_axion.experiment import read_experiment


@pytest.mark.parametrize("alpha", [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44])
def test_default_frequency_is_principal_two(alpha):
    radial_overtone, orbital_l = 0, 1
    principal = radial_overtone + orbital_l + 1
    expected = alpha * (1 - alpha**2 / (2 * principal**2))
    assert omega(alpha, principal) == expected
    # This is the old literal (2+1+1), corresponding to principal four.
    assert omega(alpha, 4) == alpha * (1 - alpha**2 / (2 * (2 + 1 + 1)**2))
    assert omega(alpha, principal) != omega(alpha, 4)


def test_experiment_uses_principal_two_and_no_hidden_defaults():
    assert read_experiment()["shared"]["principal_n"] == 2
    for function in (omega, evpa_templates):
        assert inspect.signature(function).parameters["n_principal"].default is inspect.Parameter.empty
    assert LimitOptions.__dataclass_fields__["n_principal"].default is dataclasses.MISSING


def test_templates_depend_on_principal_n():
    alpha = 0.3
    amplitude = {alpha: np.array([[0.02, 0.03], [0.001, 0.002]])}
    phase = {alpha: np.array([[0.1, -0.2], [0.02, 0.03]])}
    common = dict(phase_count=64, time_conversion=0.112411, exposure_hours=6.0, day_gap_hours=24.0,
                  pair_separation_hours=120.0, reference_norm=0.1)
    two = evpa_templates(alpha, amplitude, phase, n_principal=2, **common)
    old = evpa_templates(alpha, amplitude, phase, n_principal=4, **common)
    assert not np.allclose(two[0], old[0])
    assert not np.allclose(two[1], old[1])
