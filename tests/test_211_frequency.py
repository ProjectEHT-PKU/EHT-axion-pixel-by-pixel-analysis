"""Article 211 means principal N=2, radial overtone=0; m is not in N."""
import inspect
from pathlib import Path

import numpy as np
import pytest

from eht_axion.limits import LimitOptions
from eht_axion.model import evpa_templates, omega
from eht_axion.pipeline.config import read_config


@pytest.mark.parametrize("alpha", [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44])
def test_default_frequency_is_principal_two(alpha):
    radial_overtone, orbital_l = 0, 1
    principal = radial_overtone + orbital_l + 1
    expected = alpha * (1 - alpha**2 / (2 * principal**2))
    assert omega(alpha) == expected
    # This is the old literal (2+1+1), corresponding to principal four.
    assert omega(alpha, 4) == alpha * (1 - alpha**2 / (2 * (2 + 1 + 1)**2))
    assert omega(alpha) != omega(alpha, 4)


def test_research_and_all_public_defaults_agree_on_principal_two():
    root = Path(__file__).resolve().parents[1]
    assert read_config(root / "configs/research.yaml")["inference"]["principal_n"] == 2
    assert LimitOptions().n_principal == 2
    assert inspect.signature(omega).parameters["n_principal"].default == 2
    assert inspect.signature(evpa_templates).parameters["n_principal"].default == 2


def test_template_defaults_use_two_without_changing_other_options():
    alpha = 0.3
    amplitude = {alpha: np.array([[0.02, 0.03], [0.001, 0.002]])}
    phase = {alpha: np.array([[0.1, -0.2], [0.02, 0.03]])}
    default = evpa_templates(alpha, amplitude, phase)
    explicit = evpa_templates(alpha, amplitude, phase, n_principal=2)
    old = evpa_templates(alpha, amplitude, phase, n_principal=4)
    for actual, expected in zip(default, explicit):
        np.testing.assert_array_equal(actual, expected)
    assert not np.allclose(default[0], old[0])
    assert not np.allclose(default[1], old[1])
