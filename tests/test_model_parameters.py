"""Explicit model choices (norm, cloud model), with small synthetic fit products."""

from dataclasses import replace

import numpy as np
import pytest

from eht_axion.analysis.likelihood import LimitOptions
from eht_axion.analysis.model import cloud_amplitude_scaling, evpa_templates
from eht_axion.analysis.observations import ObservationData
from eht_axion.analysis.likelihood import calculate_observed

TEMPLATE = dict(phase_count=7, time_conversion=0.112411, exposure_hours=6.0, day_gap_hours=24.0,
                pair_separation_hours=120.0, n_principal=2, reference_norm=0.1)
CLOUD = dict(coefficient=0.03, pivot_alpha=0.05, exponent=1.75)


def options(**changes):
    return replace(LimitOptions(
        prior="log-uniform", credibility=0.95, coupling_points=13, coupling_min=0.02, coupling_max=2.0,
        grmhd_sigma_degrees=2.0, phase_count=7, time_conversion=0.112411, exposure_hours=6.0,
        day_gap_hours=24.0, pair_separation_hours=120.0, n_principal=2, reference_norm=0.1,
        cloud_coefficient=0.03, cloud_pivot_alpha=0.05, cloud_exponent=1.75), **changes)


@pytest.fixture
def products():
    amplitude = {0.1: np.array([[0.7, 1.1], [0.08, 0.12]])}
    phase = {0.1: np.array([[0.3, -0.4], [0.04, 0.09]])}
    return amplitude, phase


@pytest.fixture
def observations():
    covariance = np.array([[0.04, 0.006], [0.006, 0.03]])
    return ObservationData(
        means=tuple(np.zeros(2) for _ in range(4)),
        standard_deviations=tuple(np.full(2, 0.2) for _ in range(4)),
        covariances=tuple(covariance / 2 for _ in range(4)),
        pair_differences=(np.array([0.05, -0.02]), np.array([0.01, 0.03])),
        pair_covariances=(covariance, covariance * 1.3),
        selected_indices=np.array([0, 1]),
    )


def test_reference_norm_scales_signal_inversely(products):
    baseline = evpa_templates(0.1, *products, **TEMPLATE)
    changed = evpa_templates(0.1, *products, **dict(TEMPLATE, reference_norm=0.25))
    assert np.any(baseline != 0)
    np.testing.assert_allclose(changed, baseline * (0.1 / 0.25), rtol=2e-15)


@pytest.mark.parametrize("alpha", [0.1, 0.25, 0.44])
def test_cloud_model_formula(alpha):
    baseline = cloud_amplitude_scaling(alpha, **CLOUD)
    assert baseline == 0.03 * (alpha / 0.05) ** 1.75
    changed = cloud_amplitude_scaling(alpha, coefficient=0.06, pivot_alpha=0.08, exponent=1.2)
    assert changed / baseline == pytest.approx(2 * (alpha / 0.08) ** 1.2 / (alpha / 0.05) ** 1.75)
    assert cloud_amplitude_scaling(alpha, **dict(CLOUD, exponent=0)) == 0.03


POSITIVE_INVALID = [0.0, -0.1, np.nan, np.inf, -np.inf, True, False, np.bool_(True), np.bool_(False)]
FINITE_INVALID = [np.nan, np.inf, -np.inf, True, False, np.bool_(True), np.bool_(False)]


@pytest.mark.parametrize("value", POSITIVE_INVALID)
def test_reference_norm_rejects_invalid_values(products, value):
    with pytest.raises(ValueError, match="reference_norm"):
        evpa_templates(0.1, *products, **dict(TEMPLATE, reference_norm=value))


@pytest.mark.parametrize("name", ["coefficient", "pivot_alpha"])
@pytest.mark.parametrize("value", POSITIVE_INVALID)
def test_cloud_positive_parameter_guards(name, value):
    with pytest.raises(ValueError, match=name):
        cloud_amplitude_scaling(0.1, **dict(CLOUD, **{name: value}))


@pytest.mark.parametrize("value", FINITE_INVALID)
def test_cloud_exponent_guards(value):
    with pytest.raises(ValueError, match="exponent"):
        cloud_amplitude_scaling(0.1, **dict(CLOUD, exponent=value))


@pytest.mark.parametrize("name,values", [
    ("reference_norm", POSITIVE_INVALID), ("cloud_coefficient", POSITIVE_INVALID),
    ("cloud_pivot_alpha", POSITIVE_INVALID), ("cloud_exponent", FINITE_INVALID),
])
def test_limit_options_reject_invalid_parameters(observations, products, name, values):
    for value in values:
        with pytest.raises(ValueError):
            calculate_observed(observations, *products, [0.1], options=options(**{name: value}))


@pytest.mark.parametrize("mode", ["std", "cor"])
@pytest.mark.parametrize("changes,ratio", [
    ({"reference_norm": 0.2}, 0.5),
    ({"cloud_coefficient": 0.06}, 2.0),
    ({"cloud_pivot_alpha": 0.1}, 2 ** -1.75),
    ({"cloud_exponent": 2.75}, 2.0),
])
def test_model_amplitude_ratio_equals_coupling_ratio(observations, products, mode, changes, ratio):
    # At fixed coupling, scaling the model amplitude by r equals scaling the coupling by r.
    base = options()
    _, changed_cut, changed = calculate_observed(observations, *products, [0.1], mode=mode,
                                                 options=replace(base, **changes))
    _, shifted_cut, shifted = calculate_observed(
        observations, *products, [0.1], mode=mode,
        options=replace(base, coupling_min=base.coupling_min * ratio, coupling_max=base.coupling_max * ratio))
    np.testing.assert_allclose(changed["log_likelihood"], shifted["log_likelihood"], rtol=3e-14, atol=2e-12)
    np.testing.assert_allclose(changed_cut, shifted_cut - np.log10(ratio), atol=2e-13)
