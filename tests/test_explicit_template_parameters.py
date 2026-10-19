"""Explicit model choices, tested only with small synthetic fit products."""

from dataclasses import replace

import numpy as np
import pytest

from eht_axion.limits import LimitOptions, calculate_from_products
from eht_axion.model import AXION_NORM, cloud_amplitude_scaling, evpa_templates
from eht_axion.observations import ObservationData


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
        pair_differences=(np.zeros(2), np.zeros(2)),
        pair_covariances=(covariance, covariance * 1.3),
        selected_indices=np.array([0, 1]),
    )


@pytest.mark.parametrize("legacy", [False, True])
def test_reference_norm_scales_signal_and_sigma_inversely(products, legacy):
    baseline = evpa_templates(0.1, *products, phase_count=7, legacy=legacy)
    changed = evpa_templates(
        0.1, *products, phase_count=7, legacy=legacy, reference_norm=0.25,
    )
    for original, scaled in zip(baseline, changed):
        assert np.any(original != 0)
        np.testing.assert_allclose(scaled, original * (0.1 / 0.25), rtol=2e-15)


@pytest.mark.parametrize("alpha", [0.1, 0.25, 0.44])
def test_cloud_historical_default_and_explicit_model_ratio(alpha):
    baseline = cloud_amplitude_scaling(alpha)
    assert baseline == 0.03 * (alpha / 0.05) ** (7 / 4)
    changed = cloud_amplitude_scaling(
        alpha, coefficient=0.06, pivot_alpha=0.08, exponent=1.2,
    )
    ratio = 2 * (alpha / 0.08) ** 1.2 / (alpha / 0.05) ** 1.75
    assert changed / baseline == pytest.approx(ratio)
    assert cloud_amplitude_scaling(alpha, exponent=0) == 0.03
    assert cloud_amplitude_scaling(alpha, exponent=-1) == pytest.approx(0.03 * 0.05 / alpha)


POSITIVE_INVALID = [0.0, -0.1, np.nan, np.inf, -np.inf, True, False,
                    np.bool_(True), np.bool_(False)]
FINITE_INVALID = [np.nan, np.inf, -np.inf, True, False, np.bool_(True), np.bool_(False)]


@pytest.mark.parametrize("value", POSITIVE_INVALID)
def test_reference_norm_rejects_invalid_values(products, value):
    with pytest.raises(ValueError, match="reference_norm"):
        evpa_templates(0.1, *products, reference_norm=value)


@pytest.mark.parametrize("name", ["alpha", "coefficient", "pivot_alpha"])
@pytest.mark.parametrize("value", POSITIVE_INVALID)
def test_cloud_positive_parameter_guards(name, value):
    kwargs = {"alpha": 0.1, name: value}
    with pytest.raises(ValueError, match=name):
        cloud_amplitude_scaling(**kwargs)


@pytest.mark.parametrize("value", FINITE_INVALID)
def test_cloud_exponent_guards(value):
    with pytest.raises(ValueError, match="exponent"):
        cloud_amplitude_scaling(0.1, exponent=value)


@pytest.mark.parametrize("name,values", [
    ("reference_norm", POSITIVE_INVALID),
    ("cloud_coefficient", POSITIVE_INVALID),
    ("cloud_pivot_alpha", POSITIVE_INVALID),
    ("cloud_exponent", FINITE_INVALID),
])
def test_limit_options_forward_invalid_parameters(observations, products, name, values):
    for value in values:
        options = LimitOptions(phase_count=5, coupling_points=7, **{name: value})
        with pytest.raises(ValueError, match=name.removeprefix("cloud_")):
            calculate_from_products(observations, *products, [0.1], options=options)


@pytest.mark.parametrize("legacy", [False, True])
def test_explicit_historical_defaults_preserve_results(observations, products, legacy):
    options = LimitOptions(legacy=legacy, phase_count=7, coupling_points=13)
    assert options.reference_norm == AXION_NORM == 0.1
    assert (options.cloud_coefficient, options.cloud_pivot_alpha, options.cloud_exponent) == (
        0.03, 0.05, 1.75,
    )
    explicit = replace(options, reference_norm=0.1, cloud_coefficient=0.03,
                       cloud_pivot_alpha=0.05, cloud_exponent=7 / 4)
    default_result = calculate_from_products(observations, *products, [0.1], options=options)
    explicit_result = calculate_from_products(observations, *products, [0.1], options=explicit)
    np.testing.assert_array_equal(default_result[1], explicit_result[1])
    assert default_result[2] == explicit_result[2]


@pytest.mark.parametrize("mode", ["std", "cor_grmhd"])
@pytest.mark.parametrize("changes,ratio", [
    ({"reference_norm": 0.2}, 0.5),
    ({"cloud_coefficient": 0.06}, 2.0),
    ({"cloud_pivot_alpha": 0.1}, 2 ** -1.75),
    ({"cloud_exponent": 2.75}, 2.0),
    ({"reference_norm": 0.2, "cloud_coefficient": 0.045,
      "cloud_pivot_alpha": 0.08, "cloud_exponent": 1.2},
     0.5 * 1.5 * (0.1 / 0.08) ** 1.2 / 2 ** 1.75),
])
def test_model_amplitude_ratio_equals_coupling_ratio(
    observations, products, mode, changes, ratio,
):
    # At fixed coupling, changing the full model amplitude by r is exactly
    # equivalent to multiplying coupling by r, including propagated sigma.
    options = LimitOptions(phase_count=7, coupling_points=13,
                           coupling_min=0.02, coupling_max=2.0)
    _, changed_cut, changed = calculate_from_products(
        observations, *products, [0.1], mode=mode, options=replace(options, **changes),
    )
    _, shifted_cut, shifted = calculate_from_products(
        observations, *products, [0.1], mode=mode,
        options=replace(options, coupling_min=options.coupling_min * ratio,
                        coupling_max=options.coupling_max * ratio),
    )
    np.testing.assert_allclose(changed["log_likelihood"], shifted["log_likelihood"],
                               rtol=3e-14, atol=2e-12)
    np.testing.assert_allclose(changed_cut, shifted_cut - np.log10(ratio), atol=2e-13)
