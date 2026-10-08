import numpy as np
import pytest

from scipy.integrate import quad

from eht_axion.model import TIME_CONVERSION, evpa_templates, interpolate_fit_parameters, omega


def test_phase_interpolation_crosses_branch_cut_not_zero():
    first = np.array([[np.deg2rad(179.0)], [0.1]])
    second = np.array([[np.deg2rad(-179.0)], [0.3]])
    result = interpolate_fit_parameters(0.15, {0.1: first, 0.2: second}, circular=True)
    assert abs(abs(np.rad2deg(result[0, 0])) - 180) < 1e-10
    assert result[1, 0] == pytest.approx(0.2)


def test_linear_interpolation_for_amplitudes():
    first = np.array([[1.0], [2.0]])
    second = np.array([[3.0], [4.0]])
    np.testing.assert_allclose(
        interpolate_fit_parameters(0.15, {0.1: first, 0.2: second}),
        [[2.0], [3.0]],
    )


@pytest.mark.parametrize("node", [0.1, 0.2, 0.3])
def test_exact_phase_nodes_preserve_values_not_wrapped(node):
    values = {alpha: np.array([[5.0 + alpha, -7.0], [0.1, 0.3]])
              for alpha in (0.1, 0.2, 0.3)}
    result = interpolate_fit_parameters(node, values, circular=True)
    np.testing.assert_array_equal(result, values[node])
    result[0, 0] = 0
    assert values[node][0, 0] != 0


@pytest.mark.parametrize("legacy, count", [(False, 64), (True, 16)])
def test_default_phase_grid_and_legacy_shared_uncertainty(legacy, count):
    amplitude = {0.1: np.array([[1.0], [0.2]])}
    phase = {0.1: np.array([[0.3], [0.1]])}
    signal, error = evpa_templates(0.1, amplitude, phase, legacy=legacy)
    assert signal.shape == error.shape == (count, 2, 1)
    assert np.allclose(signal[0], signal[-1]) == legacy
    assert np.allclose(error[:, 0], error[:, 1]) == legacy


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("timing", [
    (TIME_CONVERSION, 6.0, 120.0, 24.0, 4),
    (0.112444863, 4.2, 93.0, 26.7, 2),
])
def test_two_cosine_exposure_pair_differences_oracle(legacy, timing):
    conversion, exposure, gap, separation, principal = timing
    alpha = 0.3
    amplitude = {alpha: np.array([[1.2, 0.7], [0.13, 0.08]])}
    phase = {alpha: np.array([[0.2, -0.6], [0.07, 0.11]])}
    signal, error = evpa_templates(
        alpha, amplitude, phase, legacy=legacy, phase_count=9,
        legacy_shared_uncertainty=False, time_conversion=conversion,
        exposure_hours=exposure, day_gap_hours=separation,
        pair_separation_hours=gap, n_principal=principal,
    )
    frequency = conversion * alpha * (1 - alpha**2 / (2 * principal**2))
    initial = np.linspace(0, 2 * np.pi, 9, endpoint=legacy)
    effective_gap = gap / 2 if legacy else gap
    expected = np.empty_like(signal)
    expected_error = np.empty_like(error)
    for p, initial_phase in enumerate(initial):
        for epoch, midpoint in enumerate((0, effective_gap)):
            for pixel in range(2):
                phi = initial_phase + phase[alpha][0, pixel]

                def average(center, angle):
                    return quad(lambda t: np.cos(frequency * t + angle),
                                center - exposure / 2, center + exposure / 2)[0] / exposure

                def difference(angle):
                    return (average(midpoint - separation / 2, angle)
                            - average(midpoint + separation / 2, angle)) / 0.1

                response = difference(phi)
                derivative = (difference(phi + 1e-5) - difference(phi - 1e-5)) / 2e-5
                expected[p, epoch, pixel] = amplitude[alpha][0, pixel] * response
                expected_error[p, epoch, pixel] = np.hypot(
                    amplitude[alpha][1, pixel] * response,
                    amplitude[alpha][0, pixel] * phase[alpha][1, pixel] * derivative,
                )
    np.testing.assert_allclose(signal, expected, atol=2e-13, rtol=2e-13)
    np.testing.assert_allclose(error, expected_error, atol=2e-10, rtol=2e-10)
    if legacy:
        physical, _ = evpa_templates(
            alpha, amplitude, phase, phase_count=9, duplicate_endpoint=True,
            time_conversion=conversion, exposure_hours=exposure, day_gap_hours=separation,
            pair_separation_hours=gap, n_principal=principal,
        )
        np.testing.assert_allclose(signal[:, 0], physical[:, 0])
        assert not np.allclose(signal[:, 1], physical[:, 1])


def test_omega_keeps_historical_n4_and_supports_explicit_n():
    assert omega(0.3) == 0.3 * (1 - 0.3**2 / 32)
    assert omega(0.3, 2) == 0.3 * (1 - 0.3**2 / 8)


@pytest.mark.parametrize("kwargs", [
    {"phase_count": 1}, {"phase_count": 2.5}, {"phase_count": True},
    {"time_conversion": np.nan}, {"exposure_hours": -1},
    {"day_gap_hours": np.inf}, {"day_gap_hours": 0},
    {"pair_separation_hours": -1},
    {"n_principal": 0}, {"n_principal": 2.0},
])
def test_template_parameter_guards(kwargs):
    parameters = {0.1: np.array([[1.0], [0.2]])}
    with pytest.raises(ValueError):
        evpa_templates(0.1, parameters, parameters, **kwargs)


@pytest.mark.parametrize("values", [
    {}, {0.1: np.ones(2)}, {0.1: np.ones((3, 1))},
    {np.nan: np.ones((2, 1))}, {0.1: np.array([[np.inf], [0.1]])},
    {0.1: np.array([[1.0], [-0.1]])},
    {0.1: np.ones((2, 1)), 0.2: np.ones((2, 2))},
])
def test_interpolation_shape_and_finite_guards(values):
    with pytest.raises(ValueError):
        interpolate_fit_parameters(0.15, values)


def test_amplitude_phase_shape_guard():
    with pytest.raises(ValueError, match="shapes differ"):
        evpa_templates(0.1, {0.1: np.ones((2, 2))}, {0.1: np.ones((2, 1))})


def test_zero_exposure_is_instantaneous_limit():
    parameters = {0.1: np.array([[0.3], [0.1]])}
    signal, error = evpa_templates(0.1, parameters, parameters, exposure_hours=0)
    assert np.all(np.isfinite(signal)) and np.all(np.isfinite(error))
