import numpy as np
import pytest

from scipy.integrate import quad

from eht_axion.model import evpa_templates, interpolate_fit_parameters, omega

TEMPLATE = dict(phase_count=64, time_conversion=0.112411, exposure_hours=6.0, day_gap_hours=24.0,
                pair_separation_hours=120.0, n_principal=2, reference_norm=0.1)


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
        interpolate_fit_parameters(0.15, {0.1: first, 0.2: second}, circular=False),
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


def test_phase_grid_has_no_duplicated_endpoint():
    amplitude = {0.1: np.array([[1.0], [0.2]])}
    phase = {0.1: np.array([[0.3], [0.1]])}
    signal = evpa_templates(0.1, amplitude, phase, **TEMPLATE)
    assert signal.shape == (64, 2, 1)
    assert not np.allclose(signal[0], signal[-1])


@pytest.mark.parametrize("timing", [
    (0.112411, 6.0, 120.0, 24.0, 4),
    (0.112444863, 4.2, 93.0, 26.7, 2),
])
def test_two_cosine_exposure_pair_differences_oracle(timing):
    conversion, exposure, gap, separation, principal = timing
    alpha = 0.3
    amplitude = {alpha: np.array([[1.2, 0.7], [0.13, 0.08]])}
    phase = {alpha: np.array([[0.2, -0.6], [0.07, 0.11]])}
    signal = evpa_templates(
        alpha, amplitude, phase, phase_count=9, time_conversion=conversion,
        exposure_hours=exposure, day_gap_hours=separation,
        pair_separation_hours=gap, n_principal=principal, reference_norm=0.1,
    )
    frequency = conversion * alpha * (1 - alpha**2 / (2 * principal**2))
    initial = np.linspace(0, 2 * np.pi, 9, endpoint=False)
    expected = np.empty_like(signal)
    for p, initial_phase in enumerate(initial):
        for epoch, midpoint in enumerate((0, gap)):
            for pixel in range(2):
                phi = initial_phase + phase[alpha][0, pixel]

                def average(center, angle):
                    return quad(lambda t: np.cos(frequency * t + angle),
                                center - exposure / 2, center + exposure / 2)[0] / exposure

                response = (average(midpoint - separation / 2, phi)
                            - average(midpoint + separation / 2, phi)) / 0.1
                expected[p, epoch, pixel] = amplitude[alpha][0, pixel] * response
    np.testing.assert_allclose(signal, expected, atol=2e-13, rtol=2e-13)


@pytest.mark.parametrize("alpha", [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44])
def test_omega_principal_n(alpha):
    # 211 state: radial overtone 0 + l 1 + 1 = principal N 2.
    assert omega(alpha, 2) == alpha * (1 - alpha**2 / 8)
    assert omega(alpha, 4) == alpha * (1 - alpha**2 / 32)


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
        evpa_templates(0.1, parameters, parameters, **dict(TEMPLATE, **kwargs))


@pytest.mark.parametrize("values", [
    {}, {0.1: np.ones(2)}, {0.1: np.ones((3, 1))},
    {np.nan: np.ones((2, 1))}, {0.1: np.array([[np.inf], [0.1]])},
    {0.1: np.array([[1.0], [-0.1]])},
    {0.1: np.ones((2, 1)), 0.2: np.ones((2, 2))},
])
def test_interpolation_shape_and_finite_guards(values):
    with pytest.raises(ValueError):
        interpolate_fit_parameters(0.15, values, circular=False)


def test_amplitude_phase_shape_guard():
    with pytest.raises(ValueError, match="shapes differ"):
        evpa_templates(0.1, {0.1: np.ones((2, 2))}, {0.1: np.ones((2, 1))}, **TEMPLATE)


def test_zero_exposure_is_instantaneous_limit():
    parameters = {0.1: np.array([[0.3], [0.1]])}
    signal = evpa_templates(0.1, parameters, parameters, **dict(TEMPLATE, exposure_hours=0))
    assert np.all(np.isfinite(signal))
