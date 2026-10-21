"""Focused synthetic tests; no EHT data, reconstruction, or science runs."""

import inspect
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import block_diag
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

import eht_axion.limits as null_limits
from eht_axion.limits import LimitOptions, calculate_from_products
from eht_axion.model import cloud_amplitude_scaling, evpa_templates
from eht_axion.observations import ObservationData
from eht_axion.research import calculate_observed
from eht_axion.statistics import posterior_upper_limit


def products(dimension="1d", zero=False):
    shape = (4,) if dimension == "1d" else (2, 2)
    amp = np.array([[0.7, 1.1, 0.4, 0.8], [0.3, 0.1, 0.2, 0.4]])
    phi = np.array([[0.3, -0.8, 0.5, 1.2], [0.1, 0.04, 0.2, 0.1]])
    if zero:
        amp[:] = 0
    return {0.1: amp.reshape((2,) + shape)}, {0.1: phi.reshape((2,) + shape)}


def observations(zero=False):
    covariance = np.array([[0.04, 0.006, -0.002],
                           [0.006, 0.03, 0.005], [-0.002, 0.005, 0.05]])
    return ObservationData(
        means=tuple(np.zeros(3) for _ in range(4)),
        standard_deviations=tuple(np.full(3, 0.4) for _ in range(4)),
        covariances=tuple(covariance / 2 for _ in range(4)),
        pair_differences=(np.zeros((2, 3)) if zero
                          else np.array([[0.1, -0.2, 0.05], [-0.05, 0.03, 0.2]])),
        pair_covariances=(covariance, covariance * 1.3),
        selected_indices=np.array([3, 0, 1]),  # noncontiguous and deliberately unsorted
    )


def options(**kwargs):
    return replace(LimitOptions(phase_count=9, coupling_points=23,
                                coupling_min=0.01, coupling_max=8), **kwargs)


def full_signal(amplitude, phase, opts):
    signal, _ = evpa_templates(
        0.1, amplitude, phase, phase_count=opts.phase_count,
        circular_phase_interpolation=opts.circular_phase_interpolation,
        time_conversion=opts.time_conversion, exposure_hours=opts.exposure_hours,
        day_gap_hours=opts.day_gap_hours, pair_separation_hours=opts.pair_separation_hours,
        n_principal=opts.n_principal, reference_norm=opts.reference_norm,
    )
    return signal.reshape(signal.shape[0], 2, -1) * cloud_amplitude_scaling(
        0.1, coefficient=opts.cloud_coefficient, pivot_alpha=opts.cloud_pivot_alpha,
        exponent=opts.cloud_exponent,
    )


@pytest.mark.parametrize("dimension", ["1d", "2d"])
@pytest.mark.parametrize("mode", ["std", "cor", "cor_fixed_background"])
@pytest.mark.parametrize("prior", ["uniform", "log-uniform"])
def test_dense_external_gaussian_shared_phase_oracle(dimension, mode, prior):
    obs = observations()
    amplitude, phase = products(dimension)
    opts = options(prior=prior, grmhd_sigma_degrees=3.0, reference_norm=0.15,
                   cloud_coefficient=0.04, cloud_pivot_alpha=0.07, cloud_exponent=1.2)
    alpha, cuts, diag = calculate_observed(obs, amplitude, phase, [0.1], dimension, mode, opts)
    signal = full_signal(amplitude, phase, opts)[:, :, obs.selected_indices]
    covariances = [c.copy() for c in obs.pair_covariances]
    if mode == "std":
        covariances = [np.diag(np.diag(c)) for c in covariances]
    if mode == "cor_fixed_background":
        covariances = [c + np.deg2rad(opts.grmhd_sigma_degrees)**2 * np.eye(3) for c in covariances]
    # Independent SciPy dense multivariate-normal oracle over BOTH pairs.
    covariance = block_diag(*covariances)
    data = np.asarray(obs.pair_differences).ravel()
    grid = np.asarray(diag["couplings"])
    oracle = np.array([
        [multivariate_normal.logpdf(data, mean=g * t.ravel(), cov=covariance) for t in signal]
        for g in grid
    ])
    constant = -0.5 * (data.size * np.log(2 * np.pi) + np.linalg.slogdet(covariance)[1])
    expected = logsumexp(oracle, axis=1) - np.log(opts.phase_count)
    np.testing.assert_allclose(diag["phase_log_likelihood"][0], oracle - constant,
                               rtol=2e-12, atol=1e-11)
    np.testing.assert_allclose(diag["log_likelihood"][0], expected - constant,
                               rtol=2e-12, atol=1e-11)
    assert cuts[0] == pytest.approx(np.log10(posterior_upper_limit(grid, expected, prior=prior)))
    np.testing.assert_array_equal(alpha, [0.1])
    assert diag["data_nonzero_truth"] is True
    assert diag["response_model_variance_included"] is False
    assert diag["effective_options"]["combine_epochs_by_product"] is True
    assert diag["fixed_covariance_logdet_omitted"] is True
    np.testing.assert_allclose(diag["effective_pair_covariances"], covariances)
    assert np.all(np.diff(diag["cdf"][0]) >= 0)
    assert diag["cdf"][0][-1] == 1
    json.dumps(diag, allow_nan=False)
    # A different nuisance phase per pair is NOT this model.
    independent = []
    for g in grid:
        independent.append(sum(logsumexp([
            multivariate_normal.logpdf(obs.pair_differences[e], mean=g * t[e], cov=covariances[e])
            for t in signal
        ]) - np.log(opts.phase_count) for e in range(2)))
    assert not np.allclose(expected - expected[0], np.array(independent) - independent[0])


def test_analytic_single_pixel_two_epochs_fixed_phase_quadratic_and_optimum():
    amplitude = {0.1: np.array([[0.8], [0.3]])}
    phase = {0.1: np.array([[0.6], [0.2]])}
    opts = options(phase_count=4, coupling_points=80, coupling_min=0.001, coupling_max=4)
    template = full_signal(amplitude, phase, opts)
    truth = 0.7
    data = truth * template[0]
    obs = ObservationData((), (), (), data, (np.array([[0.001]]), np.array([[0.002]])), np.array([0]))
    _, _, diag = calculate_observed(obs, amplitude, phase, [0.1], mode="cor", options=opts)
    grid = np.asarray(diag["couplings"])
    power = template[0, 0, 0]**2 / 0.001 + template[0, 1, 0]**2 / 0.002
    expected = -0.5 * power * (grid - truth)**2
    np.testing.assert_allclose(np.array(diag["phase_log_likelihood"])[0, :, 0], expected, atol=2e-12)
    fit = diag["fits"][0]
    # Analytic optimum for this fixed phase is (t^T C^-1 d)/(t^T C^-1 t).
    assert fit["data_template_cross_by_phase"][0] / fit["template_power_by_phase"][0] == pytest.approx(truth)
    assert fit["joint_grid_coupling"] == pytest.approx(grid[np.argmin(abs(grid - truth))])
    null_obs = replace(obs, pair_differences=np.zeros((2, 1)))
    _, null_cut, null_diag = calculate_observed(null_obs, amplitude, phase, [0.1], mode="cor", options=opts)
    assert null_diag["fits"][0]["joint_grid_coupling"] == pytest.approx(opts.coupling_min)
    assert diag["fits"][0]["marginal_mle_coupling"] > null_diag["fits"][0]["marginal_mle_coupling"]
    assert diag["data_nonzero_truth"] and not null_diag["data_nonzero_truth"]
    assert np.log10(truth) > null_cut[0]


@pytest.mark.parametrize("prior", ["uniform", "log-uniform"])
def test_zero_signal_analytic_finite_prior_and_actual_data_constant(prior):
    amplitude, phase = products(zero=True)
    opts = options(prior=prior, coupling_min=1, coupling_max=101, coupling_points=501)
    obs = observations()
    _, cuts, diag = calculate_observed(obs, amplitude, phase, [0.1], mode="cor", options=opts)
    chi2 = sum(d @ np.linalg.solve(c, d) for d, c in zip(obs.pair_differences, obs.pair_covariances))
    np.testing.assert_allclose(diag["log_likelihood"][0], -chi2 / 2, atol=1e-14)
    analytic = 96 if prior == "uniform" else 101**0.95
    assert 10**cuts[0] == pytest.approx(analytic, rel=2e-5)
    zero = replace(obs, pair_differences=np.zeros((2, 3)))
    _, zero_cuts, zero_diag = calculate_observed(zero, amplitude, phase, [0.1], mode="cor", options=opts)
    np.testing.assert_allclose(cuts, zero_cuts, atol=1e-14)
    np.testing.assert_allclose(zero_diag["log_likelihood"], 0, atol=1e-14)


@pytest.mark.parametrize("dimension", ["1d", "2d"])
@pytest.mark.parametrize("mode", ["std", "cor"])
@pytest.mark.parametrize("std_covariance", ["covariance", "legacy-sigma"])
def test_zero_data_matches_unchanged_corrected_null_api(dimension, mode, std_covariance):
    obs = observations(zero=True)
    amplitude, phase = products(dimension)
    opts = options(std_covariance=std_covariance)
    alpha, cuts, diag = calculate_observed(obs, amplitude, phase, [0.1], dimension, mode, opts)
    null_alpha, null_cuts, null_diag = calculate_from_products(
        obs, amplitude, phase, [0.1], dimension, mode, opts,
    )
    np.testing.assert_array_equal(alpha, null_alpha)
    np.testing.assert_allclose(cuts, null_cuts, rtol=0, atol=1e-13)
    for key in ("log_likelihood", "cdf", "posterior_density"):
        np.testing.assert_allclose(diag[key], null_diag[key], rtol=2e-12, atol=1e-11)


def test_fixed_background_loosen_in_specific_zero_data_example_only():
    amplitude, phase = products()
    obs = observations(zero=True)
    opts = options(coupling_points=101, coupling_max=100, grmhd_sigma_degrees=10)
    _, plain, _ = calculate_observed(obs, amplitude, phase, [0.1], mode="cor", options=opts)
    _, background, diag = calculate_observed(obs, amplitude, phase, [0.1],
                                            mode="cor_fixed_background", options=opts)
    # This is a particular zero-data example, not a monotonicity claim for real data.
    assert background[0] > plain[0]
    _, no_background, _ = calculate_observed(
        obs, amplitude, phase, [0.1], mode="cor_fixed_background",
        options=replace(opts, grmhd_sigma_degrees=0),
    )
    np.testing.assert_allclose(plain, no_background, atol=1e-14)
    assert diag["background_variance_radians2"] == pytest.approx(np.deg2rad(10)**2)


def test_selected_before_validation_and_response_uncertainties_excluded():
    obs = observations()
    amplitude, phase = products()
    _, cuts, diag = calculate_observed(obs, amplitude, phase, [0.1], options=options())
    # The excluded full-grid pixel is allowed to be NaN; uncertainties of used
    # pixels remain validated but must NOT add coupling-dependent variance.
    amplitude[0.1][:, 2] = np.nan
    phase[0.1][:, 2] = np.nan
    amplitude[0.1][1, obs.selected_indices] = 7
    phase[0.1][1, obs.selected_indices] = 4
    _, same_cuts, same_diag = calculate_observed(obs, amplitude, phase, [0.1], options=options())
    np.testing.assert_array_equal(cuts, same_cuts)
    np.testing.assert_array_equal(diag["log_likelihood"], same_diag["log_likelihood"])


def test_gaussian_residual_is_not_rewrapped(monkeypatch):
    import eht_axion.research.observed as observed

    signal = np.full((2, 2, 3), 2.0)
    monkeypatch.setattr(observed, "evpa_templates", lambda *a, **k: (signal, np.zeros_like(signal)))
    obs = replace(observations(), pair_differences=np.full((2, 3), -1.5))
    amplitude, phase = products()
    opts = options(cloud_coefficient=1, cloud_exponent=0, coupling_min=1, coupling_max=2)
    _, _, diag = calculate_observed(obs, amplitude, phase, [0.1], mode="cor", options=opts)
    expected = sum(multivariate_normal.logpdf(np.full(3, -3.5), cov=c)
                   - multivariate_normal.logpdf(np.zeros(3), cov=c) for c in obs.pair_covariances)
    assert diag["log_likelihood"][0][0] == pytest.approx(expected)
    assert any("not rewrapped" in warning for warning in diag["warnings"])


@pytest.mark.parametrize("change", [
    {"legacy": True}, {"legacy": None}, {"legacy": "false"},
    {"combine_epochs_by_product": False}, {"duplicate_phase_endpoint": True},
    {"include_parameter_logdet": "yes"}, {"posterior_method": "legacy"},
    {"prior": "bad"}, {"coupling_points": 1}, {"coupling_points": True},
    {"coupling_min": 0}, {"coupling_max": np.inf}, {"coupling_min": True},
    {"std_covariance": "bad"}, {"grmhd_sigma_degrees": -1},
    {"phase_count": 1}, {"phase_count": 2.5}, {"n_principal": False},
    {"time_conversion": 0}, {"exposure_hours": -1}, {"day_gap_hours": 0},
    {"pair_separation_hours": -1}, {"reference_norm": 0},
    {"cloud_coefficient": 0}, {"cloud_pivot_alpha": -1}, {"cloud_exponent": np.nan},
])
def test_invalid_options(change):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_observed(observations(), amplitude, phase, [0.1], options=options(**change))


@pytest.mark.parametrize("change", [
    {"selected_indices": np.array([])}, {"selected_indices": np.array([0.5, 1.5, 2.5])},
    {"selected_indices": np.array([0, 0, 1])}, {"selected_indices": np.array([-1, 0, 1])},
    {"selected_indices": np.array([0, 1, 4])},
    {"pair_differences": np.full((2, 3), np.nan)}, {"pair_differences": np.zeros((3, 2))},
    {"pair_differences": np.zeros((1, 3))}, {"pair_covariances": (np.eye(3),)},
    {"pair_covariances": (np.eye(2), np.eye(2))},
    {"pair_covariances": (np.full((3, 3), np.inf), np.eye(3))},
    {"pair_covariances": (np.array([[1, 0.1, 0], [0, 1, 0], [0, 0, 1]]), np.eye(3))},
    {"pair_covariances": (np.ones((3, 3)), np.eye(3))},
    {"pair_covariances": (np.array([[1, 2, 0], [2, 1, 0], [0, 0, 1]]), np.eye(3))},
])
@pytest.mark.parametrize("mode", ["std", "cor_fixed_background"])
def test_invalid_observations_including_non_spd_even_with_background(change, mode):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_observed(replace(observations(), **change), amplitude, phase, [0.1], mode=mode)


@pytest.mark.parametrize("alphas", [[], [np.nan], [0], [0.2, 0.1], [0.1, 0.1], [[0.1]]])
def test_invalid_alpha_grid(alphas):
    amplitude, phase = products()
    with pytest.raises(ValueError, match="alphas"):
        calculate_observed(observations(), amplitude, phase, alphas)


@pytest.mark.parametrize("kind", ["empty", "key", "different_keys", "ndim", "shape", "finite", "std"])
def test_invalid_fit_products(kind):
    amplitude, phase = products()
    if kind == "empty":
        amplitude = {}
    elif kind == "key":
        amplitude = {np.nan: amplitude[0.1]}
    elif kind == "different_keys":
        phase = {0.2: phase[0.1]}
    elif kind == "ndim":
        amplitude = {0.1: np.ones((2, 2, 2))}
    elif kind == "shape":
        amplitude = {0.1: np.ones((3, 4))}
    elif kind == "finite":
        amplitude[0.1][0, 0] = np.inf
    else:
        phase[0.1][1, 0] = -1
    with pytest.raises(ValueError):
        calculate_observed(observations(), amplitude, phase, [0.1])


@pytest.mark.parametrize("arguments", [{"mode": "cor_grmhd"}, {"mode": "bad"}, {"dimension": "3d"}])
def test_invalid_request(arguments):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_observed(observations(), amplitude, phase, [0.1], **arguments)


def test_multiple_alphas_interpolation_and_resolved_options():
    amplitude, phase = products("2d")
    amplitude[0.2] = amplitude[0.1] * 1.3
    phase[0.2] = phase[0.1] + np.array([[[2.0]], [[0.01]]])
    opts = options(phase_count=None, circular_phase_interpolation=np.bool_(False))
    obs = observations(zero=True)
    alpha, cuts, diag = calculate_observed(obs, amplitude, phase, [0.1, 0.15, 0.2],
                                           dimension="2d", mode="cor", options=opts)
    _, null_cuts, null_diag = calculate_from_products(obs, amplitude, phase, alpha,
                                                     dimension="2d", mode="cor", options=opts)
    np.testing.assert_allclose(cuts, null_cuts, atol=1e-13)
    np.testing.assert_allclose(diag["log_likelihood"], null_diag["log_likelihood"], atol=1e-11)
    assert np.shape(diag["phase_log_likelihood"]) == (3, 23, 64)
    assert diag["effective_options"]["phase_count"] == 64
    assert diag["effective_options"]["circular_phase_interpolation"] is False
    assert diag["effective_options"]["include_parameter_logdet"] is False
    json.dumps(diag, allow_nan=False)


@pytest.mark.parametrize("stds", [(), (np.ones(3),) * 3, (np.ones(2),) * 4,
                                      (np.full(3, -1.0),) * 4, (np.full(3, np.nan),) * 4,
                                      (np.zeros(3),) * 4])
def test_invalid_historical_daily_sigmas(stds):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_observed(replace(observations(), standard_deviations=stds), amplitude, phase,
                           [0.1], options=options(std_covariance="legacy-sigma"))


def test_option_type_and_null_observed_option_remain_rejected():
    amplitude, phase = products()
    with pytest.raises(TypeError, match="LimitOptions"):
        calculate_observed(observations(), amplitude, phase, [0.1], options={})
    with pytest.raises(TypeError, match="observed"):
        LimitOptions(observed=True)
    with pytest.raises(TypeError, match="observed"):
        calculate_from_products(observations(), amplitude, phase, [0.1], observed=True)
    assert "observed" not in inspect.signature(calculate_from_products).parameters


def test_no_overwrite_or_null_function_replacement(tmp_path, monkeypatch):
    # Snapshot NULL implementation/API/tests. This test never writes to them.
    root = Path(null_limits.__file__).resolve().parent.parent
    paths = [root / name for name in ("eht_axion/limits.py", "tests/test_limits.py",
                                      "tests/test_null_only_results.py")]
    before = {path: path.read_bytes() for path in paths}
    original_calculation = null_limits.calculate_from_products
    monkeypatch.chdir(tmp_path)
    obs = observations()
    amplitude, phase = products()
    original_data = np.array(obs.pair_differences, copy=True)
    original_covariances = np.array(obs.pair_covariances, copy=True)
    calculate_observed(obs, amplitude, phase, [0.1], options=options())
    calculate_observed(obs, amplitude, phase, [0.1], mode="cor_fixed_background", options=options())
    assert list(tmp_path.iterdir()) == []
    assert null_limits.calculate_from_products is original_calculation
    assert {path: path.read_bytes() for path in paths} == before
    np.testing.assert_array_equal(obs.pair_differences, original_data)
    np.testing.assert_array_equal(obs.pair_covariances, original_covariances)
    _, actual_null, _ = original_calculation(obs, amplitude, phase, [0.1], options=options())
    _, zero_null, _ = original_calculation(observations(zero=True), amplitude, phase,
                                           [0.1], options=options())
    np.testing.assert_array_equal(actual_null, zero_null)
