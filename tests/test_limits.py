import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.special import logsumexp

from eht_axion.limits import (
    ALPHAS, LimitOptions, _legacy_upper_limit, calculate_from_products, calculate_limits,
)
from eht_axion.model import cloud_amplitude_scaling, evpa_templates
from eht_axion.observations import ObservationData
from eht_axion.statistics import gaussian_log_likelihood, posterior_upper_limit


def products(dimension="1d", zero_signal=False, zero_errors=False):
    shape = (4,) if dimension == "1d" else (2, 2)
    amplitude = np.array([[1.0, 0.7, 1.1, 0.4], [0.1, 0.2, 0.0, 0.0]])
    phase = np.array([[0.3, -0.8, 4.0, 0.5], [0.03, 0.05, 0.0, 0.0]])
    if zero_signal:
        amplitude[:] = 0
    if zero_errors:
        amplitude[1] = phase[1] = 0
    return ({0.1: amplitude.reshape((2,) + shape)}, {0.1: phase.reshape((2,) + shape)})


def observations(diagonal=False):
    covariance = np.array([[0.04, 0.006, -0.002],
                           [0.006, 0.03, 0.005], [-0.002, 0.005, 0.05]])
    if diagonal:
        covariance = np.diag(np.diag(covariance))
    # Deliberately different historical sigmas and covariance diagonals.
    return ObservationData(
        means=tuple(np.zeros(3) for _ in range(4)),
        standard_deviations=tuple(np.full(3, 0.4) for _ in range(4)),
        covariances=tuple(covariance / 2 for _ in range(4)),
        pair_differences=(np.array([0.1, -0.2, 0.05]), np.array([-0.05, 0.03, 0.2])),
        pair_covariances=(covariance, covariance * 1.3),
        selected_indices=np.array([0, 1, 3]),
    )


def direct_null_loglikes(obs, amplitude, phase, dimension, mode, options, grid):
    """Independent null residual/Cholesky oracle, never the fast quadratic/eigensolver."""
    signal, error = evpa_templates(
        0.1, amplitude, phase, legacy=options.legacy, phase_count=options.phase_count,
        time_conversion=options.time_conversion, exposure_hours=options.exposure_hours,
        day_gap_hours=options.day_gap_hours, pair_separation_hours=options.pair_separation_hours,
        n_principal=options.n_principal,
    )
    signal = signal.reshape(signal.shape[0], 2, -1)[:, :, obs.selected_indices]
    error = error.reshape(error.shape[0], 2, -1)[:, :, obs.selected_indices]
    scale = cloud_amplitude_scaling(0.1)
    result = []
    for coupling in grid:
        phase_likes = []
        for p in range(signal.shape[0]):
            epoch_likes = []
            for epoch in range(2):
                covariance = obs.pair_covariances[epoch].copy()
                if mode == "std":
                    if options.std_covariance == "legacy-sigma":
                        covariance = np.diag(obs.standard_deviations[2 * epoch]**2
                                             + obs.standard_deviations[2 * epoch + 1]**2)
                    else:
                        covariance = np.diag(np.diag(covariance))
                if mode == "cor_grmhd":
                    covariance[np.diag_indices(3)] += (
                        np.deg2rad(options.grmhd_sigma_degrees)**2
                        + (error[p, epoch] * coupling * scale)**2
                    )
                residual = signal[p, epoch] * coupling * scale
                include = (not options.legacy if options.include_parameter_logdet is None
                           else options.include_parameter_logdet)
                epoch_likes.append(gaussian_log_likelihood(
                    residual, covariance, include_logdet=include and mode == "cor_grmhd",
                ))
            phase_likes.append(logsumexp(epoch_likes) if options.legacy else sum(epoch_likes))
        result.append(logsumexp(phase_likes) - np.log(signal.shape[0]))
    return np.array(result)


@pytest.mark.parametrize("dimension", ["1d", "2d"])
@pytest.mark.parametrize("mode", ["std", "cor", "cor_grmhd"])
@pytest.mark.parametrize("legacy", [False, True])
def test_new_api_likelihood_and_limits_match_direct_null_gaussian(dimension, mode, legacy):
    obs = observations()
    amplitude, phase = products(dimension)
    options = LimitOptions(legacy=legacy, coupling_points=19,
                           coupling_min=0.01, coupling_max=30, phase_count=7)
    alphas, cuts, diagnostics = calculate_from_products(
        obs, amplitude, phase, [0.1], dimension, mode, options,
    )
    expected = direct_null_loglikes(obs, amplitude, phase, dimension, mode, options,
                               diagnostics["couplings"])
    np.testing.assert_allclose(diagnostics["log_likelihood"][0], expected, atol=2e-11, rtol=2e-12)
    if legacy:
        cut = _legacy_upper_limit(diagnostics["couplings"], np.exp(expected - expected.max()))
    else:
        cut = np.log10(posterior_upper_limit(diagnostics["couplings"], expected))
    np.testing.assert_allclose(cuts, [cut], atol=1e-12, rtol=0)
    np.testing.assert_array_equal(alphas, [0.1])
    assert np.shape(diagnostics["cdf"]) == (1, 19)
    assert diagnostics["cdf"][0][-1] == pytest.approx(1.0)
    json.dumps(diagnostics, allow_nan=False)
    assert np.all(np.diff(diagnostics["cdf"][0]) >= 0)
    assert diagnostics["boundaries"][0]["lower_bin_mass"] >= 0


@pytest.mark.parametrize("diagonal", [False, True])
@pytest.mark.parametrize("zero_errors", [False, True])
@pytest.mark.parametrize("include_logdet", [False, True])
def test_exact_heteroscedastic_eigen_oracle_zero_eigenvalues(diagonal, zero_errors, include_logdet):
    obs = observations(diagonal)
    amplitude, phase = products(zero_errors=zero_errors)
    options = LimitOptions(phase_count=8, coupling_points=31,
                           coupling_min=1e-4, coupling_max=1e3, grmhd_sigma_degrees=0.8,
                           include_parameter_logdet=include_logdet)
    _, _, diagnostic = calculate_from_products(obs, amplitude, phase, [0.1],
                                               "1d", "cor_grmhd", options)
    expected = direct_null_loglikes(obs, amplitude, phase, "1d", "cor_grmhd", options,
                               diagnostic["couplings"])
    np.testing.assert_allclose(diagnostic["log_likelihood"][0], expected, rtol=3e-12, atol=2e-9)


@pytest.mark.parametrize("mode", ["std", "cor", "cor_grmhd"])
def test_new_api_synthetic_analytic_constant_likelihood(mode):
    obs = observations()
    amplitude, phase = products(zero_signal=True)
    options = LimitOptions(prior="uniform", coupling_min=1, coupling_max=101,
                           coupling_points=9)
    _, cuts, diagnostic = calculate_from_products(obs, amplitude, phase, [0.1],
                                                 "1d", mode, options)
    # With a zero template and zero propagated model errors, the uniform-prior
    # posterior is exactly constant; diagnostic nonzero differences are ignored.
    assert cuts[0] == pytest.approx(np.log10(96.0), abs=1e-14)
    assert diagnostic["phase_count"] == 64
    np.testing.assert_allclose(diagnostic["cdf"][0], (np.array(diagnostic["couplings"]) - 1) / 100)
    np.testing.assert_allclose(diagnostic["log_likelihood"][0],
                               diagnostic["log_likelihood"][0][0], atol=1e-14)


@pytest.mark.parametrize("dimension", ["1d", "2d"])
@pytest.mark.parametrize("mode", ["std", "cor", "cor_grmhd"])
@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("diagonal", [False, True])
def test_nonzero_pair_differences_do_not_affect_null_inference(dimension, mode, legacy, diagonal):
    obs = observations(diagonal)
    amplitude, phase = products(dimension)
    options = LimitOptions(legacy=legacy, phase_count=7, coupling_points=19)
    _, cuts, diagnostic = calculate_from_products(
        obs, amplitude, phase, [0.1], dimension, mode, options,
    )
    zero_obs = replace(obs, pair_differences=(np.zeros(3), np.zeros(3)))
    _, null_cuts, null_diagnostic = calculate_from_products(
        zero_obs, amplitude, phase, [0.1], dimension, mode, options,
    )
    np.testing.assert_array_equal(cuts, null_cuts)
    np.testing.assert_array_equal(diagnostic["log_likelihood"], null_diagnostic["log_likelihood"])
    assert diagnostic == null_diagnostic


def test_observed_option_is_removed():
    with pytest.raises(TypeError, match="observed"):
        LimitOptions(observed=True)


def test_std_defaults_to_diagonal_of_same_cor_covariance():
    amplitude, phase = products()
    obs = observations()
    options = LimitOptions(phase_count=6)
    _, cuts, diag = calculate_from_products(obs, amplitude, phase, [0.1], options=options)
    diagonal_obs = replace(obs, pair_covariances=tuple(np.diag(np.diag(c))
                                                     for c in obs.pair_covariances))
    _, corcuts, cor = calculate_from_products(diagonal_obs, amplitude, phase, [0.1],
                                            mode="cor", options=options)
    np.testing.assert_array_equal(cuts, corcuts)
    np.testing.assert_allclose(diag["log_likelihood"], cor["log_likelihood"], atol=1e-14)
    _, oldcuts, _ = calculate_from_products(obs, amplitude, phase, [0.1],
                                           options=replace(options, std_covariance="legacy-sigma"))
    assert not np.allclose(cuts, oldcuts)


@pytest.mark.parametrize("mode", ["std", "cor", "cor_grmhd"])
def test_compatibility_loader_keeps_legacy_sigma_and_return_pair(monkeypatch, mode):
    import eht_axion.limits as limits

    obs = observations()
    amplitude, phase = products()
    monkeypatch.setattr(limits, "load_1d_observations", lambda *args: obs)
    monkeypatch.setattr(limits, "load_fit_parameters", lambda *args: (amplitude, phase))
    options = LimitOptions(legacy=True, coupling_points=11)
    alphas, cuts = calculate_limits(Path("unused"), dimension="1d", mode=mode, options=options)
    expected_alphas, expected_cuts, diagnostic = calculate_from_products(
        obs, amplitude, phase, ALPHAS, "1d", mode, replace(options, std_covariance="legacy-sigma"),
    )
    np.testing.assert_array_equal(alphas, expected_alphas)
    np.testing.assert_array_equal(cuts, expected_cuts)
    assert diagnostic["phase_count"] == 16


def test_legacy_templates_and_posterior_strict_notebook_regression():
    obs = observations()
    amplitude, phase = products()
    # The historical oracle requires explicit N=4, never the corrected default.
    options = LimitOptions(legacy=True, std_covariance="legacy-sigma", n_principal=4)
    _, cuts, diagnostics = calculate_from_products(obs, amplitude, phase, [0.1], options=options)
    # Independent notebook expression, including duplicate endpoint, omega4,
    # half-speed five-day shift, likelihood addition and discrete log-grid CDF.
    frequency = 0.112411 * 0.1 * (1 - 0.1**2 / 32)
    exposure = frequency / 2 * 6
    constant = 2 * np.sin(frequency / 2 * 24) * np.sin(exposure) / exposure / 0.1
    logs = []
    for g in diagnostics["couplings"]:
        phase_likes = []
        for p in np.linspace(0, 2 * np.pi, 16):
            epoch_likes = []
            for epoch in range(2):
                template = constant * amplitude[0.1][0] * np.sin(
                    p + phase[0.1][0] + epoch * frequency / 2 * 24 * 5
                )
                residual = g * cloud_amplitude_scaling(0.1) * template[obs.selected_indices]
                epoch_likes.append(-0.5 * np.sum(residual**2 / 0.32))
            phase_likes.append(logsumexp(epoch_likes))
        logs.append(logsumexp(phase_likes) - np.log(16))
    np.testing.assert_allclose(diagnostics["log_likelihood"][0], logs, atol=2e-12, rtol=2e-14)
    expected_cut = _legacy_upper_limit(np.array(diagnostics["couplings"]),
                                       np.exp(np.array(logs) - max(logs)))
    assert cuts[0] == pytest.approx(expected_cut, abs=1e-13)


@pytest.mark.parametrize("legacy,method", [(True, None), (True, "legacy"), (False, "legacy")])
def test_legacy_discrete_cdf_rejects_uniform_prior_before_loading(monkeypatch, legacy, method):
    import eht_axion.limits as limits

    def forbidden(*args):
        pytest.fail("invalid prior must be rejected before loading products")

    monkeypatch.setattr(limits, "load_1d_observations", forbidden)
    options = LimitOptions(legacy=legacy, posterior_method=method, prior="uniform")
    amplitude, phase = products()
    with pytest.raises(ValueError, match="legacy posterior_method requires prior='log-uniform'"):
        calculate_from_products(observations(), amplitude, phase, [0.1], options=options)
    with pytest.raises(ValueError, match="legacy posterior_method requires prior='log-uniform'"):
        calculate_limits(Path("unused"), dimension="1d", mode="std", options=options)


@pytest.mark.parametrize("legacy", [False, True])
def test_explicit_trapezoid_uniform_prior_is_valid_and_metadata_is_truthful(legacy):
    amplitude, phase = products(zero_signal=True)
    options = LimitOptions(legacy=legacy, posterior_method="trapezoid", prior="uniform",
                           coupling_min=1, coupling_max=101, coupling_points=9)
    _, cuts, diagnostic = calculate_from_products(observations(), amplitude, phase, [0.1],
                                                 options=options)
    assert cuts[0] == pytest.approx(np.log10(96.0), abs=1e-14)
    assert diagnostic["prior"] == "uniform"
    assert diagnostic["posterior_method"] == "trapezoid"
    assert diagnostic["density_measure"] == "dg"


@pytest.mark.parametrize("legacy", [False, True])
def test_explicit_legacy_cdf_keeps_equal_log_grid_result_and_metadata(legacy):
    amplitude, phase = products(zero_signal=True)
    options = LimitOptions(legacy=legacy, posterior_method="legacy", coupling_points=9)
    _, cuts, diagnostic = calculate_from_products(observations(), amplitude, phase, [0.1],
                                                 options=options)
    grid = np.asarray(diagnostic["couplings"])
    assert cuts[0] == pytest.approx(_legacy_upper_limit(grid, np.ones(9)), abs=1e-14)
    np.testing.assert_allclose(diagnostic["cdf"][0], np.arange(1, 10) / 9)
    assert diagnostic["prior"] == "log-uniform"
    assert diagnostic["posterior_method"] == "legacy"
    assert diagnostic["density_measure"] == "discrete-log-grid-mass"


@pytest.mark.parametrize("change", [
    {"coupling_points": 1}, {"coupling_points": 2.5}, {"coupling_min": 0},
    {"coupling_max": np.inf}, {"grmhd_sigma_degrees": np.nan},
    {"grmhd_sigma_degrees": -1}, {"posterior_method": "oops"},
    {"std_covariance": "oops"}, {"prior": "oops"},
])
def test_option_guards(change):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_from_products(observations(), amplitude, phase, [0.1],
                                options=replace(LimitOptions(), **change))


@pytest.mark.parametrize("change", [
    {"selected_indices": np.array([])}, {"selected_indices": np.array([0.5, 1.5])},
    {"selected_indices": np.array([0, 0, 1])}, {"selected_indices": np.array([0, 1, 4])},
    {"pair_differences": (np.array([np.nan, 0, 0]), np.zeros(3))},
    {"pair_differences": (np.zeros(2), np.zeros(2))},
    {"pair_differences": (np.zeros(3),)},
    {"pair_covariances": (np.eye(2), np.eye(2))},
    {"pair_covariances": (np.array([[1, 0.1, 0], [0, 1, 0], [0, 0, 1]]), np.eye(3))},
    {"pair_covariances": (np.full((3, 3), np.inf), np.eye(3))},
])
def test_observation_guards(change):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_from_products(replace(observations(), **change), amplitude, phase, [0.1])


@pytest.mark.parametrize("alphas", [[], [np.nan], [0.0], [0.2, 0.1], [0.1, 0.1]])
def test_alpha_guards(alphas):
    amplitude, phase = products()
    with pytest.raises(ValueError):
        calculate_from_products(observations(), amplitude, phase, alphas)
