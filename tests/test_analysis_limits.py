"""Independent-day working covariance semantics; no large source files."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from eht_axion.evpa import circular_difference
from eht_axion.analysis.likelihood import LimitOptions
from eht_axion.analysis.observations import (
    branch_risk, build_observations, covariance_diagnostics, moments,
)
from eht_axion.analysis.likelihood import calculate_observed


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'eht_ensemble_limits_driver', ROOT / 'scripts/analysis_limits.py')
driver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(driver)


def samples():
    residual = np.array([[-.10, -.20], [0., .10], [.10, .10], [.02, -.03]])
    residual -= residual.mean(axis=0)
    centers = np.array([[1.54, -.2], [1.64, -.1], [1.4, .1], [1.2, .4]])
    return [c + s * residual for c, s in zip(centers, (1, 2, 1.5, .5))], centers


def test_population_scatter_not_mean_error_and_pi_invariance():
    data, _ = samples()
    mean, cov, resultant = moments(data[0])
    changed = data[0] + np.pi * np.array([[0, 1], [2, -1], [1, 0], [-3, 2]])
    m2, c2, r2 = moments(changed)
    np.testing.assert_allclose(circular_difference(mean, m2), 0, atol=1e-14)
    np.testing.assert_allclose(cov, np.cov(data[0], rowvar=False, ddof=0), atol=1e-14)
    np.testing.assert_allclose(c2, cov, atol=1e-14)
    np.testing.assert_allclose(r2, resultant)
    np.testing.assert_allclose(moments(np.repeat(data[0], 3, axis=0))[1], cov, atol=1e-14)


def test_independent_day_pair_sum_and_zero_data_control():
    data, centers = samples()
    obs, null, metadata = build_observations(data, np.array([4, 14]), resultant_min=.7)
    circular_means = np.array([.5 * np.angle(np.mean(np.exp(2j * v), axis=0)) for v in data])
    np.testing.assert_allclose(obs.means, circular_means, atol=1e-14)
    assert np.max(np.abs(circular_difference(obs.means, centers))) > .001
    np.testing.assert_allclose(obs.pair_differences,
        circular_difference(circular_means[[1, 3]], circular_means[[0, 2]]), atol=1e-14)
    offsets = np.array(metadata['local_branch_mean_minus_circular_radians'])
    np.testing.assert_allclose(offsets, circular_difference(centers, circular_means), atol=1e-14)
    for k, (a, b) in enumerate(((0, 1), (2, 3))):
        np.testing.assert_allclose(obs.pair_covariances[k], obs.covariances[a] + obs.covariances[b])
    assert null.pair_covariances is obs.pair_covariances
    assert null.covariances is obs.covariances
    np.testing.assert_array_equal(null.pair_differences, np.zeros((2, 2)))
    assert metadata['covariance_of_mean'] is False
    # Arbitrary independent draw permutations cannot change the answer.
    permuted = [v[np.random.default_rng(i).permutation(len(v))] for i, v in enumerate(data)]
    other, _, _ = build_observations(permuted, np.array([4, 14]), resultant_min=.7)
    np.testing.assert_allclose(other.pair_covariances, obs.pair_covariances, atol=1e-14)


@pytest.mark.parametrize('value', [True, np.nan, np.inf, 0, 1.1])
def test_invalid_concentration(value):
    data, _ = samples()
    with pytest.raises(ValueError, match='resultant_min'):
        build_observations(data, np.array([4, 14]), resultant_min=value)


@pytest.mark.parametrize('bad', [np.ones((1, 2)), np.ones((4, 0)), [[True], [False]],
                                 [[np.nan], [0]], [[1j], [2j]]])
def test_invalid_samples(bad):
    with pytest.raises(ValueError, match='EVPA samples'):
        moments(bad)


def test_branch_risk_reports_seams_without_selection():
    differences = (np.array([np.pi / 2 - .01, .1]), np.array([-.2, .3]))
    covariances = (np.diag([.02**2, .1**2]), np.diag([.03**2, .04**2]))
    risk = branch_risk(differences, covariances)
    assert risk[0]['count_within_one_sigma'] == 1
    assert risk[0]['count_within_three_sigma'] == 1
    assert risk[1]['count_within_three_sigma'] == 0
    assert risk[0]['distance_to_seam_in_sigma'][0] == pytest.approx(.5)
    assert risk[0]['max_abs_difference_degrees'] == pytest.approx(90 - np.rad2deg(.01))
    stats = covariance_diagnostics(covariances)
    assert stats[0]['condition_number'] == pytest.approx(25)
    assert stats[0]['eigenvalue_min'] == .02**2


def test_matched_aperture_manual_intensity_quadrature():
    config = driver.read_experiment()
    operator = driver.theory_aperture_operator(config)
    row, col = np.indices((64, 64))
    frame = np.zeros((64, 64, 4))
    # East-localized I and spatially varying Q/U expose wrong axes or weighting.
    frame[..., 0] = 1 + 5 * np.exp(-((row - 32)**2 + (col - 17)**2) / 20)
    frame[..., 1] = np.cos(.02 * row + .03 * col)
    frame[..., 2] = np.sin(.02 * row + .03 * col)
    actual = driver.matched_aperture_samples(frame[None], config, operator)[0]
    template = frame.transpose(1, 0, 2)
    qi = operator @ (template[..., 0] * template[..., 1]).ravel()
    ui = operator @ (template[..., 0] * template[..., 2]).ravel()
    manual = np.roll(.5 * np.arctan2(ui, qi), -90)
    np.testing.assert_allclose(actual, manual, atol=1e-14)
    flux_only = np.roll(.5 * np.arctan2(operator @ template[..., 2].ravel(),
                                      operator @ template[..., 1].ravel()), -90)
    assert np.max(np.abs(actual - flux_only)) > .001


def test_no_hidden_covariance_repair():
    with pytest.raises(ValueError, match='positive definite'):
        build_observations([np.zeros((10, 2))] * 4, np.array([4, 14]), resultant_min=.7)


def test_three_fixed_c_modes_share_grid_and_phase_product():
    data, _ = samples()
    obs, null, _ = build_observations(data, np.array([0, 1]), resultant_min=.7)
    amp = {.1: np.array([[1., .8], [0., 0.]])}
    phase = {.1: np.array([[.1, .4], [0., 0.]])}
    options = LimitOptions(prior='log-uniform', credibility=.95, coupling_points=12, coupling_min=1 / 137,
                           coupling_max=1000.,
                           grmhd_sigma_degrees=2., phase_count=8, time_conversion=0.112411,
                           exposure_hours=6., day_gap_hours=24., pair_separation_hours=120., n_principal=2,
                           reference_norm=.1, cloud_coefficient=.03, cloud_pivot_alpha=.05,
                           cloud_exponent=1.75)
    for mode in ('std', 'cor', 'cor_fixed_background'):
        for center in (obs, null):
            _, cuts, info = calculate_observed(center, amp, phase, alphas=[.1], mode=mode,
                                                options=options)
            assert np.isfinite(cuts).all()
            assert info['response_model_variance_included'] is False
            assert info['phase_count'] == 8
            assert info['credibility'] == .95
            assert info['data_nonzero_truth'] == (center is obs)
            if mode == 'cor_fixed_background':
                expected = (np.array(obs.pair_covariances)
                            + np.eye(2) * np.deg2rad(options.grmhd_sigma_degrees)**2)
                np.testing.assert_allclose(info['effective_pair_covariances'], expected)


def test_native_east_feature_transpose_and_no_second_aperture_roll(tmp_path, monkeypatch):
    import json

    config = driver.read_experiment()
    ensemble, theory = tmp_path / 'ensemble', tmp_path / 'theory'
    ensemble.mkdir()
    grid = (np.arange(64) - 31.5) * 1.875
    np.savez(ensemble / 'axes.npz', east_uas=-grid, north_uas=grid, angle_deg=np.arange(360))
    arrays = {}
    row, col = np.indices((64, 64))
    # Nonuniform east/north feature: omitted transpose produces a different angle.
    native_chi = .002 * row + .004 * col
    for day in driver.DAYS:
        maps = np.ones((64, 64, 4))
        maps[..., 1] = np.cos(2 * native_chi)
        maps[..., 2] = np.sin(2 * native_chi)
        directory = ensemble / f'day-{day}'
        directory.mkdir(parents=True)
        arrays[directory / 'stokes.npy'] = np.broadcast_to(maps, (1000, 64, 64, 4))
        arrays[directory / 'evpa.npy'] = np.broadcast_to(native_chi, (1000, 64, 64))
        aperture = np.broadcast_to(np.arange(360) * .001, (1000, 360))
        np.savez(directory / 'aperture.npz', samples_rad=aperture, angle_deg=np.arange(360),
                 intensity=np.ones_like(aperture), polarized_brightness=np.ones_like(aperture))
    report = dict(schema='eht-ensemble-v2', days=[str(d) for d in driver.DAYS],
        draws_per_day=1000, rmax_uas=60., aperture_full_width_deg=10.,
        angle_convention='EVPA radians modulo pi, E of N; aperture angles degrees E of N; no roll',
        map_order='draw,y,x,IQUV; native FITS E decreases along x, N increases along y')
    (ensemble / 'report.json').write_text(json.dumps(report))
    records = {ensemble: dict(stage='eht:ensemble', config=dict(draw_count=4000))}
    original_load = np.load
    monkeypatch.setattr(driver.np, 'load',
                        lambda p, **kw: arrays[p] if p in arrays else original_load(p, **kw))
    processed = theory / 'background'
    processed.mkdir(parents=True)
    np.save(processed / 'stokes.npy', np.ones((5, 1, 64, 64, 4)))
    geometry = driver.analysis_geometry(config)
    delivered = tmp_path / 'raptor/res-400/background/t1600-t1690/stokes-t1600-t1690.npy'
    records[processed] = dict(stage='analysis:background',
                             config=dict(geometry['image_preprocess'], state_indices=[0]),
                             inputs=[dict(path=str(delivered)), dict(path=str(delivered.parent / 'x.json'))])
    checked, cadences = [], []
    monkeypatch.setattr(driver, 'validate_processed_source',
                        lambda record, config, **kw: checked.append((record['stage'], kw['alpha'],
                                                                     kw['resolution'])))
    monkeypatch.setattr(driver, 'input_cadence_hours',
                        lambda paths, config: cadences.append(paths) or 8.)
    for path in records:
        (path / 'manifest.json').write_text('{}')
    monkeypatch.setattr(driver, 'verify_stage', lambda path, **kw: records[path])
    monkeypatch.setattr(driver, 'require_declared_outputs', lambda *args: None)
    valid2d = np.zeros((64, 64), dtype=bool)
    valid2d[20, 35] = True
    products = {'1d': (None, None, None, np.ones(360, dtype=bool)),
                '2d': (None, None, None, valid2d)}
    prepared, _, _, sigmas, background_inputs = driver.prepare_limit_inputs(
        ensemble, theory, config, products, verified=set())
    assert checked == [('analysis:background', None, 400)]
    assert cadences == [[delivered]]
    assert sigmas == {dim: dict(stokes=str(processed / 'stokes.npy'), cadence_hours=8.)
                      for dim in ('1d', '2d')}
    assert delivered.parent / 'generated.json' in background_inputs
    selected, daily, coords, _ = prepared['2d']
    np.testing.assert_array_equal(selected, [20 * 64 + 35])
    np.testing.assert_allclose(daily[0], .002 * 35 + .004 * 20)
    np.testing.assert_allclose(coords, [[-grid[20], grid[35]]])
    from eht_axion.analysis.templates import _integrate
    expected = _integrate(maps.transpose(1, 0, 2), driver.theory_aperture_operator(config))
    indices = prepared['1d'][0]
    np.testing.assert_allclose(prepared['1d'][1][0][0], expected[indices])
    assert not np.allclose(expected[indices], indices * .001)  # native sectors were deliberately different
    wrong = _integrate(maps, driver.theory_aperture_operator(config))
    assert not np.allclose(expected[indices], wrong[indices])
    np.savez(ensemble / 'axes.npz', east_uas=-grid, north_uas=grid,
             angle_deg=np.roll(np.arange(360), -90))
    with pytest.raises(ValueError, match='no extra roll'):
        driver.prepare_limit_inputs(ensemble, theory, config, products, verified=set())


def test_options_bind_actual_physical_scale_and_new_background():
    config = driver.read_experiment()
    options = driver.options_from_config(config, .031)
    assert options.phase_count == 64
    assert options.coupling_points == 120
    assert options.coupling_min == 1 / 137
    assert options.coupling_max == 1000
    assert options.prior == 'log-uniform'
    assert options.n_principal == 2
    assert options.reference_norm == .1
    assert options.grmhd_sigma_degrees == np.rad2deg(.031)
    assert options.time_conversion == 1 / driver.physical_scale(config).tg_hours
    for sigma in (True, np.nan, -1):
        with pytest.raises(ValueError, match='Background sigma'):
            driver.options_from_config(config, sigma)


@pytest.mark.parametrize('dimension', ['1d', '2d'])
def test_selected_background_exact_feature_domain(tmp_path, monkeypatch, dimension):
    import json

    config = driver.read_experiment()
    values = np.array([[0., 0.], [.01, .02], [.02, .03], [.05, .07], [.08, .12]])
    chi = np.zeros((5, 64, 64))
    chi[:, 20, 35], chi[:, 30, 40] = values.T
    maps = np.ones((5, 1, 64, 64, 4))
    maps[:, 0, ..., 1], maps[:, 0, ..., 2] = np.cos(2 * chi), np.sin(2 * chi)
    path = tmp_path / 'background.npy'
    np.save(path, maps)
    selected = np.array([20 * 64 + 35, 30 * 64 + 40])
    if dimension == '1d':
        selected = np.array([4, 14])
        import eht_axion.analysis.templates as pol
        monkeypatch.setattr(driver, 'theory_aperture_operator', lambda config: None)
        def integrated(frame, operator):
            out = np.zeros(360)
            out[94] = .5 * np.arctan2(frame[20, 35, 2], frame[20, 35, 1])
            out[104] = .5 * np.arctan2(frame[30, 40, 2], frame[30, 40, 1])
            return out
        monkeypatch.setattr(pol, 'integrate_stokes', integrated)
    sigma = driver.selected_background(tmp_path, dimension, selected, config,
        dict(stokes=str(path), cadence_hours=8.))
    expected = values[3:] - values[:-3]
    assert sigma == pytest.approx(expected.std(ddof=0))
    with np.load(tmp_path / f'background-selected-{dimension}.npz') as z:
        np.testing.assert_allclose(z['differences'], expected, atol=1e-14)
    meta = json.loads((tmp_path / f'background-selected-{dimension}.json').read_text())
    assert meta['lag'] == 3 and meta['lag_hours'] == 24.
    assert 'source_unmasked_sigma' not in meta
    assert meta['sqrt2_correction'] is False
    assert meta['lag_hours'] == 24


def test_save_all_six_cases_and_zero_data_control(tmp_path, monkeypatch):
    config = driver.read_experiment()
    data, _ = samples()
    alphas = np.array(config['shared']['alphas'])
    products, prepared = {}, {}
    for dim in ('1d', '2d'):
        products[dim] = (alphas, {}, {}, None)
        prepared[dim] = (np.array([4, 14]), data, np.array([[1., 2.], [3., 4.]]), {'synthetic': True})
    calls = []

    def calculate(obs, amp, phase, *, alphas, dimension, mode, options):
        calls.append((dimension, mode, obs, options))
        return alphas, np.zeros(8), dict(couplings=[.01, 1.], log_likelihood=[[0., -1.]] * 8,
            posterior_density=[[1., .1]] * 8, cdf=[[0., 1.]] * 8)

    monkeypatch.setattr(driver, 'calculate_observed', calculate)
    monkeypatch.setattr(driver, 'plot_results', lambda *args: None)
    monkeypatch.setattr(driver, 'selected_background', lambda out, dim, selected, config, source: source['s'])
    report = driver.save_results(tmp_path, products, prepared, config,
                                  {'1d': dict(s=.02), '2d': dict(s=.03)}, {'resultant_min': .7})
    assert len(calls) == 12
    assert len(report['supported_cases']) == 6
    assert len(report['results']) == 12
    assert report['response_model_variance_included'] is False
    assert len((tmp_path / 'limits.csv').read_text().splitlines()) == 97
    for i in range(0, 12, 2):
        dim, mode, obs, options = calls[i]
        null_dim, null_mode, null, null_options = calls[i + 1]
        assert (dim, mode) == (null_dim, null_mode)
        assert null.pair_covariances is obs.pair_covariances
        assert options is null_options
        np.testing.assert_array_equal(null.pair_differences, np.zeros((2, 2)))
        assert np.any(np.array(obs.pair_differences) != 0)
        assert options.grmhd_sigma_degrees == np.rad2deg(.02 if dim == '1d' else .03)
