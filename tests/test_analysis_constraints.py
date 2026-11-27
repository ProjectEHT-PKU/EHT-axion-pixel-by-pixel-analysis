"""Independent-day working covariance semantics; no large source files."""
import numpy as np
import pytest

from eht_axion.evpa import circular_difference
from eht_axion.analysis.likelihood import LimitOptions
from eht_axion.analysis.observations import (
    branch_risk, build_observations, covariance_diagnostics, moments,
)
from eht_axion.analysis.likelihood import calculate_observed
from scripts.analysis import constraints as driver
from eht_axion.analysis.results import read_results
import json


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
    obs, zero_data, metadata = build_observations(data, np.array([4, 14]), resultant_min=.7)
    circular_means = np.array([.5 * np.angle(np.mean(np.exp(2j * v), axis=0)) for v in data])
    np.testing.assert_allclose(obs.means, circular_means, atol=1e-14)
    assert np.max(np.abs(circular_difference(obs.means, centers))) > .001
    np.testing.assert_allclose(obs.pair_differences,
        circular_difference(circular_means[[1, 3]], circular_means[[0, 2]]), atol=1e-14)
    offsets = np.array(metadata['local_branch_mean_minus_circular_radians'])
    np.testing.assert_allclose(offsets, circular_difference(centers, circular_means), atol=1e-14)
    for k, (a, b) in enumerate(((0, 1), (2, 3))):
        np.testing.assert_allclose(obs.pair_covariances[k], obs.covariances[a] + obs.covariances[b])
    assert zero_data.pair_covariances is obs.pair_covariances
    assert zero_data.covariances is obs.covariances
    np.testing.assert_array_equal(zero_data.pair_differences, np.zeros((2, 2)))
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
    evpa, intensity, polarized = driver.matched_aperture_samples(frame[None], operator)
    actual = evpa[0]
    template = frame.transpose(1, 0, 2)
    qi = operator @ (template[..., 0] * template[..., 1]).ravel()
    ui = operator @ (template[..., 0] * template[..., 2]).ravel()
    manual = np.roll(.5 * np.arctan2(ui, qi), -90)
    np.testing.assert_allclose(actual, manual, atol=1e-14)
    # Selection statistics: the same weights W and the same intensity weighting.
    sector_i = operator @ template[..., 0].ravel()
    sector_p = operator @ (template[..., 0] * np.hypot(template[..., 1], template[..., 2])).ravel()
    np.testing.assert_allclose(intensity[0], np.roll(sector_i, -90), rtol=1e-13)
    np.testing.assert_allclose(polarized[0], np.roll(sector_p / sector_i, -90), rtol=1e-13)
    flux_only = np.roll(.5 * np.arctan2(operator @ template[..., 2].ravel(),
                                      operator @ template[..., 1].ravel()), -90)
    assert np.max(np.abs(actual - flux_only)) > .001


def test_no_hidden_covariance_repair():
    with pytest.raises(ValueError, match='positive definite'):
        build_observations([np.zeros((10, 2))] * 4, np.array([4, 14]), resultant_min=.7)


def test_three_fixed_c_modes_share_grid_and_phase_product():
    data, _ = samples()
    obs, zero_data, _ = build_observations(data, np.array([0, 1]), resultant_min=.7)
    amp = {.1: np.array([[1., .8], [0., 0.]])}
    phase = {.1: np.array([[.1, .4], [0., 0.]])}
    options = LimitOptions(prior='log-uniform', credibility=.95, coupling_points=12, coupling_min=1 / 137,
                           coupling_max=1000.,
                           grmhd_sigma_degrees=2., phase_count=8, time_conversion=0.112411,
                           exposure_hours=6., day_gap_hours=24., pair_separation_hours=120., n_principal=2,
                           reference_norm=.1, cloud_coefficient=.03, cloud_pivot_alpha=.05,
                           cloud_exponent=1.75)
    for mode in ('std', 'cor', 'cor_fixed_background'):
        for observation in (obs, zero_data):
            _, cuts, info = calculate_observed(observation, amp, phase, alphas=[.1], mode=mode,
                                                options=options)
            assert np.isfinite(cuts).all()
            assert info['response_model_variance_included'] is False
            assert info['phase_count'] == 8
            assert info['credibility'] == .95
            assert info['data_nonzero_truth'] == (observation is obs)
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
    np.savez(ensemble / 'axes.npz', east_uas=-grid, north_uas=grid)
    arrays = {}
    row, col = np.indices((64, 64))
    # Nonuniform east/north feature: omitted transpose produces a different angle.
    native_chi = .002 * row + .004 * col
    maps = np.ones((64, 64, 4))
    maps[..., 1] = np.cos(2 * native_chi)
    maps[..., 2] = np.sin(2 * native_chi)
    arrays[ensemble / 'stokes.npy'] = np.broadcast_to(maps, (len(driver.DAYS), 1000, 64, 64, 4))
    report = dict(schema='eht-ensemble', days=[str(d) for d in driver.DAYS], draws_per_day=1000,
        angle_convention='EVPA radians modulo pi, east of north',
        map_order='day,draw,y,x,IQUV; native FITS E decreases along x, N increases along y')
    (ensemble / 'report.json').write_text(json.dumps(report))
    records = {ensemble: dict(stage='eht:ensemble', config=dict(draw_count=4000))}
    original_load = np.load
    monkeypatch.setattr(driver.np, 'load',
                        lambda p, **kw: arrays[p] if p in arrays else original_load(p, **kw))
    theory.mkdir(parents=True)
    np.save(theory / 'background.npy', np.ones((5, 1, 64, 64, 4)))
    delivered = tmp_path / 'raptor/images/res-400/background/stokes.npy'
    records[theory] = dict(stage='analysis:expected-signal',
                           inputs=[str(tmp_path / 'raptor/images/res-100/alpha-0.10/stokes.npy'),
                                   str(delivered), str(delivered.parent / 'generated.json')])
    cadences = []
    monkeypatch.setattr(driver, 'input_cadence_hours',
                        lambda path, config: cadences.append(path) or 8.)
    for path in records:
        (path / 'manifest.json').write_text('{}')
    monkeypatch.setattr(driver, 'verify_stage', lambda path, **kw: records[path])
    monkeypatch.setattr(driver, 'require_declared_outputs', lambda *args: None)
    valid2d = np.zeros((64, 64), dtype=bool)
    valid2d[20, 35] = True
    products = {'1d': (None, None, None, np.ones(360, dtype=bool)),
                '2d': (None, None, None, valid2d)}
    prepared, _, _, sigmas, background_inputs = driver.prepare_constraint_inputs(
        ensemble, theory, config, products)
    assert cadences == [delivered]
    assert sigmas == {dim: dict(stokes=str(theory / 'background.npy'), cadence_hours=8.)
                      for dim in ('1d', '2d')}
    assert delivered.parent / 'generated.json' in background_inputs
    selected, daily, coords, _ = prepared['2d']
    np.testing.assert_array_equal(selected, [20 * 64 + 35])
    np.testing.assert_allclose(daily[0], .002 * 35 + .004 * 20)
    np.testing.assert_allclose(coords, [[-grid[20], grid[35]]])
    from eht_axion.analysis.expected_signal import _integrate
    expected = _integrate(maps.transpose(1, 0, 2), driver.theory_aperture_operator(config))
    indices = prepared['1d'][0]
    # Uniform I and |P|: every sector passes both brightness thresholds; only the sector rule acts.
    np.testing.assert_array_equal(indices, np.arange(4, 360, 10))
    np.testing.assert_allclose(prepared['1d'][1][0][0], expected[indices])
    wrong = _integrate(maps, driver.theory_aperture_operator(config))
    assert not np.allclose(expected[indices], wrong[indices])


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
        import eht_axion.analysis.expected_signal as pol
        monkeypatch.setattr(driver, 'theory_aperture_operator', lambda config: None)
        def integrated(frame, operator):
            out = np.zeros(360)
            out[94] = .5 * np.arctan2(frame[20, 35, 2], frame[20, 35, 1])
            out[104] = .5 * np.arctan2(frame[30, 40, 2], frame[30, 40, 1])
            return out
        monkeypatch.setattr(pol, 'integrate_stokes', integrated)
    sigma, meta = driver.selected_background(dimension, selected, config,
        dict(stokes=str(path), cadence_hours=8.))
    expected = values[3:] - values[:-3]
    assert sigma == pytest.approx(expected.std(ddof=0)) and meta['sigma_radians'] == sigma
    np.testing.assert_allclose(meta['differences'], expected, atol=1e-14)
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
    monkeypatch.setattr(driver, 'selected_background',
                        lambda dim, selected, config, source: (source['s'], dict(sigma_radians=source['s'])))
    report = driver.save_results(tmp_path, products, prepared, config,
                                  {'1d': dict(s=.02), '2d': dict(s=.03)}, {'resultant_min': .7})
    assert len(calls) == 12
    assert len(report['supported_cases']) == 6
    assert len(report['results']) == 12
    assert report['response_model_variance_included'] is False
    assert len((tmp_path / 'constraints.csv').read_text().splitlines()) == 97
    assert {p.name for p in tmp_path.iterdir()} == {'constraints.csv', 'report.json', 'results.npz'}
    stored, records = read_results(tmp_path, ('observations', 'background', 'cases'))
    assert stored['results'] == json.loads(json.dumps(report['results']))
    assert set(records['cases']) == set(report['results'])
    assert records['background']['2d'] == dict(sigma_radians=.03)
    case = records['cases']['observed-cor-1d']
    np.testing.assert_array_equal(case['couplings'], [.01, 1.])
    assert case['log_likelihood'].shape == (8, 2) and case['log10_c95'].shape == (8,)
    assert records['observations']['1d']['selection'] == {'synthetic': True}
    assert records['observations']['1d']['evpa_samples_radians'].shape == (4, 4, 2)
    for i in range(0, 12, 2):
        dim, mode, obs, options = calls[i]
        zero_dim, zero_mode, zero_data, zero_options = calls[i + 1]
        assert (dim, mode) == (zero_dim, zero_mode)
        assert zero_data.pair_covariances is obs.pair_covariances
        assert options is zero_options
        np.testing.assert_array_equal(zero_data.pair_differences, np.zeros((2, 2)))
        assert np.any(np.array(obs.pair_differences) != 0)
        assert options.grmhd_sigma_degrees == np.rad2deg(.02 if dim == '1d' else .03)
