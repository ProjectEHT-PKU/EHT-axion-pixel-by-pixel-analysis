"""Entry-point contracts with synthetic data only; never execute science imaging."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('reconstruct_eht_official_test',
                                               ROOT / 'scripts/reconstruct_eht_official.py')
w = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(w)


@pytest.fixture
def config():
    return yaml.safe_load((ROOT / 'configs/eht-official.yaml').read_text())


@pytest.fixture
def plan_args(tmp_path, config, monkeypatch):
    """Fully local synthetic receipts; no GRMHD inputs and no optional package import."""
    raw = tmp_path / config['inputs'][95][0]
    raw.write_bytes(b'synthetic circular UVFITS placeholder')
    config['inputs'][95][1] = w.common.sha256(raw)
    receipt = tmp_path / 'receipt.json'
    receipt.write_text(json.dumps(dict(complete=True, files=[dict(path=str(raw),
        size=raw.stat().st_size, sha256=config['inputs'][95][1])])) )
    feed, official = tmp_path / 'feed.json', tmp_path / 'official.py'
    feed.write_text('{}')
    official.write_text("ttype     = 'nfft'\nfov       = 128*eh.RADPERUAS\n")
    config.update(input_root=str(tmp_path), input_receipt=str(receipt), feed_rotation_file=str(feed),
                  official_script=str(official))
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(config))
    dist = SimpleNamespace(version='1.3.2', files=[], locate_file=lambda f: f)
    monkeypatch.setattr(w.metadata, 'distribution', lambda name: dist)
    return w.parser().parse_args(['--config', str(path), '--stage', 'I', '--day', '95',
                                  '--output-root', str(tmp_path / 'new'), '--plan'])


def test_plan_no_optional_import_no_writes(plan_args, monkeypatch):
    monkeypatch.setitem(sys.modules, 'ehtim', None)
    monkeypatch.setitem(sys.modules, 'astropy', None)
    plan = w.build_plan(plan_args)
    assert plan['day'] == 95 and plan['stage'] == 'I'
    assert plan['scientifically_converged'] is False
    assert not plan_args.output_root.exists()
    paths = {Path(r['path']) for r in plan['inputs']}
    assert ROOT / 'scripts/reconstruct_eht_official.py' in paths
    assert ROOT / 'eht_axion/pipeline/common.py' in paths
    assert plan_args.config in paths
    assert all(len(r['sha256']) == 64 for r in plan['inputs'])


def test_cli_review_stage_draw_contract(plan_args):
    plan_args.plan = False
    with pytest.raises(ValueError, match='review-file'):
        w.build_plan(plan_args)
    plan_args.plan = True
    plan_args.stage = 'fiducial'
    with pytest.raises(ValueError, match='i-root'):
        w.build_plan(plan_args)
    plan_args.i_root = Path('not-read-yet')
    plan_args.stage = 'uncertainty'
    with pytest.raises(ValueError, match='D-roots|d-roots'):
        w.build_plan(plan_args)
    plan_args.stage, plan_args.i_root, plan_args.draws = 'I', None, 25
    with pytest.raises(ValueError, match='only to uncertainty'):
        w.build_plan(plan_args)


def test_wrong_hash_existing_destination_and_output_source_protection(plan_args):
    plan_args.output_root.mkdir()
    with pytest.raises(FileExistsError, match='new stage/day'):
        w.build_plan(plan_args)
    plan_args.output_root.rmdir()
    plan_args.output_root = plan_args.config.parent
    with pytest.raises(FileExistsError):
        w.build_plan(plan_args)
    plan_args.output_root = plan_args.config.parent / 'different'
    raw = plan_args.config.parent / yaml.safe_load(plan_args.config.read_text())['inputs'][95][0]
    raw.write_bytes(b'tampered')
    with pytest.raises(ValueError, match='checksum'):
        w.build_plan(plan_args)


@pytest.mark.parametrize('text', ['', "ttype = 'nfft'\nfov       = 128*eh.RADPERUAS",
    "ttype     = 'nfft'\n" * 2 + 'fov       = 128*eh.RADPERUAS',
    "ttype     = 'nfft'\nfov = 128*eh.RADPERUAS"])
def test_adapt_requires_unique_exact_anchors(text):
    with pytest.raises(ValueError, match='anchor'):
        w.adapt_official(text)


def test_adaptation_changes_only_two_exact_strings():
    text = ("# nfft 128 remains here\nttype     = 'nfft' # retained\n"
            'fov       = 128*eh.RADPERUAS # retained\nrest = 5\n')
    assert w.adapt_official(text) == text.replace("ttype     = 'nfft'", "ttype     = 'direct'").replace(
        'fov       = 128*eh.RADPERUAS', 'fov       = 120*eh.RADPERUAS')


class Image:
    def __init__(self):
        self.imvec = np.array([.1, .2, .15, .15])
        self.qvec = self.imvec * .1
        self.uvec = np.zeros(4)
        self.vvec = np.zeros(4)
        self.xdim = self.ydim = 2
        self.rf = 227070703125.

    def copy(self):
        return deepcopy(self)

    def add_qu(self, q, u):
        self.qvec, self.uvec = q.ravel().copy(), u.ravel().copy()

    def add_pol_image(self, v, pol):
        assert pol == 'V'
        self.vvec = v.ravel().copy()

    def total_flux(self):
        return self.imvec.sum()

    def fovx(self):
        return 120.

    def blur_circ(self, *args):
        result = self.copy()
        for field in ('imvec', 'qvec', 'uvec', 'vvec'):
            a = getattr(self, field)
            setattr(result, field, .5 * a + .5 * a.mean())
        return result

    def save_fits(self, path):
        Path(path).write_text('synthetic image')


class Observation:
    def __init__(self):
        self.polrep = 'circ'
        self.rf = 227070703125.
        self.data = np.array([('AZ', 'PV', .3+1j, .03, .09, 2e8, 0.),
                              ('LM', 'PV', .4+1j, .07, .12, 3e8, 0.),
                              ('AZ', 'LM', .9+1j, .02, .08, 1e7, 0.)],
            dtype=[('t1', 'U2'), ('t2', 'U2'), ('rlvis', complex), ('rlsigma', float), ('lrsigma', float),
                   ('u', float), ('v', float)])
        self.tarr = np.array([('AZ', 1.), ('LM', -1.), ('PV', -1.)],
                             dtype=[('site', 'U2'), ('fr_elev', float)])
        self.tkey = {s: i for i, s in enumerate(self.tarr['site'])}

    def copy(self):
        return deepcopy(self)

    def flag_uvdist(self, uv_min):
        result = self.copy()
        result.data = result.data[np.hypot(result.data['u'], result.data['v']) >= uv_min].copy()
        return result

    def switch_polrep(self, pol):
        assert pol == 'stokes'
        result = self.copy()
        result.polrep = pol
        return result

    def save_uvfits(self, path):
        Path(path).write_text('synthetic native circular')


def test_i_native_before_full_pipeline_restored_model_shared_all_station_gains(tmp_path, config, monkeypatch):
    source = tmp_path / 'official.py'
    source.write_text("ttype     = 'nfft'\nfov       = 128*eh.RADPERUAS\n")
    config['official_script'] = str(source)
    native, fixed = Observation(), Image()
    core = SimpleNamespace(load_native=Mock(return_value=(native, {'synthetic': True})))
    gains = np.array([(0., 2.+1j, 2.+1j)], dtype=[('time', float), ('rscale', complex), ('lscale', complex)])
    original_data = native.data.copy()
    cal = SimpleNamespace(tarr=native.tarr[[2, 0, 1]], data={s: gains.copy() for s in ('AZ', 'LM', 'PV')})

    def apply(obs, **kwargs):
        assert obs is not native
        assert obs.polrep == 'circ'
        np.testing.assert_array_equal(obs.data, original_data[:2])
        assert np.all(np.hypot(obs.data['u'], obs.data['v']) >= config['uv_min'])
        np.testing.assert_array_equal(obs.tarr, cal.tarr)
        assert obs.tkey == {'PV': 0, 'AZ': 1, 'LM': 2}
        assert kwargs == dict(interp='nearest', extrapolate=True)
        return obs

    cal.applycal = Mock(side_effect=apply)
    eh = SimpleNamespace(selfcal=Mock(return_value=cal), RADPERUAS=1.)

    def run(path, run_name):
        core.load_native.assert_called_once()  # native captured BEFORE I conversion
        assert run_name == '__main__'
        assert '120*eh.RADPERUAS' in Path(path).read_text()
        assert sys.argv[-1] == str(tmp_path / 'final-fixedI.fits')
        Path(sys.argv[-1]).write_text('restored final image')
        return dict(im_out=fixed,
                    imgr=SimpleNamespace(out_last=Mock(side_effect=AssertionError('not raw I'))))

    monkeypatch.setattr(w.runpy, 'run_path', run)
    monkeypatch.setattr(w, 'validate_fixed', Mock())
    argv = sys.argv.copy()
    w.run_i(dict(parameters=config, raw='synthetic.uvfits'), core, eh, tmp_path)
    assert sys.argv == argv
    (fitcopy, model), kw = eh.selfcal.call_args
    assert fitcopy.polrep == 'stokes' and fitcopy is not native
    np.testing.assert_array_equal(fitcopy.data, original_data[:2])
    assert model is fixed
    assert kw['pol'] == 'I' and kw['sites'] == ['AZ', 'LM', 'PV'] and kw['caltable'] is True
    assert kw['processes'] == -1
    cal.applycal.assert_called_once()
    np.testing.assert_array_equal(native.data, original_data)
    receipt = json.loads((tmp_path / 'polarization-input.json').read_text())
    assert receipt['shared_R_equals_L'] and receipt['no_I_noise_inflation'] and receipt['no_reverse_taper']
    assert receipt['rf_hz'] == native.rf
    assert receipt['uv_min'] == 1e8 and receipt['zero_baselines_excluded_before_shared_gain_solve']
    assert receipt['official_I_zero_baseline_workflow_unchanged']
    assert (tmp_path / 'clean-gaincal-pol.uvfits').is_file()


def test_i_rejects_nonshared_gains(tmp_path, config, monkeypatch):
    source = tmp_path / 'official.py'
    source.write_text("ttype     = 'nfft'\nfov       = 128*eh.RADPERUAS\n")
    config['official_script'] = str(source)
    native = Observation()
    bad = np.array([(1.+0j, 2.+0j)], dtype=[('rscale', complex), ('lscale', complex)])
    cal = SimpleNamespace(data={s: bad for s in ('AZ', 'LM', 'PV')}, applycal=Mock())
    eh = SimpleNamespace(selfcal=Mock(return_value=cal))
    monkeypatch.setattr(w.runpy, 'run_path', Mock(return_value=dict(im_out=Image())))
    monkeypatch.setattr(w, 'validate_fixed', Mock())
    with pytest.raises(AssertionError):
        w.run_i(dict(parameters=config, raw='synthetic'),
                SimpleNamespace(load_native=lambda *a: (native, {})), eh, tmp_path)
    cal.applycal.assert_not_called()


def test_fixed_d_never_fits_and_fsys_once_before_d_on_fresh_base(config, monkeypatch):
    base, fixed = Observation(), Image()
    original = base.data.copy()
    events = []

    def apply(obs, d):
        assert obs is not base and d == {'AZ': [[.01, .02], [.03, .04]]}
        assert obs.noisy and np.all(obs.data['rlsigma'] == 3 * original['rlsigma'])
        obs.data['rlsigma'] *= 2
        events.append('D')
        return obs

    def noise(obs, fsys):
        assert fsys == .005 and obs is not base and not getattr(obs, 'noisy', False)
        np.testing.assert_array_equal(obs.data, original)
        obs.data['rlsigma'] *= 3
        obs.noisy = True
        events.append('noise')
        return obs

    fit = Mock(side_effect=AssertionError('D fit forbidden in uncertainty'))
    core = SimpleNamespace(apply_fixed_d=Mock(side_effect=apply), add_fsys=Mock(side_effect=noise),
                           make_imaging_obs=lambda obs: obs, fit_forward_d=fit)

    def image_p(obs, *args):
        assert np.all(obs.data['rlsigma'] == 6 * original['rlsigma'])
        return fixed, []

    monkeypatch.setattr(w, 'image_p', Mock(side_effect=image_p))
    for _ in range(2):
        w.fixed_d_image(base, fixed, {'AZ': [[.01, .02], [.03, .04]]}, config, core, None)
    assert events == ['noise', 'D', 'noise', 'D']
    fit.assert_not_called()
    assert w.image_p.call_count == 2
    np.testing.assert_array_equal(base.data, original)


def test_weighted_fraction_restart_keeps_exact_i_and_domain(config):
    fixed = Image()
    image = w.initial_image(fixed, config)
    image.qvec = .99 * fixed.imvec
    image.uvec = .01 * fixed.imvec
    restart = w.restart(image, fixed, config, SimpleNamespace(RADPERUAS=1.))
    np.testing.assert_array_equal(restart.imvec, fixed.imvec)
    assert np.all(np.hypot(restart.qvec, restart.uvec) < restart.imvec)
    assert np.all(restart.vvec == 0)


def test_sampler_four_day_covariance_ddof_and_reproducible_full_matrix(config):
    from eht_axion.research import eht_official as core

    days = {str(day): {s: [[.01 * k, .02 * k], [-.02 * k, .01 * k]]
                      for s in config['solve_sites']} for k, day in enumerate(w.DAYS)}
    samples = w.sample_dterms(days, config, 95, core)
    again = w.sample_dterms(days, config, 95, core)
    assert samples == again
    assert len(samples['draws']) == 1000 and len(samples['sites']) == 7
    assert samples != w.sample_dterms(days, config, 96, core)
    values = np.array([[.01 * k, .02 * k] for k in range(4)])
    np.testing.assert_allclose(samples['distributions']['AZ']['R']['covariance'], np.cov(values.T, ddof=1))
    assert samples['residual_component_sigma'] == .01
    config['residual_d_convention'] = 'complex_rms'
    rms = w.sample_dterms(days, config, 95, core)
    assert rms['residual_component_sigma'] == .01 / np.sqrt(2)


def test_pilot_saves_all_draws_and_only_executed_prefix(tmp_path, config, monkeypatch):
    reference = dict(source_kind='Published rounded values, not private posterior',
                     stored_units='dimensionless',
                     days={str(day): {s: [[.01, .02], [.03, .04]] for s in config['solve_sites']}
                           for day in w.DAYS})
    src = tmp_path / 'table.json'
    src.write_text(json.dumps(reference))
    args = SimpleNamespace(published_dterms=src, d_roots=None)
    fixed = Image()
    executor = Mock(return_value=(fixed, []))
    monkeypatch.setattr(w, 'fixed_d_image', executor)
    from eht_axion.research import eht_official as core

    w.run_uncertainty(dict(day=95, draws=2, label='engineering-partial-pilot'), args,
                      Observation(), fixed, config, core, SimpleNamespace(RADPERUAS=1.), tmp_path)
    assert executor.call_count == 2
    with np.load(tmp_path / 'all-draws.npz', allow_pickle=False) as z:
        assert z['dterms'].shape == (1000, 7, 2)
        np.testing.assert_array_equal(z['draw_id'], np.arange(1000))
    receipt = json.loads((tmp_path / 'sampling-receipt.json').read_text())
    assert receipt['executed_draws'] == 2 and receipt['configured_draws'] == 1000
    assert receipt['D_refitted'] is False
    assert not (tmp_path / 'draw-0002-raw-IQU.fits').exists()
    assert json.loads(src.read_text()) == reference  # protected D source unchanged


def test_p_imager_shift_leaves_gradient_fixed_i_weights_and_numerical_stop(config, monkeypatch):
    import scipy.optimize as opt

    fixed, obs, imager_calls = Image(), Observation(), []

    class Imager:
        def __init__(self, passed_obs, init, **kwargs):
            assert passed_obs is obs
            imager_calls.append(kwargs)
            self.image = init.copy()

        def objfunc(self, x):
            return float(np.dot(x, x))

        def objgrad(self, x):
            return 2 * x

        def make_image_P(self, **kwargs):
            x = np.array([.2, .3])
            np.testing.assert_allclose(self.objgrad(x), 2 * x)
            np.testing.assert_allclose(self.objfunc(x), np.dot(x, x) + 2 + 100 * np.log(2))
            opt.minimize(self.objfunc, x, jac=self.objgrad)
            return self.image

    def minimize(fun, x, **kwargs):
        return SimpleNamespace(fun=fun(x), jac=kwargs['jac'](x), nit=4, success=True,
                               status=0, message='synthetic optimizer exit')

    monkeypatch.setattr(opt, 'minimize', minimize)
    eh = SimpleNamespace(imager=SimpleNamespace(Imager=Imager))
    image, records = w.image_p(obs, fixed, fixed, config, eh)
    np.testing.assert_array_equal(image.imvec, fixed.imvec)
    call = imager_calls[0]
    assert call['maxit'] == 500 and call['stop'] == 1e-8
    assert call['norm_init'] is False and call['norm_reg'] is True
    assert call['pol'] == 'P' and call['transform'] == ['mcv']
    assert call['data_term'] == {'pvis': 1., 'm': 1.}
    assert call['reg_term'] == {'hw': 100., 'ptv': 1.}
    assert records[0]['nit'] == 4 and records[0]['transformed_gradient_inf'] == .6


def test_fiducial_total_d_solves_against_fresh_base_and_final_p_only(tmp_path, config, monkeypatch):
    fixed, base, fits, applications = Image(), Observation(), [], []
    config['rounds'], config['restart_rounds'] = 3, [2]
    tarr = np.array([(s, .01+.02j, -.03+.04j) for s in config['solve_sites']],
                    dtype=[('site', 'U2'), ('dr', complex), ('dl', complex)])

    def add_fsys(obs, fraction):
        assert not getattr(obs, 'noisy', False)
        obs.noisy = True
        return obs

    def solve(obs, image, sites):
        assert obs is not base and obs.noisy and not getattr(obs, 'corrected', False)
        assert sites == config['solve_sites']
        fits.append(obs)
        return SimpleNamespace(tarr=tarr)

    def apply(obs, dterms):
        assert obs is not base and obs.noisy and not getattr(obs, 'corrected', False)
        obs.corrected = True
        applications.append(obs)
        return obs

    core = SimpleNamespace(add_fsys=add_fsys, fit_forward_d=solve, apply_fixed_d=apply,
                           make_imaging_obs=lambda obs: obs)
    p = Mock(side_effect=lambda obs, fixed, init, c, eh: (init, [dict(nit=2)]))
    monkeypatch.setattr(w, 'image_p', p)
    w.run_fiducial(base, fixed, config, core, SimpleNamespace(RADPERUAS=1.), tmp_path)
    assert len(fits) == 3 and len(applications) == 3  # rounds 2/3 plus final P-only
    assert p.call_count == 4
    assert not getattr(base, 'noisy', False) and not getattr(base, 'corrected', False)
    final = json.loads((tmp_path / 'final-dterms.json').read_text())
    assert final['PV'] == [[.01, .02], [-.03, .04]]
    assert json.loads((tmp_path / 'steps.json').read_text())[-1]['stage'] == 'final-P-only'
    assert not json.loads((tmp_path / 'termination.json').read_text())['scientifically_converged']


def test_source_stage_must_declare_consumed_products(tmp_path, monkeypatch):
    verified = Mock(return_value=dict(stage='research:eht-official:I', config=dict(day=95)))
    monkeypatch.setattr(w.common, 'verify_stage', verified)
    declared = Mock(side_effect=ValueError('Consumed product is not a declared stage output'))
    monkeypatch.setattr(w.common, 'require_declared_outputs', declared)
    with pytest.raises(ValueError, match='not a declared'):
        w.stage_inputs(tmp_path, 'I', 95, w.I_PRODUCTS)
    verified.assert_called_once_with(tmp_path, recursive=True, allow_imported=False)
    declared.assert_called_once_with(tmp_path, [tmp_path / p for p in w.I_PRODUCTS])
