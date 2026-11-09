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
SPEC = importlib.util.spec_from_file_location('eht_reconstruct_test',
                                               ROOT / 'scripts/eht_reconstruct.py')
w = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(w)


LOCATIONS = dict(input_root='/data/uvfits', input_receipt='/data/receipt.json',
                 feed_rotation_file='/data/feed.json', official_script='/data/official.py')


@pytest.fixture
def config():
    """The flat eht-imaging parameter set derived from experiment.yaml."""
    return w.experiment.eht_imaging_parameters(w.experiment.read_experiment(), LOCATIONS)


@pytest.fixture
def plan_args(tmp_path, monkeypatch):
    """Fully local synthetic receipts; no GRMHD inputs and no optional package import."""
    experiment = w.experiment.read_experiment()
    raw = tmp_path / experiment['eht_imaging']['inputs'][95][0]
    raw.write_bytes(b'synthetic circular UVFITS placeholder')
    experiment['eht_imaging']['inputs'][95][1] = w.provenance.sha256(raw)
    receipt = tmp_path / 'receipt.json'
    receipt.write_text(json.dumps(dict(complete=True, files=[dict(path=str(raw),
        size=raw.stat().st_size, sha256=experiment['eht_imaging']['inputs'][95][1])])) )
    feed, official = tmp_path / 'feed.json', tmp_path / 'official.py'
    feed.write_text('{}')
    official.write_text("ttype     = 'nfft'\nfov       = 128*eh.RADPERUAS\n")
    path = tmp_path / 'experiment.yaml'
    path.write_text(yaml.safe_dump(experiment))
    dist = SimpleNamespace(version='1.3.2', files=[], locate_file=lambda f: f)
    monkeypatch.setattr(w.metadata, 'distribution', lambda name: dist)
    return w.parser().parse_args(['--config', str(path), '--stage', 'I', '--day', '95',
                                  '--input-root', str(tmp_path), '--input-receipt', str(receipt),
                                  '--feed-rotation-file', str(feed), '--official-script', str(official),
                                  '--no-env-file', '--output-root', str(tmp_path / 'new'), '--plan'])


def test_plan_no_optional_import_no_writes(plan_args, monkeypatch):
    monkeypatch.setitem(sys.modules, 'ehtim', None)
    monkeypatch.setitem(sys.modules, 'astropy', None)
    plan = w.build_plan(plan_args)
    assert plan['day'] == 95 and plan['stage'] == 'I'
    assert plan['scientifically_converged'] is False
    assert not plan_args.output_root.exists()
    paths = {Path(r['path']) for r in plan['inputs']}
    assert ROOT / 'scripts/eht_reconstruct.py' in paths
    assert ROOT / 'eht_axion/provenance.py' in paths
    assert plan_args.config in paths
    assert all(len(r['sha256']) == 64 for r in plan['inputs'])
    assert plan['parameters']['input_root'] == str(plan_args.input_root)
    assert (plan['parameters']['npix'], plan['parameters']['fov_uas']) == (64, 120.0)
    assert plan['display'] == 'additional 19uas I/Q/U kernel, NOT an effective PSF'


def test_locations_are_required_absolute_paths(plan_args):
    plan_args.official_script = None
    with pytest.raises(ValueError, match='official_script'):
        w.build_plan(plan_args)
    plan_args.official_script = Path('relative.py')
    with pytest.raises(ValueError, match='absolute'):
        w.build_plan(plan_args)


def test_cli_review_and_stage_contract(plan_args):
    plan_args.plan = False
    with pytest.raises(ValueError, match='review-file'):
        w.build_plan(plan_args)
    plan_args.plan = True
    plan_args.stage = 'fiducial'
    with pytest.raises(ValueError, match='i-root'):
        w.build_plan(plan_args)
    plan_args.stage, plan_args.i_root = 'I', Path('not-read-yet')
    with pytest.raises(ValueError, match='i-root'):
        w.build_plan(plan_args)
    with pytest.raises(SystemExit):
        w.parser().parse_args(['--stage', 'uncertainty', '--day', '95', '--output-root', 'x'])


def test_wrong_hash_existing_destination_and_output_source_protection(plan_args):
    plan_args.output_root.mkdir()
    with pytest.raises(FileExistsError, match='new stage/day'):
        w.build_plan(plan_args)
    plan_args.output_root.rmdir()
    plan_args.output_root = plan_args.config.parent
    with pytest.raises(FileExistsError):
        w.build_plan(plan_args)
    plan_args.output_root = plan_args.config.parent / 'different'
    experiment = yaml.safe_load(plan_args.config.read_text())
    raw = plan_args.config.parent / experiment['eht_imaging']['inputs'][95][0]
    raw.write_bytes(b'tampered')
    with pytest.raises(ValueError, match='checksum'):
        w.build_plan(plan_args)


@pytest.mark.parametrize('text', ['', "ttype = 'nfft'\nfov       = 128*eh.RADPERUAS",
    "ttype     = 'nfft'\n" * 2 + 'fov       = 128*eh.RADPERUAS',
    "ttype     = 'nfft'\nfov = 128*eh.RADPERUAS"])
def test_adapt_requires_unique_exact_anchors(text):
    with pytest.raises(ValueError, match='anchor'):
        w.adapt_official(text, 120.0)


def test_adaptation_changes_only_two_exact_strings():
    text = ("# nfft 128 remains here\nttype     = 'nfft' # retained\n"
            'fov       = 128*eh.RADPERUAS # retained\nrest = 5\n')
    expected = text.replace("ttype     = 'nfft'", "ttype     = 'direct'")
    assert w.adapt_official(text, 120.0) == expected.replace(
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
                SimpleNamespace(load_native=lambda *a, **k: (native, {})), eh, tmp_path)
    cal.applycal.assert_not_called()


def test_weighted_fraction_restart_keeps_exact_i_and_domain(config):
    fixed = Image()
    image = w.initial_image(fixed, config)
    image.qvec = .99 * fixed.imvec
    image.uvec = .01 * fixed.imvec
    restart = w.restart(image, fixed, config, SimpleNamespace(RADPERUAS=1.))
    np.testing.assert_array_equal(restart.imvec, fixed.imvec)
    assert np.all(np.hypot(restart.qvec, restart.uvec) < restart.imvec)
    assert np.all(restart.vvec == 0)


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

    def solve(obs, image, sites, optimizer_options):
        assert obs is not base and obs.noisy and not getattr(obs, 'corrected', False)
        assert sites == config['solve_sites']
        assert optimizer_options is None  # the fiducial stage uses the official D optimizer options
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
    verified = Mock(return_value=dict(stage='eht:i', config=dict(day=95)))
    monkeypatch.setattr(w.provenance, 'verify_stage', verified)
    declared = Mock(side_effect=ValueError('Consumed product is not a declared stage output'))
    monkeypatch.setattr(w.provenance, 'require_declared_outputs', declared)
    with pytest.raises(ValueError, match='not a declared'):
        w.stage_inputs(tmp_path, 'I', 95, w.I_PRODUCTS)
    verified.assert_called_once_with(tmp_path, recursive=True)
    declared.assert_called_once_with(tmp_path, [tmp_path / p for p in w.I_PRODUCTS])
