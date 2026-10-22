"""No real imaging: pinned planning, strict safety, and tiny upstream API doubles."""
import copy
from contextlib import contextmanager
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

import rebuild_eht_fiducial as r


@pytest.fixture
def planned(tmp_path):
    cfg = yaml.safe_load((r.ROOT / 'configs/eht-fiducial.yaml').read_text())
    uv = tmp_path / r.NAME
    uv.write_bytes(b'synthetic UVFITS receipt test, not science')
    receipt = tmp_path / 'receipt.json'
    receipt.write_text(json.dumps(dict(commit=r.COMMIT, relative_path='hops_data/April11/' + r.NAME,
        source='https://github.com/eventhorizontelescope/2023-D01-01', bytes=uv.stat().st_size,
        sha256=r.digest(uv))))
    feed = tmp_path / 'feed.json'
    feed.write_text(json.dumps({'table': [dict(site=s, fr_par=1, fr_elev=0, fr_off=0)
                                        for s in ['AA', 'AP', 'AZ', 'LM', 'PV', 'SM', 'JC']]}))
    cfg.update(input_uvfits=str(uv), input_receipt=str(receipt), input_sha256=r.digest(uv),
               feed_rotation_file=str(feed))
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(cfg))
    return r.parser().parse_args(['--config', str(path), '--output-root', str(tmp_path / 'new'),
                                 '--max-rounds', '2', '--maxit', '30'])


def test_plan_is_no_optional_import_no_writes(planned):
    plan = r.main(['--config', str(planned.config), '--output-root', str(planned.output_root),
                   '--plan', '--max-rounds', '2', '--maxit', '30'])
    assert plan['label'] == 'engineering-partial'
    assert plan['parameters']['i_maxit'] == plan['parameters']['maxit'] == 30
    assert plan['parameters']['max_rounds'] == 2
    assert plan['parameters']['ttype'] == 'direct'
    assert not planned.output_root.exists()
    assert plan['hashes'][str(planned.config)] == r.digest(planned.config)
    check = subprocess.run([sys.executable, '-S', str(Path(r.__file__)), '--help'], capture_output=True)
    assert check.returncode == 0


def test_full_defaults(planned):
    planned.max_rounds = planned.maxit = None
    p = r.build_plan(planned)
    assert p['label'] == 'official-inspired-full-budget'
    assert p['parameters']['max_rounds'] == 50
    assert p['parameters']['maxit'] == p['parameters']['i_maxit'] == 100
    assert p['i_data'] == {'amp': .2, 'cphase': 1., 'logcamp': 1.}
    assert p['p_reg'] == {'hw': 100., 'ptv': 1.}


@pytest.mark.parametrize('value', ['100', True, float('nan'), 0, 1.5])
def test_strict_config(planned, value):
    c = yaml.safe_load(planned.config.read_text())
    c['maxit'] = value
    planned.config.write_text(yaml.safe_dump(c))
    with pytest.raises(ValueError):
        r.build_plan(planned)


def test_receipt_tamper_and_limits(planned):
    planned.max_rounds = 51
    with pytest.raises(ValueError, match='override'):
        r.build_plan(planned)
    planned.max_rounds = 2
    c = yaml.safe_load(planned.config.read_text())
    Path(c['input_uvfits']).write_bytes(b'tampered')
    with pytest.raises(ValueError, match='receipt'):
        r.build_plan(planned)


def test_output_safety(planned):
    planned.output_root.mkdir()
    with pytest.raises(FileExistsError):
        r.build_plan(planned)
    planned.output_root.rmdir()
    planned.output_root.symlink_to(planned.config.parent / 'absent', target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        r.build_plan(planned)
    planned.output_root = planned.config.parent
    with pytest.raises(ValueError, match='overlap'):
        r.build_plan(planned)


class Obs:
    def __init__(self):
        # JC deliberately absent; fixed sites must be intersected, not required.
        self.data = np.zeros(6, dtype=[('t1', 'U2'), ('t2', 'U2'), ('vis', complex),
                                      ('qsigma', float), ('usigma', float)])
        self.data['t1'] = ['AA', 'AP', 'AZ', 'LM', 'PV', 'SM']
        self.data['t2'] = 'AA'
        self.data['vis'], self.data['qsigma'], self.data['usigma'] = 2, 3, 4
        self.tarr = np.zeros(6, dtype=[('site', 'U2'), ('dr', complex), ('dl', complex),
                                      ('fr_par', float), ('fr_elev', float), ('fr_off', float)])
        self.tarr['site'] = self.data['t1']
        self.mjd, self.source, self.rf, self.dcal = 57854, 'M87', 227070703125., True

    def copy(self):
        return copy.deepcopy(self)

    def switch_polrep(self, *a, **kw):
        return self.copy()

    def add_scans(self):
        pass

    def reorder_tarr_snr(self):
        pass

    def avg_coherent(self, seconds, scan_avg):
        assert seconds == 0 and scan_avg
        return self.copy()

    def flag_uvdist(self, uv_min):
        assert uv_min == 1e8
        return self.copy()

    def res(self):
        return 25.


class Image:
    blurs = []

    def __init__(self):
        self.imvec = np.ones(4096) * .6 / 4096
        self.qvec = self.uvec = np.zeros(4096)

    def copy(self):
        return copy.deepcopy(self)

    def add_gauss(self, *a):
        return self.copy()

    def blur_circ(self, *args):
        self.blurs.append(args)
        return self.copy()

    def add_qu(self, q, u):
        self.qvec, self.uvec = q.ravel(), u.ravel()

    def add_pol_image(self, *a):
        pass

    def save_fits(self, path):
        Path(path).write_text('synthetic IQU output')


def test_noise_not_cumulative():
    obs = Obs()
    once = r.pol_noise(obs, .005)
    assert np.allclose(once.data['qsigma']**2 - 9, .01**2 / 2)
    assert np.allclose(once.data['usigma']**2 - 16, .01**2 / 2)
    assert np.array_equal(r.pol_noise(obs, .005).data, once.data)
    assert np.all(obs.data['qsigma'] == 3)


def test_two_round_api_contract(planned, monkeypatch, tmp_path):
    plan = r.build_plan(planned)
    out = tmp_path / 'mock-stage'
    out.mkdir()
    seen = []

    class Imager:
        def __init__(self, obs, init, prior, flux, **kw):
            self.image, self.kw = init.copy(), kw
            seen.append(kw)

        def make_image_I(self, **kw):
            return self.image

        def make_image_P(self, **kw):
            assert self.kw['pol'] == 'P' and self.kw['transform'] == ['mcv']
            self.image.qvec *= 1.1
            return self.image

    def leakage(obs, image, **kw):
        assert np.all(obs.tarr['dr'] == 0)  # never previous calibrated obs
        assert kw == dict(sites=['AZ', 'LM', 'PV'], leakage_tol=.1,
                          ttype='direct', show_solution=False)
        result = obs.copy()
        result.tarr['dr'][2:5] = .01j
        result.dcal = True
        return result

    @contextmanager
    def records():
        yield [dict(nit=3, success=False, status=1, message='mock limit', objective=1.)]

    monkeypatch.setattr(r, 'optimizer_records', records)
    monkeypatch.setattr(r, 'pol_chi2', lambda *a: dict(pvis=1., m=2.))
    eh = SimpleNamespace(RADPERUAS=1., obsdata=SimpleNamespace(load_uvfits=lambda *a, **k: Obs()),
        image=SimpleNamespace(make_square=lambda *a: Image()), imager=SimpleNamespace(Imager=Imager),
        selfcal=lambda obs, *a, **kw: obs.copy())
    Image.blurs.clear()
    r.run_rebuild(plan, eh, SimpleNamespace(leakage_cal=leakage), out)
    termination = json.loads((out / 'termination.json').read_text())
    assert termination['rounds'] == 2 and not termination['engineering_stable']
    assert not termination['scientifically_converged']
    assert len(seen) == 6  # initial I + 3 blur cycles + two P steps
    assert all(kw['maxit'] == 30 for kw in seen)
    assert (19., 19.) in Image.blurs
    assert len(list(out.glob('round-*-dterms.json'))) == 2
    assert (out / 'final-raw-IQU.fits').exists()
    row = json.loads((out / 'round-02.json').read_text())
    assert row['optimizers'][0]['nit'] == 3 and row['chi2_after_D']['m'] == 2


def test_optimizer_records_actual_steps_and_exception():
    import scipy.optimize as opt

    original = opt.minimize
    with r.optimizer_records() as records:
        opt.minimize(lambda x: (x[0] - 1)**2, [0.], method='L-BFGS-B')
    assert records[0]['nit'] > 0 and opt.minimize is original
    with pytest.raises(RuntimeError, match='upstream'):
        with r.optimizer_records():
            raise RuntimeError('upstream')
    assert opt.minimize is original
