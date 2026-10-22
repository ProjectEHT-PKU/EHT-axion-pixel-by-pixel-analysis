"""Fast synthetic controls only: no science cases or optional ehtim import."""
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

import validate_eht_apr11 as r
import rebuild_eht_fiducial as reused
from test_eht_fiducial_plan import Obs, Image


@pytest.fixture
def planned(tmp_path, monkeypatch):
    c = yaml.safe_load((r.ROOT/'configs/eht-apr11-validation.yaml').read_text())
    source = tmp_path/'source'
    source.mkdir()
    for name in (r.RAW, r.UV):
        (source/name).write_bytes(b'synthetic only')
    feed = tmp_path/'feed.json'
    feed.write_text(json.dumps(dict(status='verified', table=[
        dict(site=s, fr_par=p, fr_elev=e, fr_off=o) for s, p, e, o in
        [('AA',1,0,0), ('AP',1,1,0), ('AZ',1,1,0), ('LM',1,-1,0),
         ('PV',1,-1,0), ('SM',1,-1,45), ('JC',1,0,0)]])))
    c.update(i_root=str(source), feed_file=str(feed))
    config = tmp_path/'config.yaml'
    config.write_text(yaml.safe_dump(c))
    monkeypatch.setattr(r, 'verified_products', lambda root, stage, names: [root/n for n in names])
    return r.parser().parse_args(['--config', str(config), '--case', 'baseline',
                                 '--output-root', str(tmp_path/'new')])


def test_planning_controls_no_optional_imports_or_output(planned):
    a = r.build_plan(planned)
    assert a['rounds'] == 100 and a['evpa_deg'] == 0 and a['restart_rounds'] == []
    assert not a['audit_template']['working_model_approved']
    for case, evpa, restarts in [('evpa45',45,[]), ('restart',0,[2,12,22,32,42])]:
        planned.case = case
        b = r.build_plan(planned)
        assert b['selection_sha256'] == a['selection_sha256']
        assert b['evpa_deg'] == evpa and b['restart_rounds'] == restarts
        assert b['parameters'] == a['parameters'] and b['hashes'] == a['hashes']
    assert 'ehtim' not in sys.modules and not planned.output_root.exists()
    help_result = subprocess.run([sys.executable, '-S', r.__file__, '--help'], capture_output=True)
    assert help_result.returncode == 0


@pytest.mark.parametrize('key,value', [('unknown',1), ('schema',2), ('maxit',True), ('fsys',.02), ('norm_reg','true')])
def test_reject_config(planned, key, value):
    c = yaml.safe_load(planned.config.read_text())
    c[key] = value
    planned.config.write_text(yaml.safe_dump(c))
    with pytest.raises(ValueError):
        r.build_plan(planned)


def test_norm_selection_and_rounds(planned):
    a = r.build_plan(planned)
    planned.norm_reg = False
    assert r.build_plan(planned)['selection_sha256'] != a['selection_sha256']
    planned.rounds = 101
    with pytest.raises(ValueError):
        r.build_plan(planned)


def test_safety(planned):
    planned.output_root.mkdir()
    with pytest.raises(FileExistsError):
        r.build_plan(planned)
    planned.output_root.rmdir()
    planned.output_root.symlink_to(planned.config.parent/'absent', target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        r.build_plan(planned)
    planned.output_root.unlink()
    planned.output_root = Path(yaml.safe_load(planned.config.read_text())['i_root'])/'nested'
    with pytest.raises(ValueError, match='overlap'):
        r.build_plan(planned)


def test_feed_degrees_reject_radians(planned):
    feed = Path(yaml.safe_load(planned.config.read_text())['feed_file'])
    assert r.feed_table(feed)['SM'][2] == 45
    data = json.loads(feed.read_text())
    data['table'][5]['fr_off'] = np.pi/4
    feed.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='DEGREES'):
        r.feed_table(feed)


def test_projection_fixed_i_domain():
    fixed, image = Image(), Image()
    image.qvec = fixed.imvec*4
    init, count = r.projected_init(image, fixed, 20.)
    assert count == 4096
    assert np.array_equal(init.imvec, fixed.imvec)
    assert np.all(np.hypot(init.qvec, init.uvec) < fixed.imvec)
    assert np.all(image.qvec == fixed.imvec*4)


def test_gate_no_self_approval(planned, monkeypatch):
    p = r.build_plan(planned)
    audit = planned.config.parent/'audit.json'
    r.dump(audit, p['audit_template'])
    args = ['--config',str(planned.config), '--case','baseline', '--output-root',str(planned.output_root)]
    with pytest.raises(ValueError, match='BEFORE'):
        r.main(args)
    with pytest.raises(ValueError, match='parent audit'):
        r.main(args+['--audit-file',str(audit)])
    assert not planned.output_root.exists()


def test_two_round_noise_once_actual_records_and_outputs(planned, monkeypatch):
    planned.case, planned.rounds = 'restart', 2
    plan = r.build_plan(planned)
    out = planned.config.parent/'synthetic-output'
    out.mkdir()
    noise_calls, starts = [], []
    original_noise = reused.pol_noise
    monkeypatch.setattr(reused, 'pol_noise', lambda obs, f: (noise_calls.append(f) or original_noise(obs,f)))
    monkeypatch.setattr(reused, 'pol_chi2', lambda *a: dict(pvis=1.,m=2.))
    monkeypatch.setattr(r, 'objective', lambda *a: dict(total=3.,regularizers={'hw':.1,'ptv':.2}))
    monkeypatch.setattr(r, 'residuals', lambda *a: {'ALL': {'rows':6}})

    class Fixed(Image):
        def regrid_image(self, fov, npix):
            assert (fov,npix) == (120.,64)
            return self.copy()

    class Imager:
        def __init__(self, obs, init, prior, flux, **kw):
            self.image = init.copy()
            starts.append(copy.deepcopy(kw))
            assert np.allclose(obs.switch_polrep('stokes').data['qsigma']**2, 9+.01**2/2)
        def make_image_P(self, **kw):
            return self.image

    @contextmanager
    def records():
        yield [dict(nit=3,status=1,success=False,objective=3.)]

    def leakage(obs, image, **kw):
        assert np.all(obs.tarr['dr'] == 0) and kw['inverse'] is False
        assert np.allclose(obs.data['qsigma']**2, 9+.01**2/2)
        return obs.copy()

    monkeypatch.setattr(reused, 'optimizer_records', records)
    eh = SimpleNamespace(RADPERUAS=1., image=SimpleNamespace(load_fits=lambda *a: Fixed()),
        obsdata=SimpleNamespace(load_uvfits=lambda *a,**kw: Obs()), imager=SimpleNamespace(Imager=Imager))
    r.run(plan, eh, SimpleNamespace(leakage_cal=leakage), out)
    assert noise_calls == [.005] and len(starts) == 2
    assert all(s['norm_reg'] and not s['norm_init'] and s['stop'] == 1e-4 for s in starts)
    row = json.loads((out/'round-002.json').read_text())
    assert row['blur_restart'] and row['optimizers'][0]['nit'] == 3
    assert all(pair == [[0.,0.],[0.,0.]] for pair in row['dterms'].values())
    assert (out/'round-002-raw-IQU.fits').exists() and not (out/'round-001-raw-IQU.fits').exists()
    assert json.loads((out/'termination.json').read_text())['rounds'] == 2
