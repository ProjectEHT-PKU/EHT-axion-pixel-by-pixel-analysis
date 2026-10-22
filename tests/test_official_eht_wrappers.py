"""Wrapper contracts only: no released sources, optional stack, or real imaging."""

import ast
import builtins
from contextlib import contextmanager
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock, call

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / 'scripts/run_official_eht_i_clean_pol.py'


def load_wrapper(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT))
    monkeypatch.syspath_prepend(str(ROOT / 'scripts'))
    spec = importlib.util.spec_from_file_location('_official_clean_pol_test', WRAPPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def arguments(tmp_path):
    return [
        '--official-script', str(tmp_path / 'official.py'),
        '--input-uvfits', str(tmp_path / 'input.uvfits'),
        '--input-receipt', str(tmp_path / 'receipt.json'),
        '--output-root', str(tmp_path / 'output'),
    ]


def test_import_does_not_load_optional_stack(monkeypatch):
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name.split('.')[0] in {'ehtim', 'astropy'}:
            pytest.fail(f'Optional dependency imported before main: {name}')
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', guarded_import)
    assert callable(load_wrapper(monkeypatch).main)


def test_existing_output_and_wrong_version_fail_before_work(monkeypatch, tmp_path):
    wrapper = load_wrapper(monkeypatch)
    version = Mock(return_value='1.3.1')
    stage = Mock(side_effect=AssertionError('Must not stage rejected runs'))
    monkeypatch.setattr(wrapper, 'version', version)
    monkeypatch.setattr(wrapper, 'stage_directory', stage)
    monkeypatch.setitem(sys.modules, 'ehtim', None)
    monkeypatch.setitem(sys.modules, 'astropy', None)
    output = tmp_path / 'output'
    output.mkdir()
    sentinel = output / 'keep.txt'
    sentinel.write_text('do not overwrite')
    # All input paths are deliberately absent: guards must precede any reads.
    with pytest.raises(FileExistsError, match='New official-I output'):
        wrapper.main(arguments(tmp_path))
    version.assert_not_called()
    assert sentinel.read_text() == 'do not overwrite'
    sentinel.unlink()
    output.rmdir()
    with pytest.raises(ValueError, match='Requires pinned ehtim1.3.2'):
        wrapper.main(arguments(tmp_path))
    version.assert_called_once_with('ehtim')
    stage.assert_not_called()
    assert not output.exists()


def test_backend_anchor_requires_exactly_one_exact_match(monkeypatch, tmp_path):
    wrapper = load_wrapper(monkeypatch)
    monkeypatch.setattr(wrapper, 'version', lambda name: '1.3.2')
    stage = Mock(side_effect=AssertionError('Must not stage invalid source'))
    monkeypatch.setattr(wrapper, 'stage_directory', stage)
    monkeypatch.setitem(sys.modules, 'ehtim', None)
    monkeypatch.setitem(sys.modules, 'astropy', None)
    for text in ('# absent\n', "ttype = 'nfft'\n", "ttype     = 'nfft'\n" * 2):
        (tmp_path / 'official.py').write_text(text)
        with pytest.raises(ValueError, match='adaptation anchor differs'):
            wrapper.main(arguments(tmp_path))
        assert (tmp_path / 'official.py').read_text() == text
        assert not (tmp_path / 'output').exists()
    stage.assert_not_called()


def test_official_execution_exports_clean_original_gaincal_checkpoint(monkeypatch, tmp_path):
    wrapper = load_wrapper(monkeypatch)
    monkeypatch.setattr(wrapper, 'version', Mock(return_value='1.3.2'))
    # Tiny stand-ins only; never scan or hash an installed optional package.
    package_code = tmp_path / 'fake-ehtim/pkg.py'
    package_code.parent.mkdir()
    package_code.write_text('# synthetic package code\n')
    ehtim = ModuleType('ehtim')
    ehtim.__file__, ehtim.RADPERUAS = str(package_code), 2.
    astropy = ModuleType('astropy')
    utils = ModuleType('astropy.utils')
    iers = ModuleType('astropy.utils.iers')
    iers.conf = SimpleNamespace(auto_download=True)
    astropy.utils, utils.iers = utils, iers
    for name, module in [('ehtim', ehtim), ('astropy', astropy),
                         ('astropy.utils', utils), ('astropy.utils.iers', iers)]:
        monkeypatch.setitem(sys.modules, name, module)

    source = "# synthetic official script\nttype     = 'nfft'\n# nfft stays in comments\n"
    official = tmp_path / 'official.py'
    official.write_text(source)
    uv = tmp_path / 'input.uvfits'
    uv.write_text('synthetic input')
    receipt = tmp_path / 'receipt.json'
    receipt.write_text('{}')
    output = tmp_path / 'output'
    staged = []

    @contextmanager
    def stage(out, kind, config, inputs):
        staged.append((kind, config, inputs))
        out.mkdir()
        yield out

    monkeypatch.setattr(wrapper, 'stage_directory', stage)

    def save(label):
        return Mock(side_effect=lambda path: Path(path).write_text(label))

    image = SimpleNamespace(save_fits=save('raw I'))
    imgr = SimpleNamespace(out_last=Mock(return_value=image))
    class Observation:
        # No taper/noise methods: only copying and station reordering are allowed.
        def __init__(self):
            self.data = np.array(
                [('AZ', 'LM', 2.+3.j, .1), ('PV', 'AZ', 4.-1.j, .3)],
                dtype=[('t1', 'U2'), ('t2', 'U2'), ('vis', complex), ('sigma', float)],
            )
            self.tarr = np.array([('AZ', 11.), ('LM', 22.), ('PV', 33.)],
                                 dtype=[('site', 'U2'), ('x', float)])

        def copy(self):
            return deepcopy(self)

    obs_orig = Observation()
    original_data, original_tarr = obs_orig.data.copy(), obs_orig.tarr.copy()
    obs_sc = SimpleNamespace(save_uvfits=save('selfcal'), chisq=Mock(return_value=3.))
    clean = SimpleNamespace(save_uvfits=save('clean polarization'))

    def applycal(obs, **kwargs):
        assert obs is not obs_orig and obs is not obs_sc
        # Copying/reordering must not retaper visibilities or inflate errors.
        np.testing.assert_array_equal(obs.data, original_data)
        assert obs.data.tobytes() == original_data.tobytes()
        np.testing.assert_array_equal(obs.tarr, original_tarr[[2, 0, 1]])
        np.testing.assert_array_equal(obs.tarr['site'], caltab.tarr['site'])
        assert kwargs == {'interp': 'nearest', 'extrapolate': True}
        return clean

    caltab = SimpleNamespace(tarr=original_tarr[[2, 0, 1]].copy(),
                             applycal=Mock(side_effect=applycal))

    def execute(path, run_name):
        assert run_name == '__main__'
        assert Path(path) == output / 'official-I-direct.py'
        adapted = Path(path).read_text()
        assert adapted == source.replace("ttype     = 'nfft'", "ttype     = 'direct'", 1)
        assert adapted.count("ttype     = 'direct'") == 1
        assert (output / 'official-I-original.py').read_text() == source
        assert sys.argv == [path, '-i', str(uv.resolve()),
                            '-o', str(output / 'official-I-restored.fits')]
        assert iers.conf.auto_download is False
        Path(sys.argv[-1]).write_text('restored I')
        return dict(imgr=imgr, obs_orig=obs_orig, obs_sc=obs_sc, caltab=caltab,
                    fov=256., npix=64)

    # main changes argv itself; register restoration even if an assertion fails.
    monkeypatch.setattr(sys, 'argv', list(sys.argv))
    runner = Mock(side_effect=execute)
    monkeypatch.setattr(wrapper.runpy, 'run_path', runner)
    wrapper.main(arguments(tmp_path))
    runner.assert_called_once_with(str(output / 'official-I-direct.py'), run_name='__main__')
    assert official.read_text() == source
    kind, identity, inputs = staged[0]
    assert len(staged) == 1 and kind == 'research:official-eht-I'
    assert identity['ehtim_version'] == '1.3.2'
    assert inputs == [WRAPPER.resolve(), official.resolve(), uv.resolve(),
                      receipt.resolve(), package_code.resolve()]
    assert caltab.applycal.call_count == 1
    np.testing.assert_array_equal(obs_orig.data, original_data)
    np.testing.assert_array_equal(obs_orig.tarr, original_tarr)
    image.save_fits.assert_called_once_with(str(output / 'official-I-pre-restoration.fits'))
    obs_sc.save_uvfits.assert_called_once_with(str(output / 'official-I-final-selfcal.uvfits'))
    clean.save_uvfits.assert_called_once_with(str(output / 'official-I-gaincal-clean-pol.uvfits'))
    assert (output / 'official-I-gaincal-clean-pol.uvfits').read_text() == 'clean polarization'
    record = json.loads((output / 'polarization-input.json').read_text())
    assert record['source'] == 'obs_orig: before reverse taper / I systematic noise'
    assert record['noise'] == 'original errors propagated by gain application; no I2% additive inflation'
    assert record['gains'] == 'final caltab solved by released I pipeline'
    assert record['remaining_D_sites'] == ['AZ', 'LM', 'PV']
    assert record['frcal_from_release'] is True
    keys = ['vis', 'amp', 'cphase', 'logcamp']
    assert obs_sc.chisq.call_args_list == [call(image, dtype=key, ttype='direct') for key in keys]
    report = json.loads((output / 'report.json').read_text())
    assert report['chi2_final_selfcal'] == dict.fromkeys(keys, 3.)
    assert report['original_fov_uas'] == 128. and report['pixels'] == 64


def test_clean_pol_rebuild_declares_raw_i_and_clean_uv_inputs():
    # Static input-plan contract only: do not import numpy/ehtim/astropy or run P/D.
    tree = ast.parse((ROOT / 'scripts/rebuild_eht_pol_clean_input.py').read_text())
    assignments = {
        node.targets[0].id: node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign) and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
    }
    for name, filename in [('raw', 'official-I-pre-restoration.fits'),
                           ('uv', 'official-I-gaincal-clean-pol.uvfits')]:
        value = assignments[name]
        assert isinstance(value, ast.BinOp) and isinstance(value.op, ast.Div)
        assert ast.unparse(value.left) == 'args.i_root'
        assert ast.literal_eval(value.right) == filename
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)]
    declared = [node for node in calls if isinstance(node.func, ast.Name)
                and node.func.id == 'require_declared_outputs']
    assert len(declared) == 1
    assert [ast.unparse(arg) for arg in declared[0].args] == ['args.i_root', '[raw, uv]']
    inputs = assignments['inputs']
    assert isinstance(inputs, ast.List)
    assert [ast.unparse(node) for node in inputs.elts] == [
        'Path(__file__).resolve()',
        "Path(__file__).with_name('rebuild_eht_fiducial.py').resolve()",
        'manifest_path(args.i_root)', 'raw', 'uv', 'args.feed_file.resolve()',
    ]
