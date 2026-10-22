"""Synthetic contracts for the continuation entry point."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('refine_test', ROOT/'scripts/refine_eht_optimization.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


def args(tmp_path):
    return w.parser().parse_args(['--mode', 'fixed-d', '--day', '100', '--draw-id', '0',
        '--source-root', str(tmp_path/'source'), '--i-root', str(tmp_path/'I'),
        '--output-root', str(tmp_path/'new')])


@pytest.mark.parametrize('key,value', [('p_maxit', 0), ('outer_rounds', -1),
    ('checkpoint_every', 0), ('p_ftol', float('nan')), ('p_gtol', 0), ('draw_id', -1)])
def test_guards(tmp_path, key, value):
    a = args(tmp_path)
    setattr(a, key, value)
    with pytest.raises(ValueError):
        w.validate_args(a)


def test_defaults(tmp_path):
    a = args(tmp_path)
    w.validate_args(a)
    assert (a.p_maxit, a.p_ftol, a.p_gtol, a.outer_rounds) == (10000, 1e-12, 1e-8, 10)


def test_shift_only_objective(tmp_path):
    class Imager:
        def __init__(self, *ar, **kw):
            self.kw = kw
        def objfunc(self, x):
            return np.dot(x, x)
        def objgrad(self, x):
            return 2*x
    eh = SimpleNamespace(imager=SimpleNamespace(Imager=Imager))
    c = dict(p_data={'pvis': 1., 'm': 1.}, p_reg={'hw': 100., 'ptv': 1.}, norm_reg=True)
    image = SimpleNamespace(total_flux=lambda: .6)
    imgr = w.make_imager(None, image, image, c, eh, args(tmp_path))
    x = np.array([1., 2.])
    assert imgr.objfunc(x) == pytest.approx(5 + 2 + 100*np.log(2))
    np.testing.assert_array_equal(imgr.objgrad(x), 2*x)
    assert imgr.kw['transform'] == ['mcv']
    assert imgr.kw['data_term'] is c['p_data']
    assert imgr.kw['reg_term'] is c['p_reg']


def test_noising_fresh_once_before_D():
    from copy import deepcopy
    class Base:
        noisy = False
        def copy(self):
            return deepcopy(self)
    base, events = Base(), []
    def noise(b, fsys):
        assert not b.noisy
        b.noisy = True
        events.append('noise')
        return b
    def apply(b, d):
        assert b.noisy
        events.append('D')
        return b
    core = SimpleNamespace(add_fsys=noise, apply_fixed_d=apply,
        make_imaging_obs=lambda b: events.append('adapter') or b)
    for _ in range(3):
        w.corrected(base, {}, {'fsys': .005}, core)
    assert not base.noisy
    assert events == ['noise', 'D', 'adapter']*3


def test_actual_imports():
    assert Path(w.original.__file__).resolve() == ROOT/'scripts/reconstruct_eht_official.py'
    assert Path(w.common.__file__).resolve() == ROOT/'eht_axion/pipeline/common.py'
