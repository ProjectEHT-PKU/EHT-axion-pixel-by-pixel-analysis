"""Verify the D-only options hook does not alter the P solver or model."""
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np

import refine_eht_d_tolerance as w


def test_tight_options_are_scoped_to_d(monkeypatch, tmp_path):
    import scipy.optimize as opt
    from eht_axion.research import eht_official as core
    from ehtim.calibrating import pol_cal_new
    actual = Mock(return_value=SimpleNamespace(fun=1., jac=np.array([1e-7]), nit=8,
        nfev=10, njev=10, success=True, status=0, message='ftol', x=np.zeros(1)))
    monkeypatch.setattr(opt, 'minimize', actual)
    def fit(base, image, sites):
        assert base == 'base' and image == 'P' and sites == ['AZ']
        return pol_cal_new.opt.minimize(lambda x: 2., np.zeros(1), options={'ftol': 1e-6})
    monkeypatch.setattr(core, 'fit_forward_d', fit)
    monkeypatch.setattr(w.runner, 'build_plan', lambda *a: {'inputs': []})
    monkeypatch.setattr(w.runner.original, 'dump', Mock())
    def main(argv):
        assert argv == ['synthetic']
        core.fit_forward_d('base', 'P', ['AZ'])
        options = actual.call_args.kwargs['options']
        assert options['ftol'] == 1e-13 and options['gtol'] == 1e-9
        assert options['maxiter'] == 10000
        # The process-wide scipy object was restored after the scoped D call.
        assert opt.minimize is actual
        w.runner.original.dump(tmp_path/'steps.json', [])
        assert w.runner.original.dump is not None
        return 'done'
    monkeypatch.setattr(w.runner, 'main', main)
    assert w.main(['synthetic']) == 'done'
    assert core.fit_forward_d is fit
