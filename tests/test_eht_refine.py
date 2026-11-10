"""scripts/eht/refine.py: settings from experiment.yaml and the scoped D-optimizer options (synthetic)."""
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from scripts.eht import refine as w
from eht_axion.experiment import read_experiment
from eht_axion.eht import calibration as core


def test_refine_settings_of_experiment_yaml():
    refine = read_experiment()['eht_imaging']['refine']
    assert refine == dict(outer_rounds=100, gap_atol=1e-6, maxiter=20000, check_every=500, time_limit=300.,
                          d_optimizer=dict(ftol=1e-13, gtol=1e-9, maxiter=10000, maxcor=50, maxls=50,
                                           maxfun=100000))


def test_command_line_names_only_roots():
    options = {a.dest for a in w.parser()._actions} - {'help'}
    assert options == {'config', 'day', 'source_root', 'i_root', 'output_root', 'data_root', 'env_file',
                       'no_env_file'}
    with pytest.raises(SystemExit):
        w.parser().parse_args(['--day', '95', '--source-root', 'a', '--i-root', 'b', '--output-root', 'c',
                               '--outer-rounds', '2'])


def test_recording_optimizer_merges_options_without_touching_scipy(monkeypatch):
    import scipy.optimize

    result = SimpleNamespace(fun=1., jac=np.array([1e-7]), nit=8, nfev=10, njev=10, success=True,
                             status=0, message='ftol', x=np.zeros(1))
    module = SimpleNamespace(minimize=Mock(return_value=result), OptimizeResult=object)
    records = []
    proxy = core._RecordingOptimizer(module, dict(ftol=1e-13, maxiter=10000), records)
    assert proxy.minimize(lambda x: 2., np.zeros(1), method='L-BFGS-B',
                          options={'maxiter': 10, 'gtol': 1e-5}) is result
    options = module.minimize.call_args.kwargs['options']
    assert options == dict(maxiter=10000, gtol=1e-5, ftol=1e-13)
    assert records[0]['call'] == 1 and records[0]['initial_objective'] == 2.
    assert records[0]['options'] == options and records[0]['g_inf'] == 1e-7
    assert proxy.OptimizeResult is object  # everything else is the real module
    assert scipy.optimize.minimize.__module__.startswith('scipy')


def test_private_solver_namespace_uses_the_proxy_only():
    pytest.importorskip('ehtim')
    from ehtim.calibrating import pol_cal_new
    sentinel = object()
    function, receipt = core._forward_solver(sentinel)
    assert function.__globals__['opt'] is sentinel
    assert pol_cal_new.opt is not sentinel
    assert core._forward_solver()[0].__globals__['opt'] is pol_cal_new.opt
    assert receipt['prior'] == 'explicit constant zero'
