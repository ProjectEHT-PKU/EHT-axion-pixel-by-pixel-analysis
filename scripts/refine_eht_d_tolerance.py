#!/usr/bin/env python3
"""Same-model numerical D-tolerance experiment using certified P subproblems.

Frozen previous drivers are not edited. Every D fit's actual options, objective,
gradient and initialization are recorded. No D prior or data term is added.
"""
from pathlib import Path
import sys
import time
from unittest.mock import patch

import refine_eht_convex as runner


def main(argv=None):
    import numpy as np
    import scipy.optimize as opt
    from eht_axion.research import eht_official as core
    from eht_axion.pipeline import common
    from ehtim.calibrating import pol_cal_new
    argv = sys.argv[1:] if argv is None else list(argv)
    original_plan = runner.build_plan
    original_fit = core.fit_forward_d
    records = []
    numerical = dict(ftol=1e-13, gtol=1e-9, maxiter=10000, maxcor=50, maxls=50, maxfun=100000)

    def build_plan(args):
        p = original_plan(args)
        p['inputs'] = common.describe_inputs([Path(r['path']) for r in p['inputs']]
                                             + [Path(__file__).resolve()])
        p['D_numerical_options'] = numerical
        p['D_goal'] = 'numeric stationarity diagnostic only, not D uncertainty or global joint optimum'
        return p

    def fit(base, image, sites):
        actual = opt.minimize
        started = time.monotonic()

        def minimize(fun, x0, *args, **kwargs):
            options = dict(kwargs.get('options', {}))
            options.update(numerical)
            kwargs['options'] = options
            initial = float(fun(x0))
            result = actual(fun, x0, *args, **kwargs)
            record = dict(call=len(records)+1, options=options, initial_objective=initial,
                initial_coordinates=np.asarray(x0).tolist(), objective=float(result.fun),
                g_inf=float(np.max(abs(result.jac))), nit=int(result.nit), nfev=int(result.nfev),
                njev=int(result.njev), success=bool(result.success), status=int(result.status),
                message=str(result.message), seconds=time.monotonic()-started,
                final_coordinates=np.asarray(result.x).tolist())
            records.append(record)
            print('D_OPTIMIZER', record, flush=True)
            return result

        # Narrow process-local scope: P optimizer is never patched by this hook.
        with patch.object(pol_cal_new.opt, 'minimize', minimize):
            return original_fit(base, image, sites=sites)

    old_dump = runner.original.dump

    def dump(path, value):
        old_dump(path, value)
        if Path(path).name == 'steps.json':
            old_dump(Path(path).parent/'D-optimizer-records.json', records)

    runner.build_plan = build_plan
    try:
        with patch.object(core, 'fit_forward_d', fit), patch.object(runner.original, 'dump', dump):
            return runner.main(argv)
    finally:
        runner.build_plan = original_plan


if __name__ == '__main__':
    main()
