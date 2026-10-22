#!/usr/bin/env python3
"""Fixed-D runner with explicit absent-station handling and P-only polish.

Reuses the frozen reconstruction entrypoint without editing its executed bytes.
Only uncertainty is supported. Sampler retains all seven stations, but a day's
Jones correction uses only stations present that day; unused draws are recorded.
Five 20uas restart/optimization passes, then an unblurred final polish. No D fit.
A pilot is not an automatically certified full ensemble.
"""
from pathlib import Path
import sys

import reconstruct_eht_official as driver

PASSES = 5


def fixed_d_image(base, fixed, dterms, c, core, eh):
    present = set(map(str, base.tarr['site']))
    available = {s: pair for s, pair in dterms.items() if s in present}
    if set(available) != present:
        raise ValueError('sampled D set does not cover every station in this observation')
    noisy = core.add_fsys(base.copy(), c['fsys'])
    corrected = core.make_imaging_obs(core.apply_fixed_d(noisy, available))
    image, records = driver.initial_image(fixed, c), []
    # These are numerical optimization passes on ONE fixed objective and D set.
    # No measurement redrawing, objective weight change, or leakage fit occurs.
    short = dict(c, maxit=100)
    for k in range(PASSES):
        init = image if k == 0 else driver.restart(image, fixed, c, eh)
        image, actual = driver.image_p(corrected, fixed, init, short, eh)
        records += [dict(r, p_only_pass=k + 1, blur_restart=k > 0,
                         unused_station_draws=sorted(set(dterms) - present)) for r in actual]
    image, actual = driver.image_p(corrected, fixed, image, c, eh)
    records += [dict(r, p_only_pass=PASSES + 1, blur_restart=False,
                     unused_station_draws=sorted(set(dterms) - present)) for r in actual]
    return image, records


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    args = driver.parser().parse_args(argv)
    if args.stage != 'uncertainty':
        raise ValueError('This entrypoint supports uncertainty only')
    original = driver.build_plan

    def build_plan(request):
        plan = original(request)
        paths = [Path(r['path']) for r in plan['inputs']] + [Path(__file__).resolve()]
        plan['inputs'] = driver.common.describe_inputs(paths)
        plan['fixed_D_schedule'] = dict(p_only_passes=PASSES, per_pass_maxit=100,
            blur_restart_uas=plan['parameters']['restart_blur_uas'],
            final_unblurred_polish_maxit=plan['parameters']['maxit'],
            fit_D=False, absent_stations='ignore only unused draws, record each optimizer row',
            validated_full_ensemble=False)
        return plan

    driver.build_plan, driver.fixed_d_image = build_plan, fixed_d_image
    try:
        return driver.main(argv)
    finally:
        driver.build_plan = original


if __name__ == '__main__':
    main()
