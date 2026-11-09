#!/usr/bin/env python3
"""EHT imaging: exact convex P refinement and D refits of a fiducial stage; one day per new root.

Produces the refined ("tightD") fiducial D terms used for the D-term Monte Carlo. Each
round solves the fixed-I/fixed-D P problem to a certified duality gap, then refits
D (official zero-prior forward fit) against the same P; each changed D defines a
fresh convex problem. All numerical settings (rounds, P solver budget and
tolerance, D optimizer options) come from experiment.yaml eht_imaging.refine; the
command line only names the input and output roots. No blur restart, sampling or
science-acceptance threshold.
"""
import argparse
from importlib import metadata
import json
from pathlib import Path
import sys

from _bootstrap import ROOT
import eht_reconstruct as original
from eht_axion import experiment
from eht_axion import provenance
from eht_axion.eht import convex as solver


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=experiment.DEFAULT_PATH, help='default: experiment.yaml')
    p.add_argument('--day', type=int, choices=experiment.EHT_DAYS, required=True)
    p.add_argument('--source-root', type=Path, required=True, help='completed fiducial stage of this day')
    p.add_argument('--i-root', type=Path, required=True, help='completed I stage of this day')
    p.add_argument('--output-root', type=Path, required=True, help='new directory')
    return p


def build_plan(a):
    config = experiment.read_experiment(a.config)
    refine = config['eht_imaging']['refine']
    source = provenance.verify_stage(a.source_root, recursive=True)
    im, files = original.stage_inputs(a.i_root, 'I', a.day, original.I_PRODUCTS)
    if source['stage'] != 'eht:fiducial' or source['config']['day'] != a.day:
        raise ValueError('source must be the fiducial stage of the same day')
    c = source['config']['parameters']
    if c != im['config']['parameters']:
        raise ValueError('source and I configurations differ')
    locations = {key: c[key] for key in experiment.EHT_LOCATIONS}
    if c != json.loads(json.dumps(experiment.eht_imaging_parameters(config, locations))):
        raise ValueError('source imaging parameters differ from experiment.yaml eht_imaging')
    recorded = {Path(r['path']).resolve() for r in provenance.describe_inputs(
        provenance.source_closure([provenance.manifest_path(a.source_root)]))}
    if any((a.i_root/n).resolve() not in recorded for n in original.I_PRODUCTS):
        raise ValueError('source does not descend from this same I stage')
    rawname, dname = 'final-raw-IQU.fits', 'final-dterms.json'
    consumed = [a.source_root/rawname, a.source_root/dname]
    provenance.require_declared_outputs(a.source_root, consumed)
    paths = [*files, provenance.manifest_path(a.source_root), *consumed, Path(__file__), a.config,
             Path(original.__file__), ROOT/'pixi.lock', Path(c['feed_rotation_file']),
             *sorted((ROOT/'eht_axion').rglob('*.py'))]
    # Full Python/package metadata closure, without hashing every native .so.
    for name in ('ehtim', 'numpy', 'scipy', 'astropy'):
        dist = metadata.distribution(name)
        if name == 'ehtim' and dist.version != '1.3.2':
            raise ValueError('requires ehtim 1.3.2')
        paths.extend(Path(dist.locate_file(f)) for f in dist.files or []
                     if str(f).endswith('.py') or Path(str(f)).name in
                     ('METADATA', 'RECORD', 'direct_url.json'))
        # Conda can install additional Python files absent from wheel RECORD.
        paths.extend(Path(dist.locate_file(name)).rglob('*.py'))
    paths = sorted(set(p.resolve(strict=True) for p in paths))
    destination = a.output_root.absolute()
    provenance.reject_directory_symlinks(destination)
    if destination.exists():
        raise FileExistsError(destination)
    for path in provenance.source_closure(paths):
        if path == destination or path.is_relative_to(destination) or destination.is_relative_to(path):
            raise ValueError('source/output overlap')
    return dict(day=a.day, parameters=c, refine=refine, rawname=rawname, dname=dname,
                inputs=provenance.describe_inputs(paths),
                cli={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                scope='floating-point fixed-I/D objective certificate; not interval/image uncertainty',
                D_budget='explicit outer budget; joint D optimality not verified',
                D_goal='numeric stationarity diagnostic only, not D uncertainty or global joint optimum')


def run(a, plan, core, eh, out):
    import numpy as np
    c, refine = plan['parameters'], plan['refine']
    fixed = eh.image.load_fits(str(a.i_root/original.I_PRODUCTS[0]))
    original.validate_fixed(fixed, c, eh)
    image = eh.image.load_fits(str(a.source_root/plan['rawname']))
    np.testing.assert_array_equal(image.imvec, fixed.imvec)
    if ((image.xdim, image.ydim, image.psize, image.rf) !=
            (fixed.xdim, fixed.ydim, fixed.psize, fixed.rf) or np.any(image.vvec != 0)):
        raise ValueError('archived image domain/V differs')
    base = original.load_gaincal_base(a.i_root/original.I_PRODUCTS[1], c, eh)
    base = base.flag_uvdist(uv_min=c['uv_min'])
    data = json.loads((a.source_root/plan['dname']).read_text())
    present = set(map(str, base.tarr['site']))
    dterms = original.d_dictionary({s: v for s, v in data.items() if s in present})
    steps, d_records = [], []

    def problem(d):
        obs = core.make_imaging_obs(core.apply_fixed_d(core.add_fsys(base.copy(), c['fsys']), d))
        return solver.exact_cartesian_problem(fixed, obs, library=eh, p_data=c['p_data'],
                                              p_reg=c['p_reg'], norm_reg=c['norm_reg'])

    count = refine['outer_rounds'] + 1
    for r in range(1, count+1):
        label = 'final-P' if r == count else f'round-{r:02d}'
        before = image.qvec+1j*image.uvec
        # New problem, B/b, steps, duals and LB history on EVERY call.
        p, report, trace = solver.solve_problem(problem(dterms), before,
            gap_atol=refine['gap_atol'], maxiter=refine['maxiter'], check_every=refine['check_every'],
            time_limit=refine['time_limit'], callback=lambda row: print(
                f'PDHG day={a.day} call={r} '+json.dumps(row, allow_nan=False), flush=True))
        image = fixed.copy()
        image.add_qu(p.real.reshape(fixed.ydim, fixed.xdim), p.imag.reshape(fixed.ydim, fixed.xdim))
        image.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
        row = dict(call=r, P=report, trace=trace,
                   raw_QU_delta_L2=float(np.linalg.norm(p-before)))
        if r < count:
            solved = core.fit_forward_d(core.add_fsys(base.copy(), c['fsys']), image, sites=c['solve_sites'],
                                        optimizer_options=refine['d_optimizer'], records=d_records)
            new = original.d_dictionary({str(t['site']): [t['dr'], t['dl']]
                                        for t in solved.tarr if str(t['site']) in c['solve_sites']})
            # SAME physical P endpoint, different declared fixed-D objectives.
            fresh = problem(new)
            after = solver.Split(fresh).cert(p/fresh.I)
            row.update(before_D_same_P=dict(primal=report['best_primal'],
                       lower=report['best_lower'], gap=report['best_gap']),
                       after_D_same_P=after,
                       D_max_delta=max(float(np.max(abs(np.asarray(new[s])-np.asarray(dterms[s]))))
                                       for s in new))
            dterms = new
        original.save_image(image, out/label, c, eh)
        original.dump(out/(label+'-dterms.json'), dterms)
        steps.append(row)
        original.dump(out/'steps.json', steps)
        original.dump(out/'D-optimizer-records.json', d_records)
    original.save_image(image, out/'final', c, eh)
    original.dump(out/'final-dterms.json', dterms)
    # Small archive: no matrices, duals, per-iteration images, or duplicate I.
    np.savez_compressed(out/'best-P.npz', P=p)
    original.dump(out/'termination.json', dict(
        all_P_numerically_certified=all(s['P']['numerically_certified'] for s in steps),
        joint_D_optimality='not verified', D_refitted=True, outer_rounds=refine['outer_rounds'],
        reason='explicit outer budget; final P aligned with last D', blur_restarts=False))
    return steps


def main(argv=None):
    a = parser().parse_args(argv)
    with provenance.verification_session():
        plan = build_plan(a)
        import ehtim as eh
        from astropy.utils import iers
        from eht_axion.eht import calibration as core
        iers.conf.auto_download = False
        declared = {r['path'] for r in plan['inputs']}
        for name in ('ehtim', 'numpy', 'scipy', 'astropy'):
            module = sys.modules[name]
            package_sources = Path(module.__file__).resolve().parent.rglob('*.py')
            if any(str(p.resolve()) not in declared for p in package_sources):
                raise ValueError('actual imported package outside declared closure')
        with provenance.stage_directory(a.output_root, 'eht:refine', plan,
                                    [Path(r['path']) for r in plan['inputs']]) as out:
            if out is None:
                raise FileExistsError(a.output_root)
            original.dump(out/'plan.json', plan)
            print(f'READY day={a.day} working={out}', flush=True)
            return run(a, plan, core, eh, out)


if __name__ == '__main__':
    main()
