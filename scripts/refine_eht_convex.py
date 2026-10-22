#!/usr/bin/env python3
"""Exact convex P refinement of verified archived fits; one day per new root.

No blur restart, sampling, or science-acceptance threshold. Coupled D uses the
original zero-prior forward fit; each changed D defines a fresh convex problem.
"""
import argparse
from importlib import metadata
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'scripts')]
import reconstruct_eht_official as original  # noqa: E402
from eht_axion.pipeline import common  # noqa: E402
from eht_axion.research import eht_optimizer as solver  # noqa: E402


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('fiducial', 'fixed-d'), required=True)
    p.add_argument('--day', type=int, choices=original.DAYS, required=True)
    for name in ('source-root', 'i-root', 'output-root'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--draw-id', type=int)
    p.add_argument('--outer-rounds', type=int, default=2)
    p.add_argument('--gap-atol', type=float, default=1e-6)
    p.add_argument('--maxiter', type=int, default=20000)
    p.add_argument('--check-every', type=int, default=500)
    p.add_argument('--time-limit', type=float, default=300)
    return p


def build_plan(a):
    import math
    if (min(a.outer_rounds, a.maxiter, a.check_every) < 1
            or not math.isfinite(a.gap_atol) or a.gap_atol < 0
            or not math.isfinite(a.time_limit) or a.time_limit <= 0
            or (a.mode == 'fiducial' and a.draw_id is not None)
            or (a.draw_id is not None and a.draw_id < 0)):
        raise ValueError('invalid numerical budget or draw selection')
    source = common.verify_stage(a.source_root, recursive=True, allow_imported=False)
    im, files = original.stage_inputs(a.i_root, 'I', a.day, original.I_PRODUCTS)
    stage = source['stage']
    allowed = ('research:eht-official:fiducial',) if a.mode == 'fiducial' else (
        'research:eht-official:fiducial', 'research:eht-official:uncertainty',
        'research:eht-fixed-d')
    if stage not in allowed or source['config']['day'] != a.day:
        raise ValueError('source mode/day mismatch')
    c = source['config']['parameters']
    if c != im['config']['parameters']:
        raise ValueError('source and I configurations differ')
    if c['p_data'] != {'pvis': 1, 'm': 1} or c['p_reg'] != {'hw': 100, 'ptv': 1} or c['norm_reg'] is not True:
        raise ValueError('only exact frozen paper weights supported')
    recorded = {Path(r['path']).resolve() for r in common.describe_inputs(
        common.source_closure([common.manifest_path(a.source_root)]))}
    if any((a.i_root/n).resolve() not in recorded for n in original.I_PRODUCTS):
        raise ValueError('source does not descend from this same I stage')
    draw = a.draw_id
    rawname = 'final-raw-IQU.fits' if draw is None else f'draw-{draw:04d}-raw-IQU.fits'
    dname = 'final-dterms.json' if draw is None else 'all-draws.json'
    consumed = [a.source_root/rawname, a.source_root/dname]
    common.require_declared_outputs(a.source_root, consumed)
    if draw is not None and draw >= len(json.loads(consumed[1].read_text())['draws']):
        raise ValueError('draw outside frozen samples')
    paths = [*files, common.manifest_path(a.source_root), *consumed, Path(__file__),
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
    common.reject_directory_symlinks(destination)
    if destination.exists():
        raise FileExistsError(destination)
    for path in common.source_closure(paths):
        if path == destination or path.is_relative_to(destination) or destination.is_relative_to(path):
            raise ValueError('source/output overlap')
    return dict(day=a.day, mode=a.mode, parameters=c, rawname=rawname, dname=dname,
                inputs=common.describe_inputs(paths),
                cli={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                scope='floating-point fixed-I/D objective certificate; not interval/image uncertainty',
                D_budget='explicit outer budget; joint D optimality not verified')


def run(a, plan, core, eh, out):
    import numpy as np
    c = plan['parameters']
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
    if a.draw_id is not None:
        data = {s: [v['R'], v['L']] for s, v in data['draws'][a.draw_id].items()}
    present = set(map(str, base.tarr['site']))
    dterms = original.d_dictionary({s: v for s, v in data.items() if s in present})
    steps = []
    def problem(d):
        obs = core.make_imaging_obs(core.apply_fixed_d(core.add_fsys(base.copy(), c['fsys']), d))
        return solver.exact_cartesian_problem(fixed, obs, eh)
    count = a.outer_rounds+1 if a.mode == 'fiducial' else 1
    for r in range(1, count+1):
        label = 'final-P' if r == count else f'round-{r:02d}'
        before = image.qvec+1j*image.uvec
        # New problem, B/b, steps, duals and LB history on EVERY call.
        p, report, trace = solver.solve_problem(problem(dterms), before,
            gap_atol=a.gap_atol, maxiter=a.maxiter, check_every=a.check_every,
            time_limit=a.time_limit, callback=lambda row: print(
                f'PDHG day={a.day} call={r} '+json.dumps(row, allow_nan=False), flush=True))
        image = fixed.copy()
        image.add_qu(p.real.reshape(fixed.ydim, fixed.xdim), p.imag.reshape(fixed.ydim, fixed.xdim))
        image.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
        row = dict(call=r, P=report, trace=trace,
                   raw_QU_delta_L2=float(np.linalg.norm(p-before)),
                   statistics=original.image_statistics(image))
        if a.mode == 'fiducial' and r < count:
            solved = core.fit_forward_d(core.add_fsys(base.copy(), c['fsys']), image, sites=c['solve_sites'])
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
    original.save_image(image, out/'final', c, eh)
    original.dump(out/'final-dterms.json', dterms)
    # Small archive: no matrices, duals, per-iteration images, or duplicate I.
    np.savez_compressed(out/'best-P.npz', P=p)
    original.dump(out/'termination.json', dict(
        all_P_numerically_certified=all(s['P']['numerically_certified'] for s in steps),
        joint_D_optimality='not verified' if a.mode == 'fiducial' else 'not applicable: frozen D',
        D_refitted=a.mode == 'fiducial', outer_rounds=a.outer_rounds if a.mode == 'fiducial' else 0,
        reason='explicit outer budget; final P aligned with last D', blur_restarts=False))
    return steps


def main(argv=None):
    a = parser().parse_args(argv)
    with common.verification_session():
        plan = build_plan(a)
        import ehtim as eh
        from astropy.utils import iers
        from eht_axion.research import eht_official as core
        iers.conf.auto_download = False
        declared = {r['path'] for r in plan['inputs']}
        for name in ('ehtim', 'numpy', 'scipy', 'astropy'):
            module = sys.modules[name]
            package_sources = Path(module.__file__).resolve().parent.rglob('*.py')
            if any(str(p.resolve()) not in declared for p in package_sources):
                raise ValueError('actual imported package outside declared closure')
        with common.stage_directory(a.output_root, 'research:eht-convex:'+a.mode, plan,
                                    [Path(r['path']) for r in plan['inputs']]) as out:
            if out is None:
                raise FileExistsError(a.output_root)
            original.dump(out/'plan.json', plan)
            print(f'READY day={a.day} working={out}', flush=True)
            return run(a, plan, core, eh, out)


if __name__ == '__main__':
    main()
