#!/usr/bin/env python3
"""Unblurred, fixed-I continuation of completed official P/D or frozen-D fits.

One day per invocation. No convergence verdict is inferred from solver success.
Callbacks checkpoint a continuous L-BFGS-B solve; outer calls lose its history.
"""
from __future__ import annotations

import argparse
from importlib import metadata
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '4'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
import reconstruct_eht_official as original  # noqa: E402
from eht_axion.pipeline import common  # noqa: E402


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('fiducial', 'fixed-d'), required=True)
    p.add_argument('--day', '--days', dest='day', type=int, choices=original.DAYS, required=True,
                   help='single day (95, 96, 100, 101); independent new root per invocation')
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--i-root', type=Path, required=True)
    p.add_argument('--draw-id', type=int)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--p-maxit', type=int, default=10000)
    p.add_argument('--p-ftol', type=float, default=1e-12)
    p.add_argument('--p-gtol', type=float, default=1e-8)
    p.add_argument('--outer-rounds', type=int, default=10)
    p.add_argument('--checkpoint-every', '--chunk', dest='checkpoint_every', type=int, default=50,
                   help='callback interval, NOT a solver restart interval')
    p.add_argument('--plan', action='store_true')
    return p


def validate_args(a):
    import math
    for key in ('p_maxit', 'outer_rounds', 'checkpoint_every'):
        if getattr(a, key) < 1:
            raise ValueError(f'{key} must be positive')
    for key in ('p_ftol', 'p_gtol'):
        if not math.isfinite(getattr(a, key)) or not 0 < getattr(a, key) < 1:
            raise ValueError(f'{key} must be finite and in (0,1)')
    if (a.mode == 'fixed-d') != (a.draw_id is not None):
        raise ValueError('--draw-id required exactly for fixed-d')
    if a.draw_id is not None and a.draw_id < 0:
        raise ValueError('negative draw-id')


def build_plan(a):
    validate_args(a)
    source = common.verify_stage(a.source_root, recursive=True, allow_imported=False)
    im, files = original.stage_inputs(a.i_root, 'I', a.day, original.I_PRODUCTS)
    expected = 'research:eht-official:' + ('fiducial' if a.mode == 'fiducial' else 'uncertainty')
    if source['stage'] != expected or source['config']['day'] != a.day:
        raise ValueError('source mode/day mismatch')
    c = source['config']['parameters']
    if c != im['config']['parameters']:
        raise ValueError('source and I configurations differ')
    # Match the actual fixed-I bytes, not just the grid/configuration.
    recorded = {Path(r['path']).resolve(): r for r in common.describe_inputs(
        common.source_closure([common.manifest_path(a.source_root)]))}
    for name in original.I_PRODUCTS:
        path = (a.i_root / name).resolve()
        if path not in recorded:
            raise ValueError('source does not descend from this same I-stage product')
    rawname = 'final-raw-IQU.fits' if a.mode == 'fiducial' else f'draw-{a.draw_id:04d}-raw-IQU.fits'
    names = [rawname, 'final-dterms.json' if a.mode == 'fiducial' else 'all-draws.json']
    common.require_declared_outputs(a.source_root, [a.source_root / n for n in names])
    if a.mode == 'fixed-d':
        samples = json.loads((a.source_root / 'all-draws.json').read_text())
        if a.draw_id >= len(samples['draws']):
            raise ValueError('draw-id outside frozen samples')
    paths = [*files, common.manifest_path(a.source_root), *[a.source_root / n for n in names],
             Path(__file__), Path(original.__file__), ROOT / 'pixi.lock',
             *sorted((ROOT / 'eht_axion').rglob('*.py'))]
    # Hash runtime package sources and native extensions, including local patches.
    for package in ('ehtim', 'scipy', 'numpy', 'astropy'):
        dist = metadata.distribution(package)
        if package == 'ehtim' and dist.version != '1.3.2':
            raise ValueError('requires ehtim==1.3.2')
        paths.extend(Path(dist.locate_file(f)) for f in dist.files or []
                     if str(f).endswith(('.py', '.so')) or Path(str(f)).name in ('METADATA', 'RECORD'))
    paths = sorted(set(p.resolve(strict=True) for p in paths))
    out = a.output_root.absolute()
    common.reject_directory_symlinks(out)
    if out.exists():
        raise FileExistsError(out)
    for path in common.source_closure(paths):
        if path == out or path.is_relative_to(out) or out.is_relative_to(path):
            raise ValueError(f'source/output overlap: {path}')
    return dict(day=a.day, mode=a.mode, parameters=c, rawname=rawname,
                inputs=common.describe_inputs(paths),
                cli={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                scientifically_converged=False, blur_restarts=False,
                history='continuous within P call; reset after each outer/D update')


def assert_image(image, fixed):
    import numpy as np
    np.testing.assert_array_equal(image.imvec, fixed.imvec)
    p = image.qvec + 1j * image.uvec
    if not np.isfinite(p).all() or np.any(abs(p) >= fixed.imvec) or np.any(image.vvec != 0):
        raise ValueError('nonfinite/unphysical P or V != 0')
    if (image.xdim, image.ydim, image.psize, image.rf) != (fixed.xdim, fixed.ydim, fixed.psize, fixed.rf):
        raise ValueError('image domain differs from fixed I')


def make_imager(obs, fixed, init, c, eh, a):
    import numpy as np
    shift = sum(c['p_data'].values()) + c['p_reg']['hw'] * np.log(2.)

    class ContinuationImager(eh.imager.Imager):
        def objfunc(self, x):
            return super().objfunc(x) + shift

    return ContinuationImager(obs, init, prior_im=init, flux=float(fixed.total_flux()),
                             pol='P', transform=['mcv'], data_term=c['p_data'],
                             reg_term=c['p_reg'], norm_reg=c['norm_reg'], norm_init=False,
                             maxit=a.p_maxit, stop=a.p_ftol, ttype='direct')


def physical(imgr, x, eh):
    return eh.imager.transform_imarr(
        eh.imager.unpack_imarr(x, imgr._xarr, imgr._which_solve),
        imgr.transform_next, imgr._which_solve)


def evaluate(imgr, x, eh):
    import numpy as np
    j, g = float(imgr.objfunc(x)), np.asarray(imgr.objgrad(x))
    if not np.isfinite(j) or not np.isfinite(g).all():
        raise ValueError('nonfinite objective/actual analytic gradient')
    chi = imgr.make_chisq_dict(physical(imgr, x, eh))
    return dict(objective=j, g_inf=float(np.max(abs(g))), g_l2=float(np.linalg.norm(g)),
                data_chi2={k: float(v) for k, v in chi.items()})


def image_p(obs, fixed, init, c, eh, a, prefix):
    import numpy as np
    import scipy.optimize as opt
    imgr = make_imager(obs, fixed, init, c, eh, a)
    records, minimize = [], opt.minimize

    def intercept(fun, x0, *args, **kwargs):
        if kwargs.get('method') != 'L-BFGS-B' or not callable(kwargs.get('jac')):
            raise ValueError('expected analytic L-BFGS-B')
        options = dict(kwargs.get('options', {}), maxiter=a.p_maxit, ftol=a.p_ftol,
                       gtol=a.p_gtol, maxls=50, maxcor=50, maxfun=a.p_maxit * 10 + 100)
        kwargs['options'] = options
        old_callback = kwargs.get('callback')
        initial = evaluate(imgr, x0, eh)
        previous, count = np.array(x0, copy=True), 0
        trace = []

        def checkpoint(x, final=False):
            nonlocal previous
            row = dict(iteration=count, **evaluate(imgr, x, eh),
                       stepnorm=float(np.linalg.norm(x - previous)),
                       step_scope='since previous recorded checkpoint', final=final)
            row['objective_minus_initial'] = row['objective'] - initial['objective']
            trace.append(row)
            np.savez(str(prefix) + '-checkpoint.npz', x=x, physical=physical(imgr, x, eh),
                     iteration=count)
            original.dump(str(prefix) + '-progress.json', dict(options=options, initial=initial, trace=trace))
            previous = np.array(x, copy=True)

        def callback(x):
            nonlocal count
            count += 1
            if old_callback:
                old_callback(x)
            if count % a.checkpoint_every == 0:
                checkpoint(x)

        kwargs['callback'] = callback
        result = minimize(fun, x0, *args, **kwargs)
        checkpoint(result.x, final=True)
        records.append(dict(options=options, initial=initial, final=trace[-1], trace=trace,
                            nit=int(result.nit), nfev=int(result.nfev), njev=int(result.njev),
                            success=bool(result.success), status=int(result.status),
                            message=str(result.message), history_reset_at_call_start=True))
        return result

    with patch.object(opt, 'minimize', intercept):
        image = imgr.make_image_P(show_updates=False)
    assert_image(image, fixed)
    return image, records


def corrected(base, dterms, c, core):
    return core.make_imaging_obs(core.apply_fixed_d(core.add_fsys(base.copy(), c['fsys']), dterms))


def compare_images(before, after, c, eh):
    import numpy as np
    p0, p1 = before.qvec + 1j * before.uvec, after.qvec + 1j * after.uvec
    result = dict(raw_relative_QU_L2=float(np.linalg.norm(p1-p0) / np.linalg.norm(p0)))
    b = c['display_blur_uas'] * eh.RADPERUAS
    x, y = before.blur_circ(b, b), after.blur_circ(b, b)
    p, q = x.qvec + 1j*x.uvec, y.qvec + 1j*y.uvec
    mask = ((x.imvec > .1*x.imvec.max()) & (y.imvec > .1*y.imvec.max())
            & (abs(p) > .2*abs(p).max()) & (abs(q) > .2*abs(q).max()))
    w = np.sqrt(abs(p[mask])*abs(q[mask]))
    result.update(bright_mask_pixels=int(mask.sum()),
                  bright_mask='intersection: additional19uas I>0.1peakI, P>0.2peakP; diagnostic only',
                  evpa_rms_rad=float(np.sqrt(np.sum(w*(.5*np.angle(q[mask]*p[mask].conj()))**2)/w.sum()))
                  if w.sum() else None,
                  before=original.image_statistics(x), after=original.image_statistics(y))
    return result


def run(a, plan, core, eh, out):
    import numpy as np
    c = plan['parameters']
    fixed = eh.image.load_fits(str(a.i_root / original.I_PRODUCTS[0]))
    original.validate_fixed(fixed, c, eh)
    image = eh.image.load_fits(str(a.source_root / plan['rawname']))
    assert_image(image, fixed)
    base = original.load_gaincal_base(a.i_root / original.I_PRODUCTS[1], c, eh)
    base = base.flag_uvdist(uv_min=c['uv_min'])
    if a.mode == 'fixed-d':
        draw = json.loads((a.source_root / 'all-draws.json').read_text())['draws'][a.draw_id]
        present = set(map(str, base.tarr['site']))
        dterms = original.d_dictionary({s: [v['R'], v['L']] for s, v in draw.items() if s in present})
    else:
        dterms = original.d_dictionary(json.loads((a.source_root / 'final-dterms.json').read_text()))
    steps = []
    for r in range(1, a.outer_rounds + 1):
        prefix = out / f'round-{r:02d}'
        before = image.copy()
        image, optimizers = image_p(corrected(base, dterms, c, core), fixed, image, c, eh, a, prefix)
        row = dict(round=r, optimizers=optimizers, **compare_images(before, image, c, eh))
        if a.mode == 'fiducial':
            solved = core.fit_forward_d(core.add_fsys(base.copy(), c['fsys']), image, sites=c['solve_sites'])
            new = original.d_dictionary({str(t['site']): [t['dr'], t['dl']] for t in solved.tarr
                                         if str(t['site']) in c['solve_sites']})
            row['D_max_delta_by_station'] = {
                s: float(np.max(np.linalg.norm(np.asarray(new[s])-np.asarray(dterms[s]), axis=1)))
                for s in new}
            # Evaluate both endpoints with the SAME newly corrected D objective.
            obs = corrected(base, new, c, core)
            row['common_new_D_endpoints'] = []
            for endpoint in (before, image):
                probe = make_imager(obs, fixed, endpoint, c, eh, a)
                row['common_new_D_endpoints'].append(evaluate(probe, probe._xinit, eh))
            row['after_D_at_same_P'] = row['common_new_D_endpoints'][1]
            row['before_D_at_same_P'] = optimizers[-1]['final']
            dterms = new
        original.save_image(image, prefix, c, eh)
        original.dump(str(prefix) + '-dterms.json', dterms)
        steps.append(row)
        original.dump(out / 'steps.json', steps)
        print(f'CHECKPOINT day={a.day} round={r} path={prefix} (no convergence verdict)', flush=True)
    # In coupled mode align delivered P with the last D, without another D fit.
    if a.mode == 'fiducial':
        image, optimizers = image_p(corrected(base, dterms, c, core), fixed, image, c, eh, a, out / 'final-P')
        steps.append(dict(stage='final-P-only', optimizers=optimizers))
    original.save_image(image, out / 'final', c, eh)
    original.dump(out / 'final-dterms.json', dterms)
    original.dump(out / 'steps.json', steps)
    original.dump(out / 'termination.json', dict(reason='explicit outer budget',
                  scientifically_converged=False, D_refitted=a.mode == 'fiducial',
                  outer_rounds=a.outer_rounds, blur_restarts=False))
    return steps


def main(argv=None):
    a = parser().parse_args(argv)
    with common.verification_session():
        plan = build_plan(a)
        if a.plan:
            print(json.dumps(plan, indent=2))
            return plan
        import ehtim as eh
        from astropy.utils import iers
        from eht_axion.research import eht_official as core
        iers.conf.auto_download = False
        declared = {r['path'] for r in plan['inputs']}
        for module in (eh, sys.modules['scipy'], sys.modules['numpy']):
            if str(Path(module.__file__).resolve()) not in declared:
                raise ValueError('actual imported library is outside hashed closure')
        with common.stage_directory(
                a.output_root, 'research:eht-optimization:' + a.mode, plan,
                [Path(r['path']) for r in plan['inputs']]) as out:
            if out is None:
                raise FileExistsError(a.output_root)
            original.dump(out / 'plan.json', plan)
            print(f'READY day={a.day} working={out} final={a.output_root.absolute()}', flush=True)
            return run(a, plan, core, eh, out)


if __name__ == '__main__':
    main()
