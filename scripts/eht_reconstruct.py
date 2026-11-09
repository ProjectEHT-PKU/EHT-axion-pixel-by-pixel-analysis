#!/usr/bin/env python3
"""Process 2: Paper IV/VII polarimetric imaging (stages I, fiducial, uncertainty), one day per root.

Scientific parameters: experiment.yaml (eht_imaging + shared.image). Machine-local
locations (uvfits directory, input receipt, feed-rotation table, official script)
come from CLI options or .env (EHT_IMAGING_*). Not a bitwise historical reproduction.

--plan only verifies/hashes inputs; no optional imaging import or output writes.
Execute only after external review. A --review-file records that external decision;
its presence is NOT an automatic approval or scientific-convergence assertion.
Each invocation publishes exactly one transactional stage/day to a new root.
"""
from __future__ import annotations

import argparse
from importlib import metadata
import json
from pathlib import Path
import runpy
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from eht_axion import experiment  # noqa: E402
from eht_axion.pipeline import common  # noqa: E402
from _environment import add_environment_options, resolve_environment  # noqa: E402

DAYS = (95, 96, 100, 101)
I_PRODUCTS = ('final-fixedI.fits', 'clean-gaincal-pol.uvfits')


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=experiment.DEFAULT_PATH, help='default: experiment.yaml')
    p.add_argument('--stage', choices=('I', 'fiducial', 'uncertainty'), required=True)
    p.add_argument('--day', type=int, choices=DAYS, required=True)
    for name, env in LOCATION_ENV.items():
        p.add_argument('--' + name.replace('_', '-'), type=Path, help=f'absolute path; overrides ${env}')
    p.add_argument('--i-root', type=Path, help='completed I stage for the SAME day')
    source = p.add_mutually_exclusive_group()
    source.add_argument('--d-roots', nargs=4, type=Path,
                        help='four completed fiducial stages (all four days)')
    source.add_argument('--published-dterms', type=Path,
                        help='Table5 JSON; published rounded, not private posterior')
    p.add_argument('--output-root', required=True, type=Path, help='new directory for this stage/day')
    p.add_argument('--review-file', type=Path,
                   help='external review bytes, required for execution; not self-approval')
    p.add_argument('--draws', type=int, help='uncertainty pilot prefix of configured full ensemble')
    p.add_argument('--plan', action='store_true')
    add_environment_options(p)
    return p


LOCATION_ENV = dict(input_root='EHT_IMAGING_INPUT_ROOT', input_receipt='EHT_IMAGING_INPUT_RECEIPT',
                    feed_rotation_file='EHT_IMAGING_FEED_ROTATION',
                    official_script='EHT_IMAGING_OFFICIAL_SCRIPT')


def locations(args):
    """CLI > exported environment > .env; every location must be an absolute path."""
    environment = resolve_environment(args)
    result = {}
    for name, env in LOCATION_ENV.items():
        value = getattr(args, name) or (Path(environment[env]) if environment.get(env) else None)
        if value is None or not Path(value).is_absolute():
            raise ValueError(f'{name}: absolute path required (--{name.replace("_", "-")} or {env})')
        result[name] = Path(value)
    return result


def adapt_official(text, fov_uas):
    """Only these two unique exact source substitutions are permitted."""
    for old, new in (("ttype     = 'nfft'", "ttype     = 'direct'"),
                     ('fov       = 128*eh.RADPERUAS', f'fov       = {fov_uas:g}*eh.RADPERUAS')):
        if text.count(old) != 1:
            raise ValueError(f'Official adaptation anchor differs: {old}')
        text = text.replace(old, new, 1)
    return text


def stage_inputs(root, stage, day=None, products=()):
    """Verify the entire source closure AND that each consumed file is declared."""
    m = common.verify_stage(root, recursive=True)
    if m['stage'] != 'research:eht-official:' + stage or (day is not None and m['config']['day'] != day):
        raise ValueError(f'Expected official {stage} stage for day {day}: {root}')
    paths = [root / name for name in products]
    common.require_declared_outputs(root, paths)
    return m, [common.manifest_path(root), *paths]


def build_plan(args):
    config = experiment.read_experiment(args.config)
    c = experiment.eht_imaging_parameters(config, locations(args))
    from eht_axion.research import eht_official as core

    if set(c['fixed_sites']) != set(core._PRECORRECTED):
        raise ValueError(f'fixed_sites must be the HOPS-precorrected stations {core._PRECORRECTED}')
    if (args.stage == 'I' and args.i_root) or (args.stage != 'I' and args.i_root is None):
        raise ValueError('--i-root required exactly for fiducial/uncertainty')
    if args.stage == 'uncertainty':
        if not (args.d_roots or args.published_dterms):
            raise ValueError('uncertainty requires --d-roots or --published-dterms')
    elif args.d_roots or args.published_dterms or args.draws is not None:
        raise ValueError('D sources and --draws apply only to uncertainty')
    draws = args.draws if args.draws is not None else c['uncertainty_draws']
    if not 1 <= draws <= c['uncertainty_draws']:
        raise ValueError(f"--draws must be in [1,{c['uncertainty_draws']}]")
    rawroot = Path(c['input_root'])
    name, expected = c['inputs'][args.day]
    raw = (rawroot / name).resolve(strict=True)
    receipt = json.loads(Path(c['input_receipt']).read_text())
    matches = [r for r in receipt['files'] if Path(r['path']).name == name]
    if (receipt.get('complete') is not True or len(matches) != 1 or matches[0]['sha256'] != expected
            or raw.stat().st_size != matches[0]['size'] or common.sha256(raw) != expected):
        raise ValueError('raw HOPS input checksum/receipt mismatch')
    adapt_official(Path(c['official_script']).read_text(), c['fov_uas'])
    inputs = [args.config.resolve(), raw, Path(c['input_receipt']), Path(c['feed_rotation_file']),
              Path(c['official_script']), *sorted((ROOT / 'eht_axion').rglob('*.py')),
              Path(__file__).resolve(), ROOT / 'pixi.lock']
    if args.i_root:
        m, files = stage_inputs(args.i_root, 'I', args.day, I_PRODUCTS)
        upstream = m['config']['parameters']
        for key in ('npix', 'fov_uas', 'flux_jy'):
            if upstream[key] != c[key]:
                raise ValueError(f'I-stage domain differs: {key}')
        if upstream['feed_rotation_file'] != c['feed_rotation_file']:
            raise ValueError('I-stage feed table differs')
        inputs.extend(files)
    if args.d_roots:
        days = []
        for root in args.d_roots:
            m, files = stage_inputs(root, 'fiducial', products=('final-dterms.json',))
            days.append(m['config']['day'])
            if m['config']['parameters'] != json.loads(json.dumps(c)):
                raise ValueError('all four fiducial D sources must use the same declared configuration')
            inputs.extend(files)
        if sorted(days) != list(DAYS):
            raise ValueError('--d-roots must cover each of the four days exactly once')
    if args.published_dterms:
        published = json.loads(args.published_dterms.read_text())
        if published.get('stored_units') != 'dimensionless' or set(published['days']) != set(map(str, DAYS)):
            raise ValueError('published D input requires four dimensionless days')
        for day in DAYS:
            values = d_dictionary(published['days'][str(day)])
            if set(values) != set(c['solve_sites']):
                raise ValueError('published D requires AZ/LM/PV, both hands, for each day')
        inputs.append(args.published_dterms)
    if args.review_file:
        if not args.review_file.read_bytes().strip():
            raise ValueError('empty external review file')
        inputs.append(args.review_file)
    if not args.plan and args.review_file is None:
        raise ValueError('execution requires --review-file (external decision; no automatic self-approval)')
    # Enumerate optional library source without importing ehtim or executing imaging.
    dist = metadata.distribution('ehtim')
    if dist.version != '1.3.2':
        raise ValueError('requires ehtim==1.3.2')
    package = Path(dist.locate_file('ehtim'))
    inputs.extend(package.rglob('*.py'))
    inputs.extend(Path(dist.locate_file(f)) for f in dist.files or []
                  if Path(str(f)).name in ('METADATA', 'RECORD', 'direct_url.json'))
    inputs = sorted(set(p.resolve(strict=True) for p in inputs))
    destination = args.output_root.absolute()
    common.reject_directory_symlinks(destination)
    if destination.exists():
        raise FileExistsError(f'new stage/day output root required: {destination}')
    for path in common.source_closure(inputs):
        if path == destination or path.is_relative_to(destination) or destination.is_relative_to(path):
            raise ValueError(f'output/source overlap: {path}')
    return dict(day=args.day, stage=args.stage, parameters=c, raw=str(raw), draws=draws,
                label=('engineering-partial-pilot'
                       if args.stage == 'uncertainty' and draws < c['uncertainty_draws'] else 'full-budget'),
                output_root=str(destination), inputs=common.describe_inputs(inputs),
                cli={k: [str(p) for p in v] if isinstance(v, list) else str(v) if isinstance(v, Path) else v
                     for k, v in vars(args).items()},
                scope='Paper IV/VII method implementation; not bitwise historical; review not inferred',
                display=f"additional {c['display_blur_uas']:g}uas I/Q/U kernel, NOT an effective PSF",
                scientifically_converged=False)


def run_i(plan, core, eh, out):
    """Preserve native circular BEFORE the official Stokes-only I workflow."""
    import numpy as np

    c = plan['parameters']
    native, native_meta = core.load_native(plan['raw'], c['feed_rotation_file'])
    text = Path(c['official_script']).read_text()
    (out / 'official-I-original.py').write_text(text)
    executable = out / 'official-I-120-direct.py'
    executable.write_text(adapt_official(text, c['fov_uas']))
    with patch.object(sys, 'argv', [str(executable), '-i', plan['raw'], '-o', str(out / I_PRODUCTS[0])]):
        ns = runpy.run_path(str(executable), run_name='__main__')
    fixed = ns['im_out']  # restored 5uas image, NOT the pre-restoration imgr checkpoint
    validate_fixed(fixed, c, eh)
    if not np.isclose(fixed.rf, native.rf, rtol=0, atol=1.):
        raise ValueError('final fixed I frequency differs from the native observation')
    # The old pipeline excludes PV in its last selfcal. This additional solve is
    # explicit, on untapered scan-averaged data, against the SAME restored fixed I.
    # Exclude extended-flux zero baselines BEFORE solving compact-model gains;
    # the unchanged official I workflow above may still use its zero baselines.
    native = native.flag_uvdist(uv_min=c['uv_min'])
    fitcopy = native.switch_polrep('stokes')
    stations = sorted(set(native.data['t1']) | set(native.data['t2']))
    cal = eh.selfcal(fitcopy, fixed, pol='I', sites=stations, method='both',
                     ttype='direct', gain_tol=[.02, .2], caltable=True, processes=-1)
    for station in stations:
        gains = cal.data[station]
        if not len(gains) or not np.isfinite(gains['rscale']).all():
            raise ValueError(f'missing/nonfinite shared gains for {station}')
        np.testing.assert_array_equal(gains['rscale'], gains['lscale'])
    clean = native.copy()
    lookup = {str(s): n for n, s in enumerate(clean.tarr['site'])}
    clean.tarr = clean.tarr[[lookup[str(s)] for s in cal.tarr['site']]].copy()
    clean.tkey = {s: n for n, s in enumerate(clean.tarr['site'])}
    clean = cal.applycal(clean, interp='nearest', extrapolate=True)
    if clean.polrep != 'circ':
        raise ValueError('gain application lost native circular representation')
    clean.save_uvfits(str(out / I_PRODUCTS[1]))
    np.savez(out / 'shared-gains.npz', **{str(s): cal.data[s] for s in stations})
    dump(out / 'polarization-input.json', dict(native_load=native_meta, stations=stations,
         shared_R_equals_L=True, uv_min=float(c['uv_min']),
         zero_baselines_excluded_before_shared_gain_solve=True,
         official_I_zero_baseline_workflow_unchanged=True,
         source='pre-I native circular, scan averaged; preserved unequal RL/LR errors',
         gain_model=f"final 5uas-restored {c['fov_uas']:g}uas fixed I; all stations including PV",
         no_I_noise_inflation=True, no_reverse_taper=True, rf_hz=float(native.rf),
         fov_uas=float(fixed.fovx() / eh.RADPERUAS), flux_jy=float(fixed.total_flux()),
         note='0.6Jy is the official compact-flux constraint, not post-hoc image renormalization'))


def load_gaincal_base(path, c, eh):
    """Only consume our verified I-stage product; never scan-average it twice.

    UVFITS does not preserve all calibration flags/feed coefficients. Restore the
    explicit I-stage contract, not the loader defaults or raw release filename.
    """
    obs = eh.obsdata.load_uvfits(str(path), polrep='circ', remove_nan=False)
    feed = {row['site']: row for row in json.loads(Path(c['feed_rotation_file']).read_text())['table']}
    for row in obs.tarr:
        for key in ('fr_par', 'fr_elev', 'fr_off'):
            row[key] = feed[str(row['site'])][key]
    obs.tarr['dr'] = obs.tarr['dl'] = 0j
    obs.ampcal = obs.phasecal = obs.frcal = True
    obs.dcal = False
    obs.official_metadata = dict(fsys_applied=False, source='verified I-stage clean-gaincal-pol.uvfits',
                                 scan_average=False, reason='already scan averaged in I stage')
    return obs


def validate_fixed(image, c, eh):
    import numpy as np

    if ((image.xdim, image.ydim) != (c['npix'], c['npix'])
            or not np.isclose(image.fovx() / eh.RADPERUAS, c['fov_uas'], rtol=0, atol=1e-8)
            or not np.isfinite(image.imvec).all() or np.any(image.imvec <= 0)):
        raise ValueError('fixed I must already be positive finite on the shared sky grid; '
                         'no regrid or renormalization')


def initial_image(fixed, c):
    import numpy as np

    image = fixed.copy()
    p = c['initial_m'] * fixed.imvec * np.exp(2j * np.deg2rad(c['initial_evpa_deg']))
    image.add_qu(p.real.reshape(fixed.ydim, fixed.xdim), p.imag.reshape(fixed.ydim, fixed.xdim))
    image.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
    return image


def restart(image, fixed, c, eh):
    """Blur fractional complex P, I-weighted: |P/I|<1 stays inside convex disk.

    K*(I p)/(K*I) is an I-weighted local complex fraction; multiply by the
    ORIGINAL I. This is NOT literal Q/U-only convolution with I restored:
    it reweights P to retain the physical domain at initialization only.
    It changes neither the objective nor the fitted output; it is not a clip.
    """
    import numpy as np

    blurred = image.blur_circ(c['restart_blur_uas'] * eh.RADPERUAS,
                              c['restart_blur_uas'] * eh.RADPERUAS)
    p = (blurred.qvec + 1j * blurred.uvec) / blurred.imvec
    if not np.isfinite(p).all() or np.any(np.abs(p) >= 1):
        raise ValueError('restart outside physical fractional polarization domain')
    init = fixed.copy()
    p *= (1 - c['restart_domain_epsilon']) * fixed.imvec
    init.add_qu(p.real.reshape(fixed.ydim, fixed.xdim), p.imag.reshape(fixed.ydim, fixed.xdim))
    init.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
    return init


def image_p(obs, fixed, init, c, eh):
    import numpy as np
    import scipy.optimize as opt

    class ShiftedImager(eh.imager.Imager):
        def objfunc(self, x):
            return super().objfunc(x) + sum(c['p_data'].values()) + c['p_reg']['hw'] * np.log(2.)

    records, original = [], opt.minimize

    def capture(fun, x0, *args, **kwargs):
        first = float(fun(x0))
        jac = kwargs.get('jac')
        if not np.isfinite(first) or (callable(jac) and not np.isfinite(jac(x0)).all()):
            raise ValueError('nonfinite initial P objective/gradient; optimization not started')
        result = original(fun, x0, *args, **kwargs)
        if not np.isfinite(result.fun) or not np.isfinite(result.jac).all():
            raise ValueError('nonfinite P optimizer result')
        records.append(dict(initial_objective=first, objective=float(result.fun), nit=int(result.nit),
                            success=bool(result.success), status=int(result.status),
                            message=str(result.message),
                            transformed_gradient_inf=float(np.max(np.abs(result.jac)))))
        return result

    imager = ShiftedImager if c['shift_objective_constant'] else eh.imager.Imager
    imgr = imager(obs, init, prior_im=init, flux=float(fixed.total_flux()), pol='P', transform=['mcv'],
                  data_term=c['p_data'], reg_term=c['p_reg'], norm_reg=c['norm_reg'], norm_init=False,
                  maxit=c['maxit'], stop=c['stop'], ttype='direct')
    with patch.object(opt, 'minimize', capture):
        image = imgr.make_image_P(show_updates=False)
    np.testing.assert_array_equal(image.imvec, fixed.imvec)
    if np.any(image.vvec != 0):
        raise ValueError('P-only solve changed V=0')
    return image, records


def d_dictionary(dterms):
    """JSON boundary: station -> [DR,DL], each represented [real,imag]."""
    import numpy as np

    result = {}
    for station, pair in dterms.items():
        a = np.asarray(pair)
        if a.shape == (2, 2) and not np.iscomplexobj(a):
            a = a[:, 0] + 1j * a[:, 1]
        if a.shape != (2,) or not np.isfinite(a).all():
            raise ValueError(f'invalid D terms for {station}')
        result[str(station)] = [[float(z.real), float(z.imag)] for z in a]
    return result


def fixed_d_image(base, fixed, dterms, c, core, eh):
    """Fresh base -> fsys once -> sampled D once -> P only. NO D fit.

    Single maxit-limited optimization is provisional pending fixed-D pilot;
    this is not yet a validated multiscale MC imaging schedule.
    """
    noisy_base = core.add_fsys(base.copy(), c['fsys'])
    corrected = core.make_imaging_obs(core.apply_fixed_d(noisy_base, dterms))
    return image_p(corrected, fixed, initial_image(fixed, c), c, eh)


def save_image(image, prefix, c, eh):
    image.save_fits(str(prefix) + '-raw-IQU.fits')
    b = c['display_blur_uas'] * eh.RADPERUAS
    display = image.blur_circ(b, b)
    display.save_fits(str(prefix) + '-display-IQU-additional19uas.fits')
    return display


def run_fiducial(base, fixed, c, core, eh, out):
    image, dterms, steps = initial_image(fixed, c), {}, []
    for r in range(1, c['rounds'] + 1):
        init = restart(image, fixed, c, eh) if r in c['restart_rounds'] else image.copy()
        noisy_base = core.add_fsys(base.copy(), c['fsys'])
        corrected = core.apply_fixed_d(noisy_base.copy(), dterms) if dterms else noisy_base.copy()
        image, optimizers = image_p(core.make_imaging_obs(corrected), fixed, init, c, eh)
        # Fit TOTAL D against the same uncorrected noisy base, fixed errors.
        solved = core.fit_forward_d(noisy_base.copy(), image, sites=c['solve_sites'])
        dterms = d_dictionary({str(row['site']): [row['dr'], row['dl']]
                              for row in solved.tarr if str(row['site']) in c['solve_sites']})
        dump(out / f'round-{r:02d}-dterms.json', dterms)
        save_image(image, out / f'round-{r:02d}', c, eh)
        steps.append(dict(round=r, blur_restart=r in c['restart_rounds'], optimizers=optimizers))
    # Align the delivered P image with the last solved D, with no further D fit.
    image, optimizers = image_p(core.make_imaging_obs(
        core.apply_fixed_d(core.add_fsys(base.copy(), c['fsys']), dterms)), fixed, image, c, eh)
    save_image(image, out / 'final', c, eh)
    dump(out / 'final-dterms.json', dterms)
    dump(out / 'steps.json', steps + [dict(stage='final-P-only', optimizers=optimizers)])
    dump(out / 'termination.json', dict(rounds=c['rounds'], reason='fixed round budget',
         scientifically_converged=False,
         maxit=c['maxit'], maxit_scope='numerical per-call cap, not a uniquely mandated Table6 value',
         note='stop is numerical ftol/gtol, not a paper-mandated science threshold'))


def sample_dterms(days, c, day, core):
    """Use the shared four-day estimator/sampler; day-specific deterministic seed."""
    import numpy as np

    daily = {str(d): {s: pair for s, pair in d_dictionary(days[str(d)]).items() if s in c['solve_sites']}
             for d in DAYS}
    if any(set(v) != set(c['solve_sites']) for v in daily.values()):
        raise ValueError('each D source day must contain AZ/LM/PV, both hands')
    distribution = core.d_distribution(daily, ddof=c['covariance_ddof'])
    sd = c['residual_d_sigma'] / (np.sqrt(2) if c['residual_d_convention'] == 'complex_rms' else 1)
    # Distinct days must not accidentally share the same random realization.
    seed = int(np.random.SeedSequence([c['seed'], day]).generate_state(1)[0])
    return core.sample_dterms(distribution, c['solve_sites'] + c['fixed_sites'], c['uncertainty_draws'],
                              seed, residual_component_sigma=sd)


def image_statistics(image):
    """Whole-image quantities only; aperture/profile/beta2 definitions await review."""
    import numpy as np

    i, p = image.imvec, image.qvec + 1j * image.uvec
    return dict(net_m=float(abs(p.sum()) / i.sum()), resolved_m=float(np.abs(p).sum() / i.sum()),
                integrated_evpa_rad=float(.5 * np.angle(p.sum())) if abs(p.sum()) else None)


def run_uncertainty(plan, args, base, fixed, c, core, eh, out):
    import numpy as np

    if args.published_dterms:
        source = json.loads(args.published_dterms.read_text())
        days = source['days']
    else:
        days = {str(common.complete_manifest(root)['config']['day']):
                json.loads((root / 'final-dterms.json').read_text()) for root in args.d_roots}
        source = dict(source_kind='four local consistent fiducial stages',
                      roots=[str(p) for p in args.d_roots])
    samples = sample_dterms(days, c, plan['day'], core)
    core.write_d_samples(samples, out / 'all-draws.json', out / 'all-draws.csv')
    sites, draws = samples['sites'], samples['draws']
    matrix = np.asarray([[[complex(*draw[s][h]) for h in ('R', 'L')] for s in sites] for draw in draws])
    np.savez(out / 'all-draws.npz', dterms=matrix, stations=np.asarray(sites), draw_id=np.arange(len(draws)))
    dump(out / 'D-source.json', dict(source=source, ddof=1, residual_sigma=c['residual_d_sigma'],
         residual_convention=c['residual_d_convention'], derived_day_seed=samples['seed']))
    stats = []
    for draw_id, draw in enumerate(draws[:plan['draws']]):
        dterms = {s: [v['R'], v['L']] for s, v in draw.items()}
        image, optimizers = fixed_d_image(base, fixed, dterms, c, core, eh)
        displayed = save_image(image, out / f'draw-{draw_id:04d}', c, eh)
        stats.append(dict(draw_id=draw_id, optimizers=optimizers, **image_statistics(displayed)))
        dump(out / 'statistics.json', stats)
    dump(out / 'sampling-receipt.json', dict(configured_draws=len(draws), executed_draws=len(stats),
         draw_ids=list(range(len(stats))), label=plan['label'], failed_draws=0,
         optimizer_unsuccessful_draw_ids=[s['draw_id'] for s in stats
                                         if any(not r['success'] for r in s['optimizers'])],
         exception_policy='abort transaction; never publish a half-complete ensemble stage',
         imaging_schedule='single P optimization; multiscale schedule/pilot validation pending',
         mc_schedule_validated=False,
         seed=samples['seed'], master_seed=c['seed'], D_refitted=False, no_thermal_bootstrap=True,
         deferred_statistics=['aperture EVPA/P', 'beta2', 'azimuth profiles'],
         scientifically_converged=False))


def main(argv=None):
    args = parser().parse_args(argv)
    with common.verification_session():
        plan = build_plan(args)
        if args.plan:
            print(json.dumps(plan, indent=2))
            return plan
        import ehtim as eh
        from astropy.utils import iers
        from eht_axion.research import eht_official as core

        iers.conf.auto_download = False
        c = plan['parameters']
        # Refuse a shadow import: all actually loaded library source must belong
        # to the package whose source bytes were declared before execution.
        declared = {r['path'] for r in plan['inputs']}
        if any(str(p.resolve()) not in declared for p in Path(eh.__file__).resolve().parent.rglob('*.py')):
            raise ValueError('loaded ehtim source is outside declared package closure')
        with common.stage_directory(args.output_root, 'research:eht-official:' + args.stage, plan,
                                    [Path(r['path']) for r in plan['inputs']]) as out:
            if out is None:
                raise FileExistsError('new stage/day root required')
            dump(out / 'plan.json', plan)
            if args.stage == 'I':
                run_i(plan, core, eh, out)
            else:
                fixed = eh.image.load_fits(str(args.i_root / I_PRODUCTS[0]))
                validate_fixed(fixed, c, eh)
                base = load_gaincal_base(args.i_root / I_PRODUCTS[1], c, eh)
                base = base.flag_uvdist(uv_min=c['uv_min'])
                if args.stage == 'fiducial':
                    run_fiducial(base, fixed, c, core, eh, out)
                else:
                    run_uncertainty(plan, args, base, fixed, c, core, eh, out)
        return plan


if __name__ == '__main__':
    main()
