#!/usr/bin/env python3
"""EHT imaging, D-term Monte Carlo: fixed-I/fixed-D P images for every D draw (4 days x draws).

Science values (draw count, seed, fsys, uv cut, D-term statistics, sky grid,
blur and the P solver settings eht_imaging.monte_carlo) come from
experiment.yaml. Inputs per day: the fixed I image and gain-calibrated data of the I stage
(eht/intensity/day-<day>) and the refined D terms of the refine stage (eht/refine/day-<day>).
Execution options (workers, threads, memory) are CLI.
Input-only: no export, SSH, D fit or thermal draw.

Receipts are this script's batch records (batch.json, taskreceipt.json, finish.json), not stage manifests.
A floating-point gap certificate is numerical, not scientific acceptance.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import redirect_stdout, redirect_stderr
import fcntl
import importlib.metadata
import json
import multiprocessing
import os
import resource
import time
import traceback

from eht_axion.environment import DATA_ROOT_KEY, add_data_root_option, data_root, fill_defaults, layout_path
from eht_axion import experiment, provenance
from eht_axion.experiment import DEFAULT_PATH, EHT_DAYS
from eht_axion.provenance import contained_path as contained

DAYS = tuple(str(day) for day in EHT_DAYS)  # bundle and task keys are strings
_CACHE = {}


def science_settings(config_path):
    """The experiment values this Monte Carlo uses (recorded in batch.json)."""
    from eht_axion.experiment import read_experiment

    config = read_experiment(config_path)
    c, image = config['eht_imaging'], config['shared']['image']
    return dict(npix=image['shape'], fov_uas=image['fov_uas'], blur_fwhm_uas=image['blur_fwhm_uas'],
                uv_min=c['uv_min'], fsys=c['fsys'], covariance_ddof=c['covariance_ddof'],
                residual_component_sigma=c['residual_d_sigma'], solve_sites=c['solve_sites'],
                fixed_sites=c['fixed_sites'], sites=sorted(c['solve_sites'] + c['fixed_sites']),
                p_data=c['p_data'], p_reg=c['p_reg'], norm_reg=c['norm_reg'],
                feed_rotation=c['feed_rotation'],
                draws_per_day=c['uncertainty_draws'], seed=c['seed'], monte_carlo=c['monte_carlo'])


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def atomic(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('wb') as f:
        f.write(encoded(value) + b'\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


I_FILES = {'I': 'final-fixedI.fits', 'obs': 'clean-gaincal-pol.uvfits'}
D_FILES = ('final-dterms.json', 'termination.json')


def _by_day(roots, stage, config):
    """Check each completed stage and index it by its recorded day."""
    result = {}
    for root in roots:
        root = Path(root).resolve(strict=True)
        manifest = provenance.verify_stage(root)
        if manifest['stage'] != stage:
            raise ValueError(f'{root}: expected stage {stage}, found {manifest["stage"]}')
        parameters = manifest['config']['parameters']
        locations = {key: parameters[key] for key in experiment.EHT_LOCATIONS}
        if parameters != json.loads(json.dumps(experiment.eht_imaging_parameters(config, locations))):
            raise ValueError(f'{root}: imaging parameters differ from experiment.yaml eht_imaging')
        day = str(manifest['config']['day'])
        if day in result:
            raise ValueError(f'two {stage} roots for day {day}')
        result[day] = (root, manifest)
    if set(result) != set(DAYS):
        raise ValueError(f'{stage}: need exactly the days {DAYS}, got {sorted(result)}')
    return result


def verify_sources(i_roots, d_roots, config):
    """{day: {I, obs, dterms}} from the I and refine stages; every refine P solve certified."""
    images = _by_day(i_roots, 'eht:i', config)
    dterms = _by_day(d_roots, 'eht:refine', config)
    sources = {}
    for day in DAYS:
        i_root, d_root = images[day][0], dterms[day][0]
        provenance.require_declared_outputs(i_root, [i_root / name for name in I_FILES.values()])
        provenance.require_declared_outputs(d_root, [d_root / name for name in D_FILES])
        closure = set(provenance.source_closure([provenance.manifest_path(d_root)]))
        if any((i_root / name).resolve() not in closure for name in I_FILES.values()):
            raise ValueError(f'day {day}: the refine stage does not descend from this I stage')
        termination = json.loads((d_root / 'termination.json').read_text())
        if termination.get('all_P_numerically_certified') is not True:
            raise ValueError(f'day {day}: not every P solve of the refine stage is certified')
        sources[day] = {key: str(i_root / name) for key, name in I_FILES.items()}
        sources[day]['dterms'] = str(d_root / 'final-dterms.json')
    return sources


def configure_threads(n):
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
        os.environ[key] = str(n)
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    os.environ['MKL_DYNAMIC'] = 'FALSE'


def code_identity():
    """Git commit of the repository (tracked files must be clean) and package versions.

    The commit pins the scripts, eht_axion and pixi.lock; installed packages are
    recorded by version only.
    """
    from eht_axion import provenance
    packages = {}
    for name in ('ehtim', 'numpy', 'scipy', 'astropy', 'threadpoolctl'):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = 'not-installed'
    return dict(provenance.code_identity(), packages=packages, python=sys.version.split()[0])


def make_draws(sources, n, seed, science):
    import numpy as np
    from eht_axion.eht import calibration as core
    solve = set(science['solve_sites'])
    daily = {d: json.loads(Path(sources[d]['dterms']).read_text()) for d in DAYS}
    daily = {d: {s: v for s, v in terms.items() if s in solve} for d, terms in daily.items()}
    if any(set(v) != solve for v in daily.values()):
        raise ValueError(f'each fiducial must contain {sorted(solve)}')
    distribution = core.d_distribution(daily, ddof=science['covariance_ddof'])
    return {d: core.sample_dterms(distribution, science['sites'], n,
                int(np.random.SeedSequence([seed, int(d)]).generate_state(1)[0]),
                residual_component_sigma=science['residual_component_sigma'],
                fixed_sites=science['fixed_sites']) for d in DAYS}


def enumerate_tasks(draws):
    return [dict(id=f'{d}-{k:04d}', day=d, draw=k, dterms=draws[d]['draws'][k])
            for k in range(len(draws[DAYS[0]]['draws'])) for d in DAYS]


def attempt_budget(config, numerical_attempt):
    if numerical_attempt == 1:
        return config['maxiter'], config['time_limit']
    return config['retry_maxiter'], config['retry_time_limit']


def verified_receipt(directory, task):
    p = directory / 'taskreceipt.json'
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    if r['task'] != task:
        raise ValueError('receipt task mismatch: ' + str(p))
    for rel in r['outputs']:
        contained(directory, rel)
    if r['status'] == 'certified':
        cert = json.loads((directory / 'certificate.json').read_text())
        report = cert['solver']
        if not (report['numerically_certified'] and report['best_gap'] <= cert['gap_atol']
                and report['best_gap'] >= -1e-10):
            raise ValueError('invalid certified receipt')
    return r


def scan_task(out, task):
    """Uncertified/error/interrupted attempts are never counted as completed."""
    receipts = []
    for directory in sorted((out / 'attempts' / task['id']).glob('attempt-*')):
        r = verified_receipt(directory, task)
        if r:
            receipts.append((directory, r))
    completed = [str(p.relative_to(out)) for p, r in receipts if r['status'] == 'certified']
    numerical = sum(r['status'] == 'uncertified' for _, r in receipts)
    errors = sum(r['status'] == 'error' for _, r in receipts)
    return completed, numerical, errors


def init_worker(threads):
    configure_threads(threads)
    import numpy  # noqa: F401 -- only after environment setup
    try:
        from threadpoolctl import threadpool_limits, threadpool_info
    except ImportError:
        return  # Optional inspection; pre-import environment is authoritative.
    global _LIMITS
    _LIMITS = threadpool_limits(limits=threads)
    if any(x['num_threads'] > threads for x in threadpool_info()):
        raise RuntimeError('BLAS thread limit not honored')


def run_attempt(sources, task, directory, config, numerical_attempt, science):
    directory = Path(directory)
    started = time.monotonic()
    status = 'error'
    with (directory / 'job.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
        try:
            import numpy as np
            from scripts.eht import reconstruct as driver
            from eht_axion.eht import calibration as core, convex as opt
            eh = core._library()
            d = task['day']
            if d not in _CACHE:
                fixed = eh.image.load_fits(sources[d]['I'])
                driver.validate_fixed(fixed, science, eh)
                base = driver.load_gaincal_base(Path(sources[d]['obs']), science, eh)
                base = base.flag_uvdist(uv_min=science['uv_min'])
                core._native(base)
                _CACHE[d] = fixed, base
            fixed, base = _CACHE[d]
            present = set(map(str, base.tarr['site']))
            if not present <= task['dterms'].keys():
                raise ValueError('draw does not cover observation stations')
            terms = {s: task['dterms'][s] for s in present}
            noisy = core.add_fsys(base.copy(), science['fsys'])
            obs = core.make_imaging_obs(core.apply_fixed_d(noisy, terms))
            problem = opt.exact_cartesian_problem(fixed, obs, p_data=science['p_data'],
                                                  p_reg=science['p_reg'], norm_reg=science['norm_reg'])
            maxiter, limit = attempt_budget(config, numerical_attempt)
            initial = config['initial_p_fraction'] * fixed.imvec.astype(complex)
            P, report, trace = opt.solve_problem(problem, initial,
                gap_atol=config['gap_atol'], maxiter=maxiter, time_limit=limit,
                check_every=config['check_every'],
                callback=lambda row: atomic(directory / 'progress.json', row))
            image = fixed.copy()
            image.add_qu(P.real.reshape(fixed.ydim, fixed.xdim), P.imag.reshape(fixed.ydim, fixed.xdim))
            image.add_pol_image(np.zeros((fixed.ydim, fixed.xdim)), 'V')
            driver.save_image(image, directory / 'best', science, eh)
            cert = dict(solver=report, gap_atol=config['gap_atol'], trace=trace,
                        numerical_attempt=numerical_attempt, budget=dict(maxiter=maxiter, time_limit=limit),
                        unused_station_draws=sorted(set(task['dterms']) - present),
                        science_acceptance=False,
                        initial_P=f"{config['initial_p_fraction']:g} I real, no warm start",
                        worker_lifetime_peak_rss_gib=(
                            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20))
            atomic(directory / 'certificate.json', cert)
            status = 'certified' if report['numerically_certified'] else 'uncertified'
        except Exception:
            traceback.print_exc()
    outputs = [p.name for p in sorted(directory.iterdir()) if p.is_file()]
    receipt = dict(schema='eht-sample-task', task=task,
                   status=status, numerical_attempt=numerical_attempt,
                   seconds=time.monotonic()-started, outputs=outputs,
                   scientifically_certified=False)
    atomic(directory / 'taskreceipt.json', receipt)
    return status


def memory_projection(usage, worker_rss, workers):
    """Reserve the entire pool before lazy allocations materialize."""
    reserve = max(0, workers * 2**30 - worker_rss)
    return dict(worker_rss_bytes=worker_rss, unmaterialized_reserve_bytes=reserve,
                projected_memory_bytes=usage + reserve)


def child_rss():
    """Direct pool children plus tiny resource tracker; no private pool APIs."""
    total = 0
    children = Path(f'/proc/self/task/{os.getpid()}/children').read_text().split()
    for pid in children:
        try:
            rows = Path(f'/proc/{pid}/status').read_text().splitlines()
            total += int(next(r for r in rows if r.startswith('VmRSS:')).split()[1]) * 1024
        except (OSError, ValueError, StopIteration):
            pass  # Exited child: zero is conservative for reservations.
    return total


def resource_wait_expired(blocked_since, now, active, pending, admitted, timeout):
    if active or not pending or admitted:
        return None, False
    since = now if blocked_since is None else blocked_since
    return since, now - since >= timeout


def memory_usage():
    """Account cgroup if available; otherwise conservative sum of this UID's RSS."""
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            if line.startswith('0::'):
                relative = Path(line[3:])
                account = next((part for part in relative.parts
                                if part.startswith('user-') and part.endswith('.slice')), None)
                if account is not None:
                    parts = relative.parts[:relative.parts.index(account)+1]
                    p = Path('/sys/fs/cgroup').joinpath(*parts[1:]) / 'memory.current'
                    return int(p.read_text()), 'account cgroup-v2 memory.current'
                p = Path('/sys/fs/cgroup') / line[3:].lstrip('/') / 'memory.current'
                return int(p.read_text()), 'cgroup-v2 memory.current'
    except (OSError, ValueError):
        pass
    rss = 0
    for p in Path('/proc').glob('[0-9]*/status'):
        try:
            rows = p.read_text().splitlines()
            uid = next(r for r in rows if r.startswith('Uid:')).split()[1]
            if int(uid) == os.getuid():
                rss += int(next(r for r in rows if r.startswith('VmRSS:')).split()[1]) * 1024
        except (OSError, ValueError, StopIteration):
            pass
    return rss, 'UID summed RSS (shared pages counted repeatedly)'


def controller(args, out):
    science = science_settings(args.config)
    sources = verify_sources(args.i_root, args.d_root, experiment.read_experiment(args.config))
    configure_threads(args.threads)
    config = dict(draws_per_day=science['draws_per_day'], seed=science['seed'], **science['monte_carlo'],
                  threads=args.threads)
    plan = dict(schema='eht-sample-batch', inputs=sources,
                code=code_identity(), config=config, science=science)
    planpath = out / 'batch.json'
    if planpath.exists():
        if json.loads(planpath.read_text()) != plan:
            raise ValueError('resume requires the same inputs, science settings and code commit')
    else:
        if any(p.name != '.controller.lock' for p in out.iterdir()):
            raise ValueError('new output root must be empty')
        atomic(planpath, plan)
    draws = make_draws(sources, science['draws_per_day'], science['seed'], science)
    p = out / 'all-draws.json'
    if p.exists():
        if json.loads(p.read_text()) != draws:
            raise ValueError('persisted draw mismatch')
    else:
        atomic(p, draws)
    tasks = enumerate_tasks(draws)
    pending, done, exhausted = deque(), {}, []
    for t in tasks:
        completed, numerical, errors = scan_task(out, t)
        if completed:
            done[t['id']] = completed[0]
        elif numerical >= 2 or errors >= 2:
            exhausted.append(t['id'])
        else:
            pending.append((t, numerical + 1))
    active = {}
    submitted = 0
    blocked_since = None
    peak_usage = peak_projected = 0
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker,
            initargs=(args.threads,), mp_context=multiprocessing.get_context('spawn')) as pool:
        while pending or active:
            if args.task_limit is not None and submitted >= args.task_limit and not active:
                break
            usage, scope = memory_usage()
            available = int(next(r.split()[1] for r in Path('/proc/meminfo').read_text().splitlines()
                                 if r.startswith('MemAvailable:'))) * 1024
            projection = memory_projection(usage, child_rss(), args.workers)
            peak_usage = max(peak_usage, usage)
            peak_projected = max(peak_projected, projection['projected_memory_bytes'])
            # Admission only: full pool reservation; never kill active workers.
            admitted = (bool(pending) and len(active) < args.workers
                        and (args.task_limit is None or submitted < args.task_limit)
                        and projection['projected_memory_bytes'] <= args.max_active_memory_gib*2**30
                        and available >= projection['unmaterialized_reserve_bytes'] + 2*2**30)
            blocked_since, expired = resource_wait_expired(
                blocked_since, time.monotonic(), len(active), len(pending), admitted,
                args.resource_wait)
            memory_state = dict(memory_bytes=usage, memory_scope=scope,
                                available_memory_bytes=available, **projection,
                                peak_memory_bytes=peak_usage, peak_projected_memory_bytes=peak_projected)
            if expired:
                state = dict(status='resource_blocked', resumable=True,
                             certified=len(done), total=len(tasks), active=0, pending=len(pending),
                             exhausted=exhausted, blocked_seconds=time.monotonic()-blocked_since,
                             resource_wait_seconds=args.resource_wait, **memory_state)
                atomic(out / 'state.json', state)
                print(json.dumps(state, sort_keys=True), file=sys.stderr, flush=True)
                return 3  # Pool context exits, freeing idle caches; no finish receipt.
            if admitted:
                t, number = pending.popleft()
                parent = out / 'attempts' / t['id']
                parent.mkdir(parents=True, exist_ok=True)
                ordinal = 1 + max([int(p.name.split('-')[1]) for p in parent.glob('attempt-*')] or [0])
                directory = parent / f'attempt-{ordinal:04d}'
                directory.mkdir()  # Interrupted attempts remain untouched forever.
                future = pool.submit(run_attempt, sources, t, directory, config, number, science)
                active[future] = t
                submitted += 1
                continue
            atomic(out / 'state.json', dict(certified=len(done), total=len(tasks),
                   active=len(active), pending=len(pending), exhausted=exhausted,
                   **memory_state, updated_unix=time.time()))
            ready, _ = wait(active, timeout=5, return_when=FIRST_COMPLETED) if active else ([], [])
            if not active:
                time.sleep(5)
            for f in ready:
                t = active.pop(f)
                try:
                    f.result()
                except Exception as exc:
                    # A dead pool is not a numerical failure; preserve attempts and exit for restart.
                    raise RuntimeError('worker/process failure; restart same command to recover') from exc
                completed, numerical, errors = scan_task(out, t)
                if completed:
                    done[t['id']] = completed[0]
                elif numerical < 2 and errors < 2:
                    pending.append((t, numerical + 1))
                else:
                    exhausted.append(t['id'])
    if len(done) == len(tasks):
        finish = dict(schema='eht-sample-finish',
                      numerically_certified=len(done), scientifically_certified=False,
                      receipts=dict(sorted(done.items())))
        if (out / 'finish.json').exists():
            if json.loads((out / 'finish.json').read_text()) != finish:
                raise ValueError('finish receipt mismatch')
        else:
            atomic(out / 'finish.json', finish)
        # state.json is the progress of a running batch; a finished batch has finish.json.
        (out / 'state.json').unlink(missing_ok=True)
        return 0
    atomic(out / 'state.json', dict(certified=len(done), total=len(tasks),
                                   active=0, pending=len(pending), exhausted=exhausted))
    return 0 if args.task_limit is not None and not exhausted else 2


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--i-root', type=Path, action='append',
                   help='completed I stage; give all four days (default: eht/intensity/day-<day>)')
    p.add_argument('--d-root', type=Path, action='append',
                   help='completed refine stage; give all four days (default: eht/refine/day-<day>)')
    p.add_argument('--output-root', type=Path, help='batch root; resumable (default: eht/sample)')
    p.add_argument('--workers', type=int, default=64)
    p.add_argument('--resource-wait', type=float, default=300.,
                   help='idle admission timeout seconds; exit 3, safe to resume')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--config', type=Path, help='experiment parameter file (default: configs/experiment.yaml)')
    p.add_argument('--task-limit', type=int,
                   help='attempt cap; full draw plan retained, no partial finish receipt')
    p.add_argument('--resume', action='store_true',
                   help='explicit resume intent; resuming the same batch is safe')
    p.add_argument('--max-active-memory-gib', type=float, default=100.)
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(p, args, {'output_root': ('eht.sample', {})})
    for name, key in (('i_root', 'eht.intensity'), ('d_root', 'eht.refine')):
        if getattr(args, name) is None:
            root = data_root(args)
            if root is None:
                p.error(f'Provide --{name.replace("_", "-")} (four times), or --data-root / ${DATA_ROOT_KEY}')
            setattr(args, name, [layout_path(root, key, day=day) for day in EHT_DAYS])
    import math
    if args.config is None:
        args.config = DEFAULT_PATH
    if (min(args.workers, args.threads) < 1
            or (args.task_limit is not None and args.task_limit < 1)
            or not all(math.isfinite(v) and v > 0 for v in (args.max_active_memory_gib, args.resource_wait))):
        p.error('positive workers, threads, memory and wait required')
    out = args.output_root.resolve()
    for root in args.i_root + args.d_root:
        root = root.resolve()
        if out.is_relative_to(root) or root.is_relative_to(out):
            p.error('input and output roots must be separate, non-nested directories')
    out.mkdir(parents=True, exist_ok=True)
    fd = os.open(out / '.controller.lock', os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return controller(args, out)
    finally:
        os.close(fd)  # Kernel releases lock even after crash; no stale-lock timeout.


if __name__ == '__main__':
    raise SystemExit(main())
