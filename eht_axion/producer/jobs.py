"""Bounded concurrent, per-job atomic RAPTOR runner; never a run transaction.

The build is produced by the core-owned setup/make Python entry point. This
module only consumes its immutable, declared products. Failed attempts remain
outside the stage temporary directory. Resume requires both verified lineage
and the success/schema receipts, never mere directory or HDF5 existence.
"""
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextvars import ContextVar, copy_context
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time

from ..pipeline import common
from .outputs import read_core_constants, validate_job_output

PARAMETER_UNITS = {
    "MBH": "(Msun)", "DISTANCE": "(kpc)", "M_UNIT": "(g)", "R_LOW": "(-)", "R_HIGH": "(-)",
    "INCLINATION": "(deg)", "IMG_WIDTH": "(pixels)", "IMG_HEIGHT": "(pixels)",
    "CAM_SIZE_X": "(Rg)", "CAM_SIZE_Y": "(Rg)", "FREQS_PER_DEC": "(-)", "FREQ_MIN": "(Hz)",
    "STEPSIZE": "(-)", "AXION_ALPHA": "(-)", "AXION_OMEGA": "(-)", "AXION_NORM": "(-)",
    "AXION_PHASE": "(-)", "AXION_TABLE": "(-)", "MAX_LEVEL": "(-)",
}
_INTEGER_FIELDS = {"IMG_WIDTH", "IMG_HEIGHT", "AXION_PHASE", "MAX_LEVEL"}


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _token(value, name):
    value = str(value)
    if not value or len(os.fsencode(value)) > 255 or any(c.isspace() or c == "\0" for c in value):
        raise ValueError(f"{name} must be one nonempty token of at most 255 bytes")
    return value


def parameter_values(parameters: dict) -> dict:
    """Accept plan's explicit {unit,value} records or already-normalized values."""
    result = {}
    for name, item in parameters.items():
        if isinstance(item, dict):
            if set(item) != {"unit", "value"} or item["unit"] != PARAMETER_UNITS.get(name):
                raise ValueError(f"Invalid explicit unit/value for {name}")
            result[name] = item["value"]
        else:
            result[name] = item
    return result


def serialize_parameters(parameters: dict) -> str:
    """Emit the exact 19 LABEL UNIT VALUE triples accepted by the new core."""
    if set(parameters) != set(PARAMETER_UNITS):
        raise ValueError("model.in requires exactly the 19 declared core parameters")
    p = parameter_values(parameters)
    for name, value in p.items():
        if name == "AXION_TABLE":
            _token(value, name)
        elif name in _INTEGER_FIELDS:
            _integer(value, name)
        elif isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
            raise ValueError(f"{name} must be a finite number")
    positive = {"MBH", "DISTANCE", "M_UNIT", "R_LOW", "R_HIGH", "CAM_SIZE_X", "CAM_SIZE_Y",
                "FREQS_PER_DEC", "FREQ_MIN", "STEPSIZE"}
    if any(p[name] <= 0 for name in positive):
        raise ValueError("Core scale/model parameters must be positive")
    if not 0 <= p["INCLINATION"] <= 180:
        raise ValueError("INCLINATION must be in 0..180 degrees")
    if (p["IMG_WIDTH"] <= 0 or p["IMG_WIDTH"] != p["IMG_HEIGHT"] or p["IMG_WIDTH"] % 10
            or (p["IMG_WIDTH"] // 10) ** 2 > 2**31 - 1):
        raise ValueError("Square image resolution must be a positive multiple of 10")
    if p["MAX_LEVEL"] != 1:
        raise ValueError("First-scope runner supports only MAX_LEVEL=1")
    if not 0 <= p["AXION_PHASE"] <= 15:
        raise ValueError("AXION_PHASE must be in 0..15")
    if p["AXION_NORM"] < 0 or p["AXION_ALPHA"] < 0 or p["AXION_OMEGA"] < 0:
        raise ValueError("Axion parameters must be nonnegative")
    if p["AXION_NORM"] > 0 and (p["AXION_ALPHA"] <= 0 or p["AXION_OMEGA"] <= 0):
        raise ValueError("Enabled alpha and independent omega must be positive")
    if not math.isfinite(p["MBH"] * 1.989e33) or not math.isfinite(p["DISTANCE"] * 3.086e21):
        raise ValueError("Core cgs conversion overflows")
    lines = []
    for name, unit in PARAMETER_UNITS.items():
        value = p[name]
        text = str(value) if name == "AXION_TABLE" or name in _INTEGER_FIELDS else format(value, ".17g")
        lines.append(f"{name}\t{unit}\t{text}\n")
    return "".join(lines)


def _build_products(build_dir, manifest):
    declared = [build_dir / e["path"] for e in manifest["outputs"]]

    def one(name, preferred=()):
        for relative in preferred:
            if build_dir / relative in declared:
                return build_dir / relative
        matches = [p for p in declared if p.name == name]
        if len(matches) != 1:
            raise ValueError(f"Build must declare one effective {name}; found {len(matches)}")
        return matches[0]

    binary = one("RAPTOR", ("RAPTOR", "run/RAPTOR"))
    grid = one("grid_mks.in", ("grid_mks.in", "run/grid_mks.in"))
    header = one("definitions.h", ("run/definitions.h", "definitions.h"))
    common.require_declared_outputs(build_dir, [binary, grid, header])
    if any(p.is_symlink() for p in (binary, grid, header)):
        raise ValueError("Generated build products must not be symlinks")
    if not os.access(binary, os.X_OK):
        raise ValueError("Declared RAPTOR binary is not executable")
    return binary, grid, header


def _prepare_task(task):
    task = dict(task)
    identifier = task.get("id")
    if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}", identifier):
        raise ValueError("Unsafe task id")
    if task.get("kind") not in {"baseline", "response"}:
        raise ValueError("Task kind must be baseline or response")
    source = Path(task["source_path"]).resolve(strict=True)
    if not source.is_file() or len(os.fsencode(source)) > 255:
        raise ValueError("Snapshot must be a regular file with a path <=255 bytes")
    _integer(task["label"], "snapshot label")
    if isinstance(task["time_tg"], bool) or not math.isfinite(task["time_tg"]):
        raise ValueError("Snapshot time_tg must be finite")
    resolution = _integer(task["resolution"], "resolution", 1)
    slot = _integer(task.get("output_slot", task["phase_index"]), "output slot")
    if slot > 16:
        raise ValueError("Output slot must be in 0..16")
    p = parameter_values(task["parameters"])
    if p["IMG_WIDTH"] != resolution or p["IMG_HEIGHT"] != resolution:
        raise ValueError("Task resolution and model dimensions differ")
    if task["norm"] != p["AXION_NORM"]:
        raise ValueError("Task norm and model norm differ")
    table = None
    if task["kind"] == "baseline":
        if (task["norm"] != 0 or slot != 0 or p["AXION_ALPHA"] != 0 or p["AXION_OMEGA"] != 0
                or p["AXION_PHASE"] != 0 or task.get("omega", 0) != 0
                or task["phase_index"] != 0 or task.get("alpha") is not None):
            raise ValueError("Baseline requires norm/alpha/omega/phase/outputslot zero")
        # Disabled jobs do not consume any radial table or cloud identity.
        task.update(alpha=None, omega=0, table_id=None, table_path=None)
        p["AXION_TABLE"] = "disabled"
    else:
        if task["norm"] <= 0 or slot < 1 or p["AXION_PHASE"] != slot - 1:
            raise ValueError("Response requires positive norm and slot=AXION_PHASE+1")
        if task["alpha"] != p["AXION_ALPHA"] or task["omega"] != p["AXION_OMEGA"]:
            raise ValueError("Task alpha/omega and model parameters differ")
        if "output_slot" in task and task["phase_index"] != p["AXION_PHASE"]:
            raise ValueError("Task phase_index differs from physical AXION_PHASE")
        table = Path(task["table_path"]).resolve(strict=True)
        if Path(p["AXION_TABLE"]).resolve(strict=True) != table:
            raise ValueError("Model AXION_TABLE differs from task table_path")
        if not table.is_file():
            raise ValueError("Axion table must be a regular file")
        p["AXION_TABLE"] = _token(table, "AXION_TABLE")
        task["table_path"] = str(table)
    task.update(source_path=str(source), parameters=p, output_slot=slot)
    serialize_parameters(p)
    # Intentionally whitelist consumed identity: baseline must not acquire a
    # hidden dependency on arbitrary planning annotations or radial tables.
    config = {key: task.get(key) for key in (
        "id", "kind", "label", "time_tg", "source_path", "resolution", "alpha", "omega", "norm",
        "phase_index", "output_slot", "parameters",
    )}
    if table is not None:
        config.update({key: task.get(key) for key in ("table_id", "table_path", "table_spin")})
    return task, config, [source] + ([] if table is None else [table])


def _validate_plan_binding(plan, prepared, manifest, build_dir, grid):
    """Reject drift between a production plan, live footer/table, and build.

    A bare tasks-only plan is the lower-level synthetic runner API. Production
    plans carry config/snapshot/table inventories and must satisfy every gate.
    """
    for task, _, paths in prepared:
        for path in paths:
            digest = common.sha256(path)
            expected = task.get("source_sha256" if path == paths[0] else "table_sha256")
            if expected is not None and expected != digest:
                raise ValueError(f"Declared plan input hash differs from live bytes: {path}")
    expected_commit = plan.get("expected_core_commit")
    if expected_commit is not None:
        if manifest["config"].get("source_commit") != expected_commit:
            raise ValueError("Build source commit differs from expected_core_commit")
        receipt_path = build_dir / "receipt.json"
        common.require_declared_outputs(build_dir, [receipt_path])
        if json.loads(receipt_path.read_text()).get("source_commit") != expected_commit:
            raise ValueError("Build receipt source commit differs from expected_core_commit")
    if "config" not in plan:
        return
    if expected_commit is None:
        raise ValueError("Production plan requires expected_core_commit")
    from ..experiment import validate_ray_tracing
    from ..model import omega
    from .plan import _parameters, canonical_alpha, read_snapshot_footer, validate_radial_table

    config = validate_ray_tracing(plan["config"])
    p = config["raptor"]
    if manifest["config"].get("settings") != p["build"]:
        raise ValueError("Build settings differ from production plan")
    # grid_mks.in's exact legacy field spellings are owned by the core.
    fields = {}
    for line in grid.read_text().splitlines():
        if not line.strip():
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key in fields:
            raise ValueError("Duplicate core grid field")
        fields[key] = float(value.strip())
    expected_grid = {f"nxlone{i+1}": v for i, v in enumerate(p["grid"]["nxlone"])}
    expected_grid.update({f"xprobminx{i+1}": v for i, v in enumerate(p["grid"]["xprobmin"])})
    expected_grid.update({f"xprobmax{i+1}": v for i, v in enumerate(p["grid"]["xprobmax"])})
    expected_grid["hslope"] = p["grid"]["hslope"]
    if fields != expected_grid:
        raise ValueError("Effective build grid differs from production plan")
    macro_path = build_dir / "macros.txt"
    common.require_declared_outputs(build_dir, [macro_path])
    macro_text = macro_path.read_text()
    from .build import _verify_macros
    _verify_macros(macro_text, p["build"])
    macros = dict(re.findall(r"^#define[ \t]+(\w+)[ \t]+([^\n]*)$", macro_text, re.MULTILINE))
    required = {"POL": "1", "num_frequencies": "1", "num_pixels_1d": "10", "NSPIN": "3",
                "SMR": "1", "AMR": "0", "UNIF": "0", "metric": "MKSBHAC",
                "DF": "TH", "IMGFILE": "1", "SPECFILE": "1"}
    if any(macros.get(key, "").strip().strip("()").strip() != value for key, value in required.items()):
        raise ValueError("Unsupported effective camera/POL/frequency/build macros")
    snapshots = {s["label"]: s for s in plan["snapshots"]}
    footer_cache, table_cache = {}, {}
    start, stop, step = (p["snapshots"][k] for k in ("start", "stop", "step"))
    for task, _, _ in prepared:
        label = task["label"]
        if not start <= label <= stop or (label - start) % step:
            raise ValueError("Task label is not on the production snapshot grid")
        snapshot = snapshots.get(label)
        if snapshot is None or Path(snapshot["source_path"]).resolve() != Path(task["source_path"]):
            raise ValueError("Task source differs from plan snapshot inventory")
        source = task["source_path"]
        if source not in footer_cache:
            footer_cache[source] = read_snapshot_footer(Path(source), p["grid"])
        if footer_cache[source]["time_tg"] != task["time_tg"] or task["time_tg"] != snapshot["time_tg"]:
            raise ValueError("Task time differs from live snapshot footer")
        if "metadata" in snapshot and footer_cache[source] != snapshot["metadata"]:
            raise ValueError("Live snapshot footer differs from plan metadata")
        if snapshot.get("source_sha256") not in (None, common.sha256(Path(source))):
            raise ValueError("Live snapshot bytes differ from plan inventory hash")
        enabled = task["kind"] == "response"
        alpha = task["alpha"] if enabled else None
        if enabled and alpha not in config["shared"]["alphas"]:
            raise ValueError("Task alpha is not in production config")
        resolution = p["response"]["resolution"]
        if not enabled and task["resolution"] not in {resolution, p["background"]["resolution"]}:
            raise ValueError("Baseline resolution differs from production config")
        if enabled and task["resolution"] != resolution:
            raise ValueError("Response resolution differs from production config")
        phase = task["parameters"]["AXION_PHASE"]
        expected_omega = float(omega(alpha, config["shared"]["principal_n"])) if enabled else 0.
        norm = p["response"]["reference_norm"] if enabled else 0.
        table_path = None
        if enabled:
            table_id = p["alpha_table"][canonical_alpha(alpha)]
            table = plan["radial_tables"][table_id]
            table_path = str(Path(table["path"]).resolve(strict=True))
            if (task["table_id"] != table_id or task["table_path"] != table_path
                    or task.get("table_spin") != p["radial_tables"][table_id]["spin"]):
                raise ValueError("Task table identity/spin differs from production plan")
            if table.get("sha256") not in (None, common.sha256(Path(table_path))):
                raise ValueError("Live radial table differs from plan inventory hash")
            if table_path not in table_cache:
                table_cache[table_path] = validate_radial_table(Path(table_path), [])
            if not table_cache[table_path]["alpha_min"] <= alpha <= table_cache[table_path]["alpha_max"]:
                raise ValueError("Live radial table does not cover task alpha")
        expected = parameter_values(_parameters(config, task["resolution"], alpha, expected_omega,
                                                norm, phase, table_path))
        if not enabled:
            expected["AXION_TABLE"] = "disabled"
        if expected != task["parameters"]:
            raise ValueError("Task model parameters differ from production plan")


def _json(path, record):
    common.require_new_output_file(path)
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def _kill_group(process):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


class _JobCancelled(BaseException):
    """Internal cancellation; the scheduler propagates the original cause."""


class _Processes:
    """Invocation-local registry, including the spawn/cancel race."""

    def __init__(self):
        self.lock = threading.Lock()
        self.processes = set()
        self.cancelled = threading.Event()
        self.cause = None

    def check(self):
        if self.cancelled.is_set():
            raise _JobCancelled()

    def cancel(self, cause):
        with self.lock:
            if self.cause is None:
                self.cause = cause
            self.cancelled.set()
            processes = list(self.processes)
        # Do not hold the registry lock while waiting for children.
        for process in processes:
            _kill_group(process)


_EXECUTION = ContextVar("producer_job_execution", default=None)


def _observed_cpu_set(pid):
    """OpenMP binds each TID separately; the leader mask may be one core."""
    cpus = set()
    for entry in (Path("/proc") / str(pid) / "task").glob("[0-9]*"):
        try:
            cpus.update(os.sched_getaffinity(int(entry.name)))
        except ProcessLookupError:
            # OpenMP threads can disappear between enumeration and observation.
            continue
    return sorted(cpus) if cpus else None


def _execute(command, cwd, env, stdout, stderr, timeout):
    """Keep the six-argument API used by measurement wrappers and fake jobs."""
    execution = _EXECUTION.get()
    registry, receipt = execution if execution is not None else (_Processes(), {})
    started = time.perf_counter()
    with stdout.open("xb") as out, stderr.open("xb") as err:
        with registry.lock:
            registry.check()
            process = subprocess.Popen(command, cwd=cwd, env=env, stdout=out, stderr=err,
                                       shell=False, start_new_session=True)
            registry.processes.add(process)
        receipt.update(pid=process.pid, process_group=process.pid)
        try:
            while True:
                registry.check()
                if hasattr(os, "sched_getaffinity"):
                    try:
                        observed = _observed_cpu_set(process.pid)
                        requested = receipt.get("cpu_set")
                        # taskset runs after Popen; OpenMP then narrows each TID
                        # mask. Ignore the inherited pre-taskset broad mask, but
                        # accept an observed union that is a requested subset.
                        if observed is not None and (requested is None or
                                                     set(observed).issubset(requested)):
                            receipt["actual_cpu_set"] = observed
                    except ProcessLookupError:
                        pass
                remaining = None if timeout is None else timeout - (time.perf_counter() - started)
                if remaining is not None and remaining <= 0:
                    _kill_group(process)
                    receipt.update(returncode=process.returncode, timed_out=True)
                    return process.returncode, True
                try:
                    code = process.wait(timeout=.1 if remaining is None else min(.1, remaining))
                    receipt["returncode"] = code
                    registry.check()
                    return code, False
                except subprocess.TimeoutExpired:
                    continue
        except BaseException as exc:
            if not isinstance(exc, Exception):
                registry.cancel(exc)
            _kill_group(process)
            receipt["returncode"] = process.returncode
            raise
        finally:
            receipt["process_wall_seconds"] = time.perf_counter() - started
            with registry.lock:
                registry.processes.discard(process)


def _assert_sources_unchanged(paths):
    # verification_session's sha256 validates cached stat identity, rather than
    # rereading 97 GB. It also validates identity before/after the initial read.
    for path in paths:
        common.sha256(path)


def _run_jobs(plan, build_dir, run_root, *, max_jobs, workers, threads, cpu_affinity,
              cpu_list, timeout_seconds):
    from .affinity import allocate_cpu_sets, launch_prefix

    workers = _integer(workers, "workers", 1)
    threads = _integer(threads, "threads", 1)
    cpu_sets = allocate_cpu_sets(workers, threads, mode=cpu_affinity, cpu_list=cpu_list)
    if len(cpu_sets) != workers:
        raise ValueError("CPU allocation must provide one set per worker")
    allocated = set()
    for cpu_set in cpu_sets:
        if cpu_set is None:
            if cpu_affinity != "none":
                raise ValueError("Explicit affinity requires a CPU set for every worker")
            continue
        if len(cpu_set) != threads or len(set(cpu_set)) != threads or allocated.intersection(cpu_set):
            raise ValueError("Worker CPU sets must have threads distinct, nonoverlapping CPUs")
        allocated.update(cpu_set)
    prefixes = [launch_prefix(cpu_set) if cpu_set is not None else [] for cpu_set in cpu_sets]
    if max_jobs is not None:
        _integer(max_jobs, "max_jobs")
    if timeout_seconds is not None and (
        isinstance(timeout_seconds, bool) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0
    ):
        raise ValueError("timeout_seconds must be positive and finite")
    if not isinstance(plan, dict) or not isinstance(plan.get("tasks"), list):
        raise ValueError("plan must contain a tasks list")
    tasks = plan["tasks"] if max_jobs is None else plan["tasks"][:max_jobs]
    prepared = [_prepare_task(task) for task in tasks]
    if len({t[0]["id"] for t in prepared}) != len(prepared):
        raise ValueError("Duplicate task ids")
    if "config" in plan and prepared:
        from .affinity import check_memory
        # Tasks-only plans are synthetic API fixtures. Real production plans
        # must also be safe when invoked directly, bypassing the CLI guard.
        check_memory(min(workers, len(prepared)))
    common.reject_directory_symlinks(build_dir)
    common.reject_directory_symlinks(run_root)
    build_dir, run_root = Path(build_dir).resolve(strict=True), Path(run_root).resolve()
    if run_root == build_dir or run_root.is_relative_to(build_dir):
        raise ValueError("Run container must not be inside the immutable build")
    for parent in [run_root, *run_root.parents]:
        for name in ("manifest.json",):
            marker = parent / name
            if not marker.exists() and not marker.is_symlink():
                continue
            # The analysis output root has a descriptive container manifest,
            # not a completed/replaceable scientific stage. Its production child
            # is an intended location, but never write in an actual stage tree.
            if marker.is_symlink():
                raise ValueError("Refusing an output stage symlink ancestor")
            record = json.loads(marker.read_text())
            container = (name == "manifest.json" and parent != run_root
                         and record.get("layout") == "canonical-run-v1"
                         and record.get("complete") is not True
                         and "schema" not in record and "stage" not in record
                         and "outputs" not in record)
            if not container:
                raise ValueError("Refusing to use an existing output stage as a run container/ancestor")
    manifest = common.verify_stage(build_dir, recursive=True)
    binary, grid, header = _build_products(build_dir, manifest)
    constants = read_core_constants(header)
    _validate_plan_binding(plan, prepared, manifest, build_dir, grid)
    base_inputs = [common.manifest_path(build_dir), binary, grid, header]
    closure = list(common.source_closure(base_inputs + [p for _, _, paths in prepared for p in paths]))
    write_roots = [run_root / "jobs", run_root / "attempts"]
    for source in closure:
        for write_root in write_roots:
            if source == write_root or source.is_relative_to(write_root) or write_root.is_relative_to(source):
                raise ValueError(f"Job/attempt root overlaps a consumed source: {source}")
        if run_root.is_relative_to(source.parent) and not source.is_relative_to(build_dir):
            raise ValueError(f"Run root is inside a consumed source directory: {source.parent}")
    for name in ("jobs", "attempts"):
        common.reject_directory_symlinks(run_root / name)
    # Reject shared core/fluid/table directory ancestry even where the concrete
    # input filename is not itself an ancestor. CLI performs broader root checks.
    for key in ("source_roots", "input_roots"):
        roots = plan.get(key, [])
        if isinstance(roots, dict):
            roots = roots.values()
        for root in roots:
            root = Path(root).resolve(strict=True)
            if root == run_root or root.is_relative_to(run_root) or run_root.is_relative_to(root):
                raise ValueError(f"Run root overlaps protected input root: {root}")
    # Prehash the shared closure on the main thread before copying its session
    # into workers. Per-job checks must not scan every other task's raw inputs.
    _assert_sources_unchanged(closure)
    common.runtime_identity()
    base_closure = list(common.source_closure(base_inputs))
    env = os.environ.copy()
    # "none" means unbound, not an unsigned policy inherited from the caller.
    # Never mutate os.environ: other invocations may run in this interpreter.
    for key in ("OMP_PROC_BIND", "OMP_PLACES", "GOMP_CPU_AFFINITY"):
        env.pop(key, None)
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    if cpu_affinity != "none":
        env.update(OMP_PROC_BIND="close", OMP_PLACES="cores")
    registry = _Processes()
    results = [None] * len(prepared)
    # Deterministic slot assignment preserves pinned cache identity on resume:
    # task i always uses slot i % workers, regardless of completion timing.
    next_indices = list(range(workers))
    active = {}
    peak = 0
    started = time.perf_counter()
    executor = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="raptor-job")

    def submit(slot):
        nonlocal peak
        registry.check()
        index = next_indices[slot]
        next_indices[slot] += workers
        # A Context cannot be entered concurrently. Each task gets a fresh copy
        # of the main session, including an independent stage-temporary stack.
        context = copy_context()
        future = executor.submit(context.run, _run_job, prepared[index], build_dir, run_root,
                                 binary, grid, header, constants, base_inputs, base_closure,
                                 threads, cpu_affinity, cpu_sets[slot], prefixes[slot], slot,
                                 env, timeout_seconds, registry)
        active[future] = (index, slot)
        peak = max(peak, len(active))

    try:
        for slot in range(min(workers, len(prepared))):
            submit(slot)
        while active:
            done, _ = wait(active, return_when=FIRST_COMPLETED)
            # Collect every observed failure/cancellation before launching more.
            free = []
            for future in sorted(done, key=lambda f: active[f][0]):
                index, slot = active.pop(future)
                results[index] = future.result()
                free.append(slot)
            for slot in sorted(free):
                if next_indices[slot] < len(prepared):
                    submit(slot)
    except BaseException as exc:
        registry.cancel(exc)
        for future in active:
            future.cancel()
        raise registry.cause
    finally:
        # Even on main-thread KeyboardInterrupt, all worker receipts and stage
        # cleanup finish before verification_session checks its shared records.
        executor.shutdown(wait=True, cancel_futures=True)
    counts = {status: sum(r["status"] == status for r in results)
              for status in ("complete", "cached", "failed")}
    return {"schema": 1, "run_root": str(run_root), "build_dir": str(build_dir), "jobs": results,
            "counts": counts, "selected_tasks": len(prepared), "planned_tasks": len(plan["tasks"]),
            "workers": workers, "threads": threads, "cpu_affinity": cpu_affinity,
            "cpu_sets": cpu_sets, "cpu_list": cpu_list,
            "scheduling": {"policy": "bounded-deterministic-slots", "peak_in_flight": peak,
                           "wall_seconds": time.perf_counter() - started},
            "complete": counts["failed"] == 0, "process_run": counts["complete"] > 0,
            "scientific_validation": False}


def _run_job(prepared, build_dir, run_root, binary, grid, header, constants, base_inputs,
             base_closure, threads, cpu_affinity, cpu_set, prefix, slot, env,
             timeout_seconds, registry):
    task, config, task_inputs = prepared
    destination = run_root / "jobs" / task["id"]
    attempt_parent = run_root / "attempts" / task["id"]
    task.update(core_constants=constants, core_constants_source=str(header))
    inputs = base_inputs + task_inputs
    closure = base_closure + list(common.source_closure(task_inputs))
    config.update(threads=threads, output_adapter_schema=1)
    if cpu_affinity != "none":
        # Worker count is scheduling only; CPU mask/binding can affect floating
        # point execution and is conservatively part of the job cache identity.
        config.update(cpu_affinity=cpu_affinity, cpu_set=cpu_set,
                      omp_proc_bind="close", omp_places="cores")
    started = time.perf_counter()
    attempt = temporary = receipt = None
    execution_token = None
    result = {"id": task["id"], "directory": str(destination), "worker_slot": slot,
              "cpu_set": cpu_set, "threads": threads}
    try:
        registry.check()
        common.reject_directory_symlinks(attempt_parent)
        with common.stage_directory(destination, "producer-job", config, inputs) as temporary:
            registry.check()
            if temporary is None:
                # A generic complete manifest alone is insufficient.
                common.require_declared_outputs(destination, [destination / "success.json",
                                                             destination / "output.json"])
                success = json.loads((destination / "success.json").read_text())
                if (success.get("schema") != 1 or success.get("returncode") != 0
                        or success.get("timed_out") is not False
                        or success.get("output_validated") is not True):
                    raise ValueError("Cached job lacks a successful process receipt")
                validated = validate_job_output(destination, task)
                common.require_declared_outputs(destination, [destination / validated["hdf5"],
                                                             destination / validated["spectrum"],
                                                             destination / "model.in",
                                                             destination / "grid_mks.in",
                                                             destination / "task.json"])
                if (destination / "model.in").read_text() != serialize_parameters(task["parameters"]):
                    raise ValueError("Cached job model.in differs from task")
                if json.loads((destination / "task.json").read_text()) != config:
                    raise ValueError("Cached task receipt differs from task config")
                if common.sha256(destination / "grid_mks.in") != common.sha256(grid):
                    raise ValueError("Cached job grid differs from declared build grid")
                if json.loads((destination / "output.json").read_text()) != validated:
                    raise ValueError("Cached output validation metadata differs")
                _assert_sources_unchanged(closure)
                result.update(status="cached")
            else:
                attempt_parent.mkdir(parents=True, exist_ok=True)
                attempt = Path(tempfile.mkdtemp(prefix="attempt-", dir=attempt_parent))
                command = prefix + [str(binary), "model.in", task["source_path"], str(task["output_slot"])]
                receipt = {"schema": 1, "task_id": task["id"], "command": command,
                           "started_unix": time.time(), "threads": threads, "worker_slot": slot,
                           "cpu_affinity": cpu_affinity, "cpu_set": cpu_set, "actual_cpu_set": None,
                           "omp_proc_bind": env.get("OMP_PROC_BIND"),
                           "omp_places": env.get("OMP_PLACES"),
                           "timeout_seconds": timeout_seconds, "inputs": common.describe_inputs(inputs),
                           "returncode": None, "timed_out": False, "complete": False}
                _json(attempt / "input.json", receipt)
                execution_token = _EXECUTION.set((registry, receipt))
                try:
                    _json(temporary / "task.json", config)
                    (temporary / "model.in").write_text(serialize_parameters(task["parameters"]))
                    shutil.copyfile(grid, temporary / "grid_mks.in")
                    _assert_sources_unchanged(closure)
                    registry.check()
                    code, timed_out = _execute(command, temporary, env, attempt / "stdout.log",
                                               attempt / "stderr.log", timeout_seconds)
                    receipt.update(returncode=code, timed_out=timed_out)
                    registry.check()
                    if timed_out:
                        raise RuntimeError(f"RAPTOR timed out after {timeout_seconds} seconds")
                    if code != 0:
                        raise RuntimeError(f"RAPTOR exited with returncode {code}")
                    validated = validate_job_output(temporary, task)
                    _json(temporary / "output.json", validated)
                    _json(temporary / "success.json", {"schema": 1, "returncode": 0, "timed_out": False,
                                                       "output_validated": True})
                    shutil.copyfile(attempt / "stdout.log", temporary / "stdout.log")
                    shutil.copyfile(attempt / "stderr.log", temporary / "stderr.log")
                    _assert_sources_unchanged(closure)
                    registry.check()
                except BaseException as exc:
                    if not isinstance(exc, Exception):
                        registry.cancel(exc)
                    # Preserve partial HDF and exact model input before common
                    # cleans the failed transaction. Never mark this tree complete.
                    shutil.copytree(temporary, attempt / "work", symlinks=True)
                    raise
        if receipt is not None:
            receipt.update(complete=True)
            _finish_receipt(attempt, receipt, started)
            result.update(status="complete", attempt=str(attempt), returncode=0,
                          actual_cpu_set=receipt["actual_cpu_set"],
                          process_wall_seconds=receipt.get("process_wall_seconds"),
                          pid=receipt.get("pid"))
    except Exception as exc:
        if receipt is not None:
            receipt.update(error=f"{type(exc).__name__}: {exc}")
            _finish_receipt(attempt, receipt, started)
        result.update(status="failed", attempt=None if attempt is None else str(attempt),
                      returncode=None if receipt is None else receipt["returncode"],
                      error=f"{type(exc).__name__}: {exc}")
    except BaseException as exc:
        registry.cancel(exc)
        if receipt is not None:
            receipt.update(error=f"{type(exc).__name__}: {exc}", cancelled=True)
            _finish_receipt(attempt, receipt, started)
        raise
    finally:
        if execution_token is not None:
            _EXECUTION.reset(execution_token)
    result["wall_seconds"] = time.perf_counter() - started
    return result


def _finish_receipt(attempt, receipt, started):
    receipt.update(finished_unix=time.time(), wall_seconds=time.perf_counter() - started)
    _json(attempt / "result.json", receipt)


def run_jobs(plan: dict, build_dir: Path, run_root: Path, *, max_jobs: int | None = None,
             workers: int = 1, threads: int = 2, cpu_affinity: str = "none",
             cpu_list: str | None = None, timeout_seconds: float | None = None) -> dict:
    """Run the authorized task prefix with bounded concurrency and shared hashes.

    No real data is copied, no shell is used, no old successful stage is rebuilt.
    ``max_jobs`` selects tasks, not workers. Operational failures are returned
    independently per task; invalid plans/root/build fail before launch.
    Cancellation kills all active process groups, saves attempts, joins workers,
    then propagates. Each logical worker slot owns a fixed disjoint CPU set.
    """
    with common.verification_session():
        return _run_jobs(plan, build_dir, run_root, max_jobs=max_jobs, workers=workers,
                         threads=threads, cpu_affinity=cpu_affinity, cpu_list=cpu_list,
                         timeout_seconds=timeout_seconds)
