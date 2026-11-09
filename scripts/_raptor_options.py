"""Command-line options of the raptor_* scripts: CLI > exports > project .env.raptor > project .env."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import _bootstrap  # noqa: F401
from eht_axion.experiment import DEFAULT_PATH, read_experiment
from eht_axion.provenance import reject_directory_symlinks
from _environment import add_environment_options, resolve_environment

ROOT = Path(__file__).resolve().parents[1]


def options(kind, argv=None):
    parser = argparse.ArgumentParser(
        description=f"{kind}: official RAPTOR build and Python task orchestration")
    parser.add_argument("--config", type=Path, default=DEFAULT_PATH,
                        help="Experiment parameter file (default: experiment.yaml)")
    paths = {"source_root": "EHT_RAPTOR_SOURCE_ROOT",
             "run_root": "EHT_RAPTOR_RUN_ROOT", "fluid_root": "EHT_GRMHD_ROOT",
             "radial_low": "EHT_RADIAL_LOW", "radial_high": "EHT_RADIAL_HIGH"}
    for key, env in paths.items():
        parser.add_argument("--" + key.replace("_", "-"), type=Path,
                            help=f"Overrides ${env}; relative paths use current working directory")
    parser.add_argument("--source-lock", type=Path, default=ROOT / "raptor-source.lock.json")
    if kind == "generate":
        parser.add_argument("--plan", action="store_true", help="Validate and list tasks; never launch C")
        parser.add_argument("--plan-output", type=Path,
                            help="With --plan: also write the frozen plan JSON (input of raptor_merge.py)")
        parser.add_argument("--convert-only", action="store_true",
                            help="Deliver already complete jobs; no launch")
        parser.add_argument("--max-jobs", type=int,
                            help="Run at most this many planned tasks; incomplete data not delivered")
        parser.add_argument("--threads", type=int, help="OpenMP threads per RAPTOR process")
        parser.add_argument("--workers", type=int, help="Concurrent RAPTOR processes, independent of threads")
        parser.add_argument("--cpu-affinity", choices=["none", "physical", "logical"],
                            help="Disjoint per-worker CPU sets; physical avoids SMT siblings")
        parser.add_argument("--cpu-list", help="Restrict allocation to CPU IDs/ranges, e.g. 0-47")
        parser.add_argument("--timeout-seconds", type=float)
    add_environment_options(parser)
    args = parser.parse_args(argv)
    try:
        environment = resolve_environment(args, raptor=True)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    for key, env in paths.items():
        if getattr(args, key) is None and environment.get(env):
            setattr(args, key, Path(environment[env]))
    for key in ("source_root", "run_root"):
        if getattr(args, key) is None:
            parser.error(f"Provide --{key.replace('_', '-')} or ${paths[key]}")
    try:
        reject_directory_symlinks(args.run_root)
        args.run_root = args.run_root.resolve()
        args.source_root = args.source_root.resolve(strict=True)
        args.config = args.config.resolve(strict=True)
        config = read_experiment(args.config)
        lock = json.loads(args.source_lock.read_text())
        if lock.get("repository") != "https://github.com/ProjectEHT-PKU/raptor":
            raise ValueError("Unexpected RAPTOR source lock repository")
        args.expected_commit = lock["commit"]
        protected = [ROOT, args.source_root]
        for name in ("fluid_root", "radial_low", "radial_high"):
            path = getattr(args, name)
            if path is not None:
                path = path.resolve(strict=True)
                setattr(args, name, path)
                protected.append(path if path.is_dir() else path.parent)
        for source in protected:
            if (source == args.run_root or source.is_relative_to(args.run_root)
                    or args.run_root.is_relative_to(source)):
                raise ValueError(f"Production run root overlaps protected input/code root: {source}")
        if (args.run_root / "manifest.json").exists():
            raise ValueError("Run root must be a container, not an existing stage")
        if kind == "generate":
            if args.plan and args.convert_only:
                raise ValueError("--plan and --convert-only are exclusive")
            for key in ("fluid_root", "radial_low", "radial_high"):
                if getattr(args, key) is None:
                    raise ValueError(f"Provide --{key.replace('_', '-')} or ${paths[key]}")
            if args.threads is None:
                args.threads = int(environment.get("EHT_RAPTOR_THREADS", "2"))
            if args.workers is None:
                args.workers = int(environment.get("EHT_RAPTOR_WORKERS", "1"))
            if args.cpu_affinity is None:
                args.cpu_affinity = environment.get("EHT_RAPTOR_CPU_AFFINITY", "none")
            if args.cpu_list is None:
                args.cpu_list = environment.get("EHT_RAPTOR_CPU_LIST") or None
            if args.timeout_seconds is None and environment.get("EHT_RAPTOR_TIMEOUT_SECONDS"):
                args.timeout_seconds = float(environment["EHT_RAPTOR_TIMEOUT_SECONDS"])
            if (args.threads < 1 or args.workers < 1
                    or (args.max_jobs is not None and args.max_jobs < 1)):
                raise ValueError("workers, threads and max-jobs must be positive")
            from eht_axion.raptor.affinity import allocate_cpu_sets

            args.cpu_sets = allocate_cpu_sets(args.workers, args.threads,
                                             args.cpu_affinity, args.cpu_list)
    except (ValueError, OSError, KeyError) as error:
        parser.error(str(error))
    args.config_data = config
    return args


def save_request(args, plan=None):
    """Append a request record; never overwrite the run's experiment file or earlier results."""
    import uuid

    args.run_root.mkdir(parents=True, exist_ok=True)
    target = args.run_root / "config.yaml"
    reject_directory_symlinks(target)
    content = args.config.read_bytes()
    if target.exists():
        if target.read_bytes() != content:
            raise ValueError("Run ID already uses a different science configuration; choose a new run")
    else:
        with target.open("xb") as stream:
            stream.write(content)
    directory = args.run_root / "requests"
    reject_directory_symlinks(directory)
    directory.mkdir(exist_ok=True)
    record = {"config": str(args.config), "run_root": str(args.run_root),
              "core_source": str(args.source_root), "core_commit": args.expected_commit}
    if plan is not None:
        record["plan"] = plan
    path = directory / f"{uuid.uuid4().hex}.json"
    with path.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return path
