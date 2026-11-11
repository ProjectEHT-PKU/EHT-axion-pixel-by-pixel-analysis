"""Command-line options of the RAPTOR scripts: CLI > exported environment > project .env.

Three roots (eht_axion/environment.py): ``--source-root`` / ``EHT_RAPTOR_SOURCE_ROOT`` (the
RAPTOR checkout), ``--input-root`` / ``EHT_INPUT_ROOT`` (snapshots and radial tables; generate
only) and ``--data-root`` / ``EHT_DATA_ROOT`` (the run is raptor/ below it: plan.json, build/,
jobs/, images/).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from ..environment import (DATA_ROOT_KEY, INPUT_ROOT_KEY, RAPTOR_SOURCE_ROOT_KEY, ROOT,
                           add_data_root_option, add_input_root_option, data_root, input_path,
                           input_root, layout_path, raptor_source_root, resolve_environment)
from ..experiment import DEFAULT_PATH, read_experiment
from ..provenance import reject_directory_symlinks

# generate.py inputs: option attribute -> INPUT_LAYOUT key below the input root
INPUTS = {"fluid_root": "grmhd", "radial_low": "radial.low", "radial_high": "radial.high"}


def options(kind, argv=None):
    parser = argparse.ArgumentParser(
        description=f"{kind}: official RAPTOR build and Python task orchestration")
    parser.add_argument("--config", type=Path, default=DEFAULT_PATH,
                        help="Experiment parameter file (default: configs/experiment.yaml)")
    parser.add_argument("--source-root", "--raptor-source-root", dest="source_root", type=Path,
                        help=f"RAPTOR Git checkout; overrides ${RAPTOR_SOURCE_ROOT_KEY}")
    if kind == "generate":
        add_input_root_option(parser)
        for key, layout in INPUTS.items():
            parser.add_argument("--" + key.replace("_", "-"), type=Path,
                                help=f"default: {layout} layout entry below --input-root")
        parser.add_argument("--plan", action="store_true", help="Validate and list tasks; never launch C")
        parser.add_argument("--convert-only", action="store_true",
                            help="Write the image cubes of already complete jobs; no launch")
        parser.add_argument("--max-jobs", type=int,
                            help="Run at most this many planned tasks; image cubes need all jobs")
        parser.add_argument("--threads", type=int, help="OpenMP threads per RAPTOR process")
        parser.add_argument("--workers", type=int, help="Concurrent RAPTOR processes, independent of threads")
        parser.add_argument("--cpu-affinity", choices=["none", "physical", "logical"],
                            help="Disjoint per-worker CPU sets; physical avoids SMT siblings")
        parser.add_argument("--cpu-list", help="Restrict allocation to CPU IDs/ranges, e.g. 0-47")
        parser.add_argument("--timeout-seconds", type=float)
    add_data_root_option(parser)
    args = parser.parse_args(argv)
    try:
        environment = resolve_environment(args)
        args.source_root = raptor_source_root(args, environment)
        root = data_root(args, environment)
        args.run_root = None if root is None else layout_path(root, "raptor")
        if kind == "generate":
            missing = [key for key in INPUTS if getattr(args, key) is None]
            root = input_root(args, environment) if missing else None
            for key in missing:
                if root is not None:
                    setattr(args, key, input_path(root, INPUTS[key]))
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if args.source_root is None:
        parser.error(f"Provide --source-root or ${RAPTOR_SOURCE_ROOT_KEY}")
    if args.run_root is None:
        parser.error(f"Provide --data-root or ${DATA_ROOT_KEY}")
    if kind == "generate":
        missing = ", ".join("--" + key.replace("_", "-") for key in INPUTS if getattr(args, key) is None)
        if missing:
            parser.error(f"Provide {missing}, or --input-root / ${INPUT_ROOT_KEY}")
    try:
        reject_directory_symlinks(args.run_root)
        args.run_root = args.run_root.resolve()
        args.source_root = args.source_root.resolve(strict=True)
        args.config = args.config.resolve(strict=True)
        config = read_experiment(args.config)
        protected = [ROOT, args.source_root]
        if kind == "generate":
            # Snapshots and radial tables only matter to generate; build never reads them.
            for name in INPUTS:
                path = getattr(args, name).resolve(strict=True)
                setattr(args, name, path)
                protected.append(path if path.is_dir() else path.parent)
        for source in protected:
            if (source == args.run_root or source.is_relative_to(args.run_root)
                    or args.run_root.is_relative_to(source)):
                raise ValueError(f"Run root overlaps protected input/code root: {source}")
        if (args.run_root / "manifest.json").exists():
            raise ValueError("Run root must be a container, not an existing stage")
        if kind == "generate":
            if args.plan and args.convert_only:
                raise ValueError("--plan and --convert-only are exclusive")
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
