"""Producer options: CLI > exports > project .env.raptor > project .env data."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import _bootstrap  # noqa: F401
from eht_axion.pipeline.common import reject_directory_symlinks
from eht_axion.producer.config import read_production_config
from _environment import add_environment_options, resolve_environment

ROOT = Path(__file__).resolve().parents[1]


def options(kind, argv=None):
    parser = argparse.ArgumentParser(
        description=f"{kind}: official RAPTOR build and Python task orchestration")
    paths = {"config": "EHT_PRODUCTION_CONFIG", "source_root": "EHT_RAPTOR_SOURCE_ROOT",
             "run_root": "EHT_PRODUCTION_ROOT", "fluid_root": "EHT_GRMHD_ROOT",
             "radial_low": "EHT_RADIAL_LOW", "radial_high": "EHT_RADIAL_HIGH"}
    for key, env in paths.items():
        parser.add_argument("--" + key.replace("_", "-"), type=Path,
                            help=f"Overrides ${env}; relative paths use current working directory")
    parser.add_argument("--run-id",
                        help="New run name under EHT_OUTPUT_ROOT/production; overrides EHT_RUN_ID")
    parser.add_argument("--source-lock", type=Path, default=ROOT / "raptor-source.lock.json")
    if kind == "generate":
        parser.add_argument("--plan", action="store_true", help="Validate and list tasks; never launch C")
        parser.add_argument("--convert-only", action="store_true",
                            help="Deliver already complete jobs; no launch")
        parser.add_argument("--max-jobs", type=int,
                            help="Run at most this many planned tasks; incomplete data not delivered")
        parser.add_argument("--threads", type=int)
        parser.add_argument("--timeout-seconds", type=float)
    add_environment_options(parser)
    args = parser.parse_args(argv)
    try:
        environment = resolve_environment(args, producer=True)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    for key, env in paths.items():
        if getattr(args, key) is None and environment.get(env):
            setattr(args, key, Path(environment[env]))
    if args.config is None:
        args.config = ROOT / "configs/production.yaml"
    if args.source_root is None:
        args.source_root = ROOT.parent / "raptor"
    if args.run_root is None:
        run_id = args.run_id or environment.get("EHT_RUN_ID")
        if not run_id or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", run_id):
            parser.error("Provide --run-root/EHT_PRODUCTION_ROOT, or a safe --run-id/EHT_RUN_ID")
        output = Path(environment.get("EHT_OUTPUT_ROOT", str(ROOT.parent / "eht-pipeline-runs")))
        args.run_root = output / "production" / run_id
    try:
        reject_directory_symlinks(args.run_root)
        args.run_root = args.run_root.resolve()
        args.source_root = args.source_root.resolve(strict=True)
        args.config = args.config.resolve(strict=True)
        config = read_production_config(args.config)
        lock = json.loads(args.source_lock.read_text())
        if lock.get("repository") != "https://github.com/ProjectEHT-PKU/raptor":
            raise ValueError("Unexpected RAPTOR source lock repository")
        args.expected_commit = lock["commit"]
        protected = [ROOT, args.source_root]
        if environment.get("EHT_INPUT_ROOT"):
            protected.append(Path(environment["EHT_INPUT_ROOT"]).resolve(strict=True))
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
        if (args.run_root / "manifest.json").exists() or (args.run_root / "stage.json").exists():
            raise ValueError("Run root must be a container, not an existing stage")
        if kind == "generate":
            if args.plan and args.convert_only:
                raise ValueError("--plan and --convert-only are exclusive")
            for key in ("fluid_root", "radial_low", "radial_high"):
                if getattr(args, key) is None:
                    raise ValueError(f"Provide --{key.replace('_', '-')} or ${paths[key]}")
            if args.threads is None:
                args.threads = int(environment.get("EHT_RAPTOR_THREADS", "2"))
            if args.timeout_seconds is None and environment.get("EHT_RAPTOR_TIMEOUT_SECONDS"):
                args.timeout_seconds = float(environment["EHT_RAPTOR_TIMEOUT_SECONDS"])
            if args.threads < 1 or (args.max_jobs is not None and args.max_jobs < 1):
                raise ValueError("threads and max-jobs must be positive")
    except (ValueError, OSError, KeyError) as error:
        parser.error(str(error))
    args.config_data = config
    return args


def save_request(args, plan=None):
    """Append requests, never overwrite science configuration or old run results."""
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
