"""Read explicit project dotenv data, never execute shell or search caller cwd."""
from __future__ import annotations

from io import StringIO
import os
from pathlib import Path

from dotenv import dotenv_values
from dotenv.parser import parse_stream

ROOT = Path(__file__).resolve().parents[1]
# Machine-local locations only; scientific values live in experiment.yaml.
PATH_KEYS = {
    # Data root of the default layout (docs/en/data-layout.md), read by every script
    "M87_AXION_DATA_ROOT",
    # RAPTOR simulation, read from .env and .env.raptor
    "EHT_RAPTOR_SOURCE_ROOT", "EHT_RAPTOR_RUN_ROOT",
    "EHT_GRMHD_ROOT", "EHT_RADIAL_LOW", "EHT_RADIAL_HIGH",
    # EHT imaging, read from .env
    "EHT_IMAGING_INPUT_ROOT", "EHT_IMAGING_INPUT_RECEIPT", "EHT_IMAGING_FEED_ROTATION",
    "EHT_IMAGING_OFFICIAL_SCRIPT",
}
VALUE_KEYS = {
    "EHT_RAPTOR_THREADS", "EHT_RAPTOR_TIMEOUT_SECONDS", "EHT_RAPTOR_WORKERS",
    "EHT_RAPTOR_CPU_AFFINITY", "EHT_RAPTOR_CPU_LIST",
}


def add_environment_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--env-file", type=Path,
                       help="Read this dotenv instead of defaults; file-relative paths use project root")
    group.add_argument("--no-env-file", action="store_true",
                       help="Disable dotenv loading; retain explicit CLI and exported environment values")


def resolve_environment(args, *, raptor=False):
    """CLI > exported environment > .env.raptor (raptor_* scripts) > .env.

    Return a local mapping, never contaminate os.environ between script calls.
    File-derived relative paths are anchored to this checkout, not caller cwd.
    There is no variable expansion, shell substitution or dotenv file discovery.
    """
    files = []
    if not args.no_env_file:
        if args.env_file is not None:
            if not args.env_file.is_file():
                raise ValueError(f"Dotenv file does not exist: {args.env_file}")
            files = [args.env_file]
        else:
            files = [ROOT / ".env"] + ([ROOT / ".env.raptor"] if raptor else [])
    values = {}
    for path in files:
        if not path.exists():
            continue
        if not path.is_file():
            raise ValueError(f"Dotenv must be a regular file: {path}")
        content = path.read_text(encoding="utf-8")
        bindings = list(parse_stream(StringIO(content)))
        seen = set()
        for binding in bindings:
            if binding.error:
                raise ValueError(f"Invalid dotenv data in {path}:{binding.original.line}; "
                                 "shell commands not allowed")
            if binding.key is None:
                continue
            if binding.key not in PATH_KEYS | VALUE_KEYS:
                raise ValueError(f"Unknown dotenv key {binding.key!r} in {path}; "
                                 "only the declared settings are allowed")
            if binding.key in seen:
                raise ValueError(f"Duplicate dotenv key {binding.key!r} in {path}")
            if binding.value is None:
                raise ValueError(f"Dotenv key requires '=value': {binding.key} in {path}")
            seen.add(binding.key)
        loaded = dotenv_values(stream=StringIO(content), interpolate=False)
        for key, value in loaded.items():
            if key in PATH_KEYS and value:
                p = Path(value)
                value = str(p if p.is_absolute() else ROOT / p)
            values[key] = value
    values.update(os.environ)
    return values
