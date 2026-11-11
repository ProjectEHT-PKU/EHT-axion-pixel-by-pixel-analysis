"""Machine-local settings of the scripts: the dotenv file and the three path roots.

Every machine-local location derives from one of three roots (CLI > exported environment >
the repository's ``.env``; docs/en/configuration.md#machine-environment-variables-env):

* ``EHT_RAPTOR_SOURCE_ROOT`` (``--source-root`` / ``--raptor-source-root``): the RAPTOR Git
  checkout pinned by ``raptor.source`` in configs/experiment.yaml.
* ``EHT_INPUT_ROOT`` (``--input-root``): read-only external inputs (GRMHD snapshots, radial
  tables, the official EHT data and script) in the layout ``INPUT_LAYOUT``.
* ``EHT_DATA_ROOT`` (``--data-root``): every generated product, in the layout ``LAYOUT``
  (docs/en/data-layout.md#default-layout).

Dotenv data are read explicitly (never executed, never searched in the caller's cwd). An
explicit per-path option (``--...-root``, ``--radial-low``, ...) always wins over the default
derived from a root. Without either, the option is an error.
"""
from __future__ import annotations

from io import StringIO
import os
from pathlib import Path

from dotenv import dotenv_values
from dotenv.parser import parse_stream

ROOT = Path(__file__).resolve().parents[1]
RAPTOR_SOURCE_ROOT_KEY = "EHT_RAPTOR_SOURCE_ROOT"
INPUT_ROOT_KEY = "EHT_INPUT_ROOT"
DATA_ROOT_KEY = "EHT_DATA_ROOT"
# Machine-local locations only; scientific values live in configs/experiment.yaml.
PATH_KEYS = {RAPTOR_SOURCE_ROOT_KEY, INPUT_ROOT_KEY, DATA_ROOT_KEY}
# Optional execution settings of scripts/raptor/generate.py.
VALUE_KEYS = {
    "EHT_RAPTOR_THREADS", "EHT_RAPTOR_TIMEOUT_SECONDS", "EHT_RAPTOR_WORKERS",
    "EHT_RAPTOR_CPU_AFFINITY", "EHT_RAPTOR_CPU_LIST",
}


def add_environment_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--env-file", type=Path,
                       help="Read this dotenv instead of .env; file-relative paths use project root")
    group.add_argument("--no-env-file", action="store_true",
                       help="Disable dotenv loading; retain explicit CLI and exported environment values")


def resolve_environment(args):
    """CLI > exported environment > the repository's .env (or --env-file).

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
            files = [ROOT / ".env"]
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


# Relative locations of the generated products below the data root; "{day}" is an EHT
# observing day.
LAYOUT = {
    "raptor": "raptor",
    "raptor.plan": "raptor/plan.json",
    "raptor.build": "raptor/build",
    "raptor.jobs": "raptor/jobs",
    "raptor.images": "raptor/images",
    "eht.intensity": "eht/intensity/day-{day}",
    "eht.fiducial": "eht/fiducial/day-{day}",
    "eht.refine": "eht/refine/day-{day}",
    "eht.sample": "eht/sample",
    "eht.ensemble": "eht/ensemble",
    "analysis.expected_signal": "analysis/expected-signal",
    "analysis.limits": "analysis/limits",
    "analysis.check_grid": "analysis/check-grid",
    "analysis.check_prior": "analysis/check-prior",
    "analysis.check_covariance": "analysis/check-covariance",
    "analysis.summary": "analysis/summary",
    "analysis.paper": "analysis/paper",
}

# Relative locations of the read-only external inputs below the input root. The radial
# tables are keyed by the raptor.radial_tables ids of configs/experiment.yaml; the uvfits
# file names and SHA-256 are eht_imaging.inputs of configs/experiment.yaml.
INPUT_LAYOUT = {
    "grmhd": "grmhd",
    "radial.low": "radial/result_a=0.9375_l1m1.txt",
    "radial.high": "radial/result_a=0.99_l1m1.txt",
    "eht.uvfits": "eht",
    "eht.official_script": "eht/eht-imaging_pipeline.py",
}


def add_data_root_option(parser, *, environment_options=True):
    parser.add_argument("--data-root", type=Path,
                        help=f"Data root of the default layout; overrides ${DATA_ROOT_KEY}")
    if environment_options:
        add_environment_options(parser)


def add_input_root_option(parser):
    parser.add_argument("--input-root", type=Path,
                        help=f"Root of the read-only external inputs; overrides ${INPUT_ROOT_KEY}")


def _root(args, attribute, key, environment):
    value = getattr(args, attribute, None)
    if value is not None:
        return value
    if environment is None:
        environment = resolve_environment(args)
    value = environment.get(key)
    return Path(value) if value else None


def data_root(args, environment=None):
    """--data-root, else EHT_DATA_ROOT from the resolved environment, else None."""
    return _root(args, "data_root", DATA_ROOT_KEY, environment)


def input_root(args, environment=None):
    """--input-root, else EHT_INPUT_ROOT from the resolved environment, else None."""
    return _root(args, "input_root", INPUT_ROOT_KEY, environment)


def raptor_source_root(args, environment=None):
    """--source-root / --raptor-source-root, else EHT_RAPTOR_SOURCE_ROOT, else None."""
    return _root(args, "source_root", RAPTOR_SOURCE_ROOT_KEY, environment)


def layout_path(root, key, **fields):
    return Path(root) / LAYOUT[key].format(**fields)


def input_path(root, key):
    return Path(root) / INPUT_LAYOUT[key]


def _fill(parser, args, defaults, root_of, option, key, path_of, environment):
    missing = [name for name in defaults if getattr(args, name) is None]
    if not missing:
        return None
    try:
        root = root_of(args, environment)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if root is None:
        options = ", ".join("--" + name.replace("_", "-") for name in missing)
        parser.error(f"Provide {options}, or {option} / ${key}")
    for name in missing:
        setattr(args, name, path_of(root, defaults[name]))
    return root


def fill_defaults(parser, args, defaults, environment=None):
    """Set each unset option ``attribute`` to ``layout_path(data root, key, **fields)``.

    ``defaults`` maps an argparse attribute to ``(layout key, fields)``. An option
    given on the command line is kept. Returns the data root (or None).
    """
    return _fill(parser, args, defaults, data_root, "--data-root", DATA_ROOT_KEY,
                 lambda root, item: layout_path(root, item[0], **item[1]), environment)


def fill_input_defaults(parser, args, defaults, environment=None):
    """Set each unset option ``attribute`` to ``input_path(input root, key)``.

    ``defaults`` maps an argparse attribute to an ``INPUT_LAYOUT`` key. An option given
    on the command line is kept. Returns the input root (or None).
    """
    return _fill(parser, args, defaults, input_root, "--input-root", INPUT_ROOT_KEY, input_path,
                 environment)
