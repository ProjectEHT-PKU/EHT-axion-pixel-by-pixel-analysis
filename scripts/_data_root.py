"""Default input and output roots below one data root (docs/en/data-layout.md#default-layout).

The data root is ``--data-root`` or ``M87_AXION_DATA_ROOT`` (exported, or in ``.env``;
CLI > exported environment > dotenv). An explicit ``--...-root`` option always wins
over the default derived from the data root. Without either, the option is an error.
"""
from __future__ import annotations

from pathlib import Path

from _environment import add_environment_options, resolve_environment

DATA_ROOT_KEY = "M87_AXION_DATA_ROOT"

# Relative locations of the default layout; "{day}" is an EHT observing day.
LAYOUT = {
    "raptor_plan": "raptor/plan.json",
    "raptor_runs": "raptor/runs",
    "raptor_run": "raptor/runs/main",
    "raptor_delivered": "raptor/delivered",
    "raptor_data": "raptor/delivered/data",
    "eht_intensity": "eht/intensity/day-{day}",
    "eht_fiducial": "eht/fiducial/day-{day}",
    "eht_refine": "eht/refine/day-{day}",
    "eht_bundle": "eht/bundle",
    "eht_sample": "eht/sample",
    "eht_ensemble": "eht/ensemble",
    "analysis_templates": "analysis/templates",
    "analysis_limits": "analysis/limits",
    "analysis_check_grid": "analysis/check-grid",
    "analysis_check_prior": "analysis/check-prior",
    "analysis_check_covariance": "analysis/check-covariance",
    "analysis_summary": "analysis/summary",
    "analysis_paper": "analysis/paper",
}


def add_data_root_option(parser, *, environment_options=True):
    parser.add_argument("--data-root", type=Path,
                        help=f"Data root of the default layout; overrides ${DATA_ROOT_KEY}")
    if environment_options:
        add_environment_options(parser)


def data_root(args, environment=None):
    """--data-root, else M87_AXION_DATA_ROOT from the resolved environment, else None."""
    if args.data_root is not None:
        return args.data_root
    if environment is None:
        environment = resolve_environment(args)
    value = environment.get(DATA_ROOT_KEY)
    return Path(value) if value else None


def layout_path(root, key, **fields):
    return Path(root) / LAYOUT[key].format(**fields)


def fill_defaults(parser, args, defaults, environment=None):
    """Set each unset option ``attribute`` to ``layout_path(data root, key, **fields)``.

    ``defaults`` maps an argparse attribute to ``(layout key, fields)``. An option
    given on the command line is kept. Returns the data root (or None).
    """
    missing = [name for name in defaults if getattr(args, name) is None]
    if not missing:
        return None
    try:
        root = data_root(args, environment)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if root is None:
        options = ", ".join("--" + name.replace("_", "-") for name in missing)
        parser.error(f"Provide {options}, or --data-root / ${DATA_ROOT_KEY}")
    for name in missing:
        key, fields = defaults[name]
        setattr(args, name, layout_path(root, key, **fields))
    return root
