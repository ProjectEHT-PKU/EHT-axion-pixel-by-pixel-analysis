"""Research-stage options: explicit arguments override exported environment values.

.env is intentionally loaded by the user's shell, not implicitly executed by Python.
"""

from pathlib import Path
import argparse
import os

import _bootstrap  # noqa: F401
from eht_axion.pipeline.workflow import run_stage


ENVIRONMENT_OPTIONS = {
    "config": "EHT_CONFIG",
    "input_root": "EHT_INPUT_ROOT",
    "output_root": "EHT_OUTPUT_ROOT",
    "trust_pickle": "EHT_TRUST_PICKLE",
    "rebuild": "EHT_REBUILD",
    "reuse_imported": "EHT_REUSE_IMPORTED",
}


def environment_boolean(name):
    value = os.environ.get(name, "").strip().lower()
    if value in {"", "0", "false", "no", "off"}:
        return False
    if value in {"1", "true", "yes", "on"}:
        return True
    raise ValueError(f"{name} must be true/false, 1/0, yes/no or on/off")


def parse_options(stage, argv=None):
    parser = argparse.ArgumentParser(
        description=f"Run the {stage} research stage. Arguments override exported EHT_* variables."
    )
    for option in ("config", "input_root", "output_root"):
        parser.add_argument(
            "--" + option.replace("_", "-"),
            type=Path,
            help=f"Overrides ${ENVIRONMENT_OPTIONS[option]}",
        )
    for option in ("trust_pickle", "rebuild", "reuse_imported"):
        parser.add_argument(
            "--" + option.replace("_", "-"),
            action=argparse.BooleanOptionalAction,
            default=None,
            help=f"Overrides ${ENVIRONMENT_OPTIONS[option]}; --no-* explicitly disables it",
        )
    args = parser.parse_args(argv)
    for option, name in ENVIRONMENT_OPTIONS.items():
        if getattr(args, option) is not None:
            continue  # Do not even parse a stale/invalid environment value if CLI overrides it.
        if option in {"config", "input_root", "output_root"}:
            value = os.environ.get(name, "").strip()
            setattr(args, option, Path(value) if value else None)
        else:
            try:
                setattr(args, option, environment_boolean(name))
            except ValueError as exc:
                parser.error(str(exc))
    for option in ("config", "output_root"):
        if getattr(args, option) is None:
            parser.error(
                f"Set {ENVIRONMENT_OPTIONS[option]} in the shell or pass --{option.replace('_', '-')}"
            )
    if stage in {"all", "raptor", "eht", "background"} and args.input_root is None:
        parser.error("Set EHT_INPUT_ROOT in the shell or pass --input-root")
    return args


def execute(stage, argv=None):
    args = parse_options(stage, argv)
    output = run_stage(
        args.config,
        args.output_root,
        input_root=args.input_root,
        stage=stage,
        trust_pickle=args.trust_pickle,
        rebuild=args.rebuild,
        reuse_imported=args.reuse_imported,
    )
    print(f"Results: {output}")
    return output
