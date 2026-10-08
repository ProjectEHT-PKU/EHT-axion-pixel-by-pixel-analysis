"""Expert RAPTOR tools: convert, preprocess, fit, summarize and background.

These low-level writers leave orchestration, transactions and caching to callers.
This script is not an installed or unified pipeline CLI. All angles are radians.
"""
from __future__ import annotations

import _bootstrap  # noqa: F401

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from eht_axion.pipeline.raptor import (
    _integer,
    convert_hdf5,
    estimate_background,
    fit_stokes,
    preprocess_stokes,
    summarize_fits,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Parse expert-tool arguments and delegate to reusable stage functions."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    def io(name):
        sub = commands.add_parser(name)
        sub.add_argument("inputs", nargs="+", type=Path, help="explicitly ordered input paths")
        sub.add_argument("--output-dir", required=True, type=Path)
        return sub

    def aperture(sub, *, required=True):
        sub.add_argument("--fov", required=required, type=float)
        sub.add_argument("--weighting", choices=("intensity", "flux"), default="intensity")
        sub.add_argument("--radial-max", type=float)
        sub.add_argument("--window-degrees", type=float, default=10)
        sub.add_argument("--angular-bins", type=int, default=360)

    convert = io("convert")
    convert.add_argument("--mapping", required=True, type=Path,
                         help="JSON with datasets, image_shape, block_origins/pixel_shape")
    convert.add_argument("--states-per-time", type=int,
                         help="reshape ordered HDF5 paths into a time/state grid")
    convert.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    preprocess = io("preprocess")
    for name in ("sigma", "init-fov", "final-fov", "rotation-angle"):
        preprocess.add_argument("--" + name, required=True, type=float)
    preprocess.add_argument("--final-shape", type=int, default=64)
    preprocess.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    preprocess.add_argument("--state-indices", nargs="+", type=int,
                            help="process only these states, retaining the state axis")
    fit = io("fit")
    aperture(fit)
    fit.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    summarize = commands.add_parser("summarize")
    summarize.add_argument("--fit", action="append", required=True, metavar="ALPHA=PATH",
                           help="repeat alpha to supply ordered fit shards")
    summarize.add_argument("--output-dir", required=True, type=Path)
    summarize.add_argument("--minimum-samples", type=int, default=2)
    background = io("background")
    background.add_argument("--lag", type=int, default=3)
    background.add_argument("--cadence-hours", type=float, default=26.68 / 3)
    background.add_argument("--noaxion-state", default="0", help="state index, or 'none' for 4D")
    background.add_argument("--mask", type=Path)
    background.add_argument("--difference", choices=("angular", "linear"), default="angular")
    background.add_argument("--method", choices=("mle", "legacy_hist"), default="mle")
    background.add_argument("--bins", type=int, default=100)
    background.add_argument("--minimum-samples", type=int, default=2)
    background.add_argument("--domain", choices=("2d", "1d"), default="2d")
    # FOV is only needed for 1D background, not pixel differences.
    aperture(background, required=False)
    background.add_argument("--legacy-threshold", type=float, default=0.2)
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    try:
        if command == "convert":
            mapping = json.loads(args.pop("mapping").read_text(encoding="utf-8"))
            states = args.pop("states_per_time")
            if states is not None:
                states = _integer(states, "states_per_time")
                paths = args["inputs"]
                if len(paths) % states:
                    raise ValueError("HDF5 input count must be divisible by states-per-time")
                args["inputs"] = [paths[i:i + states] for i in range(0, len(paths), states)]
            result = convert_hdf5(args.pop("inputs"), **args, **mapping)
        elif command == "summarize":
            grouped = {}
            for spec in args.pop("fit"):
                alpha, path = spec.split("=", 1)
                grouped.setdefault(float(alpha), []).append(Path(path))
            result = summarize_fits(grouped, **args)
        elif command == "background":
            state = args["noaxion_state"]
            args["noaxion_state"] = None if state.lower() == "none" else int(state)
            result = estimate_background(**args)
        else:
            result = {"preprocess": preprocess_stokes, "fit": fit_stokes}[command](**args)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.error(str(error))
    print(json.dumps({name: str(path) for name, path in result.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
