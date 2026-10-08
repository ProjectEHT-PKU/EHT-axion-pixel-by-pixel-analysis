"""Process trusted EHT image ensembles with fixed relocation and summaries.

The reusable stage leaves transaction and cache management to its caller.
"""
from __future__ import annotations

import _bootstrap  # noqa: F401

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from eht_axion.pipeline.eht import run_eht_stage


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--trust-pickle", action="store_true", help="allow executable trusted input")
    parser.add_argument("--topset-size", type=int, default=1572)
    parser.add_argument("--centering", choices=("ensemble",), default="ensemble")
    parser.add_argument("--weighting", choices=("intensity", "flux"), default="intensity")
    parser.add_argument("--inherit-2d-mask", action="store_true")
    parser.add_argument("--fov", type=float, default=32.6868)
    parser.add_argument("--radial-max", type=float)
    parser.add_argument("--window-degrees", type=int, default=10)
    parser.add_argument("--angular-bins", type=int, default=360)
    parser.add_argument("--angle-roll", type=int, default=-90, help="roll in bins; default 360-bin grid")
    parser.add_argument("--min-radius", type=float, default=2.0)
    parser.add_argument("--intensity-fraction", type=float, default=0.1)
    parser.add_argument("--no-historical-mask", dest="historical_mask", action="store_false")
    parser.add_argument("--exclude-columns", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--pixel-batch-size", type=int, default=128)
    args = parser.parse_args(argv)
    if not args.trust_pickle:
        parser.error("--trust-pickle is required: unpickling may execute arbitrary code")
    report = run_eht_stage(**vars(args))
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
