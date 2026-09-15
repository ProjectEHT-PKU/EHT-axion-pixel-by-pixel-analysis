#!/usr/bin/env python3
"""Reproduce or correct the six final EHT axion limit curves."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

from eht_axion.limits import LimitOptions, calculate_limits

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", choices=("1d", "2d"), required=True)
    parser.add_argument("--mode", choices=("std", "cor", "cor_grmhd"), required=True)
    parser.add_argument("--legacy", action="store_true", help="reproduce historical formulas")
    parser.add_argument(
        "--observed", action="store_true", help="use measured rather than null EVPA differences"
    )
    parser.add_argument("--prior", choices=("log-uniform", "uniform"), default="log-uniform")
    parser.add_argument("--coupling-points", type=int, default=60)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--observation-dir", type=Path)
    parser.add_argument("--fit-parameter-dir", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    variant = "legacy" if args.legacy else "corrected"
    observation = "observed" if args.observed else "expected"
    output_dir = args.output_dir or ROOT / "audit" / "results" / f"{variant}_{observation}"
    output_dir.mkdir(parents=True, exist_ok=True)
    options = LimitOptions(
        legacy=args.legacy,
        observed=args.observed,
        prior=args.prior,
        coupling_points=args.coupling_points,
    )
    started = time.perf_counter()
    alphas, cuts = calculate_limits(
        ROOT,
        dimension=args.dimension,
        mode=args.mode,
        options=options,
        observation_dir=args.observation_dir,
        fit_parameter_dir=args.fit_parameter_dir,
    )
    elapsed = time.perf_counter() - started
    stem = f"{args.mode}_{args.dimension}"
    np.savez(
        output_dir / f"{stem}.npz",
        alphas=alphas,
        cuts=cuts,
        variant=variant,
        observation=observation,
        prior=args.prior,
    )
    metadata = {
        "dimension": args.dimension,
        "mode": args.mode,
        "variant": variant,
        "observation": observation,
        "prior": args.prior,
        "coupling_points": args.coupling_points,
        "alphas": alphas.tolist(),
        "cuts": cuts.tolist(),
    }
    (output_dir / f"{stem}.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)
    print(f"elapsed_seconds={elapsed:.6f}; pid={os.getpid()}", flush=True)


if __name__ == "__main__":
    main()
