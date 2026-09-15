#!/usr/bin/env python3
"""Check 60-point final limits against a 120-point coupling grid."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from eht_axion.limits import LimitOptions, calculate_limits

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "audit" / "results" / "rebuilt_corrected_expected"


def main() -> None:
    report = {"coarse_points": 60, "fine_points": 120, "curves": {}}
    for dimension in ("1d", "2d"):
        observation_dir = ROOT / "audit" / ("rebuilt_1d" if dimension == "1d" else "rebuilt_2d")
        fit_dir = ROOT / "audit" / "rebuilt_1d" if dimension == "1d" else None
        _, fine = calculate_limits(
            ROOT,
            dimension=dimension,
            mode="cor_grmhd",
            options=LimitOptions(coupling_points=120),
            observation_dir=observation_dir,
            fit_parameter_dir=fit_dir,
        )
        with np.load(BASE / f"cor_grmhd_{dimension}.npz") as data:
            coarse = data["cuts"]
        difference = fine - coarse
        report["curves"][dimension] = {
            "coarse": coarse.tolist(),
            "fine": fine.tolist(),
            "difference": difference.tolist(),
            "maximum_absolute_difference": float(np.max(np.abs(difference))),
        }
    output = ROOT / "audit" / "grid_convergence.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
