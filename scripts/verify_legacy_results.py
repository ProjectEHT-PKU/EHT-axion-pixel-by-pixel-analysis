#!/usr/bin/env python3
"""Verify that the modular legacy path reproduces all six saved Notebook curves."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "audit" / "results" / "legacy_expected"
REFERENCE = ROOT / "Results" / "data"


def main() -> None:
    report = {"curves": {}}
    for dimension in ("1d", "2d"):
        for mode in ("std", "cor", "cor_grmhd"):
            generated_path = GENERATED / f"{mode}_{dimension}.npz"
            reference_name = f"{mode}_1d.npz" if dimension == "1d" else f"{mode}.npz"
            with (
                np.load(generated_path) as generated,
                np.load(REFERENCE / reference_name) as reference,
            ):
                alpha_difference = float(np.max(np.abs(generated["alphas"] - reference["alphas"])))
                cut_difference = float(np.max(np.abs(generated["cuts"] - reference["cuts"])))
            report["curves"][f"{mode}_{dimension}"] = {
                "maximum_alpha_difference": alpha_difference,
                "maximum_cut_difference": cut_difference,
                "passed": alpha_difference < 1e-14 and cut_difference < 1e-12,
            }
    report["all_passed"] = all(item["passed"] for item in report["curves"].values())
    output = ROOT / "audit" / "legacy_reproduction.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_passed"]:
        raise SystemExit("legacy reproduction failed")


if __name__ == "__main__":
    main()
