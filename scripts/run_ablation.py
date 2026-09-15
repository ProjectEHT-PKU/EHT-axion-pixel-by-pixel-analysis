#!/usr/bin/env python3
"""Measure the incremental impact of each verified inference correction."""

from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

from eht_axion.limits import LimitOptions, calculate_limits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit" / "ablation_2d_cor_grmhd.json"


def main() -> None:
    base = LimitOptions(legacy=True)
    configurations = [
        ("legacy", base),
        ("unique_phase_endpoint", replace(base, duplicate_phase_endpoint=False)),
        (
            "multiply_independent_epochs",
            replace(
                base,
                duplicate_phase_endpoint=False,
                combine_epochs_by_product=True,
            ),
        ),
        (
            "epoch_specific_model_uncertainty",
            replace(
                base,
                duplicate_phase_endpoint=False,
                combine_epochs_by_product=True,
                shared_epoch_uncertainty=False,
            ),
        ),
        (
            "parameter_dependent_logdet",
            replace(
                base,
                duplicate_phase_endpoint=False,
                combine_epochs_by_product=True,
                shared_epoch_uncertainty=False,
                include_parameter_logdet=True,
            ),
        ),
        (
            "explicit_posterior_quadrature",
            replace(
                base,
                duplicate_phase_endpoint=False,
                combine_epochs_by_product=True,
                shared_epoch_uncertainty=False,
                include_parameter_logdet=True,
                posterior_method="trapezoid",
            ),
        ),
    ]
    report = {"dimension": "2d", "mode": "cor_grmhd", "steps": {}}
    previous = None
    for name, options in configurations:
        begin = time.perf_counter()
        alphas, cuts = calculate_limits(ROOT, dimension="2d", mode="cor_grmhd", options=options)
        report["steps"][name] = {
            "alphas": alphas.tolist(),
            "cuts": cuts.tolist(),
            "change_from_previous": (None if previous is None else (cuts - previous).tolist()),
            "maximum_absolute_change_from_previous": (
                None if previous is None else float(max(abs(cuts - previous)))
            ),
            "elapsed_seconds": time.perf_counter() - begin,
        }
        previous = cuts
        OUT.write_text(json.dumps(report, indent=2) + "\n")
        print(f"completed ablation step {name}", flush=True)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
