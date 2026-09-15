#!/usr/bin/env python3
"""Spot-check saved GRMHD Gaussian fit products against their source arrays."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from utils.libanaly import fit_gaussian2

ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44)


def main() -> None:
    rng = np.random.default_rng(20260914)
    report = {}
    data_dir = ROOT / "GRMHD" / "data_out"
    for dimension, source_pattern, amp_file, phase_file, spatial_shape in (
        (
            "2d",
            "EVPA_cos_fit_result_alpha={alpha}.npy",
            "A_params_dict.npy",
            "phi0_params_dict.npy",
            (64, 64),
        ),
        (
            "1d",
            "EVPA_cos_fit_result_1d_10deg_alpha={alpha}.npy",
            "A_params_1d_dict.npy",
            "phi0_params_1d_dict.npy",
            (360,),
        ),
    ):
        saved_amp = np.load(data_dir / amp_file, allow_pickle=True).item()
        saved_phase = np.load(data_dir / phase_file, allow_pickle=True).item()
        checks = []
        flat_size = int(np.prod(spatial_shape))
        indices = rng.choice(flat_size, size=min(32, flat_size), replace=False)
        for alpha in ALPHAS:
            source = np.load(data_dir / source_pattern.format(alpha=alpha), mmap_mode="r")
            for flat_index in indices:
                coordinate = np.unravel_index(flat_index, spatial_shape)
                amp_values = source[(slice(None),) + coordinate + (0,)]
                phase_values = source[(slice(None),) + coordinate + (1,)]
                expected_amp = fit_gaussian2(
                    amp_values,
                    plot_failed=False,
                    threshold=np.nan if dimension == "2d" else 0.2,
                    bins=20,
                    nomean=dimension == "2d",
                )
                expected_phase = fit_gaussian2(
                    phase_values,
                    plot_failed=False,
                    angle=True,
                    lower_bound=-np.pi,
                    upper_bound=np.pi,
                    threshold=np.nan if dimension == "2d" else 0.2,
                    bins=20,
                    nomean=dimension == "2d",
                )
                actual_amp = saved_amp[alpha][(slice(None),) + coordinate]
                actual_phase = saved_phase[alpha][(slice(None),) + coordinate]
                checks.append(
                    {
                        "alpha": alpha,
                        "coordinate": [int(value) for value in coordinate],
                        "amplitude_max_abs_difference": float(
                            np.max(np.abs(expected_amp - actual_amp))
                        ),
                        "phase_max_abs_difference": float(
                            np.max(np.abs(expected_phase - actual_phase))
                        ),
                    }
                )
        report[dimension] = {
            "checks": len(checks),
            "amplitude_max_abs_difference": max(x["amplitude_max_abs_difference"] for x in checks),
            "phase_max_abs_difference": max(x["phase_max_abs_difference"] for x in checks),
            "details": checks,
        }
    output = ROOT / "audit" / "fit_product_validation.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {k: {x: y for x, y in v.items() if x != "details"} for k, v in report.items()}, indent=2
        )
    )


if __name__ == "__main__":
    main()
