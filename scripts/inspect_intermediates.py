#!/usr/bin/env python3
"""Validate the intermediate products used by the final-limit calculation.

This script deliberately memory-maps large arrays and reports hashes, shapes,
finite fractions, ranges, and internal consistency checks without rewriting data.
Run it from any directory; paths are resolved relative to the repository root.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ALPHAS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44)
DAYS = (95, 96, 100, 101)


def sha256(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def array_summary(array: np.ndarray) -> dict[str, object]:
    finite = np.isfinite(array)
    values = array[finite]
    return {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "finite_fraction": float(finite.mean()),
        "min": float(values.min()) if values.size else None,
        "max": float(values.max()) if values.size else None,
        "mean": float(values.mean()) if values.size else None,
    }


def main() -> None:
    report: dict[str, object] = {"root": str(ROOT), "grmhd": {}, "eht": {}}

    for alpha in ALPHAS:
        products = {}
        for stem in (
            "EVPA_cos_fit_result_alpha={alpha}.npy",
            "EVPA_cos_fit_result_1d_10deg_alpha={alpha}.npy",
        ):
            path = ROOT / "GRMHD" / "data_out" / stem.format(alpha=alpha)
            array = np.load(path, mmap_mode="r")
            products[path.name] = {
                **array_summary(array),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        report["grmhd"][str(alpha)] = products

    condition_path = ROOT / "EHT-Imaging" / "data_out" / "condition.npy"
    condition = np.load(condition_path)
    report["eht"]["condition"] = {
        **array_summary(condition),
        "selected_pixels": int(condition.sum()),
        "sha256": sha256(condition_path),
    }

    for day in DAYS:
        day_report = {}
        summary_path = ROOT / "EHT-Imaging" / "data_out" / f"{day}.npz"
        with np.load(summary_path) as data:
            theta_mean, theta_std, p_mean, i_mean = (
                data["arr_0"],
                data["arr_1"],
                data["arr_2"],
                data["arr_3"],
            )
            day_report["summary"] = {
                "theta_mean": array_summary(theta_mean),
                "theta_std": array_summary(theta_std),
                "p_mean": array_summary(p_mean),
                "i_mean": array_summary(i_mean),
                "sha256": sha256(summary_path),
            }

        all_path = ROOT / "EHT-Imaging" / "data_out" / f"{day}_all_data.npy"
        all_data = np.load(all_path, mmap_mode="r")
        # Check the saved summary against the saved topset stack.  NaNs at the
        # translated border are expected and therefore use nanmean.
        stack_mean = np.nanmean(all_data, axis=0)
        day_report["all_data"] = {
            **array_summary(all_data),
            "sha256": sha256(all_path),
            "max_abs_I_mean_difference": float(np.nanmax(np.abs(stack_mean[:, :, 0] - i_mean))),
            "max_abs_P_mean_difference": float(
                np.nanmax(np.abs(np.hypot(stack_mean[:, :, 1], stack_mean[:, :, 2]) - p_mean))
            ),
        }

        one_d_path = ROOT / "EHT-Imaging" / "data_out" / f"{day}_EVPA_1d.npy"
        one_d = np.load(one_d_path, mmap_mode="r")
        day_report["evpa_1d"] = {
            **array_summary(one_d),
            "sha256": sha256(one_d_path),
        }
        report["eht"][str(day)] = day_report

    output = ROOT / "audit" / "intermediate_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
