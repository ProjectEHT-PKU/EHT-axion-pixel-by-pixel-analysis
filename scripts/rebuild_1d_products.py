#!/usr/bin/env python3
"""Rebuild corrected 10-degree 1-D EHT and GRMHD EVPA products.

Corrections relative to the notebook:
* integrate Stokes Q/U directly (they are not fractional q/u),
* use an actual 10-degree circular window instead of three 1-degree bins,
* use a cached sparse linear operator instead of per-image splines.
Outputs go to audit/rebuilt_1d and never overwrite historical intermediates.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from eht_axion.angles import circular_difference
from eht_axion.polarization import circular_window, integrated_evpa, radial_integration_matrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit" / "rebuilt_1d"
DAYS = (95, 96, 100, 101)
ALPHAS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44)


def fit_cos_evpa(evpa: np.ndarray) -> np.ndarray:
    """Fit amplitude, phase and baseline to 17-phase EVPA arrays."""
    baseline = evpa[:, :1]
    unwrapped = baseline + circular_difference(evpa[:, 1:], baseline)
    # After all phases are unwrapped around the no-axion baseline, preserve the
    # historical cosine-fit definition by taking their ordinary mean.
    average = np.mean(unwrapped, axis=1)
    average = (average + np.pi / 2) % np.pi - np.pi / 2
    delta = circular_difference(unwrapped, average[:, None])
    coefficient = np.sum(delta * np.exp(-2j * np.pi * np.arange(16) / 16)[None, :, None], axis=1)
    result = np.empty((evpa.shape[0], evpa.shape[2], 3))
    result[..., 0] = 2 * np.abs(coefficient) / 16
    result[..., 1] = np.angle(coefficient)
    result[..., 2] = average
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    operator = csr_matrix(circular_window(radial_integration_matrix(), 10))
    corrected_topset_dir = ROOT / "audit" / "rebuilt_2d"
    eht_source_dir = (
        corrected_topset_dir
        if (corrected_topset_dir / "topset_indices.npy").exists()
        else ROOT / "EHT-Imaging" / "data_out"
    )
    report = {
        "operator_nonzero": int(operator.nnz),
        "eht_source": str(eht_source_dir.relative_to(ROOT)),
        "days": {},
        "alphas": {},
    }

    for day in DAYS:
        begin = time.perf_counter()
        data = np.load(eht_source_dir / f"{day}_all_data.npy", mmap_mode="r")
        corrected = integrated_evpa(
            np.nan_to_num(data[..., 1]),
            np.nan_to_num(data[..., 2]),
            window_degrees=1,
            operator=operator,
        )
        corrected = np.roll(corrected, -90, axis=-1)
        path = OUT / f"{day}_EVPA_1d.npy"
        np.save(path, corrected)
        old = np.load(ROOT / "EHT-Imaging" / "data_out" / f"{day}_EVPA_1d.npy", mmap_mode="r")
        difference = np.rad2deg(circular_difference(corrected, old))
        report["days"][str(day)] = {
            "shape": list(corrected.shape),
            "elapsed_seconds": time.perf_counter() - begin,
            "rms_change_degrees": float(np.sqrt(np.mean(difference**2))),
            "median_absolute_change_degrees": float(np.median(np.abs(difference))),
        }
        print(f"rebuilt EHT day {day}: {report['days'][str(day)]}", flush=True)

    for alpha in ALPHAS:
        begin = time.perf_counter()
        data = np.load(
            ROOT / "GRMHD" / "data_out" / f"grmhdDataAfterImageS_alpha={alpha}_res=100.npy",
            mmap_mode="r",
        )
        shape = data.shape[:2]
        corrected = integrated_evpa(
            data[..., 1], data[..., 2], window_degrees=1, operator=operator
        ).reshape(shape + (360,))
        corrected = np.roll(corrected, -90, axis=-1)
        np.save(OUT / f"EVPA_1d_alpha={alpha}.npy", corrected)
        fit = fit_cos_evpa(corrected)
        np.save(OUT / f"EVPA_cos_fit_result_1d_alpha={alpha}.npy", fit)
        old = np.load(
            ROOT / "GRMHD" / "data_out" / f"EVPA_1d_10deg_alpha={alpha}_res=100.npy",
            mmap_mode="r",
        )
        difference = np.rad2deg(circular_difference(corrected, old))
        report["alphas"][str(alpha)] = {
            "shape": list(corrected.shape),
            "elapsed_seconds": time.perf_counter() - begin,
            "rms_change_degrees": float(np.sqrt(np.mean(difference**2))),
            "median_absolute_change_degrees": float(np.median(np.abs(difference))),
        }
        print(f"rebuilt GRMHD alpha {alpha}: {report['alphas'][str(alpha)]}", flush=True)

    report["elapsed_seconds"] = time.perf_counter() - started
    (OUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
