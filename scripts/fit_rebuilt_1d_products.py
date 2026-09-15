#!/usr/bin/env python3
"""Fit corrected one-dimensional EHT and GRMHD intermediates."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from eht_axion.angles import circular_difference, circular_mean

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "audit" / "rebuilt_1d"
DAYS = (95, 96, 100, 101)
ALPHAS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44)


def fit_columns(values: np.ndarray, *, period: float | None = None) -> np.ndarray:
    """Fit every column with the unbinned Gaussian maximum likelihood estimator.

    The historical histogram+curve_fit result depends on an arbitrary bin count
    and rejects otherwise finite columns when the optimizer covariance exceeds
    a threshold.  For a Gaussian, the direct MLE is both faster and stable.
    """
    data = np.asarray(values, dtype=float)
    if period is not None:
        center = circular_mean(data, axis=0, period=period)
        data = center + circular_difference(data, center, period=period)
    return np.vstack((np.mean(data, axis=0), np.std(data, axis=0)))


def main() -> None:
    started = time.perf_counter()
    report = {"days": {}, "alphas": {}}
    for day in DAYS:
        begin = time.perf_counter()
        evpa = np.load(DATA / f"{day}_EVPA_1d.npy", mmap_mode="r")
        fit = fit_columns(evpa, period=np.pi)
        fit[0] = (fit[0] + np.pi / 2) % np.pi - np.pi / 2
        np.savez(DATA / f"{day}_1d.npz", fit[0], np.abs(fit[1]))
        report["days"][str(day)] = {
            "finite_bins": int(np.isfinite(fit).all(axis=0).sum()),
            "elapsed_seconds": time.perf_counter() - begin,
        }
        print(f"fit EHT day {day}: {report['days'][str(day)]}", flush=True)

    amplitude_parameters, phase_parameters = {}, {}
    for alpha in ALPHAS:
        begin = time.perf_counter()
        fit = np.load(DATA / f"EVPA_cos_fit_result_1d_alpha={alpha}.npy", mmap_mode="r")
        amplitude = fit_columns(np.asarray(fit[..., 0]))
        phase = fit_columns(np.asarray(fit[..., 1]), period=2 * np.pi)
        amplitude_parameters[alpha] = amplitude
        phase_parameters[alpha] = phase
        report["alphas"][str(alpha)] = {
            "finite_amplitude_bins": int(np.isfinite(amplitude).all(axis=0).sum()),
            "finite_phase_bins": int(np.isfinite(phase).all(axis=0).sum()),
            "elapsed_seconds": time.perf_counter() - begin,
        }
        print(f"fit GRMHD alpha {alpha}: {report['alphas'][str(alpha)]}", flush=True)
    np.save(DATA / "A_params_1d_dict.npy", amplitude_parameters, allow_pickle=True)
    np.save(DATA / "phi0_params_1d_dict.npy", phase_parameters, allow_pickle=True)
    report["elapsed_seconds"] = time.perf_counter() - started
    (DATA / "fit_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
