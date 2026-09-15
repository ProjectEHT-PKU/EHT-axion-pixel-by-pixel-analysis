#!/usr/bin/env python3
"""Rebuild four-day EHT products from the coordinate-aligned topset ranking.

Outputs are isolated under audit/rebuilt_2d and never overwrite historical data.
"""

from __future__ import annotations

import gc
import json
import pickle
import time
from pathlib import Path

import numpy as np

from eht_axion.angles import circular_difference, circular_mean
from utils.libanaly import fit_gaussian2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit" / "rebuilt_2d"
DAYS = (95, 96, 100, 101)
FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
OFFSETS = {95: (1, 0), 96: (0, 1), 100: (-1, 2), 101: (0, 0)}


def sorted_order(metadata: list[dict]) -> np.ndarray:
    return np.asarray(
        sorted(range(len(metadata)), key=lambda i: tuple(metadata[i][x] for x in FIELDS))
    )


def shift_spatial(array: np.ndarray, offset: tuple[int, int]) -> np.ndarray:
    xoff, yoff = offset
    pad = 9
    padded = np.full(array.shape[:-2] + (80, 80), np.nan, dtype=array.dtype)
    padded[..., pad : pad + 64, pad : pad + 64] = array
    return padded[..., pad + xoff : pad + xoff + 64, pad + yoff : pad + yoff + 64]


def fit_theta(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = circular_mean(theta, axis=0)
    centered = center + circular_difference(theta, center)
    fit = np.apply_along_axis(
        fit_gaussian2,
        0,
        centered,
        plot_failed=False,
        angle=True,
        lower_bound=-np.pi / 2,
        upper_bound=np.pi / 2,
        bins=100,
        threshold=0.2,
    )
    mean = (fit[0] + np.pi / 2) % np.pi - np.pi / 2
    return mean, np.abs(fit[1])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ranking = np.load(ROOT / "audit" / "topset_ranking_comparison.npz")
    selection = ranking["corrected_order"][:1572]
    np.save(OUT / "topset_indices.npy", selection)
    condition = np.ones((64, 64), dtype=bool)
    report = {"selection_size": int(selection.size), "days": {}}
    started = time.perf_counter()

    for day in DAYS:
        begin = time.perf_counter()
        with (ROOT / "EHT-Imaging" / "data_in" / f"230_{day}.pkl").open("rb") as stream:
            raw = pickle.load(stream)
        order = sorted_order(raw["meta_data"])
        chosen = order[selection]
        all_data = np.zeros((selection.size, 64, 64, 4), dtype=np.float64)
        all_data[..., 0] = np.rot90(
            np.transpose(raw["I_data"][::-1][:, :, chosen], (2, 0, 1)), k=-1, axes=(1, 2)
        )
        all_data[..., 1] = -np.rot90(
            np.transpose(raw["Q_data"][::-1][:, :, chosen], (2, 0, 1)), k=-1, axes=(1, 2)
        )
        all_data[..., 2] = -np.rot90(
            np.transpose(raw["U_data"][::-1][:, :, chosen], (2, 0, 1)), k=-1, axes=(1, 2)
        )
        theta = 0.5 * np.arctan2(all_data[..., 2], all_data[..., 1])
        mean, std = fit_theta(theta)
        all_data = shift_spatial(all_data.transpose(0, 3, 1, 2), OFFSETS[day]).transpose(0, 2, 3, 1)
        mean = shift_spatial(mean, OFFSETS[day])
        std = shift_spatial(std, OFFSETS[day])
        all_mean = np.mean(all_data, axis=0)
        p_mean = np.hypot(all_mean[..., 1], all_mean[..., 2])
        np.save(OUT / f"{day}_all_data.npy", all_data)
        np.save(OUT / f"{day}_all_data_mean.npy", all_mean)
        np.savez(OUT / f"{day}.npz", mean, std, p_mean, all_mean[..., 0])
        condition &= all_mean[..., 0] > np.nanmax(all_mean[..., 0]) * 0.1
        condition &= np.isfinite(all_mean[..., 0]) & np.isfinite(std)
        report["days"][str(day)] = {
            "offset": list(OFFSETS[day]),
            "finite_fit_pixels": int(np.isfinite(std).sum()),
            "elapsed_seconds": time.perf_counter() - begin,
        }
        print(f"rebuilt day {day}: {report['days'][str(day)]}", flush=True)
        del raw, all_data, theta, mean, std, all_mean
        gc.collect()

    grid = np.indices((64, 64))
    radius = np.hypot(grid[0] - 31.5, grid[1] - 31.5) * 0.51073125
    condition &= radius >= 2
    condition[:, :5] = False
    condition[48:55, 10:19] = False
    np.save(OUT / "condition.npy", condition)
    old_condition = np.load(ROOT / "EHT-Imaging" / "data_out" / "condition.npy")
    report.update(
        {
            "selected_pixels": int(condition.sum()),
            "historical_selected_pixels": int(old_condition.sum()),
            "condition_overlap": int(np.sum(condition & old_condition)),
            "elapsed_seconds": time.perf_counter() - started,
        }
    )
    (OUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
