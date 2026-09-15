#!/usr/bin/env python3
"""Recompute the EHT 4800-image ranking before and after coordinate alignment.

The historical notebook translates the fitted mean image but computes ranking
residuals from untranslated per-image EVPAs.  This script reproduces that score
and compares it with aligned residuals.  Large temporary arrays are written to
audit/ rather than retained in memory across all four days.
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
DAYS = (95, 96, 100, 101)
FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
OFFSETS = {95: (2, 0), 96: (1, 1), 100: (-4, 4), 101: (-1, 1)}
TOPSET_SIZE = 1572


def sorted_order(metadata: list[dict]) -> np.ndarray:
    return np.asarray(
        sorted(range(len(metadata)), key=lambda i: tuple(metadata[i][x] for x in FIELDS))
    )


def shift_2d(array: np.ndarray, offset: tuple[int, int]) -> np.ndarray:
    """Match the notebook's padded slicing translation for trailing 2 axes."""
    xoff, yoff = offset
    pad = 9
    padded_shape = array.shape[:-2] + (array.shape[-2] + 16, array.shape[-1] + 16)
    padded = np.full(padded_shape, np.nan, dtype=array.dtype)
    padded[..., pad : pad + 64, pad : pad + 64] = array
    return padded[..., pad + xoff : pad + xoff + 64, pad + yoff : pad + yoff + 64]


def fitted_mean_std(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # The von Mises MLE location used in the notebook equals the circular mean
    # for fixed unit scale.  This vectorized form is much faster and is then
    # followed by the same histogram Gaussian fit.
    center = circular_mean(theta, axis=0)
    centered = center + circular_difference(theta, center)
    result = np.apply_along_axis(
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
    mean = (result[0] + np.pi / 2) % np.pi - np.pi / 2
    return mean, np.abs(result[1])


def main() -> None:
    work = ROOT / "audit" / "topset_work"
    work.mkdir(parents=True, exist_ok=True)
    condition = np.ones((64, 64), dtype=bool)
    summaries = {}
    started = time.perf_counter()

    for day in DAYS:
        day_started = time.perf_counter()
        path = ROOT / "EHT-Imaging" / "data_in" / f"230_{day}.pkl"
        with path.open("rb") as stream:
            raw = pickle.load(stream)
        order = sorted_order(raw["meta_data"])
        np.save(work / f"order_{day}.npy", order)
        stokes = np.empty((4800, 64, 64, 3), dtype=np.float64)
        stokes[..., 0] = np.transpose(raw["I_data"][::-1][:, :, order], (2, 0, 1))
        stokes[..., 1] = -np.rot90(
            np.transpose(raw["Q_data"][::-1][:, :, order], (2, 0, 1)), k=-1, axes=(1, 2)
        )
        stokes[..., 2] = -np.rot90(
            np.transpose(raw["U_data"][::-1][:, :, order], (2, 0, 1)), k=-1, axes=(1, 2)
        )
        # I follows the same spatial rotation but no Stokes sign change.
        stokes[..., 0] = np.rot90(stokes[..., 0], k=-1, axes=(1, 2))
        theta = 0.5 * np.arctan2(stokes[..., 2], stokes[..., 1])
        mean, std = fitted_mean_std(theta)
        shifted_mean = shift_2d(mean, OFFSETS[day])
        shifted_std = shift_2d(std, OFFSETS[day])
        shifted_i_mean = shift_2d(np.mean(stokes[..., 0], axis=0), OFFSETS[day])
        condition &= shifted_i_mean > np.nanmax(shifted_i_mean) * 0.1
        condition &= np.isfinite(shifted_i_mean) & np.isfinite(shifted_std)
        np.save(work / f"theta_{day}.npy", theta)
        np.savez(work / f"fit_{day}.npz", mean=shifted_mean, std=shifted_std)
        summaries[str(day)] = {
            "offset": list(OFFSETS[day]),
            "fit_finite_pixels": int(np.isfinite(shifted_std).sum()),
            "elapsed_seconds": time.perf_counter() - day_started,
        }
        del raw, stokes, theta, mean, std
        gc.collect()
        print(f"prepared day {day}: {summaries[str(day)]}", flush=True)

    grid = np.indices((64, 64))
    radius = np.hypot(grid[0] - 31.5, grid[1] - 31.5) * 0.51073125
    condition &= radius >= 2
    condition[:, :5] = False
    condition[48:55, 10:19] = False
    selected_pixels = np.flatnonzero(condition.ravel())
    old_score = np.zeros(4800)
    fixed_score = np.zeros(4800)

    for day in DAYS:
        theta = np.load(work / f"theta_{day}.npy", mmap_mode="r")
        with np.load(work / f"fit_{day}.npz") as fit:
            mean, std = fit["mean"].ravel(), fit["std"].ravel()
        x, y = np.unravel_index(selected_pixels, (64, 64))
        xoff, yoff = OFFSETS[day]
        legacy_theta = theta[:, x, y]
        aligned_theta = theta[:, x + xoff, y + yoff]
        old_residual = circular_difference(legacy_theta, mean[selected_pixels])
        fixed_residual = circular_difference(aligned_theta, mean[selected_pixels])
        old_score += np.sum((old_residual / std[selected_pixels]) ** 2, axis=1)
        fixed_score += np.sum((fixed_residual / std[selected_pixels]) ** 2, axis=1)
        print(f"scored day {day}", flush=True)

    old_order = np.argsort(old_score)
    fixed_order = np.argsort(fixed_score)
    old_top = old_order[:TOPSET_SIZE]
    fixed_top = fixed_order[:TOPSET_SIZE]
    overlap = np.intersect1d(old_top, fixed_top)
    rank_old = np.empty(4800, dtype=int)
    rank_fixed = np.empty(4800, dtype=int)
    rank_old[old_order] = np.arange(4800)
    rank_fixed[fixed_order] = np.arange(4800)
    report = {
        "selected_pixels": int(selected_pixels.size),
        "topset_size": TOPSET_SIZE,
        "topset_overlap": int(overlap.size),
        "topset_changed": int(TOPSET_SIZE - overlap.size),
        "spearman_rank_correlation": float(np.corrcoef(rank_old, rank_fixed)[0, 1]),
        "maximum_absolute_rank_change": int(np.max(np.abs(rank_old - rank_fixed))),
        "elapsed_seconds": time.perf_counter() - started,
        "days": summaries,
    }
    np.savez(
        ROOT / "audit" / "topset_ranking_comparison.npz",
        legacy_score=old_score,
        corrected_score=fixed_score,
        legacy_order=old_order,
        corrected_order=fixed_order,
        selected_pixels=selected_pixels,
    )
    (ROOT / "audit" / "topset_ranking_comparison.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
