#!/usr/bin/env python3
"""Infer historical topset indices from saved Stokes stacks for all four days."""

from __future__ import annotations

import gc
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DAYS = (95, 96, 100, 101)
FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
SIGNATURE_COORDINATES = ((31, 31), (30, 33), (36, 28))
OFFSETS = {95: (1, 0), 96: (0, 1), 100: (-1, 2), 101: (0, 0)}


def sorted_order(metadata: list[dict]) -> np.ndarray:
    return np.asarray(
        sorted(
            range(len(metadata)),
            key=lambda index: tuple(metadata[index][field] for field in FIELDS),
        )
    )


def infer_indices(day: int) -> tuple[np.ndarray, float]:
    with (ROOT / "EHT-Imaging" / "data_in" / f"230_{day}.pkl").open("rb") as stream:
        raw = pickle.load(stream)
    order = sorted_order(raw["meta_data"])
    # Saved[x,y] = rotated[x+xoff,y+yoff] after the padded translation, and
    # np.rot90(k=-1): rotated[x,y] = input[n-1-y,x].
    xoff, yoff = OFFSETS[day]
    candidates = np.stack(
        [raw["I_data"][::-1][63 - (y + yoff), x + xoff, order] for x, y in SIGNATURE_COORDINATES],
        axis=1,
    )
    saved = np.load(ROOT / "EHT-Imaging" / "data_out" / f"{day}_all_data.npy", mmap_mode="r")
    targets = np.stack([saved[:, x, y, 0] for x, y in SIGNATURE_COORDINATES], axis=1)
    lookup = {tuple(row): index for index, row in enumerate(candidates)}
    inferred = np.array([lookup.get(tuple(row), -1) for row in targets])
    for target_index in np.flatnonzero(inferred < 0):
        inferred[target_index] = np.argmin(
            np.max(np.abs(candidates - targets[target_index]), axis=1)
        )
    error = float(np.max(np.abs(candidates[inferred] - targets)))
    del raw, saved, candidates, targets
    gc.collect()
    return inferred, error


def main() -> None:
    ranking = np.load(ROOT / "audit" / "topset_ranking_comparison.npz")
    legacy = ranking["legacy_order"][:1572]
    corrected = ranking["corrected_order"][:1572]
    report = {"days": {}}
    inferred_by_day = {}
    for day in DAYS:
        inferred, error = infer_indices(day)
        inferred_by_day[day] = inferred
        report["days"][str(day)] = {
            "unique_indices": int(np.unique(inferred).size),
            "legacy_overlap": int(np.intersect1d(inferred, legacy).size),
            "corrected_overlap": int(np.intersect1d(inferred, corrected).size),
            "maximum_signature_error": error,
        }
        print(f"validated day {day}: {report['days'][str(day)]}", flush=True)
    first = inferred_by_day[DAYS[0]]
    report["all_days_same_index_set"] = all(
        np.array_equal(np.sort(first), np.sort(inferred_by_day[day])) for day in DAYS[1:]
    )
    report["all_days_reproduce_legacy_selection"] = all(
        values["legacy_overlap"] == 1572 and values["maximum_signature_error"] == 0
        for values in report["days"].values()
    )
    output = ROOT / "audit" / "historical_selection_validation.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_days_reproduce_legacy_selection"]:
        raise SystemExit("historical selection did not reproduce for every day")


if __name__ == "__main__":
    main()
