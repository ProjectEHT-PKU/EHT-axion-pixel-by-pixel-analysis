#!/usr/bin/env python3
"""Create the final six-curve result and legacy/corrected comparison figures."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "audit" / "results" / "legacy_expected"
CORRECTED = ROOT / "audit" / "results" / "rebuilt_corrected_expected"
OUT = ROOT / "audit" / "figures"
MODES = ("std", "cor", "cor_grmhd")
DIMENSIONS = ("1d", "2d")
LABELS = {
    "std": "Diagonal observation uncertainty",
    "cor": "Observation covariance",
    "cor_grmhd": "Observation covariance + astrophysical/model uncertainty",
}
COLORS = {"std": "#0072B2", "cor": "#E69F00", "cor_grmhd": "#CC79A7"}
LINESTYLES = {"1d": "--", "2d": "-"}
MARKERS = {"1d": "s", "2d": "o"}


def load(directory: Path, mode: str, dimension: str):
    with np.load(directory / f"{mode}_{dimension}.npz") as data:
        return data["alphas"], data["cuts"]


def style_axes(ax):
    ax.set_xlabel(r"$\alpha = M_{\rm BH} m_a$")
    ax.set_ylabel(r"$\log_{10}(c_{95})$")
    ax.grid(True, alpha=0.22)
    ax.spines[["top", "right"]].set_visible(False)


def final_figure():
    fig, ax = plt.subplots(figsize=(8.2, 5.6), constrained_layout=True)
    for mode in MODES:
        for dimension in DIMENSIONS:
            alpha, cuts = load(CORRECTED, mode, dimension)
            ax.plot(
                alpha,
                cuts,
                color=COLORS[mode],
                linestyle=LINESTYLES[dimension],
                marker=MARKERS[dimension],
                markersize=4,
                linewidth=2,
                label=f"{LABELS[mode]} — {dimension.upper()}",
            )
    style_axes(ax)
    ax.set_title("M87* axion coupling limits: three uncertainty models")
    ax.legend(fontsize=8.2, ncol=1)
    fig.savefig(OUT / "final_limits_corrected.pdf")
    fig.savefig(OUT / "final_limits_corrected.png", dpi=220)
    plt.close(fig)


def comparison_figure():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), sharey=True, constrained_layout=True)
    for ax, dimension in zip(axes, DIMENSIONS):
        for mode in MODES:
            alpha, old = load(LEGACY, mode, dimension)
            _, new = load(CORRECTED, mode, dimension)
            ax.plot(
                alpha,
                old,
                color=COLORS[mode],
                linestyle=":",
                linewidth=1.6,
                label=f"{LABELS[mode]} — legacy",
            )
            ax.plot(
                alpha,
                new,
                color=COLORS[mode],
                linestyle="-",
                marker="o",
                markersize=3,
                linewidth=2,
                label=f"{LABELS[mode]} — corrected",
            )
        style_axes(ax)
        ax.set_title(dimension.upper())
    axes[1].legend(fontsize=7.2, loc="best")
    fig.suptitle("Impact of verified implementation fixes")
    fig.savefig(OUT / "legacy_vs_corrected.pdf")
    fig.savefig(OUT / "legacy_vs_corrected.png", dpi=220)
    plt.close(fig)


def write_tables():
    rows = []
    for mode in MODES:
        for dimension in DIMENSIONS:
            alpha, legacy = load(LEGACY, mode, dimension)
            _, corrected = load(CORRECTED, mode, dimension)
            for a, old, new in zip(alpha, legacy, corrected):
                rows.append(
                    {
                        "alpha": float(a),
                        "mode": mode,
                        "dimension": dimension,
                        "legacy_log10_c95": float(old),
                        "corrected_log10_c95": float(new),
                        "change": float(new - old),
                    }
                )
    with (OUT / "final_limits.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "maximum_absolute_change": max(abs(row["change"]) for row in rows),
        "curves": 6,
        "alpha_points_per_curve": 8,
        "rows": rows,
    }
    (OUT / "final_limits.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    final_figure()
    comparison_figure()
    write_tables()
    print(f"Wrote final figures and tables to {OUT}")


if __name__ == "__main__":
    main()
