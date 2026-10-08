#!/usr/bin/env python3
"""Compare original saved six curves with the script pipeline (no recomputation)."""

from pathlib import Path
import argparse
import csv
import hashlib
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from eht_axion.pipeline.layout import load_curves

MODES = ("std", "cor", "cor_grmhd")
COLORS = {"std": "#2679b5", "cor": "#e58f24", "cor_grmhd": "#8458a5"}


def load_curve(path):
    with np.load(path, allow_pickle=False) as data:
        alpha, cuts = data["alphas"].copy(), data["cuts"].copy()
    if (
        alpha.ndim != 1
        or cuts.shape != alpha.shape
        or len(alpha) < 2
        or not np.isfinite(alpha).all()
        or not np.isfinite(cuts).all()
        or np.any(np.diff(alpha) <= 0)
    ):
        raise ValueError(f"Invalid curve: {path}")
    return alpha, cuts


def compare(original_dir, pipeline_dir, output_dir):
    original, pipeline, output = map(Path, (original_dir, pipeline_dir, output_dir))
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Use a new empty output directory: {output}")
    curves, rows, summary, inputs = {}, [], {}, []
    for dimension in ("1d", "2d"):
        for mode in MODES:
            old_path = original / f"{mode}{'_1d' if dimension == '1d' else ''}.npz"
            new_path = pipeline / f"{mode}_{dimension}.npz"
            a, old = load_curve(old_path)
            if (pipeline / f"limits-{dimension}.npz").is_file():
                b, all_cuts, identity = load_curves(pipeline, dimension)
                new = all_cuts[MODES.index(mode)]
                new_path = pipeline / f"limits-{dimension}.npz"
            else:
                b, new = load_curve(new_path)
                with np.load(new_path, allow_pickle=False) as data:
                    identity = {"observation": str(data["observation"].item())}
            if a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-14):
                raise ValueError(f"Alpha grids differ: {old_path}, {new_path}")
            if identity["observation"] not in ("null", "expected"):
                raise ValueError("Original null curves must not be compared with observed data")
            difference, ratio = new - old, 10 ** (new - old)
            curves[dimension, mode] = (a, old, new)
            for alpha, previous, current, delta, factor in zip(a, old, new, difference, ratio):
                rows.append(
                    dict(
                        dimension=dimension,
                        mode=mode,
                        alpha=float(alpha),
                        original_log10_c95=float(previous),
                        pipeline_log10_c95=float(current),
                        delta_dex=float(delta),
                        pipeline_over_original=float(factor),
                    )
                )
            summary[f"{mode}_{dimension}"] = {
                "delta_dex_range": [float(difference.min()), float(difference.max())],
                "pipeline_over_original_range": [float(ratio.min()), float(ratio.max())],
                "max_abs_delta_dex": float(np.max(abs(difference))),
            }
            for path in (old_path, new_path):
                inputs.append(
                    dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                )
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(
        2, 2, figsize=(11.2, 7.8), sharex=True, layout="constrained", gridspec_kw={"height_ratios": [2.2, 1]}
    )
    for column, dimension in enumerate(("1d", "2d")):
        for mode in MODES:
            alpha, old, new = curves[dimension, mode]
            axes[0, column].plot(
                alpha, old, color=COLORS[mode], ls="--", marker="s", ms=4, markerfacecolor="white", lw=1.6
            )
            axes[0, column].plot(alpha, new, color=COLORS[mode], marker="o", ms=4, lw=2.1)
            axes[1, column].plot(alpha, new - old, color=COLORS[mode], marker="o", ms=3.5, lw=1.8)
        axes[0, column].set_title(
            "1D angular apertures" if dimension == "1d" else "2D pixels",
            fontsize=13,
            fontweight="semibold",
            loc="left",
        )
        axes[1, column].axhline(0, color="#7b8790", lw=0.8)
        axes[1, column].set_xlabel(r"$\alpha$")
        for ax in axes[:, column]:
            ax.grid(color="#e4e8ed", lw=0.7)
            ax.set_axisbelow(True)
    axes[0, 0].set_ylabel(r"$\log_{10} c_{95}$")
    axes[1, 0].set_ylabel(r"New $-$ original (dex)")
    handles = [Line2D([0], [0], color=COLORS[m], lw=2, label=m) for m in MODES]
    handles += [
        Line2D(
            [0],
            [0],
            color="#454b52",
            ls="--",
            marker="s",
            markerfacecolor="white",
            label="Original saved results",
        ),
        Line2D([0], [0], color="#454b52", marker="o", label="Script research profile"),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncols=3, frameon=False)
    fig.suptitle("M87* | Original vs script pipeline", fontsize=17, fontweight="semibold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(output / f"original_vs_pipeline.{ext}", dpi=200)
    plt.close(fig)

    # Six small panels preserve each mode's identity rather than merging twelve lines.
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.7), sharex=True, layout="constrained")
    for row, dim in enumerate(("1d", "2d")):
        for col, mode in enumerate(MODES):
            alpha, old, new = curves[dim, mode]
            ax = axes[row, col]
            ax.plot(
                alpha,
                old,
                color="#687580",
                ls="--",
                marker="s",
                ms=4,
                markerfacecolor="white",
                label="Original",
            )
            ax.plot(alpha, new, color=COLORS[mode], marker="o", ms=4, label="Script pipeline")
            ax.set_title(f"{dim.upper()} · {mode}", loc="left", fontweight="semibold")
            ax.grid(color="#e4e8ed", lw=0.7)
            if row == 1:
                ax.set_xlabel(r"$\alpha$")
            if col == 0:
                ax.set_ylabel(r"$\log_{10} c_{95}$")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Null-data comparison · method changes are not a one-bug effect", fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(output / f"six_curves_comparison.{ext}", dpi=200)
    plt.close(fig)
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "comparison": "original saved null results vs script research-profile null results",
        "warning": (
            "Includes observable, estimator, selection, background and formula changes; "
            "not isolated bug effect."
        ),
        "sign": "negative delta / ratio < 1 means a smaller numerical upper limit",
        "curves": summary,
        "inputs": inputs,
    }
    (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-dir", type=Path, required=True)
    parser.add_argument("--pipeline-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.original_dir, args.pipeline_dir, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
