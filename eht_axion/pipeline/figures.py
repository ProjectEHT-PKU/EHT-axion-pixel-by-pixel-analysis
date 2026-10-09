"""Publication-style plots from explicit, validated numerical inputs only."""

from pathlib import Path
import csv
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .layout import load_curves
from .common import (
    verify_stage, reject_directory_symlinks, require_new_output_file, require_declared_outputs,
)

MODES = ("std", "cor", "cor_grmhd")
COLORS = {"std": "#2679b5", "cor": "#e58f24", "cor_grmhd": "#8458a5"}
LABELS = {
    "std": "Diagonal observation noise",
    "cor": "Correlated observation noise",
    "cor_grmhd": "Correlated + GRMHD/model noise",
}


def load_limits(input_dir, dimension, *, allow_legacy=False):
    """Read one dimension, preserving the fixed mode order and scalar identity.

    Legacy curves are accepted only for explicit direct plotting, not as inputs
    to a canonical run. A partial canonical result never falls back to old files.
    """
    source = Path(input_dir)
    path = source / f"limits-{dimension}.npz"
    canonical = any((source / f"limits-{dim}.npz").exists() for dim in ("1d", "2d"))
    if (canonical or not allow_legacy) and not path.is_file():
        raise FileNotFoundError(f"Missing canonical limits: {path}")
    alpha, cuts, identity = load_curves(source, dimension)
    if (
        alpha.ndim != 1
        or not len(alpha)
        or not np.isfinite(alpha).all()
        or np.any(np.diff(alpha) <= 0)
        or cuts.shape != (len(MODES), len(alpha))
        or not np.isfinite(cuts).all()
    ):
        raise ValueError(f"Invalid curves/grid: {path}")
    return alpha, cuts, identity


def plot_results(input_dir, output_dir, *, reuse_imported=False):
    source, output = Path(input_dir), Path(output_dir)
    reject_directory_symlinks(output)
    for name in ("limits.png", "limits.pdf", "limits.svg", "limits-by-dimension.png",
                 "limits.csv", "figure.json"):
        require_new_output_file(output / name)
    results, reference, rows = {}, None, []
    metadata = []
    for dim in ("1d", "2d"):
        alpha, cuts, identity = load_limits(source, dim, allow_legacy=True)
        if reference is None:
            reference = alpha
        if not np.array_equal(reference, alpha):
            raise ValueError(f"Alpha grid mismatch: {source / f'limits-{dim}.npz'}")
        if metadata and identity != metadata[0]:
            raise ValueError("Mixed observations, priors or algorithm variants")
        metadata.append(identity)
        for mode, cut in zip(MODES, cuts):
            results[dim, mode] = cut
            rows.extend(
                {
                    "alpha": float(a),
                    "dimension": dim,
                    "mode": mode,
                    "log10_c95": float(c),
                    "c95": float(10**c),
                }
                for a, c in zip(alpha, cut)
            )
    canonical = any((source / f"limits-{dim}.npz").is_file() for dim in ("1d", "2d"))
    has_manifest = any((source / name).is_file() for name in ("manifest.json", "stage.json"))
    if canonical or has_manifest:
        record = verify_stage(source, recursive=True, allow_imported=reuse_imported)
        consumed = ([source / f"limits-{dim}.npz" for dim in ("1d", "2d")] if canonical else
                    [source / f"{mode}_{dim}.npz" for mode in MODES for dim in ("1d", "2d")])
        require_declared_outputs(source, consumed)
        if record["stage"] == "inference:six-curves" and canonical:
            from .inference import validate_run_lineage

            validate_run_lineage(source.parent, reference, allow_imported=reuse_imported)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "semibold",
            "savefig.facecolor": "white",
        }
    )
    fig, ax = plt.subplots(figsize=(8.4, 5.8), layout="constrained")
    for mode in MODES:
        for dim in ("1d", "2d"):
            ax.plot(
                reference,
                results[dim, mode],
                color=COLORS[mode],
                lw=2.2 if dim == "2d" else 1.7,
                ls="-" if dim == "2d" else "--",
                marker="o" if dim == "2d" else "s",
                ms=4.5,
                markerfacecolor=COLORS[mode] if dim == "2d" else "white",
            )
    from matplotlib.lines import Line2D

    handles = [Line2D([0], [0], color=COLORS[m], lw=2.2, label=LABELS[m]) for m in MODES]
    handles += [
        Line2D(
            [0],
            [0],
            color="#42454a",
            ls="--",
            marker="s",
            markerfacecolor="white",
            label="1D apertures",
        ),
        Line2D([0], [0], color="#42454a", ls="-", marker="o", label="2D pixels"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="best")
    ax.set(
        xlabel=r"Axion gravitational coupling $\alpha$",
        ylabel=r"$\log_{10} c_{95}$",
        xlim=(reference.min() - 0.01, reference.max() + 0.01),
    )
    observation = metadata[0]["observation"]
    fig.suptitle(
        "M87*  |  Axion–photon coupling limits",
        x=0.08,
        ha="left",
        fontsize=17,
        fontweight="semibold",
    )
    ax.set_title(
        f"{observation.title()} data · 95% posterior upper limits · {metadata[0]['prior']} prior",
        loc="left",
        fontsize=10,
        color="#66717d",
        pad=15,
    )
    ax.grid(color="#dfe4e8", alpha=0.75, lw=0.7)
    ax.set_axisbelow(True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(output / f"limits.{extension}", dpi=200)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True, layout="constrained")
    for ax, dim in zip(axes, ("1d", "2d")):
        for mode in MODES:
            ax.plot(
                reference,
                results[dim, mode],
                color=COLORS[mode],
                marker="o",
                lw=2,
                ms=4,
                label=LABELS[mode],
            )
        ax.set(
            title="1D angular apertures" if dim == "1d" else "2D selected pixels",
            xlabel=r"$\alpha$",
        )
        ax.grid(color="#e5e9ee", lw=0.7)
    axes[0].set_ylabel(r"$\log_{10} c_{95}$")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(output / "limits-by-dimension.png", dpi=200)
    plt.close(fig)
    with (output / "limits.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "figure.json").write_text(
        json.dumps(
            dict(metadata[0], curves=6, points=len(reference), input_directory=str(source.resolve()),
                 source_validation="recursive stage manifest" if canonical or has_manifest
                 else "legacy numerical files only; scientific lineage unavailable"),
            indent=2,
        )
        + "\n"
    )
    return output
