"""Publication-style plots from explicit, validated numerical inputs only."""

from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODES = ("std", "cor", "cor_grmhd")
COLORS = {"std": "#2679b5", "cor": "#e58f24", "cor_grmhd": "#8458a5"}
LABELS = {
    "std": "Diagonal observation noise",
    "cor": "Correlated observation noise",
    "cor_grmhd": "Correlated + GRMHD/model noise",
}


def plot_results(input_dir, output_dir):
    source, output = Path(input_dir), Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    results, reference, rows = {}, None, []
    metadata = []
    for dim in ("1d", "2d"):
        for mode in MODES:
            path = source / f"{mode}_{dim}.npz"
            with np.load(path, allow_pickle=False) as z:
                alpha, cut = z["alphas"].copy(), z["cuts"].copy()
                identity = {k: str(z[k].item()) for k in ("observation", "prior", "variant")}
            if alpha.ndim != 1 or cut.shape != alpha.shape or not np.isfinite(cut).all():
                raise ValueError(f"Invalid curve {path}")
            if reference is None:
                reference = alpha
            if not np.array_equal(reference, alpha):
                raise ValueError(f"Alpha grid mismatch: {path}")
            if metadata and identity != metadata[0]:
                raise ValueError("Mixed observations, priors or algorithm variants")
            metadata.append(identity)
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
    fig.savefig(output / "limits_by_dimension.png", dpi=200)
    plt.close(fig)
    with (output / "limits.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "figure.json").write_text(
        json.dumps(
            dict(
                metadata[0], curves=6, points=len(reference), input_directory=str(source.resolve())
            ),
            indent=2,
        )
        + "\n"
    )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    plot_results(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
