"""Visual checks of intermediate physical products (no fitting or inference)."""

from pathlib import Path
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_diagnostics(run_root, output_dir):
    run, output = Path(run_root), Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with np.load(run / "04_eht" / "summary_2d.npz") as z:
        intensity = z["intensity_mean"].copy()
    mask = np.load(run / "04_eht" / "mask.npy")
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.4), layout="constrained")
    for ax, image, day in zip(axes, intensity, [5, 6, 10, 11]):
        normalized = image / np.nanmax(image)
        ax.imshow(normalized.T, origin="lower", cmap="magma", vmin=0, vmax=1)
        ax.contour(mask.T, levels=[0.5], colors=["#7ce6d3"], linewidths=0.7)
        ax.set_title(f"April {day} · selected region", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"EHT topset mean Stokes I · common mask: {int(mask.sum())} pixels", fontsize=14)
    fig.savefig(output / "eht_selection.png", dpi=180)
    plt.close(fig)
    background = json.loads((run / "05_background" / "background.json").read_text())
    with np.load(run / "05_background" / "background_diagnostics.npz") as z:
        edges, counts = np.rad2deg(z["hist_edges"]), z["hist_counts"]
    density = counts / (counts.sum() * np.diff(edges))
    mu, sigma = np.rad2deg(background["mean"]), np.rad2deg(background["sigma"])
    fig, ax = plt.subplots(figsize=(6.5, 4), layout="constrained")
    ax.stairs(density, edges, color="#2679b5", fill=True, alpha=0.4, label="GRMHD differences")
    if sigma > 0:
        x = np.linspace(mu - 5 * sigma, mu + 5 * sigma, 600)
        ax.plot(
            x,
            np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi)),
            color="#bf5a46",
            lw=2,
            label=f"Gaussian MLE σ={sigma:.3f}°",
        )
    ax.set(
        xlabel="No-axion EVPA difference (degree)",
        ylabel="Probability density",
        title=f"Intrinsic background · lag {background['lag_hours']:.2f} h",
    )
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "grmhd_background.png", dpi=180)
    plt.close(fig)


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    plot_diagnostics(args.run_root, args.output_dir)


if __name__ == "__main__":
    main()
