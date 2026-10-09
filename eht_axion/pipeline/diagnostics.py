"""Visual checks of intermediate physical products (no fitting or inference)."""

from pathlib import Path
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .eht import DAYS
from .layout import RunLayout, date_label
from .common import verify_stage, require_new_output_file, manifest_path, require_declared_outputs
from .inference import _stage_files, _check_lineage


def plot_diagnostics(run_root, output_dir, *, reuse_imported=False):
    layout, output = RunLayout(run_root), Path(output_dir)
    for name in ("eht_selection.png", "grmhd_background.png", "diagnostics.json"):
        require_new_output_file(output / name)
    for stage in (layout.eht, layout.background):
        verify_stage(stage, recursive=True, allow_imported=reuse_imported)
    require_declared_outputs(layout.eht, [layout.eht / "mask.npy"]
                             + [layout.eht / date_label(day) / "summary.npz" for day in DAYS])
    require_declared_outputs(layout.background,
                             [layout.background / "summary.json", layout.background / "evpa-difference.npz"])
    _check_lineage(
        verify_stage(layout.background),
        _stage_files(layout.background_processed) + [manifest_path(layout.eht), layout.eht / "mask.npy"],
    )
    output.mkdir(parents=True, exist_ok=True)
    intensity = []
    for day in DAYS:
        with np.load(layout.eht / date_label(day) / "summary.npz", allow_pickle=False) as z:
            intensity.append(z["intensity_mean"].copy())
    mask = np.load(layout.eht / "mask.npy", allow_pickle=False)
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.4), layout="constrained")
    for ax, image, day in zip(axes, intensity, DAYS):
        normalized = image / np.nanmax(image)
        ax.imshow(normalized.T, origin="lower", cmap="magma", vmin=0, vmax=1)
        ax.contour(mask.T, levels=[0.5], colors=["#7ce6d3"], linewidths=0.7)
        ax.set_title(f"{date_label(day)} · selected region", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"EHT topset mean Stokes I · common mask: {int(mask.sum())} pixels", fontsize=14)
    fig.savefig(output / "eht_selection.png", dpi=180)
    plt.close(fig)
    background = json.loads((layout.background / "summary.json").read_text())
    with np.load(layout.background / "evpa-difference.npz", allow_pickle=False) as z:
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
    (output / "diagnostics.json").write_text(json.dumps({
        "run_root": str(layout.root.resolve()), "source_validation": "recursive stage manifests",
        "reuse_imported": reuse_imported,
    }, indent=2) + "\n")
