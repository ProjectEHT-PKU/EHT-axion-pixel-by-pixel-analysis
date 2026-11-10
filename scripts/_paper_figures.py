"""Figures of the paper (vector PDF): style and one function per figure.

Every function receives arrays that analysis_paper.py has already read and checked,
writes one PDF and returns the numbers it derived for the figure. Sky maps are stored
as template arrays [x, y] with x = -East and y = North; they are shown transposed with
East to the left and North up, without any further rotation.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LogNorm, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

COLORS = ("#0072B2", "#D55E00", "#009E73")
PAIRS = ("Apr 6 − 5", "Apr 11 − 10")
MODES = ("std", "cor", "cor_fixed_background")
MODE_LABELS = ("std: diagonal D-term", "cor: full D-term", "cor+bg: full D-term + background")
STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.3,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "axes.spines.top": False,
    "axes.spines.right": False,
}
FIGURES = (
    "radial_function",
    "two_photon",
    "signal_maps",
    "template_maps",
    "background",
    "observational_uncertainty",
    "data_vs_model",
    "limits",
)


def save(fig, out, name):
    """Deterministic vector PDF (no creation dates)."""
    fig.savefig(
        Path(out) / f"{name}.pdf",
        metadata={"Creator": "analysis_paper.py", "Producer": None, "CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


def grid(ax):
    ax.grid(alpha=0.18, linewidth=0.5)


def sky(ax, image, fov, **kwargs):
    """Stored [x=-East, y=North] map as vector cells, East to the left; no data rotation."""
    edges = np.linspace(-fov / 2, fov / 2, image.shape[0] + 1)
    artist = ax.pcolormesh(
        -edges,
        edges,
        image.T,
        shading="flat",
        rasterized=False,
        edgecolors="face",
        linewidth=0.35,
        antialiased=False,
        **kwargs,
    )
    ax.set(xlim=(fov / 2, -fov / 2), ylim=(-fov / 2, fov / 2), xlabel="East offset (μas)", aspect="equal")
    ax.set_xticks([-40, 0, 40])
    ax.set_yticks([-40, 0, 40])
    return artist


def colorbar(fig, ax, artist, label):
    cb = fig.colorbar(artist, ax=ax, orientation="horizontal", pad=0.18, fraction=0.07, aspect=22)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor("face")
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=8)
    return cb


def radial_function(out, profiles, window, horizon_rg, switch_alpha, sample_radii):
    """profiles: [(alpha, table spin, high-spin table?, r, R/Rmax)]; window: (radius, exponent)."""
    fig, axes = plt.subplots(
        1, 2, figsize=(6.5, 3.1), layout="constrained", gridspec_kw={"width_ratios": [1.35, 1]}
    )
    colors = plt.colormaps["viridis"](np.linspace(0.05, 0.92, len(profiles)))
    samples = []
    for (alpha, _spin, high, r, response), color in zip(profiles, colors):
        axes[0].plot(r, response, color=color, label=f"{alpha:.2f}" + ("*" if high else ""))
        samples.append(np.interp(sample_radii, r, response))
    radius, exponent = window
    r = np.geomspace(1, 60, 600)
    axes[0].plot(r, np.exp(-((r / radius) ** exponent)), "k--", linewidth=1.1, label="$W(r)$")
    axes[0].axvspan(1, horizon_rg, color=".75", alpha=0.7, zorder=-2)
    axes[0].set(
        xscale="log",
        xlim=(1, 60),
        ylim=(0, 1.04),
        xlabel="$r$ ($r_g$)",
        ylabel=r"$R_{211}/R_{\max}$",
        title="(a) Radial profile and field window",
    )
    axes[0].set_xticks([1, 2, 5, 10, 20, 60], ["1", "2", "5", "10", "20", "60"])
    low = next(spin for _, spin, high, _, _ in profiles if not high)
    highs = sorted({spin for _, spin, high, _, _ in profiles if high})
    axes[0].legend(
        ncol=3,
        title=rf"$\alpha$; $a={low:g}$" + "".join(f", * $a={s:g}$" for s in highs),
        title_fontsize=8,
        fontsize=8,
        loc="upper left",
        columnspacing=0.6,
        handlelength=1.2,
    )
    samples = np.asarray(samples)
    alphas = [p[0] for p in profiles]
    for index, color in zip(range(min(3, len(sample_radii))), COLORS):
        axes[1].plot(
            alphas,
            samples[:, index],
            "o-",
            color=color,
            markersize=3,
            label=f"$r={sample_radii[index]:g}\\,r_g$",
        )
    if switch_alpha is not None:
        axes[1].axvline(switch_alpha, color=".7", linewidth=0.8, linestyle=":")
    axes[1].set(
        xlabel=r"$\alpha$",
        ylabel=r"$R_{211}/R_{\max}$",
        title="(b) Inner-cloud radial amplitude",
        ylim=(0, 1.04),
    )
    axes[1].legend(loc="upper left", frameon=False)
    for ax in axes:
        grid(ax)
    save(fig, out, "radial_function")
    return samples


def kappa(ratio, chi12, phi12):
    """Two-emitter combination factor; ratio=|P2/P1|, chi12 and phi12 in radians."""
    p1, p2 = 1.0, ratio * np.exp(2j * chi12)
    net = np.abs(p1 + p2) ** 2
    contrast = np.where(net > 1e-12, (1 - ratio**2) / np.where(net > 1e-12, net, 1.0), np.nan)
    return np.sqrt(np.cos(phi12 / 2) ** 2 + np.sin(phi12 / 2) ** 2 * contrast**2)


def kappa_direct(ratio, chi12, phi12, amplitude, phases):
    """First Fourier component of arg(P1 e^{2i dchi1} + P2 e^{2i dchi2})/2 at a small amplitude."""
    psi = 2 * np.pi * np.arange(phases) / phases
    p1, p2 = 1.0, ratio * np.exp(2j * chi12)
    total = p1 * np.exp(2j * amplitude * np.cos(psi)) + p2 * np.exp(2j * amplitude * np.cos(psi + phi12))
    base = 0.5 * np.angle(p1 + p2)
    delta = 0.5 * np.angle(total * np.exp(-2j * base))
    return 2 * np.abs(np.mean(delta * np.exp(-1j * psi))) / amplitude


def two_photon(out, left_ratio, right_phi_degrees):
    chi = np.linspace(0, np.pi / 2, 401)
    phi = np.linspace(0, np.pi, 401)
    ratio = np.linspace(0.01, 1.0, 400)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8), constrained_layout=True)
    c, f = np.meshgrid(chi, phi)
    r, c2 = np.meshgrid(ratio, chi)
    norm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=3)
    for ax, x, y, k in (
        (axes[0], np.degrees(chi), np.degrees(phi), kappa(left_ratio, c, f)),
        (axes[1], ratio, np.degrees(chi), kappa(r, c2, np.deg2rad(right_phi_degrees))),
    ):
        mesh = ax.pcolormesh(
            x, y, np.clip(k, 0, 3), cmap="RdBu_r", norm=norm, shading="auto", rasterized=True
        )
        ax.contour(x, y, k, levels=[1.0], colors="k", linewidths=0.8)
    axes[0].set_xlabel(r"EVPA difference $|\chi_{12}|$ (deg)")
    axes[0].set_ylabel(r"axion phase difference $|\phi_{12}|$ (deg)")
    axes[0].set_title(rf"$|P_2/P_1|={left_ratio:g}$", fontsize=9)
    axes[1].set_xlabel(r"$|P_2/P_1|$")
    axes[1].set_ylabel(r"EVPA difference $|\chi_{12}|$ (deg)")
    axes[1].set_title(rf"$|\phi_{{12}}|={right_phi_degrees:g}^\circ$", fontsize=9)
    cbar = fig.colorbar(mesh, ax=axes, shrink=0.9, extend="max")
    cbar.set_label(r"$\kappa=\mathcal{A}^{\rm tot}/\mathcal{A}$")
    save(fig, out, "two_photon")


def signal_maps(out, frame, amplitude_deg, phase, fov, tick):
    """frame: no-axion processed Stokes [x, y, IQUV]; tick: stride, I and P fractions, length."""
    i, q, u = (frame[..., k] for k in range(3))
    p = np.hypot(q, u)
    chi = 0.5 * np.arctan2(u, q)
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.25), layout="constrained")
    intensity = sky(axes[0], i / i.max(), fov, cmap="afmhot", vmin=0, vmax=1)
    sparse = np.zeros(i.shape, bool)
    sparse[:: tick["stride"], :: tick["stride"]] = True
    keep = sparse & (i > tick["intensity_fraction"] * i.max()) & (p > tick["polarization_fraction"] * p.max())
    n = i.shape[0]
    coords = (np.arange(n) - (n - 1) / 2) * fov / n
    east, north = np.meshgrid(-coords, coords, indexing="ij")
    length = tick["length_uas"] * p[keep] / p.max()
    dx, dy = 0.5 * length * np.sin(chi[keep]), 0.5 * length * np.cos(chi[keep])
    segments = np.stack(
        (
            np.column_stack((east[keep] - dx, north[keep] - dy)),
            np.column_stack((east[keep] + dx, north[keep] + dy)),
        ),
        axis=1,
    )
    axes[0].add_collection(LineCollection(segments, colors="#41C9F5", linewidths=0.9))
    a_plot = sky(axes[1], amplitude_deg, fov, cmap="viridis", vmin=0)
    phi_plot = sky(axes[2], phase, fov, cmap="twilight", vmin=-np.pi, vmax=np.pi)
    axes[0].set_ylabel("North offset (μas)")
    for ax, title in zip(
        axes, ("(a) No-axion intensity + EVPA", "(b) Response amplitude", "(c) Response phase")
    ):
        ax.set_title(title)
    colorbar(fig, axes[0], intensity, r"$I/I_{\max}$")
    colorbar(fig, axes[1], a_plot, "$A/c$ (deg)")
    cb = colorbar(fig, axes[2], phi_plot, r"$\phi$ (rad)")
    cb.set_ticks([-np.pi, 0, np.pi], labels=[r"$-\pi$", "0", r"$\pi$"])
    save(fig, out, "signal_maps")
    return int(keep.sum())


def template_maps(out, mean_deg, fraction, curve_deg, spread_deg, fov):
    fig, axes = plt.subplots(
        1, 3, figsize=(6.5, 3.25), layout="constrained", gridspec_kw={"width_ratios": [1, 1, 1.15]}
    )
    m_plot = sky(axes[0], mean_deg, fov, cmap="viridis", vmin=0)
    s_plot = sky(axes[1], fraction, fov, cmap="magma", norm=LogNorm(vmin=fraction.min(), vmax=fraction.max()))
    axes[0].set_ylabel("North offset (μas)")
    axes[0].set_title("(a) Mean response amplitude")
    axes[1].set_title("(b) Fractional time scatter")
    colorbar(fig, axes[0], m_plot, r"$\langle A\rangle/c$ (deg)")
    cb = colorbar(fig, axes[1], s_plot, r"$\sigma_A/\langle A\rangle$ (log scale)")
    ticks = [t for t in (0.1, 0.3, 1, 3) if fraction.min() <= t <= fraction.max()]
    cb.set_ticks(ticks, labels=[f"{t:g}" for t in ticks])
    cb.ax.minorticks_off()
    theta = np.arange(len(curve_deg)) * 360 / len(curve_deg)
    axes[2].plot(theta, curve_deg, color=COLORS[0])
    axes[2].fill_between(
        theta,
        curve_deg - spread_deg,
        curve_deg + spread_deg,
        color=COLORS[0],
        alpha=0.25,
        label="±1σ over time",
    )
    axes[2].set(
        xlabel="Position angle (deg E of N)",
        ylabel=r"$\langle A\rangle/c$ (deg)",
        title="(c) 1D aperture template",
        xlim=(0, 360),
        ylim=(0, None),
    )
    axes[2].set_xticks([0, 90, 180, 270, 360])
    axes[2].tick_params(axis="x", labelsize=8, rotation=45)
    axes[2].legend(frameon=False, loc="lower right")
    grid(axes[2])
    save(fig, out, "template_maps")


def background(out, selected, unmasked):
    """selected/unmasked: {dim: differences in degrees}; selected also carries mean and sigma."""
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.0), layout="constrained")
    for ax, dim, panel in zip(axes, ("1d", "2d"), ("a", "b")):
        data, mu, sigma = selected[dim]
        everything = unmasked[dim]
        edges = np.linspace(everything.min(), everything.max(), 101)
        counts, _ = np.histogram(everything, edges)
        ax.stairs(
            counts / (everything.size * np.diff(edges)), edges, color=".5", linewidth=1, label="All features"
        )
        ax.hist(
            data.ravel(),
            bins=40,
            density=True,
            color=COLORS[0],
            alpha=0.28,
            edgecolor=COLORS[0],
            linewidth=0.6,
            label="Selected inference features",
        )
        x = np.linspace(min(-12, data.min()), max(12, data.max()), 500)
        ax.plot(
            x,
            np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma),
            color=COLORS[1],
            label="Gaussian, selected pooled σ",
        )
        ax.set(
            xlabel="No-axion ΔEVPA (deg)",
            ylabel="Probability density (deg$^{-1}$)",
            title=f"({panel}) {dim.upper()}: selected σ = {sigma:.3f}°",
        )
        ax.set_xlim(min(-10, data.min() - 0.5), max(10, data.max() + 0.5))
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        grid(ax)
    save(fig, out, "background")


def observational_uncertainty(out, entries):
    """entries[dim]: dict(x, sigma_deg [2, n], correlation [2, n, n], separation [n, n] or None)."""
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.1), layout="constrained")
    for ax, dim, xlabel in (
        (axes[0], "1d", "Position angle (deg E of N)"),
        (axes[1], "2d", "Selected pixel index"),
    ):
        entry = entries[dim]
        for pair, label, color in zip(range(2), PAIRS, COLORS):
            ax.plot(entry["x"], entry["sigma_deg"][pair], "o-", markersize=2.5, color=color, label=label)
        ax.set(xlabel=xlabel, ylabel="D-term pair σ (deg)")
        ax.legend(frameon=False)
        grid(ax)
    entry = entries["2d"]
    iu = np.triu_indices(len(entry["x"]), 1)
    for pair, label, color in zip(range(2), PAIRS, COLORS):
        axes[2].scatter(
            entry["separation"][iu],
            np.abs(entry["correlation"][pair][iu]),
            s=7,
            facecolors="none",
            edgecolors=color,
            linewidths=0.55,
            alpha=0.55,
            label=label,
            rasterized=False,
        )
    axes[2].set(xlabel="Pixel separation (μas)", ylabel=r"$|\rho_{ij}|$", ylim=(-0.025, 1.025))
    grid(axes[2])
    for ax, title in zip(
        axes, ("(a) 1D selected apertures", "(b) 2D selected pixels", "(c) 2D pair correlations")
    ):
        ax.set_title(title)
    save(fig, out, "observational_uncertainty")


def data_vs_model(out, entries, alpha, sector_spacing):
    """entries[dim]: dict(x, data_deg [2, n], sigma_deg [2, n], prediction_deg [2, n])."""
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.5), layout="constrained", sharey=True)
    fit_label = f"Saved grid fit, α = {alpha:.2f}"
    for column, dim in enumerate(("1d", "2d")):
        entry = entries[dim]
        x, n = entry["x"], len(entry["x"])
        for pair, label, ax in zip(range(2), PAIRS, axes[:, column]):
            ax.axhline(0, color=".5", linewidth=0.8)
            ax.errorbar(
                x,
                entry["data_deg"][pair],
                yerr=entry["sigma_deg"][pair],
                fmt="o",
                markersize=3,
                capsize=1.8,
                linewidth=0.8,
                color=COLORS[0],
                label="Observed ± marginal σ",
            )
            if dim == "1d":
                # Do not bridge gaps between non-adjacent selected sectors.
                gaps = np.flatnonzero(np.diff(x) > sector_spacing) + 1
                for k, indices in enumerate(np.split(np.arange(n), gaps)):
                    ax.plot(
                        x[indices],
                        entry["prediction_deg"][pair, indices],
                        "s-",
                        markersize=2.5,
                        color=COLORS[1],
                        label=fit_label if k == 0 else None,
                    )
            else:
                ax.plot(x, entry["prediction_deg"][pair], "s", markersize=3, color=COLORS[1], label=fit_label)
            ax.set_title(f"{dim.upper()}, {label} ({n} features)")
            ax.set_ylabel("ΔEVPA (deg)")
            grid(ax)
            if pair == 1:
                ax.set_xlabel("Position angle (deg E of N)" if dim == "1d" else "Selected pixel index")
        axes[0, column].legend(frameon=False, fontsize=8)
    save(fig, out, "data_vs_model")


def limits(out, results, features):
    """results[(data, mode, dim)] = (alphas, c95); data in {'observed', 'zero-data'}."""
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.3), layout="constrained", sharey=True)
    for ax, dim in zip(axes, ("1d", "2d")):
        for mode, color in zip(MODES, COLORS):
            for data, style in (("observed", "-"), ("zero-data", "--")):
                alphas, c95 = results[(data, mode, dim)]
                ax.plot(
                    alphas,
                    c95,
                    linestyle=style,
                    color=color,
                    marker="o" if data == "observed" else None,
                    markersize=3,
                )
        ax.set(yscale="log", xlabel=r"$\alpha$", title=f"{dim.upper()}: {features[dim]} selected features")
        ax.set_xticks([0.1, 0.2, 0.3, 0.4])
        grid(ax)
    axes[0].set_ylabel("$c_{95}$ (dimensionless)")
    axes[0].legend(
        [Line2D([], [], color=c) for c in COLORS], MODE_LABELS, loc="upper right", frameon=False, fontsize=8
    )
    axes[1].legend(
        [Line2D([], [], color=".2", linestyle=s) for s in ("-", "--")],
        ("Observed day differences", "Zero-data control"),
        loc="upper right",
        frameon=False,
    )
    save(fig, out, "limits")
