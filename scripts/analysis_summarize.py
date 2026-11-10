#!/usr/bin/env python3
"""Analysis: summary tables/figures of the limits and their checks.

Inputs are completed small stages only (limits of analysis_limits.py, the grid and
prior checks of analysis_check_grid.py, optionally the covariance check); raw
inputs are never reopened and their provenance is inherited, not re-certified.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import _bootstrap  # noqa: F401
from _data_root import add_data_root_option, data_root, fill_defaults, layout_path
import numpy as np

from eht_axion import provenance
from eht_axion.analysis.likelihood import MODES


def comparisons(primary, other):
    """Compare identical alpha/case grids, without attributing changed selections."""
    rows = []
    if set(primary["results"]) != set(other["results"]):
        raise ValueError("Case grids differ")
    for key, base in primary["results"].items():
        test = other["results"][key]
        if not np.array_equal(base["alphas"], test["alphas"]):
            raise ValueError("Alpha grids differ")
        delta = np.asarray(test["log10_c95"]) - np.asarray(base["log10_c95"])
        if not np.isfinite(delta).all():
            raise ValueError("Nonfinite limits")
        for alpha, value in zip(base["alphas"], delta):
            rows.append(
                {
                    "case": key,
                    "alpha": alpha,
                    "delta_dex": float(value),
                    "ratio": float(10**value),
                    "features_reference": base["features"],
                    "features_test": test["features"],
                }
            )
    return rows


def plot_limits(out, report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), layout="constrained")
    for ax, dim in zip(axes, ("1d", "2d")):
        for mode, color in zip(MODES, ("C0", "C1", "C2")):
            for data, style in (("observed", "-"), ("zero-data", "--")):
                row = report["results"][f"{data}-{mode}-{dim}"]
                ax.plot(
                    row["alphas"],
                    row["c95"],
                    style,
                    color=color,
                    marker="o",
                    label=f"{mode}: {data}",
                )
        ax.set(
            xlabel=r"$\alpha$",
            ylabel=r"$c_{95}$ (model coupling)",
            yscale="log",
            title=dim.upper(),
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle("Credible limits: observed data and zero-data control")
    fig.savefig(out / "limits.png", dpi=160)
    plt.close(fig)


def fit_evidence(source):
    """Saved joint-grid best fits of every observed case and alpha; no refit."""
    rows, diagnostics = [], {}
    for dim in ("1d", "2d"):
        for mode in MODES:
            key = f"observed-{mode}-{dim}"
            diagnostic = json.loads((source / (key + ".json")).read_text())
            if len(diagnostic["fits"]) != len(diagnostic["alphas"]):
                raise ValueError(f"{key}: one saved fit per alpha required")
            for alpha, fit in zip(diagnostic["alphas"], diagnostic["fits"]):
                rows.append(
                    {
                        "dimension": dim,
                        "mode": mode,
                        "alpha": alpha,
                        "components": 2 * len(diagnostic["selected_indices"]),
                        "chi2_zero_coupling": fit["data_quadratic"],
                        "chi2_best": fit["joint_grid_chi2"],
                        "best_c": fit["joint_grid_coupling"],
                        "best_phase_radians": fit["joint_grid_phase_radians"],
                    }
                )
            diagnostics[key] = diagnostic
    return rows, diagnostics


def fit_figure_name(alpha):
    return f"fits-alpha-{alpha:.2f}.png"


def plot_fits(out, diagnostics):
    """One figure per alpha: data, saved best fit and whitened residual RMS of every case."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    alphas = diagnostics[f"observed-{MODES[0]}-1d"]["alphas"]
    names = []
    for index, alpha in enumerate(alphas):
        fig, axes = plt.subplots(len(MODES), 2, figsize=(16, 4 * len(MODES)), layout="constrained",
                                 squeeze=False)
        for column, dim in enumerate(("1d", "2d")):
            for row, mode in enumerate(MODES):
                diagnostic = diagnostics[f"observed-{mode}-{dim}"]
                fit = diagnostic["fits"][index]
                data = np.asarray(diagnostic["data_radians"])
                prediction = np.asarray(fit["prediction_radians"])
                covariance = np.asarray(diagnostic["effective_pair_covariances"])
                white = np.concatenate([np.linalg.solve(np.linalg.cholesky(c), r)
                                        for c, r in zip(covariance, data - prediction)])
                ax = axes[row, column]
                n = data.shape[1]
                for pair, label, color in ((0, "96-95", "C0"), (1, "101-100", "C1")):
                    x = np.arange(n) + pair * (n + 2)
                    ax.errorbar(x, np.rad2deg(data[pair]),
                                yerr=np.rad2deg(np.sqrt(np.diag(covariance[pair]))),
                                fmt=".", color=color, label=f"{label}: data, pair scatter")
                    ax.plot(x, np.rad2deg(prediction[pair]), "-", color=color,
                            label=f"{label}: saved best fit")
                ax.set(title=f"{dim.upper()} {mode}; whitened residual RMS={np.sqrt(np.mean(white**2)):.2f}",
                       xlabel="Selected feature index (pairs stacked; not spatial distance)",
                       ylabel="Day EVPA difference (degrees)")
                ax.legend(fontsize=7)
                ax.grid(alpha=0.2)
        fig.suptitle(f"alpha = {alpha:.2f}: saved joint-grid best fits (no refit)")
        name = fit_figure_name(alpha)
        fig.savefig(out / name, dpi=120)
        plt.close(fig)
        names.append(name)
    return names


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def results_text(primary, comparisons_all, fits, covcheck, figures):
    """RESULTS.md: the configured credibility and the numbers of this run; no fixed judgements."""
    credibility = primary["credibility"]
    lines = [
        "# Results summary",
        "",
        f"Upper limits c95 are the {100 * credibility:g}% posterior quantiles of the coupling c "
        "(analysis.inference.credibility), conditional on the grid and prior of analysis.inference.",
        "",
        "## Limits",
        "",
        "| Case | Features | " + " | ".join(f"alpha {a:.2f}" for a in
                                            primary["results"]["observed-std-1d"]["alphas"]) + " |",
        "|---|---:|" + "---:|" * len(primary["results"]["observed-std-1d"]["alphas"]),
    ]
    lines += [f"| {key} | {result['features']} | " + " | ".join(f"{c:.4g}" for c in result["c95"]) + " |"
              for key, result in primary["results"].items()]
    lines += ["", "## GRMHD background", ""]
    lines += [f"- {dim}: sigma_bg = {sigma:.6g} rad ({np.rad2deg(sigma):.4g} deg)"
              for dim, sigma in primary["background_sigma_radians"].items()]
    lines += ["", "## Numerical checks", ""]
    for label, rows in comparisons_all.items():
        lines.append(f"- {label}: max |delta log10 c95| = {max(abs(r['delta_dex']) for r in rows):.6g} dex")
    if covcheck is not None:
        lines.append(f"- covariance halves: max |delta log10 c95| = {covcheck['max_abs_delta_dex']:.6g} dex")
    lines += [
        "",
        "## Saved joint-grid best fits",
        "",
        "chi2 at c = 0 and at the best grid point (c, phase); components = 2 x features.",
        "",
        "| Dimension | Mode | alpha | Components | chi2 (c = 0) | chi2 (best) | best c |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {r['dimension']} | {r['mode']} | {r['alpha']:.2f} | {r['components']} | "
        f"{r['chi2_zero_coupling']:.3f} | {r['chi2_best']:.3f} | {r['best_c']:.4g} |"
        for r in fits
    ]
    lines += ["", "![Limits](limits.png)", ""]
    lines += [f"![Saved fits]({name})" for name in figures]
    return lines


def run(args):
    destination = args.output_root.resolve()
    if destination.exists():
        raise FileExistsError("Use a new summary output root")
    sources = {"primary": args.limits.resolve()}
    sources.update({name: getattr(args, name).resolve() for name in ("gridcheck", "priorcheck")})
    if getattr(args, "covcheck", None) is not None:
        sources["covcheck"] = args.covcheck.resolve()
    reports, consumed, records = {}, [], {}
    for name, source in sources.items():
        record = provenance.verify_stage(source, recursive=False)
        expected = (
            "analysis:check-grid"
            if name in ("gridcheck", "priorcheck")
            else "analysis:limits"
        )
        if name != "covcheck" and record["stage"] != expected:
            raise ValueError(f"Unexpected {name} stage")
        records[name] = record
        paths = [source / "report.json"]
        if name == "primary":
            paths += [
                source / f"observed-{mode}-{dim}.json"
                for dim in ("1d", "2d")
                for mode in MODES
            ]
            paths += [source / f"observations-{dim}.json" for dim in ("1d", "2d")]
        provenance.require_declared_outputs(source, paths)
        consumed += paths
        reports[name] = json.loads(paths[0].read_text())
    primary = reports["primary"]
    comparisons_all = {}
    for label, reference, test in (
        ("grid_vs_primary", "primary", "gridcheck"),
        ("prior_vs_dense_grid", "gridcheck", "priorcheck"),
    ):
        comparisons_all[label] = comparisons(reports[reference], reports[test])
    # Grid and prior changes must not silently include a reselection.
    for label in ("grid_vs_primary", "prior_vs_dense_grid"):
        if any(r["features_reference"] != r["features_test"] for r in comparisons_all[label]):
            raise ValueError("Numerical check changed selection size")
    fits, diagnostics = fit_evidence(sources["primary"])
    limits = []
    for key, result in primary["results"].items():
        for a, c, logc in zip(result["alphas"], result["c95"], result["log10_c95"]):
            limits.append(
                {
                    "case": key,
                    "alpha": a,
                    "c95": c,
                    "log10_c95": logc,
                    "features": result["features"],
                }
            )
    # The source manifests are recorded inputs; their lineage is verified when
    # those stages are made and is not re-hashed here.
    consumed += [provenance.manifest_path(source) for source in sources.values()]
    consumed.append(Path(__file__).resolve())
    identity = {
        "schema": "analysis-summary-v2",
        "sources": {k: str(v) for k, v in sources.items()},
        "verification": (
            "completed small stage outputs only; raw lineage inherited, not rehashed"
        ),
        "calibration_certification": False,
    }
    with provenance.stage_directory(
        destination, "analysis:summary", identity, consumed, verify_lineage=False
    ) as out:
        write_csv(out / "limits.csv", limits)
        write_csv(out / "fit-chi2.csv", fits)
        for label, rows in comparisons_all.items():
            write_csv(out / (label + ".csv"), rows)
        evidence = dict(
            identity,
            credibility=primary["credibility"],
            results=primary["results"],
            comparisons=comparisons_all,
            fit_chi2=fits,
            background_sigma_radians=primary["background_sigma_radians"],
            covariance_check=reports.get("covcheck"),
            diagnostics={
                k: {
                    field: v[field]
                    for field in (
                        "branch_risk",
                        "working_covariance_diagnostics",
                        "boundaries",
                    )
                    if field in v
                }
                for k, v in diagnostics.items()
            },
            observations={
                dim: json.loads(
                    (sources["primary"] / f"observations-{dim}.json").read_text()
                )
                for dim in ("1d", "2d")
            },
        )
        (out / "report.json").write_text(
            json.dumps(evidence, indent=2, allow_nan=False) + "\n"
        )
        plot_limits(out, primary)
        figures = plot_fits(out, diagnostics)
        lines = results_text(primary, comparisons_all, fits, reports.get("covcheck"), figures)
        (out / "RESULTS.md").write_text("\n".join(lines) + "\n")
    return destination


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limits", type=Path,
                        help="analysis_limits.py output root (default: analysis/limits)")
    parser.add_argument("--gridcheck", type=Path, help="default: analysis/check-grid")
    parser.add_argument("--priorcheck", type=Path, help="default: analysis/check-prior")
    parser.add_argument("--output-root", type=Path, help="new directory (default: analysis/summary)")
    parser.add_argument(
        "--covcheck",
        type=Path,
        help="Completed covariance sensitivity stage; optional without a data root "
             "(default with a data root: analysis/check-covariance)",
    )
    add_data_root_option(parser)
    args = parser.parse_args(argv)
    defaults = {"limits": ("analysis_limits", {}), "gridcheck": ("analysis_check_grid", {}),
                "priorcheck": ("analysis_check_prior", {}), "output_root": ("analysis_summary", {})}
    root = data_root(args)
    if args.covcheck is None and root is not None:
        args.covcheck = layout_path(root, "analysis_check_covariance")
    fill_defaults(parser, args, defaults)
    run(args)


if __name__ == "__main__":
    main()
