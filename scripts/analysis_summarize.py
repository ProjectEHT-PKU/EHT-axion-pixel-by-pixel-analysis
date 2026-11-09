#!/usr/bin/env python3
"""Process 3, step 5: summary tables/figures of the limits and their checks.

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
import numpy as np

from eht_axion.pipeline import common

MODES = ("std", "cor", "cor_fixed_background")


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
            for center, style in (("observed", "-"), ("null-comparison", "--")):
                row = report["results"][f"{center}-{mode}-{dim}"]
                ax.plot(
                    row["alphas"],
                    row["c95"],
                    style,
                    color=color,
                    marker="o",
                    label=f"{mode}: {center}",
                )
        ax.set(
            xlabel=r"$\alpha$",
            ylabel=r"$c_{95}$ (model coupling)",
            yscale="log",
            title=dim.upper(),
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle("Conditional observed credible limits / matched NULL controls")
    fig.savefig(out / "limits-observed-null.png", dpi=160)
    plt.close(fig)


def fit_evidence(source):
    """Extract saved fits; no optimization or noise regularization here."""
    rows, diagnostics = [], {}
    for dim in ("1d", "2d"):
        for mode in MODES:
            key = f"observed-{mode}-{dim}"
            diagnostic = json.loads((source / (key + ".json")).read_text())
            index = np.flatnonzero(
                np.isclose(diagnostic["alphas"], 0.3, rtol=0, atol=1e-12)
            )
            if len(index) != 1:
                raise ValueError("Exactly one alpha=0.3 fit required")
            fit = diagnostic["fits"][int(index[0])]
            rows.append(
                {
                    "dimension": dim,
                    "mode": mode,
                    "alpha": 0.3,
                    "components": 2 * len(diagnostic["selected_indices"]),
                    "null_chi2": fit["data_quadratic"],
                    "best_chi2": fit["joint_grid_chi2"],
                    "best_c": fit["joint_grid_coupling"],
                    "phase_radians": fit["joint_grid_phase_radians"],
                }
            )
            diagnostics[key] = diagnostic
    return rows, diagnostics


def plot_fits(out, diagnostics):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(16, 12), layout="constrained")
    for column, dim in enumerate(("1d", "2d")):
        for row, mode in enumerate(MODES):
            diagnostic = diagnostics[f"observed-{mode}-{dim}"]
            index = int(np.flatnonzero(np.isclose(diagnostic["alphas"], 0.3))[0])
            fit = diagnostic["fits"][index]
            data = np.asarray(diagnostic["data_radians"])
            prediction = np.asarray(fit["prediction_radians"])
            covariance = np.asarray(diagnostic["effective_pair_covariances"])
            residual = data - prediction
            white = np.concatenate(
                [
                    np.linalg.solve(np.linalg.cholesky(c), r)
                    for c, r in zip(covariance, residual)
                ]
            )
            ax = axes[row, column]
            n = data.shape[1]
            for pair, label, color in ((0, "96-95", "C0"), (1, "101-100", "C1")):
                x = np.arange(n) + pair * (n + 2)
                ax.errorbar(
                    x,
                    np.rad2deg(data[pair]),
                    yerr=np.rad2deg(np.sqrt(np.diag(covariance[pair]))),
                    fmt=".",
                    color=color,
                    label=f"{label}: data / pair scatter",
                )
                ax.plot(
                    x,
                    np.rad2deg(prediction[pair]),
                    "-",
                    color=color,
                    label=f"{label}: saved shared-phase best fit",
                )
            ax.set(
                title=(
                    f"{dim.upper()} {mode}; "
                    f"whitened residual RMS={np.sqrt(np.mean(white**2)):.2f}"
                ),
                xlabel="Selected feature index (pairs stacked; not spatial distance)",
                ylabel="Day EVPA difference (degrees)",
            )
            ax.legend(fontsize=7)
            ax.grid(alpha=0.2)
    fig.suptitle(
        "Alpha=0.30: saved joint-grid fits, no refit; effective pair covariance"
    )
    fig.savefig(out / "alpha030-fits.png", dpi=160)
    plt.close(fig)


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    destination = args.output_root.resolve()
    if destination.exists():
        raise FileExistsError("Use a new summary output root")
    sources = {"primary": args.limits.resolve()}
    sources.update({name: getattr(args, name).resolve() for name in ("gridcheck", "priorcheck")})
    if getattr(args, "covcheck", None) is not None:
        sources["covcheck"] = args.covcheck.resolve()
    # Fail fast if a source has not committed its manifest. No waiting/monitor loop.
    for source in sources.values():
        common.complete_manifest(source)
    reports, consumed, records = {}, [], {}
    for name, source in sources.items():
        record = common.verify_stage(source, recursive=False)
        expected = (
            "research:eht120-limit-numerics"
            if name in ("gridcheck", "priorcheck")
            else "research:eht120-observed-limits"
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
        common.require_declared_outputs(source, paths)
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
    for name in ("gridcheck", "priorcheck"):
        rows = comparisons(primary, reports[name])
        if any(r["features_reference"] != r["features_test"] for r in rows):
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
    # Byte-identical receipt snapshots, not canonical stage rebinding. Their
    # filenames intentionally avoid common's recursive manifest-name dispatch.
    receipts = destination.with_name(destination.name + "-source-receipts")
    receipts.mkdir(parents=True, exist_ok=False)
    snapshots = []
    for name, source in sources.items():
        original = common.manifest_path(source)
        snapshot = receipts / (name + "-receipt.json")
        snapshot.write_bytes(original.read_bytes())
        snapshots.append(
            {
                "original": str(original),
                "snapshot": str(snapshot),
                "sha256": common.sha256(snapshot),
            }
        )
        consumed.append(snapshot)
    consumed.append(Path(__file__).resolve())
    identity = {
        "schema": "eht120-results-summary-v1",
        "sources": {k: str(v) for k, v in sources.items()},
        "receipt_snapshots": snapshots,
        "verification": (
            "completed small stage outputs only; raw lineage inherited, not rehashed"
        ),
        "canonical_NULL_modified": False,
        "calibration_certification": False,
    }
    with common.stage_directory(
        destination, "research:eht120-results-summary", identity, consumed
    ) as out:
        write_csv(out / "limits.csv", limits)
        write_csv(out / "alpha030-chi2.csv", fits)
        for label, rows in comparisons_all.items():
            write_csv(out / (label + ".csv"), rows)
        evidence = dict(
            identity,
            results=primary["results"],
            comparisons=comparisons_all,
            alpha030_chi2=fits,
            background_sigma_radians=primary["background_sigma_radians"],
            finite_mc_covariance_check=reports.get("covcheck"),
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
        plot_fits(out, diagnostics)
        lines = [
            "# EHT120 completed result comparison",
            "",
            ("Conditional 95% credible limits in model coupling c; "
            "not GeV conversion or certified exclusions."),
            "",
            "## Numerical comparisons",
            "",
            "Prior check is compared to the dense grid, not the primary coarse grid.",
        ]
        for label, rows in comparisons_all.items():
            lines.append(
                f"- {label}: max |delta log10 c95| = "
                f"{max(abs(r['delta_dex']) for r in rows):.6g} dex."
            )
        lines += [
            "",
            "## Alpha=0.30 goodness of fit",
            "",
            "| Dimension | Mode | Components | Null chi2 | Best joint-grid chi2 |",
            "|---|---|---:|---:|---:|",
        ]
        lines += [
            f"| {r['dimension']} | {r['mode']} | {r['components']} | "
            f"{r['null_chi2']:.3f} | {r['best_chi2']:.3f} |"
            for r in fits
        ]
        lines += [
            "",
            ("Correlated D-only fits are poor: "
            "no formal new-physics detection/exclusion follows."),
            ("Fixed-background mode adds selected-feature pooled pair variance, "
            "not signal-response variance."),
            ("No eigenvalue repair, diagonal ridge, refit or EVPA alignment "
            "was applied by this helper."),
            ("Observed centers away from branch seams do not validate "
            "every prior-predictive Gaussian tail."),
            "",
            "![Limits](limits-observed-null.png)",
            "",
            "![Saved fits](alpha030-fits.png)",
        ]
        (out / "RESULTS.md").write_text("\n".join(lines) + "\n")
    return destination


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limits", type=Path, required=True, help="analysis_limits.py output root")
    for name in ("gridcheck", "priorcheck", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument(
        "--covcheck",
        type=Path,
        help="Optional completed finite-MC covariance sensitivity stage",
    )
    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
