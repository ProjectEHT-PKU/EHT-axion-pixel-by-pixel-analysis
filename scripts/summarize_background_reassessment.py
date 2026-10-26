"""Postprocess completed background reassessment only; never run inference.

No automatic winner, raw ancestry traversal, or waiting for optional experiments.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import shutil
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.pipeline.common import (
    manifest_path,
    require_declared_outputs,
    sha256,
    stage_directory,
    verification_session,
    verify_stage,
)
from eht_axion.statistics import normalized_posterior_cdf, posterior_upper_limit

BASE = Path("/home/canoziia/GRMHD_DATA/analysis")
ALPHAS = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]
MODELS = ("iso", "temporal_iso", "partial_diag", "shrink_full_0.25")
RESPONSES = ("complex_mean", "frozen_mixture")
MODES = ("observed", "matched_zero")
PLOTS = ("cv-nll.png", "alpha030-c95.png", "allalpha-c95.png", "raw-full-failure.png")


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def read_rows(path):
    rows = read_csv(path)
    for r in rows:
        for key in ("alpha", "c95", "fixed_cov_condition"):
            r[key] = float(r[key])
        for key in ("contrasts", "state_count"):
            r[key] = int(r[key])
    return rows


def index_rows(rows, alphas=None):
    idx = {}
    for r in rows:
        k = tuple(r[x] for x in ("dimension", "alpha", "contrasts", "background", "response", "mode"))
        if k in idx:
            raise ValueError(f"duplicate case {k}")
        if not np.isfinite(r["c95"]) or r["c95"] <= 0:
            raise ValueError("invalid c95")
        expected_key = f"{k[0]}-a{k[1]:.6f}-h{k[2]}-{k[3]}-{k[4]}-{k[5]}"
        if r["key"] != expected_key:
            raise ValueError("case filename/key mismatch")
        if r["state_count"] != (141 if r["response"] == "frozen_mixture" else 1):
            raise ValueError("response state count mismatch")
        idx[k] = r
    if alphas is not None:
        expected = set(itertools.product(("1d", "2d"), alphas, (2, 3), MODELS, RESPONSES, MODES))
        if set(idx) != expected:
            raise ValueError("case counts/keys differ")
    return idx


def validate_posterior(path, row, metadata, *, require_phase=True):
    with np.load(path, allow_pickle=False) as z:
        c, ll, density, cdf = [z[k] for k in ("couplings", "loglike", "density", "cdf")]
        if c.ndim != 1 or any(x.shape != c.shape for x in (ll, density, cdf)):
            raise ValueError("posterior axis mismatch")
        if (
            len(c) != metadata["coupling_points"]
            or not np.isfinite(c).all()
            or np.any(c <= 0)
            or np.any(np.diff(c) <= 0)
            or not np.isfinite(ll).all()
        ):
            raise ValueError("invalid coupling/loglike grid")
        if require_phase and z["phase_loglike"].shape != (len(c), metadata["phase_count"]):
            raise ValueError("phase axis mismatch")
        if (
            not np.isfinite(density).all()
            or np.any(density < 0)
            or not np.isfinite(cdf).all()
            or np.any(np.diff(cdf) < -1e-12)
            or not np.allclose(cdf[[0, -1]], [0, 1], atol=1e-12)
        ):
            raise ValueError("invalid density/CDF")
        np.testing.assert_allclose(c[[0, -1]], metadata["prior_bounds"], rtol=1e-12)
        pd, pc = normalized_posterior_cdf(c, ll, prior="log-uniform")
        np.testing.assert_allclose(density, pd, rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(cdf, pc, rtol=1e-11, atol=1e-13)
        mass = np.sum(np.diff(c) * (density[1:] + density[:-1]) / 2)
        if not np.isclose(mass, 1, rtol=1e-11):
            raise ValueError("posterior not normalized")
        limit = posterior_upper_limit(c, ll, prior="log-uniform")
        if not np.allclose([row["c95"], float(z["c95"])], limit, rtol=1e-11, atol=0):
            raise ValueError("CSV/NPZ c95 mismatch")
        peak = int(np.argmax(ll))
        return dict(
            key=row["key"],
            sha256=sha256(path),
            normalized_mass=float(mass),
            lower_bin_mass=float(cdf[1]),
            upper_bin_mass=float(1 - cdf[-2]),
            likelihood_peak_side="lower" if peak == 0 else "upper" if peak == len(c) - 1 else "interior",
            coupling_points=len(c),
            phase_count=metadata["phase_count"],
            phase_axis_checked=require_phase,
        )


def cv_records(path):
    records = read_csv(path)
    for r in records:
        r["scoreable"] = r["scoreable"] == "True"
        for k in ("contrasts", "n_splits", "fold", "n_train", "n_test", "dimension", "rank"):
            r[k] = int(r[k])
        for k in ("mean_nll_per_dimension", "condition_number", "train_mean_rms"):
            r[k] = float(r[k]) if r.get(k) else None
        if r.get("train_mean"):
            r["train_mean_norm"] = float(np.linalg.norm(json.loads(r["train_mean"])))
        if (
            r.get("covariance_center") != "training_mean"
            or r.get("predictive_mean") != "zero"
            or r.get("ddof") != "0"
        ):
            raise ValueError("CV centering policy mismatch")
    return records


def aggregate_cv(records):
    """All folds required; never aggregate a scoreable subset of raw-full folds."""
    result = []
    groups = sorted({(r["dimension_label"], r["contrasts"], r["n_splits"]) for r in records})
    for d, h, count in groups:
        group = [r for r in records if (r["dimension_label"], r["contrasts"], r["n_splits"]) == (d, h, count)]
        for model in MODELS + ("diagonal", "full"):
            rr = sorted([r for r in group if r["model"] == model], key=lambda r: r["fold"])
            if [r["fold"] for r in rr] != list(range(count)):
                raise ValueError("missing/duplicate CV fold")
            valid = all(r["scoreable"] for r in rr)
            item = dict(
                dimension=d,
                contrasts=h,
                n_splits=count,
                model=model,
                scoreable_folds=sum(r["scoreable"] for r in rr),
                all_folds_scoreable=valid,
                weighted_nll=None,
                delta_to_iso=None,
                wins_vs_iso=None,
                worst_fold_nll=None,
            )
            if valid:
                scores = np.array([r["mean_nll_per_dimension"] for r in rr], dtype=float)
                iso = sorted([r for r in group if r["model"] == "iso"], key=lambda r: r["fold"])
                if not np.isfinite(scores).all():
                    raise ValueError("scoreable CV lacks finite NLL")
                item.update(
                    weighted_nll=float(np.average(scores, weights=[r["n_test"] for r in rr])),
                    worst_fold_nll=float(scores.max()),
                )
                if len(iso) == count and all(r["scoreable"] for r in iso):
                    delta = scores - [r["mean_nll_per_dimension"] for r in iso]
                    item.update(
                        delta_to_iso=float(np.average(delta, weights=[r["n_test"] for r in rr])),
                        wins_vs_iso=int(np.sum(delta < 0)),
                    )
            result.append(item)
    return result


def validate_cv(records, summaries, splits):
    expected = {(d, h, n) for d in ("1d", "2d") for h in (2, 3) for n in (5, 3)}
    if {(r["dimension_label"], r["contrasts"], r["n_splits"]) for r in records} != expected or len(
        records
    ) != 192:
        raise ValueError("CV grid mismatch")
    result = aggregate_cv(records)
    for d, h, count in expected:
        key = f"{d}-h{h}-folds{count}"
        folds = splits[key]
        if len(folds) != count:
            raise ValueError("split count mismatch")
        tests = []
        for i, f in enumerate(folds):
            train, test, purged = (np.array(f[x], dtype=int) for x in ("train", "test", "purged"))
            if sorted(np.concatenate((train, test, purged))) != list(range(124)):
                raise ValueError("split partition mismatch")
            support = np.array(f["support_indices"])
            lo, hi = f["test_support"]
            if np.any((support[train, 1] >= lo) & (support[train, 0] <= hi)):
                raise ValueError("shared raw support after purge")
            tests.extend(test.tolist())
            for r in records:
                if (r["dimension_label"], r["contrasts"], r["n_splits"], r["fold"]) == (d, h, count, i):
                    if (r["n_train"], r["n_test"]) != (len(train), len(test)):
                        raise ValueError("CV/splits size mismatch")
        if sorted(tests) != list(range(124)):
            raise ValueError("held-out windows mismatch")
        for r in result:
            if (r["dimension"], r["contrasts"], r["n_splits"]) != (d, h, count):
                continue
            original = summaries[key]["summary"][r["model"]]
            if original["all_folds_scoreable"] != r["all_folds_scoreable"]:
                raise ValueError("CV summary availability mismatch")
            if r["all_folds_scoreable"]:
                np.testing.assert_allclose(
                    r["weighted_nll"], original["weighted_mean_nll_per_dimension"], rtol=1e-12
                )
    return result


def make_plots(rows, records, aggregates, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    idx = index_rows(rows)
    colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for i, d in enumerate(("1d", "2d")):
        for j, h in enumerate((2, 3)):
            ax = axes[i, j]
            for count, offset, marker in ((5, -0.1, "o"), (3, 0.1, "s")):
                for x, model in enumerate(MODELS + ("diagonal",)):
                    rr = [
                        r
                        for r in records
                        if (r["dimension_label"], r["contrasts"], r["n_splits"], r["model"])
                        == (d, h, count, model)
                    ]
                    yy = [r["mean_nll_per_dimension"] for r in rr if r["scoreable"]]
                    ax.scatter(np.full(len(yy), x + offset), yy, s=15, alpha=0.4, color=colors[j])
                    a = next(
                        r
                        for r in aggregates
                        if (r["dimension"], r["contrasts"], r["n_splits"], r["model"]) == (d, h, count, model)
                    )
                    if a["all_folds_scoreable"]:
                        ax.plot(
                            x + offset,
                            a["weighted_nll"],
                            marker,
                            color="black",
                            label=f"{count}-fold weighted" if x == 0 else None,
                        )
            ax.set(title=f"{d} H{h} (compare within panel only)", ylabel="GRMHD NLL / dimension")
            ax.set_xticks(range(5), ["iso", "time iso", "partial diag", "shrink .25", "diag*"], rotation=20)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out / PLOTS[0], dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for i, d in enumerate(("1d", "2d")):
        for j, h in enumerate((2, 3)):
            ax = axes[i, j]
            for response, style in zip(RESPONSES, ("o-", "s--")):
                ax.plot(
                    range(4),
                    [idx[d, 0.3, h, b, response, "observed"]["c95"] for b in MODELS],
                    style,
                    label=response,
                )
            ax.set(title=f"{d} H{h}, alpha=.3", ylabel="conditional observed c95", yscale="log")
            ax.set_xticks(range(4), MODELS, rotation=20)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out / PLOTS[1], dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for i, d in enumerate(("1d", "2d")):
        for j, h in enumerate((2, 3)):
            ax = axes[i, j]
            for b, color in zip(("iso", "temporal_iso", "shrink_full_0.25"), colors):
                for response, marker in zip(RESPONSES, ("o", "s")):
                    for mode, style in zip(MODES, ("-", "--")):
                        ax.plot(
                            ALPHAS,
                            [idx[d, a, h, b, response, mode]["c95"] for a in ALPHAS],
                            color=color,
                            linestyle=style,
                            marker=marker,
                            markersize=3,
                            label=(
                                f"{b} "
                                + ("mean" if response == "complex_mean" else "frozen")
                                + " "
                                + ("obs" if mode == "observed" else "NULL")
                            ),
                        )
            ax.set(title=f"{d} H{h}", xlabel="alpha", ylabel="conditional c95", yscale="log")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(out / PLOTS[2], dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for i, d in enumerate(("1d", "2d")):
        for j, h in enumerate((2, 3)):
            ax = axes[i, j]
            plotted = 0
            for count, marker in ((5, "o"), (3, "s")):
                rr = [
                    r
                    for r in records
                    if (r["dimension_label"], r["contrasts"], r["n_splits"], r["model"])
                    == (d, h, count, "full")
                ]
                valid = [
                    r
                    for r in rr
                    if r["scoreable"]
                    and r["condition_number"] is not None
                    and np.isfinite(r["condition_number"])
                    and r["condition_number"] > 0
                ]
                plotted += len(valid)
                ax.scatter(
                    [r["condition_number"] for r in valid],
                    [r["mean_nll_per_dimension"] for r in valid],
                    marker=marker,
                    label=f"{count}-fold: unavailable {sum(not r['scoreable'] for r in rr)}/{count}",
                )
            ax.set(
                title=f"{d} H{h}: raw full only",
                xscale="log",
                yscale="symlog",
                xlabel="training covariance condition number",
                ylabel="held-out NLL/d (symlog)",
            )
            if not plotted:
                # An entirely unscoreable panel has no data-defined log domain.
                # Set positive display limits, not synthetic points or zero scores.
                ax.set_xlim(1, 10)
                ax.text(0.5, 0.5, "No scoreable folds", ha="center", transform=ax.transAxes)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out / PLOTS[3], dpi=160)
    plt.close(fig)


def table(headers, rows):
    return "\n".join(
        ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
        + ["| " + " | ".join(map(str, r)) + " |" for r in rows]
    )


def report_tables(rows, aggregates):
    idx = index_rows(rows)
    sections = []
    for response in RESPONSES:
        sections.append(
            f"### α=.3 observed：{response}\n\n"
            + table(
                ["背景", "1D H2", "1D H3", "2D H2", "2D H3"],
                [
                    [b]
                    + [
                        f"{idx[d, 0.3, h, b, response, 'observed']['c95']:.7g}"
                        for d in ("1d", "2d")
                        for h in (2, 3)
                    ]
                    for b in MODELS
                ],
            )
        )
    sections.append(
        "### 纯 GRMHD 预测：加权 NLL/d（越低越好；只在同一维度/H/折数内比较）\n\n"
        + table(
            ["维度/H/折数"] + list(MODELS) + ["diagonal*", "raw full*"],
            [
                [f"{d} H{h} / {n}"]
                + [
                    next(
                        (
                            f"{r['weighted_nll']:.6g}"
                            if r["all_folds_scoreable"]
                            else f"不可用 ({r['scoreable_folds']}/{n})"
                        )
                        for r in aggregates
                        if (r["dimension"], r["contrasts"], r["n_splits"], r["model"]) == (d, h, n, b)
                    )
                    for b in MODELS + ("diagonal", "full")
                ]
                for d in ("1d", "2d")
                for h in (2, 3)
                for n in (5, 3)
            ],
        )
    )
    return "\n\n".join(sections)


def run(args):
    output = args.output_root.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError("new noncanonical output root required")
    root = args.inputs.resolve()
    roots = [("main", root)] + [
        (name, getattr(args, name).resolve()) for name in ("embargo", "grid") if getattr(args, name, None)
    ]
    receipts, provenance = [], []
    receipt_root = output.with_name(output.name + "-source-receipts")
    # Verify all completed roots before creating receipts. Never poll or execute them.
    with verification_session():
        for name, source in roots:
            verify_stage(source, recursive=False, allow_imported=False)
        meta = json.loads((root / "metadata.json").read_text())
        if (
            meta["schema"] != "background-reassessment-v1"
            or meta["validate_only"]
            or meta["selected_model"] is not None
            or meta["alphas"] != ALPHAS
            or tuple(meta["inference_models"]) != MODELS
            or meta["posterior_case_count"] != 512
            or meta["prior"] != "log-uniform"
            or meta["coupling_points"] != 120
            or meta["phase_count"] != 64
        ):
            raise ValueError("unexpected completed reassessment metadata")
        rows = read_rows(root / "summary.csv")
        index_rows(rows, ALPHAS)
        records = cv_records(root / "validation.csv")
        summaries = json.loads((root / "validation-summary.json").read_text())
        splits = json.loads((root / "splits.json").read_text())
        aggregates = validate_cv(records, summaries, splits)
        consumed = [
            root / x
            for x in (
                "metadata.json",
                "summary.csv",
                "validation.csv",
                "validation-summary.json",
                "splits.json",
            )
        ]
        consumed += [root / f"{r['key']}.npz" for r in rows]
        require_declared_outputs(root, consumed)
        checks = [validate_posterior(root / f"{r['key']}.npz", r, meta) for r in rows]
        optional = {}
        for name, source in roots[1:]:
            if name == "embargo":
                files = [source / "report.json"]
                require_declared_outputs(source, files)
                data = json.loads(files[0].read_text())
                if data["EHT_used"] or data["hyperparameter_tuning"] or data["embargo_hours"] != 24:
                    raise ValueError("unexpected embargo policy")
                rr = [dict(r, dimension_label=r["dimension"]) for r in data["records"]]
                if len(rr) != 192:
                    raise ValueError("embargo record count mismatch")
                optional[name] = aggregate_cv(rr)
            else:
                files = [source / "report.json", source / "comparison.csv"]
                require_declared_outputs(source, files)
                gm = json.loads(files[0].read_text())
                if (
                    gm.get("alpha") != 0.3
                    or gm.get("no_model_changes") is not True
                    or gm.get("coupling_points") != 240
                    or gm.get("phase_count") != 128
                ):
                    raise ValueError("dense grid policy/grid mismatch")
                coarse = {r["key"]: r for r in rows if r["alpha"] == 0.3}
                gr = read_csv(files[1])
                if len(gr) != 64 or {r["key"] for r in gr} != set(coarse):
                    raise ValueError("dense grid case keys/count mismatch")
                grid_meta = dict(meta, coupling_points=240, phase_count=128)
                optional[name] = []
                for r in gr:
                    old = coarse[r["key"]]
                    if any(
                        str(r[k]) != str(old[k])
                        for k in ("dimension", "contrasts", "background", "response", "mode")
                    ):
                        raise ValueError("dense grid case fields mismatch")
                    if not np.isclose(float(r["coarse_c95"]), old["c95"], rtol=1e-12):
                        raise ValueError("dense grid reference differs from main")
                    dense = float(r["c95"])
                    delta = float(np.log10(dense / old["c95"]))
                    np.testing.assert_allclose(delta, float(r["delta_dex"]), atol=1e-12)
                    file = source / f"{r['key']}.npz"
                    files.append(file)
                    require_declared_outputs(source, [file])
                    # Independent grid stage did not save phase_loglike; do not invent a check.
                    checks.append(
                        dict(
                            source="grid",
                            **validate_posterior(file, dict(old, c95=dense), grid_meta, require_phase=False),
                        )
                    )
                    optional[name].append(
                        dict(key=r["key"], coarse_c95=old["c95"], dense_c95=dense, delta_dex=delta)
                    )
            consumed.extend(files)
        receipt_root.mkdir(parents=True, exist_ok=False)
        for i, (name, source) in enumerate(roots):
            receipt = receipt_root / f"source-receipt-{i:02d}.json"
            receipt.write_bytes(manifest_path(source).read_bytes())
            if sha256(receipt) != sha256(manifest_path(source)):
                raise ValueError("receipt byte mismatch")
            receipts.append(receipt)
            provenance.append(
                dict(
                    name=name,
                    source_manifest_path=str(manifest_path(source)),
                    receipt=receipt.name,
                    sha256=sha256(receipt),
                    recursive=False,
                )
            )
        template = args.report_template.resolve()
        text = template.read_text()
        text += "\n\n## 已校验产物自动附表\n\n" + report_tables(rows, aggregates)
        text += (
            f"\n\n主表 {len(rows)} 项；逐项 posterior 自检 {len(checks)} 项"
            "（若接入 grid，额外 64 项单列）。\n"
        )
        text += "\n".join(f"\n![{p}]({output / p})" for p in PLOTS)
        if "embargo" in optional:
            text += "\n\n### 已完成额外 24 h embargo\n\n" + table(
                ["维度", "H", "折数", "模型", "NLL/d", "相对 ISO Δ", "胜折数（非独立试验）"],
                [
                    [
                        r["dimension"],
                        r["contrasts"],
                        r["n_splits"],
                        r["model"],
                        r["weighted_nll"],
                        r["delta_to_iso"],
                        r["wins_vs_iso"],
                    ]
                    for r in optional["embargo"]
                ],
            )
        if "grid" in optional:
            worst = max(abs(r["delta_dex"]) for r in optional["grid"])
            text += (
                "\n\n已完成 α=.3 dense 240/128 对照：64 项"
                f"最大 |Δlog10 c95|={worst:.8g}；完整 grid-comparison.csv。"
                "不是单轴误差分解。\n"
            )
        config = dict(
            schema="background-reassessment-results-v1",
            canonical=False,
            scientific_recomputation=False,
            selected_model=None,
            row_count=len(rows),
            posterior_check_count=len(checks),
            source_provenance=provenance,
            optional_completed=list(optional),
            recursive_source_verification=False,
        )
        # Only receipt byte files enter manifest inputs, never canonical source manifests.
        with stage_directory(
            output,
            "research:background-reassessment-summary",
            config,
            consumed + receipts + [Path(__file__).resolve(), template],
        ) as out:
            for receipt in receipts:
                shutil.copyfile(receipt, out / receipt.name)
            for name in ("metadata.json", "validation.csv", "validation-summary.json", "splits.json"):
                shutil.copyfile(root / name, out / name)
            write_csv(out / "combined-summary.csv", rows)
            write_csv(out / "cv-comparison.csv", aggregates)
            write_csv(out / "posterior-checks.csv", checks)
            write_csv(out / "mean-response-control.csv", [r for r in rows if r["response"] == "complex_mean"])
            write_csv(
                out / "training-mean-drift.csv",
                [
                    {
                        k: r.get(k)
                        for k in (
                            "dimension_label",
                            "contrasts",
                            "n_splits",
                            "fold",
                            "model",
                            "n_train",
                            "train_mean_rms",
                            "train_mean_norm",
                        )
                    }
                    for r in records
                ],
            )
            for name, data in optional.items():
                write_csv(out / ("grid-comparison.csv" if name == "grid" else "embargo-comparison.csv"), data)
            (out / "provenance.json").write_text(json.dumps(config, indent=2) + "\n")
            (out / "BACKGROUND_REASSESSMENT_ZH.md").write_text(text)
            make_plots(rows, records, aggregates, out)
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", type=Path, default=BASE / "eht-background-reassessment-v2")
    p.add_argument("--output-root", type=Path, default=BASE / "eht-background-reassessment-results-v1")
    p.add_argument("--embargo", type=Path, help="Optional COMPLETED 24h embargo stage only")
    p.add_argument("--grid", type=Path, help="Optional COMPLETED alpha=.3 240/128 stage only")
    p.add_argument(
        "--report-template",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "reports/script-pipeline/BACKGROUND_REASSESSMENT_ZH.md",
    )
    print(run(p.parse_args(argv)), flush=True)


if __name__ == "__main__":
    main()
