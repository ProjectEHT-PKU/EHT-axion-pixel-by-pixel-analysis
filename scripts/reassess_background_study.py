"""Background-only purged validation and limited conditional posterior reassessment.

Consumes completed prepared products, never raw RT/EHT data. CV is descriptive,
not a rule for selecting the strongest limit. Outputs are noncanonical research
products. No evolving or paired joint-background/response models are fitted.
"""

from __future__ import annotations

# Environment setup must precede numerical imports.
# ruff: noqa: E402

import argparse
import csv
import inspect
import json
import os
import shutil
from pathlib import Path

# Script-local BLAS control, applied before importing numerical libraries.
for _variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_variable] = "2"

import _bootstrap  # noqa: E402,F401
import numpy as np  # noqa: E402

from eht_axion.pipeline.common import (
    describe_inputs,
    manifest_path,
    require_declared_outputs,
    stage_directory,
    verification_session,
    verify_stage,
)
from eht_axion.research import background_validation as validation
from eht_axion.research import joint_days
from run_joint_day_study import H2, H3, alpha_index, load_npz

# Fixed before any data are scored; coefficient is not tuned by this runner.
INFERENCE_MODELS = ("iso", "temporal_iso", "partial_diag", "shrink_full_0.25")
ALL_MODELS = INFERENCE_MODELS + ("diagonal", "full")
RESPONSES = ("complex_mean", "frozen_mixture")
DAY_HOURS = np.array([0.0, 24.0, 120.0, 144.0])


def jsonable(value):
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None  # undefined condition numbers stay undefined, never repaired
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def write_json(path, value):
    path.write_text(json.dumps(jsonable(value), indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                k: json.dumps(jsonable(v))
                if isinstance(v, (dict, list, np.ndarray))
                else jsonable(v)
                for k, v in row.items()
            }
            for row in rows
        )


def validate_inputs(states, obs):
    """Check frozen branches, support/time units, whole-state weights and sizes."""
    b = states["background_h3_principal"]
    means = states["observed_daily_local_unwrapped"]
    if b.ndim != 3 or b.shape[1] != 3 or means.shape != (4, b.shape[2]):
        raise ValueError("background must be (M,3,N), daily means (4,N)")
    for x in (b, means, obs["covariances"], states["frozen_daily"]):
        if not np.isfinite(x).all():
            raise ValueError("nonfinite prepared data")
    if not np.array_equal(states["selected_indices"], obs["selected_indices"]):
        raise ValueError("selected feature ordering differs")
    if not np.array_equal(states["H2"], H2) or not np.array_equal(states["H3"], H3):
        raise ValueError("contrast definitions differ")
    if not np.allclose(H2 @ means, obs["pair_differences"], rtol=0, atol=1e-12):
        raise ValueError("prepared pair branch differs")
    if not np.allclose(H3 @ means, states["observed_h3_principal"], rtol=0, atol=1e-12):
        raise ValueError("prepared H3 branch differs")
    if not np.array_equal(states["day_hours"], DAY_HOURS):
        raise ValueError("expected day times 0,24,120,144 hours")
    starts = states["window_starts_hours"]
    centers = states["window_centers_hours"]
    if len(starts) != len(b) or not np.allclose(
        centers, starts[:, None] + DAY_HOURS, rtol=0, atol=1e-12
    ):
        raise ValueError("window synchronization differs")
    validation.support_indices(
        starts, states["snapshot_times_hours"], DAY_HOURS, 6.0, window_centers=centers
    )
    f, w = states["frozen_daily"], states["frozen_weights"]
    if (
        f.shape != (len(states["alphas"]), len(w), 4, b.shape[2])
        or not np.iscomplexobj(f)
        or np.any(w <= 0)
        or not np.isfinite(w).all()
        or not np.allclose(w, 1 / len(w), rtol=0, atol=1e-14)
    ):
        raise ValueError("frozen response must have uniform whole-state weights")
    joint_days.contrast_observation(means, obs["covariances"], H3)


def check_covariance(cov, meta, name, n_train, dimension):
    """Only dimensionally undefined raw-full scoring may be skipped."""
    if not np.isfinite(cov).all() or np.trace(cov) <= 0:
        raise ValueError(f"{name}: zero background sigma or nonfinite covariance")
    if name == "full" and n_train <= dimension:
        return False
    if not meta["scoreable"]:
        raise ValueError(f"{name}: unexpected non-SPD covariance")
    np.linalg.cholesky(cov)
    return True


def evaluate_background(states, contrast_rows, n_splits):
    result = validation.evaluate_folds(
        states["background_h3_principal"],
        states["window_starts_hours"],
        states["snapshot_times_hours"],
        contrast_rows,
        n_splits=n_splits,
        day_hours=DAY_HOURS,
        width=6.0,
        window_centers=states["window_centers_hours"],
    )
    d = len(contrast_rows) * states["background_h3_principal"].shape[2]
    for record in result["records"]:
        if record["n_train"] < 2:
            continue  # undefined, retained in helper summary, not discarded
        if record["model"] == "full" and record["n_train"] <= d:
            record.update(scoreable=False, failure="n_train_le_dimension_raw_full")
            for key in (
                "nll_per_dimension",
                "mean_nll_per_dimension",
                "mean_quadratic_per_dimension",
            ):
                record.pop(key, None)
        elif not record["scoreable"]:
            raise ValueError(f"unexpected non-SPD validation fold: {record['model']}")
    # Guard even if a helper happened to score a numerically nonsingular raw full.
    full = [r for r in result["records"] if r["model"] == "full"]
    if not all(r["scoreable"] for r in full):
        result["summary"]["full"] = dict(
            all_folds_scoreable=False,
            scoreable_folds=sum(r["scoreable"] for r in full),
            diagnostic_only=True,
        )
    return result


def fit_backgrounds(states, rows):
    b = states["background_h3_principal"][:, rows, :]
    fits = {}
    for name in ALL_MODELS:
        cov, meta = validation.model_covariance(b, name)
        valid = check_covariance(cov, meta, name, len(b), b.shape[1] * b.shape[2])
        if not valid:
            meta = dict(meta, scoreable=False, failure="n_train_le_dimension_raw_full")
        fits[name] = dict(covariance=cov, metadata=meta, scoreable=valid)
    return fits


def build_cases(states, alpha, h, fits):
    """141 whole states, Cartesian mean or uniform mixture; no B/T pairing."""
    ai = alpha_index(states["alphas"], alpha)
    f = np.einsum("kd,mdn->mkn", h, states["frozen_daily"][ai])
    w = states["frozen_weights"]
    mean = np.einsum("m,mkn->kn", w, f)[None]
    for name in INFERENCE_MODELS:
        for response in RESPONSES:
            kwargs = dict(
                signal_coeff=mean if response == "complex_mean" else f,
                background_cov=fits[name]["covariance"],
            )
            if response == "frozen_mixture":
                kwargs["weights"] = w
            yield name, response, kwargs


def compare_old(rows, path):
    """H2 ISO only: H3 now pools all three contrasts and is not a recovery test."""
    with path.open(newline="") as stream:
        old = list(csv.DictReader(stream))
    comparisons = []
    for row in rows:
        if (row["contrasts"], row["background"], row["response"]) != (
            2,
            "iso",
            "complex_mean",
        ):
            continue
        matches = [
            r
            for r in old
            if r["dimension"] == row["dimension"]
            and abs(float(r["alpha"]) - row["alpha"]) < 1e-12
            and int(r["contrasts"]) == 2
            and r["mode"] == row["mode"]
            and r["case"] == "mean_background_iso"
        ]
        if len(matches) != 1:
            raise ValueError("old H2 ISO comparison must match exactly once")
        delta = float(np.log10(row["c95"] / float(matches[0]["c95"])))
        if abs(delta) > 1e-9:
            raise ValueError("H2 ISO c95 failed old-study recovery (check quadrature)")
        comparisons.append(dict(key=row["key"], log10_c95_ratio=delta))
    return comparisons


def run(args):
    root = args.inputs.resolve()
    if args.output_root.exists() or args.output_root.is_symlink():
        raise FileExistsError("new noncanonical output root required; no append/resume")
    verify_stage(root, recursive=False, allow_imported=False)
    paths = [root / "report.json"] + [
        root / f"{kind}-{dim}.npz"
        for dim in args.dimensions
        for kind in ("states", "observations")
    ]
    require_declared_outputs(root, paths)
    report = json.loads(paths[0].read_text())
    if report["schema"] != "joint-day-inputs-v1":
        raise ValueError("unsupported completed input schema")
    bundles = {}
    for dim in args.dimensions:
        s, o = (
            load_npz(root / f"{kind}-{dim}.npz") for kind in ("states", "observations")
        )
        validate_inputs(s, o)
        if s["background_h3_principal"].shape != (124, 3, {"1d": 23, "2d": 20}[dim]):
            raise ValueError("expected completed 124-window background bundle")
        if len(s["frozen_weights"]) != 141:
            raise ValueError("expected completed 141-state response bundle")
        bundles[dim] = s, o
    alphas = list(
        bundles[args.dimensions[0]][0]["alphas"] if args.alphas is None else args.alphas
    )
    if not alphas or len(set(alphas)) != len(alphas) or not np.isfinite(alphas).all():
        raise ValueError("alphas must be finite, nonempty and unique")
    for s, _ in bundles.values():
        for alpha in alphas:
            alpha_index(s["alphas"], alpha)
    receipts = [manifest_path(root)]
    sources = paths + [
        Path(__file__).resolve(),
        Path(inspect.getfile(validation)),
        Path(inspect.getfile(joint_days)),
        Path(inspect.getfile(load_npz)),
    ]
    if args.old_summary:
        old_root = args.old_summary.resolve().parent
        verify_stage(old_root, recursive=False, allow_imported=False)
        require_declared_outputs(old_root, [args.old_summary.resolve()])
        sources.append(args.old_summary.resolve())
        receipts.append(manifest_path(old_root))
    receipt_root = args.output_root.parent / (args.output_root.name + '-source-receipts')
    receipt_root.mkdir(exist_ok=False)
    for i, receipt in enumerate(receipts):
        snapshot = receipt_root / f'source-receipt-{i:02d}.json'
        snapshot.write_bytes(receipt.read_bytes())
        sources.append(snapshot)
    identity = dict(
        schema="background-reassessment-v1",
        canonical=False,
        alphas=alphas,
        dimensions=args.dimensions,
        inference_models=INFERENCE_MODELS,
        validation_models=ALL_MODELS,
        selected_model=None,
        validate_only=args.validate_only,
        phase_count=64,
        coupling_points=120,
        prior="log-uniform",
        prior_bounds=[1 / 137, 1000],
        units=dict(
            background="radians",
            covariance="radians^2",
            times="hours",
            coupling="same physical unit-coupling convention as prepared frozen_daily",
        ),
        sample_counts=dict(
            background_windows=124,
            frozen_response_states=141,
            counts_are_not_iid_effective_sample_sizes=True,
        ),
        validation_comparisons="within dimension/H only; no H2/H3 score ranking",
        policy=[
            "CV predicts zero mean using training-centered MLE covariance (ddof=0).",
            "Purging removes shared raw support, not all dependence; no iid errors.",
            "No EHT chi2, posterior strength, or c95 used to select/tune candidates.",
            "All four predefined stable models are inferred regardless of CV scores.",
            "H2 ISO uses its natural two-edge pool; H3 ISO repools all three edges.",
            "H3 changes data and background target: not pure information isolation.",
            "Old full and joint-empirical models: not preferred with sparse support.",
            "141 states define a support prior, not 141 independent observations.",
            "Frozen responses only; evolving and paired B/T temporarily omitted.",
            "Gaussian predictive scatter: no extra sqrt(N) or fitted offset.",
            "Conditional observed/matched-zero limits are not calibrated coverage.",
        ],
        source_report=report,
        receipt_sources=[str(p) for p in receipts],
    )
    with stage_directory(
        args.output_root.absolute(),
        "research:background-reassessment",
        jsonable(identity),
        sources,
    ) as out:
        shutil.copyfile(manifest_path(root), out / "prepared-source-receipt.json")
        if args.old_summary:
            shutil.copyfile(
                manifest_path(args.old_summary.resolve().parent),
                out / "old-study-source-receipt.json",
            )
        write_json(out / "input-sha256.json", describe_inputs(sources))
        records, summaries, splits, fits_by_dim = [], {}, {}, {}
        for dim, (s, _) in bundles.items():
            fits_by_dim[dim] = {}
            for k, rr in ((2, [0, 2]), (3, [0, 1, 2])):
                fits = fit_backgrounds(s, rr)
                fits_by_dim[dim][k] = fits
                summaries[f"{dim}-h{k}-full_training"] = fits
                for count in (5, 3):
                    result = evaluate_background(s, rr, count)
                    key = f"{dim}-h{k}-folds{count}"
                    splits[key] = result["folds"]
                    summaries[key] = {
                        key: value
                        for key, value in result.items()
                        if key not in ("folds", "records")
                    }
                    # Preserve per-window held-out scores and helper metadata.
                    write_json(out / f"calibration-{key}.json", result)
                    records.extend(
                        dict(r, dimension_label=dim, contrasts=k, n_splits=count)
                        for r in result["records"]
                    )
        write_csv(out / "validation.csv", records)
        write_json(out / "splits.json", splits)
        write_json(out / "validation-summary.json", summaries)
        rows = []
        if not args.validate_only:
            couplings = np.logspace(np.log10(1 / 137), np.log10(1000), 120)
            phases = np.arange(64) * 2 * np.pi / 64
            for dim, (s, obs) in bundles.items():
                for alpha in alphas:
                    for h in (H2, H3):
                        for bg, response, kwargs in build_cases(
                            s, alpha, h, fits_by_dim[dim][len(h)]
                        ):
                            for target in ("observed", "matched_zero"):
                                means = (
                                    s["observed_daily_local_unwrapped"]
                                    if target == "observed"
                                    else np.zeros_like(obs["means"])
                                )
                                result = joint_days.calculate_joint_days(
                                    means,
                                    obs["covariances"],
                                    h,
                                    couplings=couplings,
                                    phases=phases,
                                    **kwargs,
                                )
                                key = (
                                    f"{dim}-a{alpha:.6f}-h{len(h)}-"
                                    f"{bg}-{response}-{target}"
                                )
                                np.savez(
                                    out / f"{key}.npz",
                                    **{
                                        k: v
                                        for k, v in result.items()
                                        if k != "diagnostics"
                                    },
                                )
                                write_json(out / f"{key}.json", result["diagnostics"])
                                rows.append(
                                    dict(
                                        key=key,
                                        dimension=dim,
                                        alpha=float(alpha),
                                        contrasts=len(h),
                                        background=bg,
                                        response=response,
                                        mode=target,
                                        c95=result["c95"],
                                        state_count=len(kwargs["signal_coeff"]),
                                        fixed_cov_condition=result["diagnostics"][
                                            "condition_number"
                                        ],
                                    )
                                )
            write_csv(out / "summary.csv", rows)
        identity["legacy_h2_iso_recovery"] = (
            compare_old(rows, args.old_summary) if args.old_summary and rows else None
        )
        identity["posterior_case_count"] = len(rows)
        write_json(out / "metadata.json", identity)
    return args.output_root


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    analysis = Path("/home/canoziia/GRMHD_DATA/analysis")
    p.add_argument("--inputs", type=Path, default=analysis / "eht-joint-days-inputs-v1")
    p.add_argument(
        "--output-root", type=Path, default=analysis / "eht-background-reassessment-v1"
    )
    p.add_argument("--alphas", nargs="+", type=float, default=None)
    p.add_argument(
        "--dimensions", nargs="+", choices=("1d", "2d"), default=["1d", "2d"]
    )
    p.add_argument("--validate-only", action="store_true")
    p.add_argument(
        "--old-summary",
        type=Path,
        help="optional completed old study summary.csv; strict H2 ISO recovery",
    )
    args = p.parse_args(argv)
    if len(set(args.dimensions)) != len(args.dimensions):
        p.error("duplicate dimensions")
    with verification_session():
        print(run(args), flush=True)


if __name__ == "__main__":
    main()
