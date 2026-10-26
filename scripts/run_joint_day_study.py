"""Frozen four-day conditional research study; consumes completed prepared inputs only.

Legacy recovery reads existing declared baseline diagnostics and mean-template
products, retaining the original discrete phase reference and physical scaling.
No RT, raw ancestry verification, or reconstruction.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import shutil
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from research_observed_limits import options_from_config

from eht_axion.model import cloud_amplitude_scaling, evpa_templates
from eht_axion.pipeline.common import (
    manifest_path,
    require_declared_outputs,
    stage_directory,
    verify_stage,
)
from eht_axion.pipeline.products import load_templates
from eht_axion.research import joint_days
from eht_axion.research.joint_days import background_moments, calculate_joint_days

H3 = np.array([[-1.0, 1.0, 0.0, 0.0], [0.0, -1.0, 1.0, 0.0], [0.0, 0.0, -1.0, 1.0]])
H2 = H3[[0, 2]]


def load_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key].copy() for key in z.files}


def alpha_index(grid, alpha):
    matches = np.flatnonzero(np.isclose(grid, alpha, rtol=0, atol=1e-12))
    if len(matches) != 1:
        raise ValueError(f"alpha {alpha} must occur exactly once")
    return int(matches[0])


def load_baseline(baseline, templates, dim, obs):
    """Adapt existing products in memory, never invent a prepared legacy file."""
    report_path = baseline / "report.json"
    diagnostic_path = baseline / f"observed-cor_fixed_background-{dim}.json"
    observation_path = baseline / f"observations-{dim}.npz"
    declared = [report_path, diagnostic_path, observation_path]
    verify_stage(baseline, recursive=False, allow_imported=False)
    require_declared_outputs(baseline, declared)
    report = json.loads(report_path.read_text())
    diagnostic = json.loads(diagnostic_path.read_text())
    original = load_npz(observation_path)
    for key in ("selected_indices", "means", "covariances", "pair_differences"):
        if not np.array_equal(original[key], obs[key]):
            raise ValueError(f"prepared observation differs from frozen primary: {key}")
    options = options_from_config(report["scientific_profile"])
    for key in (
        "phase_count",
        "coupling_points",
        "coupling_min",
        "coupling_max",
        "prior",
        "reference_norm",
        "cloud_coefficient",
        "cloud_pivot_alpha",
        "cloud_exponent",
        "time_conversion",
        "exposure_hours",
        "day_gap_hours",
        "pair_separation_hours",
        "n_principal",
    ):
        if getattr(options, key) != diagnostic["effective_options"][key]:
            raise ValueError(f"primary configuration/diagnostic mismatch: {key}")
    directory = templates / "raptor/templates" / dim
    template_path = directory / "templates.npz"
    verify_stage(directory, recursive=False, allow_imported=False)
    require_declared_outputs(directory, [template_path])
    alphas, amp, phase, _ = load_templates(template_path, dim)
    count = diagnostic["phase_count"]
    if count % 4:
        raise ValueError("legacy phase grid must include pi/2")
    phases = np.arange(count) * 2 * np.pi / count
    coeffs = []
    for alpha in diagnostic["alphas"]:
        alpha_index(alphas, alpha)
        signal, _ = evpa_templates(
            float(alpha),
            amp,
            phase,
            phase_count=count,
            legacy=False,
            duplicate_endpoint=False,
            circular_phase_interpolation=options.circular_phase_interpolation,
            legacy_shared_uncertainty=options.shared_epoch_uncertainty,
            time_conversion=options.time_conversion,
            exposure_hours=options.exposure_hours,
            day_gap_hours=options.day_gap_hours,
            pair_separation_hours=options.pair_separation_hours,
            n_principal=options.n_principal,
            reference_norm=options.reference_norm,
        )
        signal = signal.reshape(count, 2, -1)[..., obs["selected_indices"]]
        signal *= cloud_amplitude_scaling(
            float(alpha),
            coefficient=options.cloud_coefficient,
            pivot_alpha=options.cloud_pivot_alpha,
            exponent=options.cloud_exponent,
        )
        coeff = signal[0] - 1j * signal[count // 4]
        reconstructed = np.real(np.exp(1j * phases[:, None, None]) * coeff)
        if not np.allclose(reconstructed, signal, rtol=1e-11, atol=1e-13):
            raise ValueError("legacy signal is not reproduced on the original phase grid")
        coeffs.append(coeff)
    cuts = report["results"][f"observed-cor_fixed_background-{dim}"]["c95"]
    sigma = report["background_sigma_radians"][dim]
    if not np.isclose(float(sigma) ** 2, diagnostic["background_variance_radians2"], rtol=1e-12):
        raise ValueError("primary pooled background variance differs from diagnostics")
    legacy = {
        "alphas": np.asarray(diagnostic["alphas"]),
        "signal_coeff": np.asarray(coeffs),
        "background_variance": np.asarray(diagnostic["background_variance_radians2"]),
        "recovery_couplings": np.asarray(diagnostic["couplings"]),
        "recovery_phases": phases,
        "recovery_loglike": np.asarray(diagnostic["log_likelihood"]),
        "recovery_c95": np.asarray(cuts),
    }
    return (
        legacy,
        declared + [template_path],
        [manifest_path(baseline), manifest_path(directory)],
    )


def validate_bundle(states, obs, legacy, expected_n=None):
    means = states["observed_daily_local_unwrapped"]
    n = means.shape[1] if means.ndim == 2 else -1
    if means.shape != (4, n) or (expected_n is not None and n != expected_n):
        raise ValueError("invalid selected daily means")
    if not np.array_equal(states["selected_indices"], obs["selected_indices"]):
        raise ValueError("selected feature order differs")
    if not np.array_equal(states["H2"], H2) or not np.array_equal(states["H3"], H3):
        raise ValueError("contrast order must be first, middle, last")
    if not np.allclose(H2 @ means, obs["pair_differences"], rtol=0, atol=1e-12):
        raise ValueError("daily branch changes primary pair observations")
    if not np.allclose(H3 @ means, states["observed_h3_principal"], rtol=0, atol=1e-12):
        raise ValueError("inconsistent H3 angular branch")
    a = len(states["alphas"])
    for name, weights in [
        ("frozen_daily", "frozen_weights"),
        ("evolving_daily", "evolving_weights"),
    ]:
        x = states[name]
        w = states[weights]
        if x.shape != (a, len(w), 4, n) or not np.iscomplexobj(x):
            raise ValueError(f"{name} must be complex (A,M,4,N)")
        if not np.isfinite(x).all() or np.any(w <= 0) or not np.isclose(w.sum(), 1):
            raise ValueError("invalid responses or state weights")
    m = len(states["evolving_weights"])
    if states["background_h3_principal"].shape != (m, 3, n):
        raise ValueError("background/evolving window synchronization mismatch")
    if legacy["signal_coeff"].shape != (len(legacy["alphas"]), 2, n):
        raise ValueError("invalid legacy coefficient dimensions")
    variance = legacy["background_variance"]
    if np.size(variance) != 1 or not np.isfinite(variance).all() or float(variance) < 0:
        raise ValueError("legacy background variance must be scalar nonnegative")
    joint_days.contrast_observation(means, obs["covariances"], H3)


def independent_gaussian_states(b, t, w):
    """Exact Gaussian block moments with B/T cross block zero, no MC pairing.

    Two equally weighted blocks scaled by sqrt(2) retain Cov B and Cov T;
    their cross covariance is identically zero. Only for model='gaussian'.
    """
    bm = np.einsum("m,mkn->kn", w, b)
    tm = np.einsum("m,mkn->kn", w, t)
    bb = np.concatenate([bm + np.sqrt(2) * (b - bm), np.broadcast_to(bm, b.shape)])
    tt = np.concatenate([np.broadcast_to(tm, t.shape), tm + np.sqrt(2) * (t - tm)])
    return bb, tt, np.concatenate([w / 2, w / 2])


def build_cases(states, legacy, alpha, h):
    """Same phenomenological cases for H2/H3; legacy-only controls use H2."""
    i = alpha_index(states["alphas"], alpha)
    rows = [0, 2] if len(h) == 2 else [0, 1, 2]
    b = states["background_h3_principal"][:, rows]
    w = states["evolving_weights"]
    bm, full = background_moments(b, weights=w)
    b = b - bm
    f = np.einsum("kd,mdn->mkn", h, states["frozen_daily"][i])
    e = np.einsum("kd,mdn->mkn", h, states["evolving_daily"][i])
    fw = states["frozen_weights"]
    mean = np.einsum("m,mkn->kn", fw, f)[None]
    # Fix the ISO pool to the SAME original two contrasts for both H bases.
    # Re-pooling with the middle edge would change the model in the H3 comparison.
    _, pool = background_moments(states["background_h3_principal"][:, [0, 2]], weights=w)
    scalar = np.eye(mean[0].size) * np.trace(pool) / pool.shape[0]
    cases = []

    def add(name, **kwargs):
        cases.append((name, kwargs))

    if len(h) == 2:
        old = np.eye(mean[0].size) * float(legacy["background_variance"])
        coeff = legacy["signal_coeff"][alpha_index(legacy["alphas"], alpha)][None]
        add("legacy_recovery", signal_coeff=coeff, background_cov=old)
        add("cartesian_mean_old_background", signal_coeff=mean, background_cov=old)
        add(
            "frozen_mixture_old_background",
            signal_coeff=f,
            weights=fw,
            background_cov=old,
        )
    for mode, cov in [
        ("iso", scalar),
        ("diag", np.diag(np.diag(full))),
        ("full", full),
    ]:
        add("mean_background_" + mode, signal_coeff=mean, background_cov=cov)
    add(
        "mean_background_full_offset",
        signal_coeff=mean,
        background_cov=full,
        background_mean=bm,
    )
    add("frozen_mixture_full", signal_coeff=f, weights=fw, background_cov=full)
    add("evolving_mixture_full", signal_coeff=e, weights=w, background_cov=full)
    add("joint_background_evolving_mixture", signal_coeff=e, backgrounds=b, weights=w)
    add(
        "joint_gaussian_moments",
        signal_coeff=e,
        backgrounds=b,
        weights=w,
        model="gaussian",
    )
    bi, ti, wi = independent_gaussian_states(b, e, w)
    add(
        "independent_gaussian_moments",
        signal_coeff=ti,
        backgrounds=bi,
        weights=wi,
        model="gaussian",
    )
    return cases


def recovery_check(result, legacy, alpha):
    i = alpha_index(legacy["alphas"], alpha)
    if not np.array_equal(result["couplings"], legacy["recovery_couplings"]):
        raise ValueError("recovery coupling grid mismatch")
    if not np.array_equal(result["phases"], legacy["recovery_phases"]):
        raise ValueError("recovery phase grid mismatch")
    expected = legacy["recovery_loglike"][i]
    actual = result["loglike"]
    if not np.allclose(actual - actual.max(), expected - expected.max(), rtol=1e-9, atol=1e-8):
        raise ValueError("legacy posterior shape did not recover")
    if abs(np.log10(result["c95"] / legacy["recovery_c95"][i])) > 1e-10:
        raise ValueError("legacy c95 did not recover")


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    return value


def run(args):
    root = args.inputs.resolve()
    if args.output_root.exists() or args.output_root.is_symlink():
        raise FileExistsError("new output root required")
    verify_stage(root, recursive=False, allow_imported=False)
    paths = [root / "report.json"] + [
        root / f"{kind}-{dim}.npz" for dim in args.dimensions for kind in ("states", "observations")
    ]
    require_declared_outputs(root, paths)
    report = json.loads((root / "report.json").read_text())
    if report["schema"] != "joint-day-inputs-v1":
        raise ValueError("unsupported completed input schema")
    bundles = {}
    receipts = [manifest_path(root)]
    for dim in args.dimensions:
        states, obs = [load_npz(root / f"{kind}-{dim}.npz") for kind in ("states", "observations")]
        legacy, consumed, source_receipts = load_baseline(
            args.baseline_root.resolve(), args.templates_root.resolve(), dim, obs
        )
        paths.extend(consumed)
        receipts.extend(source_receipts)
        validate_bundle(states, obs, legacy, expected_n={"1d": 23, "2d": 20}[dim])
        bundles[dim] = states, obs, legacy
    # Match the original logspace quadrature bit-for-bit at the default grid.
    couplings = np.logspace(np.log10(1 / 137), np.log10(1000), args.coupling_points)
    phases = np.arange(args.phase_count) * 2 * np.pi / args.phase_count
    identity = {
        "schema": "joint-day-study-v1",
        "alphas": args.alphas,
        "coupling_points": args.coupling_points,
        "phase_count": args.phase_count,
        "prior": "log-uniform",
        "source_report": report,
        "warnings": [
            "Conditional frozen EHT covariance, not calibrated coverage.",
            "H3 uses the same temporal state model as H2; gains are not solely new data.",
            "Component minimum chi2 is NOT mixture goodness of fit.",
            "Overlapping windows are model states, not independent observations.",
            "No covariance shrinkage, chi2 tuning, or response rescaling.",
            "ISO background pool is fixed to H2 feature variances for both H2/H3.",
            "Frozen 141-state and evolving 124-window envelopes bracket temporal response assumptions.",
            "Joint B/T mixture vs Gaussian B plus T mixture changes distribution AND cross dependence.",
            "Joint vs independent Gaussian moments isolates B/T cross covariance; marginal moments match.",
        ],
    }
    receipt_root = args.output_root.parent / (args.output_root.name + '-source-receipts')
    receipt_root.mkdir(exist_ok=False)
    receipt_copies = []
    for i, receipt in enumerate(dict.fromkeys(receipts)):
        snapshot = receipt_root / f'source-receipt-{i:02d}.json'
        shutil.copyfile(receipt, snapshot)
        receipt_copies.append(snapshot)
    identity['receipt_sources'] = [str(p) for p in dict.fromkeys(receipts)]
    sources = (
        paths
        + receipt_copies
        + [
            Path(__file__).resolve(),
            Path(inspect.getfile(joint_days)),
            Path(inspect.getfile(options_from_config)),
            Path(inspect.getfile(evpa_templates)),
            Path(inspect.getfile(load_templates)),
        ]
    )
    rows, diagnostics = [], {}
    with stage_directory(
        args.output_root.absolute(), "research:joint-day-study", identity, sources
    ) as out:
        for receipt in receipt_copies:
            shutil.copyfile(receipt, out / receipt.name)
        for dim, (states, obs, legacy) in bundles.items():
            for alpha in args.alphas:
                for h in (H2, H3):
                    for name, kwargs in build_cases(states, legacy, alpha, h):
                        for target in ("observed", "matched_zero"):
                            means = (
                                states["observed_daily_local_unwrapped"]
                                if target == "observed"
                                else np.zeros_like(obs["means"])
                            )
                            result = calculate_joint_days(
                                means,
                                obs["covariances"],
                                h,
                                couplings=couplings,
                                phases=phases,
                                **kwargs,
                            )
                            if name == "legacy_recovery" and target == "observed":
                                reference = result
                                if not (
                                    np.array_equal(couplings, legacy["recovery_couplings"])
                                    and np.array_equal(phases, legacy["recovery_phases"])
                                ):
                                    reference = calculate_joint_days(
                                        means,
                                        obs["covariances"],
                                        h,
                                        couplings=legacy["recovery_couplings"],
                                        phases=legacy["recovery_phases"],
                                        **kwargs,
                                    )
                                recovery_check(reference, legacy, alpha)
                                np.savez(
                                    out / f"recovery-reference-{dim}-a{alpha:.6f}.npz",
                                    **{k: v for k, v in reference.items() if k != "diagnostics"},
                                )
                            key = f"{dim}-a{alpha:.6f}-h{len(h)}-{name}-{target}"
                            np.savez(
                                out / f"{key}.npz",
                                **{k: v for k, v in result.items() if k != "diagnostics"},
                            )
                            d = result["diagnostics"]
                            diagnostics[key] = jsonable(
                                dict(
                                    d,
                                    observed_edge_winding=states.get("observed_winding"),
                                    observed_closure_winding=states.get("observed_closure_winding"),
                                    background_edge_winding=states.get("background_winding"),
                                    background_closure_winding=states.get(
                                        "background_closure_winding"
                                    ),
                                )
                            )
                            if "backgrounds" in kwargs:
                                bmean, cov = background_moments(
                                    kwargs["backgrounds"], weights=kwargs["weights"]
                                )
                            else:
                                cov = kwargs.get(
                                    "background_cov",
                                    np.zeros_like(result["covariance"]),
                                )
                                bmean = kwargs.get("background_mean", np.zeros_like(result["data"]))
                            eig = np.linalg.eigvalsh(cov)
                            fixed_eig = np.linalg.eigvalsh(result["covariance"])
                            coeff = kwargs["signal_coeff"]
                            weights = kwargs.get("weights", np.full(len(coeff), 1 / len(coeff)))
                            tm = np.einsum("m,mkn->kn", weights, coeff)
                            tv = float(
                                np.einsum(
                                    "m,mkn,mkn->",
                                    weights,
                                    (coeff - tm).conj(),
                                    coeff - tm,
                                ).real
                            )
                            cross_norm = 0.0
                            if "backgrounds" in kwargs:
                                bc = kwargs["backgrounds"].reshape(len(coeff), -1)
                                tc = np.concatenate(
                                    [
                                        coeff.real.reshape(len(coeff), -1),
                                        coeff.imag.reshape(len(coeff), -1),
                                    ],
                                    axis=1,
                                )
                                cross = ((bc - weights @ bc).T * weights) @ (tc - weights @ tc)
                                cross_norm = float(np.linalg.norm(cross))
                            rows.append(
                                {
                                    "key": key,
                                    "dimension": dim,
                                    "alpha": alpha,
                                    "contrasts": len(h),
                                    "mode": target,
                                    "case": name,
                                    "phase_count": len(phases),
                                    "n": means.shape[1],
                                    "state_count": len(coeff),
                                    "state_weights": "explicit whole-state",
                                    "c95": result["c95"],
                                    "cbg_rank": int(np.linalg.matrix_rank(cov)),
                                    "cbg_trace": float(np.trace(cov)),
                                    "cbg_eigen_min": float(eig[0]),
                                    "background_mean_norm": float(np.linalg.norm(bmean)),
                                    "fixed_cov_condition": d["condition_number"],
                                    "fixed_cov_eigen_min": float(fixed_eig[0]),
                                    "background_cov_role": "latent_moment_not_added"
                                    if "backgrounds" in kwargs
                                    else "added_fixed_gaussian",
                                    "phase_likelihood_model": kwargs.get("model", "empirical"),
                                    "response_variance_trace": tv,
                                    "mean_response_norm": float(np.linalg.norm(tm)),
                                    "bt_cross_cov_frobenius": cross_norm,
                                    "response_variance_over_mean_power": tv
                                    / max(
                                        float(np.linalg.norm(tm) ** 2),
                                        np.finfo(float).tiny,
                                    ),
                                    "lower_bin_mass": d["edge_interval_mass"][0],
                                    "upper_bin_mass": d["edge_interval_mass"][1],
                                    "likelihood_peak_at_boundary": bool(
                                        np.argmax(result["loglike"]) in (0, len(couplings) - 1)
                                    ),
                                    "minimum_component_chi2": d["minimum_component_chi2"],
                                }
                            )
        with (out / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (out / "metadata.json").write_text(
            json.dumps(
                jsonable(dict(identity, diagnostics=diagnostics)),
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        # Fixed constructor order, never ranked by constraint strength.
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for dim in args.dimensions:
            for alpha in args.alphas:
                fig, ax = plt.subplots(figsize=(11, 6))
                subset = [r for r in rows if r["dimension"] == dim and r["alpha"] == alpha]
                names = list(dict.fromkeys(r["case"] for r in subset))
                for k, target, marker in [
                    (2, "observed", "o"),
                    (3, "observed", "s"),
                    (2, "matched_zero", "x"),
                    (3, "matched_zero", "+"),
                ]:
                    rr = [r for r in subset if r["contrasts"] == k and r["mode"] == target]
                    ax.plot(
                        [names.index(r["case"]) for r in rr],
                        [r["c95"] for r in rr],
                        marker,
                        label=f"H{k} {target}",
                    )
                ax.set(
                    yscale="log",
                    ylabel="conditional c95",
                    title=f"{dim}, alpha={alpha:g}",
                )
                ax.set_xticks(range(len(names)), names, rotation=65, ha="right")
                ax.legend()
                fig.tight_layout()
                fig.savefig(out / f"case-controls-{dim}-a{alpha:.6f}.png", dpi=150)
                plt.close(fig)
    return args.output_root


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", type=Path, required=True)
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument(
        "--baseline-root",
        type=Path,
        default=Path("/home/canoziia/GRMHD_DATA/analysis/eht-120uas-observed-limits-v1"),
    )
    p.add_argument(
        "--templates-root",
        type=Path,
        default=Path("/home/canoziia/GRMHD_DATA/analysis/eht-120uas-observed-templates-v1"),
    )
    p.add_argument("--alphas", type=float, nargs="+", default=[0.3])
    p.add_argument("--dimensions", nargs="+", choices=["1d", "2d"], default=["1d", "2d"])
    p.add_argument("--coupling-points", type=int, default=120)
    p.add_argument("--phase-count", type=int, default=64)
    args = p.parse_args(argv)
    if (
        args.coupling_points < 2
        or args.phase_count < 4
        or not args.alphas
        or not np.isfinite(args.alphas).all()
    ):
        p.error("invalid grid")
    print(run(args), flush=True)


if __name__ == "__main__":
    main()
