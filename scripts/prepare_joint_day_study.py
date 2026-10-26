"""Prepare NEW four-day research inputs from verified completed products; no RT.

Verification is intentionally output-only (recursive=False). Exact source
manifest copies are receipts, NOT canonical upstream manifests. The new stage
hashes actual consumed processed products and implementation source modules.
"""

from __future__ import annotations

import argparse
import inspect
import json
import shutil
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.model import cloud_amplitude_scaling, omega
from eht_axion.pipeline.common import (
    manifest_path,
    require_declared_outputs,
    stage_directory,
    verify_stage,
)
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.constants import PhysicalScale, analysis_geometry
from eht_axion.pipeline.layout import alpha_label
from eht_axion.pipeline.raptor import _operator
from eht_axion.research import state_ensemble, time_sampling
from eht_axion.research.state_ensemble import (
    local_branches,
    response_days,
    selected_background,
)
from eht_axion.research.time_sampling import (
    DAY_HOURS,
    H2,
    H3,
    contrast_covariance,
    four_day_stokes,
)


def run(args):
    config = read_config(args.config)
    primary, templates = args.primary.resolve(), args.templates.resolve()
    target = args.output_root.absolute()
    if target.exists() or target.is_symlink():
        raise FileExistsError("new output root required")
    geometry = analysis_geometry(config)
    scale = PhysicalScale(**config["physical"])
    width = config["inference"]["exposure_hours"]
    if (
        width != 6
        or config["inference"]["day_gap_hours"] != 24
        or config["inference"]["pair_separation_hours"] != 120
    ):
        raise ValueError(
            "this study requires six-hour exposures and frozen four-day schedule"
        )
    consumed, receipts, verified = [], [], {}

    def source(directory, names):
        directory = Path(directory)
        if directory not in verified:
            record = verify_stage(directory, recursive=False, allow_imported=False)
            verified[directory] = record
            receipts.append(manifest_path(directory))
        paths = [directory / name for name in names]
        require_declared_outputs(directory, paths)
        consumed.extend(paths)
        return verified[directory]

    source(primary, [f"observations-{dim}.npz" for dim in ("1d", "2d")])
    bg = templates / "raptor/background"
    bg_record = source(bg, ["stokes.npy", "processing.json"])
    if {
        k: v for k, v in bg_record["config"].items() if k != "state_indices"
    } != geometry["image_preprocess"] or bg_record["config"].get("state_indices") != [
        0
    ]:
        raise ValueError("background preprocessing differs from configuration")
    summarydir = templates / "background-unmasked/1d"
    source(summarydir, ["summary.json"])
    summary = json.loads((summarydir / "summary.json").read_text())
    cadence = float(summary["cadence_hours"])
    labels = None
    responses = {dim: [] for dim in ("1d", "2d")}
    for alpha in config["data"]["alphas"]:
        directory = templates / "raptor/response" / alpha_label(alpha)
        record = source(directory, ["response-1d.npz", "response-2d.npz"])
        if {
            k: v for k, v in record["config"].items() if k != "snapshot_labels"
        } != geometry["aperturefit"]:
            raise ValueError("response aperture differs from configuration")
        for dim, dimension_responses in responses.items():
            with np.load(directory / f"response-{dim}.npz", allow_pickle=False) as z:
                current = z["snapshot_labels"]
                if labels is None:
                    labels = current.copy()
                if not np.array_equal(labels, current):
                    raise ValueError("all alpha/dimension snapshot labels must match")
                dimension_responses.append((z["amplitude"], z["phase"]))
    if labels.ndim != 1 or len(labels) < 2 or np.any(np.diff(labels) <= 0):
        raise ValueError("invalid snapshot labels")
    times = (labels - labels[0]) * scale.label_step_tg * scale.tg_hours
    if not np.allclose(np.diff(times), cadence, rtol=1e-12, atol=1e-10):
        raise ValueError(
            "physical label time and processed background cadence disagree"
        )
    processing = json.loads((bg / "processing.json").read_text())
    stokes = np.load(bg / "stokes.npy", mmap_mode="r", allow_pickle=False)
    if stokes.shape != (
        len(labels),
        1,
        config["image"]["final_shape"],
        config["image"]["final_shape"],
        4,
    ):
        raise ValueError("background and response time/grid mismatch")
    # Background native labels are recorded in the processed source shard names;
    # require the inclusive producer schedule and explicit source ordering too.
    snap = config["producer"]["snapshots"]
    expected = np.arange(snap["start"], snap["stop"] + 1, snap["step"])
    if not np.array_equal(labels, expected):
        raise ValueError("response labels differ from physical producer schedule")
    import re

    shardlabels = []
    for path in processing["inputs"]:
        match = re.search(r"stokes-t(\d+)-t(\d+)\.npy$", path)
        if not match:
            raise ValueError("cannot establish background shard labels")
        shardlabels.extend(range(int(match[1]), int(match[2]) + 1, snap["step"]))
    if not np.array_equal(shardlabels, labels):
        raise ValueError("background source shard order differs from response labels")
    starts = times + width / 2
    starts = starts[starts + DAY_HOURS[-1] + width / 2 <= times[-1]]
    if not len(starts):
        raise ValueError("no complete four-day windows")
    observations = {}
    for dim, count in (("1d", 23), ("2d", 20)):
        with np.load(primary / f"observations-{dim}.npz", allow_pickle=False) as z:
            observations[dim] = {k: z[k] for k in z.files}
        obs = observations[dim]
        indices = obs["selected_indices"]
        if (
            indices.shape != (count,)
            or indices.dtype.kind not in "iu"
            or len(np.unique(indices)) != count
        ):
            raise ValueError("frozen final selected feature contract differs")
        if obs["means"].shape != (4, count) or obs["covariances"].shape != (
            4,
            count,
            count,
        ):
            raise ValueError("invalid frozen daily moments")
        lifted, _ = local_branches(obs["means"])
        if not np.allclose(H2 @ lifted, obs["pair_differences"], atol=1e-12):
            raise ValueError(
                "local branch does not preserve primary two-pair observations"
            )
    fit = geometry["aperturefit"]
    operator = _operator(
        config["image"]["final_shape"],
        fit["fov"],
        None,
        fit["radial_max"],
        fit["window_degrees"],
        fit["angular_bins"],
    )
    background = {dim: [] for dim in responses}
    for _, days in four_day_stokes(times, stokes, width, starts):
        for dim, dimension_background in background.items():
            dimension_background.append(
                selected_background(
                    days,
                    observations[dim]["selected_indices"],
                    dim,
                    operator,
                    fit["weighting"],
                )
            )
    inputs = [
        args.config.resolve(),
        Path(__file__).resolve(),
        Path(inspect.getfile(time_sampling)),
        Path(inspect.getfile(state_ensemble)),
        Path(inspect.getfile(omega)),
        Path(inspect.getfile(_operator)),
    ] + consumed
    identity = {
        "schema": "joint-day-inputs-v1",
        "scientific_profile": config,
        "verification": (
            "source stage outputs verified; raw ancestry NOT recursively reverified"
        ),
        "response_units": "radian complex coefficient per unit coupling c",
        "phase": "delta_chi=Re[z exp(i psi)]; psi at first day center",
        "branch": (
            "raw daily angles retained; chronological minimal-gap lift is LOCAL ONLY"
        ),
        "state_weighting": "uniform overlapping model windows, NOT iid measurements",
    }
    with stage_directory(target, "research:joint-day-inputs", identity, inputs) as out:
        for i, path in enumerate(receipts):
            shutil.copyfile(path, out / f"source-receipt-{i:02d}.json")
        report = dict(
            identity,
            cadence_hours=cadence,
            exposure_hours=width,
            snapshot_count=len(labels),
            span_hours=float(times[-1]),
            window_count=len(starts),
            receipt_sources=[str(p) for p in receipts],
            dimensions={},
            warnings=[
                (
                    "Six-hour exposures are underresolved at the native "
                    "8.893-hour cadence."
                ),
                "Integration is exact only for a piecewise-linear interpolant.",
                (
                    "Response fits interpolated as complex envelopes, not refitted "
                    "from averaged phase Stokes."
                ),
                (
                    "Background native resolution 400 vs response 100; final grid 64 "
                    "shared, not identical RT discretization."
                ),
                (
                    "Frozen EHT Monte Carlo daily covariance is conditional and held "
                    "fixed; no calibration claim."
                ),
                (
                    "Minimal-gap local branch does not establish physical winding "
                    "across the 96-hour gap."
                ),
                (
                    "Overlapping windows are correlated model states; no sqrt(N) "
                    "uncertainty reduction."
                ),
            ],
        )
        for dim, dimension_responses in responses.items():
            obs = observations[dim]
            indices = obs["selected_indices"]
            daily = np.asarray(background[dim])
            lifted, branch = local_branches(daily)
            observed_lifted, observed_branch = local_branches(obs["means"])
            frozen, evolving = [], []
            cloud = config["cloud"]
            for alpha, (amp, phase) in zip(
                config["data"]["alphas"], dimension_responses
            ):
                amp, phase = (
                    amp.reshape(len(labels), -1)[:, indices],
                    phase.reshape(len(labels), -1)[:, indices],
                )
                coeff = amp * np.exp(1j * phase)
                coeff *= (
                    cloud_amplitude_scaling(
                        alpha,
                        coefficient=cloud["coefficient"],
                        pivot_alpha=cloud["pivot_alpha"],
                        exponent=cloud["exponent"],
                    )
                    / cloud["reference_norm"]
                )
                f, e = response_days(
                    times,
                    coeff,
                    starts,
                    omega(alpha, config["inference"]["principal_n"]) / scale.tg_hours,
                    width,
                )
                frozen.append(f)
                evolving.append(e)
            contrasts = branch["principal"]
            shutil.copyfile(
                primary / f"observations-{dim}.npz", out / f"observations-{dim}.npz"
            )
            np.savez(
                out / f"states-{dim}.npz",
                alphas=config["data"]["alphas"],
                selected_indices=indices,
                snapshot_labels=labels,
                snapshot_times_hours=times,
                window_centers_hours=starts[:, None] + DAY_HOURS,
                window_starts_hours=starts,
                day_hours=DAY_HOURS,
                H2=H2,
                H3=H3,
                frozen_daily=np.asarray(frozen),
                evolving_daily=np.asarray(evolving),
                frozen_weights=np.full(len(labels), 1 / len(labels)),
                evolving_weights=np.full(len(starts), 1 / len(starts)),
                background_daily=daily,
                background_daily_local_unwrapped=lifted,
                background_h3_raw=branch["raw"],
                background_h3_principal=contrasts,
                background_h3_centered=contrasts - contrasts.mean(axis=0),
                background_h3_mean=contrasts.mean(axis=0),
                background_winding=branch["winding"],
                background_closure_winding=branch["closure_winding"],
                observed_daily_local_unwrapped=observed_lifted,
                observed_h3_principal=observed_branch["principal"],
                observed_winding=observed_branch["winding"],
                observed_closure_winding=observed_branch["closure_winding"],
                covariance_h2=contrast_covariance(obs["covariances"], H2),
                covariance_h3=contrast_covariance(obs["covariances"], H3),
            )
            report["dimensions"][dim] = {
                "features": len(indices),
                "frozen_shape": list(np.shape(frozen)),
                "evolving_shape": list(np.shape(evolving)),
                "background_shape": list(daily.shape),
                "background_nonzero_edge_windings": int(
                    np.count_nonzero(branch["winding"])
                ),
                "background_nonzero_closure_windings": int(
                    np.count_nonzero(branch["closure_winding"])
                ),
                "observed_nonzero_closure_windings": int(
                    np.count_nonzero(observed_branch["closure_winding"])
                ),
                "reference_differences_outside_principal": int(
                    np.count_nonzero(abs(lifted - lifted[:, :1]) > np.pi / 2)
                ),
            }
        (out / "report.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
    return target


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--primary", type=Path, required=True)
    p.add_argument("--templates", type=Path, required=True)
    p.add_argument("--output-root", type=Path, required=True)
    args = p.parse_args(argv)
    print(run(args), flush=True)


if __name__ == "__main__":
    main()
