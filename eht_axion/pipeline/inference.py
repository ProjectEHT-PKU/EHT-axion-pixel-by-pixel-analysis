"""Six-curve inference from compact, named pipeline products."""

from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from ..limits import LimitOptions, calculate_from_products
from .constants import PhysicalScale
from .common import stage_directory, verify_stage, sha256, manifest_path
from .layout import RunLayout
from .products import load_templates, load_observations


def _stage_files(directory):
    directory = Path(directory)
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [directory / x["path"] for x in manifest["outputs"]]


def _check_lineage(manifest, expected):
    entries = manifest["inputs"]
    actual = {Path(e["path"]).resolve() for e in entries}
    required = {Path(p).resolve() for p in expected}
    if actual != required:
        raise ValueError("Upstream lineage belongs to a different run/root")
    for entry in entries:
        path = Path(entry["path"])
        if not path.is_file() or sha256(path) != entry["sha256"]:
            raise ValueError(f"Stale upstream lineage: {path}")


def run_inference(output_root, config, *, rebuild=False, reuse_imported=False):
    layout = RunLayout(output_root)
    settings = config["inference"]
    physical = PhysicalScale(**config["physical"])
    background_dir = layout.background
    verify_stage(background_dir, recursive=True, allow_imported=reuse_imported)
    background = json.loads((background_dir / "summary.json").read_text())
    if background["units"] != "rad":
        raise ValueError("Background product must be in radians")
    options = LimitOptions(
        observed=settings["observation"] == "observed",
        prior=settings["prior"],
        phase_count=settings["phase_count"],
        coupling_points=settings["coupling_points"],
        coupling_min=settings["coupling_min"],
        coupling_max=settings["coupling_max"],
        grmhd_sigma_degrees=float(np.rad2deg(background["sigma"])),
        time_conversion=1 / physical.tg_hours,
        exposure_hours=settings["exposure_hours"],
        day_gap_hours=settings["day_gap_hours"],
        pair_separation_hours=settings["pair_separation_hours"],
        n_principal=settings["principal_n"],
        std_covariance="covariance",
    )
    if settings["observation"] not in ("null", "observed"):
        raise ValueError("observation must be null or observed")
    inputs = [manifest_path(background_dir), background_dir / "summary.json"]
    for dim in ("1d", "2d"):
        for stage in [layout.observations(dim), layout.templates(dim)]:
            verify_stage(stage, recursive=True, allow_imported=reuse_imported)
            inputs.append(manifest_path(stage))
        inputs += [
            layout.observations(dim) / "observations.npz",
            layout.templates(dim) / "templates.npz",
        ]
    # Refuse intact-but-stale products or claims about settings not used upstream.
    eht_manifest = verify_stage(layout.eht)
    expected_eht = {
        "topset_size": config["eht"]["topset_size"],
        "weighting": config["aperture"]["weighting"],
        "window_degrees": config["aperture"]["window_degrees"],
        "angular_bins": config["aperture"]["angular_bins"],
        "radial_max": config["aperture"]["radial_max_rg"],
        "fov": config["image"]["final_fov_rg"],
        "inherit_2d_mask": config["eht"]["inherit_2d_mask_in_1d"],
        "min_radius": config["eht"]["min_radius_rg"],
        "intensity_fraction": config["eht"]["intensity_fraction"],
        "historical_mask": config["eht"]["historical_mask"],
    }
    if any(eht_manifest["config"].get(k) != v for k, v in expected_eht.items()):
        raise ValueError("Current config differs from completed EHT stage; rebuild upstream")
    expected_preprocess = dict(
        sigma=config["image"]["gaussian_sigma_rg"],
        init_fov=config["image"]["initial_fov_rg"],
        final_fov=config["image"]["final_fov_rg"],
        final_shape=config["image"]["final_shape"],
        rotation_angle=float(np.deg2rad(config["image"]["rotation_degrees"])),
    )
    expected_fit = dict(
        fov=config["image"]["final_fov_rg"],
        weighting=config["aperture"]["weighting"],
        radial_max=config["aperture"]["radial_max_rg"],
        window_degrees=config["aperture"]["window_degrees"],
        angular_bins=360,
    )
    for alpha in config["data"]["alphas"]:
        processed = verify_stage(layout.processed(alpha))
        if processed["config"] != expected_preprocess:
            raise ValueError("Current image config differs from RAPTOR preprocessing")
        fitted = verify_stage(layout.response(alpha))
        if fitted["config"] != expected_fit:
            raise ValueError("Current aperture config differs from RAPTOR response fits")
        _check_lineage(fitted, _stage_files(layout.processed(alpha)))
    if (
        background["lag"] != config["background"]["lag_frames"]
        or background["method"] != config["background"]["method"]
        or not np.isclose(background["cadence_hours"], 10 * physical.label_step_tg * physical.tg_hours)
    ):
        raise ValueError("Background product differs from current time/model config")
    background_manifest = verify_stage(background_dir)
    background_config = background_manifest["config"]
    if (
        background_config.get("source_alpha") != config["background"]["alpha"]
        or background_config.get("source_resolution") != config["background"]["source_resolution"]
    ):
        raise ValueError("Current background source differs from completed stage")
    bg_processed = verify_stage(layout.background_processed)
    if bg_processed["config"] != dict(expected_preprocess, state_indices=[0]):
        raise ValueError("Current image config differs from background preprocessing")
    _check_lineage(background_manifest, _stage_files(layout.background_processed) + [layout.eht / "mask.npy"])
    for dim in ("1d", "2d"):
        observations_manifest = verify_stage(layout.observations(dim))
        if (
            observations_manifest["config"]["angular_offset"] != config["eht"]["angular_offset"]
            or observations_manifest["config"]["covariance_ridge"] != config["eht"]["covariance_ridge"]
        ):
            raise ValueError("Current config differs from observation selection/covariance")
        _check_lineage(
            observations_manifest,
            _stage_files(layout.eht) + [layout.templates(dim) / "templates.npz"],
        )
        _check_lineage(
            verify_stage(layout.templates(dim)),
            [layout.response(a) / f"response-{dim}.npz" for a in config["data"]["alphas"]],
        )
        template_alphas, _, _, _ = load_templates(layout.templates(dim) / "templates.npz", dim)
        if not np.array_equal(template_alphas, config["data"]["alphas"]):
            raise ValueError("Current alpha grid differs from template products")
    effective = {
        "options": asdict(options),
        "physical": physical.metadata(),
        "background": background,
        "aperture": config["aperture"],
        "eht": config["eht"],
        "alphas": config["data"]["alphas"],
    }
    with stage_directory(
        layout.inference,
        "inference:six-curves",
        effective,
        inputs,
        rebuild=rebuild,
        reuse_imported=reuse_imported,
    ) as temporary:
        if temporary is None:
            return layout.inference
        report = dict(
            effective,
            results={},
            caveats=[
                "Null data are zero-difference/Asimov-like, not median expected limits.",
                "Model errors are diagonal in space and assumed independent between day-pairs.",
                "GRMHD background is scalar at a discrete lag, not exact 24h.",
                "Principal n=4 and image angular calibration need author confirmation.",
                "Posterior limits depend on explicit finite coupling bounds and prior.",
            ],
        )
        for dim in ("1d", "2d"):
            alphas, amplitude, phase, _ = load_templates(layout.templates(dim) / "templates.npz", dim)
            if not np.array_equal(alphas, config["data"]["alphas"]):
                raise ValueError("Current alpha grid differs from template products")
            observations = load_observations(layout.observations(dim) / "observations.npz")
            modes = ("std", "cor", "cor_grmhd")
            dimension_cuts = []
            for mode in modes:
                started = time.perf_counter()
                alpha, cuts, diagnostics = calculate_from_products(
                    observations,
                    amplitude,
                    phase,
                    alphas=alphas,
                    dimension=dim,
                    mode=mode,
                    options=options,
                )
                stem = f"{mode}_{dim}"
                if not np.array_equal(alpha, alphas) or np.shape(cuts) != alphas.shape:
                    raise ValueError(f"Inference alpha grid/curve mismatch: {stem}")
                dimension_cuts.append(cuts)
                # Full numerical diagnostics allow grid/phase/truncation review without reruns.
                (temporary / f"posterior-{mode}-{dim}.json").write_text(
                    json.dumps(diagnostics, indent=2, allow_nan=False) + "\n"
                )
                report["results"][stem] = dict(
                    alphas=alpha.tolist(),
                    cuts=cuts.tolist(),
                    elapsed_seconds=time.perf_counter() - started,
                    selected_points=len(observations.selected_indices),
                )
                print(stem, cuts, flush=True)
            np.savez(
                temporary / f"limits-{dim}.npz",
                alphas=alphas,
                modes=np.asarray(modes),
                cuts=np.stack(dimension_cuts),
                variant="script-pipeline",
                observation=settings["observation"],
                prior=options.prior,
            )
        (temporary / "run.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return layout.inference
