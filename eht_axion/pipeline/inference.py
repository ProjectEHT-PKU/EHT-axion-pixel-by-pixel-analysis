"""Six-curve inference from compact, named pipeline products."""

from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from ..limits import LimitOptions, calculate_from_products
from .constants import PhysicalScale
from .common import stage_directory, verify_stage, sha256
from .products import load_templates, load_observations


def _stage_files(directory):
    directory = Path(directory)
    manifest = verify_stage(directory)
    return [directory / "stage.json"] + [directory / x["path"] for x in manifest["outputs"]]


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


def run_inference(output_root, config, *, rebuild=False):
    output = Path(output_root)
    settings = config["inference"]
    physical = PhysicalScale(**config["physical"])
    background_dir = output / "05_background"
    verify_stage(background_dir)
    background = json.loads((background_dir / "background.json").read_text())
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
    inputs = [background_dir / "stage.json", background_dir / "background.json"]
    for dim in ("1d", "2d"):
        for stage in [output / "06_observations" / dim, output / "03_templates" / dim]:
            verify_stage(stage)
            inputs.append(stage / "stage.json")
        inputs += [
            output / "06_observations" / dim / "observations.npz",
            output / "03_templates" / dim / "templates.npz",
        ]
    # Refuse intact-but-stale products or claims about settings not used upstream.
    eht_manifest = verify_stage(output / "04_eht")
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
        processed = verify_stage(output / "01_raptor" / f"alpha-{alpha}")
        if processed["config"] != expected_preprocess:
            raise ValueError("Current image config differs from RAPTOR preprocessing")
        fitted = verify_stage(output / "02_response" / f"alpha-{alpha}")
        if fitted["config"] != expected_fit:
            raise ValueError("Current aperture config differs from RAPTOR response fits")
        _check_lineage(fitted, _stage_files(output / "01_raptor" / f"alpha-{alpha}"))
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
    bg_processed = verify_stage(output / "01_raptor_background")
    if bg_processed["config"] != dict(expected_preprocess, state_indices=[0]):
        raise ValueError("Current image config differs from background preprocessing")
    _check_lineage(
        background_manifest, _stage_files(output / "01_raptor_background") + [output / "04_eht" / "mask.npy"]
    )
    for dim in ("1d", "2d"):
        observations_manifest = verify_stage(output / "06_observations" / dim)
        if (
            observations_manifest["config"]["angular_offset"] != config["eht"]["angular_offset"]
            or observations_manifest["config"]["covariance_ridge"] != config["eht"]["covariance_ridge"]
        ):
            raise ValueError("Current config differs from observation selection/covariance")
        _check_lineage(
            observations_manifest,
            _stage_files(output / "04_eht") + [output / "03_templates" / dim / "templates.npz"],
        )
        _check_lineage(
            verify_stage(output / "03_templates" / dim),
            [output / "02_response" / f"alpha-{a}" / f"fits_{dim}.npy" for a in config["data"]["alphas"]],
        )
        template_alphas, _, _, _ = load_templates(output / "03_templates" / dim / "templates.npz", dim)
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
        output / "07_inference", "inference:six-curves", effective, inputs, rebuild=rebuild
    ) as temporary:
        if temporary is None:
            return output / "07_inference"
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
            alphas, amplitude, phase, _ = load_templates(output / "03_templates" / dim / "templates.npz", dim)
            if not np.array_equal(alphas, config["data"]["alphas"]):
                raise ValueError("Current alpha grid differs from template products")
            observations = load_observations(output / "06_observations" / dim / "observations.npz")
            for mode in ("std", "cor", "cor_grmhd"):
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
                np.savez(
                    temporary / f"{stem}.npz",
                    alphas=alpha,
                    cuts=cuts,
                    variant="script-pipeline",
                    observation=settings["observation"],
                    prior=options.prior,
                )
                # Full numerical diagnostics allow grid/phase/truncation review without reruns.
                (temporary / f"{stem}.json").write_text(
                    json.dumps(diagnostics, indent=2, allow_nan=False) + "\n"
                )
                report["results"][stem] = dict(
                    alphas=alpha.tolist(),
                    cuts=cuts.tolist(),
                    elapsed_seconds=time.perf_counter() - started,
                    selected_points=len(observations.selected_indices),
                )
                print(stem, cuts, flush=True)
        (temporary / "run.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return output / "07_inference"
