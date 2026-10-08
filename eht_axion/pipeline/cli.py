"""Orchestrate the script-first pipeline with reusable physical stage products."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

try:
    import tomllib
except ImportError:  # Python 3.10
    import tomli as tomllib
from .common import stage_directory, verify_stage
from .constants import PhysicalScale


def read_config(path):
    with Path(path).open("rb") as stream:
        config = tomllib.load(stream)
    # The published angle conventions/offset4 are defined on a one-degree grid.
    if config["aperture"]["angular_bins"] != 360:
        raise ValueError("Research profile requires 360 angular bins (one-degree grid)")
    return config


def stage_inputs(directory):
    manifest = verify_stage(directory)
    return [Path(directory) / "stage.json"] + [Path(directory) / p["path"] for p in manifest["outputs"]]


def run(args):
    from .raptor import preprocess_cube, fit_processed, summarize_fits, estimate_background
    from .eht import run_eht_stage
    from .products import build_observations

    config = read_config(args.config)
    scale = PhysicalScale(**config["physical"])
    image, aperture = config["image"], config["aperture"]
    output = args.output_root.resolve()
    source = args.input_root.resolve() if args.input_root else None
    if source is not None and (source == output or source.is_relative_to(output)):
        raise ValueError("Output root must not contain the raw input root")
    alphas = config["data"]["alphas"]
    selected = {"all"} if args.stage == "all" else {args.stage}

    def do(stage):
        return "all" in selected or stage in selected

    if args.rebuild:
        targets = {
            "raptor": ["01_raptor", "02_response", "03_templates"],
            "eht": ["04_eht"],
            "background": ["01_raptor_background", "05_background"],
            "observations": ["06_observations"],
            "inference": ["07_inference"],
            "figures": ["08_figures"],
        }
        for name, directories in targets.items():
            if do(name) and any((output / d).is_symlink() for d in directories):
                raise ValueError("Refusing to rebuild through a shared stage symlink")

    process_params = dict(
        sigma=image["gaussian_sigma_rg"],
        init_fov=image["initial_fov_rg"],
        final_fov=image["final_fov_rg"],
        final_shape=image["final_shape"],
        rotation_angle=np.deg2rad(image["rotation_degrees"]),
    )
    fit_params = dict(
        fov=image["final_fov_rg"],
        weighting=aperture["weighting"],
        radial_max=aperture["radial_max_rg"],
        window_degrees=aperture["window_degrees"],
        angular_bins=aperture["angular_bins"],
    )
    if do("raptor"):
        if source is None:
            raise ValueError("--input-root required for RAPTOR stages")
        for alpha in alphas:
            inputs = [source / "raptor-stokes" / "res-100" / f"alpha-{alpha}.npy"]
            processed = output / "01_raptor" / f"alpha-{alpha}"
            with stage_directory(
                processed, "raptor:preprocess", process_params, inputs, rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    preprocess_cube(inputs, temporary, **process_params)
            fitted = output / "02_response" / f"alpha-{alpha}"
            with stage_directory(
                fitted, "raptor:fit", fit_params, stage_inputs(processed), rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    fit_processed(processed, temporary, **fit_params)
        for dim in ("1d", "2d"):
            fits = {a: output / "02_response" / f"alpha-{a}" / f"fits_{dim}.npy" for a in alphas}
            dest = output / "03_templates" / dim
            with stage_directory(
                dest,
                "raptor:summarize",
                {"dimension": dim, "method": "mle"},
                list(fits.values()),
                rebuild=args.rebuild,
            ) as temporary:
                if temporary is not None:
                    summarize_fits(fits, temporary)
    if do("eht"):
        if source is None:
            raise ValueError("--input-root required for EHT stage")
        params = dict(
            trust_pickle=args.trust_pickle,
            topset_size=config["eht"]["topset_size"],
            weighting=aperture["weighting"],
            fov=image["final_fov_rg"],
            radial_max=aperture["radial_max_rg"],
            window_degrees=aperture["window_degrees"],
            angular_bins=aperture["angular_bins"],
            inherit_2d_mask=config["eht"]["inherit_2d_mask_in_1d"],
            min_radius=config["eht"]["min_radius_rg"],
            intensity_fraction=config["eht"]["intensity_fraction"],
            historical_mask=config["eht"]["historical_mask"],
        )
        raw = source / "eht-image-ensemble"
        with stage_directory(
            output / "04_eht",
            "eht:ensemble",
            params,
            [raw / f"230_{d}.pkl" for d in (95, 96, 100, 101)],
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                run_eht_stage(raw, temporary, **params)
    if do("background"):
        if source is None:
            raise ValueError("--input-root required for background stage")
        cfg = config["background"]
        alpha, res = cfg["alpha"], cfg["source_resolution"]
        if res == 400:
            inputs = sorted((source / "raptor-stokes" / "res-400" / f"alpha-{alpha}").glob("stokes-*.npy"))
        else:
            inputs = [source / "raptor-stokes" / "res-100" / f"alpha-{alpha}.npy"]
        if not inputs:
            raise ValueError("No background Stokes input found")
        processed = output / "01_raptor_background"
        with stage_directory(
            processed,
            "raptor:background-preprocess",
            dict(process_params, state_indices=[0]),
            inputs,
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                preprocess_cube(inputs, temporary, state_indices=[0], **process_params)
        mask = output / "04_eht" / "mask.npy"
        params = dict(
            lag=cfg["lag_frames"],
            cadence_hours=10 * scale.label_step_tg * scale.tg_hours,
            method=cfg["method"],
        )
        with stage_directory(
            output / "05_background",
            "raptor:background",
            dict(params, physical=scale.metadata(), source_alpha=alpha, source_resolution=res),
            stage_inputs(processed) + [mask],
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                estimate_background(processed / "stokes.npy", temporary, mask=mask, **params)
    if do("observations"):
        for dim in ("1d", "2d"):
            templates = output / "03_templates" / dim / "templates.npz"
            eht_dir = output / "04_eht"
            params = dict(
                angular_offset=config["eht"]["angular_offset"],
                covariance_ridge=config["eht"]["covariance_ridge"],
            )
            with stage_directory(
                output / "06_observations" / dim,
                "observations:covariance",
                dict(params, dimension=dim),
                stage_inputs(eht_dir) + [templates],
                rebuild=args.rebuild,
            ) as temporary:
                if temporary is not None:
                    build_observations(eht_dir, templates, dim, temporary / "observations.npz", **params)
    if do("inference"):
        from .inference import run_inference

        run_inference(output, config, rebuild=args.rebuild)
    if do("figures"):
        from .figures import plot_results

        inference_dir = output / "07_inference"
        with stage_directory(
            output / "08_figures",
            "figures:limits",
            {},
            stage_inputs(inference_dir),
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                plot_results(inference_dir, temporary)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    execute = sub.add_parser("run", help="run one physical stage or the whole chain")
    execute.add_argument("--config", type=Path, required=True)
    execute.add_argument("--input-root", type=Path, help="analysis-inputs/draft (read-only)")
    execute.add_argument("--output-root", type=Path, required=True)
    execute.add_argument(
        "--stage",
        choices=["all", "raptor", "eht", "background", "observations", "inference", "figures"],
        default="all",
    )
    execute.add_argument("--trust-pickle", action="store_true", help="allow trusted EHT pickle inputs")
    execute.add_argument("--rebuild", action="store_true", help="replace stages transactionally")
    verify = sub.add_parser("verify", help="verify a completed stage checksum manifest")
    verify.add_argument("directory", type=Path)
    args = parser.parse_args(argv)
    if args.command == "verify":
        print(json.dumps(verify_stage(args.directory), indent=2))
    else:
        print(f"Pipeline output: {run(args)}")


if __name__ == "__main__":
    main()
