"""Research-specific stage orchestration; no installed command or argument parser."""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np

from .config import read_config
from .common import stage_directory, verify_stage, manifest_path, reject_directory_symlinks
from .layout import RunLayout, input_raptor_files, input_eht_file
from .constants import PhysicalScale


def stage_inputs(directory):
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [Path(directory) / p["path"] for p in manifest["outputs"]]


def run_stage(
    config,
    output_root,
    *,
    input_root=None,
    stage="all",
    trust_pickle=False,
    rebuild=False,
    reuse_imported=False,
):
    """Execute a research stage using explicit paths and verified caches.

    Shared orchestration for the direct scripts; not a general public pipeline API.
    Numerical functions remain independently reusable in their own modules.
    """
    from types import SimpleNamespace

    if stage not in {"all", "raptor", "eht", "background", "observations", "inference", "figures"}:
        raise ValueError(f"Unknown research stage: {stage}")
    args = SimpleNamespace(
        config=Path(config),
        output_root=Path(output_root),
        input_root=Path(input_root) if input_root is not None else None,
        stage=stage,
        trust_pickle=trust_pickle,
        rebuild=rebuild,
        reuse_imported=reuse_imported,
    )
    from .raptor import preprocess_cube, fit_processed, summarize_fits, estimate_background
    from .eht import run_eht_stage
    from .products import build_observations

    config = read_config(args.config)
    scale = PhysicalScale(**config["physical"])
    image, aperture = config["image"], config["aperture"]
    reject_directory_symlinks(args.output_root)
    output = args.output_root.resolve()
    layout = RunLayout(output)
    reuse = getattr(args, "reuse_imported", False)
    validated_lineage = set()

    def stage(*values, **keywords):
        inputs = values[3] if len(values) > 3 else keywords.get("input_paths", ())
        for path in inputs:
            path = Path(path)
            if path.name in {"stage.json", "manifest.json"}:
                record = json.loads(path.read_text())
                if record.get("complete"):
                    verify_stage(path.parent, recursive=True, allow_imported=reuse, _seen=validated_lineage)
        return stage_directory(*values, reuse_imported=reuse, **keywords)

    source = args.input_root.resolve() if args.input_root else None
    if source is not None and (source == output or source.is_relative_to(output)):
        raise ValueError("Output root must not contain the raw input root")
    alphas = config["data"]["alphas"]
    selected = {"all"} if args.stage == "all" else {args.stage}

    def do(stage):
        return "all" in selected or stage in selected

    if args.rebuild:
        targets = {
            "raptor": ["raptor/processed", "raptor/response", "raptor/templates"],
            "eht": ["eht"],
            "background": ["raptor/background", "background"],
            "observations": ["observations"],
            "inference": ["inference"],
            "figures": ["figures"],
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
            inputs = input_raptor_files(source, alpha, 100)
            processed = layout.processed(alpha)
            with stage(
                processed, "raptor:preprocess", process_params, inputs, rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    preprocess_cube(inputs, temporary, **process_params)
            fitted = layout.response(alpha)
            with stage(
                fitted, "raptor:fit", fit_params, stage_inputs(processed), rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    fit_processed(processed, temporary, **fit_params)
        for dim in ("1d", "2d"):
            fits = {a: layout.response(a) / f"response-{dim}.npz" for a in alphas}
            dest = layout.templates(dim)
            with stage(
                dest,
                "raptor:summarize",
                {"dimension": dim, "method": "mle"},
                [p for a in alphas for p in (manifest_path(layout.response(a)), fits[a])],
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
        raw = source
        with stage(
            layout.eht,
            "eht:ensemble",
            params,
            [input_eht_file(raw, d) for d in (95, 96, 100, 101)],
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                run_eht_stage(raw, temporary, **params)
    if do("background"):
        if source is None:
            raise ValueError("--input-root required for background stage")
        cfg = config["background"]
        alpha, res = cfg["alpha"], cfg["source_resolution"]
        inputs = input_raptor_files(source, alpha, res)
        processed = layout.background_processed
        with stage(
            processed,
            "raptor:background-preprocess",
            dict(process_params, state_indices=[0]),
            inputs,
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                preprocess_cube(inputs, temporary, state_indices=[0], **process_params)
        mask = layout.eht / "mask.npy"
        params = dict(
            lag=cfg["lag_frames"],
            cadence_hours=10 * scale.label_step_tg * scale.tg_hours,
            method=cfg["method"],
        )
        with stage(
            layout.background,
            "raptor:background",
            dict(params, physical=scale.metadata(), source_alpha=alpha, source_resolution=res),
            stage_inputs(processed) + [manifest_path(layout.eht), mask],
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                estimate_background(processed / "stokes.npy", temporary, mask=mask, **params)
    if do("observations"):
        for dim in ("1d", "2d"):
            templates = layout.templates(dim) / "templates.npz"
            eht_dir = layout.eht
            params = dict(
                angular_offset=config["eht"]["angular_offset"],
                covariance_ridge=config["eht"]["covariance_ridge"],
            )
            with stage(
                layout.observations(dim),
                "observations:covariance",
                dict(params, dimension=dim),
                stage_inputs(eht_dir) + [manifest_path(layout.templates(dim)), templates],
                rebuild=args.rebuild,
            ) as temporary:
                if temporary is not None:
                    build_observations(eht_dir, templates, dim, temporary / "observations.npz", **params)
    if do("inference"):
        from .inference import run_inference

        run_inference(output, config, rebuild=args.rebuild, reuse_imported=reuse)
    if do("figures"):
        from .figures import plot_results

        inference_dir = layout.inference
        with stage(
            layout.figures,
            "figures:limits",
            {},
            stage_inputs(inference_dir),
            rebuild=args.rebuild,
        ) as temporary:
            if temporary is not None:
                plot_results(inference_dir, temporary, reuse_imported=reuse)
    output.mkdir(parents=True, exist_ok=True)
    snapshot = output / "config.yaml"
    reject_directory_symlinks(snapshot)
    if snapshot.exists() and snapshot.read_bytes() != args.config.read_bytes():
        import hashlib

        digest = hashlib.sha256(args.config.read_bytes()).hexdigest()[:12]
        snapshot = output / "configs" / f"{digest}.yaml"
        reject_directory_symlinks(snapshot)
        snapshot.parent.mkdir(exist_ok=True)
    if snapshot.exists():
        if snapshot.read_bytes() != args.config.read_bytes():
            raise ValueError(f"Config snapshot identity collision: {snapshot}")
    else:
        snapshot.write_bytes(args.config.read_bytes())
    import platform
    from importlib.metadata import version

    root_record = {
        "layout": "canonical-run-v1",
        "requested_config_path": str(snapshot.relative_to(output)),
        "requested_config": config,
        "last_requested_stage": args.stage,
        "reuse_imported_requested": reuse,
        "completion": "See individual stage manifests; this is not a completion marker",
        "source_root": str(source) if source else None,
        "environment": {
            "python": platform.python_version(),
            "numpy": version("numpy"),
            "scipy": version("scipy"),
        },
    }
    import uuid

    invocations = output / "invocations"
    reject_directory_symlinks(invocations)
    invocations.mkdir(exist_ok=True)
    invocation = invocations / f"{uuid.uuid4().hex}.json"
    invocation.write_text(json.dumps(root_record, indent=2) + "\n")
    root_manifest = output / "manifest.json"
    reject_directory_symlinks(root_manifest)
    if not root_manifest.exists():
        root_manifest.write_text(
            json.dumps(
                {
                    "layout": "canonical-run-v1",
                    "initial_config": str(snapshot.relative_to(output)),
                    "source_root": str(source) if source else None,
                    "completion": "See individual stage manifests",
                    "invocations": "Append-only records in invocations/",
                },
                indent=2,
            )
            + "\n"
        )
    return output
