"""Research-specific stage orchestration; no installed command or argument parser."""

from __future__ import annotations
import json
from pathlib import Path

from .config import read_config
from .common import (
    stage_directory, verify_stage, manifest_path, reject_directory_symlinks, verification_session,
)
from .layout import RunLayout, input_raptor_files
from .constants import analysis_geometry


def stage_inputs(directory):
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [Path(directory) / p["path"] for p in manifest["outputs"]]


@verification_session()
def run_stage(
    config,
    output_root,
    *,
    input_root=None,
    raptor_input_root=None,
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

    if stage not in {"all", "raptor"}:
        raise ValueError(f"Unknown research stage: {stage}")
    args = SimpleNamespace(
        config=Path(config),
        output_root=Path(output_root),
        input_root=Path(input_root) if input_root is not None else None,
        raptor_input_root=Path(raptor_input_root) if raptor_input_root is not None else None,
        stage=stage,
        trust_pickle=trust_pickle,
        rebuild=rebuild,
        reuse_imported=reuse_imported,
    )
    from .raptor import preprocess_cube, fit_processed, summarize_fits

    config = read_config(args.config)
    if "final_fov_uas" in config["image"] and not {"cloud", "eht_input"} <= config.keys():
        raise ValueError("Angular analysis requires explicit cloud and eht_input contracts")
    geometry = analysis_geometry(config)
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
    raptor_source = args.raptor_input_root.resolve() if args.raptor_input_root else source
    for raw_root in (source, raptor_source):
        if raw_root is not None and (raw_root == output or raw_root.is_relative_to(output)
                                     or output.is_relative_to(raw_root)):
            raise ValueError("Output root must not overlap a raw input root")
    alphas = config["data"]["alphas"]
    selected = {"all"} if args.stage == "all" else {args.stage}

    def do(stage):
        return "all" in selected or stage in selected

    if args.rebuild:
        targets = {
            "raptor": ["raptor/processed", "raptor/response", "raptor/templates"],
        }
        for name, directories in targets.items():
            if do(name) and any((output / d).is_symlink() for d in directories):
                raise ValueError("Refusing to rebuild through a shared stage symlink")

    process_params = geometry["image_preprocess"]
    fit_params = geometry["aperturefit"]
    if do("raptor"):
        if raptor_source is None:
            raise ValueError("--input-root or --raptor-input-root required for RAPTOR stages")
        from ..producer.delivery import generated_inputs

        for alpha in alphas:
            cubes = input_raptor_files(raptor_source, alpha, 100)
            inputs, labels = generated_inputs(cubes, config, alpha=alpha, resolution=100,
                                             verified=validated_lineage)
            processed = layout.processed(alpha)
            with stage(
                processed, "raptor:preprocess", process_params, inputs, rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    preprocess_cube(cubes, temporary, **process_params)
            fitted = layout.response(alpha)
            effective_fit = dict(fit_params)
            if labels is not None:
                effective_fit["snapshot_labels"] = labels
            with stage(
                fitted, "raptor:fit", effective_fit, stage_inputs(processed), rebuild=args.rebuild
            ) as temporary:
                if temporary is not None:
                    fit_processed(processed, temporary, **effective_fit)
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
        "derived_geometry": geometry,
        "last_requested_stage": args.stage,
        "reuse_imported_requested": reuse,
        "completion": "See individual stage manifests; this is not a completion marker",
        "source_root": str(source) if source else None,
        "raptor_source_root": str(raptor_source) if raptor_source else None,
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
