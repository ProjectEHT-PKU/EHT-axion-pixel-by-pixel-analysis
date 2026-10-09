"""Explicit, non-destructive naming migration; no scientific recomputation.

Origins are live, checksummed dependencies, not unverified metadata receipts.
A run is built privately on the destination filesystem and published once.
"""

import json
import os
import shutil
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from .common import (
    describe_inputs,
    manifest_path,
    reject_directory_symlinks,
    sha256,
    signature,
    stage_directory,
    verify_stage,
)
from .layout import (
    RunLayout,
    alpha_label,
    date_label,
    input_eht_file,
    input_raptor_files,
)

DIMS = ("1d", "2d")
MODES = ("std", "cor", "cor_grmhd")
DAYS = (95, 96, 100, 101)


def _json(path):
    return json.loads(Path(path).read_text())


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _implementation():
    return {"kind": "migration", "module": __name__, "sha256": sha256(Path(__file__))}


def _trees(source, destination):
    # Check lexical paths BEFORE resolve: an alias must never bypass write safety.
    reject_directory_symlinks(destination)
    source, destination = (
        Path(source).resolve(strict=True),
        Path(destination).absolute(),
    )
    if not source.is_dir():
        raise ValueError("Migration source must be a directory")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Migration destination must be new")
    destination = destination.resolve()
    if (
        source == destination
        or source.is_relative_to(destination)
        or destination.is_relative_to(source)
    ):
        raise ValueError("Source and destination must be separate directory trees")
    return source, destination


@contextmanager
def _atomic_root(destination):
    """No destination is visible until every stage and its lineage are ready."""
    reject_directory_symlinks(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.with_name(destination.name + ".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(
            f"Migration locked: {lock}; inspect owner before clearing stale lock"
        ) from exc
    root = None
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump({"pid": os.getpid(), "started": time.time()}, stream)
        if destination.exists() or destination.is_symlink():
            raise FileExistsError("Migration destination must be new")
        root = Path(
            tempfile.mkdtemp(
                prefix=destination.name + ".building-", dir=destination.parent
            )
        )
        yield root
        reject_directory_symlinks(destination)
        if destination.exists() or destination.is_symlink():
            raise FileExistsError("Migration destination appeared during construction")
        root.rename(destination)
        root = None
    finally:
        try:
            if root is not None:
                shutil.rmtree(root)
        finally:
            lock.unlink(missing_ok=True)


def link(source, target):
    source = Path(source).resolve(strict=True)
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(source)


def stage_files(directory):
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [
        Path(directory) / p["path"] for p in manifest["outputs"]
    ]


def _unchanged(records):
    if describe_inputs([entry["path"] for entry in records]) != records:
        raise ValueError("Migration source changed during construction")


def _source_snapshot(directories):
    """Capture the full verified source graph, including raw and origin inputs."""
    paths, seen = set(), set()

    def visit(directory):
        path = manifest_path(directory).resolve(strict=True)
        if path in seen:
            return
        seen.add(path)
        manifest = verify_stage(directory, recursive=True)
        paths.update(p.resolve(strict=True) for p in stage_files(directory))
        for entry in manifest["inputs"] + manifest.get("origin_inputs", []):
            upstream = Path(entry["path"]).resolve(strict=True)
            paths.add(upstream)
            if upstream.name in {"stage.json", "manifest.json"}:
                record = _json(upstream)
                if record.get("complete") and record.get("schema") == 1:
                    visit(upstream.parent)

    for directory in directories:
        visit(directory)
    return describe_inputs(sorted(paths))


def _rebind(root, destination, stages):
    """Rehash in topological order, reading physical bytes but storing final paths.

    Only exact temporary-root prefixes are mapped; origins remain source paths.
    Rebuild ancestry from upstream manifests rather than editing nested receipts
    blindly (their hashes change when an upstream manifest is rebound).
    """

    def logical(path):
        path = Path(path)
        return (
            str(destination / path.relative_to(root))
            if path.is_relative_to(root)
            else str(path)
        )

    for relative in stages:
        path = root / relative / "manifest.json"
        manifest = _json(path)
        ancestors = []
        for field in ("inputs", "origin_inputs"):
            records = []
            for entry in manifest.get(field, []):
                physical = Path(entry["path"])
                record = describe_inputs([physical])[0]
                record["path"] = logical(physical)
                records.append(record)
                if physical.name in {"stage.json", "manifest.json"}:
                    upstream = _json(physical)
                    if upstream.get("imported_from") or upstream.get(
                        "imported_ancestors"
                    ):
                        ancestors.append(
                            {
                                "manifest": record["path"],
                                "sha256": record["sha256"],
                                "original": upstream.get("imported_from"),
                                "ancestors": upstream.get("imported_ancestors", []),
                            }
                        )
            if field in manifest:
                manifest[field] = records
        payload, key = signature(
            manifest["stage"], manifest["config"], manifest["inputs"]
        )
        manifest.update(payload, signature=key)
        manifest.pop("imported_ancestors", None)
        if ancestors:
            manifest["imported_ancestors"] = ancestors
        for entry in manifest["outputs"]:
            if "external_symlink" in entry:
                entry["external_symlink"] = logical(entry["external_symlink"])
        _write_json(path, manifest)

    # The final paths do not exist yet. Validate every rebound edge using the
    # reverse mapping, without publishing a partially verified run.
    for relative in stages:
        manifest = _json(root / relative / "manifest.json")
        for entry in manifest["inputs"] + manifest.get("origin_inputs", []):
            logical_path = Path(entry["path"])
            physical = (
                root / logical_path.relative_to(destination)
                if logical_path.is_relative_to(destination)
                else logical_path
            )
            actual = describe_inputs([physical])[0]
            if (actual["size"], actual["sha256"]) != (entry["size"], entry["sha256"]):
                raise ValueError(f"Invalid rebound lineage: {logical_path}")


def migrate_inputs(source, destination):
    """Create a canonical input view, preserving all original shared data bytes."""
    source, destination = _trees(source, destination)
    source_manifest = _json(source / "manifest.json")
    entries = source_manifest["files"]
    source_files = {}
    for entry in entries:
        rel = Path(entry["path"])
        if rel.is_absolute() or ".." in rel.parts or str(rel) in source_files:
            raise ValueError("Unsafe/duplicate shared input identity")
        source_files[str(rel)] = entry
    pairs = []
    for alpha in (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44):
        pairs.append(
            (
                input_raptor_files(source, alpha, 100)[0],
                Path("raptor/res-100") / alpha_label(alpha) / "stokes.npy",
            )
        )
    for raw in input_raptor_files(source, 0.3, 400):
        pairs.append((raw, Path("raptor/res-400") / alpha_label(0.3) / raw.name))
    for day in DAYS:
        pairs.append(
            (
                input_eht_file(source, day),
                Path("eht") / date_label(day) / "image-ensemble.pkl",
            )
        )
    for raw in sorted(
        (source / "raptor-transport-diagnostics/snapshot-2110-res-400").glob("*.npy")
    ):
        pairs.append((raw, Path("transport/snapshot-2110-res-400") / raw.name))
    identities = describe_inputs([p for p, _ in pairs])
    for (raw, _), record in zip(pairs, identities):
        expected = source_files.get(str(raw.relative_to(source)))
        if (
            expected is None
            or expected["sha256"] != record["sha256"]
            or expected["bytes"] != record["size"]
        ):
            raise ValueError(f"Shared input checksum mismatch: {raw}")
    snapshot = describe_inputs([source / "manifest.json", *[p for p, _ in pairs]])
    with _atomic_root(destination) as root:
        # stage_directory requires a new leaf; keep its own staging inside root.
        view = root / "view"
        with stage_directory(
            view,
            "inputs:canonical-view",
            {},
            [source / "manifest.json"],
            origin_paths=[p for p, _ in pairs],
            imported_from={
                "source": str(source),
                "purpose": "read-only symlink view",
                "implementation": _implementation(),
            },
        ) as temp:
            for raw, target in pairs:
                link(raw, temp / target)
            _write_json(
                temp / "dataset.json",
                {
                    "layout": "canonical-inputs-v1",
                    "publication_state": "not part of runtime path",
                    "units": "Original Stokes intensity units preserved; not independently verified",
                    "axes": {
                        "raptor": ["time", "axion_state", "x", "y", "stokes"],
                        "stokes": ["I", "Q", "U", "V"],
                    },
                    "source_manifest": str(source / "manifest.json"),
                    "files": [
                        dict(identity, path=str(target))
                        for (_, target), identity in zip(pairs, identities)
                    ],
                },
            )
        verify_stage(view, recursive=True)
        _unchanged(snapshot)
        for product in view.iterdir():
            product.rename(root / product.name)
        view.rmdir()
    return destination


def _response_labels(directory):
    metadata = _json(directory / "fit.json")
    arrays = [
        np.load(directory / f"fits_{dim}.npy", mmap_mode="r", allow_pickle=False)
        for dim in DIMS
    ]
    for dim, array in zip(DIMS, arrays):
        if (
            array.ndim != (3 if dim == "1d" else 4)
            or array.shape[-1] != 3
            or any(n == 0 for n in array.shape)
            or array.dtype.kind != "f"
        ):
            raise ValueError(f"Invalid {dim} response dimensions/dtype: {array.shape}")
    if arrays[0].shape[0] != arrays[1].shape[0]:
        raise ValueError("Response dimensions have inconsistent snapshot counts")
    length = arrays[0].shape[0]
    provided = metadata.get("snapshot_labels")
    if provided is not None:
        labels = np.asarray(provided)
        if (
            labels.shape != (length,)
            or labels.dtype.kind not in "iu"
            or np.any(labels[1:] <= labels[:-1])
        ):
            raise ValueError("Invalid source snapshot_labels dimensions/type/order")
    trusted = (
        provided is not None
        and metadata.get("snapshot_labels_verified") is True
        and metadata.get("snapshot_labels_source") not in (None, "default_convention")
    )
    if trusted:
        return labels, "original_fit_metadata", True
    if provided is not None:
        return labels, "original_fit_metadata_unverified", False
    return 1600 + 10 * np.arange(length), "default_convention", False


def _inference_data(directory, alphas):
    reference, cuts = None, {}
    for dim in DIMS:
        values = []
        for mode in MODES:
            with np.load(directory / f"{mode}_{dim}.npz", allow_pickle=False) as z:
                current = {
                    k: z[k].copy()
                    for k in ("alphas", "cuts", "variant", "observation", "prior")
                }
            identity = {
                k: str(current[k].item()) for k in ("variant", "observation", "prior")
            }
            if identity["observation"] not in {"null", "expected"}:
                raise ValueError("Only zero-day-difference results can be migrated")
            if (
                identity["prior"] not in {"log-uniform", "uniform"}
                or not identity["variant"]
            ):
                raise ValueError("Invalid source inference identity")
            if (
                not np.array_equal(current["alphas"], alphas)
                or current["cuts"].shape != (len(alphas),)
                or not np.isfinite(current["cuts"]).all()
            ):
                raise ValueError("Invalid source inference alpha grid/cuts")
            if reference is not None:
                for key in ("alphas", "variant", "observation", "prior"):
                    if not np.array_equal(reference[key], current[key]):
                        raise ValueError("Mixed source inference identities")
            reference = current
            values.append(current["cuts"])
        cuts[dim] = np.stack(values)
    return reference, cuts


def _preflight(source):
    old_run = _json(source / "07_inference/run.json")
    if old_run.get("options", {}).get("observed", False) is not False:
        raise ValueError("Measured-day-difference runs are not supported")
    results = old_run["results"]
    alphas = np.asarray(next(iter(results.values()))["alphas"])
    if (
        alphas.ndim != 1
        or len(alphas) == 0
        or alphas.dtype.kind not in "fi"
        or not np.isfinite(alphas).all()
        or np.any(np.diff(alphas) <= 0)
    ):
        raise ValueError("Invalid source alpha grid")
    for alpha in alphas:
        alpha_label(alpha)
    for value in results.values():
        if not np.array_equal(value["alphas"], alphas):
            raise ValueError("Mixed source run alpha grids")
    directories = []

    def check(relative, stage, required, expected=None):
        directory = source / relative
        reject_directory_symlinks(directory)
        manifest = verify_stage(directory, recursive=True)
        if not manifest["inputs"]:
            raise ValueError(
                f"Unverifiable source stage (no declared inputs): {directory}"
            )
        if manifest["stage"] != stage:
            raise ValueError(f"Invalid source stage identity: {directory}")
        declared = {entry["path"] for entry in manifest["outputs"]}
        if not set(required).issubset(declared):
            raise ValueError(f"Missing declared source products: {directory}")
        if expected is not None:
            recorded = {Path(e["path"]).resolve() for e in manifest["inputs"]}
            if recorded != {Path(p).resolve() for p in expected}:
                raise ValueError(
                    f"Original lineage belongs to a different run: {directory}"
                )
        directories.append(directory)
        return directory

    labels = {}
    for alpha in alphas:
        processed = check(
            f"01_raptor/alpha-{alpha:g}",
            "raptor:preprocess",
            ["stokes.npy", "preprocess.json"],
        )
        response = check(
            f"02_response/alpha-{alpha:g}",
            "raptor:fit",
            ["fits_1d.npy", "fits_2d.npy", "fit.json"],
            stage_files(processed),
        )
        labels[float(alpha)] = _response_labels(response)
    for dim in DIMS:
        check(
            f"03_templates/{dim}",
            "raptor:summarize",
            ["templates.npz", "summarize.json"],
            [
                source / "02_response" / f"alpha-{a:g}" / f"fits_{dim}.npy"
                for a in alphas
            ],
        )
    required = [
        "mask.npy",
        "selection.npz",
        "report.json",
        "summary_1d.npz",
        "summary_2d.npz",
    ]
    required += [f"{kind}_{day}.npy" for day in DAYS for kind in ("stokes", "evpa_1d")]
    eht = check("04_eht", "eht:ensemble", required)
    for dim in DIMS:
        with np.load(eht / f"summary_{dim}.npz", allow_pickle=False) as z:
            keys = [
                "mean",
                "std",
                "center",
                "unwrapped_mean",
                "concentration",
                "circular_std",
                "count",
            ]
            shapes = [z[key].shape for key in keys]
            if (
                len(set(shapes)) != 1
                or len(shapes[0]) != (2 if dim == "1d" else 3)
                or shapes[0][0] != len(DAYS)
            ):
                raise ValueError("Invalid source EHT summary dimensions")
            if dim == "1d" and z["covariance"].shape != (
                len(DAYS),
                shapes[0][1],
                shapes[0][1],
            ):
                raise ValueError("Invalid source EHT covariance dimensions")
            if dim == "2d" and z["intensity_mean"].shape != shapes[0]:
                raise ValueError("Invalid source EHT intensity dimensions")
    processed = check(
        "01_raptor_background",
        "raptor:background-preprocess",
        ["stokes.npy", "preprocess.json"],
    )
    check(
        "05_background",
        "raptor:background",
        ["background.json", "background_diagnostics.npz"],
        stage_files(processed) + [eht / "mask.npy"],
    )
    for dim in DIMS:
        check(
            f"06_observations/{dim}",
            "observations:covariance",
            ["observations.npz", "observations.json"],
            stage_files(eht) + [source / "03_templates" / dim / "templates.npz"],
        )
    inputs = [
        source / "05_background/stage.json",
        source / "05_background/background.json",
    ]
    for dim in DIMS:
        inputs += [
            source / "06_observations" / dim / "stage.json",
            source / "03_templates" / dim / "stage.json",
            source / "06_observations" / dim / "observations.npz",
            source / "03_templates" / dim / "templates.npz",
        ]
    inference = check(
        "07_inference",
        "inference:six-curves",
        ["run.json"]
        + [f"{m}_{d}.{ext}" for m in MODES for d in DIMS for ext in ("npz", "json")],
        inputs,
    )
    _inference_data(inference, alphas)
    check("08_figures", "figures:limits", [], stage_files(inference))
    return old_run, alphas, labels, _source_snapshot(directories)


def migrate_run(source, destination):
    source, destination = _trees(source, destination)
    old_run, alphas, labels, snapshot = _preflight(source)
    with _atomic_root(destination) as root:
        layout, stages = RunLayout(root), []

        def publish(old, new, inputs, writer):
            manifest = verify_stage(old, recursive=True)
            origin = {
                "path": str(old),
                "manifest_sha256": sha256(manifest_path(old)),
                "original_implementation": manifest["implementation"],
                "original_signature": manifest["signature"],
                "implementation": _implementation(),
                "contract": "live-verified-declared-original-lineage",
            }
            config = json.loads(json.dumps(manifest["config"]))
            if manifest["stage"] == "inference:six-curves" and "observed" in config.get("options", {}):
                if config["options"].pop("observed") is not False:
                    raise ValueError("Measured-day-difference runs are not supported")
            with stage_directory(
                new,
                manifest["stage"],
                config,
                inputs,
                imported_from=origin,
                origin_paths=stage_files(old),
            ) as temp:
                writer(old, temp)
            stages.append(str(new.relative_to(root)))

        def processed_writer(src, dst):
            link(src / "stokes.npy", dst / "stokes.npy")
            shutil.copy2(src / "preprocess.json", dst / "processing.json")

        for alpha in alphas:
            old = source / "01_raptor" / f"alpha-{alpha:g}"
            publish(
                old,
                layout.processed(alpha),
                [Path(e["path"]) for e in verify_stage(old)["inputs"]],
                processed_writer,
            )

            def response_writer(src, dst, alpha=alpha):
                snapshot_labels, label_source, verified = labels[float(alpha)]
                for dim in DIMS:
                    x = np.load(
                        src / f"fits_{dim}.npy", mmap_mode="r", allow_pickle=False
                    )
                    np.savez(
                        dst / f"response-{dim}.npz",
                        amplitude=x[..., 0],
                        phase=x[..., 1],
                        baseline=x[..., 2],
                        snapshot_labels=snapshot_labels,
                    )
                    with np.load(dst / f"response-{dim}.npz", allow_pickle=False) as z:
                        for i, key in enumerate(("amplitude", "phase", "baseline")):
                            np.testing.assert_array_equal(z[key], x[..., i])
                metadata = _json(src / "fit.json")
                metadata.update(
                    snapshot_labels=snapshot_labels.tolist(),
                    snapshot_labels_source=label_source,
                    snapshot_labels_verified=verified,
                    migration_implementation=_implementation(),
                )
                _write_json(dst / "response.json", metadata)

            publish(
                source / "02_response" / f"alpha-{alpha:g}",
                layout.response(alpha),
                stage_files(layout.processed(alpha)),
                response_writer,
            )
        for dim in DIMS:

            def template_writer(src, dst):
                link(src / "templates.npz", dst / "templates.npz")
                shutil.copy2(src / "summarize.json", dst / "summary.json")

            inputs = [
                p
                for a in alphas
                for p in (
                    manifest_path(layout.response(a)),
                    layout.response(a) / f"response-{dim}.npz",
                )
            ]
            publish(
                source / "03_templates" / dim,
                layout.templates(dim),
                inputs,
                template_writer,
            )

        def eht_writer(src, dst):
            for name in ("selection.npz", "mask.npy", "report.json"):
                if name.endswith(".npy"):
                    link(src / name, dst / name)
                else:
                    shutil.copy2(src / name, dst / name)
            with (
                np.load(src / "summary_1d.npz", allow_pickle=False) as one,
                np.load(src / "summary_2d.npz", allow_pickle=False) as two,
            ):
                for k, day in enumerate(DAYS):
                    daydir = dst / date_label(day)
                    daydir.mkdir()
                    link(src / f"stokes_{day}.npy", daydir / "stokes.npy")
                    link(src / f"evpa_1d_{day}.npy", daydir / "evpa-1d.npy")
                    fields = {
                        f"evpa_{key}_{dim}": z[key][k]
                        for dim, z in (("1d", one), ("2d", two))
                        for key in (
                            "mean",
                            "std",
                            "center",
                            "unwrapped_mean",
                            "concentration",
                            "circular_std",
                            "count",
                        )
                    }
                    fields.update(
                        intensity_mean=two["intensity_mean"][k],
                        covariance_1d=one["covariance"][k],
                    )
                    np.savez(daydir / "summary.npz", **fields)

        old_eht = source / "04_eht"
        publish(
            old_eht,
            layout.eht,
            [Path(e["path"]) for e in verify_stage(old_eht)["inputs"]],
            eht_writer,
        )
        old_bg = source / "01_raptor_background"
        publish(
            old_bg,
            layout.background_processed,
            [Path(e["path"]) for e in verify_stage(old_bg)["inputs"]],
            processed_writer,
        )

        def background_writer(src, dst):
            shutil.copy2(src / "background.json", dst / "summary.json")
            shutil.copy2(
                src / "background_diagnostics.npz", dst / "evpa-difference.npz"
            )

        publish(
            source / "05_background",
            layout.background,
            stage_files(layout.background_processed)
            + [manifest_path(layout.eht), layout.eht / "mask.npy"],
            background_writer,
        )
        for dim in DIMS:

            def observation_writer(src, dst):
                shutil.copy2(src / "observations.npz", dst / "observations.npz")
                shutil.copy2(src / "observations.json", dst / "observations.json")

            publish(
                source / "06_observations" / dim,
                layout.observations(dim),
                stage_files(layout.eht)
                + [
                    manifest_path(layout.templates(dim)),
                    layout.templates(dim) / "templates.npz",
                ],
                observation_writer,
            )
        inference_inputs = [
            manifest_path(layout.background),
            layout.background / "summary.json",
        ]
        for dim in DIMS:
            inference_inputs += [
                manifest_path(layout.observations(dim)),
                manifest_path(layout.templates(dim)),
                layout.observations(dim) / "observations.npz",
                layout.templates(dim) / "templates.npz",
            ]

        def inference_writer(src, dst):
            reference, cuts = _inference_data(src, alphas)
            for dim in DIMS:
                for mode in MODES:
                    shutil.copy2(
                        src / f"{mode}_{dim}.json", dst / f"posterior-{mode}-{dim}.json"
                    )
                np.savez(
                    dst / f"limits-{dim}.npz",
                    alphas=reference["alphas"],
                    cuts=cuts[dim],
                    modes=MODES,
                    **{
                        key: reference[key]
                        for key in ("variant", "observation", "prior")
                    },
                )
            shutil.copy2(src / "run.json", dst / "run.json")

        publish(
            source / "07_inference",
            layout.inference,
            inference_inputs,
            inference_writer,
        )

        def figures_writer(src, dst):
            # Copy declared outputs only, never undeclared files beside a stage.
            for entry in verify_stage(src)["outputs"]:
                relative = Path(entry["path"])
                name = (
                    Path("limits-by-dimension.png")
                    if str(relative) == "limits_by_dimension.png"
                    else relative
                )
                target = dst / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src / relative, target)

        publish(
            source / "08_figures",
            layout.figures,
            stage_files(layout.inference),
            figures_writer,
        )
        # Check both converted products and original live ancestry before rebinding.
        for relative in stages:
            verify_stage(root / relative, recursive=True)
        _unchanged(snapshot)
        _rebind(root, destination, stages)
        _write_json(
            root / "migration.json",
            {
                "source": str(source),
                "layout": "canonical-run-v1",
                "scientific_recomputation": False,
                "implementation": _implementation(),
                "origin_contract": "live-verified-declared-original-lineage",
                "provenance_scope": (
                    "Verifies recorded original dependencies and expected legacy stage edges; "
                    "cannot establish provenance omitted from the original manifests"
                ),
                "large_arrays": "symlinked; source files unchanged",
                "converted_response_arrays": "exact equality verified",
                "completed_stages": stages,
                "source_run_metadata": old_run,
            },
        )
        _unchanged(snapshot)
    return destination
