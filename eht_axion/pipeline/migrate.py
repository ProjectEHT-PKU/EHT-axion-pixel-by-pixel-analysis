"""Explicit, non-destructive naming migration; no scientific recomputation."""

from pathlib import Path
import json
import shutil
import numpy as np
from .common import stage_directory, verify_stage, manifest_path, describe_inputs
from .layout import RunLayout, date_label, alpha_label, input_raptor_files, input_eht_file


def link(source, target):
    source = Path(source).resolve(strict=True)
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.symlink_to(source)


def stage_files(directory):
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [Path(directory) / p["path"] for p in manifest["outputs"]]


def migrate_inputs(source, destination):
    """Create a canonical input view, preserving all original shared data bytes."""
    source, destination = Path(source).resolve(strict=True), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Input-view destination must be new")
    destination = destination.resolve()
    if source == destination or source.is_relative_to(destination) or destination.is_relative_to(source):
        raise ValueError("Source and destination must be separate directory trees")
    source_manifest = json.loads((source / "manifest.json").read_text())
    source_files = {entry["path"]: entry for entry in source_manifest["files"]}
    pairs = []
    for alpha in (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44):
        raw = input_raptor_files(source, alpha, 100)[0]
        pairs.append((raw, Path("raptor/res-100") / alpha_label(alpha) / "stokes.npy"))
    for raw in input_raptor_files(source, 0.3, 400):
        pairs.append((raw, Path("raptor/res-400") / alpha_label(0.3) / raw.name))
    for day in (95, 96, 100, 101):
        pairs.append((input_eht_file(source, day), Path("eht") / date_label(day) / "image-ensemble.pkl"))
    old_transport = source / "raptor-transport-diagnostics/snapshot-2110-res-400"
    for raw in sorted(old_transport.glob("*.npy")):
        pairs.append((raw, Path("transport/snapshot-2110-res-400") / raw.name))
    # Verify the source manifest BEFORE creating a view of potentially invalid data.
    identities = describe_inputs([p for p, _ in pairs])
    for (raw, _), record in zip(pairs, identities):
        rel = str(raw.relative_to(source))
        expected = source_files.get(rel)
        if expected is None or expected["sha256"] != record["sha256"] or expected["bytes"] != record["size"]:
            raise ValueError(f"Shared input checksum mismatch: {raw}")
    with stage_directory(
        destination,
        "inputs:canonical-view",
        {},
        [source / "manifest.json"],
        imported_from={"source": str(source), "purpose": "read-only symlink view"},
    ) as temp:
        for (raw, target), identity in zip(pairs, identities):
            link(raw, temp / target)
        (temp / "dataset.json").write_text(
            json.dumps(
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
                        dict(identity, path=str(target)) for (_, target), identity in zip(pairs, identities)
                    ],
                },
                indent=2,
            )
            + "\n"
        )
    return destination


def migrate_run(source, destination):
    source, destination = Path(source).resolve(strict=True), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Run migration destination must be new")
    destination = destination.resolve()
    if destination == source or source.is_relative_to(destination) or destination.is_relative_to(source):
        raise ValueError("Source and destination must be separate directory trees")
    old_run = json.loads((source / "07_inference/run.json").read_text())
    if old_run.get("options", {}).get("observed", False) is not False:
        raise ValueError("Measured-day-difference runs are not supported")
    alphas = next(iter(old_run["results"].values()))["alphas"]
    layout = RunLayout(destination)
    stages = []
    verified_sources = set()

    def publish(old, new, inputs, writer, expected_old_inputs=None):
        manifest = verify_stage(old, recursive=True, _seen=verified_sources)
        if expected_old_inputs is not None:
            recorded = {Path(e["path"]).resolve() for e in manifest["inputs"]}
            expected = {Path(p).resolve() for p in expected_old_inputs}
            if recorded != expected:
                raise ValueError(f"Original lineage belongs to a different run: {old}")
        origin = {
            "path": str(old),
            "manifest_sha256": describe_inputs([manifest_path(old)])[0]["sha256"],
            "original_implementation": manifest["implementation"],
            "original_signature": manifest["signature"],
        }
        with stage_directory(
            new, manifest["stage"], manifest["config"], inputs, imported_from=origin
        ) as temp:
            writer(Path(old), temp)
        stages.append(str(new.relative_to(destination)))

    for alpha in alphas:
        old = source / "01_raptor" / f"alpha-{alpha:g}"

        def processed_writer(src, dst):
            link(src / "stokes.npy", dst / "stokes.npy")
            shutil.copy2(src / "preprocess.json", dst / "processing.json")

        publish(
            old,
            layout.processed(alpha),
            [Path(e["path"]) for e in verify_stage(old)["inputs"]],
            processed_writer,
        )
        old_response = source / "02_response" / f"alpha-{alpha:g}"

        def response_writer(src, dst):
            for dim in ("1d", "2d"):
                x = np.load(src / f"fits_{dim}.npy", mmap_mode="r", allow_pickle=False)
                np.savez(
                    dst / f"response-{dim}.npz",
                    amplitude=x[..., 0],
                    phase=x[..., 1],
                    baseline=x[..., 2],
                    snapshot_labels=1600 + 10 * np.arange(len(x)),
                )
                with np.load(dst / f"response-{dim}.npz", allow_pickle=False) as z:
                    for i, key in enumerate(("amplitude", "phase", "baseline")):
                        np.testing.assert_array_equal(z[key], x[..., i])
            shutil.copy2(src / "fit.json", dst / "response.json")

        publish(
            old_response,
            layout.response(alpha),
            stage_files(layout.processed(alpha)),
            response_writer,
            expected_old_inputs=stage_files(old),
        )
    for dim in ("1d", "2d"):

        def template_writer(src, dst):
            link(src / "templates.npz", dst / "templates.npz")
            shutil.copy2(src / "summarize.json", dst / "summary.json")

        publish(
            source / "03_templates" / dim,
            layout.templates(dim),
            [layout.response(a) / f"response-{dim}.npz" for a in alphas],
            template_writer,
            expected_old_inputs=[source / "02_response" / f"alpha-{a:g}" / f"fits_{dim}.npy" for a in alphas],
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
            for k, day in enumerate((95, 96, 100, 101)):
                daydir = dst / date_label(day)
                daydir.mkdir()
                link(src / f"stokes_{day}.npy", daydir / "stokes.npy")
                link(src / f"evpa_1d_{day}.npy", daydir / "evpa-1d.npy")
                fields = {}
                for dim, z in [("1d", one), ("2d", two)]:
                    for key in (
                        "mean",
                        "std",
                        "center",
                        "unwrapped_mean",
                        "concentration",
                        "circular_std",
                        "count",
                    ):
                        fields[f"evpa_{key}_{dim}"] = z[key][k]
                fields["intensity_mean"] = two["intensity_mean"][k]
                fields["covariance_1d"] = one["covariance"][k]
                np.savez(daydir / "summary.npz", **fields)

    old_eht = source / "04_eht"
    publish(old_eht, layout.eht, [Path(e["path"]) for e in verify_stage(old_eht)["inputs"]], eht_writer)
    old_bg = source / "01_raptor_background"
    publish(
        old_bg,
        layout.background_processed,
        [Path(e["path"]) for e in verify_stage(old_bg)["inputs"]],
        processed_writer,
    )

    def background_writer(src, dst):
        shutil.copy2(src / "background.json", dst / "summary.json")
        shutil.copy2(src / "background_diagnostics.npz", dst / "evpa-difference.npz")

    publish(
        source / "05_background",
        layout.background,
        stage_files(layout.background_processed) + [layout.eht / "mask.npy"],
        background_writer,
        expected_old_inputs=stage_files(source / "01_raptor_background") + [source / "04_eht/mask.npy"],
    )
    for dim in ("1d", "2d"):

        def observation_writer(src, dst):
            shutil.copy2(src / "observations.npz", dst / "observations.npz")
            shutil.copy2(src / "observations.json", dst / "observations.json")

        publish(
            source / "06_observations" / dim,
            layout.observations(dim),
            stage_files(layout.eht) + [layout.templates(dim) / "templates.npz"],
            observation_writer,
            expected_old_inputs=stage_files(source / "04_eht")
            + [source / "03_templates" / dim / "templates.npz"],
        )
    inference_inputs = [manifest_path(layout.background), layout.background / "summary.json"]
    for dim in ("1d", "2d"):
        inference_inputs += [
            manifest_path(layout.observations(dim)),
            manifest_path(layout.templates(dim)),
            layout.observations(dim) / "observations.npz",
            layout.templates(dim) / "templates.npz",
        ]

    def inference_writer(src, dst):
        modes = ["std", "cor", "cor_grmhd"]
        for dim in ("1d", "2d"):
            values = []
            reference = None
            for mode in modes:
                with np.load(src / f"{mode}_{dim}.npz", allow_pickle=False) as z:
                    current = {k: z[k].copy() for k in ("alphas", "cuts", "variant", "observation", "prior")}
                if str(current["observation"].item()) not in {"null", "expected"}:
                    raise ValueError("Only zero-day-difference results can be migrated")
                if reference is not None:
                    for key in ("alphas", "variant", "observation", "prior"):
                        np.testing.assert_array_equal(reference[key], current[key])
                reference = current
                values.append(current["cuts"])
                shutil.copy2(src / f"{mode}_{dim}.json", dst / f"posterior-{mode}-{dim}.json")
            np.savez(
                dst / f"limits-{dim}.npz",
                alphas=reference["alphas"],
                cuts=np.stack(values),
                modes=modes,
                **{key: reference[key] for key in ("variant", "observation", "prior")},
            )
        shutil.copy2(src / "run.json", dst / "run.json")

    old_inference_inputs = [source / "05_background/stage.json", source / "05_background/background.json"]
    for dim in ("1d", "2d"):
        old_inference_inputs += [
            source / "06_observations" / dim / "stage.json",
            source / "03_templates" / dim / "stage.json",
            source / "06_observations" / dim / "observations.npz",
            source / "03_templates" / dim / "templates.npz",
        ]
    publish(
        source / "07_inference",
        layout.inference,
        inference_inputs,
        inference_writer,
        expected_old_inputs=old_inference_inputs,
    )

    def figures_writer(src, dst):
        for path in src.iterdir():
            if path.is_file() and path.name != "stage.json":
                name = "limits-by-dimension.png" if path.name == "limits_by_dimension.png" else path.name
                shutil.copy2(path, dst / name)

    publish(
        source / "08_figures",
        layout.figures,
        stage_files(layout.inference),
        figures_writer,
        expected_old_inputs=stage_files(source / "07_inference"),
    )
    (destination / "migration.json").write_text(
        json.dumps(
            {
                "source": str(source),
                "layout": "canonical-run-v1",
                "scientific_recomputation": False,
                "large_arrays": "symlinked; source files unchanged",
                "converted_response_arrays": "exact equality verified",
                "completed_stages": stages,
                "source_run_metadata": old_run,
            },
            indent=2,
        )
        + "\n"
    )
    return destination
