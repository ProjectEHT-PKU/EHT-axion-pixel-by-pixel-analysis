"""Verified job HDF -> existing canonical mmap products, with live ancestry."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..model import omega
from ..pipeline.common import (
    manifest_path, require_declared_outputs, stage_directory, verification_session, verify_stage,
)
from ..pipeline.layout import alpha_label
from ..pipeline.raptor import convert_hdf5


def _job_files(run_root, task, build_inputs, verified=None):
    directory = Path(run_root) / "jobs" / task["id"]
    record = verify_stage(directory, recursive=True, allow_imported=False, _seen=verified)
    from .jobs import _prepare_task

    prepared, expected, raw_inputs = _prepare_task(task)
    if record["stage"] != "producer-job" or any(record["config"].get(k) != v for k, v in expected.items()):
        raise ValueError("Completed job belongs to a different generation task")
    require_declared_outputs(directory, [directory / "success.json", directory / "output.json"])
    success = json.loads((directory / "success.json").read_text())
    if (success.get("schema") != 1 or success.get("returncode") != 0
            or success.get("timed_out") is not False or success.get("output_validated") is not True):
        raise ValueError("Job lacks a successful validated execution receipt")
    if {Path(e["path"]).resolve() for e in record["inputs"]} != {
        Path(p).resolve() for p in build_inputs + raw_inputs
    }:
        raise ValueError("Job lineage does not declare the selected build and source inputs")
    from .outputs import read_core_constants, validate_job_output

    prepared["core_constants"] = read_core_constants(build_inputs[-1])
    prepared["core_constants_source"] = str(build_inputs[-1])
    if json.loads((directory / "output.json").read_text()) != validate_job_output(directory, prepared):
        raise ValueError("Completed job output validation differs from its receipt")
    hdfs = [directory / entry["path"] for entry in record["outputs"] if entry["path"].endswith(".h5")]
    if len(hdfs) != 1:
        raise ValueError(f"Expected one declared job HDF: {directory}")
    return hdfs[0], manifest_path(directory)


def _mapping(resolution, frequency):
    blocks = resolution // 10
    return dict(
        datasets={name: f"{name}{frequency:e}" for name in "IQUV"},
        image_shape=(resolution, resolution), pixel_shape=(10, 10),
        block_origins=[(10 * (b // blocks), 10 * (b % blocks)) for b in range(blocks**2)],
    )


def _identity(config, *, alpha, resolution, labels, times, state_count):
    p = config["producer"]
    return dict(
        schema="raptor-generated-v1", physical=config["physical"],
        initial_fov_rg=config["image"]["initial_fov_rg"], alpha=alpha,
        resolution=resolution, snapshot_labels=labels, times_tg=times,
        state_count=state_count, reference_norm=p["response"]["reference_norm"],
        phase_count=16, principal_n=config["inference"]["principal_n"],
        omega=None if alpha is None else omega(alpha, config["inference"]["principal_n"]),
        frequency_hz=p["plasma"]["frequency_hz"],
        producer=p, axes=["time", "state", "x", "y", "stokes"],
        channels=list("IQUV"), intensity_units="Jy/sr", frame_time_source="GRMHD footer",
        phase_convention="state0 baseline; state k+1 has phase k*pi/8",
    )


def deliver_stokes(plan, run_root, *, job_roots=None, output_root=None):
    """Publish verified groups, preserving each job's real build/source ancestry.

    Optional job_roots merges independent runs without copying or re-signing
    their job receipts. Every planned task must exist in exactly one source.
    output_root is a new delivery container, never a source run or raw directory.
    """
    from ..pipeline.common import reject_directory_symlinks, source_closure

    root = Path(run_root).resolve()
    roots = [root] if job_roots is None else [Path(p).resolve(strict=True) for p in job_roots]
    if not roots or len(set(roots)) != len(roots):
        raise ValueError("Delivery requires distinct nonempty job roots")
    target = root if output_root is None else Path(output_root).absolute()
    reject_directory_symlinks(target)
    target = target.resolve()
    if output_root is not None:
        for source in roots:
            if target == source or target.is_relative_to(source) or source.is_relative_to(target):
                raise ValueError("Delivery output must be disjoint from source job roots")
    config = plan["config"]
    tasks = plan["tasks"]
    outputs = []
    with verification_session():
        from .jobs import _build_products, _prepare_task, _validate_plan_binding

        lookup, task_sources, assignments = {}, {}, {source: [] for source in roots}
        for task in tasks:
            key = (task["kind"], task["resolution"], task["label"], task.get("alpha"),
                   task["phase_index"])
            if key in lookup or task["id"] in task_sources:
                raise ValueError("Duplicate planned task")
            matches = [source for source in roots
                       if (source / "jobs" / task["id"]).exists()
                       or (source / "jobs" / task["id"]).is_symlink()]
            if len(matches) != 1:
                raise ValueError(f"Expected exactly one source for task {task['id']}, found {len(matches)}")
            source = matches[0]
            lookup[key] = task
            task_sources[task["id"]] = source
            assignments[source].append(task)
        build_inputs, verified = {}, set()
        for source, selected in assignments.items():
            if not selected:
                continue
            build_dir = source / "build"
            build_manifest = verify_stage(build_dir, recursive=True, allow_imported=False, _seen=verified)
            binary, grid_path, header = _build_products(build_dir, build_manifest)
            _validate_plan_binding(plan, [_prepare_task(t) for t in selected],
                                   build_manifest, build_dir, grid_path)
            build_inputs[source] = [manifest_path(build_dir), binary, grid_path, header]
        # Reject delivery roots inside any actual consumed raw/code source tree.
        if output_root is not None:
            sources = [Path(t["source_path"]) for t in tasks]
            sources += [Path(t["table_path"]) for t in tasks if t.get("table_path")]
            sources += [p for values in build_inputs.values() for p in values]
            for source in source_closure(sources):
                if (source == target or source.is_relative_to(target)
                        or target.is_relative_to(source.parent)):
                    raise ValueError(f"Delivery output overlaps a consumed source: {source}")
        source_records = [dict(task_id=t["id"], root=str(task_sources[t["id"]])) for t in tasks]
        response_res = config["producer"]["response"]["resolution"]
        background_res = config["producer"]["background"]["resolution"]
        snap = config["producer"]["snapshots"]
        labels = list(range(snap["start"], snap["stop"] + 1, snap["step"]))
        # Preserve current res100 whole-file input contract. High-res delivery is
        # per 10-time shard, compatible with the existing inclusive-step layout.
        groups = [(a, response_res, labels, 17) for a in config["data"]["alphas"]]
        groups += [(None, background_res, labels[i:i + 10], 1) for i in range(0, len(labels), 10)]
        for alpha, resolution, group_labels, states in groups:
            grid, inputs, times = [], [], []
            for label in group_labels:
                baseline = lookup[("baseline", resolution, label, None, 0)]
                selected = [baseline]
                if alpha is not None:
                    selected += [lookup[("response", resolution, label, alpha, k)] for k in range(16)]
                row = []
                for task in selected:
                    source = task_sources[task["id"]]
                    hdf, receipt = _job_files(source, task, build_inputs[source], verified)
                    row.append(hdf)
                    inputs += [receipt, hdf]
                grid.append(row)
                times.append(baseline["time_tg"])
            alias = config["background"]["alpha"] if alpha is None else alpha
            base = target / "data" / "raptor" / f"res-{resolution}" / alpha_label(alias)
            if alpha is None:
                destination = base / f"t{group_labels[0]}-t{group_labels[-1]}"
                filename = f"stokes-t{group_labels[0]}-t{group_labels[-1]}.npy"
            else:
                destination, filename = base, "stokes.npy"
            identity = _identity(config, alpha=alpha, resolution=resolution, labels=group_labels,
                                 times=times, state_count=states)
            with stage_directory(destination, "producer:delivery", identity,
                                 list(dict.fromkeys(inputs))) as temp:
                if temp is not None:
                    result = convert_hdf5(grid, temp, **_mapping(resolution, identity["frequency_hz"]))
                    if filename != "stokes.npy":
                        result["stokes"].rename(temp / filename)
                    (temp / "generated.json").write_text(json.dumps(identity, indent=2) + "\n")
                    used = {t["id"] for label in group_labels for t in tasks
                            if t["label"] == label and t["resolution"] == resolution
                            and (t["kind"] == "baseline" or t.get("alpha") == alpha)}
                    (temp / "sources.json").write_text(json.dumps(
                        [r for r in source_records if r["task_id"] in used], indent=2) + "\n")
            outputs.append(destination / filename)
    return outputs


def generated_inputs(paths, config, *, alpha, resolution):
    """Validate generated source identity, return declared edges and real labels.

    Legacy raw data remains usable without producer configuration. A production
    profile must never consume unlabeled legacy raw data by silent fallback.
    """
    paths = [Path(p) for p in paths]
    inputs, labels = list(paths), []
    generated = [(p.parent / "generated.json").is_file() for p in paths]
    if not any(generated):
        for path in paths:
            for name in ("manifest.json", "stage.json"):
                marker = path.parent / name
                if marker.is_file():
                    record = json.loads(marker.read_text())
                    if (str(record.get("stage", "")).startswith("producer")
                            or record.get("config", {}).get("schema") == "raptor-generated-v1"
                            or any(e.get("path") == "generated.json" for e in record.get("outputs", []))):
                        raise ValueError("Generated delivery metadata is missing; refusing legacy downgrade")
        if "producer" in config:
            raise ValueError("Production profile requires verified generated RAPTOR inputs")
        return inputs, None
    if not all(generated):
        raise ValueError("Mixed generated and legacy RAPTOR inputs")
    if "producer" not in config:
        raise ValueError("Generated RAPTOR inputs require their production science profile")
    for path in paths:
        directory = path.parent
        verify_stage(directory, recursive=True, allow_imported=False)
        metadata = directory / "generated.json"
        require_declared_outputs(directory, [path, metadata])
        record = json.loads(metadata.read_text())
        state_count = 1 if resolution == config["producer"]["background"]["resolution"] else 17
        actual_alpha = None if state_count == 1 else alpha
        expected = _identity(config, alpha=actual_alpha, resolution=resolution,
                             labels=record["snapshot_labels"], times=record["times_tg"],
                             state_count=state_count)
        if record != expected:
            raise ValueError("Generated RAPTOR science identity differs from analysis configuration")
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        if a.shape != (len(record["snapshot_labels"]), state_count, resolution, resolution, 4):
            raise ValueError("Generated RAPTOR shape/identity mismatch")
        del a
        labels += record["snapshot_labels"]
        inputs += [manifest_path(directory), metadata]
    snap = config["producer"]["snapshots"]
    if labels != list(range(snap["start"], snap["stop"] + 1, snap["step"])):
        raise ValueError("Generated labels differ from selected snapshot interval")
    return list(dict.fromkeys(inputs)), labels


def validate_processed_source(manifest, config, *, alpha, resolution):
    """Recheck producer contracts from independent inference/verify entrypoints."""
    raw = [Path(e["path"]) for e in manifest["inputs"] if Path(e["path"]).suffix == ".npy"]
    if not raw:
        raise ValueError("Preprocess stage declares no raw Stokes input")
    required, _ = generated_inputs(raw, config, alpha=alpha, resolution=resolution)
    actual = {Path(e["path"]).resolve() for e in manifest["inputs"]}
    if actual != {p.resolve() for p in required}:
        raise ValueError("Preprocess source lineage omits generated identity or contains foreign inputs")
