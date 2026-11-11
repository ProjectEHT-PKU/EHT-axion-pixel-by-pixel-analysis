"""Delivery: RAPTOR job HDF5 -> Stokes cubes (delivered/res-*/{alpha-*,background}/), and checks of
delivered cubes against this experiment (generated_inputs, input_cadence_hours, validate_processed_source)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..physics import omega
from ..provenance import (
    manifest_path, require_declared_outputs, stage_directory, verify_stage,
)
from .convert import convert_hdf5
from .layout import product_directory
from .plan import PHASE_COUNT


def _job_files(run_root, task, build_inputs):
    directory = Path(run_root) / "jobs" / task["id"]
    record = verify_stage(directory)
    from .jobs import _prepare_task

    _, expected, raw_inputs = _prepare_task(task)
    if record["stage"] != "raptor:job" or any(record["config"].get(k) != v for k, v in expected.items()):
        raise ValueError("Completed job belongs to a different generation task")
    require_declared_outputs(directory, [directory / "success.json", directory / "output.json"])
    success = json.loads((directory / "success.json").read_text())
    if (success.get("returncode") != 0
            or success.get("timed_out") is not False or success.get("output_validated") is not True):
        raise ValueError("Job lacks a successful validated execution receipt")
    if {Path(p).resolve() for p in record["inputs"]} != {
        Path(p).resolve() for p in build_inputs + raw_inputs
    }:
        raise ValueError("Job lineage does not declare the selected build and source inputs")
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


def _identity(config, *, alpha, resolution, labels, times):
    """generated.json: the ray-tracing configuration plus this product's coordinates."""
    from ..experiment import ray_tracing_config

    config = ray_tracing_config(config)
    principal_n = config["shared"]["principal_n"]
    return dict(
        schema="raptor-generated", ray_tracing=config, alpha=alpha, resolution=resolution,
        snapshot_labels=labels, times_tg=times,
        state_count=1 if alpha is None else PHASE_COUNT + 1,
        omega=None if alpha is None else omega(alpha, principal_n),
        axes=["time", "state", "x", "y", "stokes"], channels=list("IQUV"), intensity_units="Jy/sr",
        frame_time_source="GRMHD footer",
        phase_convention="state0 axion-off; state k+1 has axion phase k*pi/8",
    )


def deliver_stokes(plan, run_root):
    """Publish the Stokes cubes of the completed jobs of run_root into run_root/delivered.

    run_root holds build/, jobs/<task id>/ and delivered/res-<n>/{alpha-<a>,background/tA-tB}/.
    """
    root = Path(run_root).resolve()
    config = plan["config"]
    tasks = plan["tasks"]
    outputs = []
    from .jobs import _build_products, _prepare_task, _validate_plan_binding

    lookup = {}
    for task in tasks:
        key = (task["kind"], task["resolution"], task["label"], task.get("alpha"), task["phase_index"])
        if key in lookup:
            raise ValueError("Duplicate planned task")
        lookup[key] = task
    build_dir = root / "build"
    build_manifest = verify_stage(build_dir)
    binary, grid_path, header = _build_products(build_dir, build_manifest)
    _validate_plan_binding(plan, [_prepare_task(t) for t in tasks], build_manifest, build_dir, grid_path)
    build_inputs = [manifest_path(build_dir), binary, grid_path, header]
    response_res = config["raptor"]["response"]["resolution"]
    background_res = config["raptor"]["background"]["resolution"]
    snap = config["raptor"]["snapshots"]
    labels = list(range(snap["start"], snap["stop"] + 1, snap["step"]))
    # Response products are one cube per alpha; background products are
    # shards of at most 10 snapshots on the inclusive label grid.
    groups = [(a, response_res, labels) for a in config["shared"]["alphas"]]
    groups += [(None, background_res, labels[i:i + 10]) for i in range(0, len(labels), 10)]
    for alpha, resolution, group_labels in groups:
        grid, inputs, times = [], [], []
        for label in group_labels:
            baseline = lookup[("baseline", resolution, label, None, 0)]
            selected = [baseline]
            if alpha is not None:
                selected += [lookup[("response", resolution, label, alpha, k)]
                             for k in range(PHASE_COUNT)]
            row = []
            for task in selected:
                hdf, receipt = _job_files(root, task, build_inputs)
                row.append(hdf)
                inputs += [receipt, hdf]
            grid.append(row)
            times.append(baseline["time_tg"])
        base = product_directory(root / "delivered", alpha, resolution)
        if alpha is None:
            destination = base / f"t{group_labels[0]}-t{group_labels[-1]}"
            filename = f"stokes-t{group_labels[0]}-t{group_labels[-1]}.npy"
        else:
            destination, filename = base, "stokes.npy"
        identity = _identity(config, alpha=alpha, resolution=resolution, labels=group_labels,
                             times=times)
        with stage_directory(destination, "raptor:delivery", identity,
                             list(dict.fromkeys(inputs))) as temp:
            if temp is not None:
                frequency = config["raptor"]["plasma"]["frequency_hz"]
                result = convert_hdf5(grid, temp, **_mapping(resolution, frequency))
                if filename != "stokes.npy":
                    result["stokes"].rename(temp / filename)
                (temp / "generated.json").write_text(json.dumps(identity, indent=2) + "\n")
        outputs.append(destination / filename)
    return outputs


def generated_inputs(paths, config, *, alpha, resolution):
    """Verify delivered RAPTOR products against this experiment; return inputs and labels.

    ``alpha`` is None for the axion-off background products.

    Each product must be a complete raptor:delivery stage whose generated.json
    equals the identity this experiment would deliver (ray-tracing settings,
    alpha, resolution); the labels must cover raptor.snapshots exactly.
    """
    paths = [Path(p) for p in paths]
    inputs, labels = list(paths), []
    for path in paths:
        directory = path.parent
        verify_stage(directory)
        metadata = directory / "generated.json"
        require_declared_outputs(directory, [path, metadata])
        record = json.loads(metadata.read_text())
        expected = _identity(config, alpha=alpha, resolution=resolution,
                             labels=record.get("snapshot_labels"), times=record.get("times_tg"))
        if record != expected:
            raise ValueError(f"Generated RAPTOR identity differs from this experiment: {metadata}")
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        if a.shape != (len(record["snapshot_labels"]), record["state_count"], resolution, resolution, 4):
            raise ValueError("Generated RAPTOR shape/identity mismatch")
        del a
        labels += record["snapshot_labels"]
        inputs += [manifest_path(directory), metadata]
    snap = config["raptor"]["snapshots"]
    if labels != list(range(snap["start"], snap["stop"] + 1, snap["step"])):
        raise ValueError("Generated labels differ from raptor.snapshots")
    return list(dict.fromkeys(inputs)), labels


def input_cadence_hours(paths, config):
    """Frame spacing in hours from the recorded GRMHD footer times (must be uniform)."""
    from ..experiment import physical_scale

    times = []
    for path in paths:
        times.extend(json.loads((Path(path).parent / "generated.json").read_text())["times_tg"])
    delta = np.diff(np.asarray(times, dtype=float))
    if (not len(delta) or not np.isfinite(delta).all() or np.any(delta <= 0)
            or not np.allclose(delta, delta[0], rtol=1e-12, atol=0)):
        raise ValueError("Background requires a uniform finite footer-time grid")
    return float(delta[0]) * physical_scale(config).tg_hours


def validate_processed_source(manifest, config, *, alpha, resolution):
    """Recheck that a preprocess stage consumed exactly the delivered products."""
    raw = [Path(p) for p in manifest["inputs"] if Path(p).suffix == ".npy"]
    if not raw:
        raise ValueError("Preprocess stage declares no raw Stokes input")
    required, _ = generated_inputs(raw, config, alpha=alpha, resolution=resolution)
    actual = {Path(p).resolve() for p in manifest["inputs"]}
    if actual != {p.resolve() for p in required}:
        raise ValueError("Preprocess source lineage omits generated identity or contains foreign inputs")
