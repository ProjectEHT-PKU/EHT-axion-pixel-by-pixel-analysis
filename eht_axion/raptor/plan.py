"""RAPTOR task planning: read-only BHAC snapshot/radial-table checks and the JSON task list.

No build, subprocess, output directory, discovery of earlier outputs or simulation is
performed here. Execution must hash live source/table/build bytes and bind those
hashes to its own session records. Only the small BHAC footers and the radial
tables are read here.

Task ids are SHA-256 digests of exactly what determines one RAPTOR image: the
snapshot (name + footer metadata), the 19 model.in values except the table path,
the radial-table content hash and spin, the core grid and the build settings.
Analysis, sky-grid and eht-imaging settings never enter a task id.
"""

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import struct

import numpy as np

from ..experiment import ray_tracing_config, validate_ray_tracing
from .. import physics
from .layout import alpha_label

SCHEMA_VERSION = 2
# The core applies AXION_PHASE * pi/8: 16 phases cover one period (state k+1 = phase k).
PHASE_COUNT = 16
# model/bhac/model.c: EOF-40 is eight 32-bit ints then a double; immediately
# preceding it are ndimini ints nx[] and neqparini doubles. NSPIN is 3.
FOOTER = struct.Struct("<8id")


def alpha_key(alpha):
    """raptor.alpha_table keys omit the filesystem layout's ``alpha-`` prefix."""
    return alpha_label(alpha).removeprefix("alpha-")


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode("utf-8")).hexdigest()


def _path(path, name):
    path = Path(path).resolve()
    text = str(path)
    if len(text.encode("utf-8")) > 255 or any(c.isspace() or c == "\0" for c in text):
        raise ValueError(f"{name} must be a single C path token of at most 255 UTF-8 bytes: {path}")
    return path


def file_sha256(path: Path) -> str:
    """Stream once per caller-selected unique input; never load a fluid into RAM."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_snapshot_footer(path: Path, grid: dict, fluid_spin: float) -> dict:
    """Read only the little-endian Linux BHAC footer, not its fluid or forest.

    Payload length is checked against a necessary lower bound, not certified as
    a full AMR-tree or body validation. This is deliberately a light preflight.
    """
    path = Path(path)
    with path.open("rb") as source:
        size = source.seek(0, 2)
        if size < FOOTER.size:
            raise ValueError(f"Truncated BHAC footer: {path}")
        source.seek(-FOOTER.size, 2)
        nleafs, levmax, ndim, ndir, nw, nws, neq, iteration, time = FOOTER.unpack(source.read(40))
        if (nleafs <= 0 or not 1 <= levmax <= 30 or ndim != 3 or ndir != 3
                or not 10 <= nw <= 64 or not 0 <= nws <= 64 or not 4 <= neq <= 64 or iteration < 0):
            raise ValueError(f"Unsupported or corrupt BHAC footer layout (little-endian 3D): {path}")
        if not math.isfinite(time):
            raise ValueError(f"Nonfinite BHAC footer time: {path}")
        meta_size = ndim * 4 + neq * 8
        if size < FOOTER.size + meta_size:
            raise ValueError(f"Truncated BHAC nx/neqpar metadata: {path}")
        source.seek(-(FOOTER.size + meta_size), 2)
        nx = list(struct.unpack("<3i", source.read(12)))
        neqpar = list(struct.unpack(f"<{neq}d", source.read(neq * 8)))
    if any(n <= 0 or total % n for n, total in zip(nx, grid["nxlone"])):
        raise ValueError(f"BHAC block nx must divide the configured core nxlone grid: {path}")
    if not all(math.isfinite(value) for value in neqpar) or neqpar[0] <= 1:
        raise ValueError(f"Invalid BHAC neqpar/gamma: {path}")
    if not math.isclose(neqpar[3], fluid_spin, rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"BHAC footer spin must be raptor.fluid_spin={fluid_spin}, got {neqpar[3]}: {path}")
    cells = math.prod(nx)
    payload = nleafs * (cells * nw + math.prod(n + 1 for n in nx) * nws) * 8
    minimum_size = payload + nleafs * 4 + FOOTER.size + meta_size
    if size < minimum_size:
        raise ValueError(f"Truncated BHAC payload by footer-derived size lower bound: {path}")
    return {
        "byte_order": "little", "size_bytes": size, "nleafs": nleafs,
        "levmax": levmax, "ndim": ndim, "ndir": ndir, "nw": nw, "nws": nws,
        "neqpar_count": neq, "iteration": iteration, "nx": nx, "neqpar": neqpar,
        "time_tg": time, "spin": neqpar[3], "gamma": neqpar[0],
        "minimum_size_bytes": minimum_size,
        "body_verified": False,
    }


def validate_radial_table(path: Path, required_alphas) -> dict:
    """Load one small headerless alpha,r,R table and mirror core grid checks.

    CSV cannot establish its spin or normalization provenance. Spin is declared
    by the profile; amplitudes are validated but never rescaled or normalized.
    Radius coverage during actual ray tracing remains the core's domain check.
    """
    try:
        values = np.loadtxt(path, delimiter=",", comments=None, ndmin=2)
    except (ValueError, OSError) as exc:
        raise ValueError(f"Invalid headerless alpha,r,R CSV: {path}: {exc}") from exc
    if values.ndim != 2 or values.shape[1] != 3 or values.shape[0] < 2:
        raise ValueError(f"Radial CSV must have exactly three columns and at least two rows: {path}")
    if not np.all(np.isfinite(values)) or np.any(values[:, 2] < 0):
        raise ValueError(f"Radial CSV requires finite alpha,r,R with R >= 0: {path}")
    # Keep original row order: sorting would conceal a core-incompatible file.
    boundaries = np.r_[0, np.flatnonzero(np.diff(values[:, 0]) != 0) + 1, len(values)]
    lengths = np.diff(boundaries)
    nr = int(lengths[0])
    alphas = values[boundaries[:-1], 0]
    if nr < 2 or np.any(lengths != nr) or np.any(np.diff(alphas) <= 0):
        raise ValueError(f"Radial CSV requires a rectangular increasing alpha grid with nr >= 2: {path}")
    radius = values[:nr, 1]
    if np.any(np.diff(radius) <= 0):
        raise ValueError(f"Radial CSV radius grid must be strictly increasing: {path}")
    for start in boundaries[:-1]:
        if not np.array_equal(values[start:start + nr, 1], radius):
            raise ValueError(f"Radial CSV alpha blocks must share identical radius grids: {path}")
    for alpha in required_alphas:
        if not alphas[0] <= alpha <= alphas[-1]:
            raise ValueError(f"Radial CSV does not cover alpha={alpha}: {path}")
    return {
        "row_count": int(len(values)), "alpha_count": int(len(alphas)), "radius_count": nr,
        "alpha_min": float(alphas[0]), "alpha_max": float(alphas[-1]),
        "radius_min": float(radius[0]), "radius_max": float(radius[-1]),
        "amplitude_min": float(values[:, 2].min()), "amplitude_max": float(values[:, 2].max()),
        "normalization": "unchanged", "provenance_verified": False,
    }


def _parameters(config, resolution, alpha, omega, norm, phase, table_path):
    """Exact 19-field core/run/model.in order. Each entry is LABEL -> UNIT, VALUE."""
    shared, raptor = config["shared"], config["raptor"]
    p = raptor["plasma"]
    entries = [
        ("MBH", "(Msun)", shared["mass_solar"]),
        ("DISTANCE", "(kpc)", shared["distance_mpc"] * 1000),
        ("M_UNIT", "(g)", p["mass_unit_g"]),
        ("R_LOW", "(-)", p["r_low"]), ("R_HIGH", "(-)", p["r_high"]),
        ("INCLINATION", "(deg)", p["inclination_degrees"]),
        ("IMG_WIDTH", "(pixels)", resolution), ("IMG_HEIGHT", "(pixels)", resolution),
        ("CAM_SIZE_X", "(Rg)", raptor["camera_fov_rg"]),
        ("CAM_SIZE_Y", "(Rg)", raptor["camera_fov_rg"]),
        ("FREQS_PER_DEC", "(-)", 1), ("FREQ_MIN", "(Hz)", p["frequency_hz"]),
        ("STEPSIZE", "(-)", p["step_size"]),
        ("AXION_ALPHA", "(-)", 0.0 if alpha is None else alpha),
        ("AXION_OMEGA", "(-)", omega), ("AXION_NORM", "(-)", norm),
        ("AXION_PHASE", "(-)", phase),
        ("AXION_TABLE", "(-)", "disabled.csv" if table_path is None else table_path),
        ("MAX_LEVEL", "(-)", 1),
    ]
    return {name: {"unit": unit, "value": value} for name, unit, value in entries}


def task_identity(snapshot, parameters, table, config):
    """Digest of everything that determines one RAPTOR image (see module docstring)."""
    raptor = config["raptor"]
    values = {name: item["value"] if isinstance(item, dict) else item for name, item in parameters.items()}
    values.pop("AXION_TABLE")
    return _digest({
        "schema_version": SCHEMA_VERSION,
        "snapshot": {"name": snapshot["snapshot_name"], "metadata": snapshot["metadata"]},
        "parameters": values,
        "table": None if table is None else {"spin": table["spin"], "sha256": table["sha256"]},
        "grid": raptor["grid"], "build": raptor["build"],
    })


def plan_generation(config: dict, fluid_root: Path, radial_paths: dict[str, Path],
                    build_path: Path | None = None, *, hash_inputs: bool = False) -> dict:
    """Plan inclusive configured labels; no lexicographic listing or hidden slice.

    ``config`` is the experiment (or its ray-tracing subset). ``build_path`` is an
    opaque caller-selected build location, recorded only. Radial tables are small
    and always hashed (their content defines task ids); ``hash_inputs`` also
    streams every snapshot once. Execution still has to hash live inputs.
    """
    config = deepcopy(validate_ray_tracing(ray_tracing_config(config)))
    shared, p = config["shared"], config["raptor"]
    if type(hash_inputs) is not bool:
        raise ValueError("hash_inputs must be a boolean")
    if not isinstance(radial_paths, dict) or set(radial_paths) != set(p["radial_tables"]):
        raise ValueError(f"radial_paths requires exactly {sorted(p['radial_tables'])}")
    config_sha = _digest(config)
    fluid_root = Path(fluid_root).resolve()
    build = None if build_path is None else str(_path(build_path, "build_path"))

    snapshots = []
    start, stop, step = (p["snapshots"][key] for key in ("start", "stop", "step"))
    for label in range(start, stop + 1, step):
        name = f"data{label:06d}.dat"
        path = _path(fluid_root / name, "source_path")
        metadata = read_snapshot_footer(path, p["grid"], p["fluid_spin"])
        if snapshots:
            previous = snapshots[-1]
            expected = shared["label_step_tg"] * (label - previous["label"])
            if not math.isclose(metadata["time_tg"] - previous["time_tg"], expected,
                                rel_tol=0, abs_tol=1e-8):
                raise ValueError(f"BHAC footer time cadence disagrees with labels: {name}")
            if metadata["gamma"] != previous["metadata"]["gamma"]:
                raise ValueError(f"BHAC gamma changes across snapshots: {name}")
        snapshots.append({
            "label": label, "snapshot_name": name, "time_tg": metadata["time_tg"],
            "source_path": str(path), "source_sha256": file_sha256(path) if hash_inputs else None,
            "metadata": metadata,
        })

    tables, table_cache = {}, {}
    for table_id in p["radial_tables"]:
        path = _path(radial_paths[table_id], f"radial_paths.{table_id}")
        required = [alpha for alpha in shared["alphas"]
                    if p["alpha_table"][alpha_key(alpha)] == table_id]
        # Unique physical files are loaded once, even if explicitly mapped twice.
        key = str(path)
        if key not in table_cache:
            table_cache[key] = (validate_radial_table(path, []), file_sha256(path))
        metadata, digest = deepcopy(table_cache[key][0]), table_cache[key][1]
        if any(not metadata["alpha_min"] <= alpha <= metadata["alpha_max"] for alpha in required):
            raise ValueError(f"Radial CSV {table_id} does not cover required alphas {required}: {path}")
        tables[table_id] = {"path": key, "spin": p["radial_tables"][table_id]["spin"],
                            "sha256": digest, "metadata": metadata}

    tasks = []

    def append_task(snapshot, resolution, alpha=None, phase=0):
        enabled = alpha is not None
        table_id = p["alpha_table"][alpha_key(alpha)] if enabled else None
        table = tables[table_id] if enabled else None
        omega = float(physics.omega(alpha, shared["principal_n"])) if enabled else 0.0
        norm = p["response"]["reference_norm"] if enabled else 0.0
        table_path = None if table is None else table["path"]
        parameters = _parameters(config, resolution, alpha, omega, norm, phase, table_path)
        tasks.append({
            "id": task_identity(snapshot, parameters, table, config), "label": snapshot["label"],
            "snapshot_name": snapshot["snapshot_name"], "time_tg": snapshot["time_tg"],
            "source_path": snapshot["source_path"], "source_sha256": snapshot["source_sha256"],
            "kind": "response" if enabled else "baseline", "resolution": resolution,
            "alpha": alpha, "omega": omega, "norm": norm,
            "phase_index": phase, "output_slot": phase + 1 if enabled else 0,
            "table_id": table_id, "table_path": table_path,
            "table_spin": None if table is None else table["spin"],
            "table_sha256": None if table is None else table["sha256"],
            "parameters": parameters,
        })

    response_resolution, background_resolution = p["response"]["resolution"], p["background"]["resolution"]
    for snapshot in snapshots:
        # Deduplicate baseline per snapshot/resolution, never per alpha.
        for resolution in dict.fromkeys((response_resolution, background_resolution)):
            append_task(snapshot, resolution)
        for alpha in shared["alphas"]:
            for phase in range(PHASE_COUNT):
                append_task(snapshot, response_resolution, alpha, phase)
    if len({task["id"] for task in tasks}) != len(tasks):
        raise ValueError("Duplicate task identities")

    baseline_count = sum(task["kind"] == "baseline" for task in tasks)
    response_count = len(tasks) - baseline_count
    # Only unavoidable numerical buffers are counted, not HDF5/container,
    # allocator/pointer, AMR-tree, plot or checkpoint overhead.
    stokes_bytes = sum(task["resolution"]**2 * 4 * 8 for task in tasks)
    memory_bytes = max(s["metadata"]["nleafs"] * math.prod(s["metadata"]["nx"]) * (8 + 6) * 8
                       for s in snapshots)
    from .build import _arguments, _settings

    return {
        "schema_version": SCHEMA_VERSION, "config_sha256": config_sha, "config": config,
        "build_path": build, "hash_inputs": hash_inputs,
        "verification": "input_hashes_only" if hash_inputs else "unverified",
        "snapshots": snapshots, "radial_tables": tables, "tasks": tasks,
        "counts": {"snapshots": len(snapshots), "baseline": baseline_count,
                   "response": response_count, "background": len(snapshots), "total": len(tasks)},
        "requirements": {
            "source_sha256_required_at_execution": True, "table_sha256_required_at_execution": True,
            "build_sha256_required_at_execution": True, "live_rehash_required_at_execution": True,
            "free_space_check_owner": "caller", "output_stokes_lower_bound_bytes": stokes_bytes,
            "per_process_memory_lower_bound_bytes": memory_bytes,
            "unique_source_bytes": sum(s["metadata"]["size_bytes"] for s in snapshots),
            "official_setup_args": _arguments(_settings(p["build"])), "build_command": ["make"],
            "label_time_offset_tg": snapshots[0]["time_tg"] - start * shared["label_step_tg"],
        },
        "warnings": [
            "Task ids describe planned ray-tracing inputs; "
            "execution must bind live source/table/build SHA256.",
            "Footer/payload lower-bound checks do not validate fluid values or AMR forest semantics.",
            "Radial CSV spin/provenance cannot be inferred from alpha,r,R columns; "
            "configured spins are declarations.",
            "Declared radial-table spins " + ", ".join(
                f"{k}={v['spin']}" for k, v in sorted(p["radial_tables"].items()))
            + f"; fluid spin raptor.fluid_spin={p['fluid_spin']}.",
            "No radial amplitudes were renormalized; ray-radius coverage remains a core runtime check.",
            "Resource estimates are lower bounds only; "
            "caller must check disk/RAM with overhead and concurrency.",
            "This plan builds and launches nothing and reads no earlier scientific outputs.",
        ],
    }
