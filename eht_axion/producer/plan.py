"""Read-only BHAC producer preflight and JSON-serializable task planning.

No build, subprocess, output directory, old-output discovery or simulation is
performed here. Plan IDs describe intended science, NOT verified input content.
Execution must hash live source/table/build bytes and bind those hashes to its
own session records. The default only reads the small BHAC footers and tables.
"""

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import struct

import numpy as np

from .. import model
from .config import canonical_alpha, validate_production_config

SCHEMA_VERSION = 1
# model/bhac/model.c: EOF-40 is eight 32-bit ints then a double; immediately
# preceding it are ndimini ints nx[] and neqparini doubles. NSPIN is 3.
FOOTER = struct.Struct("<8id")
EXPECTED_SPIN = 0.9375


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


def read_snapshot_footer(path: Path, grid: dict) -> dict:
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
    if not math.isclose(neqpar[3], EXPECTED_SPIN, rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"BHAC actual spin must be {EXPECTED_SPIN}, got {neqpar[3]}: {path}")
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
    p, physical = config["producer"]["plasma"], config["physical"]
    # Exact 19-field core/run/model.in order. Each entry is LABEL -> UNIT, VALUE.
    entries = [
        ("MBH", "(Msun)", physical["mass_solar"]),
        ("DISTANCE", "(kpc)", physical["distance_mpc"] * 1000),
        ("M_UNIT", "(g)", p["mass_unit_g"]),
        ("R_LOW", "(-)", p["r_low"]), ("R_HIGH", "(-)", p["r_high"]),
        ("INCLINATION", "(deg)", p["inclination_degrees"]),
        ("IMG_WIDTH", "(pixels)", resolution), ("IMG_HEIGHT", "(pixels)", resolution),
        ("CAM_SIZE_X", "(Rg)", config["image"]["initial_fov_rg"]),
        ("CAM_SIZE_Y", "(Rg)", config["image"]["initial_fov_rg"]),
        ("FREQS_PER_DEC", "(-)", 1), ("FREQ_MIN", "(Hz)", p["frequency_hz"]),
        ("STEPSIZE", "(-)", p["step_size"]),
        ("AXION_ALPHA", "(-)", 0.0 if alpha is None else alpha),
        ("AXION_OMEGA", "(-)", omega), ("AXION_NORM", "(-)", norm),
        ("AXION_PHASE", "(-)", phase),
        ("AXION_TABLE", "(-)", "disabled.csv" if table_path is None else table_path),
        ("MAX_LEVEL", "(-)", 1),
    ]
    return {name: {"unit": unit, "value": value} for name, unit, value in entries}


def plan_generation(config: dict, fluid_root: Path, radial_paths: dict[str, Path],
                    build_path: Path | None = None, *, hash_inputs: bool = False) -> dict:
    """Plan inclusive configured labels; no lexicographic listing or hidden slice.

    ``build_path`` is an opaque caller-selected build location, recorded only.
    Optional hashing streams each unique source/table once and does NOT verify
    body semantics or build bytes. Execution still has to hash live inputs.
    """
    validate_production_config(config)
    if "producer" not in config:
        raise ValueError("Generation planning requires producer section")
    if type(hash_inputs) is not bool:
        raise ValueError("hash_inputs must be a boolean")
    if not isinstance(radial_paths, dict) or set(radial_paths) != {"low", "high"}:
        raise ValueError("radial_paths requires exactly low and high paths")
    config = deepcopy(config)
    p = config["producer"]
    config_sha = _digest(config)
    fluid_root = Path(fluid_root).resolve()
    build = None if build_path is None else str(_path(build_path, "build_path"))
    hashes = {}

    def input_hash(path):
        if not hash_inputs:
            return None
        key = str(path)
        if key not in hashes:
            hashes[key] = file_sha256(path)
        return hashes[key]

    snapshots = []
    start, stop, step = (p["snapshots"][key] for key in ("start", "stop", "step"))
    for label in range(start, stop + 1, step):
        name = f"data{label:06d}.dat"
        path = _path(fluid_root / name, "source_path")
        metadata = read_snapshot_footer(path, p["grid"])
        if snapshots:
            previous = snapshots[-1]
            expected = config["physical"]["label_step_tg"] * (label - previous["label"])
            if not math.isclose(metadata["time_tg"] - previous["time_tg"], expected,
                                rel_tol=0, abs_tol=1e-8):
                raise ValueError(f"BHAC footer time cadence disagrees with labels: {name}")
            if metadata["gamma"] != previous["metadata"]["gamma"]:
                raise ValueError(f"BHAC gamma changes across snapshots: {name}")
        snapshots.append({
            "label": label, "snapshot_name": name, "time_tg": metadata["time_tg"],
            "source_path": str(path), "source_sha256": input_hash(path), "metadata": metadata,
        })

    tables, table_cache = {}, {}
    for table_id in ("low", "high"):
        path = _path(radial_paths[table_id], f"radial_paths.{table_id}")
        required = [alpha for alpha in config["data"]["alphas"]
                    if p["alpha_table"][canonical_alpha(alpha)] == table_id]
        # Unique physical files are loaded once, even if explicitly mapped twice.
        key = str(path)
        if key not in table_cache:
            table_cache[key] = validate_radial_table(path, [])
        metadata = deepcopy(table_cache[key])
        if any(not metadata["alpha_min"] <= alpha <= metadata["alpha_max"] for alpha in required):
            raise ValueError(f"Radial CSV {table_id} does not cover required alphas {required}: {path}")
        tables[table_id] = {
            "path": key, "spin": p["radial_tables"][table_id]["spin"],
            "sha256": input_hash(path), "metadata": metadata,
        }

    tasks = []

    def append_task(snapshot, resolution, alpha=None, phase=0):
        enabled = alpha is not None
        table_id = p["alpha_table"][canonical_alpha(alpha)] if enabled else None
        table = tables[table_id] if enabled else None
        omega = float(model.omega(alpha, config["inference"]["principal_n"])) if enabled else 0.0
        norm = p["response"]["reference_norm"] if enabled else 0.0
        kind = "response" if enabled else "baseline"
        identity = {
            "schema_version": SCHEMA_VERSION, "config_sha256": config_sha,
            "snapshot_name": snapshot["snapshot_name"], "time_tg": snapshot["time_tg"],
            "kind": kind, "resolution": resolution, "alpha": alpha,
            "omega": omega, "norm": norm, "phase_index": phase, "table_id": table_id,
        }
        table_path = None if table is None else table["path"]
        parameters = _parameters(config, resolution, alpha, omega, norm, phase, table_path)
        if not enabled:
            # A disabled axion model consumes no table, alpha, omega, norm or
            # phase. Analysis/inference settings likewise do not alter its rays.
            identity = {
                "schema_version": SCHEMA_VERSION, "kind": "baseline", "label": snapshot["label"],
                "snapshot_name": snapshot["snapshot_name"], "time_tg": snapshot["time_tg"],
                "spin": snapshot["metadata"]["spin"], "gamma": snapshot["metadata"]["gamma"],
                "parameters": {key: value for key, value in parameters.items()
                               if not key.startswith("AXION_")},
                "grid": p["grid"], "build": p["build"],
            }
        tasks.append({
            "id": _digest(identity), "label": snapshot["label"],
            "snapshot_name": snapshot["snapshot_name"], "time_tg": snapshot["time_tg"],
            "source_path": snapshot["source_path"], "source_sha256": snapshot["source_sha256"],
            "kind": kind, "resolution": resolution, "alpha": alpha, "omega": omega, "norm": norm,
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
        for alpha in config["data"]["alphas"]:
            for phase in range(p["response"]["phase_count"]):
                append_task(snapshot, response_resolution, alpha, phase)

    baseline_count = sum(task["kind"] == "baseline" for task in tasks)
    response_count = len(tasks) - baseline_count
    # Only unavoidable numerical buffers are counted, not HDF5/container,
    # allocator/pointer, AMR-tree, plot or checkpoint overhead.
    stokes_bytes = sum(task["resolution"]**2 * 4 * 8 for task in tasks)
    memory_bytes = max(s["metadata"]["nleafs"] * math.prod(s["metadata"]["nx"]) * (8 + 6) * 8
                       for s in snapshots)
    build_settings = p["build"]
    setup_args = [
        "--code=bhac", "--metric=mks", "--int=rk4", "--grid=SMR", "--rad=pol", "--sf=sfc",
        "--bflip=plus", f"--sigma-cut={build_settings['sigma_cut']}",
        f"--thetae-min={build_settings['thetae_min']}", f"--thetae-max={build_settings['thetae_max']}",
        "--uniform=0",
    ]
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
            "official_setup_args": setup_args, "build_command": ["make"],
            "label_time_offset_tg": snapshots[0]["time_tg"] - start * config["physical"]["label_step_tg"],
        },
        "warnings": [
            "Plan IDs do not verify input content; execution must bind live source/table/build SHA256.",
            "Footer/payload lower-bound checks do not validate fluid values or AMR forest semantics.",
            "Radial CSV spin/provenance cannot be inferred from alpha,r,R columns; "
            "configured spins are declarations.",
            "High table spin=0.99 intentionally differs from actual fluid spin=0.9375 "
            "(explicit migration benchmark).",
            "No radial amplitudes were renormalized; ray-radius coverage remains a core runtime check.",
            "Resource estimates are lower bounds only; "
            "caller must check disk/RAM with overhead and concurrency.",
            "This plan builds and launches nothing and uses no old scientific outputs.",
        ],
    }
