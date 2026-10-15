"""Structural adapter for the core's POL=1, one-frequency, uniform camera.

IQUV are Jy/sr (NOT Jy/pixel); axis0=x/alpha, axis1=y/beta. Neither
polarization signs nor coordinate order are repaired here. This validates
serialization/units, not a scientific model or a target source flux.
"""
from __future__ import annotations

from pathlib import Path
import re

import h5py
import numpy as np

# model/bhac/definitions.h in the pinned core, not modern IAU constants.
CORE_CONSTANTS = {"GGRAV": 6.674e-8, "MSUN": 1.989e33,
                  "SPEED_OF_LIGHT": 2.99792458e10, "KPCTOCM": 3.086e21}


def read_core_constants(path: Path) -> dict:
    """Read literal conversion constants from the effective build header.

    Deliberately do not eval C expressions or silently substitute modern values.
    A changed header with nonliteral definitions requires a new adapter.
    """
    text = Path(path).read_text()
    result = {}
    for name in CORE_CONSTANTS:
        matches = re.findall(rf"^\s*#define\s+{name}\s+\(?([+\-\d.eE]+)\)?\s*(?://.*)?$",
                             text, re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f"Cannot identify core constant {name} in {path}")
        value = float(matches[0])
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"Invalid core constant {name}")
        result[name] = value
    return result


def hdf_mapping(resolution: int, frequency: float) -> dict:
    """Explicit kwargs for pipeline.raptor.convert_hdf5; block and pixel x-major."""
    if isinstance(resolution, bool) or not isinstance(resolution, (int, np.integer)):
        raise ValueError("resolution must be a positive integer multiple of 10")
    if resolution <= 0 or resolution % 10:
        raise ValueError("resolution must be a positive integer multiple of 10")
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("frequency must be positive and finite")
    blocks = resolution // 10
    return {"datasets": {c: f"{c}{frequency:e}" for c in "IQUV"},
            "image_shape": [int(resolution), int(resolution)],
            "block_origins": [[10 * x, 10 * y] for x in range(blocks) for y in range(blocks)],
            "pixel_shape": [10, 10]}


def _output_file(directory, relative):
    path = directory / relative
    if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(directory.resolve()):
        raise ValueError(f"Missing/unsafe output file: {path}")
    if any(p.is_symlink() for p in path.parents if p != directory.parent):
        raise ValueError(f"Output path contains a symlink: {path}")
    return path


def validate_job_output(directory: Path, task: dict) -> dict:
    """Validate one successful process's products and return JSON-safe metadata.

    task['parameters'] contains the 19 core values or explicit {unit,value} records.
    output_slot is the output slot, while AXION_PHASE is the physical 0..15
    phase (legacy flat task dictionaries may use phase_index as the slot).
    Optional core_constants must come from a declared effective header;
    direct callers without a build use the documented pinned-core constants.
    """
    directory = Path(directory)
    p = {key: value["value"] if isinstance(value, dict) else value
         for key, value in task["parameters"].items()}
    if p["MAX_LEVEL"] != 1:
        raise ValueError("Only uniform MAX_LEVEL=1 output is supported")
    resolution = task["resolution"]
    mapping = hdf_mapping(resolution, p["FREQ_MIN"])
    if p["IMG_WIDTH"] != resolution or p["IMG_HEIGHT"] != resolution:
        raise ValueError("Task resolution differs from model image dimensions")
    slot = task.get("output_slot", task["phase_index"])
    if isinstance(slot, bool) or not isinstance(slot, int) or not 0 <= slot <= 16:
        raise ValueError("output slot must be an integer in 0..16")
    hdf_relative = f"output/img_data_{slot}.h5"
    spec_relative = f"output/spectrum_{slot}_{p['INCLINATION']:.2f}.dat"
    hdf = _output_file(directory, hdf_relative)
    spectrum = _output_file(directory, spec_relative)
    constants = dict(task.get("core_constants", CORE_CONSTANTS))
    if set(constants) != set(CORE_CONSTANTS) or any(
        not np.isfinite(v) or v <= 0 for v in constants.values()
    ):
        raise ValueError("Invalid core conversion constants")
    for key in ("CAM_SIZE_X", "CAM_SIZE_Y", "MBH", "DISTANCE"):
        if not np.isfinite(p[key]) or p[key] <= 0:
            raise ValueError(f"Invalid geometry parameter {key}")
    rg = constants["GGRAV"] * (p["MBH"] * constants["MSUN"])
    rg = rg / constants["SPEED_OF_LIGHT"] / constants["SPEED_OF_LIGHT"]
    distance = p["DISTANCE"] * constants["KPCTOCM"]
    d_omega = (p["CAM_SIZE_X"] / resolution) * (p["CAM_SIZE_Y"] / resolution)
    d_omega *= (rg / distance) ** 2
    if not np.isfinite(d_omega) or d_omega <= 0:
        raise ValueError("Nonfinite/zero angular pixel area")
    shape = ((resolution // 10) ** 2, 100)
    frequency = p["FREQ_MIN"]
    names = dict(mapping["datasets"], tau=f"tau{frequency:e}", tauF=f"tauF{frequency:e}",
                 alpha="alpha", beta="beta")
    integrated = []
    negative_i = overpolarized = 0
    with h5py.File(hdf, "r") as source:
        if set(source.keys()) != set(names.values()):
            raise ValueError("Expected exactly I/Q/U/V/tau/tauF/alpha/beta at one frequency")
        for name in names.values():
            # External/soft HDF links are not a self-contained output.
            if not isinstance(source.get(name, getlink=True), h5py.HardLink):
                raise ValueError(f"Unsafe HDF5 link: {name}")
            data = source[name]
            if not isinstance(data, h5py.Dataset) or data.shape != shape or data.dtype.kind not in "iuf":
                raise ValueError(f"Invalid HDF5 dataset shape/type: {name}, expected {shape}")
            if data.is_virtual or data.external:
                raise ValueError(f"External HDF5 storage is unsupported: {name}")
        totals = np.zeros(4, dtype=np.float64)
        scale = np.array([p["CAM_SIZE_X"], p["CAM_SIZE_Y"]], dtype=float)
        # Bounded memory: at most one block of eight channels at a time.
        pixel_x, pixel_y = np.arange(100) // 10, np.arange(100) % 10
        for block, (origin_x, origin_y) in enumerate(mapping["block_origins"]):
            arrays = {key: np.asarray(source[name][block], dtype=np.float64)
                      for key, name in names.items()}
            if not all(np.isfinite(a).all() for a in arrays.values()):
                raise ValueError(f"Nonfinite HDF5 channel in block {block}")
            expected_x = -scale[0] / 2 + (origin_x + pixel_x + .5) * scale[0] / resolution
            expected_y = -scale[1] / 2 + (origin_y + pixel_y + .5) * scale[1] / resolution
            for axis, expected, extent in (("alpha", expected_x, scale[0]),
                                            ("beta", expected_y, scale[1])):
                if not np.allclose(arrays[axis], expected, rtol=0,
                                   atol=64 * np.finfo(float).eps * extent):
                    raise ValueError(f"Incorrect {axis} pixel centers/order in block {block}")
            stokes = np.stack([arrays[c] for c in "IQUV"])
            totals += np.sum(stokes, axis=1, dtype=np.float64)
            negative_i += int(np.count_nonzero(stokes[0] < 0))
            overpolarized += int(np.count_nonzero(np.linalg.norm(stokes[1:], axis=0) > stokes[0]))
        integrated = totals * d_omega
    try:
        rows = np.loadtxt(spectrum, ndmin=2)
    except (OSError, ValueError) as exc:
        raise ValueError("Unreadable spectrum") from exc
    if rows.shape != (1, 5) or not np.isfinite(rows).all():
        raise ValueError("Spectrum requires one finite row: frequency I Q U V")
    if not np.isclose(rows[0, 0], frequency, rtol=1e-14, atol=0):
        raise ValueError("Spectrum frequency mismatch")
    # The file uses %.15e, unlike diagnostic stderr (%.5e). No target-flux cut.
    # A scale-aware absolute bound also handles cancellation in Q/U/V.
    tolerance = 1e-10 * max(float(np.max(np.abs(integrated))), np.finfo(float).tiny)
    if not np.allclose(rows[0, 1:], integrated, rtol=1e-10, atol=tolerance):
        raise ValueError("Spectrum flux mismatch with HDF5 Jy/sr times angular pixel area")
    return {"schema": 1, "validated": True, "hdf5": hdf_relative, "spectrum": spec_relative,
            "mapping": mapping, "datasets": names, "hdf_shape": list(shape),
            "axes": ["x/alpha", "y/beta", "Stokes"], "stokes": list("IQUV"),
            "block_order": "x-major; block=x*B+y", "pixel_order": "x-major; pixel=x*10+y",
            "intensity_units": "Jy/sr", "spectrum_units": "Jy", "coordinate_units": "Rg",
            "spectrum_columns": ["frequency_hz", "I_Jy", "Q_Jy", "U_Jy", "V_Jy"],
            "optical_depth_channels": ["tau", "tauF"], "optical_depth_units": "dimensionless",
            "angular_pixel_area_sr": float(d_omega), "core_constants": constants,
            "core_constants_source": task.get("core_constants_source", "pinned BHAC definitions.h"),
            "hdf_integrated_flux_jy": integrated.tolist(), "spectrum_flux_jy": rows[0, 1:].tolist(),
            "frequency_hz": float(frequency), "snapshot_label": task.get("label"),
            "time_tg": task.get("time_tg"), "kind": task.get("kind"), "state_slot": slot,
            "alpha": task.get("alpha"), "omega": task.get("omega"), "norm": task.get("norm"),
            "physical_phase_index": p["AXION_PHASE"],
            "diagnostics": {"negative_I_pixels": negative_i, "polarization_exceeds_I_pixels": overpolarized},
            "scientific_validation": False}
