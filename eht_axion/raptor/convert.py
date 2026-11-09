"""Delivered RAPTOR HDF5 -> float64 Stokes cube (stokes.npy), and small array-file helpers.

convert_hdf5 reads the core's tiled IQUV datasets of a time x state grid of job
outputs. The helpers (dtype/number checks, new output files that never alias an
input, JSON metadata) are shared with the analysis image stages
(analysis/templates.py). Inputs are read-only; outputs are new files.
"""
from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap

from ..provenance import require_new_output_file


def _real_dtype(dtype):
    dtype = np.dtype(dtype)
    if dtype.kind not in "iuf":
        raise ValueError(f"expected real numeric data, got {dtype}")
    return dtype


def _positive(value, name):
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or int(value) != value or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _outputs(output_dir, names, inputs=()):
    paths = {name: Path(output_dir) / filename for name, filename in names.items()}
    sources = [Path(p).resolve() for p in inputs]
    for path in paths.values():
        if path.resolve() in sources or any(
            path.exists() and path.samefile(source) for source in sources
        ):
            raise ValueError(f"output aliases a read-only input: {path}")
        require_new_output_file(path)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    return paths


def _metadata(path, **values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _hdf_sources(sources):
    rows = [list(row) for row in sources]
    if not rows:
        raise ValueError("at least one HDF5 source is required")
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("HDF5 source grid must be rectangular (time, state)")
    return [Path(p) for row in rows for p in row], (len(rows), len(rows[0]))


def convert_hdf5(
    sources, output_dir, *, datasets: Mapping[str, str], image_shape: tuple[int, int],
    block_origins, pixel_shape: tuple[int, int],
) -> dict[str, Path]:
    """Convert the core's tiled IQUV datasets to a float64 (T,S,H,W,4) stokes.npy.

    ``sources`` is a rectangular time x state grid of HDF5 paths. block_origins[b]
    is the (axis0,axis1) pixel origin of source block b; flattened block pixels
    use C order, never an inferred order.
    """
    import h5py  # Optional dependency: only conversion needs h5py.

    if set(datasets) != set("IQUV") or len(set(datasets.values())) != 4:
        raise ValueError("datasets must map I,Q,U,V to four distinct explicit dataset names")
    if len(image_shape) != 2 or len(pixel_shape) != 2:
        raise ValueError("image_shape and pixel_shape must have two dimensions")
    image_shape = tuple(_integer(v, "image_shape") for v in image_shape)
    pixel_shape = tuple(_integer(v, "pixel_shape") for v in pixel_shape)
    paths, leading = _hdf_sources(sources)
    raw = np.asarray(block_origins)
    if raw.ndim != 2 or raw.shape[1] != 2 or raw.dtype.kind not in "iu":
        raise ValueError("block_origins must be an integer (blocks,2) array")
    origins = raw.astype(np.int64)
    coverage = np.zeros(image_shape, dtype=bool)
    for x, y in origins:
        if x < 0 or y < 0 or x + pixel_shape[0] > image_shape[0] or (
            y + pixel_shape[1] > image_shape[1]
        ):
            raise ValueError("block lies outside image_shape")
        region = coverage[x:x + pixel_shape[0], y:y + pixel_shape[1]]
        if region.any():
            raise ValueError("overlapping HDF5 blocks")
        region[:] = True
    if not coverage.all():
        raise ValueError("block mapping does not cover the entire image")

    # Validate every header before creating output; inspect no large data here.
    expected = [(len(origins), int(np.prod(pixel_shape))), (len(origins),) + pixel_shape]
    for path in paths:
        with h5py.File(path, "r") as source:
            for name in datasets.values():
                if name not in source or not isinstance(source[name], h5py.Dataset):
                    raise ValueError(f"missing HDF5 dataset {name!r} in {path}")
                data = source[name]
                _real_dtype(data.dtype)
                if data.shape not in expected:
                    raise ValueError(f"dataset {name!r} shape {data.shape}, expected {expected}")
    result = _outputs(output_dir, {"stokes": "stokes.npy", "metadata": "convert.json"}, paths)
    shape = leading + image_shape + (4,)
    output = open_memmap(result["stokes"], mode="w+", dtype=np.float64, shape=shape)
    frames = output.reshape((len(paths),) + image_shape + (4,))
    for i, path in enumerate(paths):
        with h5py.File(path, "r") as source:
            for channel, label in enumerate("IQUV"):
                data = source[datasets[label]]
                for b, (x, y) in enumerate(origins):
                    frames[i, x:x + pixel_shape[0], y:y + pixel_shape[1], channel] = (
                        data[b].reshape(pixel_shape)
                    )
    output.flush()
    del frames, output
    _metadata(result["metadata"], stage="convert", inputs=[str(p) for p in paths],
              datasets=dict(datasets), shape=list(shape), dtype="float64",
              block_origins=origins.tolist(),
              pixel_shape=pixel_shape, pixel_order="C", stokes_axis=list("IQUV"))
    return result
