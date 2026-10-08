"""Notebook-free, frame-streamed RAPTOR stages (all angles are radians).

Stage functions write only to output_dir, return a mapping of product names to
Paths, and intentionally leave transactions/caching to their caller. Input .npy
files are always read-only memmaps. See docs/raptor-stage.md for conventions.
"""
from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import curve_fit

__all__ = [
    "convert_hdf5", "preprocess_frame", "preprocess_stokes", "preprocess_cube",
    "fit_phase_dft", "fit_stokes", "fit_processed", "summarize_fits",
    "estimate_background", "main",
]


def _real_dtype(dtype, *, output=False):
    dtype = np.dtype(dtype)
    if output:
        if dtype not in (np.dtype("float32"), np.dtype("float64")):
            raise ValueError("output dtype must be float32 or float64")
    elif dtype.kind not in "iuf":
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
        if path.exists():
            raise FileExistsError(f"refusing to overwrite stage output: {path}")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    return paths


def _metadata(path, **values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")


class _Series:
    """Logical time series over explicitly ordered shards; never concatenate."""

    def __init__(self, inputs, *, ndim=(4, 5), stokes=True):
        if isinstance(inputs, np.ndarray):
            self.paths = []
            self.arrays = [inputs]
        else:
            if isinstance(inputs, (str, Path)):
                inputs = [inputs]
            self.paths = [Path(p) for p in inputs]
            if not self.paths:
                raise ValueError("at least one input shard is required")
            if len({p.resolve() for p in self.paths}) != len(self.paths):
                raise ValueError("duplicate input shards")
            self.arrays = []
            for path in self.paths:
                array = self._load(path)
                _real_dtype(array.dtype)
                self.arrays.append(array)
        first = self.arrays[0]
        for array in self.arrays:
            _real_dtype(array.dtype)
            if array.ndim not in ndim or array.shape[0] == 0:
                raise ValueError(f"invalid input time-series shape: {array.shape}")
            if array.shape[1:] != first.shape[1:] or array.dtype != first.dtype:
                raise ValueError("shards must have identical non-time shape and dtype")
            if stokes and (
                array.shape[-1] != 4 or array.shape[-2] != array.shape[-3]
                or array.shape[-2] < 2
                or (array.ndim == 5 and array.shape[1] == 0)
            ):
                raise ValueError("Stokes input must be (T,[S,]N,N,4), N >= 2")
        self.ends = np.cumsum([a.shape[0] for a in self.arrays])
        self.shape = (int(self.ends[-1]),) + first.shape[1:]
        self.dtype = first.dtype

    def _load(self, path):
        if path.suffix != ".npy":
            raise ValueError(f"expected a .npy input: {path}")
        return np.load(path, mmap_mode="r", allow_pickle=False)

    def frame(self, index):
        shard = int(np.searchsorted(self.ends, index, side="right"))
        start = 0 if shard == 0 else int(self.ends[shard - 1])
        return self.arrays[shard][index - start]

    def describe(self):
        return {"inputs": [str(p) for p in self.paths], "shape": list(self.shape),
                "input_dtype": str(self.dtype), "shard_order": "explicit caller order"}


def _snapshot_labels(labels, length):
    if labels is None:
        return 1600 + 10 * np.arange(length, dtype=np.int64)
    labels = np.asarray(labels)
    if labels.shape != (length,) or labels.dtype.kind not in "iu":
        raise ValueError("snapshot_labels must be a one-dimensional integer array of length T")
    return labels


class _FitSeries(_Series):
    """Read named responses or legacy columns with identical frame semantics."""

    def _load(self, path):
        if path.suffix == ".npy":
            return super()._load(path)
        if path.suffix != ".npz":
            raise ValueError(f"expected a .npy or .npz fit input: {path}")
        with np.load(path, allow_pickle=False) as response:
            keys = ("amplitude", "phase", "baseline")
            if not set((*keys, "snapshot_labels")).issubset(response.files):
                raise ValueError(f"response requires amplitude, phase, baseline, snapshot_labels: {path}")
            arrays = [response[key] for key in keys]
            first = arrays[0]
            if first.ndim not in (2, 3) or first.shape[0] == 0:
                raise ValueError(f"invalid response shape: {first.shape}")
            for array in arrays:
                _real_dtype(array.dtype)
                if array.shape != first.shape or array.dtype != first.dtype:
                    raise ValueError("response arrays must have identical shape and dtype")
            _snapshot_labels(response["snapshot_labels"], first.shape[0])
            # Moderate products only; stacking preserves each legacy column exactly.
            return np.stack(arrays, axis=-1)


def _hdf_sources(sources):
    if isinstance(sources, (str, Path)):
        return [Path(sources)], ()
    sources = list(sources)
    if not sources:
        raise ValueError("at least one HDF5 source is required")
    if all(isinstance(p, (str, Path)) for p in sources):
        return [Path(p) for p in sources], (len(sources),)
    rows = [list(row) for row in sources]
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("HDF5 source grid must be rectangular (time, state)")
    return [Path(p) for row in rows for p in row], (len(rows), len(rows[0]))


def convert_hdf5(
    sources, output_dir, *, datasets: Mapping[str, str], image_shape: tuple[int, int],
    block_origins=None, pixel_shape: tuple[int, int] | None = None, dtype="float64",
) -> dict[str, Path]:
    """Convert named IQUV datasets, optionally tiled, to stokes.npy.

    A path gives (H,W,4), a path list (T,H,W,4), a nested time/state path
    list (T,S,H,W,4). block_origins[b] is the (axis0,axis1) pixel origin of
    source block b. Flattened block pixels use C order, never inferred order.
    """
    import h5py  # Optional dependency: only conversion needs h5py.

    if set(datasets) != set("IQUV") or len(set(datasets.values())) != 4:
        raise ValueError("datasets must map I,Q,U,V to four distinct explicit dataset names")
    if len(image_shape) != 2:
        raise ValueError("image_shape must have two dimensions")
    image_shape = tuple(_integer(v, "image_shape") for v in image_shape)
    dtype = _real_dtype(dtype, output=True)
    paths, leading = _hdf_sources(sources)
    origins = None
    if block_origins is not None:
        if pixel_shape is None or len(pixel_shape) != 2:
            raise ValueError("pixel_shape is required with block_origins")
        pixel_shape = tuple(_integer(v, "pixel_shape") for v in pixel_shape)
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
    elif pixel_shape is not None:
        raise ValueError("pixel_shape requires block_origins")

    # Validate every header before creating output; inspect no large data here.
    for path in paths:
        with h5py.File(path, "r") as source:
            for name in datasets.values():
                if name not in source or not isinstance(source[name], h5py.Dataset):
                    raise ValueError(f"missing HDF5 dataset {name!r} in {path}")
                data = source[name]
                _real_dtype(data.dtype)
                expected = [image_shape] if origins is None else [
                    (len(origins), int(np.prod(pixel_shape))), (len(origins),) + pixel_shape,
                ]
                if data.shape not in expected:
                    raise ValueError(f"dataset {name!r} shape {data.shape}, expected {expected}")
    result = _outputs(output_dir, {"stokes": "stokes.npy", "metadata": "convert.json"}, paths)
    shape = leading + image_shape + (4,)
    output = open_memmap(result["stokes"], mode="w+", dtype=dtype, shape=shape)
    frames = output.reshape((len(paths),) + image_shape + (4,))
    for i, path in enumerate(paths):
        with h5py.File(path, "r") as source:
            for channel, label in enumerate("IQUV"):
                data = source[datasets[label]]
                if origins is None:
                    frames[i, ..., channel] = data[...]
                else:
                    for b, (x, y) in enumerate(origins):
                        frames[i, x:x + pixel_shape[0], y:y + pixel_shape[1], channel] = (
                            data[b].reshape(pixel_shape)
                        )
    output.flush()
    del frames, output
    _metadata(result["metadata"], stage="convert", inputs=[str(p) for p in paths],
              datasets=dict(datasets), shape=list(shape), dtype=str(dtype),
              block_origins=None if origins is None else origins.tolist(),
              pixel_shape=pixel_shape, pixel_order="C", stokes_axis=list("IQUV"))
    return result


def preprocess_frame(
    stokes, *, sigma: float, init_fov: float, final_fov: float,
    rotation_angle: float, final_shape: int = 64,
) -> np.ndarray:
    """Blur then inverse-map a square frame; actively rotate Q+iU by exp(2iθ).

    sigma is a Gaussian standard deviation in the same physical units as FOV,
    not FWHM or pixels. Axis0 is x, axis1 is y. Nearest boundary extension and
    bilinear resampling match the old ImageS blur=True path. No flux renorm.
    """
    stokes = np.asarray(stokes)
    _real_dtype(stokes.dtype)
    if stokes.ndim != 3 or stokes.shape[-1] != 4 or (
        stokes.shape[0] != stokes.shape[1] or stokes.shape[0] < 2
    ):
        raise ValueError("frame must have shape (N,N,4), N >= 2")
    for value, name in ((init_fov, "init_fov"), (final_fov, "final_fov")):
        _positive(value, name)
    if not np.isfinite(sigma) or sigma < 0 or not np.isfinite(rotation_angle):
        raise ValueError("sigma must be finite and nonnegative; rotation_angle must be finite")
    final_shape = _integer(final_shape, "final_shape", 2)
    size = stokes.shape[0]
    grid = (np.arange(final_shape) - (final_shape - 1) / 2) * final_fov / final_shape
    x, y = np.meshgrid(grid, grid, indexing="ij")
    c, s = np.cos(rotation_angle), np.sin(rotation_angle)
    coords = np.array([(c * x + s * y) * size / init_fov + (size - 1) / 2,
                       (-s * x + c * y) * size / init_fov + (size - 1) / 2])
    output = np.empty((final_shape, final_shape, 4), dtype=np.float64)
    for channel in range(4):
        blurred = gaussian_filter(np.asarray(stokes[..., channel], dtype=np.float64),
                                  sigma * size / init_fov, mode="nearest")
        output[..., channel] = map_coordinates(blurred, coords, order=1, mode="nearest",
                                              prefilter=False)
    q, u = output[..., 1].copy(), output[..., 2].copy()
    output[..., 1] = q * np.cos(2 * rotation_angle) - u * np.sin(2 * rotation_angle)
    output[..., 2] = q * np.sin(2 * rotation_angle) + u * np.cos(2 * rotation_angle)
    return output


def preprocess_stokes(
    inputs, output_dir, *, sigma: float, init_fov: float, final_fov: float,
    rotation_angle: float, final_shape: int = 64, dtype="float64", state_indices=None,
) -> dict[str, Path]:
    """Process one image at a time; state_indices selects states retaining axis1."""
    series = _Series(inputs)
    dtype = _real_dtype(dtype, output=True)
    # Validate geometry before creating products, using only a tiny dummy frame.
    params = dict(sigma=sigma, init_fov=init_fov, final_fov=final_fov,
                  rotation_angle=rotation_angle, final_shape=final_shape)
    preprocess_frame(np.zeros((2, 2, 4)), **{**params, "final_shape": 2})
    final_shape = _integer(final_shape, "final_shape", 2)
    if state_indices is not None:
        if len(series.shape) != 5:
            raise ValueError("state_indices requires a 5D input with a state axis")
        state_indices = [_integer(i, "state index", 0) for i in state_indices]
        if not state_indices or len(set(state_indices)) != len(state_indices):
            raise ValueError("state_indices must be nonempty and unique")
        if max(state_indices) >= series.shape[1]:
            raise ValueError("state_indices is outside the state axis")
    states = (list(range(series.shape[1])) if state_indices is None else state_indices
              ) if len(series.shape) == 5 else None
    leading = (series.shape[0], len(states)) if states is not None else (series.shape[0],)
    shape = leading + (final_shape, final_shape, 4)
    result = _outputs(output_dir, {"stokes": "stokes.npy", "metadata": "processing.json"},
                      series.paths)
    output = open_memmap(result["stokes"], mode="w+", dtype=dtype, shape=shape)
    for t in range(series.shape[0]):
        if len(series.shape) == 4:
            output[t] = preprocess_frame(series.frame(t), **params)
        else:
            for destination_state, source_state in enumerate(states):
                output[t, destination_state] = preprocess_frame(series.frame(t)[source_state], **params)
    output.flush()
    del output
    _metadata(result["metadata"], stage="preprocess", **series.describe(), **params,
              output_shape=list(shape), dtype=str(dtype), state_indices=state_indices,
              sigma_units="FOV units (stddev)",
              interpolation="bilinear", boundary="nearest", spin2="Q+iU -> exp(2iθ)(Q+iU)")
    return result


def preprocess_cube(
    sources, output_dir, *, sigma: float, init_fov: float, final_fov: float,
    rotation_angle: float, final_shape: int = 64, dtype="float64", state_indices=None,
) -> dict[str, Path]:
    """Parent-orchestrator entry point; see preprocess_stokes for the contract."""
    return preprocess_stokes(sources, output_dir, sigma=sigma, init_fov=init_fov,
                             final_fov=final_fov, rotation_angle=rotation_angle,
                             final_shape=final_shape, dtype=dtype, state_indices=state_indices)


def _wrap(angle, period):
    return (angle + period / 2) % period - period / 2


def _evpa(stokes):
    q, u = stokes[..., 1], stokes[..., 2]
    return np.where((q == 0) & (u == 0), np.nan, 0.5 * np.arctan2(u, q))


def fit_phase_dft(evpa, *, phase_axis: int = 0) -> np.ndarray:
    """Return [..., amplitude, phi, mean_EVPA] for mean+A*cos(ψ+phi).

    Exactly 16 equally spaced ψ=2πk/16 are used. A 17-state input has state0
    discarded, including for the circular branch reference. EVPA has period π;
    center on its circular mean, lift onto that branch, then take the DFT.
    A missing phase invalidates that location. Zero amplitude has undefined phi.
    """
    values = np.moveaxis(np.asarray(evpa, dtype=np.float64), phase_axis, 0)
    if values.shape[0] == 17:
        values = values[1:]
    if values.shape[0] != 16:
        raise ValueError("expected 16 phases or noaxion state0 + 16 phases")
    finite = np.all(np.isfinite(values), axis=0)
    vector = np.mean(np.exp(2j * values), axis=0)
    center = np.angle(vector) / 2
    valid = finite & (np.abs(vector) > 1e-12)
    lifted = center + _wrap(values - center, np.pi)
    mean = np.mean(lifted, axis=0)
    delta = lifted - mean
    weights = np.exp(-2j * np.pi * np.arange(16) / 16).reshape(
        (16,) + (1,) * (values.ndim - 1)
    )
    coefficient = np.sum(delta * weights, axis=0) / 8
    amplitude = np.abs(coefficient)
    # Numerical cancellation at constant EVPA must not produce arbitrary phases.
    zero = amplitude <= 64 * np.finfo(float).eps
    amplitude = np.where(zero, 0.0, amplitude)
    phase = np.where(zero, np.nan, _wrap(np.angle(coefficient), 2 * np.pi))
    output = np.stack((amplitude, phase, _wrap(mean, np.pi)), axis=-1)
    return np.where(valid[..., None], output, np.nan)


def _operator(shape, fov, operator, radial_max, window_degrees, angular_bins):
    _positive(fov, "fov")
    angular_bins = _integer(angular_bins, "angular_bins")
    if angular_bins % 4:
        raise ValueError("EHT 90-degree alignment requires angular_bins divisible by 4")
    if operator is None:
        from eht_axion.polarization import aperture_operator

        operator = aperture_operator(shape=shape, fov=fov, radial_max=radial_max,
                                     window_degrees=window_degrees, angular_bins=angular_bins)
    if operator.shape != (angular_bins, shape * shape):
        raise ValueError(f"operator must have shape ({angular_bins}, {shape * shape})")
    return operator


def _integrate(stokes, operator, weighting):
    # No fallback to the old flux-only operator: that would silently change science.
    from eht_axion.polarization import integrate_stokes

    evpa = integrate_stokes(stokes, operator, weighting=weighting)
    # RAPTOR aperture zero is +x. EHT bin zero corresponds to RAPTOR +y:
    # output[k] = aperture[k+90] on the standard 360-bin grid.
    shift = -evpa.shape[-1] // 4
    if evpa.shape[-1] % 4:
        raise ValueError("EHT 90-degree alignment requires angular_bins divisible by 4")
    return np.roll(evpa, shift, axis=-1)


def fit_stokes(
    inputs, output_dir, *, fov: float, operator=None, weighting="intensity",
    radial_max=None, window_degrees=10, angular_bins=360, dtype="float64",
    snapshot_labels=None,
) -> dict[str, Path]:
    """Fit processed (T,16|17,N,N,4) Stokes, streaming one time frame.

    response-{2d,1d}.npz keys: amplitude, phase, baseline (all radians),
    snapshot_labels. Labels are caller-provided or the unverified convention
    1600+10*i, even for short synthetic inputs; they are not inferred from data.
    """
    series = _Series(inputs, ndim=(5,))
    if series.shape[1] not in (16, 17):
        raise ValueError("fit requires 16 phase states or state0 + 16 phases")
    if weighting not in ("intensity", "flux"):
        raise ValueError("weighting must be intensity or flux")
    dtype = _real_dtype(dtype, output=True)
    labels = _snapshot_labels(snapshot_labels, series.shape[0])
    n = series.shape[-2]
    operator = _operator(n, fov, operator, radial_max, window_degrees, angular_bins)
    result = _outputs(output_dir, {"fits_2d": "response-2d.npz", "fits_1d": "response-1d.npz",
                                  "metadata": "response.json"}, series.paths)
    two = np.empty((series.shape[0], n, n, 3), dtype=dtype)
    one = np.empty((series.shape[0], angular_bins, 3), dtype=dtype)
    for t in range(series.shape[0]):
        # Exclude noaxion before any integration or circular statistics.
        states = series.frame(t)[-16:]
        two[t] = fit_phase_dft(_evpa(states))
        one[t] = fit_phase_dft(_integrate(states, operator, weighting))
    for name, data in (("fits_2d", two), ("fits_1d", one)):
        np.savez(result[name], amplitude=data[..., 0], phase=data[..., 1],
                 baseline=data[..., 2], snapshot_labels=labels)
    del two, one
    _metadata(result["metadata"], stage="fit", **series.describe(), dtype=str(dtype), fov=fov,
              weighting=weighting, radial_max=radial_max, window_degrees=window_degrees,
              angular_bins=angular_bins, state0_used=False, phase_count=16,
              phase_convention="mean + amplitude*cos(2*pi*k/16 + phi)", units="rad",
              response_keys=["amplitude", "phase", "baseline", "snapshot_labels"],
              baseline_convention="mean EVPA on fitted circular branch",
              snapshot_labels=labels.tolist(),
              snapshot_labels_source="default_convention" if snapshot_labels is None else "caller",
              snapshot_labels_verified=False,
              angular_alignment="EHT: np.roll(aperture_evpa, -angular_bins//4, axis=-1)",
              angular_roll_bins=-angular_bins // 4,
              operator_source="caller or polarization.aperture_operator")
    return result


def fit_processed(
    input_dir, output_dir, *, fov: float, operator=None, weighting="intensity",
    radial_max=None, window_degrees=10, angular_bins=360, dtype="float64",
    snapshot_labels=None,
) -> dict[str, Path]:
    """Fit input_dir/stokes.npy; parent owns transactions and cache handling."""
    return fit_stokes(Path(input_dir) / "stokes.npy", output_dir, fov=fov, operator=operator,
                      weighting=weighting, radial_max=radial_max, window_degrees=window_degrees,
                      angular_bins=angular_bins, dtype=dtype, snapshot_labels=snapshot_labels)


def _moments_update(count, mean, m2, values):
    valid = np.isfinite(values)
    new_count = count + valid
    delta = np.where(valid, values - mean, 0.0)
    new_mean = mean + delta / np.maximum(new_count, 1)
    new_m2 = m2 + delta * np.where(valid, values - new_mean, 0.0)
    return new_count, new_mean, new_m2


def summarize_fits(
    fits_by_alpha: Mapping[float, object], output_dir, *, minimum_samples=2,
) -> dict[str, Path]:
    """Gaussian MLE time distributions; circularly lift phases before fitting.

    Supply one dimension at a time: {alpha: response_path}, legacy 3-column
    NPY/arrays, or ordered mixed-format shards. Labels do not enter statistics.
    templates.npz arrays have shape (alpha, spatial...).
    No pickle, normality gate, or histogram-fit rejection is used.
    """
    minimum_samples = _integer(minimum_samples, "minimum_samples")
    if not fits_by_alpha:
        raise ValueError("at least one alpha is required")
    alphas = sorted(fits_by_alpha)
    if not np.all(np.isfinite(alphas)):
        raise ValueError("alphas must be finite numeric values")
    series_list = [_FitSeries(fits_by_alpha[a], ndim=(3, 4), stokes=False) for a in alphas]
    spatial = series_list[0].shape[1:-1]
    if any(s.shape[-1] != 3 or s.shape[1:-1] != spatial for s in series_list):
        raise ValueError("fit products must have matching spatial shape and 3 columns")
    result = _outputs(output_dir, {"templates": "templates.npz", "metadata": "summary.json"},
                      [p for s in series_list for p in s.paths])
    products = {name: np.full((len(alphas),) + spatial, np.nan) for name in (
        "amplitude_mean", "amplitude_std", "phase_mean", "phase_std",
    )}
    counts = {name: np.zeros((len(alphas),) + spatial, dtype=np.int64)
              for name in ("amplitude_count", "phase_count")}
    for a, series in enumerate(series_list):
        count = np.zeros(spatial, dtype=np.int64)
        mean, m2 = np.zeros(spatial), np.zeros(spatial)
        phasor = np.zeros(spatial, dtype=complex)
        for t in range(series.shape[0]):
            frame = series.frame(t)
            count, mean, m2 = _moments_update(count, mean, m2, frame[..., 0])
            phase = frame[..., 1]
            phasor += np.where(np.isfinite(phase), np.exp(1j * phase), 0)
        good = count >= minimum_samples
        products["amplitude_mean"][a] = np.where(good, mean, np.nan)
        products["amplitude_std"][a] = np.where(good, np.sqrt(m2 / np.maximum(count, 1)), np.nan)
        counts["amplitude_count"][a] = count
        center = np.angle(phasor)
        count = np.zeros(spatial, dtype=np.int64)
        mean, m2 = np.zeros(spatial), np.zeros(spatial)
        for t in range(series.shape[0]):
            delta = _wrap(series.frame(t)[..., 1] - center, 2 * np.pi)
            count, mean, m2 = _moments_update(count, mean, m2, delta)
        good = (count >= minimum_samples) & (np.abs(phasor) > 1e-12 * count)
        products["phase_mean"][a] = np.where(good, _wrap(center + mean, 2 * np.pi), np.nan)
        products["phase_std"][a] = np.where(good, np.sqrt(m2 / np.maximum(count, 1)), np.nan)
        counts["phase_count"][a] = count
    np.savez(result["templates"], alphas=np.asarray(alphas, dtype=np.float64), **products, **counts)
    _metadata(result["metadata"], stage="summarize", alphas=alphas,
              inputs=[s.describe() for s in series_list], minimum_samples=minimum_samples,
              method="Gaussian MLE, ddof=0; phases lifted about circular mean", units="rad",
              normality_gate=False, spatial_shape=list(spatial))
    return result


def estimate_background(
    inputs, output_dir, *, lag=3, cadence_hours=26.68 / 3, noaxion_state=0,
    mask=None, difference="angular", method="mle", bins=100, minimum_samples=2,
    domain="2d", fov=None, operator=None, weighting="intensity", radial_max=None,
    window_degrees=10, angular_bins=360, legacy_threshold=0.2,
) -> dict[str, Path]:
    """Estimate pooled noaxion EVPA(t+lag)-EVPA(t) Gaussian scatter in radians.

    mask is a static boolean (N,N) selection (True = included). domain='1d'
    integrates the masked Stokes first. For (T,N,N,4) inputs set noaxion_state=None.
    Zero differences remain valid; zero polarization is undefined. Default lag3
    is 26.68 hours, not an exact day. No sqrt(2) single-epoch correction is made.
    """
    series = _Series(inputs)
    lag = _integer(lag, "lag")
    bins = _integer(bins, "bins", 2)
    minimum_samples = _integer(minimum_samples, "minimum_samples")
    _positive(cadence_hours, "cadence_hours")
    if series.shape[0] <= lag:
        raise ValueError("background needs more time frames than lag")
    if len(series.shape) == 4:
        if noaxion_state is not None:
            raise ValueError("4D noaxion input requires noaxion_state=None")
    else:
        if noaxion_state is None:
            raise ValueError("5D input requires an explicit noaxion_state index")
        noaxion_state = _integer(noaxion_state, "noaxion_state", 0)
        if noaxion_state >= series.shape[1]:
            raise ValueError("noaxion_state is outside the state axis")
    if difference not in ("angular", "linear") or method not in ("mle", "legacy_hist"):
        raise ValueError("difference must be angular|linear and method mle|legacy_hist")
    if domain not in ("2d", "1d") or weighting not in ("intensity", "flux"):
        raise ValueError("domain must be 2d|1d and weighting intensity|flux")
    n = series.shape[-2]
    mask_paths = []
    if isinstance(mask, (str, Path)):
        mask_paths = [Path(mask)]
        mask = np.load(mask, mmap_mode="r", allow_pickle=False)
    if mask is not None:
        mask = np.asarray(mask)
        if mask.shape != (n, n) or mask.dtype.kind != "b":
            raise ValueError("mask must be a boolean array with shape (N,N)")
    if domain == "1d":
        if fov is None:
            raise ValueError("fov is required for integrated background")
        operator = _operator(n, fov, operator, radial_max, window_degrees, angular_bins)
    if not np.isfinite(legacy_threshold) or legacy_threshold <= 0:
        raise ValueError("legacy_threshold must be finite and positive")

    def angle(t):
        frame = series.frame(t)
        if noaxion_state is not None:
            frame = frame[noaxion_state]
        if domain == "2d":
            evpa = _evpa(frame)
            return evpa if mask is None else evpa[mask]
        if mask is not None:
            frame = np.where(mask[..., None], frame, 0)
        return _integrate(frame, operator, weighting)

    def chunks():
        for t in range(series.shape[0] - lag):
            delta = angle(t + lag) - angle(t)
            if difference == "angular":
                delta = _wrap(delta, np.pi)
            yield delta[np.isfinite(delta)].ravel()

    # Merge per-frame moments, keeping memory bounded by a single image.
    count, mean, m2 = 0, 0.0, 0.0
    lower, upper = np.inf, -np.inf
    for values in chunks():
        if not values.size:
            continue
        k = values.size
        chunk_mean = float(np.mean(values))
        delta = chunk_mean - mean
        m2 += float(np.sum((values - chunk_mean) ** 2)) + delta**2 * count * k / (count + k)
        mean += delta * k / (count + k)
        count += k
        lower, upper = min(lower, float(values.min())), max(upper, float(values.max()))
    if count < minimum_samples:
        raise ValueError(f"only {count} finite background differences; need {minimum_samples}")
    sigma_mle = float(np.sqrt(max(m2, 0) / count))
    # The main histogram includes all finite differences, even outliers.
    edges = np.histogram_bin_edges(np.array([lower, upper]), bins=bins)
    hist = np.zeros(bins, dtype=np.int64)
    for values in chunks():
        hist += np.histogram(values, bins=edges)[0]
    sigma, fitted_mean = sigma_mle, float(mean)
    extra = {}
    if method == "legacy_hist" and sigma_mle > 0:
        if count < max(40, minimum_samples):
            raise ValueError("legacy_hist requires at least 40 finite differences")
        fit_edges = np.linspace(mean - 5 * sigma_mle, mean + 5 * sigma_mle, bins + 1)
        fit_counts = np.zeros(bins, dtype=np.int64)
        for values in chunks():
            fit_counts += np.histogram(values, bins=fit_edges)[0]
        density = fit_counts / (fit_counts.sum() * np.diff(fit_edges))
        centers = (fit_edges[:-1] + fit_edges[1:]) / 2

        def normal_pdf(x, mu, sd):
            return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (np.sqrt(2 * np.pi) * sd)

        params, covariance = curve_fit(normal_pdf, centers, density, p0=[mean, sigma_mle],
                                       bounds=([-np.inf, np.finfo(float).tiny], [np.inf, np.inf]),
                                       maxfev=50000)
        fitted_mean, sigma = map(float, params)
        error = float(np.sqrt(abs(covariance[1, 1])))
        if not np.isfinite(error) or error > sigma * legacy_threshold:
            raise ValueError("legacy histogram fit failed sigma uncertainty threshold")
        extra = {"fit_hist_counts": fit_counts, "fit_hist_edges": fit_edges,
                 "fit_sigma_error": np.asarray(error)}
    result = _outputs(output_dir, {"metadata": "summary.json",
                                  "diagnostics": "evpa-difference.npz"},
                      series.paths + mask_paths)
    np.savez(result["diagnostics"], hist_counts=hist, hist_edges=edges,
             count=np.asarray(count), mean=np.asarray(fitted_mean), sigma=np.asarray(sigma),
             mean_mle=np.asarray(mean), sigma_mle=np.asarray(sigma_mle), **extra)
    _metadata(result["metadata"], stage="background", **series.describe(), sigma=sigma,
              mean=fitted_mean, sigma_mle=sigma_mle, mean_mle=float(mean), finite_count=int(count),
              units="rad", lag=lag, cadence_hours=cadence_hours,
              lag_hours=lag * cadence_hours, exact_day=False, noaxion_state=noaxion_state,
              difference=difference, method=method, domain=domain, weighting=weighting,
              fov=fov, radial_max=radial_max, window_degrees=window_degrees,
              angular_bins=angular_bins, minimum_samples=minimum_samples, histogram_bins=bins,
              mask=None if mask is None else {"included_pixels": int(mask.sum()),
                                             "path": str(mask_paths[0]) if mask_paths else None},
              legacy_threshold=legacy_threshold, zero_differences_retained=True,
              single_epoch_correction=False,
              angular_roll_bins=-angular_bins // 4 if domain == "1d" else None)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Independent CLI; orchestration/transactions may wrap the same functions."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    def io(name):
        sub = commands.add_parser(name)
        sub.add_argument("inputs", nargs="+", type=Path, help="explicitly ordered input paths")
        sub.add_argument("--output-dir", required=True, type=Path)
        return sub

    def aperture(sub, *, required=True):
        sub.add_argument("--fov", required=required, type=float)
        sub.add_argument("--weighting", choices=("intensity", "flux"), default="intensity")
        sub.add_argument("--radial-max", type=float)
        sub.add_argument("--window-degrees", type=float, default=10)
        sub.add_argument("--angular-bins", type=int, default=360)

    convert = io("convert")
    convert.add_argument("--mapping", required=True, type=Path,
                         help="JSON with datasets, image_shape, block_origins/pixel_shape")
    convert.add_argument("--states-per-time", type=int,
                         help="reshape ordered HDF5 paths into a time/state grid")
    convert.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    preprocess = io("preprocess")
    for name in ("sigma", "init-fov", "final-fov", "rotation-angle"):
        preprocess.add_argument("--" + name, required=True, type=float)
    preprocess.add_argument("--final-shape", type=int, default=64)
    preprocess.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    preprocess.add_argument("--state-indices", nargs="+", type=int,
                            help="process only these states, retaining the state axis")
    fit = io("fit")
    aperture(fit)
    fit.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    summarize = commands.add_parser("summarize")
    summarize.add_argument("--fit", action="append", required=True, metavar="ALPHA=PATH",
                           help="repeat alpha to supply ordered fit shards")
    summarize.add_argument("--output-dir", required=True, type=Path)
    summarize.add_argument("--minimum-samples", type=int, default=2)
    background = io("background")
    background.add_argument("--lag", type=int, default=3)
    background.add_argument("--cadence-hours", type=float, default=26.68 / 3)
    background.add_argument("--noaxion-state", default="0", help="state index, or 'none' for 4D")
    background.add_argument("--mask", type=Path)
    background.add_argument("--difference", choices=("angular", "linear"), default="angular")
    background.add_argument("--method", choices=("mle", "legacy_hist"), default="mle")
    background.add_argument("--bins", type=int, default=100)
    background.add_argument("--minimum-samples", type=int, default=2)
    background.add_argument("--domain", choices=("2d", "1d"), default="2d")
    # FOV is only needed for 1D background, not pixel differences.
    aperture(background, required=False)
    background.add_argument("--legacy-threshold", type=float, default=0.2)
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    try:
        if command == "convert":
            mapping = json.loads(args.pop("mapping").read_text(encoding="utf-8"))
            states = args.pop("states_per_time")
            if states is not None:
                states = _integer(states, "states_per_time")
                paths = args["inputs"]
                if len(paths) % states:
                    raise ValueError("HDF5 input count must be divisible by states-per-time")
                args["inputs"] = [paths[i:i + states] for i in range(0, len(paths), states)]
            result = convert_hdf5(args.pop("inputs"), **args, **mapping)
        elif command == "summarize":
            grouped = {}
            for spec in args.pop("fit"):
                alpha, path = spec.split("=", 1)
                grouped.setdefault(float(alpha), []).append(Path(path))
            result = summarize_fits(grouped, **args)
        elif command == "background":
            state = args["noaxion_state"]
            args["noaxion_state"] = None if state.lower() == "none" else int(state)
            result = estimate_background(**args)
        else:
            result = {"preprocess": preprocess_stokes, "fit": fit_stokes}[command](**args)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.error(str(error))
    print(json.dumps({name: str(path) for name, path in result.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
