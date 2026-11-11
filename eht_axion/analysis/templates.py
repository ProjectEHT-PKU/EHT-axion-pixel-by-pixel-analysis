"""Simulated images -> response templates (all angles in radians).

preprocess_cube blurs, rotates and resamples the delivered RAPTOR Stokes cubes
onto the shared sky grid; fit_processed fits the first-harmonic EVPA response
over the 16 axion phases per pixel (2D) and per aperture sector (1D);
summarize_fits forms the template mean/std over the snapshots; load_templates
reads them back. RunLayout names the stage directories of scripts/analysis/templates.py.

Stage functions write only to output_dir, return a mapping of product names to
Paths and leave transactions/caching to their caller (provenance.stage_directory).
Input .npy files are always read-only memmaps.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap
from scipy.ndimage import gaussian_filter, map_coordinates

from ..evpa import finite_evpa
from ..raptor.convert import _integer, _metadata, _outputs, _positive, _real_dtype
from ..raptor.layout import alpha_label
from .aperture import aperture_operator, integrate_stokes, sector_brightness

# Numerical constants of the stages (not scientific parameters).
_MINIMUM_SAMPLES = 2  # a standard deviation needs at least two snapshots


class _Series:
    """One (T, ...) time series read through mmap, one image at a time."""

    def __init__(self, path, *, ndim=(4, 5), stokes=True):
        self.path = Path(path)
        array = self._load(self.path)
        _real_dtype(array.dtype)
        if array.ndim not in ndim or array.shape[0] == 0:
            raise ValueError(f"invalid input time-series shape: {array.shape}")
        if stokes and (
            array.shape[-1] != 4 or array.shape[-2] != array.shape[-3]
            or array.shape[-2] < 2
            or (array.ndim == 5 and array.shape[1] == 0)
        ):
            raise ValueError("Stokes input must be (T,[S,]N,N,4), N >= 2")
        self.array = array
        self.shape = array.shape
        self.dtype = array.dtype

    def _load(self, path):
        if path.suffix != ".npy":
            raise ValueError(f"expected a .npy input: {path}")
        return np.load(path, mmap_mode="r", allow_pickle=False)

    def frame(self, index):
        return self.array[index]

    def describe(self):
        return {"input": str(self.path), "shape": list(self.shape), "input_dtype": str(self.dtype)}


def _snapshot_labels(labels, length):
    labels = np.asarray(labels)
    if labels.shape != (length,) or labels.dtype.kind not in "iu":
        raise ValueError("snapshot_labels must be a one-dimensional integer array of length T")
    return labels


class _FitSeries(_Series):
    """Read a named response file (response-*.npz of fit_processed) as (T, spatial..., 3)."""

    def _load(self, path):
        if path.suffix != ".npz":
            raise ValueError(f"expected a response .npz fit input: {path}")
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
            return np.stack(arrays, axis=-1)


def preprocess_frame(
    stokes, *, sigma: float, init_fov: float, final_fov: float,
    rotation_angle: float, final_shape: int,
) -> np.ndarray:
    """Blur then inverse-map a square frame; actively rotate Q+iU by exp(2iθ).

    sigma is a Gaussian standard deviation in the same physical units as FOV,
    not FWHM or pixels. Axis0 is x, axis1 is y. Gaussian blur with nearest-value
    boundary extension, then bilinear resampling. No flux renormalization.
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


def preprocess_cube(
    path, output_dir, *, sigma: float, init_fov: float, final_fov: float,
    rotation_angle: float, final_shape: int, state_indices=None,
) -> dict[str, Path]:
    """Preprocess a delivered (T,S,N,N,4) cube one image at a time into float64 stokes.npy.

    ``state_indices`` selects states (all if None) and retains the state axis.
    """
    series = _Series(path, ndim=(5,))
    # Validate geometry before creating products, using only a tiny dummy frame.
    params = dict(sigma=sigma, init_fov=init_fov, final_fov=final_fov,
                  rotation_angle=rotation_angle, final_shape=final_shape)
    preprocess_frame(np.zeros((2, 2, 4)), **{**params, "final_shape": 2})
    final_shape = _integer(final_shape, "final_shape", 2)
    if state_indices is not None:
        state_indices = [_integer(i, "state index", 0) for i in state_indices]
        if not state_indices or len(set(state_indices)) != len(state_indices):
            raise ValueError("state_indices must be nonempty and unique")
        if max(state_indices) >= series.shape[1]:
            raise ValueError("state_indices is outside the state axis")
    states = list(range(series.shape[1])) if state_indices is None else state_indices
    shape = (series.shape[0], len(states), final_shape, final_shape, 4)
    result = _outputs(output_dir, {"stokes": "stokes.npy", "metadata": "processing.json"},
                      [series.path])
    output = open_memmap(result["stokes"], mode="w+", dtype=np.float64, shape=shape)
    for t in range(series.shape[0]):
        for destination_state, source_state in enumerate(states):
            output[t, destination_state] = preprocess_frame(series.frame(t)[source_state], **params)
    output.flush()
    del output
    _metadata(result["metadata"], stage="preprocess", **series.describe(), **params,
              output_shape=list(shape), dtype="float64", state_indices=state_indices,
              sigma_units="FOV units (stddev)",
              interpolation="bilinear", boundary="nearest", spin2="Q+iU -> exp(2iθ)(Q+iU)")
    return result


def _wrap(angle, period):
    return (angle + period / 2) % period - period / 2


def _evpa(stokes):
    return finite_evpa(stokes[..., 1], stokes[..., 2])


def fit_phase_dft(evpa, *, phase_axis: int = 0) -> np.ndarray:
    """Return [..., amplitude, phi, mean_EVPA] for mean+A*cos(ψ+phi).

    Exactly 16 equally spaced ψ=2πk/16 are required. EVPA has period π;
    center on its circular mean, lift onto that branch, then take the DFT.
    A missing phase invalidates that location. Zero amplitude has undefined phi.
    """
    values = np.moveaxis(np.asarray(evpa, dtype=np.float64), phase_axis, 0)
    if values.shape[0] != 16:
        raise ValueError("expected 16 phases")
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


def _operator(shape, fov, radial_max, window_degrees, angular_bins):
    _positive(fov, "fov")
    angular_bins = _integer(angular_bins, "angular_bins")
    if angular_bins % 4:
        raise ValueError("EHT 90-degree alignment requires angular_bins divisible by 4")
    return aperture_operator(shape=shape, fov=fov, radial_max=radial_max,
                             window_degrees=window_degrees, angular_bins=angular_bins)


def _sky_bins(values):
    """Operator angle zero is +axis0 (west); bin k of the result is centred k degrees east of
    north: output[k] = values[k + bins/4] (axis 1 is north)."""
    if values.shape[-1] % 4:
        raise ValueError("EHT 90-degree alignment requires angular_bins divisible by 4")
    return np.roll(values, -values.shape[-1] // 4, axis=-1)


def _integrate(stokes, operator):
    """Sector EVPA (integrate_stokes) in sky bins."""
    return _sky_bins(integrate_stokes(stokes, operator))


def _brightness(stokes, operator):
    """Sector intensity and I-weighted polarized brightness (sector_brightness) in sky bins."""
    intensity, polarized = sector_brightness(stokes, operator)
    return _sky_bins(intensity), _sky_bins(polarized)


def fit_processed(
    input_dir, output_dir, *, fov: float, radial_max: float, window_degrees: float,
    angular_bins: int, snapshot_labels,
) -> dict[str, Path]:
    """Fit input_dir/stokes.npy, processed (T,17,N,N,4), streaming one time frame.

    State 0 (axion off) is excluded; the 16 phase states are fitted per pixel (2D)
    and per intensity-weighted aperture bin (1D). response-{2d,1d}.npz keys:
    amplitude, phase, baseline (all radians), snapshot_labels (from the caller,
    i.e. from the verified delivered data).
    """
    series = _Series(Path(input_dir) / "stokes.npy", ndim=(5,))
    if series.shape[1] != 17:
        raise ValueError("fit requires state0 + 16 phase states")
    labels = _snapshot_labels(snapshot_labels, series.shape[0])
    n = series.shape[-2]
    operator = _operator(n, fov, radial_max, window_degrees, angular_bins)
    result = _outputs(output_dir, {"fits_2d": "response-2d.npz", "fits_1d": "response-1d.npz",
                                  "metadata": "response.json"}, [series.path])
    two = np.empty((series.shape[0], n, n, 3))
    one = np.empty((series.shape[0], angular_bins, 3))
    for t in range(series.shape[0]):
        # Exclude noaxion before any integration or circular statistics.
        states = series.frame(t)[-16:]
        two[t] = fit_phase_dft(_evpa(states))
        one[t] = fit_phase_dft(_integrate(states, operator))
    for name, data in (("fits_2d", two), ("fits_1d", one)):
        np.savez(result[name], amplitude=data[..., 0], phase=data[..., 1],
                 baseline=data[..., 2], snapshot_labels=labels)
    del two, one
    _metadata(result["metadata"], stage="fit", **series.describe(), dtype="float64", fov=fov,
              weighting="intensity", radial_max=radial_max, window_degrees=window_degrees,
              angular_bins=angular_bins, state0_used=False, phase_count=16,
              phase_convention="mean + amplitude*cos(2*pi*k/16 + phi)", units="rad",
              response_keys=["amplitude", "phase", "baseline", "snapshot_labels"],
              baseline_convention="mean EVPA on fitted circular branch",
              snapshot_labels=labels.tolist(),
              angular_alignment="EHT: np.roll(aperture_evpa, -angular_bins//4, axis=-1)",
              angular_roll_bins=-angular_bins // 4,
              operator_source="eht_axion.analysis.aperture.aperture_operator")
    return result


def _moments_update(count, mean, m2, values):
    valid = np.isfinite(values)
    new_count = count + valid
    delta = np.where(valid, values - mean, 0.0)
    new_mean = mean + delta / np.maximum(new_count, 1)
    new_m2 = m2 + delta * np.where(valid, values - new_mean, 0.0)
    return new_count, new_mean, new_m2


def summarize_fits(fits_by_dimension: Mapping[str, Mapping[float, object]], output_dir) -> dict[str, Path]:
    """Gaussian MLE time distributions; circularly lift phases before fitting.

    ``fits_by_dimension`` is {"1d": {alpha: response-1d.npz}, "2d": {alpha: response-2d.npz}}.
    templates.npz holds "<dimension>/<name>" arrays of shape (alpha, spatial...). Labels do
    not enter statistics. A location needs at least two finite snapshots. No pickle,
    normality gate, or histogram-fit rejection is used.
    """
    if not fits_by_dimension:
        raise ValueError("at least one dimension is required")
    series = {dim: [_FitSeries(fits[a], ndim=(3, 4), stokes=False) for a in sorted(fits)]
              for dim, fits in fits_by_dimension.items()}
    result = _outputs(output_dir, {"templates": "templates.npz", "metadata": "summary.json"},
                      [s.path for dim in series for s in series[dim]])
    arrays, metadata = {}, {}
    for dim, fits in fits_by_dimension.items():
        alphas, products, spatial = _summarize_dimension(fits, series[dim])
        arrays.update({f"{dim}/{name}": value for name, value in products.items()})
        metadata[dim] = dict(alphas=alphas, inputs=[s.describe() for s in series[dim]],
                             spatial_shape=list(spatial))
    np.savez(result["templates"], **arrays)
    _metadata(result["metadata"], stage="summarize", dimensions=metadata,
              minimum_samples=_MINIMUM_SAMPLES,
              method="Gaussian MLE, ddof=0; phases lifted about circular mean", units="rad",
              normality_gate=False)
    return result


def _summarize_dimension(fits_by_alpha, series_list):
    minimum_samples = _MINIMUM_SAMPLES
    if not fits_by_alpha:
        raise ValueError("at least one alpha is required")
    alphas = sorted(fits_by_alpha)
    if not np.all(np.isfinite(alphas)):
        raise ValueError("alphas must be finite numeric values")
    spatial = series_list[0].shape[1:-1]
    if any(s.shape[-1] != 3 or s.shape[1:-1] != spatial for s in series_list):
        raise ValueError("fit products must have matching spatial shape and 3 columns")
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
    return alphas, dict(alphas=np.asarray(alphas, dtype=np.float64), **products, **counts), spatial


def load_templates(path, dimension):
    """(alphas, amplitude, phase, valid) of one dimension of templates.npz."""
    with np.load(path, allow_pickle=False) as z:
        alphas = z[f"{dimension}/alphas"].copy()
        keys = ["amplitude_mean", "amplitude_std", "phase_mean", "phase_std"]
        values = [z[f"{dimension}/{k}"].copy() for k in keys]
    if dimension not in {"1d", "2d"} or alphas.ndim != 1 or not len(alphas):
        raise ValueError("Invalid template dimension/alpha grid")
    if not np.isfinite(alphas).all() or np.any(np.diff(alphas) <= 0):
        raise ValueError("Template alpha grid must be finite and increasing")
    expected_ndim = 2 if dimension == "1d" else 3
    if any(v.shape != values[0].shape or v.ndim != expected_ndim for v in values):
        raise ValueError("Template summary shapes differ")
    if values[0].shape[0] != len(alphas):
        raise ValueError("Template alpha axis mismatch")
    amplitude = {float(a): np.stack((values[0][k], values[1][k])) for k, a in enumerate(alphas)}
    phase = {float(a): np.stack((values[2][k], values[3][k])) for k, a in enumerate(alphas)}
    valid = np.ones(values[0].shape[1:], dtype=bool)
    for v in values:
        valid &= np.isfinite(v).all(axis=0)
    valid &= (values[1] >= 0).all(axis=0) & (values[3] >= 0).all(axis=0)
    return alphas, amplitude, phase, valid


@dataclass(frozen=True)
class RunLayout:
    """Stage directories written by scripts/analysis/templates.py (template root)."""

    root: Path

    def __post_init__(self):
        object.__setattr__(self, "root", Path(self.root))

    def processed(self, alpha):
        return self.root / "processed" / alpha_label(alpha)

    def response(self, alpha):
        return self.root / "response" / alpha_label(alpha)

    @property
    def templates(self):
        return self.root / "templates"

    @property
    def background_processed(self):
        return self.root / "background"
