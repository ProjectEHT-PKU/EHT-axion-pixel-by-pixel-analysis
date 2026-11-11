"""Simulated images -> expected axion signal (all angles in radians).

For every alpha, snapshot by snapshot: preprocess_frame blurs, rotates and resamples the
16 axion-phase images of the RAPTOR image cube onto the shared sky grid; snapshot_response
fits the first-harmonic EVPA response over the 16 phases per pixel (2D) and per aperture
sector (1D); expected_signal forms its mean/std over the snapshots, the expected signal
amplitude A(x) and phase phi(x) of the paper, per unit RAPTOR normalisation. Nothing in
between is stored. processed_background puts the axion-off background images on the same
grid; load_expected_signal reads one dimension of expected.npz back.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

from ..evpa import finite_evpa
from ..raptor.convert import _integer, _positive, _real_dtype
from .aperture import aperture_operator, integrate_stokes, sector_brightness

# Numerical constants (not scientific parameters).
_MINIMUM_SAMPLES = 2  # a standard deviation needs at least two snapshots


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


def snapshot_response(states, operator):
    """(pixel, sector) first-harmonic fits [..., amplitude, phase, baseline] of the 16 axion
    phase images (16, N, N, 4) of one snapshot on the sky grid; the axion-off image is not used."""
    return fit_phase_dft(_evpa(states)), fit_phase_dft(_integrate(states, operator))


def alpha_responses(cube, preprocess, operator):
    """Per-snapshot responses (T, N, N, 3) and (T, bins, 3) of one RAPTOR image cube (T, 17, M, M, 4).

    One snapshot at a time: states 1..16 (axion phases k pi/8) are put on the sky grid and fitted.
    """
    if cube.ndim != 5 or cube.shape[1] != 17:
        raise ValueError("response requires (T, state0 + 16 phase states, N, N, 4) images")
    pixels, sectors = [], []
    for t in range(cube.shape[0]):
        states = np.stack([preprocess_frame(cube[t, s], **preprocess) for s in range(1, 17)])
        pixel, sector = snapshot_response(states, operator)
        pixels.append(pixel)
        sectors.append(sector)
    return np.array(pixels), np.array(sectors)


def processed_background(cube, preprocess):
    """Axion-off images (T, 1, N, N, 4) of the background cube (T, 1, M, M, 4) on the sky grid."""
    if cube.ndim != 5 or cube.shape[1] != 1:
        raise ValueError("background requires (T, 1, N, N, 4) axion-off images")
    return np.array([[preprocess_frame(cube[t, 0], **preprocess)] for t in range(cube.shape[0])])


def _moments_update(count, mean, m2, values):
    valid = np.isfinite(values)
    new_count = count + valid
    delta = np.where(valid, values - mean, 0.0)
    new_mean = mean + delta / np.maximum(new_count, 1)
    new_m2 = m2 + delta * np.where(valid, values - new_mean, 0.0)
    return new_count, new_mean, new_m2


def expected_signal(responses_by_alpha: Mapping[float, np.ndarray]):
    """Gaussian MLE over the snapshots of {alpha: (T, spatial..., 3) responses}; phases are lifted
    about their circular mean. Returns (alphas, arrays) with arrays of shape (alpha, spatial...):
    alphas, amplitude_mean, amplitude_std, phase_mean, phase_std, amplitude_count, phase_count.
    A location needs at least two finite snapshots; no normality gate or fit rejection."""
    minimum_samples = _MINIMUM_SAMPLES
    if not responses_by_alpha:
        raise ValueError("at least one alpha is required")
    alphas = sorted(responses_by_alpha)
    if not np.all(np.isfinite(alphas)):
        raise ValueError("alphas must be finite numeric values")
    series_list = [np.asarray(responses_by_alpha[a], dtype=np.float64) for a in alphas]
    spatial = series_list[0].shape[1:-1]
    if any(s.ndim < 3 or s.shape[-1] != 3 or s.shape[1:-1] != spatial or not len(s) for s in series_list):
        raise ValueError("responses must have matching (T, spatial..., 3) shapes")
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
            frame = series[t]
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
            delta = _wrap(series[t][..., 1] - center, 2 * np.pi)
            count, mean, m2 = _moments_update(count, mean, m2, delta)
        good = (count >= minimum_samples) & (np.abs(phasor) > 1e-12 * count)
        products["phase_mean"][a] = np.where(good, _wrap(center + mean, 2 * np.pi), np.nan)
        products["phase_std"][a] = np.where(good, np.sqrt(m2 / np.maximum(count, 1)), np.nan)
        counts["phase_count"][a] = count
    return alphas, dict(alphas=np.asarray(alphas, dtype=np.float64), **products, **counts)




def load_expected_signal(path, dimension):
    """(alphas, amplitude, phase, valid) of one dimension ("1d" sectors, "2d" pixels) of expected.npz."""
    with np.load(path, allow_pickle=False) as z:
        alphas = z[f"{dimension}/alphas"].copy()
        keys = ["amplitude_mean", "amplitude_std", "phase_mean", "phase_std"]
        values = [z[f"{dimension}/{k}"].copy() for k in keys]
    if dimension not in {"1d", "2d"} or alphas.ndim != 1 or not len(alphas):
        raise ValueError("Invalid expected-signal dimension/alpha grid")
    if not np.isfinite(alphas).all() or np.any(np.diff(alphas) <= 0):
        raise ValueError("Expected-signal alpha grid must be finite and increasing")
    expected_ndim = 2 if dimension == "1d" else 3
    if any(v.shape != values[0].shape or v.ndim != expected_ndim for v in values):
        raise ValueError("Expected-signal shapes differ")
    if values[0].shape[0] != len(alphas):
        raise ValueError("Expected-signal alpha axis mismatch")
    amplitude = {float(a): np.stack((values[0][k], values[1][k])) for k, a in enumerate(alphas)}
    phase = {float(a): np.stack((values[2][k], values[3][k])) for k, a in enumerate(alphas)}
    valid = np.ones(values[0].shape[1:], dtype=bool)
    for v in values:
        valid &= np.isfinite(v).all(axis=0)
    valid &= (values[1] >= 0).all(axis=0) & (values[3] >= 0).all(axis=0)
    return alphas, amplitude, phase, valid
