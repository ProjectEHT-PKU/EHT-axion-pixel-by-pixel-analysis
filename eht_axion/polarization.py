"""Fast linear operators for radial integration of Stokes images."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .stokes import finite_evpa


def _linear_pixel_weights(coordinate: NDArray, grid_min: float, spacing: float, size: int):
    # RectBivariateSpline.ev clamps out-of-domain coordinates to the boundary.
    position = np.clip((coordinate - grid_min) / spacing, 0, size - 1)
    lower = np.floor(position).astype(int)
    lower = np.clip(lower, 0, size - 2)
    fraction = position - lower
    return lower, lower + 1, 1 - fraction, fraction


@lru_cache(maxsize=8)
def radial_integration_matrix(
    shape: int = 64,
    fov: float = 32.6868,
    angular_bins: int = 360,
    radial_max: float = 16.343,
    radial_step: float = 0.2,
    radial_subsamples: int = 10,
    angular_subsamples: int = 10,
) -> NDArray[np.floating]:
    """Build a matrix mapping flattened images to one-degree radial flux bins.

    The quadrature grid matches the historical notebook but applies it directly
    to Stokes Q and U.  Bilinear weights make repeated integration a matrix
    multiplication rather than thousands of spline constructions.
    """
    spacing = fov / shape
    grid_min = -fov / 2 + spacing / 2
    dr = radial_step / radial_subsamples
    dphi = 2 * np.pi / angular_bins / angular_subsamples
    radii = np.arange(round(radial_max / radial_step) * radial_subsamples + 1) * dr
    matrix = np.zeros((angular_bins, shape * shape), dtype=float)
    angular_window = np.ones(angular_subsamples + 1)
    angular_window[[0, -1]] = 0.5

    for angular_bin in range(angular_bins):
        # Historical bins are centered on integer degrees.
        phi = (
            angular_bin * angular_subsamples
            + np.arange(angular_subsamples + 1)
            - angular_subsamples / 2
        ) * dphi
        x = radii[:, None] * np.cos(phi[None, :])
        y = radii[:, None] * np.sin(phi[None, :])
        x0, x1, wx0, wx1 = _linear_pixel_weights(x, grid_min, spacing, shape)
        y0, y1, wy0, wy1 = _linear_pixel_weights(y, grid_min, spacing, shape)
        quadrature = radii[:, None] * dr * dphi * angular_window[None, :]
        for xi, wx in ((x0, wx0), (x1, wx1)):
            for yi, wy in ((y0, wy0), (y1, wy1)):
                np.add.at(
                    matrix[angular_bin], (xi * shape + yi).ravel(), (quadrature * wx * wy).ravel()
                )
    return matrix


def circular_window(matrix: NDArray, width: int) -> NDArray:
    """Sum adjacent angular bins using an exactly ``width``-degree window."""
    if width <= 0 or width > matrix.shape[0]:
        raise ValueError("window width must be within the angular grid")
    if width == 1:
        return matrix.copy()
    if hasattr(matrix, "tocsr"):
        from scipy.sparse import csr_matrix

        left = width // 2
        indices = np.arange(matrix.shape[0])
        return csr_matrix(
            sum(
                matrix[(indices + offset) % matrix.shape[0]]
                for offset in range(-left, width - left)
            )
        )
    left = width // 2
    offsets = range(-left, width - left)
    return sum(np.roll(matrix, -offset, axis=0) for offset in offsets)


def integrated_evpa(
    stokes_q: ArrayLike,
    stokes_u: ArrayLike,
    *,
    window_degrees: int = 10,
    operator: NDArray | None = None,
) -> NDArray[np.floating]:
    """Radially integrate Stokes flux and return EVPA in radians.

    Input arrays have shape ``(..., height, width)``.  Output has shape
    ``(..., 360)`` and EVPA period pi.
    """
    q = np.asarray(stokes_q, dtype=np.float64)
    u = np.asarray(stokes_u, dtype=np.float64)
    if q.shape != u.shape or q.ndim < 2 or q.shape[-1] != q.shape[-2]:
        raise ValueError("Q and U must have matching square image dimensions")
    one_degree = operator if operator is not None else radial_integration_matrix(shape=q.shape[-1])
    windowed = circular_window(one_degree, window_degrees)
    q_flat = q.reshape((-1, q.shape[-2] * q.shape[-1]))
    u_flat = u.reshape((-1, u.shape[-2] * u.shape[-1]))
    # scipy sparse matrices expose ``tocsr``; using sparse dot is about an
    # order of magnitude cheaper for this operator (roughly 4% nonzero).
    if hasattr(windowed, "tocsr"):
        q_integrated = windowed.dot(q_flat.T).T
        u_integrated = windowed.dot(u_flat.T).T
    else:
        q_integrated = q_flat @ windowed.T
        u_integrated = u_flat @ windowed.T
    evpa = finite_evpa(q_integrated, u_integrated)
    return evpa.reshape(q.shape[:-2] + (one_degree.shape[0],))


@lru_cache(maxsize=8)
def aperture_operator(
    shape: int = 64,
    fov: float = 32.6868,
    radial_max: float | None = None,
    window_degrees: float = 10,
    angular_bins: int = 360,
    radial_resolution: float = 0.02,
    angular_resolution_degrees: float = 0.2,
):
    """Bilinear integration of r dr dphi over exactly centered pizza apertures.

    Angle zero is +x, increasing toward +y. Pixel centers are half a pixel
    inside the FOV, matching ImageS. Edge clamping retains the historical spline
    boundary convention. Both radial endpoints have trapezoidal half-weights;
    the final radius is exact. This new operator does not alter legacy helpers.
    """
    from scipy.sparse import csr_matrix

    radius = fov / 2 if radial_max is None else radial_max
    if (
        shape < 2
        or angular_bins < 1
        or not np.isfinite(fov)
        or fov <= 0
        or not np.isfinite(radius)
        or not 0 < radius <= fov / 2
        or not np.isfinite(window_degrees)
        or not 0 < window_degrees <= 360
        or radial_resolution <= 0
        or angular_resolution_degrees <= 0
    ):
        raise ValueError("Invalid aperture grid or radius")
    nr = max(1, int(np.ceil(radius / radial_resolution)))
    na = max(1, int(np.ceil(window_degrees / angular_resolution_degrees)))
    radii = np.linspace(0, radius, nr + 1)
    dr = radius / nr
    span = np.deg2rad(window_degrees)
    angular_offsets = np.linspace(-span / 2, span / 2, na + 1)
    rw = radii * dr
    rw[[0, -1]] *= 0.5
    aw = np.full(na + 1, span / na)
    aw[[0, -1]] *= 0.5
    weight = rw[:, None] * aw[None, :]
    spacing = fov / shape
    grid_min = -fov / 2 + spacing / 2
    rows = []
    for k in range(angular_bins):
        phi = 2 * np.pi * k / angular_bins + angular_offsets
        x = radii[:, None] * np.cos(phi)
        y = radii[:, None] * np.sin(phi)
        x0, x1, wx0, wx1 = _linear_pixel_weights(x, grid_min, spacing, shape)
        y0, y1, wy0, wy1 = _linear_pixel_weights(y, grid_min, spacing, shape)
        row = np.zeros(shape * shape)
        for xi, wx in ((x0, wx0), (x1, wx1)):
            for yi, wy in ((y0, wy0), (y1, wy1)):
                np.add.at(row, (xi * shape + yi).ravel(), (weight * wx * wy).ravel())
        rows.append(csr_matrix(row[None, :]))
    from scipy.sparse import vstack

    return vstack(rows, format="csr")


def integrate_stokes(stokes, operator, weighting: str = "intensity"):
    """Coherent aperture EVPA in radians, preserving explicit weight semantics.

    intensity: arg integral I(Q+iU)/2, EHT-VII intensity-weighted statistic.
    flux: arg integral (Q+iU)/2, total polarized flux angle.
    Only all-four-channel NaN pixels are spatial padding (zero flux). Other
    nonfinite I/Q/U values are rejected, even in flux mode; unused V values do
    not affect this observable. Products are float64 and must remain finite.
    Zero or nonfinite integrated polarization is undefined (NaN), not zero EVPA.
    """
    data = np.asarray(stokes, dtype=np.float64)
    if data.ndim < 3 or data.shape[-1] != 4 or data.shape[-2] != data.shape[-3]:
        raise ValueError("Stokes must have square trailing (H,W,4) dimensions")
    if len(operator.shape) != 2 or operator.shape[1] != data.shape[-2] ** 2:
        raise ValueError("Invalid operator shape for Stokes image")
    if weighting not in ("intensity", "flux"):
        raise ValueError("weighting must be intensity or flux")
    padding = np.isnan(data).all(axis=-1)
    if np.any(~np.isfinite(data[..., :3]) & ~padding[..., None]):
        raise ValueError("Nonfinite I/Q/U outside all-channel NaN spatial padding")
    q = np.where(padding, 0.0, data[..., 1])
    u = np.where(padding, 0.0, data[..., 2])
    if weighting == "intensity":
        intensity = np.where(padding, 0.0, data[..., 0])
        with np.errstate(over="ignore", invalid="ignore"):
            q = q * intensity
            u = u * intensity
        if not np.isfinite(q).all() or not np.isfinite(u).all():
            raise ValueError("Nonfinite float64 intensity-weighted Stokes product")
    flat_q = q.reshape(-1, operator.shape[1])
    flat_u = u.reshape(-1, operator.shape[1])
    qi = operator.dot(flat_q.T).T
    ui = operator.dot(flat_u.T).T
    angle = finite_evpa(qi, ui)
    return angle.reshape(data.shape[:-3] + (operator.shape[0],))
