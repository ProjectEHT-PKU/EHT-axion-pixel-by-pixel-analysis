"""Fast linear operators for radial integration of Stokes images."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray


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
    q = np.asarray(stokes_q)
    u = np.asarray(stokes_u)
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
    evpa = 0.5 * np.arctan2(u_integrated, q_integrated)
    return evpa.reshape(q.shape[:-2] + (one_degree.shape[0],))
