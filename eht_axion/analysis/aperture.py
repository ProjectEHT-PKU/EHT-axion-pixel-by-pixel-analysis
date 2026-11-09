"""1D observable: intensity-weighted EVPA in pizza-slice apertures (sparse quadrature operator)."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.typing import NDArray

from ..evpa import finite_evpa


def _linear_pixel_weights(coordinate: NDArray, grid_min: float, spacing: float, size: int):
    # RectBivariateSpline.ev clamps out-of-domain coordinates to the boundary.
    position = np.clip((coordinate - grid_min) / spacing, 0, size - 1)
    lower = np.floor(position).astype(int)
    lower = np.clip(lower, 0, size - 2)
    fraction = position - lower
    return lower, lower + 1, 1 - fraction, fraction


# Quadrature resolution of the aperture integral (numerical, not scientific):
# radial step in the FOV unit and angular step in degrees.
RADIAL_RESOLUTION = 0.02
ANGULAR_RESOLUTION_DEGREES = 0.2


@lru_cache(maxsize=8)
def aperture_operator(*, shape: int, fov: float, radial_max: float, window_degrees: float,
                      angular_bins: int):
    """Bilinear integration of r dr dphi over exactly centered pizza apertures.

    Angle zero is +x, increasing toward +y. Pixel centers are half a pixel
    inside the FOV, matching ImageS. Edge clamping retains the historical spline
    boundary convention. Both radial endpoints have trapezoidal half-weights;
    the final radius is exact.
    """
    from scipy.sparse import csr_matrix

    radius = radial_max
    radial_resolution, angular_resolution_degrees = RADIAL_RESOLUTION, ANGULAR_RESOLUTION_DEGREES
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


def integrate_stokes(stokes, operator):
    """Intensity-weighted coherent aperture EVPA, arg integral I(Q+iU)/2 (EHT VII), radians.

    Only all-four-channel NaN pixels are spatial padding (zero flux). Other
    nonfinite I/Q/U values are rejected; unused V values do
    not affect this observable. Products are float64 and must remain finite.
    Zero or nonfinite integrated polarization is undefined (NaN), not zero EVPA.
    """
    data = np.asarray(stokes, dtype=np.float64)
    if data.ndim < 3 or data.shape[-1] != 4 or data.shape[-2] != data.shape[-3]:
        raise ValueError("Stokes must have square trailing (H,W,4) dimensions")
    if len(operator.shape) != 2 or operator.shape[1] != data.shape[-2] ** 2:
        raise ValueError("Invalid operator shape for Stokes image")
    padding = np.isnan(data).all(axis=-1)
    if np.any(~np.isfinite(data[..., :3]) & ~padding[..., None]):
        raise ValueError("Nonfinite I/Q/U outside all-channel NaN spatial padding")
    q = np.where(padding, 0.0, data[..., 1])
    u = np.where(padding, 0.0, data[..., 2])
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
