"""Per-day summaries of the fixed-D Monte Carlo images in native FITS sky coordinates.

Only what the analysis consumes is written: per-draw EVPA maps, 360-sector
aperture brightness (feature-selection masks) and whole-field diagnostics for a
sanity check against EHT Paper VII. Means, covariances and day differences are
computed by the analysis on its selected features.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .polarimetry import aperture_curves, evpa

DAYS = ('95', '96', '100', '101')
GLOBAL_COLUMNS = ('total_I', 'm_net', 'm_resolved', 'evpa_rad', 'beta2_real', 'beta2_imag')


def global_diagnostics(image, east, north):
    """Whole-field m_net, m_resolved, EVPA and beta_2 (E-of-N azimuth).

    beta_2 = sum[(Q+iU) exp(-2i phi)]/sum I. No annular paper-specific cut.
    phi origin is the grid sky zero (no re-centring, Paper VII Sec. 5.2).
    """
    i, q, u = (image[..., k] for k in range(3))
    p, total = q + 1j*u, i.sum()
    ee, nn = np.meshgrid(east, north)
    phi = np.arctan2(ee, nn)
    net = p.sum()
    if total <= 0:
        return [np.nan] * 6
    beta = np.sum(p * np.exp(-2j*phi)) / total
    return [total, abs(net)/total, np.abs(p).sum()/total,
            float(evpa(net.real, net.imag)), beta.real, beta.imag]


def summarize(samples, east, north, directory, rmax_uas, window):
    """Write evpa.npy, aperture.npz and global.npz for one day (pickle-free)."""
    directory = Path(directory)
    n = len(samples)
    angles = np.arange(360.)
    chi = np.lib.format.open_memmap(directory / 'evpa.npy', mode='w+', dtype='float64',
                                   shape=samples.shape[:-1])
    intensity, brightness = np.empty((n, 360)), np.empty((n, 360))
    diagnostics = np.empty((n, len(GLOBAL_COLUMNS)))
    for k, image in enumerate(samples):
        chi[k] = evpa(image[..., 1], image[..., 2])
        chi[k][image[..., 0] <= 0] = np.nan
        curve = aperture_curves(image, east, north, angles, window=window, rmax_uas=rmax_uas)
        intensity[k], brightness[k] = curve['intensity'], curve['polarized_brightness']
        diagnostics[k] = global_diagnostics(image, east, north)
    chi.flush()
    np.savez(directory / 'aperture.npz', angle_deg=angles, intensity=intensity,
             polarized_brightness=brightness)
    np.savez(directory / 'global.npz', values=diagnostics, columns=np.array(GLOBAL_COLUMNS))
