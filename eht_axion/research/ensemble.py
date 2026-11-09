"""Conditional fixed-D Monte Carlo summaries in native FITS sky coordinates.

No theory-axis adapter, angle roll, thermal bootstrap or joint day posterior.
All covariance is equal-weight population (ddof=0) on a local circular branch.
"""
from __future__ import annotations

import numpy as np

from .polarimetry import aperture_curves, evpa, shift_stokes

DAYS = ('95', '96', '100', '101')
PAIRS = (('96', '95'), ('101', '100'))


def branch_difference(a, b):
    """Signed principal EVPA difference, radians, in [-pi/2,pi/2)."""
    return (np.asarray(a) - np.asarray(b) + np.pi / 2) % np.pi - np.pi / 2


def circular_summary(samples, rmin=.7):
    """Require all draws finite at a location; never silently impute samples.

    std_rad is RMS of local branch residuals about the circular mean. Covariance
    additionally subtracts the arithmetic mean of those residuals (ddof=0).
    """
    a = np.asarray(samples, dtype=float)
    if a.ndim < 2 or a.shape[0] < 2 or not np.isfinite(rmin) or not 0 <= rmin <= 1:
        raise ValueError('need >=2 draws and 0 <= rmin <= 1')
    finite = np.isfinite(a).all(axis=0)
    z = np.mean(np.exp(2j * np.where(np.isfinite(a), a, 0)), axis=0)
    mean = np.mod(np.angle(z) / 2, np.pi)
    resultant = np.where(finite, np.abs(z), np.nan)
    mean = np.where(finite, mean, np.nan)
    residual = branch_difference(a, mean)
    std = np.sqrt(np.mean(residual**2, axis=0))
    return dict(mean_rad=mean, resultant=resultant, std_rad=std,
                valid=finite & (resultant >= rmin))


def selected_covariance(samples, mean_rad, indices=None):
    """Reconstruct covariance of selected flattened spatial/aperture columns.

    Invalid columns remain NaN. No n-by-n spatial matrix unless explicitly asked.
    Use each day's covariance separately, then add for independent day errors.
    """
    a = np.asarray(samples, dtype=float)
    mean = np.asarray(mean_rad, dtype=float)
    if a.ndim < 2 or a.shape[1:] != mean.shape or len(a) < 2:
        raise ValueError('sample/mean shape mismatch or too few draws')
    flat, center = a.reshape(len(a), -1), mean.ravel()
    if indices is not None:
        flat, center = flat[:, indices], center[indices]
    residual = branch_difference(flat, center)
    valid = np.isfinite(residual).all(axis=0)
    out = np.full((len(center), len(center)), np.nan)
    x = residual[:, valid]
    x = x - x.mean(axis=0)
    out[np.ix_(valid, valid)] = x.T @ x / len(a)
    return out


def independent_difference(later, earlier, covariance_later, covariance_earlier):
    """Not a matched-index sample difference: days have independent MC seeds."""
    if np.shape(later['mean_rad']) != np.shape(earlier['mean_rad']):
        raise ValueError('day summaries must have identical geometry')
    return dict(difference_rad=branch_difference(later['mean_rad'], earlier['mean_rad']),
                covariance_rad2=np.asarray(covariance_later) + np.asarray(covariance_earlier),
                valid=later['valid'] & earlier['valid'])


def centroid_offset(mean_i, east, north):
    """I-only fractional (dy,dx) moving the centroid to sky offset (0,0).

    Signed I weights, no positivity clipping or +/-8 pixel search clipping.
    Fail if total I is nonpositive or measured centroid is outside the image.
    """
    i, e, n = np.asarray(mean_i), np.asarray(east), np.asarray(north)
    if i.shape != (len(n), len(e)) or not np.isfinite(i).all() or i.sum() <= 0:
        raise ValueError('finite I with positive total and matching axes required')
    ce, cn = np.sum(i * e[None, :]) / i.sum(), np.sum(i * n[:, None]) / i.sum()
    if not min(e) <= ce <= max(e) or not min(n) <= cn <= max(n):
        raise ValueError('I centroid outside field')
    return np.array([-cn / (n[1] - n[0]), -ce / (e[1] - e[0])])


def translate_ensemble(samples, offset, destination):
    """Same nonperiodic bilinear IQUV translation for every draw."""
    for k, image in enumerate(samples):
        destination[k] = shift_stokes(image, offset)
    destination.flush()


def global_diagnostics(image, east, north):
    """Whole-field m_net, m_resolved, EVPA and beta_2 (E-of-N azimuth).

    beta_2 = sum[(Q+iU) exp(-2i phi)]/sum I. No annular paper-specific cut.
    phi origin is grid sky zero; centered variant is I-centroid referenced.
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


def summarize(samples, east, north, directory, rmin, rmax_uas, window):
    """Write reusable, pickle-free map and 360-sector summaries for one day."""
    from pathlib import Path
    directory = Path(directory)
    n = len(samples)
    angles = np.arange(360.)
    chi = np.lib.format.open_memmap(directory / 'evpa.npy', mode='w+', dtype='float64',
                                   shape=samples.shape[:-1])
    aperture = np.empty((n, 360))
    intensity, brightness = np.empty_like(aperture), np.empty_like(aperture)
    diagnostics = np.empty((n, 6))
    for k, image in enumerate(samples):
        chi[k] = evpa(image[..., 1], image[..., 2])
        chi[k][image[..., 0] <= 0] = np.nan
        curve = aperture_curves(image, east, north, angles, window=window, rmax_uas=rmax_uas)
        aperture[k] = curve['evpa_rad']
        intensity[k], brightness[k] = curve['intensity'], curve['polarized_brightness']
        diagnostics[k] = global_diagnostics(image, east, north)
    chi.flush()
    maps = circular_summary(chi, rmin)
    mean_stokes = np.mean(samples, axis=0)
    maps['valid'] &= (mean_stokes[..., 0] > 0) & (np.hypot(mean_stokes[..., 1], mean_stokes[..., 2]) > 0)
    np.savez(directory / 'map-summary.npz', **maps, mean_stokes=mean_stokes)
    summary = circular_summary(aperture, rmin)
    covariance = selected_covariance(aperture, summary['mean_rad'])
    np.savez(directory / 'aperture.npz', samples_rad=aperture, **summary,
             covariance_rad2=covariance, angle_deg=angles, intensity=intensity,
             polarized_brightness=brightness)
    np.savez(directory / 'global.npz', values=diagnostics,
             columns=np.array(['total_I', 'm_net', 'm_resolved', 'evpa_rad', 'beta2_real', 'beta2_imag']))
    return summary, covariance
