"""Native-sky research utilities; no RAPTOR mapping, beam deconvolution or calibration.

I,Q,U,V are never flipped or sign-converted. EVPA is east of north (radians,
modulo pi); sky tick components are (sin(chi), cos(chi)). Registration measures
relative image translation, NOT absolute/official consensus astrometric calibration.
"""

import numpy as np
from scipy.ndimage import gaussian_filter, shift
from scipy.signal import correlate

from ..angles import circular_difference
from ..stokes import finite_evpa


def _image(image):
    a = np.asarray(image, dtype=np.float64)
    if a.ndim != 3 or a.shape[-1] != 4 or 0 in a.shape or not np.isfinite(a).all():
        raise ValueError("image must be finite nonempty [ny,nx,4] I,Q,U,V")
    return a


def load_stokes_fits(path):
    """Read exactly four 2-D I,Q,U,V HDUs with optional Astropy (memmap=False).

    Return dict: stokes [ny,nx,4], header (primary copy), east_uas [nx],
    north_uas [ny], pixel_uas (dy,dx positive), frequency_hz, bunit.
    Primary WCS is required; absent extension WCS inherits it. Present extension
    geometry/units/frequency must agree. CRVAL is retained, not added to offsets.
    Only unrotated native RA---SIN/DEC--SIN axes in degrees are supported.
    """
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdus:
        if len(hdus) != 4:
            raise ValueError("expected exactly four I,Q,U,V HDUs")
        h = hdus[0].header.copy()
        # Official eht-imaging FITS writer stores the absolute source position
        # in OBSRA/OBSDEC rather than CRVAL. Relative offsets need only CRPIX.
        geometry = ("CDELT1", "CDELT2", "CRPIX1", "CRPIX2")
        if any(k not in h or not np.isfinite(h[k]) for k in geometry):
            raise ValueError("missing/nonfinite primary WCS geometry")
        if h["CDELT1"] == 0 or h["CDELT2"] == 0:
            raise ValueError("pixel scales must be nonzero")
        if (h.get("CTYPE1"), h.get("CTYPE2")) != ("RA---SIN", "DEC--SIN"):
            raise ValueError("expected native RA---SIN/DEC--SIN axes")
        frequency_keys = ("FREQ", "RESTFRQ", "RESTFREQ")
        frequencies = [float(h[k]) for k in frequency_keys if k in h]
        if not frequencies or not np.isfinite(frequencies).all() or min(frequencies) <= 0:
            raise ValueError("positive finite FREQ/RESTFRQ/RESTFREQ required (Hz)")
        frequency = frequencies[0]
        if not np.allclose(frequencies, frequency, rtol=1e-12, atol=0):
            raise ValueError("inconsistent frequency metadata")
        bunit = str(h.get("BUNIT", "")).strip()
        if not bunit:
            raise ValueError("BUNIT required; no brightness-unit conversion is performed")
        planes = []
        for channel, hdu in zip('IQUV', hdus):
            p, hdr = np.asarray(hdu.data), hdu.header
            if 'STOKES' in hdr and str(hdr['STOKES']).strip().upper() != channel:
                raise ValueError('HDU Stokes order must be I,Q,U,V')
            if (p.ndim != 2 or 0 in p.shape or hdr.get("NAXIS") != 2
                    or (hdr.get("NAXIS2"), hdr.get("NAXIS1")) != p.shape
                    or (planes and p.shape != planes[0].shape)):
                raise ValueError("HDU shapes/NAXIS must agree and be 2-D")
            for k in geometry + ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "OBSRA", "OBSDEC"):
                if k in hdr and (k not in h or hdr[k] != h[k]):
                    raise ValueError(f"inconsistent HDU {k}")
            for axis in (1, 2):
                if str(hdr.get(f"CUNIT{axis}", "deg")).lower() != "deg":
                    raise ValueError("spatial units must be deg")
                if hdr.get(f"CROTA{axis}", 0) != 0:
                    raise ValueError("rotated axes unsupported")
                for other in (1, 2):
                    if (f"CD{axis}_{other}" in hdr
                            or hdr.get(f"PC{axis}_{other}", int(axis == other)) != int(axis == other)):
                        raise ValueError("CD/rotated PC axes unsupported")
            if str(hdr.get("BUNIT", bunit)).strip() != bunit:
                raise ValueError("inconsistent HDU BUNIT")
            for k in frequency_keys:
                if k in hdr and not np.isclose(hdr[k], frequency, rtol=1e-12, atol=0):
                    raise ValueError("inconsistent HDU frequency")
            planes.append(np.array(p, dtype=np.float64, copy=True))
    image = _image(np.stack(planes, axis=-1))
    east = (np.arange(image.shape[1]) + 1 - h["CRPIX1"]) * h["CDELT1"] * 3.6e9
    north = (np.arange(image.shape[0]) + 1 - h["CRPIX2"]) * h["CDELT2"] * 3.6e9
    return dict(stokes=image, header=h, east_uas=east, north_uas=north,
                pixel_uas=np.abs([h["CDELT2"], h["CDELT1"]]) * 3.6e9,
                frequency_hz=frequency, bunit=bunit)


def blur_stokes_fwhm(image, kernel_uas, pixel):
    """Additional Gaussian FWHM, NOT effective PSF or beam deconvolution.

    pixel is positive uas/pixel, scalar or (dy,dx). An input 5-uas restoration
    is not an effective PSF estimate. Zero padding applies independently to IQUV.
    """
    a = _image(image)
    p = np.broadcast_to(np.asarray(pixel, dtype=float), (2,))
    if not np.isfinite(p).all() or np.any(p <= 0) or not np.isfinite(kernel_uas) or kernel_uas < 0:
        raise ValueError("finite nonnegative FWHM and positive pixel scales required")
    sigma = kernel_uas / (np.sqrt(8 * np.log(2)) * p)
    return gaussian_filter(a, sigma=(*sigma, 0), mode="constant", cval=0)


def shift_stokes(image, offset_yx):
    """Translate every channel identically in pixels; no polarization rotation."""
    a, offset = _image(image), np.asarray(offset_yx, dtype=float)
    if offset.shape != (2,) or not np.isfinite(offset).all():
        raise ValueError("offset_yx must be a finite (dy,dx)")
    return shift(a, (*offset, 0), order=1, mode="constant", cval=0, prefilter=False)


def estimate_translation(imageI, referenceI, maxshiftpx=8):
    """Return integer (dy,dx) to APPLY to imageI, fitting I only, not Q/U.

    Whole-map mean-subtracted cross-correlation / whole-map L2 norms; zero
    outside map, FFT, search +/-maxshiftpx per axis. No overlap renormalization,
    periodic wrapping, subpixel optimization or absolute calibration is implied.
    """
    a, b = np.asarray(imageI, dtype=float), np.asarray(referenceI, dtype=float)
    if a.ndim != 2 or a.shape != b.shape or not a.size or not np.isfinite([a, b]).all():
        raise ValueError("I maps must be finite, nonempty, matching 2-D arrays")
    if isinstance(maxshiftpx, bool) or not isinstance(maxshiftpx, (int, np.integer)) or maxshiftpx < 0:
        raise ValueError("maxshiftpx must be a nonnegative integer")
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        raise ValueError("flat I maps cannot determine translation")
    a, b = a - a.mean(), b - b.mean()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    scores = correlate(b / nb, a / na, mode="full", method="fft")
    y, x = np.indices(scores.shape) - (np.array(a.shape) - 1)[:, None, None]
    scores[(np.abs(y) > maxshiftpx) | (np.abs(x) > maxshiftpx)] = -np.inf
    index = np.unravel_index(np.argmax(scores), scores.shape)
    return int(y[index]), int(x[index])


def evpa(q, u):
    """Finite nonzero Q+iU only, radians in [0,pi); invalid pairs give NaN."""
    return np.mod(finite_evpa(q, u), np.pi)


def evpa_tick_vectors(q, u):
    """Return (dx_east,dy_north) unit directors, not array-index directions."""
    chi = evpa(q, u)
    return np.sin(chi), np.cos(chi)


def complex_stokes_mean(image, weights=None):
    """Return weighted mean Q+iU, never mean angles; skip invalid Q/U pairs.

    Weights default to uniform; supplied weights must match [ny,nx], finite and
    nonnegative. No valid positive-weight samples returns complex NaN.
    """
    a = np.asarray(image, dtype=float)
    if a.ndim != 3 or a.shape[-1] != 4:
        raise ValueError("expected [ny,nx,4]")
    w = np.ones(a.shape[:2]) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != a.shape[:2] or not np.isfinite(w).all() or np.any(w < 0):
        raise ValueError("weights must match image and be finite nonnegative")
    z = a[..., 1] + 1j * a[..., 2]
    valid = np.isfinite(z) & (w > 0)
    return np.sum(w[valid] * z[valid]) / w[valid].sum() if valid.any() else complex(np.nan, np.nan)


def aperture_curves(image, east, north, binangles, window=10, rmax_uas=None):
    """All-radius sectors centered at binangles (degrees E of N), full width window.

    east/north are 1-D uas axes. Equal-area pixels; optional rmax_uas only.
    Return angle_deg, intensity=sum I, polarized_brightness=sum(I*|Q+iU|)/sum I,
    evpa_rad=0.5 arg sum(I*(Q+iU)) modulo pi, count. Empty/nonpositive-I bins
    yield NaN observables; zero coherent polarization has undefined EVPA.
    """
    a = _image(image)
    e, n, angles = (np.asarray(v, dtype=float) for v in (east, north, binangles))
    if e.shape != (a.shape[1],) or n.shape != (a.shape[0],) or angles.ndim != 1:
        raise ValueError("axis/bin shapes incompatible")
    if not all(np.isfinite(v).all() for v in (e, n, angles)) or not 0 < window <= 360:
        raise ValueError("finite axes/angles and 0 < window <= 360 required")
    for axis in (e, n):
        d = np.diff(axis)
        if d.size and (d[0] == 0 or not np.allclose(d, d[0], rtol=1e-10, atol=0)):
            raise ValueError("axes must be uniformly spaced nonzero pixel grids")
    if rmax_uas is not None and (not np.isfinite(rmax_uas) or rmax_uas <= 0):
        raise ValueError("rmax_uas must be positive finite")
    ee, nn = np.meshgrid(e, n)
    phi = np.arctan2(ee, nn)
    radial = np.ones(ee.shape, bool) if rmax_uas is None else np.hypot(ee, nn) <= rmax_uas
    intensity, brightness, chi, count = [], [], [], []
    for angle in angles:
        distance = np.abs(circular_difference(phi, np.deg2rad(angle), period=2*np.pi))
        mask = radial & (distance <= np.deg2rad(window/2))
        i, q, u, _ = a[mask].T
        total, z = i.sum(), np.sum(i * (q + 1j*u))
        intensity.append(total)
        count.append(mask.sum())
        brightness.append(np.sum(i * np.hypot(q, u)) / total if total > 0 else np.nan)
        chi.append(evpa(z.real, z.imag) if total > 0 else np.nan)
    return dict(angle_deg=angles.copy(), intensity=np.array(intensity),
                polarized_brightness=np.array(brightness), evpa_rad=np.array(chi), count=np.array(count))
