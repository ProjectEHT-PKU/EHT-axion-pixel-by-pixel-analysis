"""EHT images: native-sky FITS reader and per-day summaries of the D-term Monte Carlo images.

I,Q,U,V are never flipped or sign-converted; EVPA is east of north (radians,
modulo pi). No RAPTOR mapping, beam deconvolution or re-centring (the Paper VII
Sec. 5.2 alignment is a zero shift here).

Per day the stage writes the per-draw Stokes maps (by scripts/eht/summarize.py), the per-draw
EVPA maps and the whole-field image-averaged quantities of Paper VII Sec. 5.3 (m_net,
<|m|>, EVPA, beta_2) for comparison with the paper. Sector apertures, means,
covariances and day differences are computed by the analysis.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..evpa import finite_evpa
from ..experiment import EHT_DAYS


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


def evpa(q, u):
    """Finite nonzero Q+iU only, radians in [0,pi); invalid pairs give NaN."""
    return np.mod(finite_evpa(q, u), np.pi)


DAYS = tuple(str(day) for day in EHT_DAYS)  # day keys of the Monte Carlo batch
GLOBAL_COLUMNS = ('total_I', 'm_net', 'm_resolved', 'evpa_rad', 'beta2_real', 'beta2_imag')


def global_diagnostics(image, east, north):
    """Image-averaged quantities of EHT Paper VII Sec. 5.3, as pixel sums over the whole field.

    m_net = |sum P| / sum I (eq. 12); m_resolved = sum |P| / sum I (<|m|>, eq. 13);
    EVPA of sum P; beta_2 = sum P exp(-2i phi) / sum I (eq. 14 with rho_min = 0 and
    rho_max covering the full field of view, so I_ring is the total I, as in the paper).
    P = Q + iU; phi is the position angle east of north about the grid centre. The
    paper centres each image by cross-correlation with the Paper IV consensus image;
    here the centre is the grid centre (zero shift). The images are the blurred ones
    (shared.image.blur_fwhm_uas = 19 uas, the eht-imaging kernel of Paper VII Sec. 5.1).
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


def summarize(samples, east, north, directory):
    """Write evpa.npy and global.npz for one day (pickle-free)."""
    directory = Path(directory)
    chi = np.lib.format.open_memmap(directory / 'evpa.npy', mode='w+', dtype='float64',
                                   shape=samples.shape[:-1])
    diagnostics = np.empty((len(samples), len(GLOBAL_COLUMNS)))
    for k, image in enumerate(samples):
        chi[k] = evpa(image[..., 1], image[..., 2])
        chi[k][image[..., 0] <= 0] = np.nan
        diagnostics[k] = global_diagnostics(image, east, north)
    chi.flush()
    np.savez(directory / 'global.npz', values=diagnostics, columns=np.array(GLOBAL_COLUMNS))
