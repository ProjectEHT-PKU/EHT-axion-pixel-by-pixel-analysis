"""Synthetic-only research polarimetry contract tests (no observational files)."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.ndimage import gaussian_filter

from eht_axion.research.polarimetry import (
    aperture_curves, blur_stokes_fwhm, complex_stokes_mean, estimate_translation,
    evpa, evpa_tick_vectors, load_stokes_fits, shift_stokes,
)


def test_native_fits(tmp_path):
    fits = pytest.importorskip("astropy.io.fits")
    h = fits.Header(dict(CTYPE1="RA---SIN", CTYPE2="DEC--SIN", CDELT1=-2/3.6e9,
                         CDELT2=3/3.6e9, CRPIX1=2.5, CRPIX2=1.5,
                         CRVAL1=187.7, CRVAL2=12.4, FREQ=230e9, BUNIT="Jy/pixel"))
    planes = [np.arange(12.).reshape(3, 4) + 20*k for k in range(4)]
    hdus = fits.HDUList([fits.PrimaryHDU(planes[0], h)] +
                        [fits.ImageHDU(p) for p in planes[1:]])
    path = tmp_path / "synthetic.fits"
    hdus.writeto(path)
    result = load_stokes_fits(path)
    assert_allclose(result["stokes"], np.stack(planes, axis=-1))
    assert result["stokes"].dtype.isnative
    assert_allclose(result["east_uas"], [3, 1, -1, -3])
    assert_allclose(result["north_uas"], [-1.5, 1.5, 4.5])
    assert_allclose(result["pixel_uas"], [3, 2])
    assert result["frequency_hz"] == 230e9
    assert result["bunit"] == "Jy/pixel"
    assert result["header"]["CRVAL1"] == 187.7
    mutations = [lambda: hdus[1].header.update(CDELT1=1),
                 lambda: hdus[1].header.update(CTYPE1="RA---TAN"),
                 lambda: hdus[1].header.update(BUNIT="K"),
                 lambda: hdus[1].header.update(FREQ=220e9),
                 lambda: hdus[1].header.update(PC1_2=0.1),
                 lambda: hdus[1].header.update(CUNIT2="rad"),
                 lambda: setattr(hdus[1], "data", np.zeros((2, 4))),
                 lambda: hdus[1].data.__setitem__((0, 0), np.nan),
                 lambda: hdus[0].header.update(CDELT1=0)]
    original = [hdu.copy() for hdu in hdus]
    for mutate in mutations:
        hdus = fits.HDUList([hdu.copy() for hdu in original])
        mutate()
        hdus.writeto(path, overwrite=True)
        with pytest.raises(ValueError):
            load_stokes_fits(path)
    fits.HDUList(original[:3]).writeto(path, overwrite=True)
    with pytest.raises(ValueError):
        load_stokes_fits(path)


def test_evpa_and_ticks():
    assert_allclose(evpa([1, 0, -1, 0], [0, 1, 0, -1]), np.arange(4)*np.pi/4)
    dx, dy = evpa_tick_vectors([1, 0, -1], [0, 1, 0])
    assert_allclose(dx, [0, np.sqrt(.5), 1])
    assert_allclose(dy, [1, np.sqrt(.5), 0], atol=1e-15)
    assert np.isnan(evpa([0, np.nan, 1], [0, 1, np.inf])).all()
    image = np.zeros((1, 2, 4))
    image[0, :, 1:3] = [[-1, .1], [-1, -.1]]
    assert complex_stokes_mean(image) == -1 + 0j
    image[0, 0, 1] = np.nan
    assert complex_stokes_mean(image) == -1 - .1j


def test_gaussian_fwhm_and_boundaries():
    a = np.zeros((51, 51, 4))
    a[25, 25] = [1, 2, -3, 4]
    fwhm, pixel = 4., (2., 1.)
    out = blur_stokes_fwhm(a, fwhm, pixel)
    sigma = fwhm / (np.sqrt(8*np.log(2))*np.array(pixel))
    for k in range(4):
        assert_allclose(out[..., k], gaussian_filter(a[..., k], sigma, mode="constant", cval=0))
    assert_allclose(out[25, 26, 0]/out[25, 25, 0], np.exp(-1/(2*sigma[1]**2)))
    assert_allclose(blur_stokes_fwhm(a, 0, 1), a)
    a[:] = 1
    assert blur_stokes_fwhm(a, 4, 1)[0, 0, 0] < .5


def test_translation_sign_coherence_and_flat_rejection():
    a = np.zeros((41, 43, 4))
    a[15:20, 17:21] = [1, 2, -3, 4]
    a[23, 26] = [2, 4, -6, 8]
    moved = shift_stokes(a, (3, -4))
    assert estimate_translation(moved[..., 0], a[..., 0], 5) == (-3, 4)
    assert estimate_translation(a[..., 0], a[..., 0]) == (0, 0)
    assert_allclose(shift_stokes(moved, (-3, 4)), a)
    half = shift_stokes(a, (.5, -.5))
    assert_allclose(half[..., 1], 2*half[..., 0])
    assert_allclose(half[..., 2], -3*half[..., 0])
    assert_allclose(half[..., 3], 4*half[..., 0])
    for flat in (np.zeros((4, 4)), np.full((4, 4), .1)):
        with pytest.raises(ValueError):
            estimate_translation(flat, flat)


def test_aperture_pa_weighting_wrap_and_radius():
    # North at x=1,y=2; east at x=0,y=1 (east is displayed left).
    a = np.zeros((3, 3, 4))
    a[2, 1] = [1, 1, 0, 0]
    a[1, 0] = [2, 0, 3, 0]
    result = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [0, 90, 360])
    assert_allclose(result["polarized_brightness"], [1, 3, 1])
    assert_allclose(result["evpa_rad"], [0, np.pi/4, 0])
    assert_allclose(result["intensity"], [1, 2, 1])
    wide = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [45], window=100)
    assert_allclose(wide["polarized_brightness"], [7/3])
    assert_allclose(wide["evpa_rad"], [.5*np.arctan2(6, 1)])
    cut = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [0], rmax_uas=.5)
    assert np.isnan(cut["evpa_rad"]).all()


def test_invalid_inputs():
    a = np.zeros((3, 4, 4))
    for operation in (lambda: blur_stokes_fwhm(a, -1, 1),
                      lambda: blur_stokes_fwhm(a, 1, 0),
                      lambda: shift_stokes(a, (np.nan, 0)),
                      lambda: estimate_translation(np.ones((2, 2)), np.ones((3, 2))),
                      lambda: aperture_curves(a, [0], [0], [0])):
        with pytest.raises(ValueError):
            operation()
    a[0, 0, 2] = np.nan
    with pytest.raises(ValueError):
        shift_stokes(a, (0, 0))
