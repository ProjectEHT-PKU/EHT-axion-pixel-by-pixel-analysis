"""Synthetic-only research polarimetry contract tests (no observational files)."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from eht_axion.research.polarimetry import aperture_curves, evpa, load_stokes_fits


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


def test_evpa():
    assert_allclose(evpa([1, 0, -1, 0], [0, 1, 0, -1]), np.arange(4)*np.pi/4)
    assert np.isnan(evpa([0, np.nan, 1], [0, 1, np.inf])).all()


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
    with pytest.raises(ValueError):
        aperture_curves(a, [0], [0], [0])
