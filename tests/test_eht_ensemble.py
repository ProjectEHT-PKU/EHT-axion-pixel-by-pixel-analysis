"""eht.ensemble: native FITS reader, sector apertures, spin-2 EVPA and per-day summaries (synthetic)."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from eht_axion.eht import ensemble as en
from eht_axion.eht.ensemble import aperture_curves, evpa, load_stokes_fits


def test_east_point_and_template_transpose():
    east = -(np.arange(64) - 31.5) * 1.875
    north = (np.arange(64) - 31.5) * 1.875
    image = np.zeros((64, 64, 4))
    # Slightly north of exact east on the even grid: well inside 10-degree sector.
    image[32, 16] = [2., 0., 1., 0.]
    curve = aperture_curves(image, east, north, np.arange(360.), window=10, rmax_uas=60.)
    assert curve['intensity'][90] == 2
    assert np.isnan(curve['evpa_rad'][0])
    assert curve['evpa_rad'][90] == pytest.approx(np.pi / 4)
    template = image.transpose(1, 0, 2)
    assert np.unravel_index(np.argmax(template[..., 0]), (64, 64)) == (16, 32)
    assert -east[16] < 0


def test_written_summary_contract(tmp_path):
    east = -(np.arange(8) - 3.5) * 1.875
    north = (np.arange(8) - 3.5) * 1.875
    samples = np.zeros((3, 8, 8, 4))
    samples[..., 0] = 2.
    samples[..., 1] = np.cos(2 * np.arange(3)[:, None, None] * .02)
    samples[..., 2] = np.sin(2 * np.arange(3)[:, None, None] * .02)
    samples[1, 0, 0, 0] = 0.  # nonpositive I: EVPA undefined
    en.summarize(samples, east, north, tmp_path, 60., 10)
    assert {p.name for p in tmp_path.iterdir()} == {'evpa.npy', 'aperture.npz', 'global.npz'}
    with np.load(tmp_path / 'aperture.npz', allow_pickle=False) as z:
        assert set(z.files) == {'angle_deg', 'intensity', 'polarized_brightness'}
        assert z['intensity'].shape == z['polarized_brightness'].shape == (3, 360)
        np.testing.assert_array_equal(z['angle_deg'], np.arange(360))
    chi = np.load(tmp_path / 'evpa.npy', mmap_mode='r')
    assert chi.shape == (3, 8, 8)
    assert np.isnan(chi[1, 0, 0])
    np.testing.assert_allclose(chi[2, 3, 3], .04, atol=1e-15)
    with np.load(tmp_path / 'global.npz', allow_pickle=False) as z:
        assert list(z['columns']) == list(en.GLOBAL_COLUMNS)
        assert z['values'].shape == (3, 6)
        np.testing.assert_allclose(z['values'][0, 1:3], .5)  # uniform field, |P|/I = 1/2


def test_evpa_spin2_convention():
    assert evpa(.4, .8) == pytest.approx(.5 * np.arctan2(.8, .4))
    assert np.isnan(evpa(0., 0.))
    # Q+iU spin-2 rotation changes chi by half its complex phase.
    p = (.4 + .8j) * np.exp(2j * .3)
    assert evpa(p.real, p.imag) - evpa(.4, .8) == pytest.approx(.3)


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
    result = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [0, 90, 360], window=10)
    assert_allclose(result["polarized_brightness"], [1, 3, 1])
    assert_allclose(result["evpa_rad"], [0, np.pi/4, 0])
    assert_allclose(result["intensity"], [1, 2, 1])
    wide = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [45], window=100)
    assert_allclose(wide["polarized_brightness"], [7/3])
    assert_allclose(wide["evpa_rad"], [.5*np.arctan2(6, 1)])
    cut = aperture_curves(a, [1, 0, -1], [-1, 0, 1], [0], window=10, rmax_uas=.5)
    assert np.isnan(cut["evpa_rad"]).all()


def test_invalid_inputs():
    a = np.zeros((3, 4, 4))
    with pytest.raises(ValueError):
        aperture_curves(a, [0], [0], [0], window=10)
