"""eht.ensemble: native FITS reader, spin-2 EVPA, Paper VII image averages, per-day summaries (synthetic)."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from eht_axion.eht import ensemble as en
from eht_axion.eht.ensemble import evpa, load_stokes_fits


def test_global_diagnostics_paper_vii_definitions():
    """Paper VII eqs. 12-14: pure m=2 pattern P = 0.2 I exp(2i phi), phi east of north."""
    east = -(np.arange(16) - 7.5) * 1.875
    north = (np.arange(16) - 7.5) * 1.875
    ee, nn = np.meshgrid(east, north)
    phi = np.arctan2(ee, nn)
    image = np.zeros((16, 16, 4))
    image[..., 0] = np.exp(-(ee**2 + nn**2) / 200)
    p = .2 * image[..., 0] * np.exp(2j * phi)
    image[..., 1], image[..., 2] = p.real, p.imag
    total, m_net, m_resolved, chi, beta_re, beta_im = en.global_diagnostics(image, east, north)
    assert total == pytest.approx(image[..., 0].sum())
    assert m_net == pytest.approx(0, abs=1e-14)       # symmetric pattern: no net polarization
    assert m_resolved == pytest.approx(.2)            # <|m|> = sum |P| / sum I
    assert (beta_re, beta_im) == (pytest.approx(.2), pytest.approx(0, abs=1e-14))
    # A position-angle rotation by delta rotates beta_2 by exp(2i delta) only through P.
    rotated = image.copy()
    q = p * np.exp(2j * .3)
    rotated[..., 1], rotated[..., 2] = q.real, q.imag
    beta = complex(*en.global_diagnostics(rotated, east, north)[4:])
    assert beta == pytest.approx(.2 * np.exp(2j * .3))


def test_written_summary_contract(tmp_path):
    east = -(np.arange(8) - 3.5) * 1.875
    north = (np.arange(8) - 3.5) * 1.875
    samples = np.zeros((3, 8, 8, 4))
    samples[..., 0] = 2.
    samples[..., 1] = np.cos(2 * np.arange(3)[:, None, None] * .02)
    samples[..., 2] = np.sin(2 * np.arange(3)[:, None, None] * .02)
    samples[1, 0, 0, 0] = 0.  # nonpositive I: EVPA undefined
    en.summarize(samples, east, north, tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {'evpa.npy', 'global.npz'}
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
