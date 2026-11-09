"""Native sky/spin-2 regression tests of the per-day EHT ensemble summaries."""
import numpy as np
import pytest

from eht_axion.research import ensemble as en
from eht_axion.research.polarimetry import aperture_curves, evpa


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
