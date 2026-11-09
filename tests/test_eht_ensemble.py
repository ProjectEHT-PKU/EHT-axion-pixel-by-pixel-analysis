"""Native sky/spin-2 and independent-day MC regression tests."""
import numpy as np
import pytest

from eht_axion.research import ensemble as en
from eht_axion.research.polarimetry import aperture_curves, evpa, shift_stokes


def test_seam_and_population_covariance():
    samples = np.deg2rad([[179., 10.], [1., 12.], [0., 11.]])
    summary = en.circular_summary(samples)
    assert abs(en.branch_difference(summary['mean_rad'][0], 0)) < 1e-12
    covariance = en.selected_covariance(samples, summary['mean_rad'])
    residual = en.branch_difference(samples, summary['mean_rad'])
    np.testing.assert_allclose(covariance, np.cov(residual.T, ddof=0), atol=1e-15)
    assert np.linalg.eigvalsh(covariance).min() >= -1e-15
    assert summary['valid'].all()


def test_independent_difference_is_not_paired_draw_covariance():
    a = np.array([[.1, .2], [.2, .3], [.3, .4]])
    b = a + .2
    sa, sb = en.circular_summary(a), en.circular_summary(b)
    ca = en.selected_covariance(a, sa['mean_rad'])
    cb = en.selected_covariance(b, sb['mean_rad'])
    result = en.independent_difference(sb, sa, cb, ca)
    np.testing.assert_allclose(result['difference_rad'], .2)
    np.testing.assert_allclose(result['covariance_rad2'], ca + cb)
    assert result['covariance_rad2'][0, 0] > 0
    assert np.var(b[:, 0] - a[:, 0]) < 1e-30


def test_masks_and_selection_shapes():
    a = np.zeros((4, 2, 2))
    a[:, 0, 0] = np.nan
    a[:, 1, 1] = np.arange(4) * np.pi / 4
    summary = en.circular_summary(a, .7)
    assert summary['valid'].shape == (2, 2)
    assert not summary['valid'][0, 0]
    assert not summary['valid'][1, 1]
    c = en.selected_covariance(a, summary['mean_rad'], [1, 2])
    assert c.shape == (2, 2)
    np.testing.assert_array_equal(c, 0)
    with pytest.raises(ValueError):
        en.circular_summary(a, 1.1)


def test_east_point_and_template_transpose():
    east = -(np.arange(64) - 31.5) * 1.875
    north = (np.arange(64) - 31.5) * 1.875
    image = np.zeros((64, 64, 4))
    # Slightly north of exact east on the even grid: well inside 10-degree sector.
    image[32, 16] = [2., 0., 1., 0.]
    curve = aperture_curves(image, east, north, np.arange(360.), rmax_uas=60.)
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
    summary, covariance = en.summarize(samples, east, north, tmp_path, .7, 60., 10)
    assert covariance.shape == (360, 360)
    assert summary['valid'].shape == (360,)
    with np.load(tmp_path / 'aperture.npz', allow_pickle=False) as z:
        assert z['samples_rad'].shape == (3, 360)
        np.testing.assert_array_equal(z['angle_deg'], np.arange(360))
    with np.load(tmp_path / 'map-summary.npz', allow_pickle=False) as z:
        assert z['valid'].shape == (8, 8)
        assert z['valid'].all()
        np.testing.assert_allclose(z['mean_rad'], .02, atol=1e-15)
    assert np.load(tmp_path / 'evpa.npy', mmap_mode='r').shape == (3, 8, 8)


def test_centroid_shift_no_spin_rotation():
    east = -(np.arange(65) - 32) * 1.875
    north = (np.arange(65) - 32) * 1.875
    image = np.zeros((65, 65, 4))
    image[35, 28] = [2., .4, .8, .1]
    offset = en.centroid_offset(image[..., 0], east, north)
    np.testing.assert_allclose(offset, [-3., 4.])
    shifted = shift_stokes(image, offset)
    np.testing.assert_allclose(shifted[32, 32], image[35, 28])
    assert evpa(shifted[32, 32, 1], shifted[32, 32, 2]) == evpa(.4, .8)
    assert np.isnan(evpa(0., 0.))
    # Q+iU spin-2 rotation changes chi by half its complex phase.
    p = (.4 + .8j) * np.exp(2j * .3)
    assert en.branch_difference(evpa(p.real, p.imag), evpa(.4, .8)) == pytest.approx(.3)
