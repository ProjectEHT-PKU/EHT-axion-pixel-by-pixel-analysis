"""Physical angle comparisons, not science acceptance thresholds."""
import numpy as np
import pytest

from check_eht_outer_stability import compare, summary


def image(angle, fraction=.1):
    result = np.zeros((2, 3, 4))
    result[..., 0] = 1
    result[..., 1] = fraction*np.cos(2*np.deg2rad(angle))
    result[..., 2] = fraction*np.sin(2*np.deg2rad(angle))
    return result


def test_identical_images():
    a = image(25)
    r = compare(a, a.copy(), np.ones((2, 3), bool))
    assert r['evpa_rms_deg'] == pytest.approx(0, abs=1e-12)
    assert r['qu_relative_l2'] == 0
    assert r['integrated_evpa_difference_deg'] == pytest.approx(0, abs=1e-12)


def test_modulo_pi_seam():
    r = compare(image(89), image(-89), np.ones((2, 3), bool))
    assert r['evpa_rms_deg'] == pytest.approx(2.)
    assert r['integrated_evpa_difference_deg'] == pytest.approx(2.)


def test_complex_cancellation_is_not_mean_pixel_angle():
    a = image(0)
    a[0, :, 1] *= -1
    r = summary(a, np.zeros((2, 3)))
    assert r['net_m'] == pytest.approx(0)
    assert r['resolved_m'] == pytest.approx(.1)
    assert r['beta2_abs'] == pytest.approx(0)


def test_mask_does_not_modify_whole_image_net_statistic():
    a, b = image(0), image(5)
    mask = np.zeros((2, 3), bool)
    mask[0, 0] = True
    r = compare(a, b, mask)
    assert r['pixels'] == 1 and r['evpa_rms_deg'] == pytest.approx(5)
    assert r['net_m_a'] == pytest.approx(.1)
