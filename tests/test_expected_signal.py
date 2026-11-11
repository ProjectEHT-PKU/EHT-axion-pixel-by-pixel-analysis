"""analysis.templates: preprocessing, response fit and template summary (synthetic cubes only)."""
import hashlib

import numpy as np
import pytest
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter
from scipy.sparse import csr_matrix

from eht_axion.analysis import expected_signal as stages


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, array):
    np.save(path, array)
    return path


def operator(n, bins=360):
    return csr_matrix(np.ones((bins, n * n)))


def use_operator(monkeypatch, op):
    """Replace the aperture quadrature by a synthetic operator of the right shape."""
    monkeypatch.setattr(stages, "_operator", lambda *args: op)


def phase_cube(times=3, n=3, *, states=17, amplitude=0.12, phase=0.7, mean=1.51):
    """Constant spatial images with an exact wrapped EVPA cosine oracle."""
    angles = mean + amplitude * np.cos(2 * np.pi * np.arange(16) / 16 + phase)
    if states == 17:
        angles = np.r_[-0.8, angles]  # Unrelated noaxion must not affect the fit.
    cube = np.ones((times, states, n, n, 4), dtype=np.float64)
    cube[..., 1] = np.cos(2 * angles)[None, :, None, None]
    cube[..., 2] = np.sin(2 * angles)[None, :, None, None]
    cube[..., 3] = -0.2
    return cube


def identity(n):
    """preprocess_frame settings that leave an n x n frame unchanged."""
    return dict(sigma=0, init_fov=3, final_fov=3, rotation_angle=0, final_shape=n)


FIT = dict(fov=3, radial_max=1.5, window_degrees=10, angular_bins=360)


@pytest.mark.parametrize("theta", [0, np.pi / 4, np.pi / 2, 0.6 * np.pi])
def test_preprocess_constant_spin2(theta):
    frame = np.broadcast_to([5.0, 2.0, -1.0, 0.3], (7, 7, 4)).copy()
    original = frame.copy()
    output = stages.preprocess_frame(frame, sigma=0.7, init_fov=7, final_fov=4,
                                     rotation_angle=theta, final_shape=4)
    assert output.shape == (4, 4, 4)
    assert output.dtype == np.float64
    np.testing.assert_allclose(output[..., 0], 5)
    np.testing.assert_allclose(output[..., 1], 2 * np.cos(2 * theta) + np.sin(2 * theta))
    np.testing.assert_allclose(output[..., 2], 2 * np.sin(2 * theta) - np.cos(2 * theta))
    np.testing.assert_allclose(output[..., 3], 0.3)
    np.testing.assert_array_equal(frame, original)


def test_preprocess_matches_independent_spline_blur_path():
    # Independent bilinear spline implementation of ImageS.Filter(blur=True).
    rng = np.random.default_rng(17)
    frame = rng.normal(size=(9, 9, 4))
    theta, sigma, initial, final, n = 0.43, 0.2, 9.0, 5.0, 6
    grid = (np.arange(9) - 4) * initial / 9
    out_grid = (np.arange(n) - (n - 1) / 2) * final / n
    x, y = np.meshgrid(out_grid, out_grid, indexing="ij")
    xv, yv = x * np.cos(theta) + y * np.sin(theta), -x * np.sin(theta) + y * np.cos(theta)
    sampled = np.empty((n, n, 4))
    for c in range(4):
        blurred = gaussian_filter(frame[..., c], sigma * 9 / initial, mode="nearest")
        sampled[..., c] = RectBivariateSpline(grid, grid, blurred, kx=1, ky=1).ev(xv, yv)
    q, u = sampled[..., 1].copy(), sampled[..., 2].copy()
    sampled[..., 1] = q * np.cos(2 * theta) - u * np.sin(2 * theta)
    sampled[..., 2] = q * np.sin(2 * theta) + u * np.cos(2 * theta)
    output = stages.preprocess_frame(frame, sigma=sigma, init_fov=initial, final_fov=final,
                                     rotation_angle=theta, final_shape=n)
    np.testing.assert_allclose(output, sampled, atol=2e-15)


def test_preprocess_rotation_crop_coordinates():
    x, y = np.meshgrid(np.arange(7) - 3, np.arange(7) - 3, indexing="ij")
    frame = np.zeros((7, 7, 4))
    frame[..., 0] = x + 10 * y
    output = stages.preprocess_frame(frame, sigma=0, init_fov=7, final_fov=3,
                                     rotation_angle=np.pi / 2, final_shape=3)
    xx, yy = np.meshgrid(np.arange(3) - 1, np.arange(3) - 1, indexing="ij")
    np.testing.assert_allclose(output[..., 0], yy - 10 * xx, atol=1e-14)


@pytest.mark.parametrize("phase", [-3.0, -0.4, 0.7, 3.0])
def test_dft_phase_oracle(phase):
    cube = phase_cube(times=1, states=16, phase=phase)
    angles = 0.5 * np.arctan2(cube[0, ..., 2], cube[0, ..., 1])
    fitted = stages.fit_phase_dft(angles)
    np.testing.assert_allclose(fitted, np.broadcast_to([0.12, phase, 1.51], (3, 3, 3)),
                               atol=2e-14)
    np.testing.assert_allclose(stages.fit_phase_dft(np.moveaxis(angles, 0, -1), phase_axis=-1),
                               fitted)


def test_dft_zero_amplitude_and_missing_phase():
    angles = np.full((16, 4), 0.2)
    angles[2, 1] = np.nan
    result = stages.fit_phase_dft(angles)
    np.testing.assert_array_equal(result[[0, 2, 3], 0], 0)
    assert np.isnan(result[:, 1]).all()  # No arbitrary phase at zero amplitude.
    assert np.isnan(result[1]).all()
    np.testing.assert_allclose(result[[0, 2, 3], 2], 0.2)
    for count in (15, 17):
        with pytest.raises(ValueError, match="16 phases"):
            stages.fit_phase_dft(np.ones((count, 2)))


def test_alpha_responses_exclude_state0():
    cube = phase_cube(times=4)
    cube[:, 0] = np.nan  # State0 even invalid must be entirely irrelevant.
    pixels, sectors = stages.alpha_responses(cube, identity(3), operator(3))
    assert pixels.shape == (4, 3, 3, 3) and sectors.shape == (4, 360, 3)
    for data in (pixels, sectors):
        assert data.dtype == np.float64
        np.testing.assert_allclose(data[..., 0], 0.12, atol=2e-14)
        np.testing.assert_allclose(data[..., 1], 0.7, atol=2e-14)
        np.testing.assert_allclose(data[..., 2], 1.51, atol=2e-14)
    with pytest.raises(ValueError, match="16 phase states"):
        stages.alpha_responses(cube[:, 1:], identity(3), operator(3))


def test_alpha_responses_are_the_snapshot_fits_of_the_sky_grid_images():
    op = operator(2, bins=4)
    cube = phase_cube(times=141, n=2)
    cube[0, 3, 0, 0, 1:3] = 0  # An invalid phase gives NaN, not a dropped snapshot.
    params = dict(sigma=0.3, init_fov=3, final_fov=2, rotation_angle=0.4, final_shape=2)
    pixels, sectors = stages.alpha_responses(cube, params, op)
    for t in (0, 70, 140):
        states = np.stack([stages.preprocess_frame(cube[t, s], **params) for s in range(1, 17)])
        pixel, sector = stages.snapshot_response(states, op)
        np.testing.assert_array_equal(pixels[t], pixel)
        np.testing.assert_array_equal(sectors[t], sector)
        np.testing.assert_array_equal(pixel, stages.fit_phase_dft(stages._evpa(states)))
        np.testing.assert_array_equal(sector, stages.fit_phase_dft(stages._integrate(states, op)))


def test_processed_background_keeps_the_state_axis():
    cube = phase_cube(times=3)[:, :1]
    background = stages.processed_background(cube, identity(3))
    assert background.shape == (3, 1, 3, 3, 4)
    np.testing.assert_allclose(background, cube, atol=1e-14)
    with pytest.raises(ValueError, match="axion-off"):
        stages.processed_background(phase_cube(times=3), identity(3))


def test_sector_response_roll_minus90_matches_eht():
    n = 20
    amplitude = np.linspace(0.01, 0.15, n * n).reshape(n, n)
    phases = np.linspace(-2, 2, n * n).reshape(n, n)
    psi = (2 * np.pi * np.arange(16) / 16)[:, None, None]
    angles = 0.2 + amplitude * np.cos(psi + phases)
    states = np.ones((16, n, n, 4))
    states[..., 1] = np.cos(2 * angles)
    states[..., 2] = np.sin(2 * angles)
    op = csr_matrix((np.ones(360), (np.arange(360), np.arange(360))), shape=(360, n * n))
    _, sector = stages.snapshot_response(states, op)
    np.testing.assert_allclose(sector[:, 0], np.roll(amplitude.ravel()[:360], -90), atol=1e-14)
    np.testing.assert_allclose(sector[:, 1], np.roll(phases.ravel()[:360], -90), atol=1e-14)


@pytest.mark.parametrize("spatial", [(2, 2), (360,)])
def test_expected_signal_multi_alpha_and_circular_phase(spatial):
    phases = np.array([np.pi - 0.1, -np.pi + 0.1, np.pi - 0.2, -np.pi + 0.2])
    amplitudes = np.array([1, 2, 3, 4])
    fits = np.zeros((4,) + spatial + (3,))
    fits[..., 0] = amplitudes.reshape((4,) + (1,) * len(spatial))
    fits[..., 1] = phases.reshape((4,) + (1,) * len(spatial))
    alphas, arrays = stages.expected_signal({0.3: fits, 0.1: fits.copy()})
    assert alphas == [0.1, 0.3]
    for name in ("alphas", "amplitude_mean", "amplitude_std", "phase_mean", "phase_std"):
        assert arrays[name].dtype == np.float64
    np.testing.assert_array_equal(arrays["alphas"], [0.1, 0.3])
    assert arrays["amplitude_mean"].shape == (2,) + spatial
    np.testing.assert_allclose(arrays["amplitude_mean"], 2.5)
    np.testing.assert_allclose(arrays["amplitude_std"], np.std(amplitudes))
    np.testing.assert_allclose(np.abs(arrays["phase_mean"]), np.pi)
    np.testing.assert_allclose(arrays["phase_std"], np.sqrt(0.025))
    np.testing.assert_array_equal(arrays["phase_count"], 4)


def test_expected_signal_missing_samples_and_shapes():
    fits = np.ones((3, 2, 3))
    fits[:2, 0, :2] = np.nan
    _, arrays = stages.expected_signal({0.2: fits})
    assert np.isnan(arrays["amplitude_mean"][0, 0])
    assert np.isnan(arrays["phase_std"][0, 0])
    assert arrays["amplitude_count"][0, 0] == 1
    assert arrays["amplitude_mean"][0, 1] == 1
    for bad in ({0.2: fits, 0.3: np.ones((3, 3, 3))}, {0.2: np.ones((3, 2, 2))}, {}):
        with pytest.raises(ValueError):
            stages.expected_signal(bad)
