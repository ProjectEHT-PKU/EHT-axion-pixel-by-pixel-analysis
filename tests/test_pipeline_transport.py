"""Small synthetic oracles only; no notebook execution or shared-data computation."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

from eht_axion.pipeline import transport as t
from scripts import transport_diagnostics


def grid(n, fov):
    g = (np.arange(n) - (n - 1) / 2) * fov / n
    return g, np.meshgrid(g, g, indexing="ij")


def make_inputs(root, n=41):
    root.mkdir()
    fov = 30.0
    _, (x, y) = grid(n, fov)
    radius = np.hypot(x, y)
    psi = 2 * np.pi * np.arange(16) / 16
    field = np.broadcast_to(0.08 * np.cos(psi[:, None, None] + 0.7), (16, n, n)).copy()
    baseline = np.zeros((n, n, 4))
    baseline[..., 0] = np.exp(-(radius / 8)**2)
    baseline[..., 1] = 0.3 * baseline[..., 0] * np.cos(3.0)
    baseline[..., 2] = 0.3 * baseline[..., 0] * np.sin(3.0)
    baseline[..., 3] = 0.02 * baseline[..., 0]
    np.save(root / "r_coord_400.npy", radius)
    np.save(root / "axion_field_coord_400_16.npy", field)
    np.save(root / "common_noaxion.npy", baseline)
    for name, factor in zip(t.VARIANTS, (1.0, 1.1, 0.9, 0.7)):
        cube = t.analytical_stokes(baseline, factor * field)
        np.save(root / f"{name}.npy", cube)
    return fov, field, baseline


def test_crossing_bilinear_residual_orientation_and_bounds():
    n, fov, r_s = 81, 30, 9.0
    g, (x, y) = grid(n, fov)
    radius = np.hypot(x, y)
    before = radius.copy()
    coords, info = t.get_360_x_y(radius, r_s, fov=fov, angular_bins=36,
                                return_diagnostics=True)
    assert info["valid"].all()
    np.testing.assert_array_equal(info["crossing_count"], 1)
    np.testing.assert_allclose(info["residual"], 0, atol=1e-10)
    np.testing.assert_allclose(coords[[0, 9, 18, 27]], [[0, 9], [-9, 0], [0, -9], [9, 0]],
                               atol=1e-10)
    np.testing.assert_allclose(RegularGridInterpolator((g, g), radius)(coords), r_s, atol=1e-10)
    np.testing.assert_array_equal(radius, before)
    assert np.isnan(t.get_360_x_y(radius, 100, fov=fov, angular_bins=8)).all()
    old = t.get_360_x_y(radius, r_s, fov=fov, angular_bins=36, crossing="legacy-bracket")
    np.testing.assert_allclose(coords, old, atol=1e-10)


def test_multiple_crossings_scanned_not_single_endpoint_bracket():
    _, (x, y) = grid(101, 24)
    radius = (np.hypot(x, y) - 3)**2 + 1
    coords, info = t.get_360_x_y(radius, 5, fov=24, angular_bins=8,
                                return_diagnostics=True)
    np.testing.assert_array_equal(info["crossing_count"], 2)
    np.testing.assert_allclose(np.hypot(*coords.T), 1, atol=0.02)
    assert np.isnan(t.get_360_x_y(radius, 5, fov=24, angular_bins=8,
                                  crossing="legacy-bracket")).all()


def test_mask_identical_all_channels_and_processing_order():
    n = 9
    source = np.zeros((n, n, 4))
    source[4, 4] = [1, 2, 3, 4]
    mask = np.ones((n, n), dtype=bool)
    mask[4, 4] = False
    before = source.copy()
    a = t.mask_and_blur_frame(source, mask, sigma=1, fov=n, mask_order=t.ORDERS[0])
    b = t.mask_and_blur_frame(source, mask, sigma=1, fov=n, mask_order=t.ORDERS[1])
    np.testing.assert_array_equal(a, 0)  # All IQUV removed before convolution.
    np.testing.assert_array_equal(b[4, 4], 0)  # All IQUV masked after convolution.
    assert np.all(b[4, 5] > 0)  # Contaminating unselected source has leaked outward.
    for channel in range(4):
        oracle = gaussian_filter(source[..., channel], 1, mode="nearest") * mask
        np.testing.assert_allclose(b[..., channel], oracle)
    np.testing.assert_array_equal(source, before)
    with pytest.raises(ValueError, match="boolean"):
        t.mask_and_blur_frame(source, mask.astype(float), sigma=1, fov=n, mask_order=t.ORDERS[0])


def test_legacy_rms_independent_oracle_and_dft_harmonic_distinction():
    psi = 2 * np.pi * np.arange(16) / 16
    evpa = 1.5 + 0.12 * np.cos(psi + 0.7) + 0.03 * np.cos(2 * psi - 0.4)
    wrapped = (evpa + np.pi / 2) % np.pi - np.pi / 2
    seventeen = np.r_[1.5, wrapped]
    np.testing.assert_allclose(t.legacy_rms_amplitude(seventeen), np.hypot(0.12, 0.03), atol=1e-14)
    np.testing.assert_allclose(t.fit_phase_dft(seventeen), [0.12, 0.7, 1.5], atol=1e-14)
    # Verify the old 16-state off-by-one centering instead of silently "fixing" it.
    mean16 = evpa[1:].mean()
    oracle = np.sqrt(np.sum((evpa - mean16)**2) / 8)
    np.testing.assert_allclose(t.legacy_rms_amplitude(wrapped), oracle, atol=1e-14)


def test_analytical_synthesis_and_phase_convention(tmp_path):
    fov, field, baseline = make_inputs(tmp_path / "inputs")
    before = baseline.copy()
    cube = t.analytical_stokes(baseline, field)
    assert cube.shape == (17, 41, 41, 4)
    np.testing.assert_allclose(cube[..., 0], np.broadcast_to(baseline[..., 0], (17, 41, 41)))
    np.testing.assert_allclose(np.hypot(cube[..., 1], cube[..., 2]),
                               np.broadcast_to(np.hypot(baseline[..., 1], baseline[..., 2]),
                                               (17, 41, 41)))
    np.testing.assert_array_equal(cube[..., 3], 0)  # Matches notebook analytical model.
    np.testing.assert_array_equal(baseline, before)
    radius = np.load(tmp_path / "inputs" / "r_coord_400.npy")
    coords = t.get_360_x_y(radius, 9, fov=fov, angular_bins=36)
    evpa = t.get_analytical_evpa(baseline, field, oneD_coord=coords, init_fov=fov)
    np.testing.assert_allclose(t.fit_phase_dft(evpa), np.tile([0.08, 0.7, 1.5], (36, 1)), atol=2e-14)
    # Deliberately unrelated baseline state must not alter the DFT fit.
    evpa[0] = np.nan
    np.testing.assert_allclose(t.fit_phase_dft(evpa), np.tile([0.08, 0.7, 1.5], (36, 1)), atol=2e-14)
    a = t.get_analytical_A_1d(field, radius, 9, fov=fov, angular_bins=36,
                              rotation_angle=0, amplitude_method="dft")
    np.testing.assert_allclose(a, 0.08, atol=1e-14)


def test_native_sampler_bilinear_qu_and_active_rotation():
    n, fov = 21, 20
    g, (x, y) = grid(n, fov)
    cube = np.zeros((16, n, n, 4))
    cube[..., 0] = 1
    cube[..., 1] = 1 + 0.01 * x
    cube[..., 2] = 0.01 * y
    phi = 2 * np.pi * np.arange(36) / 36
    coords = np.array([-5 * np.sin(phi), 5 * np.cos(phi)]).T
    rotation = np.pi / 2
    evpa = t.get_raptor_evpa(cube, oneD_coord=coords, init_fov=fov, rotation_angle=rotation)
    qu = RegularGridInterpolator((g, g), cube[0, ..., 1:3])(coords)
    native = np.arctan2(qu[:, 1], qu[:, 0]) / 2
    expected = (np.roll(native, 9) + rotation + np.pi / 2) % np.pi - np.pi / 2
    np.testing.assert_allclose(evpa, np.broadcast_to(expected, evpa.shape), atol=1e-14)


def test_spin2_periodic_interpolation_avoids_branch_seam():
    values = np.array([np.pi / 2 - 0.01, -np.pi / 2 + 0.01, 0.2, 0.3])
    theta = -np.pi / 4
    corrected = t._angular_sample(values, theta, "spin2-periodic")
    old = t._angular_sample(values, theta, "legacy-angle")
    # First target is halfway between nearly identical period-pi EVPA values.
    delta_correct = (corrected[0] - theta + np.pi / 2) % np.pi - np.pi / 2
    delta_old = (old[0] - theta + np.pi / 2) % np.pi - np.pi / 2
    assert abs(abs(delta_correct) - np.pi / 2) < 1e-12
    assert abs(delta_old) < 1e-12


def test_stage_products_readonly_inputs_and_no_runtime_notebook(tmp_path, monkeypatch):
    root = tmp_path / "inputs"
    fov, _, _ = make_inputs(root)
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob("*.npy")}
    original_open = Path.open
    original_load = np.load

    def no_notebook(path, *args, **kwargs):
        assert path.suffix != ".ipynb", "production must not open any notebook"
        return original_open(path, *args, **kwargs)

    def readonly_load(path, **kwargs):
        if Path(path) in hashes:
            assert kwargs == {"mmap_mode": "r", "allow_pickle": False}
        return original_load(path, **kwargs)

    monkeypatch.setattr(Path, "open", no_notebook)
    monkeypatch.setattr(np, "load", readonly_load)
    result = t.run_transport_diagnostics(root, tmp_path / "out", fov=fov, sigma=0.5,
                                         radii=(4, 9), r_width=2, angular_bins=36,
                                         amplitude_method="dft", rotation_angle=0.37)
    assert set(result) == {"diagnostics", "metadata", "amplitudes", "analytic_ratios",
                           "mask_order_comparison", "source_rings"}
    for name, path in result.items():
        assert path.is_file()
        if name not in ("diagnostics", "metadata"):
            assert path.read_bytes().startswith(b"\x89PNG")
            assert path.parent.name == "figures"
    meta = json.loads(result["metadata"].read_text())
    assert "inferred" in meta["variant_semantics_status"]
    assert "not verified" in meta["validation_status"]
    assert meta["method_choice"]["mask_channels"] == list("IQUV")
    assert meta["analytic_and_transport_rotation_matched"] is True
    assert meta["missing_crossings"] == [0, 0]
    with original_load(result["diagnostics"], allow_pickle=False) as data:
        assert data["fits_dft"].shape == (5, 2, 2, 2, 36, 3)
        assert data["analytic_ratio"].shape == (5, 2, 2, 2, 36)
        assert data["variants"].tolist()[:4] == list(t.VARIANTS)
        np.testing.assert_allclose(data["analytic_ratio"][:, :, :, 0],
                                   np.broadcast_to(np.array([1, 1.1, 0.9, 0.7, 1])[:, None, None, None],
                                                   (5, 2, 2, 36)), atol=1e-13)
        np.testing.assert_allclose(data["source_color_ratio"], 1.1, atol=1e-13)
        np.testing.assert_allclose(data["fits_dft"][0, ..., 1], 0.7, atol=1e-13)
        np.testing.assert_allclose(data["common_state0_minus_noaxion"][..., :3], 0, atol=1e-15)
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in hashes} == hashes
    with pytest.raises(FileExistsError):
        t.run_transport_diagnostics(root, tmp_path / "out", fov=fov, radii=(4, 9), angular_bins=8)


def test_script_outside_checkout_and_pointer_error(tmp_path):
    root = Path(__file__).resolve().parents[1]
    source = tmp_path / "inputs"
    fov, _, _ = make_inputs(source, n=21)
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    command = [sys.executable, str(root / "scripts" / "transport_diagnostics.py"), str(source),
               "--output-dir", str(tmp_path / "out"), "--fov", str(fov), "--sigma", "0.5",
               "--radii", "4", "9", "--r-width", "3", "--angular-bins", "8"]
    run = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert Path(json.loads(run.stdout)["metadata"]).is_file()
    (source / "common.npy").write_text("version https://git-lfs.github.com/spec/v1\n")
    run = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True)
    assert run.returncode != 0
    assert "Git LFS pointer" in run.stderr


def test_script_main_defaults_and_error_translation(tmp_path, monkeypatch, capsys):
    seen = {}

    def stage(**kwargs):
        seen.update(kwargs)
        return {"metadata": tmp_path / "transport.json"}

    monkeypatch.setattr(transport_diagnostics, "run_transport_diagnostics", stage)
    args = ["input", "--output-dir", str(tmp_path)]
    assert transport_diagnostics.main(args) == 0
    assert json.loads(capsys.readouterr().out) == {"metadata": str(tmp_path / "transport.json")}
    assert seen == dict(
        input_dir=Path("input"), output_dir=tmp_path, radii=[4, 5, 6, 7, 8, 9, 10],
        r_width=0.5, fov=46.2262, sigma=2.16871, rotation_angle=0, angular_bins=360,
        mask_order="mask-before-blur", amplitude_method="legacy-rms",
        angular_interpolation="spin2-periodic", crossing="first-outward", spin=0.9375,
        denominator_floor=1e-12,
    )

    def fail(**kwargs):
        raise ValueError("invalid synthetic input")

    monkeypatch.setattr(transport_diagnostics, "run_transport_diagnostics", fail)
    with pytest.raises(SystemExit) as error:
        transport_diagnostics.main(args)
    assert error.value.code == 2
    assert "invalid synthetic input" in capsys.readouterr().err
    assert not hasattr(t, "main")
    assert not hasattr(t, "argparse")


def test_missing_polarization_crossings_and_validation():
    cube = np.zeros((16, 5, 5, 4))
    coords = np.zeros((8, 2))
    coords[0] = np.nan
    assert np.isnan(t.get_raptor_evpa(cube, oneD_coord=coords, init_fov=5)).all()
    with pytest.raises(ValueError, match="16\\|17"):
        t.get_raptor_evpa(cube[:15], oneD_coord=coords)
    with pytest.raises(ValueError, match="r_s"):
        t.get_360_x_y(np.ones((5, 5)), -1)
    with pytest.raises(ValueError, match="angular_bins"):
        t.get_360_x_y(np.ones((5, 5)), 9, angular_bins=3)
    with pytest.raises(ValueError, match="including 4 and 9"):
        t.run_transport_diagnostics("unused", "unused", radii=[9])
    np.testing.assert_array_equal(t._ratio(np.array([1., 2.]), np.array([0., 1.]), 1e-12),
                                   [np.nan, 2.])
