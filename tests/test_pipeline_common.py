import numpy as np
import pytest
from eht_axion.pipeline.common import stage_directory, verify_stage
from eht_axion.pipeline.constants import PhysicalScale
from eht_axion.polarization import aperture_operator, integrate_stokes


def test_atomic_stage_cache_and_invalidations(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("data")
    out = tmp_path / "stage"
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        (temporary / "value.txt").write_text("product")
    verify_stage(out)
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        assert temporary is None
    with pytest.raises(ValueError, match="Stale cache"):
        with stage_directory(out, "test", {"a": 2}, [source]):
            pass
    with pytest.raises(RuntimeError, match="injected"):
        with stage_directory(out, "test", {"a": 2}, [source], rebuild=True) as temporary:
            (temporary / "value.txt").write_text("wrong")
            raise RuntimeError("injected")
    assert (out / "value.txt").read_text() == "product"
    assert not out.with_suffix(".lock").exists()
    (out / "value.txt").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        verify_stage(out)


def test_physical_scale_and_validation():
    scale = PhysicalScale(mass_solar=6.5e9, distance_mpc=16.8, label_step_tg=0.1)
    assert scale.tg_hours == pytest.approx(8.89324754)
    assert scale.angular_rg_uas == pytest.approx(3.81899325)
    with pytest.raises(ValueError):
        PhysicalScale(mass_solar=float("nan"), distance_mpc=16.8, label_step_tg=0.1)


def test_aperture_area_exact_and_weighting():
    operator = aperture_operator(
        shape=4, fov=4, radial_max=1.9, angular_bins=4, window_degrees=10, radial_resolution=0.2
    )
    area = 0.5 * 1.9**2 * np.deg2rad(10)
    np.testing.assert_allclose(np.asarray(operator.sum(axis=1)).ravel(), area, rtol=1e-12)
    image = np.zeros((4, 4, 4))
    image[..., 0] = 3
    image[..., 1] = 1
    image[..., 2] = 1
    np.testing.assert_allclose(integrate_stokes(image, operator), np.pi / 8)
    np.testing.assert_allclose(integrate_stokes(image, operator, "flux"), np.pi / 8)
    image[...] = 0
    assert np.isnan(integrate_stokes(image, operator)).all()


def test_window_center_symmetry():
    op = aperture_operator(
        shape=16, fov=4, radial_max=1.5, angular_bins=360, window_degrees=10, radial_resolution=0.1
    )
    # x-positive aperture is exactly symmetric under reflection y -> -y.
    row = op[0].toarray().reshape(16, 16)
    np.testing.assert_allclose(row, row[:, ::-1], atol=1e-14)
