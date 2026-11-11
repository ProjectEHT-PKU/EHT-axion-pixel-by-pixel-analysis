import numpy as np
import pytest

from eht_axion.raptor.layout import alpha_label, image_cube, product_directory


def test_alpha_identity_format():
    assert alpha_label(0.3) == "alpha-0.30"
    for alpha in [float("nan"), float("inf"), 0.305, 0, -1, True, 1e100]:
        with pytest.raises(ValueError):
            alpha_label(alpha)


def cube(path, times=1):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.zeros((times, 1, 2, 2, 4)))
    return path


def test_image_cube_paths(tmp_path):
    whole = cube(tmp_path / "res-100/alpha-0.20/stokes.npy")
    assert image_cube(tmp_path, 0.2, 100) == whole
    assert product_directory(tmp_path, None, 400) == tmp_path / "res-400/background"
    background = cube(tmp_path / "res-400/background/stokes.npy")
    assert image_cube(tmp_path, None, 400) == background
    with pytest.raises(FileNotFoundError):
        image_cube(tmp_path, 0.3, 400)
