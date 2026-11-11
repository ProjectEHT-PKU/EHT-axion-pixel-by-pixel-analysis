import numpy as np
import pytest

from eht_axion.analysis.templates import RunLayout
from eht_axion.raptor.layout import alpha_label, delivered_cube, product_directory


def test_alpha_identity_format():
    assert alpha_label(0.3) == "alpha-0.30"
    for alpha in [float("nan"), float("inf"), 0.305, 0, -1, True, 1e100]:
        with pytest.raises(ValueError):
            alpha_label(alpha)


def cube(path, times=1):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.zeros((times, 1, 2, 2, 4)))
    return path


def test_delivered_cube_paths(tmp_path):
    whole = cube(tmp_path / "res-100/alpha-0.20/stokes.npy")
    assert delivered_cube(tmp_path, 0.2, 100) == whole
    assert product_directory(tmp_path, None, 400) == tmp_path / "res-400/background"
    background = cube(tmp_path / "res-400/background/stokes.npy")
    assert delivered_cube(tmp_path, None, 400) == background
    with pytest.raises(FileNotFoundError):
        delivered_cube(tmp_path, 0.3, 400)


def test_run_layout_directories(tmp_path):
    layout = RunLayout(tmp_path)
    assert layout.processed(0.1) == tmp_path / "processed/alpha-0.10"
    assert layout.response(0.44) == tmp_path / "response/alpha-0.44"
    assert layout.templates("2d") == tmp_path / "templates/2d"
    assert layout.background_processed == tmp_path / "background"
