import numpy as np
import pytest

from eht_axion.provenance import manifest_path, stage_directory, verify_stage
from eht_axion.analysis.templates import RunLayout
from eht_axion.raptor.layout import BACKGROUND_DIRECTORY_ALPHA, alpha_label, input_raptor_files


def test_alpha_identity_format():
    assert alpha_label(0.3) == "alpha-0.30"
    for alpha in [float("nan"), float("inf"), 0.305, 0, -1, True, 1e100]:
        with pytest.raises(ValueError):
            alpha_label(alpha)


def cube(path, times=1):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.zeros((times, 1, 2, 2, 4)))
    return path


def test_response_whole_cube_and_background_shards(tmp_path):
    whole = cube(tmp_path / "raptor/res-100/alpha-0.20/stokes.npy")
    assert input_raptor_files(tmp_path, 0.2, 100) == [whole]
    base = tmp_path / "raptor/res-400" / alpha_label(BACKGROUND_DIRECTORY_ALPHA)
    late = cube(base / "t1700-t1790/stokes-t1700-t1790.npy")
    early = cube(base / "t1600-t1690/stokes-t1600-t1690.npy")
    assert input_raptor_files(tmp_path, BACKGROUND_DIRECTORY_ALPHA, 400) == [early, late]


def test_missing_ambiguous_and_misnamed_products(tmp_path):
    with pytest.raises(FileNotFoundError):
        input_raptor_files(tmp_path, 0.3, 400)
    base = tmp_path / "raptor/res-400/alpha-0.30"
    cube(base / "stokes.npy")
    cube(base / "t1600-t1690/stokes-t1600-t1690.npy")
    with pytest.raises(ValueError, match="Ambiguous"):
        input_raptor_files(tmp_path, 0.3, 400)
    (base / "stokes.npy").unlink()
    cube(base / "t1700-t1790/stokes-t1700-t1780.npy")
    with pytest.raises(ValueError, match="shard name"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_run_layout_directories(tmp_path):
    layout = RunLayout(tmp_path)
    assert layout.processed(0.1) == tmp_path / "raptor/processed/alpha-0.10"
    assert layout.response(0.44) == tmp_path / "raptor/response/alpha-0.44"
    assert layout.templates("2d") == tmp_path / "raptor/templates/2d"
    assert layout.background_processed == tmp_path / "raptor/background"


def test_recursive_lineage_refuses_changed_source(tmp_path):
    source = tmp_path / "raw.txt"
    source.write_text("original")
    stage = tmp_path / "processed"
    with stage_directory(stage, "test", {}, [source]) as temp:
        (temp / "product.txt").write_text("derived")
    assert manifest_path(stage).name == "manifest.json"
    source.write_text("changed")
    assert verify_stage(stage)["complete"]
    with pytest.raises(ValueError, match="Stale upstream lineage"):
        verify_stage(stage, recursive=True)


def test_stage_products_must_not_be_symlinks(tmp_path):
    source = tmp_path / "source.npy"
    np.save(source, np.arange(5))
    with pytest.raises(ValueError, match="regular files"):
        with stage_directory(tmp_path / "stage", "test", {}, [source]) as temp:
            (temp / "value.npy").symlink_to(source)
    assert not (tmp_path / "stage").exists()


def test_stage_rejects_leaf_directory_link(tmp_path):
    target = tmp_path / "shared"
    target.mkdir()
    (target / "keep.txt").write_text("original")
    leaf = tmp_path / "run/alpha-0.30"
    leaf.parent.mkdir()
    leaf.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="shared directory symlink"):
        with stage_directory(leaf, "test", {}, []):
            pass
    assert (target / "keep.txt").read_text() == "original"
