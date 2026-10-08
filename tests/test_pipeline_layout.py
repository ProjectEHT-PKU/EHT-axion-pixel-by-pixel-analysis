import numpy as np
import pytest
from eht_axion.pipeline.layout import alpha_label, date_label, input_raptor_files, input_eht_file
from eht_axion.pipeline.common import stage_directory, verify_stage, manifest_path


def test_dataset_identity_format():
    assert alpha_label(0.3) == "alpha-0.30"
    assert date_label(95) == "2017-04-05"
    assert date_label(101) == "2017-04-11"
    for alpha in [float("nan"), 0.305, 0, -1]:
        with pytest.raises(ValueError):
            alpha_label(alpha)


def test_canonical_and_legacy_input_resolution(tmp_path):
    p = tmp_path / "raptor/res-100/alpha-0.30/stokes.npy"
    p.parent.mkdir(parents=True)
    np.save(p, np.zeros((1, 17, 2, 2, 4)))
    assert input_raptor_files(tmp_path, 0.3, 100) == [p]
    q = tmp_path / "eht/2017-04-05/image-ensemble.pkl"
    q.parent.mkdir(parents=True)
    q.write_bytes(b"trusted-test")
    assert input_eht_file(tmp_path, 95) == q


def test_imported_manifest_preserves_link_target_and_requires_optin(tmp_path):
    source = tmp_path / "source.npy"
    np.save(source, np.arange(5))
    stage = tmp_path / "imported"
    config = {"preserved_setting": 3}
    with stage_directory(stage, "test", config, [source], imported_from={"original": "source"}) as temp:
        (temp / "value.npy").symlink_to(source)
    manifest = verify_stage(stage)
    assert manifest_path(stage).name == "manifest.json"
    assert manifest["scientific_recomputation"] is False
    assert manifest["outputs"][0]["shape"] == [5]
    with pytest.raises(ValueError, match="--reuse-imported"):
        with stage_directory(stage, "test", config, [source]):
            pass
    with stage_directory(stage, "test", config, [source], reuse_imported=True) as temp:
        assert temp is None
    with pytest.raises(ValueError, match="settings/lineage"):
        with stage_directory(stage, "test", {"preserved_setting": 4}, [source], reuse_imported=True):
            pass
    target = stage / "value.npy"
    target.unlink()
    target.symlink_to(tmp_path / "missing.npy")
    with pytest.raises(ValueError, match="Changed imported"):
        verify_stage(stage)
    assert source.is_file()


def test_recursive_lineage_refuses_changed_source(tmp_path):
    source = tmp_path / "raw.txt"
    source.write_text("original")
    stage = tmp_path / "processed"
    with stage_directory(stage, "test", {}, [source]) as temp:
        (temp / "product.txt").write_text("derived")
    source.write_text("changed")
    assert verify_stage(stage)["complete"]
    with pytest.raises(ValueError, match="Stale upstream lineage"):
        verify_stage(stage, recursive=True)


def test_rebuild_rejects_leaf_directory_link(tmp_path):
    target = tmp_path / "shared"
    target.mkdir()
    (target / "keep.txt").write_text("original")
    leaf = tmp_path / "run/alpha-0.30"
    leaf.parent.mkdir()
    leaf.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="shared directory symlink"):
        with stage_directory(leaf, "test", {}, [], rebuild=True):
            pass
    assert (target / "keep.txt").read_text() == "original"


def test_imported_upstream_requires_permission(tmp_path):
    source = tmp_path / "source.npy"
    np.save(source, np.arange(2))
    stage = tmp_path / "imported"
    with stage_directory(stage, "test", {}, [source], imported_from={"original": "source"}) as temp:
        (temp / "data.npy").symlink_to(source)
    with pytest.raises(ValueError, match="Consuming imported"):
        verify_stage(stage, recursive=True, allow_imported=False)
    assert verify_stage(stage, recursive=True, allow_imported=True)["complete"]


def test_legacy_manifest_is_readable_but_dual_is_ambiguous(tmp_path):
    stage = tmp_path / "stage"
    with stage_directory(stage, "test", {}, []) as temp:
        (temp / "data.txt").write_text("value")
    (stage / "manifest.json").rename(stage / "stage.json")
    assert verify_stage(stage)["complete"]
    (stage / "manifest.json").write_bytes((stage / "stage.json").read_bytes())
    with pytest.raises(ValueError, match="Ambiguous"):
        verify_stage(stage)
