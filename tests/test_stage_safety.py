"""Regression tests for output links, ancestry and environment-aware reuse."""

import json
import numpy as np
import pytest

from eht_axion.pipeline import common
from eht_axion.pipeline.raptor import preprocess_stokes


@pytest.mark.parametrize("rebuild", [False, True])
@pytest.mark.parametrize("kind", ["leaf", "ancestor", "existing"])
def test_stage_rejects_all_directory_links_before_write(tmp_path, rebuild, kind):
    external = tmp_path / "external"
    link = tmp_path / "link"
    if kind == "existing":
        external.mkdir()
    link.symlink_to(external, target_is_directory=True)
    destination = link / "stage" if kind == "ancestor" else link
    with pytest.raises(ValueError, match="shared directory symlink"):
        with common.stage_directory(destination, "test", {}, rebuild=rebuild):
            pytest.fail("must reject before writing")
    assert not (external / "stage").exists()
    if kind != "existing":
        assert not external.exists()


def test_lowlevel_rejects_dangling_file_link(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    external = tmp_path / "external.npy"
    (out / "stokes.npy").symlink_to(external)
    with pytest.raises(ValueError, match="output file symlink"):
        preprocess_stokes(np.ones((1, 2, 2, 4)), out, sigma=0, init_fov=2,
                          final_fov=2, rotation_angle=0, final_shape=2)
    assert not external.exists()


def test_origin_inputs_are_validated_recursively(tmp_path):
    raw = tmp_path / "raw.txt"
    raw.write_text("v1")
    origin = tmp_path / "origin"
    with common.stage_directory(origin, "test", {}, [raw]) as t:
        (t / "value.txt").write_text("original science")
    imported = tmp_path / "imported"
    with common.stage_directory(imported, "test", {}, [],
                                imported_from={"path": str(origin)},
                                origin_paths=[origin / "manifest.json", origin / "value.txt"]) as t:
        (t / "value.txt").write_text("converted")
    common.verify_stage(imported, recursive=True)
    raw.write_text("v2")
    with pytest.raises(ValueError, match="Stale upstream lineage"):
        common.verify_stage(imported, recursive=True)


def test_imported_ancestor_cannot_be_laundered(tmp_path):
    imported = tmp_path / "imported"
    with common.stage_directory(imported, "test", {}, [], imported_from={"path": "archive"}) as t:
        (t / "value.txt").write_text("old")
    child = tmp_path / "child"
    inputs = [imported / "manifest.json", imported / "value.txt"]
    with pytest.raises(ValueError, match="reuse-imported"):
        with common.stage_directory(child, "test", {}, inputs):
            pytest.fail("must reject imported ancestor")
    assert not child.exists()
    with common.stage_directory(child, "test", {}, inputs, reuse_imported=True) as t:
        (t / "value.txt").write_text("new calculation using old inputs")
    record = common.verify_stage(child, recursive=True)
    assert record["imported_ancestors"]
    with pytest.raises(ValueError, match="reuse-imported"):
        common.verify_stage(child, recursive=True, allow_imported=False)


def test_legacy_migration_receipt_is_not_verified_origin(tmp_path):
    origin = tmp_path / "origin"
    with common.stage_directory(origin, "test", {}, []) as t:
        (t / "x.txt").write_text("original")
    target = tmp_path / "migration"
    receipt = {"path": str(origin), "manifest_sha256": common.sha256(origin / "manifest.json")}
    with common.stage_directory(target, "test", {}, [], imported_from=receipt) as t:
        (t / "x.txt").write_text("converted")
    with pytest.raises(ValueError, match="Unverifiable legacy migration origin"):
        common.verify_stage(target, recursive=True)


def test_runtime_change_invalidates_cache(tmp_path, monkeypatch):
    target = tmp_path / "stage"
    with common.stage_directory(target, "test", {}, []) as t:
        (t / "x.txt").write_text("v1")
    original = common.runtime_identity()
    monkeypatch.setattr(common, "runtime_identity", lambda: dict(original, python="different"))
    with pytest.raises(ValueError, match="Stale cache"):
        with common.stage_directory(target, "test", {}, []):
            pass


def test_rebuild_cannot_delete_recursive_raw_source(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    raw = target / "raw.txt"
    raw.write_text("must survive")
    external = tmp_path / "external-stage"
    with common.stage_directory(external, "test", {}, [raw]) as t:
        (t / "value.txt").write_text("derived")
    with pytest.raises(ValueError, match="contains an input"):
        with common.stage_directory(target, "test", {}, [external / "manifest.json"], rebuild=True):
            pytest.fail("must reject indirect raw containment before replacement")
    assert raw.read_text() == "must survive"
    assert not (tmp_path / "target.previous").exists()


def test_dangling_backup_never_replaced(tmp_path):
    target = tmp_path / "stage"
    backup = tmp_path / "stage.previous"
    backup.symlink_to(tmp_path / "missing", target_is_directory=True)
    with pytest.raises(ValueError, match="Backup already exists"):
        with common.stage_directory(target, "test", {}, []) as t:
            (t / "value.txt").write_text("new")
    assert backup.is_symlink()
    assert not target.exists()


def test_stale_manifest_input_refused_before_child_write(tmp_path):
    raw = tmp_path / "raw.txt"
    raw.write_text("v1")
    parent = tmp_path / "parent"
    with common.stage_directory(parent, "test", {}, [raw]) as t:
        (t / "value.txt").write_text("p")
    raw.write_text("v2")
    child = tmp_path / "child"
    with pytest.raises(ValueError, match="Stale upstream lineage"):
        with common.stage_directory(child, "test", {}, [parent / "manifest.json"]):
            pass
    assert not child.exists()


def test_runtime_and_origins_are_recorded(tmp_path):
    raw = tmp_path / "raw.txt"
    raw.write_text("v1")
    out = tmp_path / "stage"
    with common.stage_directory(out, "test", {}, [], imported_from={"path": "archive"},
                                origin_paths=[raw]) as t:
        (t / "value.txt").write_text("converted")
    record = json.loads((out / "manifest.json").read_text())
    assert record["runtime"]["packages"]["numpy"]
    assert record["runtime"]["pixi_lock_sha256"]
    assert record["origin_inputs"][0]["path"] == str(raw)
    assert record["scientific_recomputation"] is False
