"""Regression tests for output links, ancestry and environment-aware reuse."""

import json
import numpy as np
import pytest

from eht_axion.pipeline import common
from eht_axion.pipeline.raptor import preprocess_cube


@pytest.mark.parametrize("kind", ["leaf", "ancestor", "existing"])
def test_stage_rejects_all_directory_links_before_write(tmp_path, kind):
    external = tmp_path / "external"
    link = tmp_path / "link"
    if kind == "existing":
        external.mkdir()
    link.symlink_to(external, target_is_directory=True)
    destination = link / "stage" if kind == "ancestor" else link
    with pytest.raises(ValueError, match="shared directory symlink"):
        with common.stage_directory(destination, "test", {}):
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
        np.save(tmp_path / "cube.npy", np.ones((1, 1, 2, 2, 4)))
        preprocess_cube(tmp_path / "cube.npy", out, sigma=0, init_fov=2,
                        final_fov=2, rotation_angle=0, final_shape=2)
    assert not external.exists()


def test_runtime_change_invalidates_cache(tmp_path, monkeypatch):
    target = tmp_path / "stage"
    with common.stage_directory(target, "test", {}, []) as t:
        (t / "x.txt").write_text("v1")
    original = common.runtime_identity()
    monkeypatch.setattr(common, "runtime_identity", lambda: dict(original, python="different"))
    with pytest.raises(ValueError, match="Stale cache"):
        with common.stage_directory(target, "test", {}, []):
            pass


def test_existing_target_cannot_contain_recursive_raw_source(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    raw = target / "raw.txt"
    raw.write_text("must survive")
    external = tmp_path / "external-stage"
    with common.stage_directory(external, "test", {}, [raw]) as t:
        (t / "value.txt").write_text("derived")
    with pytest.raises(ValueError, match="contains an input"):
        with common.stage_directory(target, "test", {}, [external / "manifest.json"]):
            pytest.fail("must reject indirect raw containment before any write")
    assert raw.read_text() == "must survive"


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


def test_runtime_is_recorded_and_outputs_verified(tmp_path):
    raw = tmp_path / "raw.txt"
    raw.write_text("v1")
    out = tmp_path / "stage"
    with common.stage_directory(out, "test", {}, [raw]) as t:
        (t / "value.txt").write_text("derived")
    record = json.loads((out / "manifest.json").read_text())
    assert record["runtime"]["packages"]["numpy"]
    assert record["runtime"]["pixi_lock_sha256"]
    assert record["inputs"][0]["path"] == str(raw)
    assert not out.with_suffix(".lock").exists()
    (out / "value.txt").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        common.verify_stage(out)
