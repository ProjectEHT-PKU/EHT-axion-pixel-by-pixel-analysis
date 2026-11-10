"""provenance: atomic stage directories, reuse of equal stages, code identity, verification."""
import json
import subprocess

import pytest

from eht_axion import provenance
from eht_axion.provenance import stage_directory, verify_stage

REAL_CODE_IDENTITY = provenance.code_identity  # captured before the autouse fixture replaces it


def test_atomic_stage_cache_and_invalidations(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("data")
    out = tmp_path / "stage"
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        (temporary / "value.txt").write_text("product")
    verify_stage(out)
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        assert temporary is None
    with pytest.raises(ValueError, match="A different stage exists"):
        with stage_directory(out, "test", {"a": 2}, [source]):
            pass
    with pytest.raises(RuntimeError, match="injected"):
        with stage_directory(tmp_path / "failed", "test", {"a": 2}, [source]) as temporary:
            (temporary / "value.txt").write_text("wrong")
            raise RuntimeError("injected")
    assert not (tmp_path / "failed").exists()
    assert (out / "value.txt").read_text() == "product"
    assert not out.with_suffix(".lock").exists()
    (out / "value.txt").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        verify_stage(out)


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True, text=True).stdout


def test_code_identity_is_the_commit_of_a_clean_tracked_tree(tmp_path):
    repository = tmp_path / "repo"
    repository.mkdir()
    git(repository, "init", "-q")
    (repository / "code.py").write_text("x = 1\n")
    git(repository, "add", "code.py")
    git(repository, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "c")
    commit = git(repository, "rev-parse", "HEAD").strip()
    (repository / "untracked.txt").write_text("ignored by the check")
    assert REAL_CODE_IDENTITY(repository) == {"git_commit": commit}
    (repository / "code.py").write_text("x = 2\n")
    with pytest.raises(ValueError, match="Uncommitted changes to tracked files"):
        REAL_CODE_IDENTITY(repository)
    with pytest.raises(ValueError, match="git checkout"):
        REAL_CODE_IDENTITY(tmp_path / "missing")


def test_manifest_records_commit_and_versions_only(tmp_path):
    raw = tmp_path / "raw.txt"
    raw.write_text("data")
    with stage_directory(tmp_path / "stage", "test", {}, [raw]) as out:
        (out / "value.txt").write_text("v")
    record = json.loads((tmp_path / "stage" / "manifest.json").read_text())
    assert record["implementation"] == {"git_commit": "0" * 40}
    assert "signature" not in record and "pixi_lock_sha256" not in record["runtime"]
    assert record["inputs"][0]["sha256"] == provenance.sha256(raw)


def test_dirty_tree_refuses_a_new_stage_before_writing(tmp_path, monkeypatch):
    def dirty(root=None):
        raise ValueError("Uncommitted changes to tracked files")
    monkeypatch.setattr(provenance, "code_identity", dirty)
    with pytest.raises(ValueError, match="Uncommitted"):
        with stage_directory(tmp_path / "stage", "test", {}, []) as out:
            (out / "value.txt").write_text("v")
    assert not (tmp_path / "stage").exists()
