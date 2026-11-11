"""provenance: atomic stage directories, reuse of equal stages, code identity, unsafe paths."""
import json
import subprocess

import numpy as np
import pytest

from eht_axion import provenance
from eht_axion.provenance import stage_directory, verify_stage

REAL_CODE_IDENTITY = provenance.code_identity  # captured before the autouse fixture replaces it


def test_stage_is_atomic_reused_only_when_equal_and_records_paths(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("data")
    out = tmp_path / "stage"
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        (temporary / "value.txt").write_text("product")
        np.save(temporary / "array.npy", np.zeros((2, 3)))
    record = verify_stage(out)
    assert list(record) == ["stage", "config", "inputs", "implementation", "runtime", "complete",
                            "outputs", "elapsed_seconds"]
    assert record["inputs"] == [str(source)]
    assert record["implementation"] == {"git_commit": "0" * 40}
    assert {e["path"] for e in record["outputs"]} == {"value.txt", "array.npy"}
    assert next(e for e in record["outputs"] if e["path"] == "array.npy")["shape"] == [2, 3]
    assert "sha256" not in json.dumps(record)
    with stage_directory(out, "test", {"a": 1}, [source]) as temporary:
        assert temporary is None
    with pytest.raises(ValueError, match="A different stage exists"):
        with stage_directory(out, "test", {"a": 2}, [source]):
            pass
    with pytest.raises(RuntimeError, match="injected"):
        with stage_directory(tmp_path / "failed", "test", {}, [source]) as temporary:
            (temporary / "value.txt").write_text("wrong")
            raise RuntimeError("injected")
    assert not (tmp_path / "failed").exists()
    assert not out.with_suffix(".lock").exists()
    (out / "value.txt").unlink()
    with pytest.raises(ValueError, match="Missing stage output"):
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


def test_dirty_tree_refuses_a_new_stage_before_writing(tmp_path, monkeypatch):
    def dirty(root=None):
        raise ValueError("Uncommitted changes to tracked files")
    monkeypatch.setattr(provenance, "code_identity", dirty)
    with pytest.raises(ValueError, match="Uncommitted"):
        with stage_directory(tmp_path / "stage", "test", {}, []) as out:
            (out / "value.txt").write_text("v")
    assert not (tmp_path / "stage").exists()


@pytest.mark.parametrize("kind", ["leaf", "ancestor"])
def test_never_writes_through_directory_symlinks(tmp_path, kind):
    external = tmp_path / "external"
    external.mkdir()
    link = tmp_path / "link"
    link.symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match="shared directory symlink"):
        with stage_directory(link / "stage" if kind == "ancestor" else link, "test", {}):
            pytest.fail("must reject before writing")
    assert list(external.iterdir()) == []


def test_products_must_be_regular_files_and_outputs_never_contain_inputs(tmp_path):
    source = tmp_path / "source.npy"
    np.save(source, np.arange(5))
    with pytest.raises(ValueError, match="regular files"):
        with stage_directory(tmp_path / "stage", "test", {}, [source]) as temp:
            (temp / "value.npy").symlink_to(source)
    assert not (tmp_path / "stage").exists()
    with pytest.raises(ValueError, match="contains an input"):
        with stage_directory(tmp_path, "test", {}, [source]):
            pass
    assert source.is_file()
