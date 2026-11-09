"""provenance: atomic stage directories, cache signature and output verification."""
import pytest

from eht_axion.provenance import stage_directory, verify_stage


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
        with stage_directory(tmp_path / "failed", "test", {"a": 2}, [source]) as temporary:
            (temporary / "value.txt").write_text("wrong")
            raise RuntimeError("injected")
    assert not (tmp_path / "failed").exists()
    assert (out / "value.txt").read_text() == "product"
    assert not out.with_suffix(".lock").exists()
    (out / "value.txt").write_text("tampered")
    with pytest.raises(ValueError, match="checksum"):
        verify_stage(out)
