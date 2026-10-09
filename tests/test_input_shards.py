"""Automatic res400 selection: tiny real cubes, metadata checks only."""

import numpy as np
import pytest

from eht_axion.pipeline.layout import input_raptor_files


DIRECTORIES = ("raptor/res-400/alpha-0.30", "raptor-stokes/res-400/alpha-0.3")


def cube(path, times=1, *, shape=None, dtype="float32"):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.zeros(shape or (times, 1, 2, 2, 4), dtype=dtype))
    return path


@pytest.mark.parametrize("directory", DIRECTORIES)
def test_numeric_sort_inclusive_endpoints_and_tiny_shards(tmp_path, directory):
    root = tmp_path / directory
    late = cube(root / "stokes-t1000-t1010.npy", times=2)
    early = cube(root / "stokes-t980-t990.npy", times=2)
    single = cube(root / "stokes-t1020-t1020.npy")
    assert input_raptor_files(tmp_path, 0.3, 400) == [early, late, single]


@pytest.mark.parametrize("directory", DIRECTORIES)
@pytest.mark.parametrize("times", [1, 3, 9])
def test_whole_file_has_no_production_length_requirement(tmp_path, directory, times):
    path = cube(tmp_path / directory / "stokes.npy", times=times)
    assert input_raptor_files(tmp_path, 0.3, 400) == [path]


@pytest.mark.parametrize("first,second,match", [
    ((900, 910, 2), (930, 940, 2), "gap"),
    ((900, 920, 3), (920, 930, 2), "overlap"),
    ((900, 930, 4), (910, 920, 2), "overlap"),
    ((900, 910, 2), (900, 910, 2), "overlap"),
])
def test_gaps_overlaps_and_duplicate_numeric_labels(tmp_path, first, second, match):
    root = tmp_path / DIRECTORIES[0]
    a, b, times = first
    cube(root / f"stokes-t{a}-t{b}.npy", times)
    a, b, times = second
    # Zero padding is legal spelling, but cannot create a second time identity.
    cube(root / f"stokes-t{a:05d}-t{b:05d}.npy", times)
    with pytest.raises(ValueError, match=match):
        input_raptor_files(tmp_path, 0.3, 400)


@pytest.mark.parametrize("name,times,match", [
    ("stokes-t900-t910.npy", 1, "time length.*frame count"),
    ("stokes-t900-t910.npy", 3, "time length.*frame count"),
    ("stokes-t900-t905.npy", 2, "step 10"),
    ("stokes-t905-t915.npy", 2, "step 10"),
    ("stokes-t920-t900.npy", 3, "ordered endpoints"),
    ("stokes-t900.0-t910.npy", 2, "shard name"),
    ("stokes-t900-t910-extra.npy", 2, "shard name"),
    ("stokes-tbad.npy", 1, "shard name"),
    ("stokes-t-10-t0.npy", 2, "shard name"),
])
def test_bad_frame_counts_steps_and_names(tmp_path, name, times, match):
    cube(tmp_path / DIRECTORIES[0] / name, times)
    with pytest.raises(ValueError, match=match):
        input_raptor_files(tmp_path, 0.3, 400)


@pytest.mark.parametrize("shape,dtype", [
    ((2,), "float32"), ((1, 2, 2, 4), "float32"),
    ((0, 1, 2, 2, 4), "float32"), ((1, 0, 2, 2, 4), "float32"),
    ((1, 1, 1, 1, 4), "float32"), ((1, 1, 2, 3, 4), "float32"),
    ((1, 1, 2, 2, 3), "float32"), ((1, 1, 2, 2, 4), "bool"),
    ((1, 1, 2, 2, 4), "complex64"), ((1, 1, 2, 2, 4), "object"),
])
@pytest.mark.parametrize("name", ["stokes.npy", "stokes-t900-t900.npy"])
def test_invalid_cube_headers_are_not_path_only_fixtures(tmp_path, shape, dtype, name):
    cube(tmp_path / DIRECTORIES[0] / name, shape=shape, dtype=dtype)
    with pytest.raises(ValueError, match="RAPTOR"):
        input_raptor_files(tmp_path, 0.3, 400)


@pytest.mark.parametrize("shape,dtype", [
    ((1, 2, 2, 2, 4), "float32"), ((1, 1, 3, 3, 4), "float32"),
    ((1, 1, 2, 2, 4), "float64"),
])
def test_shard_non_time_shape_and_dtype_must_match(tmp_path, shape, dtype):
    root = tmp_path / DIRECTORIES[0]
    cube(root / "stokes-t900-t900.npy")
    cube(root / "stokes-t910-t910.npy", shape=shape, dtype=dtype)
    with pytest.raises(ValueError, match="identical non-time shape and dtype"):
        input_raptor_files(tmp_path, 0.3, 400)


@pytest.mark.parametrize("left,right", [
    ((0, "stokes.npy"), (0, "stokes-t900-t900.npy")),
    ((1, "stokes.npy"), (1, "stokes-t900-t900.npy")),
    ((0, "stokes.npy"), (1, "stokes-t900-t900.npy")),
    ((0, "stokes-t900-t900.npy"), (1, "stokes-t900-t900.npy")),
    ((0, "stokes-t900-t900.npy"), (1, "stokes.npy")),
    ((0, "stokes.npy"), (1, "stokes.npy")),
])
def test_multiple_source_groups_fail_before_header_loading(tmp_path, monkeypatch, left, right):
    for index, name in (left, right):
        cube(tmp_path / DIRECTORIES[index] / name)

    def no_loading(*args, **kwargs):
        pytest.fail("ambiguity must be rejected before loading either source")

    monkeypatch.setattr(np, "load", no_loading)
    with pytest.raises(ValueError, match="Ambiguous res400"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_aliases_do_not_override_ambiguity_policy(tmp_path):
    path = cube(tmp_path / DIRECTORIES[0] / "stokes.npy")
    alias = tmp_path / DIRECTORIES[1] / "stokes.npy"
    alias.parent.mkdir(parents=True)
    alias.symlink_to(path)
    with pytest.raises(ValueError, match="Ambiguous res400"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_repeated_shard_inode_is_rejected(tmp_path):
    first = cube(tmp_path / DIRECTORIES[0] / "stokes-t900-t900.npy")
    alias = first.with_name("stokes-t910-t910.npy")
    alias.hardlink_to(first)
    with pytest.raises(ValueError, match="Duplicate RAPTOR"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_validation_uses_only_readonly_mmap_metadata(tmp_path, monkeypatch):
    root = tmp_path / DIRECTORIES[0]
    paths = [cube(root / f"stokes-t{t}-t{t}.npy") for t in (900, 910)]
    before = {p: (p.stat().st_size, p.stat().st_mtime_ns, p.read_bytes()) for p in paths}
    original_load = np.load
    loaded = []

    class HeaderOnly:
        def __init__(self, array):
            self.shape, self.dtype = array.shape, array.dtype

        def __array__(self, *args, **kwargs):
            pytest.fail("must not materialize pixels")

        def __getitem__(self, key):
            pytest.fail("must not read pixels")

    def load(path, **kwargs):
        assert kwargs == {"mmap_mode": "r", "allow_pickle": False}
        array = original_load(path, **kwargs)
        assert isinstance(array, np.memmap) and not array.flags.writeable
        loaded.append(path)
        return HeaderOnly(array)

    monkeypatch.setattr(np, "load", load)
    assert input_raptor_files(tmp_path, 0.3, 400) == paths
    assert loaded == paths
    assert {p: (p.stat().st_size, p.stat().st_mtime_ns, p.read_bytes()) for p in paths} == before


def test_no_source_and_corrupt_source(tmp_path):
    with pytest.raises(FileNotFoundError, match="No RAPTOR input"):
        input_raptor_files(tmp_path, 0.3, 400)
    path = tmp_path / DIRECTORIES[0] / "stokes.npy"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"not a numpy header")
    with pytest.raises(ValueError, match="Invalid RAPTOR input header"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_disguised_npz_is_not_a_single_cube(tmp_path):
    path = tmp_path / DIRECTORIES[0] / "stokes.npy"
    path.parent.mkdir(parents=True)
    with path.open("wb") as stream:
        np.savez(stream, stokes=np.zeros((1, 1, 2, 2, 4)))
    with pytest.raises(ValueError, match="single NPY cube"):
        input_raptor_files(tmp_path, 0.3, 400)


def test_truncated_npy_is_rejected_without_loading_pixels(tmp_path):
    path = cube(tmp_path / DIRECTORIES[0] / "stokes-t900-t900.npy")
    with path.open("r+b") as stream:
        stream.truncate(path.stat().st_size - 1)
    with pytest.raises(ValueError, match="Invalid RAPTOR input header"):
        input_raptor_files(tmp_path, 0.3, 400)
