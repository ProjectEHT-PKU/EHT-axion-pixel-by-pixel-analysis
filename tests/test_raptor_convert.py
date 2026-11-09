"""raptor.convert: core HDF5 blocks -> Stokes cube (synthetic files only)."""
import hashlib

import numpy as np
import pytest

from eht_axion.raptor import convert


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_hdf5_explicit_names_and_block_mapping(tmp_path):
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "data.h5"
    names = {"I": "group/intensity", "Q": "group/q", "U": "group/u", "V": "group/v"}
    origins = [[2, 2], [0, 0], [2, 0], [0, 2]]
    expected = np.zeros((4, 4, 4))
    with h5py.File(p, "w") as file:
        file.create_dataset("aaa_unrelated", data=np.zeros((40, 2)))
        for c, label in enumerate("IQUV"):
            blocks = np.arange(16, dtype="float32").reshape(4, 4) + c * 100
            file.create_dataset(names[label], data=blocks)
            for b, (x, y) in enumerate(origins):
                expected[x:x + 2, y:y + 2, c] = blocks[b].reshape(2, 2)
    before = digest(p)
    kwargs = dict(datasets=names, image_shape=(4, 4), block_origins=origins, pixel_shape=(2, 2))
    grid = convert.convert_hdf5([[p, p], [p, p]], tmp_path / "grid", **kwargs)
    array = np.load(grid["stokes"], allow_pickle=False)
    assert array.dtype == np.float64
    np.testing.assert_array_equal(array, np.broadcast_to(expected, (2, 2, 4, 4, 4)))
    assert digest(p) == before


def test_hdf5_bad_mappings(tmp_path):
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "blocks.h5"
    with h5py.File(p, "w") as file:
        for c, label in enumerate("IQUV"):
            file.create_dataset(label, data=np.full((1, 4), c))
    names = {c: c for c in "IQUV"}
    one = dict(datasets=names, image_shape=(2, 2), block_origins=[[0, 0]], pixel_shape=(2, 2))
    assert np.load(convert.convert_hdf5([[p]], tmp_path / "ok", **one)["stokes"]).shape == (1, 1, 2, 2, 4)
    with pytest.raises(ValueError, match="missing HDF5"):
        convert.convert_hdf5([[p]], tmp_path / "bad", **dict(one, datasets={**names, "Q": "missing"}))
    with pytest.raises(ValueError, match="overlapping"):
        convert.convert_hdf5([[p]], tmp_path / "bad", **dict(one, block_origins=[[0, 0], [0, 0]]))
    with pytest.raises(ValueError, match="cover"):
        convert.convert_hdf5([[p]], tmp_path / "bad", **dict(one, image_shape=(4, 4)))
    with pytest.raises(ValueError, match="rectangular"):
        convert.convert_hdf5([[p, p], [p]], tmp_path / "bad", **one)
    assert not (tmp_path / "bad").exists()
