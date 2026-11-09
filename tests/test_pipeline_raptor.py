"""Synthetic-only tests: no notebook execution and no repository data cubes."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter
from scipy.sparse import csr_matrix

from eht_axion.pipeline import raptor


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, array):
    np.save(path, array)
    return path


def operator(n, bins=360):
    return csr_matrix(np.ones((bins, n * n)))


def use_operator(monkeypatch, op):
    """Replace the aperture quadrature by a synthetic operator of the right shape."""
    monkeypatch.setattr(raptor, "_operator", lambda *args: op)


def phase_cube(times=3, n=3, *, states=17, amplitude=0.12, phase=0.7, mean=1.51):
    """Constant spatial images with an exact wrapped EVPA cosine oracle."""
    angles = mean + amplitude * np.cos(2 * np.pi * np.arange(16) / 16 + phase)
    if states == 17:
        angles = np.r_[-0.8, angles]  # Unrelated noaxion must not affect the fit.
    cube = np.ones((times, states, n, n, 4), dtype=np.float64)
    cube[..., 1] = np.cos(2 * angles)[None, :, None, None]
    cube[..., 2] = np.sin(2 * angles)[None, :, None, None]
    cube[..., 3] = -0.2
    return cube


def processed(tmp_path, cube, name="processed"):
    directory = tmp_path / name
    directory.mkdir()
    save(directory / "stokes.npy", cube)
    return directory


def response(path, fits, labels):
    np.savez(path, amplitude=fits[..., 0], phase=fits[..., 1], baseline=fits[..., 2],
             snapshot_labels=np.asarray(labels))
    return path


FIT = dict(fov=3, radial_max=1.5, window_degrees=10, angular_bins=360)


@pytest.mark.parametrize("theta", [0, np.pi / 4, np.pi / 2, 0.6 * np.pi])
def test_preprocess_constant_spin2(theta):
    frame = np.broadcast_to([5.0, 2.0, -1.0, 0.3], (7, 7, 4)).copy()
    original = frame.copy()
    output = raptor.preprocess_frame(frame, sigma=0.7, init_fov=7, final_fov=4,
                                     rotation_angle=theta, final_shape=4)
    assert output.shape == (4, 4, 4)
    assert output.dtype == np.float64
    np.testing.assert_allclose(output[..., 0], 5)
    np.testing.assert_allclose(output[..., 1], 2 * np.cos(2 * theta) + np.sin(2 * theta))
    np.testing.assert_allclose(output[..., 2], 2 * np.sin(2 * theta) - np.cos(2 * theta))
    np.testing.assert_allclose(output[..., 3], 0.3)
    np.testing.assert_array_equal(frame, original)


def test_preprocess_matches_independent_spline_blur_path():
    # Independent bilinear spline implementation of ImageS.Filter(blur=True).
    rng = np.random.default_rng(17)
    frame = rng.normal(size=(9, 9, 4))
    theta, sigma, initial, final, n = 0.43, 0.2, 9.0, 5.0, 6
    grid = (np.arange(9) - 4) * initial / 9
    out_grid = (np.arange(n) - (n - 1) / 2) * final / n
    x, y = np.meshgrid(out_grid, out_grid, indexing="ij")
    xv, yv = x * np.cos(theta) + y * np.sin(theta), -x * np.sin(theta) + y * np.cos(theta)
    sampled = np.empty((n, n, 4))
    for c in range(4):
        blurred = gaussian_filter(frame[..., c], sigma * 9 / initial, mode="nearest")
        sampled[..., c] = RectBivariateSpline(grid, grid, blurred, kx=1, ky=1).ev(xv, yv)
    q, u = sampled[..., 1].copy(), sampled[..., 2].copy()
    sampled[..., 1] = q * np.cos(2 * theta) - u * np.sin(2 * theta)
    sampled[..., 2] = q * np.sin(2 * theta) + u * np.cos(2 * theta)
    output = raptor.preprocess_frame(frame, sigma=sigma, init_fov=initial, final_fov=final,
                                     rotation_angle=theta, final_shape=n)
    np.testing.assert_allclose(output, sampled, atol=2e-15)


def test_preprocess_rotation_crop_coordinates():
    x, y = np.meshgrid(np.arange(7) - 3, np.arange(7) - 3, indexing="ij")
    frame = np.zeros((7, 7, 4))
    frame[..., 0] = x + 10 * y
    output = raptor.preprocess_frame(frame, sigma=0, init_fov=7, final_fov=3,
                                     rotation_angle=np.pi / 2, final_shape=3)
    xx, yy = np.meshgrid(np.arange(3) - 1, np.arange(3) - 1, indexing="ij")
    np.testing.assert_allclose(output[..., 0], yy - 10 * xx, atol=1e-14)


def test_preprocess_shards_memmap_and_paths(tmp_path, monkeypatch):
    cube = phase_cube(times=5)
    whole = save(tmp_path / "whole.npy", cube)
    shards = [save(tmp_path / "late-name.npy", cube[:2]),
              save(tmp_path / "early-name.npy", cube[2:])]
    before = {p: digest(p) for p in [whole, *shards]}
    original_load = np.load
    loads = []

    def checked_load(path, **kwargs):
        array = original_load(path, **kwargs)
        if Path(path) in before:
            assert kwargs == {"mmap_mode": "r", "allow_pickle": False}
            assert isinstance(array, np.memmap)
            assert not array.flags.writeable
            loads.append(Path(path))
        return array

    monkeypatch.setattr(np, "load", checked_load)
    params = dict(sigma=0.1, init_fov=3, final_fov=2, rotation_angle=0.3, final_shape=2)
    one = raptor.preprocess_cube(whole, tmp_path / "a", **params)
    many = raptor.preprocess_cube(shards, tmp_path / "new" / "b", **params)
    first, second = original_load(one["stokes"]), original_load(many["stokes"])
    assert first.shape == (5, 17, 2, 2, 4)
    assert first.dtype == np.float64
    np.testing.assert_array_equal(first, second)
    assert set(loads) == set(before)
    assert {p: digest(p) for p in before} == before
    metadata = json.loads(many["metadata"].read_text())
    assert metadata["inputs"] == [str(p) for p in shards]
    assert many["stokes"] == tmp_path / "new" / "b" / "stokes.npy"
    assert many["metadata"].name == "processing.json"


@pytest.mark.parametrize("phase", [-3.0, -0.4, 0.7, 3.0])
def test_dft_phase_oracle(phase):
    cube = phase_cube(times=1, states=16, phase=phase)
    angles = 0.5 * np.arctan2(cube[0, ..., 2], cube[0, ..., 1])
    fitted = raptor.fit_phase_dft(angles)
    np.testing.assert_allclose(fitted, np.broadcast_to([0.12, phase, 1.51], (3, 3, 3)),
                               atol=2e-14)
    np.testing.assert_allclose(raptor.fit_phase_dft(np.moveaxis(angles, 0, -1), phase_axis=-1),
                               fitted)


def test_dft_zero_amplitude_and_missing_phase():
    angles = np.full((16, 4), 0.2)
    angles[2, 1] = np.nan
    result = raptor.fit_phase_dft(angles)
    np.testing.assert_array_equal(result[[0, 2, 3], 0], 0)
    assert np.isnan(result[:, 1]).all()  # No arbitrary phase at zero amplitude.
    assert np.isnan(result[1]).all()
    np.testing.assert_allclose(result[[0, 2, 3], 2], 0.2)
    for count in (15, 17):
        with pytest.raises(ValueError, match="16 phases"):
            raptor.fit_phase_dft(np.ones((count, 2)))


def test_fit_processed_excludes_state0(tmp_path, monkeypatch):
    use_operator(monkeypatch, operator(3))
    cube = phase_cube(times=4)
    cube[:, 0] = np.nan  # State0 even invalid must be entirely irrelevant.
    outputs = raptor.fit_processed(processed(tmp_path, cube), tmp_path / "fits", **FIT,
                                   snapshot_labels=[1600, 1610, 1620, 1630])
    for name, shape in [("fits_2d", (4, 3, 3)), ("fits_1d", (4, 360))]:
        assert outputs[name].name == name.replace("fits_", "response-") + ".npz"
        with np.load(outputs[name], allow_pickle=False) as data:
            assert set(data.files) == {"amplitude", "phase", "baseline", "snapshot_labels"}
            for key, value in [("amplitude", 0.12), ("phase", 0.7), ("baseline", 1.51)]:
                assert data[key].shape == shape
                assert data[key].dtype == np.float64
                np.testing.assert_allclose(data[key], value, atol=2e-14)
            np.testing.assert_array_equal(data["snapshot_labels"], [1600, 1610, 1620, 1630])
    assert outputs["metadata"].name == "response.json"
    meta = json.loads(outputs["metadata"].read_text())
    assert meta["weighting"] == "intensity"
    assert meta["state0_used"] is False


def test_response_exact_kernel_columns_and_snapshot_labels(tmp_path, monkeypatch):
    # 141 tiny frames exercise the production label convention, not production data.
    op = operator(2, bins=4)
    use_operator(monkeypatch, op)
    cube = phase_cube(times=141, n=2)
    cube[0, 3, 0, 0, 1:3] = 0  # Invalid phase: NaNs must also survive serialization.
    labels = np.arange(2110, 2251)
    directory = processed(tmp_path, cube)
    before = digest(directory / "stokes.npy")
    result = raptor.fit_processed(directory, tmp_path / "fit", **dict(FIT, fov=2, angular_bins=4),
                                  snapshot_labels=labels)
    for name in ("fits_2d", "fits_1d"):
        series = raptor._FitSeries(result[name], ndim=(3, 4), stokes=False)
        with np.load(result[name], allow_pickle=False) as data:
            np.testing.assert_array_equal(data["snapshot_labels"], labels)
            for t in range(141):
                states = cube[t, -16:]
                evpa = raptor._evpa(states) if name == "fits_2d" else raptor._integrate(states, op)
                expected = raptor.fit_phase_dft(evpa)
                np.testing.assert_array_equal(series.frame(t), expected)
                for column, key in enumerate(("amplitude", "phase", "baseline")):
                    np.testing.assert_array_equal(data[key][t], expected[..., column])
    meta = json.loads(result["metadata"].read_text())
    assert meta["snapshot_labels"] == labels.tolist()
    assert digest(directory / "stokes.npy") == before


@pytest.mark.parametrize("labels", [None, [1600], [1600.0, 1610.0], [[1600, 1610]], [True, False]])
def test_fit_rejects_invalid_snapshot_labels_before_writing(tmp_path, monkeypatch, labels):
    use_operator(monkeypatch, operator(2, bins=4))
    with pytest.raises(ValueError, match="snapshot_labels"):
        raptor.fit_processed(processed(tmp_path, phase_cube(times=2, n=2)), tmp_path / "bad",
                             **dict(FIT, fov=2, angular_bins=4), snapshot_labels=labels)
    assert not (tmp_path / "bad").exists()


def test_preprocess_only_requested_states_retains_axis(tmp_path, monkeypatch):
    cube = phase_cube(times=3)
    original = raptor.preprocess_frame
    calls = []

    def record(frame, **kwargs):
        calls.append(frame)
        return original(frame, **kwargs)

    monkeypatch.setattr(raptor, "preprocess_frame", record)
    result = raptor.preprocess_cube(save(tmp_path / "cube.npy", cube), tmp_path / "bg", sigma=0,
                                    init_fov=3, final_fov=3, rotation_angle=0, final_shape=3,
                                    state_indices=[0])
    data = np.load(result["stokes"])
    assert data.shape == (3, 1, 3, 3, 4)
    np.testing.assert_allclose(data, cube[:, :1], atol=1e-14)
    assert len(calls) == 4  # One tiny geometry validator + exactly T state0 images.
    assert json.loads(result["metadata"].read_text())["state_indices"] == [0]


def test_fit_1d_roll_minus90_matches_eht(tmp_path, monkeypatch):
    n = 20
    amplitude = np.linspace(0.01, 0.15, n * n).reshape(n, n)
    phases = np.linspace(-2, 2, n * n).reshape(n, n)
    psi = (2 * np.pi * np.arange(16) / 16)[:, None, None]
    angles = 0.2 + amplitude * np.cos(psi + phases)
    cube = np.ones((1, 17, n, n, 4))
    cube[0, 1:, ..., 1] = np.cos(2 * angles)
    cube[0, 1:, ..., 2] = np.sin(2 * angles)
    use_operator(monkeypatch, csr_matrix((np.ones(360), (np.arange(360), np.arange(360))),
                                         shape=(360, n * n)))
    result = raptor.fit_processed(processed(tmp_path, cube), tmp_path / "fits", **dict(FIT, fov=4),
                                  snapshot_labels=[1600])
    with np.load(result["fits_1d"], allow_pickle=False) as one:
        np.testing.assert_allclose(one["amplitude"][0], np.roll(amplitude.ravel()[:360], -90),
                                   atol=1e-14)
        np.testing.assert_allclose(one["phase"][0], np.roll(phases.ravel()[:360], -90), atol=1e-14)
    assert json.loads(result["metadata"].read_text())["angular_roll_bins"] == -90


@pytest.mark.parametrize("spatial", [(2, 2), (360,)])
def test_summary_multi_alpha_named_arrays_and_circular_phase(tmp_path, spatial):
    phases = np.array([np.pi - 0.1, -np.pi + 0.1, np.pi - 0.2, -np.pi + 0.2])
    amplitudes = np.array([1, 2, 3, 4])
    fits = np.zeros((4,) + spatial + (3,))
    fits[..., 0] = amplitudes.reshape((4,) + (1,) * len(spatial))
    fits[..., 1] = phases.reshape((4,) + (1,) * len(spatial))
    first = response(tmp_path / "a.npz", fits, np.arange(4))
    second = response(tmp_path / "b.npz", fits, np.arange(4))
    before = {p: digest(p) for p in (first, second)}
    result = raptor.summarize_fits({0.3: first, 0.1: second}, tmp_path / "summary")
    assert result["metadata"].name == "summary.json"
    with np.load(result["templates"], allow_pickle=False) as templates:
        for name in ("alphas", "amplitude_mean", "amplitude_std", "phase_mean", "phase_std"):
            assert name in templates
            assert templates[name].dtype == np.float64
        np.testing.assert_array_equal(templates["alphas"], [0.1, 0.3])
        assert templates["amplitude_mean"].shape == (2,) + spatial
        np.testing.assert_allclose(templates["amplitude_mean"], 2.5)
        np.testing.assert_allclose(templates["amplitude_std"], np.std(amplitudes))
        np.testing.assert_allclose(np.abs(templates["phase_mean"]), np.pi)
        np.testing.assert_allclose(templates["phase_std"], np.sqrt(0.025))
        np.testing.assert_array_equal(templates["phase_count"], 4)
    assert {p: digest(p) for p in before} == before


@pytest.mark.parametrize("problem", ["missing", "shape", "dtype", "labels", "object", "npy"])
def test_summary_rejects_malformed_responses(tmp_path, problem):
    values = {key: np.ones((3, 4)) for key in ("amplitude", "phase", "baseline")}
    values["snapshot_labels"] = np.array([1600, 1610, 1620])
    if problem == "missing":
        del values["baseline"]
    elif problem == "shape":
        values["phase"] = np.ones((3, 5))
    elif problem == "dtype":
        values["phase"] = values["phase"].astype("float32")
    elif problem == "labels":
        values["snapshot_labels"] = np.array([1600])
    elif problem == "object":
        values["baseline"] = values["baseline"].astype(object)
    path = tmp_path / "bad.npz"
    np.savez(path, **values)
    if problem == "npy":  # Bare 3-column arrays are not response files.
        path = save(tmp_path / "bad.npy", np.ones((3, 4, 3)))
    with pytest.raises(ValueError):
        raptor.summarize_fits({0.3: path}, tmp_path / "summary")
    assert not (tmp_path / "summary").exists()


def test_summary_missing_samples(tmp_path):
    fits = np.ones((3, 2, 3))
    fits[:2, 0, :2] = np.nan
    result = raptor.summarize_fits({0.2: response(tmp_path / "r.npz", fits, [1, 2, 3])},
                                   tmp_path / "summary")
    with np.load(result["templates"]) as data:
        assert np.isnan(data["amplitude_mean"][0, 0])
        assert np.isnan(data["phase_std"][0, 0])
        assert data["amplitude_count"][0, 0] == 1
        assert data["amplitude_mean"][0, 1] == 1


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
    grid = raptor.convert_hdf5([[p, p], [p, p]], tmp_path / "grid", **kwargs)
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
    assert np.load(raptor.convert_hdf5([[p]], tmp_path / "ok", **one)["stokes"]).shape == (1, 1, 2, 2, 4)
    with pytest.raises(ValueError, match="missing HDF5"):
        raptor.convert_hdf5([[p]], tmp_path / "bad", **dict(one, datasets={**names, "Q": "missing"}))
    with pytest.raises(ValueError, match="overlapping"):
        raptor.convert_hdf5([[p]], tmp_path / "bad", **dict(one, block_origins=[[0, 0], [0, 0]]))
    with pytest.raises(ValueError, match="cover"):
        raptor.convert_hdf5([[p]], tmp_path / "bad", **dict(one, image_shape=(4, 4)))
    with pytest.raises(ValueError, match="rectangular"):
        raptor.convert_hdf5([[p, p], [p]], tmp_path / "bad", **one)
    assert not (tmp_path / "bad").exists()


def test_validation_paths_dtype_and_shards(tmp_path):
    cube = phase_cube(times=2)
    kwargs = dict(sigma=0, init_fov=3, final_fov=2, rotation_angle=0, final_shape=2)
    path = save(tmp_path / "stokes.npy", cube)
    before = digest(path)
    with pytest.raises(ValueError, match="aliases"):
        raptor.preprocess_cube(path, tmp_path, **kwargs)
    alias = tmp_path / "alias"
    alias.mkdir()
    (alias / "stokes.npy").symlink_to(path)
    with pytest.raises(ValueError, match="aliases"):
        raptor.preprocess_cube(path, alias, **kwargs)
    with pytest.raises(ValueError, match="real numeric"):
        raptor.preprocess_cube(save(tmp_path / "complex.npy", cube.astype(complex)), tmp_path / "bad",
                               **kwargs)
    wrong = save(tmp_path / "wrong.npy", cube.astype("float32"))
    with pytest.raises(ValueError, match="identical"):
        raptor.preprocess_cube([path, wrong], tmp_path / "bad", **kwargs)
    with pytest.raises(ValueError, match="duplicate"):
        raptor.preprocess_cube([path, path], tmp_path / "bad", **kwargs)
    with pytest.raises(ValueError, match="sigma"):
        raptor.preprocess_cube(path, tmp_path / "bad", **{**kwargs, "sigma": -1})
    assert digest(path) == before
    assert not (tmp_path / "bad").exists()


def test_existing_products_are_not_overwritten(tmp_path):
    kwargs = dict(sigma=0, init_fov=3, final_fov=2, rotation_angle=0, final_shape=2)
    source = save(tmp_path / "cube.npy", phase_cube(times=1))
    result = raptor.preprocess_cube(source, tmp_path / "stage", **kwargs)
    before = digest(result["stokes"])
    with pytest.raises(FileExistsError):
        raptor.preprocess_cube(source, tmp_path / "stage", **kwargs)
    assert digest(result["stokes"]) == before
