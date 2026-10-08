"""Synthetic-only tests: no notebook execution and no repository data cubes."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter
from scipy.sparse import csr_matrix

from eht_axion.pipeline import raptor
from scripts import raptor_tools


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, array):
    np.save(path, array)
    return path


def operator(n, bins=360):
    return csr_matrix(np.ones((bins, n * n)))


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


def test_preprocess_matches_scipy_legacy_blur_path():
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


def test_preprocess_shards_memmap_paths_and_dtype(tmp_path, monkeypatch):
    cube = phase_cube(times=5).astype("float32")
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
    params = dict(sigma=0.1, init_fov=3, final_fov=2, rotation_angle=0.3,
                  final_shape=2, dtype="float32")
    one = raptor.preprocess_cube(whole, tmp_path / "a", **params)
    many = raptor.preprocess_cube(shards, tmp_path / "new" / "b", **params)
    first, second = original_load(one["stokes"]), original_load(many["stokes"])
    assert first.shape == (5, 17, 2, 2, 4)
    assert first.dtype == np.float32
    np.testing.assert_array_equal(first, second)
    assert set(loads) == set(before)
    assert {p: digest(p) for p in before} == before
    metadata = json.loads(many["metadata"].read_text())
    assert metadata["inputs"] == [str(p) for p in shards]
    assert many["stokes"] == tmp_path / "new" / "b" / "stokes.npy"
    assert many["metadata"].name == "processing.json"


@pytest.mark.parametrize("states", [16, 17])
@pytest.mark.parametrize("phase", [-3.0, -0.4, 0.7, 3.0])
def test_dft_phase_oracle_and_state0_exclusion(states, phase):
    cube = phase_cube(times=1, states=states, phase=phase)
    angles = 0.5 * np.arctan2(cube[0, ..., 2], cube[0, ..., 1])
    if states == 17:
        angles[0] = np.nan  # State0 even invalid must be entirely irrelevant.
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
    with pytest.raises(ValueError, match="16 phases"):
        raptor.fit_phase_dft(np.ones((15, 2)))


def test_fit_processed_and_fit_shards(tmp_path):
    cube = phase_cube(times=4)
    processed = tmp_path / "processed"
    processed.mkdir()
    save(processed / "stokes.npy", cube)
    shards = [save(tmp_path / "s1.npy", cube[:1]), save(tmp_path / "s2.npy", cube[1:])]
    kwargs = dict(fov=3, operator=operator(3))
    outputs = raptor.fit_processed(processed, tmp_path / "fits", **kwargs)
    sharded = raptor.fit_stokes(shards, tmp_path / "shardfits", **kwargs)
    for name, shape in [("fits_2d", (4, 3, 3)), ("fits_1d", (4, 360))]:
        assert outputs[name].name == name.replace("fits_", "response-") + ".npz"
        with np.load(outputs[name], allow_pickle=False) as data, np.load(sharded[name]) as other:
            assert set(data.files) == {"amplitude", "phase", "baseline", "snapshot_labels"}
            for key, value in [("amplitude", 0.12), ("phase", 0.7), ("baseline", 1.51)]:
                assert data[key].shape == shape
                assert data[key].dtype == np.float64
                np.testing.assert_allclose(data[key], value, atol=2e-14)
                np.testing.assert_array_equal(data[key], other[key])
            np.testing.assert_array_equal(data["snapshot_labels"], [1600, 1610, 1620, 1630])
    assert outputs["metadata"].name == "response.json"
    meta = json.loads(outputs["metadata"].read_text())
    assert meta["weighting"] == "intensity"
    assert meta["state0_used"] is False
    assert meta["snapshot_labels_source"] == "default_convention"
    assert meta["snapshot_labels_verified"] is False


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("provided", [False, True])
def test_response_exact_kernel_columns_and_snapshot_labels(tmp_path, dtype, provided):
    # 141 tiny frames exercise the production label convention, not production data.
    cube = phase_cube(times=141, n=2)
    cube[0, 3, 0, 0, 1:3] = 0  # Invalid phase: NaNs must also survive serialization.
    labels = np.arange(2110, 2251) if provided else None
    processed = tmp_path / "processed"
    processed.mkdir()
    path = save(processed / "stokes.npy", cube)
    before = digest(path)
    result = raptor.fit_processed(processed, tmp_path / "fit", fov=2,
                                  operator=operator(2, bins=4), angular_bins=4,
                                  dtype=dtype, snapshot_labels=labels)
    expected_labels = labels if provided else 1600 + 10 * np.arange(141)
    for name in ("fits_2d", "fits_1d"):
        series = raptor._FitSeries(result[name], ndim=(3, 4), stokes=False)
        with np.load(result[name], allow_pickle=False) as response:
            np.testing.assert_array_equal(response["snapshot_labels"], expected_labels)
            for t in range(141):
                states = cube[t, -16:]
                if name == "fits_2d":
                    evpa = raptor._evpa(states)
                else:
                    evpa = raptor._integrate(states, operator(2, bins=4), "intensity")
                expected = raptor.fit_phase_dft(evpa).astype(dtype)
                np.testing.assert_array_equal(series.frame(t), expected)
                for column, key in enumerate(("amplitude", "phase", "baseline")):
                    assert response[key].dtype == np.dtype(dtype)
                    np.testing.assert_array_equal(response[key][t], expected[..., column])
    meta = json.loads(result["metadata"].read_text())
    assert meta["snapshot_labels"] == expected_labels.tolist()
    assert meta["snapshot_labels_source"] == ("caller" if provided else "default_convention")
    assert meta["snapshot_labels_verified"] is False
    assert digest(path) == before
    assert not (tmp_path / "fit" / "fits_2d.npy").exists()


@pytest.mark.parametrize("labels", [[1600], [1600.0, 1610.0], [[1600, 1610]], [True, False]])
def test_fit_rejects_invalid_snapshot_labels_before_writing(tmp_path, labels):
    with pytest.raises(ValueError, match="snapshot_labels"):
        raptor.fit_stokes(phase_cube(times=2, n=2), tmp_path / "bad", fov=2,
                          operator=operator(2, bins=4), angular_bins=4, snapshot_labels=labels)
    assert not (tmp_path / "bad").exists()


def test_preprocess_only_requested_states_retains_axis(tmp_path, monkeypatch):
    cube = phase_cube(times=3)
    original = raptor.preprocess_frame
    calls = []

    def record(frame, **kwargs):
        calls.append(frame)
        return original(frame, **kwargs)

    monkeypatch.setattr(raptor, "preprocess_frame", record)
    result = raptor.preprocess_cube(cube, tmp_path / "bg", sigma=0, init_fov=3,
                                    final_fov=3, rotation_angle=0, final_shape=3,
                                    state_indices=[0])
    data = np.load(result["stokes"])
    assert data.shape == (3, 1, 3, 3, 4)
    np.testing.assert_allclose(data, cube[:, :1], atol=1e-14)
    assert len(calls) == 4  # One tiny geometry validator + exactly T state0 images.
    assert json.loads(result["metadata"].read_text())["state_indices"] == [0]


def test_fit_1d_roll_minus90_matches_eht(tmp_path):
    n = 20
    amplitude = np.linspace(0.01, 0.15, n * n).reshape(n, n)
    phases = np.linspace(-2, 2, n * n).reshape(n, n)
    psi = (2 * np.pi * np.arange(16) / 16)[:, None, None]
    angles = 0.2 + amplitude * np.cos(psi + phases)
    cube = np.ones((1, 16, n, n, 4))
    cube[0, ..., 1] = np.cos(2 * angles)
    cube[0, ..., 2] = np.sin(2 * angles)
    op = csr_matrix((np.ones(360), (np.arange(360), np.arange(360))), shape=(360, n * n))
    result = raptor.fit_stokes(cube, tmp_path / "fits", fov=4, operator=op)
    with np.load(result["fits_1d"], allow_pickle=False) as one:
        np.testing.assert_allclose(one["amplitude"][0], np.roll(amplitude.ravel()[:360], -90),
                                   atol=1e-14)
        np.testing.assert_allclose(one["phase"][0], np.roll(phases.ravel()[:360], -90), atol=1e-14)
    assert json.loads(result["metadata"].read_text())["angular_roll_bins"] == -90


def test_fit_uses_parent_weighting_and_operator(tmp_path):
    # Spatially varying polarization produces genuinely different definitions.
    cube = phase_cube(times=1, n=2)
    cube[:, :, 0, 0, 0] = 10
    cube[:, :, 0, 0, 2] *= -1
    op = operator(2, bins=4)
    a = raptor.fit_stokes(cube, tmp_path / "intensity", fov=2, operator=op, angular_bins=4)
    b = raptor.fit_stokes(cube, tmp_path / "flux", fov=2, operator=op,
                          angular_bins=4, weighting="flux")
    with np.load(a["fits_1d"]) as first, np.load(b["fits_1d"]) as second:
        assert not np.allclose(first["amplitude"], second["amplitude"])


@pytest.mark.parametrize("spatial", [(2, 2), (360,)])
def test_summary_multi_alpha_named_arrays_and_circular_phase(tmp_path, spatial):
    phases = np.array([np.pi - 0.1, -np.pi + 0.1, np.pi - 0.2, -np.pi + 0.2])
    amplitudes = np.array([1, 2, 3, 4])
    fits = np.zeros((4,) + spatial + (3,))
    fits[..., 0] = amplitudes.reshape((4,) + (1,) * len(spatial))
    fits[..., 1] = phases.reshape((4,) + (1,) * len(spatial))
    p = save(tmp_path / "fits.npy", fits)
    s = [save(tmp_path / "shard1.npy", fits[:1]), save(tmp_path / "shard2.npy", fits[1:])]
    result = raptor.summarize_fits({0.3: s, 0.1: p}, tmp_path / "summary")
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


@pytest.mark.parametrize("spatial", [(2, 2), (4,)])
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_summary_npz_legacy_and_mixed_shards_element_identical(tmp_path, spatial, dtype):
    rng = np.random.default_rng(174)
    fits = rng.normal(size=(5,) + spatial + (3,)).astype(dtype)
    fits[0, ..., 1] = np.nan
    fits[2, ..., 0] = np.nan
    fits[-1, ..., 2] = np.nan  # Baseline is loaded faithfully, not used by summaries.
    legacy = save(tmp_path / "fits.npy", fits)

    def response(path, values, labels):
        np.savez(path, amplitude=values[..., 0], phase=values[..., 1],
                 baseline=values[..., 2], snapshot_labels=labels)
        return path

    named = response(tmp_path / "response.npz", fits, np.arange(5))
    first = response(tmp_path / "first.npz", fits[:2], [42, 43])
    last = save(tmp_path / "last.npy", fits[2:])
    before = {p: digest(p) for p in (legacy, named, first, last)}
    for path in (legacy, named):
        series = raptor._FitSeries(path, ndim=(3, 4), stokes=False)
        for t in range(5):
            np.testing.assert_array_equal(series.frame(t), fits[t])
        if path == legacy:
            assert isinstance(series.arrays[0], np.memmap)
            assert not series.arrays[0].flags.writeable
    outputs = [raptor.summarize_fits({0.3: source}, tmp_path / str(i))
               for i, source in enumerate((legacy, named, [first, last], fits))]
    with np.load(outputs[0]["templates"], allow_pickle=False) as expected:
        for output in outputs[1:]:
            with np.load(output["templates"], allow_pickle=False) as actual:
                assert set(actual.files) == set(expected.files)
                for key in expected.files:
                    np.testing.assert_array_equal(actual[key], expected[key])
    assert {p: digest(p) for p in before} == before


@pytest.mark.parametrize("problem", ["missing", "shape", "dtype", "labels", "object"])
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
    else:
        values["baseline"] = values["baseline"].astype(object)
    path = tmp_path / "bad.npz"
    np.savez(path, **values)
    with pytest.raises(ValueError):
        raptor.summarize_fits({0.3: path}, tmp_path / "summary")
    assert not (tmp_path / "summary").exists()


def test_summary_missing_samples(tmp_path):
    fits = np.ones((3, 2, 3))
    fits[:2, 0, :2] = np.nan
    result = raptor.summarize_fits({0.2: fits}, tmp_path / "summary")
    with np.load(result["templates"]) as data:
        assert np.isnan(data["amplitude_mean"][0, 0])
        assert np.isnan(data["phase_std"][0, 0])
        assert data["amplitude_count"][0, 0] == 1
        assert data["amplitude_mean"][0, 1] == 1


def background_cube(angles):
    angles = np.asarray(angles)
    cube = np.ones((len(angles), 1, 2, 2, 4))
    cube[..., 1] = np.cos(2 * angles)[:, None, None, None]
    cube[..., 2] = np.sin(2 * angles)[:, None, None, None]
    return cube


def test_background_angular_mle_and_diagnostics(tmp_path):
    angles = np.array([1.50, -1.50, -1.45, -1.45, 1.50, 1.48])
    cube = background_cube(angles)
    # Include a truly zero difference, and a zero-polarization pixel to exclude.
    cube[:, :, 0, 0, 1:3] = 0
    paths = [save(tmp_path / "a.npy", cube[:2]), save(tmp_path / "b.npy", cube[2:])]
    result = raptor.estimate_background(paths, tmp_path / "bg", lag=1, bins=9)
    assert result["metadata"].name == "summary.json"
    assert result["diagnostics"].name == "evpa-difference.npz"
    differences = (np.diff(angles) + np.pi / 2) % np.pi - np.pi / 2
    meta = json.loads(result["metadata"].read_text())
    assert isinstance(meta["sigma"], float)
    assert meta["sigma"] == pytest.approx(np.std(differences))
    assert meta["mean"] == pytest.approx(np.mean(differences))
    assert meta["finite_count"] == 15
    assert meta["zero_differences_retained"] is True
    with np.load(result["diagnostics"], allow_pickle=False) as data:
        assert set(data.files) == {"hist_counts", "hist_edges", "count", "mean", "sigma",
                                   "mean_mle", "sigma_mle"}
        assert data["hist_counts"].sum() == 15
        assert data["hist_edges"].shape == (10,)
        assert data["sigma"].shape == ()


def test_background_lag3_constant_zero_and_4d_input(tmp_path):
    cube = background_cube(np.zeros(6))[:, 0]
    result = raptor.estimate_background(cube, tmp_path / "bg", noaxion_state=None)
    meta = json.loads(result["metadata"].read_text())
    assert meta["sigma"] == 0
    assert meta["mean"] == 0
    assert meta["lag_hours"] == pytest.approx(26.68)
    assert meta["exact_day"] is False
    assert meta["finite_count"] == 12
    with np.load(result["diagnostics"]) as data:
        assert data["hist_counts"].sum() == 12


def test_background_mask_linear_and_integrated(tmp_path):
    cube = background_cube([1.5, -1.5, 1.5, -1.5])
    mask = np.array([[True, False], [False, False]])
    mask_path = save(tmp_path / "mask.npy", mask)
    linear = raptor.estimate_background(cube, tmp_path / "linear", lag=1,
                                        mask=mask_path, difference="linear")
    meta = json.loads(linear["metadata"].read_text())
    assert meta["sigma"] == pytest.approx(np.std([-3, 3, -3]))
    assert meta["finite_count"] == 3
    angular = raptor.estimate_background(cube, tmp_path / "angular", lag=1,
                                         mask=mask, domain="1d", fov=2,
                                         operator=operator(2, bins=4), angular_bins=4)
    meta = json.loads(angular["metadata"].read_text())
    assert meta["finite_count"] == 12
    assert meta["sigma"] == pytest.approx(np.std([np.pi - 3, 3 - np.pi, np.pi - 3]))


def test_background_legacy_hist_algorithm(tmp_path):
    rng = np.random.default_rng(83)
    angles = rng.normal(0, 0.05, size=500)
    cube = background_cube(angles)
    result = raptor.estimate_background(cube, tmp_path / "bg", lag=1, method="legacy_hist",
                                        bins=30, legacy_threshold=1)
    meta = json.loads(result["metadata"].read_text())
    assert 0.04 < meta["sigma"] < 0.10
    assert meta["method"] == "legacy_hist"
    with np.load(result["diagnostics"]) as data:
        assert "fit_hist_counts" in data
        assert np.isfinite(data["fit_sigma_error"])


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
    kwargs = dict(datasets=names, image_shape=(4, 4), block_origins=origins,
                  pixel_shape=(2, 2), dtype="float32")
    result = raptor.convert_hdf5(p, tmp_path / "single", **kwargs)
    array = np.load(result["stokes"], allow_pickle=False)
    assert array.dtype == np.float32
    np.testing.assert_array_equal(array, expected)
    grid = raptor.convert_hdf5([[p, p], [p, p]], tmp_path / "grid", **kwargs)
    np.testing.assert_array_equal(np.load(grid["stokes"]), np.broadcast_to(expected, (2, 2, 4, 4, 4)))
    assert digest(p) == before


def test_hdf5_dense_and_bad_mappings(tmp_path):
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "dense.h5"
    with h5py.File(p, "w") as file:
        for c, label in enumerate("IQUV"):
            file.create_dataset(label, data=np.full((3, 3), c))
    names = {c: c for c in "IQUV"}
    result = raptor.convert_hdf5([p, p], tmp_path / "dense", datasets=names, image_shape=(3, 3))
    assert np.load(result["stokes"]).shape == (2, 3, 3, 4)
    with pytest.raises(ValueError, match="missing HDF5"):
        raptor.convert_hdf5(p, tmp_path / "bad", datasets={**names, "Q": "missing"},
                            image_shape=(3, 3))
    with pytest.raises(ValueError, match="overlapping"):
        raptor.convert_hdf5(p, tmp_path / "bad", datasets=names, image_shape=(2, 2),
                            block_origins=[[0, 0], [0, 0]], pixel_shape=(2, 2))
    with pytest.raises(ValueError, match="cover"):
        raptor.convert_hdf5(p, tmp_path / "bad", datasets=names, image_shape=(4, 4),
                            block_origins=[[0, 0]], pixel_shape=(2, 2))
    assert not (tmp_path / "bad").exists()


def test_validation_paths_dtype_and_shards(tmp_path):
    cube = phase_cube(times=2)
    kwargs = dict(sigma=0, init_fov=3, final_fov=2, rotation_angle=0)
    path = save(tmp_path / "stokes.npy", cube)
    before = digest(path)
    with pytest.raises(ValueError, match="aliases"):
        raptor.preprocess_cube(path, tmp_path, **kwargs)
    alias = tmp_path / "alias"
    alias.mkdir()
    (alias / "stokes.npy").symlink_to(path)
    with pytest.raises(ValueError, match="aliases"):
        raptor.preprocess_cube(path, alias, **kwargs)
    with pytest.raises(ValueError, match="output dtype"):
        raptor.preprocess_cube(path, tmp_path / "bad", **kwargs, dtype="int16")
    with pytest.raises(ValueError, match="real numeric"):
        raptor.preprocess_cube(cube.astype(complex), tmp_path / "bad", **kwargs)
    wrong = save(tmp_path / "wrong.npy", cube.astype("float32"))
    with pytest.raises(ValueError, match="identical"):
        raptor.preprocess_cube([path, wrong], tmp_path / "bad", **kwargs)
    with pytest.raises(ValueError, match="duplicate"):
        raptor.preprocess_cube([path, path], tmp_path / "bad", **kwargs)
    with pytest.raises(ValueError, match="sigma"):
        raptor.preprocess_cube(path, tmp_path / "bad", **{**kwargs, "sigma": -1})
    with pytest.raises(ValueError, match="operator"):
        raptor.fit_stokes(path, tmp_path / "bad", fov=2, operator=operator(2))
    with pytest.raises(ValueError, match="boolean"):
        raptor.estimate_background(cube, tmp_path / "bad", lag=1, mask=np.ones((3, 3)))
    assert digest(path) == before
    assert not (tmp_path / "bad").exists()


def test_existing_products_are_not_overwritten(tmp_path):
    kwargs = dict(sigma=0, init_fov=3, final_fov=2, rotation_angle=0, final_shape=2)
    result = raptor.preprocess_cube(phase_cube(times=1), tmp_path / "stage", **kwargs)
    before = digest(result["stokes"])
    with pytest.raises(FileExistsError):
        raptor.preprocess_cube(phase_cube(times=1), tmp_path / "stage", **kwargs)
    assert digest(result["stokes"]) == before


def test_script_subcommands_small_end_to_end(tmp_path):
    script = Path(__file__).resolve().parents[1] / "scripts" / "raptor_tools.py"
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)

    def run(*args):
        process = subprocess.run([sys.executable, str(script), *map(str, args)],
                                 cwd=tmp_path, env=env, capture_output=True, text=True, check=True)
        return json.loads(process.stdout)

    path = save(tmp_path / "input.npy", phase_cube(times=4, n=2))
    processed, fitted, summary, bg = [tmp_path / name for name in ("processed", "fit", "sum", "bg")]
    assert "stokes" in run("preprocess", path, "--output-dir", processed, "--sigma", 0,
                            "--init-fov", 2, "--final-fov", 2, "--rotation-angle", 0,
                            "--final-shape", 2)
    assert "fits_1d" in run("fit", processed / "stokes.npy", "--output-dir", fitted,
                             "--fov", 2, "--angular-bins", 4)
    assert "templates" in run("summarize", "--fit", f"0.3={fitted / 'response-1d.npz'}",
                              "--output-dir", summary)
    assert "diagnostics" in run("background", processed / "stokes.npy", "--output-dir", bg)
    h5py = pytest.importorskip("h5py")
    hdf = tmp_path / "image.h5"
    with h5py.File(hdf, "w") as file:
        for c in "IQUV":
            file.create_dataset(c, data=np.ones((2, 2)))
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"datasets": {c: c for c in "IQUV"}, "image_shape": [2, 2]}))
    assert "stokes" in run("convert", hdf, "--output-dir", tmp_path / "converted",
                            "--mapping", mapping, "--states-per-time", 1)


def test_script_main_background_defaults_and_validation(tmp_path, monkeypatch, capsys):
    seen = {}

    def background(**kwargs):
        seen.update(kwargs)
        return {"metadata": tmp_path / "summary.json"}

    monkeypatch.setattr(raptor_tools, "estimate_background", background)
    args = ["background", "input.npy", "--output-dir", str(tmp_path)]
    assert raptor_tools.main(args) == 0
    assert json.loads(capsys.readouterr().out) == {"metadata": str(tmp_path / "summary.json")}
    assert seen == dict(
        inputs=[Path("input.npy")], output_dir=tmp_path, lag=3, cadence_hours=26.68 / 3,
        noaxion_state=0, mask=None, difference="angular", method="mle", bins=100,
        minimum_samples=2, domain="2d", fov=None, weighting="intensity", radial_max=None,
        window_degrees=10, angular_bins=360, legacy_threshold=0.2,
    )
    assert raptor_tools.main([*args, "--noaxion-state", "none"]) == 0
    assert seen["noaxion_state"] is None
    with pytest.raises(SystemExit) as error:
        raptor_tools.main([*args, "--noaxion-state", "invalid"])
    assert error.value.code == 2
    assert not hasattr(raptor, "main")
    assert not hasattr(raptor, "argparse")
