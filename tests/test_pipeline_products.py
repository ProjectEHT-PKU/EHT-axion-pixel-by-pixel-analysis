import numpy as np
import pytest
from eht_axion.pipeline.eht import summarize_evpa
from eht_axion.pipeline.products import build_observations, load_observations, load_templates
from eht_axion.pipeline.common import stage_directory, describe_inputs
from eht_axion.pipeline.inference import _check_lineage
from eht_axion.observations import ObservationData
from eht_axion.limits import calculate_from_products, LimitOptions


def test_lineage_must_bind_current_run_paths(tmp_path):
    a = tmp_path / "run-a.npy"
    b = tmp_path / "run-b.npy"
    np.save(a, np.arange(3))
    np.save(b, np.arange(3))
    # Even equal bytes from a different live run may not silently replace lineage.
    manifest = {"inputs": describe_inputs([b])}
    with pytest.raises(ValueError, match="different run/root"):
        _check_lineage(manifest, [a])
    _check_lineage(manifest, [b])


def test_output_root_cannot_contain_input(tmp_path):
    raw = tmp_path / "raw.npy"
    np.save(raw, np.arange(5))
    with pytest.raises(ValueError, match="contains an input"):
        with stage_directory(tmp_path, "test", {}, [raw], rebuild=True):
            pass
    assert raw.is_file()


@pytest.mark.parametrize("layout", ["legacy", "canonical", "mixed"])
@pytest.mark.parametrize("dimension", ["1d", "2d"])
def test_saved_angular_branch_is_used_for_covariance(tmp_path, layout, dimension, monkeypatch):
    # Broad asymmetric ensemble would change covariance if rewrapped about mean.
    data = np.deg2rad([0, 0, 0, 0, 80, 80, -80])[:, None]
    theta = data if dimension == "1d" else data[:, :, None]
    summary = summarize_evpa(theta)
    stokes = np.zeros((len(data), 1, 1, 4))
    stokes[..., 0] = 1
    stokes[..., 1] = np.cos(2 * data[:, :, None])
    stokes[..., 2] = np.sin(2 * data[:, :, None])
    dates = ("2017-04-05", "2017-04-06", "2017-04-10", "2017-04-11")
    if layout != "canonical":
        np.savez(tmp_path / f"summary_{dimension}.npz",
                 **{k: np.stack([v] * 4) for k, v in summary.items()})
    for k, (day, date) in enumerate(zip((95, 96, 100, 101), dates)):
        if layout == "canonical" or (layout == "mixed" and k % 2 == 0):
            day_dir = tmp_path / date
            day_dir.mkdir()
            np.savez(day_dir / "summary.npz",
                     **{f"evpa_{key}_{dimension}": value for key, value in summary.items()})
            path = day_dir / ("evpa-1d.npy" if dimension == "1d" else "stokes.npy")
        else:
            filename = f"evpa_1d_{day}.npy" if dimension == "1d" else f"stokes_{day}.npy"
            path = tmp_path / filename
        np.save(path, data if dimension == "1d" else stokes)
    if dimension == "2d":
        np.save(tmp_path / "mask.npy", np.ones((1, 1), dtype=bool))
    template = tmp_path / "templates.npz"
    template_shape = (1, *summary["mean"].shape)
    np.savez(
        template,
        alphas=[0.1],
        amplitude_mean=np.full(template_shape, 0.1),
        amplitude_std=np.full(template_shape, 0.01),
        phase_mean=np.full(template_shape, 0.3),
        phase_std=np.full(template_shape, 0.1),
    )
    alphas, amplitude, phase, valid = load_templates(template, dimension)
    np.testing.assert_array_equal(alphas, [0.1])
    np.testing.assert_allclose(amplitude[0.1][0], 0.1)
    np.testing.assert_allclose(phase[0.1][0], 0.3)
    assert valid.all() and valid.shape == summary["mean"].shape
    # Reading/migrating existing products must not copy or rewrite input files.
    inputs_before = {p: (p.stat().st_ino, p.stat().st_size, p.stat().st_mtime_ns)
                     for p in tmp_path.rglob("*") if p.is_file()}
    original_load = np.load
    ensemble_reads = []

    def tracking_load(path, *args, **kwargs):
        result = original_load(path, *args, **kwargs)
        if str(path).endswith(".npy") and not str(path).endswith("mask.npy"):
            assert kwargs.get("mmap_mode") == "r"
            assert isinstance(result, np.memmap)
            ensemble_reads.append(path)
        return result

    monkeypatch.setattr(np, "load", tracking_load)
    out = tmp_path / "observations.npz"
    report = build_observations(tmp_path, template, dimension, out, angular_offset=0)
    assert len(ensemble_reads) == 4
    assert report["sample_counts"] == [7] * 4
    for path, identity in inputs_before.items():
        assert (path.stat().st_ino, path.stat().st_size, path.stat().st_mtime_ns) == identity
    with np.load(out, allow_pickle=False) as z:
        assert set(z.files) == {
            "means", "standard_deviations", "covariances", "pair_differences",
            "pair_covariances", "selected_indices", "dimension", "spatial_shape",
        }
        variance = summary["std"].ravel()[0] ** 2
        np.testing.assert_allclose(z["means"], np.tile(summary["mean"].ravel(), (4, 1)), atol=1e-14)
        np.testing.assert_allclose(z["covariances"], variance, atol=1e-14)
        np.testing.assert_allclose(z["standard_deviations"] ** 2, variance, atol=1e-14)
        np.testing.assert_allclose(z["pair_covariances"], 2 * variance, atol=1e-14)
        np.testing.assert_array_equal(z["pair_differences"], 0)
        np.testing.assert_array_equal(z["selected_indices"], [0])
        np.testing.assert_array_equal(z["spatial_shape"], valid.shape)
        assert z["dimension"] == dimension
    loaded = load_observations(out)
    np.testing.assert_allclose(loaded.covariances, variance, atol=1e-14)
    np.testing.assert_array_equal(loaded.selected_indices, [0])


def test_excluded_nan_template_does_not_invalidate_selected_pixels():
    covariance = np.eye(1) * 0.01
    obs = ObservationData(
        tuple(np.zeros(1) for _ in range(4)),
        tuple(np.ones(1) * 0.1 for _ in range(4)),
        tuple(covariance for _ in range(4)),
        (np.zeros(1), np.zeros(1)),
        (covariance, covariance),
        np.array([0]),
    )
    amplitude = {0.1: np.array([[0.1, np.nan], [0.01, np.nan]])}
    phase = {0.1: np.array([[0.3, np.nan], [0.1, np.nan]])}
    _, cut, _ = calculate_from_products(obs, amplitude, phase, [0.1], options=LimitOptions())
    assert np.isfinite(cut).all()
