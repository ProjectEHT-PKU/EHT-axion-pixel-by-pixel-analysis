import numpy as np
import pytest
from eht_axion.pipeline.eht import summarize_evpa
from eht_axion.pipeline.products import build_observations
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


def test_saved_angular_branch_is_used_for_covariance(tmp_path):
    data = np.deg2rad([0, 0, 0, 0, 80, 80, -80])[:, None]
    summary = summarize_evpa(data)
    np.savez(tmp_path / "summary_1d.npz", **{k: np.stack([v] * 4) for k, v in summary.items()})
    for day in (95, 96, 100, 101):
        np.save(tmp_path / f"evpa_1d_{day}.npy", data)
    template = tmp_path / "templates.npz"
    np.savez(
        template,
        alphas=[0.1],
        amplitude_mean=[[0.1]],
        amplitude_std=[[0.01]],
        phase_mean=[[0.3]],
        phase_std=[[0.1]],
    )
    out = tmp_path / "observations.npz"
    build_observations(tmp_path, template, "1d", out, angular_offset=0)
    with np.load(out) as z:
        np.testing.assert_allclose(z["means"][0], summary["mean"], atol=1e-14)
        np.testing.assert_allclose(z["covariances"][0, 0, 0], summary["std"][0] ** 2)


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
