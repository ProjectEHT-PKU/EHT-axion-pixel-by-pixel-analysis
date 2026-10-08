"""Synthetic-only EHT stage regression tests; no historical data dependency."""

import json
import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from eht_axion.angles import circular_difference
from eht_axion.pipeline.eht import (
    DAYS,
    FIELDS,
    FINAL_OFFSETS,
    RANK_OFFSETS,
    align_metadata,
    build_mask,
    ensemble_covariance,
    ensemble_residuals,
    run_eht_stage,
    shift_spatial,
    stokes_evpa,
    summarize_evpa,
    transform_stokes,
)
from eht_axion.polarization import integrate_stokes


def metadata(n):
    return [dict(zip(FIELDS, (i, 1, 0, 1, 1, 1, 1))) for i in range(n)]


def raw_from_stokes(stokes, order):
    # Inverse of flip axis 0 then clockwise rotation; raw layout is (H,W,N).
    raw = {"meta_data": [metadata(len(stokes))[i] for i in order]}
    for channel, key in enumerate(("I_data", "Q_data", "U_data")):
        values = stokes[order, ..., channel]
        if channel in (1, 2):
            values = -values
        raw[key] = np.rot90(values, k=1, axes=(1, 2)).transpose(1, 2, 0)[::-1]
    return raw


def synthetic_inputs(directory, *, n=7, shape=10):
    directory.mkdir()
    canonical = {}
    for k, day in enumerate(DAYS):
        x, y = np.indices((shape, shape))
        stokes = np.zeros((n, shape, shape, 4))
        theta = (0.2 + 0.018 * x + 0.011 * y + 0.03 * k
                 + np.array([-0.09, 0.12, -0.03, 0.04, -0.18, 0.07, 0.01])[:n, None, None]
                 * (1 + 0.015 * x + 0.02 * y))
        stokes[..., 0] = 1 + 0.04 * x + 0.06 * y + 0.01 * np.arange(n)[:, None, None]
        stokes[..., 1] = np.cos(2 * theta)
        stokes[..., 2] = np.sin(2 * theta)
        order = np.random.default_rng(day).permutation(n)
        raw = raw_from_stokes(stokes, order)
        with (directory / f"230_{day}.pkl").open("wb") as stream:
            pickle.dump(raw, stream)
        canonical[day] = stokes
    return canonical


def tiny_operator(shape, bins=360):
    # Different angular bins sample different pixels; 360-bin roll is testable.
    rows = np.arange(bins)
    columns = (rows * 17 + 3) % (shape * shape)
    return csr_matrix((np.ones(bins), (rows, columns)), shape=(bins, shape * shape))


@pytest.mark.parametrize("offset", [(2, 0), (1, 1), (-4, 4), (-1, 1), (0, 0), (9, 0)])
def test_translation_is_padded_source_sampling_not_roll(offset):
    original = np.arange(2 * 6 * 7, dtype=float).reshape(2, 6, 7)
    shifted = shift_spatial(original, offset)
    for x in range(6):
        for y in range(7):
            sx, sy = x + offset[0], y + offset[1]
            if 0 <= sx < 6 and 0 <= sy < 7:
                np.testing.assert_array_equal(shifted[:, x, y], original[:, sx, sy])
            else:
                assert np.isnan(shifted[:, x, y]).all()
    np.testing.assert_array_equal(original, np.arange(84).reshape(2, 6, 7))


def test_transform_matches_legacy_spatial_map_and_stokes_signs():
    image = np.arange(3 * 3 * 4, dtype=float).reshape(3, 3, 4) + 1
    raw = {key: image + i for i, key in enumerate(("I_data", "Q_data", "U_data", "V_data"))}
    order = np.array([3, 0])
    actual = transform_stokes(raw, order)
    for c, key in enumerate(raw):
        expected = np.stack([np.rot90(raw[key][::-1, :, i], -1) for i in order])
        np.testing.assert_array_equal(actual[..., c], expected * (-1 if c in (1, 2) else 1))
    assert (actual[..., 0] > 0).all()  # I MUST NOT change sign.
    del raw["V_data"]
    np.testing.assert_array_equal(transform_stokes(raw, order)[..., 3], 0)


def test_metadata_matches_common_tuples_not_independent_blind_sorts():
    rows = metadata(4)
    common, order = align_metadata([rows[i] for i in (3, 0, 2, 1)])
    assert common == [tuple(row[field] for field in FIELDS) for row in rows]
    np.testing.assert_array_equal(order, [1, 3, 2, 0])
    _, order = align_metadata([rows[i] for i in (1, 2, 0, 3)], common)
    np.testing.assert_array_equal(order, [2, 0, 1, 3])
    changed = metadata(4)
    changed[-1]["TV"] = 8
    with pytest.raises(ValueError, match="identity mismatch"):
        align_metadata(changed, common)
    with pytest.raises(ValueError, match="unique"):
        align_metadata([rows[0], rows[0]])
    with pytest.raises(ValueError, match="missing"):
        align_metadata([{}])
    rows[0]["TV"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        align_metadata(rows)


def test_mean_branch_then_ordinary_mean_and_population_covariance():
    # Deliberately broad asymmetric data: circular center != ordinary mean.
    theta = np.deg2rad(np.array([[-89, 10], [88, 25], [85, 60], [-75, -65], [25, 30]]))
    summary = summarize_evpa(theta, pixel_batch_size=1)
    center = np.angle(np.exp(2j * theta).mean(axis=0)) / 2
    unwrapped = center + circular_difference(theta, center)
    mean = unwrapped.mean(axis=0)
    np.testing.assert_allclose(summary["center"], center)
    np.testing.assert_allclose(summary["mean"], circular_difference(mean, 0))
    np.testing.assert_allclose(summary["std"], unwrapped.std(axis=0, ddof=0))
    np.testing.assert_allclose(summary["concentration"], np.abs(np.exp(2j * theta).mean(axis=0)))
    residual = ensemble_residuals(theta, summary)
    np.testing.assert_allclose(residual, unwrapped - mean)
    covariance = ensemble_covariance(theta, summary, batch_size=2)
    np.testing.assert_allclose(covariance, residual.T @ residual / len(theta))
    np.testing.assert_allclose(np.diag(covariance), summary["std"] ** 2)
    assert np.max(np.abs(circular_difference(summary["mean"], center))) > 0.01
    # Narrow distribution straddling +-pi/2 must not average toward zero.
    narrow = summarize_evpa(np.deg2rad([[89], [-89], [88], [-88]]))
    assert abs(abs(narrow["mean"][0]) - np.pi / 2) < 1e-12
    assert narrow["std"][0] == pytest.approx(np.deg2rad(np.sqrt(2.5)))
    with pytest.raises(ValueError, match="no historical fit fallback"):
        summarize_evpa(theta, centering="histGaussian")


def test_zero_polarization_and_missing_samples_are_not_finite_summaries():
    stokes = np.zeros((3, 2, 2, 4))
    stokes[..., 1] = 1
    stokes[0, 0, 0, 1] = 0
    stokes[:, 1, 1, 1] = np.nan
    theta = stokes_evpa(stokes)
    summary = summarize_evpa(theta)
    assert np.isnan(theta[0, 0, 0])
    assert np.isnan(summary["mean"][0, 0])
    assert np.isnan(summary["std"][0, 0])
    assert summary["count"][0, 0] == 2
    assert summary["count"][1, 1] == 0
    covariance = ensemble_covariance(theta, summary)
    assert np.isnan(covariance[0]).all()
    assert np.isnan(covariance[-1]).all()


def test_historical_mask_geometry_and_tiny_clipped_bounds():
    intensity = np.ones((4, 64, 64))
    spread = np.ones_like(intensity)
    intensity[2, 12, 12] = 0.1  # strict > 0.1peak
    spread[0, 15, 15] = np.nan
    mask = build_mask(intensity, spread)
    assert not mask[12, 12] and not mask[15, 15]
    assert not mask[:, :5].any()
    assert not mask[48:55, 10:19].any()
    assert not mask[31:33, 31:33].any()
    assert mask[5, 5]
    tiny = build_mask(np.ones((4, 3, 3)), np.ones((4, 3, 3)), min_radius=0,
                      exclude_columns=1, exclusion_rectangles=((2, 55, 2, 19),))
    assert not tiny[:, 0].any() and not tiny[2, 2]
    assert tiny[1, 1]


def test_mini_full_stage_fixed_alignment_selection_outputs_and_roll(tmp_path):
    canonical = synthetic_inputs(tmp_path / "input")
    output = tmp_path / "output"
    operator = tiny_operator(10)
    report = run_eht_stage(tmp_path / "input", output, trust_pickle=True, topset_size=4,
                           historical_mask=False, min_radius=0, batch_size=2,
                           pixel_batch_size=7, operator=operator)
    assert report["relocation"] == "fixed; never refitted"
    assert report["rank_offsets"] == RANK_OFFSETS
    assert report["final_offsets"] == FINAL_OFFSETS
    assert report["weighting"] == "intensity"
    assert not report["inherit_2d_mask"]
    assert report["ensemble_size"] == 7
    assert len(list(output.iterdir())) == 13
    assert not list(output.glob(".eht-work-*"))
    # Independently score the aligned (not original-location) samples.
    summaries, intensities, aligned = [], [], []
    for day in DAYS:
        theta = shift_spatial(stokes_evpa(canonical[day]), RANK_OFFSETS[day])
        aligned.append(theta)
        summaries.append(summarize_evpa(theta))
        intensities.append(shift_spatial(canonical[day][..., 0].mean(axis=0), RANK_OFFSETS[day]))
    rank_mask = build_mask(np.stack(intensities), np.stack([s["std"] for s in summaries]),
                           min_radius=0, historical_mask=False, require_positive_std=True)
    expected_score = sum(np.sum((ensemble_residuals(t, s)[:, rank_mask]
                                / s["std"][rank_mask]) ** 2, axis=1)
                         for t, s in zip(aligned, summaries))
    with np.load(output / "selection.npz", allow_pickle=False) as selection:
        np.testing.assert_allclose(selection["score"], expected_score)
        expected_indices = np.argsort(expected_score, kind="stable")
        np.testing.assert_array_equal(selection["sorted_indices"], expected_indices)
        selected = selection["selected_indices"]
        np.testing.assert_array_equal(selected, expected_indices[:4])
        assert selection["identity"].dtype.kind == "U"
        assert selection["raw_sorted_indices"].shape == (4, 7)
        assert selection["selected_raw_indices"].shape == (4, 4)
        assert json.loads(selection["identity"][3]) == [3, 1, 0, 1, 1, 1, 1]
    mask = np.load(output / "mask.npy")
    assert mask.dtype == bool and mask.shape == (10, 10) and mask.any()
    with np.load(output / "summary_2d.npz") as s2, np.load(output / "summary_1d.npz") as s1:
        assert s2["mean"].shape == s2["std"].shape == (4, 10, 10)
        assert s1["mean"].shape == s1["std"].shape == (4, 360)
        assert s1["covariance"].shape == (4, 360, 360)
        for k, day in enumerate(DAYS):
            stokes = np.load(output / f"stokes_{day}.npy", mmap_mode="r")
            assert isinstance(stokes, np.memmap)
            assert stokes.shape == (4, 10, 10, 4)
            expected = shift_spatial(canonical[day][selected].transpose(0, 3, 1, 2),
                                     FINAL_OFFSETS[day]).transpose(0, 2, 3, 1)
            np.testing.assert_allclose(stokes, expected, equal_nan=True)
            evpa = np.load(output / f"evpa_1d_{day}.npy", mmap_mode="r")
            expected_evpa = np.roll(integrate_stokes(expected, operator, weighting="intensity"),
                                    -90, axis=-1)
            np.testing.assert_allclose(evpa, expected_evpa, equal_nan=True)
            summary = summarize_evpa(stokes_evpa(expected))
            np.testing.assert_allclose(s2["std"][k], summary["std"], equal_nan=True)
            np.testing.assert_allclose(s1["std"][k] ** 2, np.diag(s1["covariance"][k]),
                                       atol=1e-14, equal_nan=True)
    assert json.loads((output / "report.json").read_text())["ddof"] == 0
    with pytest.raises(FileExistsError, match="must be empty"):
        run_eht_stage(tmp_path / "input", output, trust_pickle=True, operator=operator)


def test_stage_metadata_mismatch_is_explicit_and_has_no_complete_report(tmp_path):
    synthetic_inputs(tmp_path / "input")
    path = tmp_path / "input" / "230_101.pkl"
    with path.open("rb") as stream:
        raw = pickle.load(stream)
    raw["meta_data"][0]["TV"] = 99
    with path.open("wb") as stream:
        pickle.dump(raw, stream)
    with pytest.raises(ValueError, match="identity mismatch"):
        run_eht_stage(tmp_path / "input", tmp_path / "out", trust_pickle=True, topset_size=3,
                      historical_mask=False, min_radius=0, operator=tiny_operator(10))
    assert not (tmp_path / "out" / "report.json").exists()
    assert not list((tmp_path / "out").glob(".eht-work-*"))


def test_explicit_mask_inheritance_and_weighting_options(tmp_path):
    synthetic_inputs(tmp_path / "input")
    operator = csr_matrix(np.ones((360, 100)))
    seen = []

    def observing_integrator(stokes, op, weighting):
        seen.append((stokes.copy(), weighting))
        return integrate_stokes(stokes, op, weighting)

    run_eht_stage(tmp_path / "input", tmp_path / "out", trust_pickle=True, topset_size=3,
                  min_radius=0, exclude_columns=2, inherit_2d_mask=True, weighting="flux",
                  operator=operator, integrator=observing_integrator)
    mask = np.load(tmp_path / "out" / "mask.npy")
    assert not mask[:, :2].any()
    assert seen and all(mode == "flux" for _, mode in seen)
    for batch, _ in seen:
        assert (batch[:, ~mask, :] == 0).all()


def test_batch_sizes_do_not_change_scientific_results(tmp_path):
    synthetic_inputs(tmp_path / "input")
    operator = tiny_operator(10)
    for name, batch_size, pixel_batch_size in (("a", 1, 3), ("b", 7, 100)):
        run_eht_stage(tmp_path / "input", tmp_path / name, trust_pickle=True, topset_size=4,
                      historical_mask=False, min_radius=0, operator=operator,
                      batch_size=batch_size, pixel_batch_size=pixel_batch_size)
    for filename in ("selection.npz", "summary_1d.npz", "summary_2d.npz"):
        with np.load(tmp_path / "a" / filename) as a, np.load(tmp_path / "b" / filename) as b:
            assert a.files == b.files
            for key in a.files:
                if a[key].dtype.kind in "USbi":
                    np.testing.assert_array_equal(a[key], b[key])
                else:
                    np.testing.assert_allclose(a[key], b[key], atol=1e-14, equal_nan=True)


def test_zero_polarization_does_not_silently_fallback_to_a_ranking(tmp_path):
    synthetic_inputs(tmp_path / "input")
    for day in DAYS:
        path = tmp_path / "input" / f"230_{day}.pkl"
        with path.open("rb") as stream:
            raw = pickle.load(stream)
        raw["Q_data"][:] = 0
        raw["U_data"][:] = 0
        with path.open("wb") as stream:
            pickle.dump(raw, stream)
    with pytest.raises(ValueError, match="ranking mask is empty"):
        run_eht_stage(tmp_path / "input", tmp_path / "out", trust_pickle=True, topset_size=3,
                      historical_mask=False, min_radius=0, operator=tiny_operator(10))
    assert not (tmp_path / "out" / "report.json").exists()
    assert not (tmp_path / "out" / "selection.npz").exists()


def test_trust_pickle_is_required_before_any_input_or_output_access(tmp_path):
    with pytest.raises(ValueError, match="trust_pickle"):
        run_eht_stage(tmp_path / "nonexistent", tmp_path / "output")
    assert not (tmp_path / "output").exists()
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    result = subprocess.run(
        [sys.executable, "-m", "eht_axion.pipeline.eht", "--input-dir", str(tmp_path / "bad"),
         "--output-dir", str(tmp_path / "output")],
        env=env, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 2
    assert "--trust-pickle is required" in result.stderr
    assert not (tmp_path / "output").exists()


def test_cli_mini_stage_with_real_aperture_operator(tmp_path):
    synthetic_inputs(tmp_path / "input", shape=10)
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    # Very small radial extent makes real quadrature inexpensive. No notebooks.
    result = subprocess.run(
        [sys.executable, "-m", "eht_axion.pipeline.eht", "--input-dir", str(tmp_path / "input"),
         "--output-dir", str(tmp_path / "out"), "--trust-pickle", "--topset-size", "3",
         "--fov", "10", "--radial-max", "0.2", "--no-historical-mask", "--min-radius", "0",
         "--batch-size", "2"],
        env=env, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["topset_size"] == 3
    assert np.load(tmp_path / "out" / "evpa_1d_101.npy").shape == (3, 360)
