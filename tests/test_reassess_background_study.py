"""Synthetic-only checks: no completed data products or science run required."""

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import reassess_background_study as study


def synthetic():
    rng = np.random.default_rng(903)
    m, n, a, f = 40, 2, 2, 7
    starts = np.arange(m) * 60.0 + 3
    times = np.arange(0.0, starts[-1] + 148.0, 3.0)
    means = rng.normal(size=(4, n))
    states = dict(
        background_h3_principal=rng.normal(size=(m, 3, n)),
        observed_daily_local_unwrapped=means,
        observed_h3_principal=study.H3 @ means,
        selected_indices=np.arange(n),
        H2=study.H2.copy(),
        H3=study.H3.copy(),
        day_hours=study.DAY_HOURS.copy(),
        window_starts_hours=starts,
        window_centers_hours=starts[:, None] + study.DAY_HOURS,
        snapshot_times_hours=times,
        alphas=np.array([0.2, 0.3]),
        frozen_daily=rng.normal(size=(a, f, 4, n)) + 1j * rng.normal(size=(a, f, 4, n)),
        frozen_weights=np.full(f, 1 / f),
    )
    obs = dict(
        means=means,
        covariances=np.tile(np.eye(n), (4, 1, 1)),
        pair_differences=study.H2 @ means,
        selected_indices=np.arange(n),
    )
    return states, obs


def test_bundle_and_branch_guard():
    s, o = synthetic()
    study.validate_inputs(s, o)
    o["pair_differences"] = o["pair_differences"] + 0.1
    with pytest.raises(ValueError, match="pair branch"):
        study.validate_inputs(s, o)


def test_predefined_cases_whole_states_and_cartesian_mean():
    s, _ = synthetic()
    fits = study.fit_backgrounds(s, [0, 2])
    cases = list(study.build_cases(s, 0.3, study.H2, fits))
    assert len(cases) == 8
    assert tuple(dict.fromkeys(bg for bg, _, _ in cases)) == study.INFERENCE_MODELS
    f = np.einsum("kd,mdn->mkn", study.H2, s["frozen_daily"][1])
    for bg, response, kwargs in cases:
        assert "backgrounds" not in kwargs and "background_mean" not in kwargs
        np.testing.assert_array_equal(kwargs["background_cov"], fits[bg]["covariance"])
        if response == "complex_mean":
            np.testing.assert_allclose(kwargs["signal_coeff"], f.mean(axis=0)[None])
        else:
            np.testing.assert_array_equal(kwargs["signal_coeff"], f)
            np.testing.assert_allclose(kwargs["weights"], 1 / len(f))


def test_iso_h2_old_pool_h3_natural_pool():
    s, _ = synthetic()
    s["background_h3_principal"][:, 1] *= 7
    fits2 = study.fit_backgrounds(s, [0, 2])
    fits3 = study.fit_backgrounds(s, [0, 1, 2])
    _, old = study.joint_days.background_moments(
        s["background_h3_principal"][:, [0, 2]], covariance="scalar"
    )
    np.testing.assert_allclose(fits2["iso"]["covariance"], old)
    assert fits3["iso"]["covariance"][0, 0] > old[0, 0]


def test_purged_validation_train_only_and_counts():
    s, _ = synthetic()
    for count in (5, 3):
        report = study.evaluate_background(s, [0, 2], count)
        assert len(report["records"]) == count * 6
        assert report["selected_model"] is None
        for i, fold in enumerate(report["folds"]):
            support = fold["support_indices"]
            train, test = fold["train"], fold["test"]
            for t in train:
                assert all(
                    support[t, 1] < support[j, 0] or support[t, 0] > support[j, 1]
                    for j in test
                )
            record = next(
                r for r in report["records"] if r["fold"] == i and r["model"] == "iso"
            )
            assert record["n_train"] == len(train)
            assert record["dimension"] == 4
            np.testing.assert_allclose(
                record["train_mean"],
                s["background_h3_principal"][train][:, [0, 2]].mean(axis=0),
            )


def test_no_covariance_rescue_and_rawfull_undefined():
    with pytest.raises(ValueError, match="zero background sigma"):
        study.check_covariance(np.zeros((2, 2)), {"scoreable": False}, "iso", 5, 2)
    assert not study.check_covariance(np.eye(3), {"scoreable": True}, "full", 3, 3)
    with pytest.raises(ValueError, match="non-SPD"):
        study.check_covariance(
            np.diag([1.0, 0]), {"scoreable": False}, "diagonal", 5, 2
        )


def test_positive_inference_with_fixed_background():
    s, obs = synthetic()
    fits = study.fit_backgrounds(s, [0, 2])
    _, _, kwargs = next(study.build_cases(s, 0.3, study.H2, fits))
    result = study.joint_days.calculate_joint_days(
        obs["means"],
        obs["covariances"],
        study.H2,
        couplings=np.geomspace(0.01, 10, 12),
        phases=np.arange(8) * 2 * np.pi / 8,
        **kwargs,
    )
    assert 0.01 <= result["c95"] <= 10
    assert np.all(np.linalg.eigvalsh(result["covariance"]) > 0)
    assert result["diagnostics"]["state_count"] == 1


def test_old_comparison_and_json_undefined(tmp_path):
    path = tmp_path / "old.csv"
    row = dict(
        dimension="1d",
        alpha=0.3,
        contrasts=2,
        mode="observed",
        case="mean_background_iso",
        c95=2.0,
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    new = dict(row, key="example", background="iso", response="complex_mean")
    assert study.compare_old([new], path)[0]["log10_c95_ratio"] == 0
    with pytest.raises(ValueError, match="recovery"):
        study.compare_old([dict(new, c95=3.0)], path)
    assert study.jsonable({"condition": np.inf, "rank": np.int64(2)}) == {
        "condition": None,
        "rank": 2,
    }
