"""Driver controls use synthetic inputs only; no science products are executed."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "run_joint_day_study", SCRIPTS / "run_joint_day_study.py"
)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def bundle():
    rng = np.random.default_rng(31)
    n, m, f = 2, 7, 9
    means = rng.normal(0, 0.03, (4, n))
    cov = np.tile(np.eye(n)[None] * 0.01, (4, 1, 1))
    states = {
        "alphas": np.array([0.3]),
        "selected_indices": np.arange(n),
        "H2": study.H2,
        "H3": study.H3,
        "observed_daily_local_unwrapped": means,
        "observed_h3_principal": study.H3 @ means,
        "background_h3_principal": rng.normal(0, 0.1, (m, 3, n)),
        "frozen_daily": rng.normal(size=(1, f, 4, n)) + 1j * rng.normal(size=(1, f, 4, n)),
        "evolving_daily": rng.normal(size=(1, m, 4, n)) + 1j * rng.normal(size=(1, m, 4, n)),
        "frozen_weights": np.full(f, 1 / f),
        "evolving_weights": np.full(m, 1 / m),
    }
    obs = {
        "means": means,
        "covariances": cov,
        "selected_indices": np.arange(n),
        "pair_differences": study.H2 @ means,
    }
    legacy = {
        "alphas": np.array([0.3]),
        "signal_coeff": rng.normal(size=(1, 2, n)).astype(complex),
        "background_variance": np.array(0.02),
    }
    return states, obs, legacy


def test_validate_and_shape_rejects():
    s, o, legacy = bundle()
    study.validate_bundle(s, o, legacy)
    s["evolving_daily"] = s["evolving_daily"][..., :1]
    with pytest.raises(ValueError, match="complex"):
        study.validate_bundle(s, o, legacy)


def test_all_cases_use_full_vector_and_no_double_counting():
    s, o, legacy = bundle()
    for h in (study.H2, study.H3):
        for name, kwargs in study.build_cases(s, legacy, 0.3, h):
            assert not ("backgrounds" in kwargs and "background_cov" in kwargs), name
            result = study.calculate_joint_days(
                s["observed_daily_local_unwrapped"],
                o["covariances"],
                h,
                couplings=np.geomspace(1 / 137, 1000, 6),
                phases=np.arange(8) * 2 * np.pi / 8,
                **kwargs,
            )
            assert result["data"].shape == (len(h), 2)
            assert np.isfinite(result["c95"])
            if name == "joint_background_evolving_mixture":
                assert len(kwargs["weights"]) == len(s["evolving_weights"])
                assert "background_cov" not in kwargs


def test_h3_projection_exact_data_covariance_and_cases():
    s, o, legacy = bundle()
    y3, c3 = study.joint_days.contrast_observation(
        s["observed_daily_local_unwrapped"], o["covariances"], study.H3
    )
    y2, c2 = study.joint_days.contrast_observation(
        s["observed_daily_local_unwrapped"], o["covariances"], study.H2
    )
    ix = np.array([0, 1, 4, 5])
    np.testing.assert_array_equal(y3[[0, 2]], y2)
    np.testing.assert_array_equal(c3[np.ix_(ix, ix)], c2)
    assert np.any(c3[:2, 2:4] != 0)
    c2cases = dict(study.build_cases(s, legacy, 0.3, study.H2))
    for name, kwargs in study.build_cases(s, legacy, 0.3, study.H3):
        np.testing.assert_allclose(kwargs["signal_coeff"][:, [0, 2]], c2cases[name]["signal_coeff"])
        if "background_cov" in kwargs:
            np.testing.assert_allclose(
                kwargs["background_cov"][np.ix_(ix, ix)],
                c2cases[name]["background_cov"],
            )


def test_independent_moments_preserve_blocks_and_zero_cross():
    s, _, _ = bundle()
    b = s["background_h3_principal"]
    t = s["evolving_daily"][0, :, 1:]
    w = s["evolving_weights"]
    bi, ti, wi = study.independent_gaussian_states(b, t, w)

    def moments(b, t, w):
        x = np.concatenate(
            [
                b.reshape(len(b), -1),
                t.real.reshape(len(t), -1),
                t.imag.reshape(len(t), -1),
            ],
            axis=1,
        )
        mu = w @ x
        return mu, ((x - mu).T * w) @ (x - mu)

    mu, c = moments(b, t, w)
    mui, ci = moments(bi, ti, wi)
    d = b[0].size
    np.testing.assert_allclose(mu, mui, atol=1e-14)
    np.testing.assert_allclose(c[:d, :d], ci[:d, :d], atol=1e-14)
    np.testing.assert_allclose(c[d:, d:], ci[d:, d:], atol=1e-14)
    np.testing.assert_allclose(ci[:d, d:], 0, atol=1e-14)


def test_legacy_recovery_whole_shared_phase_product():
    from scipy.special import logsumexp

    s, o, legacy = bundle()
    kwargs = dict(study.build_cases(s, legacy, 0.3, study.H2))["legacy_recovery"]
    c = np.geomspace(1 / 137, 1000, 120)
    p = np.arange(64) * 2 * np.pi / 64
    result = study.calculate_joint_days(
        o["means"], o["covariances"], study.H2, couplings=c, phases=p, **kwargs
    )
    # Independent pair Gaussian evaluations, multiplied BEFORE phase average.
    signal = np.real(np.exp(1j * p[:, None, None]) * kwargs["signal_coeff"][0])
    q = np.zeros((len(c), len(p)))
    for k in range(2):
        cov = (
            o["covariances"][2 * k]
            + o["covariances"][2 * k + 1]
            + np.eye(2) * float(legacy["background_variance"])
        )
        residual = o["pair_differences"][k] - c[:, None, None] * signal[None, :, k]
        q += np.einsum("cpi,ij,cpj->cp", residual, np.linalg.inv(cov), residual)
    expected = logsumexp(-0.5 * q, axis=1) - np.log(len(p))
    np.testing.assert_allclose(
        result["loglike"] - result["loglike"].max(),
        expected - expected.max(),
        rtol=1e-12,
        atol=1e-8,
    )
    legacy.update(
        recovery_couplings=c,
        recovery_phases=p,
        recovery_loglike=expected[None],
        recovery_c95=np.array([result["c95"]]),
    )
    study.recovery_check(result, legacy, 0.3)
    legacy["recovery_c95"] *= 1.1
    with pytest.raises(ValueError, match="c95"):
        study.recovery_check(result, legacy, 0.3)


def test_baseline_adapter_reads_existing_products_no_legacy_file(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace

    _, obs, _ = bundle()
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    np.savez(baseline / "observations-1d.npz", **obs)
    opts = {
        "phase_count": 8,
        "coupling_points": 6,
        "coupling_min": 1 / 137,
        "coupling_max": 1000.0,
        "prior": "log-uniform",
        "reference_norm": 0.1,
        "cloud_coefficient": 0.03,
        "cloud_pivot_alpha": 0.05,
        "cloud_exponent": 1.75,
        "time_conversion": 0.1,
        "exposure_hours": 6.0,
        "day_gap_hours": 24.0,
        "pair_separation_hours": 120.0,
        "n_principal": 2,
        "circular_phase_interpolation": True,
        "shared_epoch_uncertainty": False,
    }
    diagnostic = {
        "effective_options": opts,
        "phase_count": 8,
        "alphas": [0.3],
        "couplings": np.geomspace(1 / 137, 1000, 6).tolist(),
        "log_likelihood": [[0.0] * 6],
        "background_variance_radians2": 0.04,
    }
    report = {
        "scientific_profile": {},
        "background_sigma_radians": {"1d": 0.2},
        "results": {"observed-cor_fixed_background-1d": {"c95": [2.0]}},
    }
    (baseline / "report.json").write_text(json.dumps(report))
    (baseline / "observed-cor_fixed_background-1d.json").write_text(json.dumps(diagnostic))
    monkeypatch.setattr(study, "manifest_path", lambda directory: directory / "manifest.json")
    monkeypatch.setattr(study, "verify_stage", lambda *a, **kw: None)
    monkeypatch.setattr(study, "require_declared_outputs", lambda *a: None)
    monkeypatch.setattr(study, "options_from_config", lambda c: SimpleNamespace(**opts))
    monkeypatch.setattr(study, "load_templates", lambda *a: (np.array([0.3]), {}, {}, None))
    coeff = np.array([[1.0 + 2j, 3.0 + 4j], [5.0 + 6j, 7.0 + 8j]])
    phases = np.arange(8) * 2 * np.pi / 8
    signal = np.real(np.exp(1j * phases[:, None, None]) * coeff)
    monkeypatch.setattr(study, "evpa_templates", lambda *a, **kw: (signal.copy(), None))
    monkeypatch.setattr(study, "cloud_amplitude_scaling", lambda *a, **kw: 2.0)
    legacy, consumed, receipts = study.load_baseline(baseline, tmp_path / "templates", "1d", obs)
    np.testing.assert_allclose(legacy["signal_coeff"][0], 2 * coeff)
    assert len(consumed) == 4 and len(receipts) == 2
    assert not any("legacy-" in str(p) for p in consumed)
    obs["means"] = obs["means"] + 0.001
    with pytest.raises(ValueError, match="frozen primary"):
        study.load_baseline(baseline, tmp_path / "templates", "1d", obs)
