"""Small saved-product summary tests: no raw lineage or inference required."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from scripts.analysis.summarize import (
    MODES,
    comparisons,
    fit_evidence,
    plot_fits,
    plot_limits,
    results_text,
    run,
)


def toy_report():
    return {
        "credibility": 0.9,
        "background_sigma_radians": {"1d": 0.02, "2d": 0.03},
        "results": {
            f"{data}-{mode}-{dim}": {
                "alphas": [0.1, 0.3],
                "c95": [10.0, 1.0],
                "log10_c95": [1.0, 0.0],
                "features": 2,
            }
            for data in ("observed", "zero-data")
            for mode in MODES
            for dim in ("1d", "2d")
        }
    }


def test_comparisons_dense_reference_and_selection():
    base = toy_report()
    test = toy_report()
    for row in test["results"].values():
        row["log10_c95"] = [1.004, 0.004]
        row["features"] = 3
    rows = comparisons(base, test)
    assert len(rows) == 24
    assert max(abs(r["delta_dex"]) for r in rows) == pytest.approx(0.004)
    assert rows[0]["ratio"] == pytest.approx(10**0.004)
    assert rows[0]["features_test"] == 3  # changed selections are reported, not hidden
    test["results"]["observed-std-1d"]["alphas"] = [0.2, 0.3]
    with pytest.raises(ValueError, match="Alpha"):
        comparisons(base, test)


def test_plots_and_saved_fit_extraction(tmp_path):
    for dim in ("1d", "2d"):
        for mode in MODES:
            diagnostic = {
                "alphas": [0.1, 0.3],
                "selected_indices": [3, 5],
                "data_radians": [[0.01, -0.02], [0.03, 0.01]],
                "effective_pair_covariances": [np.eye(2).tolist()] * 2,
                "fits": [
                    {
                        "data_quadratic": 0.1,
                        "joint_grid_chi2": chi2,
                        "joint_grid_coupling": 0.2,
                        "joint_grid_phase_radians": 0.4,
                        "prediction_radians": [[0.0, 0.0], [0.0, 0.0]],
                    }
                    for chi2 in (0.09, 0.08)
                ],
            }
            (tmp_path / f"observed-{mode}-{dim}.json").write_text(
                json.dumps(diagnostic)
            )
    rows, diagnostics = fit_evidence(tmp_path)
    assert len(rows) == 12  # every case and every alpha, no fixed alpha
    assert [r["alpha"] for r in rows[:2]] == [0.1, 0.3]
    assert rows[0]["components"] == 4
    assert (rows[0]["chi2_zero_coupling"], rows[1]["chi2_best"]) == (0.1, 0.08)
    plot_limits(tmp_path, toy_report())
    figures = plot_fits(tmp_path, diagnostics)
    assert figures == ["fits-alpha-0.10.png", "fits-alpha-0.30.png"]
    for name in ("limits.png", *figures):
        assert (tmp_path / name).stat().st_size > 1000
    comparison = {"grid_vs_primary": [{"delta_dex": 0.002}, {"delta_dex": -0.003}]}
    text = "\n".join(results_text(toy_report(), comparison, rows, None, figures))
    assert "90% posterior quantiles" in text and "95%" not in text  # configured credibility
    assert "max |delta log10 c95| = 0.003 dex" in text
    assert "| 1d | std | 0.30 | 4 | 0.100 | 0.080 | 0.2 |" in text


def test_pending_source_fails_without_creating_summary(tmp_path):
    args = SimpleNamespace(
        **{key: tmp_path / key for key in ("limits", "gridcheck", "priorcheck", "output_root")},
        covcheck=None,
    )
    with pytest.raises((FileNotFoundError, ValueError)):
        run(args)
    assert not args.output_root.exists()
    assert not (tmp_path / "output_root-source-receipts").exists()
