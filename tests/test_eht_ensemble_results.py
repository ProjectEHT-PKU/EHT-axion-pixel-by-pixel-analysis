"""Small saved-product summary tests: no raw lineage or inference required."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from summarize_eht_ensemble_results import (
    MODES,
    comparisons,
    fit_evidence,
    plot_fits,
    plot_limits,
    run,
)


def toy_report():
    return {
        "results": {
            f"{center}-{mode}-{dim}": {
                "alphas": [0.1, 0.3],
                "c95": [10.0, 1.0],
                "log10_c95": [1.0, 0.0],
                "features": 2,
            }
            for center in ("observed", "null-comparison")
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
    assert rows[0]["features_test"] == 3  # Native differences retain changed N.
    test["results"]["observed-std-1d"]["alphas"] = [0.2, 0.3]
    with pytest.raises(ValueError, match="Alpha"):
        comparisons(base, test)


def test_plots_and_saved_fit_extraction(tmp_path):
    for dim in ("1d", "2d"):
        for mode in MODES:
            diagnostic = {
                "alphas": [0.3],
                "selected_indices": [3, 5],
                "data_radians": [[0.01, -0.02], [0.03, 0.01]],
                "effective_pair_covariances": [np.eye(2).tolist()] * 2,
                "fits": [
                    {
                        "data_quadratic": 0.1,
                        "joint_grid_chi2": 0.08,
                        "joint_grid_coupling": 0.2,
                        "joint_grid_phase_radians": 0.4,
                        "prediction_radians": [[0.0, 0.0], [0.0, 0.0]],
                    }
                ],
            }
            (tmp_path / f"observed-{mode}-{dim}.json").write_text(
                json.dumps(diagnostic)
            )
    rows, diagnostics = fit_evidence(tmp_path)
    assert len(rows) == 6
    assert rows[0]["components"] == 4
    assert rows[0]["null_chi2"] == 0.1
    plot_limits(tmp_path, toy_report())
    plot_fits(tmp_path, diagnostics)
    for name in ("limits-observed-null.png", "alpha030-fits.png"):
        assert (tmp_path / name).stat().st_size > 1000


def test_pending_native_fails_without_creating_summary(tmp_path):
    args = SimpleNamespace(
        **{
            key: tmp_path / key
            for key in ("primary", "gridcheck", "priorcheck", "native", "output_root")
        }
    )
    with pytest.raises((FileNotFoundError, ValueError)):
        run(args)
    assert not args.output_root.exists()
    assert not (tmp_path / "output_root-source-receipts").exists()
