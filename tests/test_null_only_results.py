"""Only the fixed zero-day-difference result type is accepted by current readers."""

import numpy as np
import pytest
from eht_axion.pipeline.layout import load_curves


def write_limits(directory, tag):
    for dimension in ("1d", "2d"):
        np.savez(
            directory / f"limits-{dimension}.npz",
            alphas=[0.1, 0.3],
            cuts=np.zeros((3, 2)),
            modes=["std", "cor", "cor_grmhd"],
            variant="test",
            prior="log-uniform",
            observation=tag,
        )


def test_null_results_are_accepted(tmp_path):
    write_limits(tmp_path, "null")
    alpha, cuts, metadata = load_curves(tmp_path, "1d")
    assert cuts.shape == (3, 2)
    assert metadata["observation"] == "null"


def test_measured_day_difference_results_are_rejected(tmp_path):
    write_limits(tmp_path, "observed")
    with pytest.raises(ValueError, match="Only zero-day-difference"):
        load_curves(tmp_path, "1d")
