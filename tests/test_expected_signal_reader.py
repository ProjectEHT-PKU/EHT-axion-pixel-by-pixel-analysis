import numpy as np
import pytest

from eht_axion.analysis.expected_signal import load_expected_signal


@pytest.mark.parametrize("dimension, shape", [("1d", (3,)), ("2d", (2, 2))])
def test_load_expected_signal_named_fields(tmp_path, dimension, shape):
    template = tmp_path / "expected.npz"
    full = (2, *shape)
    amplitude_std = np.full(full, 0.01)
    amplitude_std[1].flat[0] = np.nan
    values = dict(alphas=[0.1, 0.2], amplitude_mean=np.full(full, 0.1),
                  amplitude_std=amplitude_std, phase_mean=np.full(full, 0.3), phase_std=np.full(full, 0.1))
    np.savez(template, **{f"{dimension}/{k}": v for k, v in values.items()})
    alphas, amplitude, phase, valid = load_expected_signal(template, dimension)
    np.testing.assert_array_equal(alphas, [0.1, 0.2])
    np.testing.assert_allclose(amplitude[0.1][0], 0.1)
    np.testing.assert_allclose(phase[0.2][1], 0.1)
    assert valid.shape == shape and not valid.flat[0] and valid.sum() == valid.size - 1


def test_load_expected_signal_rejects_bad_alpha_grid(tmp_path):
    template = tmp_path / "expected.npz"
    values = np.zeros((2, 3))
    np.savez(template, **{"1d/alphas": [0.2, 0.1], "1d/amplitude_mean": values, "1d/amplitude_std": values,
                          "1d/phase_mean": values, "1d/phase_std": values})
    with pytest.raises(ValueError, match="increasing"):
        load_expected_signal(template, "1d")
