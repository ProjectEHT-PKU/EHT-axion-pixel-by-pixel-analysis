import numpy as np
import pytest

from eht_axion.pipeline.common import stage_directory
from eht_axion.pipeline.products import load_templates


def test_output_root_cannot_contain_input(tmp_path):
    raw = tmp_path / "raw.npy"
    np.save(raw, np.arange(5))
    with pytest.raises(ValueError, match="contains an input"):
        with stage_directory(tmp_path, "test", {}, [raw]):
            pass
    assert raw.is_file()


@pytest.mark.parametrize("dimension, shape", [("1d", (3,)), ("2d", (2, 2))])
def test_load_templates_named_fields(tmp_path, dimension, shape):
    template = tmp_path / "templates.npz"
    full = (2, *shape)
    amplitude_std = np.full(full, 0.01)
    amplitude_std[1].flat[0] = np.nan
    np.savez(template, alphas=[0.1, 0.2], amplitude_mean=np.full(full, 0.1),
             amplitude_std=amplitude_std, phase_mean=np.full(full, 0.3), phase_std=np.full(full, 0.1))
    alphas, amplitude, phase, valid = load_templates(template, dimension)
    np.testing.assert_array_equal(alphas, [0.1, 0.2])
    np.testing.assert_allclose(amplitude[0.1][0], 0.1)
    np.testing.assert_allclose(phase[0.2][1], 0.1)
    assert valid.shape == shape and not valid.flat[0] and valid.sum() == valid.size - 1


def test_load_templates_rejects_bad_alpha_grid(tmp_path):
    template = tmp_path / "templates.npz"
    values = np.zeros((2, 3))
    np.savez(template, alphas=[0.2, 0.1], amplitude_mean=values, amplitude_std=values,
             phase_mean=values, phase_std=values)
    with pytest.raises(ValueError, match="increasing"):
        load_templates(template, "1d")
