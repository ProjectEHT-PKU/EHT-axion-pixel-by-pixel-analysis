"""Saved observation reuse preserves the observed data, the covariance and the zero-data control."""
from pathlib import Path
import numpy as np
import pytest

from eht_axion.analysis.results import read_results, write_results
from scripts.analysis.check_grid import load_observation


def test_saved_observation_roundtrip(tmp_path):
    means = np.arange(12.).reshape(4, 3)/100
    cov = np.repeat(np.eye(3)[None], 4, axis=0)
    pairs = np.array([cov[0]+cov[1], cov[2]+cov[3]])
    diff = np.array([[.1, -.2, .3], [-.1, .2, -.3]])
    selected = np.array([0, 5, 9])
    write_results(tmp_path, {}, dict(observations={'1d': dict(
        means=means, standard_deviations=np.ones((4, 3)), covariances=cov, pair_differences=diff,
        pair_covariances=pairs, selected_indices=selected, selection={'stride': 1})}))
    _, records = read_results(tmp_path, ('observations',))
    assert records['observations']['1d']['selection'] == {'stride': 1}
    observed, zero_data = load_observation(records['observations']['1d'])
    np.testing.assert_array_equal(observed.pair_differences, diff)
    np.testing.assert_array_equal(zero_data.pair_differences, np.zeros_like(diff))
    np.testing.assert_array_equal(observed.pair_covariances, zero_data.pair_covariances)
    np.testing.assert_array_equal(observed.selected_indices, selected)
    np.testing.assert_array_equal(observed.means, means)


def test_missing_product_is_not_reconstructed(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_results(Path(tmp_path), ('observations',))
