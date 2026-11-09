"""Footer-time cadence contracts with tiny JSON fixtures, no scientific runs.

Identity/ancestry validation belongs to generated_inputs and its delivery tests;
these tests isolate input_cadence_hours, which only reads footer-time metadata.
"""

from copy import deepcopy
import json

import pytest

from eht_axion.experiment import read_experiment
from eht_axion.raptor.delivery import input_cadence_hours


@pytest.fixture
def config():
    return read_experiment()


def tg_hours(config):
    # Independent SI conversion, rather than the helper's PhysicalScale result.
    return 1.3271244e20 * config["shared"]["mass_solar"] / 299792458.0**3 / 3600


def footer_input(tmp_path, name, times):
    directory = tmp_path / name
    directory.mkdir()
    path = directory / "stokes.npy"
    # Deliberately not an array: cadence must not load the Stokes science data.
    path.write_bytes(b"cadence fixture, not scientific Stokes data")
    (directory / "generated.json").write_text(json.dumps({
        "snapshot_labels": [1600 + 10 * i for i in range(len(times))],
        "times_tg": times,
    }))
    return path


def test_three_footer_frames_and_changed_spacing_propagate(tmp_path, config):
    path = footer_input(tmp_path, "generated", [15010, 15012.5, 15015])
    before = deepcopy(config)
    original = input_cadence_hours([path], config)
    assert original == pytest.approx(2.5 * tg_hours(config))
    # Keep labels identical; only real footer times change, not the label grid.
    metadata = path.parent / "generated.json"
    record = json.loads(metadata.read_text())
    record["times_tg"] = [15010, 15015, 15020]
    metadata.write_text(json.dumps(record))
    changed = input_cadence_hours([path], config)
    assert changed == pytest.approx(5 * tg_hours(config))
    assert changed == pytest.approx(2 * original)
    assert config == before


def test_generated_cadence_uses_mass_not_label_step_or_distance(tmp_path, config):
    path = footer_input(tmp_path, "generated", [15010, 15013, 15016])
    original = input_cadence_hours([path], config)
    config["shared"]["mass_solar"] *= 1.25
    # The cadence helper is not a producer validator: label conversion is never
    # used for footer times, even if this analysis-only value differs from 0.1.
    config["shared"]["label_step_tg"] = 0.2
    config["shared"]["distance_mpc"] *= 0.8
    assert input_cadence_hours([path], config) == pytest.approx(1.25 * original)


def test_footer_grid_spans_multiple_shards(tmp_path, config):
    paths = [footer_input(tmp_path, "first", [15010, 15012]),
             footer_input(tmp_path, "second", [15014])]
    assert input_cadence_hours(paths, config) == pytest.approx(2 * tg_hours(config))


@pytest.mark.parametrize("times", [
    [15010, 15012, 15015],  # Nonuniform.
    [15010, 15010, 15010],  # Zero spacing.
    [15010, 15009, 15008],  # Reversed.
    [15010, float("nan"), 15012],
    [15010, float("inf"), 15012],
    [15010, -float("inf"), 15012],
    [], [15010],  # No interval can be inferred.
])
def test_invalid_footer_time_grid_fails(tmp_path, config, times):
    path = footer_input(tmp_path, "generated", times)
    with pytest.raises(ValueError, match="uniform finite footer-time grid"):
        input_cadence_hours([path], config)


def test_nonuniform_shard_boundary_fails(tmp_path, config):
    paths = [footer_input(tmp_path, "first", [15010, 15012]),
             footer_input(tmp_path, "second", [15015, 15017])]
    with pytest.raises(ValueError, match="uniform finite footer-time grid"):
        input_cadence_hours(paths, config)


def test_missing_footer_times_fails_without_label_fallback(tmp_path, config):
    path = footer_input(tmp_path, "generated", [15010, 15012, 15014])
    metadata = path.parent / "generated.json"
    record = json.loads(metadata.read_text())
    del record["times_tg"]
    metadata.write_text(json.dumps(record))
    with pytest.raises((KeyError, ValueError), match="times_tg|footer"):
        input_cadence_hours([path], config)


def test_missing_generated_metadata_fails(tmp_path, config):
    paths = [footer_input(tmp_path, "first", [15010, 15012]),
             footer_input(tmp_path, "second", [15014])]
    (paths[1].parent / "generated.json").unlink()
    with pytest.raises(FileNotFoundError):
        input_cadence_hours(paths, config)
