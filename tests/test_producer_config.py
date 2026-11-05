"""Production validation shares analysis science and rejects execution settings."""

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from eht_axion.pipeline.config import read_config
from eht_axion.producer.config import read_production_config, validate_production_config

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def profile():
    return read_production_config(ROOT / "configs/production.yaml")


def save(tmp_path, profile):
    path = tmp_path / "profile.yaml"
    path.write_text(yaml.safe_dump(profile), encoding="utf-8")
    return path


def test_new_benchmark_and_unchanged_analysis_values(profile):
    research = read_config(ROOT / "configs/research.yaml")
    assert research["physical"]["mass_solar"] == 6.5e9
    assert research["image"]["initial_fov_rg"] == 46.2262
    expected = deepcopy(research)
    expected["physical"]["mass_solar"] = 6.5e9
    expected["image"]["initial_fov_rg"] = 46.2261
    assert {k: v for k, v in profile.items() if k != "producer"} == expected
    assert profile["inference"]["principal_n"] == 2
    p = profile["producer"]
    assert p["plasma"]["mass_unit_g"] == 3e25
    assert p["plasma"]["frequency_hz"] == 230e9
    assert p["response"] == {"resolution": 100, "phase_count": 16, "reference_norm": 0.1}
    assert p["alpha_table"] == {"0.10": "low", "0.15": "low", "0.20": "low", "0.25": "low",
                                "0.30": "low", "0.35": "high", "0.40": "high", "0.44": "high"}


def test_reader_roundtrip_and_integration(tmp_path, profile):
    path = save(tmp_path, profile)
    assert read_production_config(path) == profile
    assert read_config(path) == profile
    before = deepcopy(profile)
    validate_production_config(profile)
    assert profile == before


def test_analysis_only_validation_optional_but_production_reader_requires_producer(tmp_path, profile):
    del profile["producer"]
    validate_production_config(profile)
    assert read_config(save(tmp_path, profile)) == profile
    with pytest.raises(ValueError, match="requires producer"):
        read_production_config(save(tmp_path, profile))


def test_no_recursive_analysis_read(monkeypatch, profile):
    def forbidden(*args, **kwargs):
        pytest.fail("in-memory validation must not call read_config")
    monkeypatch.setattr("eht_axion.producer.config.read_config", forbidden)
    validate_production_config(profile)


@pytest.mark.parametrize("section", ["build", "plasma", "snapshots", "response", "background",
                                     "radial_tables", "alpha_table", "grid"])
def test_every_producer_section_and_field_required(tmp_path, profile, section):
    original = deepcopy(profile["producer"][section])
    for field in original:
        profile["producer"][section] = {k: v for k, v in original.items() if k != field}
        with pytest.raises(ValueError):
            read_production_config(save(tmp_path, profile))
    profile["producer"][section] = dict(original, arbitrary_plugin=1)
    with pytest.raises(ValueError, match="unknown"):
        read_production_config(save(tmp_path, profile))
    del profile["producer"][section]
    with pytest.raises(ValueError):
        read_production_config(save(tmp_path, profile))


@pytest.mark.parametrize("section,field,value", [
    ("build", "model", "harm3d"), ("build", "metric", "cks"), ("build", "integrator", "rk2"),
    ("build", "grid", "AMR"), ("build", "polarized", False), ("build", "polarized", 1),
    ("build", "sfc", "true"), ("build", "magnetic_polarity", "minus"),
    ("build", "uniform_output", 0), ("build", "uniform_output", True),
    ("build", "sigma_cut", 0), ("build", "thetae_min", 100),
    ("plasma", "inclination_degrees", -1), ("plasma", "inclination_degrees", 181),
    ("plasma", "r_low", 0), ("plasma", "frequency_hz", -1),
    ("snapshots", "start", -1), ("snapshots", "start", True), ("snapshots", "stop", 1590),
    ("snapshots", "stop", 3001), ("snapshots", "stop", 1000000), ("snapshots", "step", 0),
    ("snapshots", "step", 10.0), ("snapshots", "step", 20), ("snapshots", "start", 1601),
    ("response", "resolution", 101), ("response", "resolution", 20),
    ("response", "resolution", 100.0), ("response", "phase_count", 64),
    ("response", "reference_norm", 0),
    ("background", "resolution", 100),
    ("grid", "nxlone", [96, 48, True]), ("grid", "xprobmin", [0.17, 0, 0]),
    ("grid", "xprobmax", [7.824046010856292, 0.5, 2]), ("grid", "hslope", 0.3),
])
def test_invalid_supported_settings(tmp_path, profile, section, field, value):
    profile["producer"][section][field] = value
    with pytest.raises(ValueError):
        read_production_config(save(tmp_path, profile))


@pytest.mark.parametrize("value", [True, "1.0", None, float("nan"), float("inf"), -float("inf")])
@pytest.mark.parametrize("section,field", [
    ("build", "sigma_cut"), ("build", "thetae_min"), ("build", "thetae_max"),
    ("plasma", "mass_unit_g"), ("plasma", "r_low"), ("plasma", "r_high"),
    ("plasma", "inclination_degrees"), ("plasma", "frequency_hz"), ("plasma", "step_size"),
    ("response", "reference_norm"),
])
def test_finite_numbers_not_coerced(tmp_path, profile, section, field, value):
    profile["producer"][section][field] = value
    with pytest.raises(ValueError, match="finite"):
        read_production_config(save(tmp_path, profile))


@pytest.mark.parametrize("mutation", ["n4", "cadence", "analysis_missing", "root_path", "macros",
                                      "alpha_numeric", "alpha_precision", "alpha_assignment", "spin"])
def test_scope_and_cross_contracts(tmp_path, profile, mutation):
    if mutation == "n4":
        profile["inference"]["principal_n"] = 4
    elif mutation == "cadence":
        profile["physical"]["label_step_tg"] = 1.0
    elif mutation == "analysis_missing":
        del profile["aperture"]["window_degrees"]
    elif mutation == "root_path":
        profile["fluid_root"] = "/tmp/data"
    elif mutation == "macros":
        profile["producer"]["build"]["macros"] = {"SIGMA_CUT": 100}
    elif mutation == "alpha_numeric":
        profile["producer"]["alpha_table"][0.1] = profile["producer"]["alpha_table"].pop("0.10")
    elif mutation == "alpha_precision":
        profile["producer"]["alpha_table"]["0.1"] = profile["producer"]["alpha_table"].pop("0.10")
    elif mutation == "alpha_assignment":
        profile["producer"]["alpha_table"]["0.35"] = "low"
    else:
        profile["producer"]["radial_tables"]["low"]["spin"] = 0.99
    with pytest.raises(ValueError):
        read_production_config(save(tmp_path, profile))


@pytest.mark.parametrize("text", [
    "producer: {}\nproducer: {}\n", "physical: {mass_solar: 1, mass_solar: 2}\n",
    "!!python/object/apply:os.system ['false']", "producer: [\n", "- item\n", "", "null\n",
])
def test_safe_unambiguous_yaml(tmp_path, text):
    path = tmp_path / "bad.yaml"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_production_config(path)


def test_small_subset_is_explicit_and_no_defaults(profile):
    profile["data"]["alphas"] = [0.1, 0.4]
    profile["producer"]["alpha_table"] = {"0.10": "low", "0.40": "high"}
    profile["producer"]["snapshots"] = {"start": 1600, "stop": 1610, "step": 10}
    validate_production_config(profile)
