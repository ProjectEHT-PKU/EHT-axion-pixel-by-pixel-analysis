"""Scientific YAML preserves numeric settings and rejects ambiguous/unsafe input."""

from pathlib import Path
import pytest
import yaml
from eht_axion.pipeline.workflow import read_config

ROOT = Path(__file__).resolve().parents[1]


def test_scientific_profiles_and_numeric_types():
    research = read_config(ROOT / "configs/research.yaml")
    assert not (ROOT / "configs/observed.yaml").exists()
    assert research["physical"]["mass_solar"] == 6.5e9
    assert isinstance(research["physical"]["mass_solar"], float)
    assert "observation" not in research["inference"]
    assert research["inference"]["coupling_min"] == 1 / 137
    assert research["inference"]["phase_count"] == 64
    assert research["eht"]["historical_mask"] is True


@pytest.mark.parametrize(
    "content,match",
    [
        ("", "must contain a mapping"),
        ("- invalid\n- root\n", "must contain a mapping"),
        ("physical: {}\nphysical: {}\n", "Duplicate YAML"),
        ("physical: {}\n", "Missing or invalid"),
        ('!!python/object/apply:os.system ["echo unsafe"]', "Invalid scientific YAML"),
    ],
)
def test_reject_invalid_duplicate_or_executable_yaml(tmp_path, content, match):
    path = tmp_path / "config.yaml"
    path.write_text(content)
    with pytest.raises(ValueError, match=match):
        read_config(path)


@pytest.mark.parametrize(
    "key,value", [("observation", "observed"), ("observation", "null"), ("observed", True)]
)
def test_removed_observation_option_is_rejected(tmp_path, key, value):
    source = read_config(ROOT / "configs/research.yaml")
    source["inference"][key] = value
    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump(source))
    with pytest.raises(ValueError, match="Only zero EVPA day differences"):
        read_config(p)


def test_yml_suffix_and_snapshot_are_supported(tmp_path):
    source = read_config(ROOT / "configs/research.yaml")
    p = tmp_path / "config.yml"
    p.write_text(yaml.safe_dump(source))
    assert read_config(p) == source
