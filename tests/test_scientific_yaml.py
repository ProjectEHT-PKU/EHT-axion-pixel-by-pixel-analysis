"""Scientific YAML preserves numeric settings and rejects ambiguous/unsafe input."""

from pathlib import Path
import copy
import pytest
import yaml
from eht_axion.pipeline.workflow import read_config

ROOT = Path(__file__).resolve().parents[1]


def test_scientific_profiles_and_numeric_types():
    research = read_config(ROOT / "configs/research.yaml")
    observed = read_config(ROOT / "configs/observed.yaml")
    expected = copy.deepcopy(research)
    expected["inference"]["observation"] = "observed"
    assert observed == expected
    assert research["physical"]["mass_solar"] == 6.5e9
    assert isinstance(research["physical"]["mass_solar"], float)
    assert research["inference"]["observation"] == "null"
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


def test_unquoted_null_is_not_silently_used(tmp_path):
    source = (ROOT / "configs/research.yaml").read_text()
    assert 'observation: "null"' in source
    p = tmp_path / "config.yaml"
    p.write_text(source.replace('observation: "null"', "observation: null"))
    with pytest.raises(ValueError, match="quoted in YAML"):
        read_config(p)


def test_yml_suffix_and_snapshot_are_supported(tmp_path):
    source = read_config(ROOT / "configs/research.yaml")
    p = tmp_path / "config.yml"
    p.write_text(yaml.safe_dump(source))
    assert read_config(p) == source
