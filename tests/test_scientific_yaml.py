"""Scientific YAML preserves numeric settings and rejects ambiguous/unsafe input."""

from pathlib import Path
import pytest
import yaml
from eht_axion.pipeline.config import read_config

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
        ("physical: [\n", "Invalid scientific YAML"),
        ("1: invalid\n", "keys must be strings"),
        ("physical: {mass_solar: 1, mass_solar: 2}\n", "Duplicate YAML"),
        ("physical: &p {mass_solar: 1}\ndata: {<<: *p, mass_solar: 2}\n", "Duplicate YAML"),
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


@pytest.fixture
def profile():
    return read_config(ROOT / "configs/research.yaml")


def save_config(tmp_path, config):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


@pytest.mark.parametrize(
    "section", ["physical", "data", "image", "aperture", "eht", "background", "inference"]
)
def test_every_section_and_field_is_required_and_unknown_fields_rejected(tmp_path, profile, section):
    original = profile[section]
    for field in original:
        profile[section] = {k: v for k, v in original.items() if k != field}
        with pytest.raises(ValueError, match="Missing or invalid"):
            read_config(save_config(tmp_path, profile))
    profile[section] = dict(original, typo=1)
    with pytest.raises(ValueError, match="Unknown"):
        read_config(save_config(tmp_path, profile))
    profile[section] = None
    with pytest.raises(ValueError, match="Missing or invalid"):
        read_config(save_config(tmp_path, profile))
    del profile[section]
    with pytest.raises(ValueError, match="Missing or invalid"):
        read_config(save_config(tmp_path, profile))


def test_unknown_root_and_non_yaml_suffix(tmp_path, profile):
    profile["input_root"] = "/not/a/scientific/parameter"
    with pytest.raises(ValueError, match="Unknown"):
        read_config(save_config(tmp_path, profile))
    path = tmp_path / "config.toml"
    path.write_text("physical = {}")
    with pytest.raises(ValueError, match="must be YAML"):
        read_config(path)


BOOL_FIELDS = [("eht", "historical_mask"), ("eht", "inherit_2d_mask_in_1d")]
INTEGER_FIELDS = [
    ("image", "final_shape"), ("aperture", "angular_bins"), ("eht", "topset_size"),
    ("eht", "angular_offset"), ("background", "source_resolution"), ("background", "lag_frames"),
    ("inference", "phase_count"), ("inference", "coupling_points"), ("inference", "principal_n"),
]
NUMBER_FIELDS = [
    ("physical", "mass_solar"), ("physical", "distance_mpc"), ("physical", "label_step_tg"),
    ("image", "initial_fov_rg"), ("image", "final_fov_rg"), ("image", "rotation_degrees"),
    ("image", "gaussian_sigma_rg"), ("aperture", "window_degrees"), ("aperture", "radial_max_rg"),
    ("eht", "min_radius_rg"), ("eht", "intensity_fraction"), ("eht", "covariance_ridge"),
    ("inference", "coupling_min"), ("inference", "coupling_max"), ("inference", "exposure_hours"),
    ("inference", "day_gap_hours"), ("inference", "pair_separation_hours"),
]


@pytest.mark.parametrize("section,field", BOOL_FIELDS)
@pytest.mark.parametrize("value", ["false", "true", 0, 1, None, [], {}])
def test_booleans_are_not_coerced(tmp_path, profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match=rf"{section}\.{field}.*boolean"):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("section,field", INTEGER_FIELDS)
@pytest.mark.parametrize("value", [True, False, 4.0, "4", None])
def test_integer_fields_exclude_bool_float_and_strings(tmp_path, profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match=rf"{section}\.{field}.*integer"):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("section,field", NUMBER_FIELDS)
@pytest.mark.parametrize("value", [True, "1.0", None, float("nan"), float("inf"), -float("inf")])
def test_numeric_fields_are_finite_real_numbers(tmp_path, profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match=rf"{section}\.{field}.*finite"):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("section,field,value", [
    ("physical", "mass_solar", 0), ("physical", "distance_mpc", -1),
    ("physical", "label_step_tg", 0), ("image", "initial_fov_rg", 0),
    ("image", "final_fov_rg", 0), ("image", "final_fov_rg", 50),
    ("image", "final_shape", 1), ("image", "gaussian_sigma_rg", -0.1),
    ("aperture", "window_degrees", 0), ("aperture", "window_degrees", 361),
    ("aperture", "radial_max_rg", 0), ("aperture", "radial_max_rg", 17),
    ("aperture", "angular_bins", 4), ("eht", "topset_size", 0),
    ("eht", "min_radius_rg", -1), ("eht", "min_radius_rg", 16.3434),
    ("eht", "intensity_fraction", -0.1), ("eht", "intensity_fraction", 1),
    ("eht", "covariance_ridge", -0.1), ("eht", "covariance_ridge", 1),
    ("eht", "angular_offset", -1), ("eht", "angular_offset", 10),
    ("background", "source_resolution", 200), ("background", "lag_frames", 0),
    ("inference", "phase_count", 1), ("inference", "coupling_points", 1),
    ("inference", "coupling_min", 0), ("inference", "coupling_max", 1 / 137),
    ("inference", "principal_n", 0), ("inference", "exposure_hours", -1),
    ("inference", "day_gap_hours", 0), ("inference", "pair_separation_hours", -1),
])
def test_ranges_and_geometry_cross_constraints(tmp_path, profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("alphas", [
    [], None, "0.1", [True], ["0.1"], [None], [float("nan")], [float("inf")],
    [0], [-0.1], [0.305], [0.10000000000000002], [0.3, 0.1], [0.1, 0.10], [0.1, 0.3, 0.3],
])
def test_alpha_grid_requires_canonical_unique_increasing_numeric_identities(tmp_path, profile, alphas):
    profile["data"]["alphas"] = alphas
    with pytest.raises(ValueError, match="data.alphas"):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("alpha", [True, "0.3", float("nan"), float("inf"), 0, -1, 0.305])
def test_background_alpha_uses_same_identity_contract(tmp_path, profile, alpha):
    profile["background"]["alpha"] = alpha
    with pytest.raises(ValueError, match="background.alpha"):
        read_config(save_config(tmp_path, profile))


@pytest.mark.parametrize("section,field", [
    ("aperture", "weighting"), ("background", "method"), ("inference", "prior"),
])
@pytest.mark.parametrize("value", ["typo", True, None, [], {}])
def test_choices_reject_unknown_and_wrong_types(tmp_path, profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match=rf"{section}\.{field}"):
        read_config(save_config(tmp_path, profile))


def test_existing_mini_profile_is_legal_without_science_changes(tmp_path, profile):
    # Same reduced settings as test_pipeline_end_to_end; no full chain execution.
    profile["data"]["alphas"] = [0.1, 0.3]
    profile["image"]["final_shape"] = 8
    profile["eht"].update(topset_size=40, min_radius_rg=0.0, historical_mask=False, covariance_ridge=0.001)
    profile["background"]["source_resolution"] = 100
    profile["inference"].update(phase_count=16, coupling_points=30)
    assert read_config(save_config(tmp_path, profile)) == profile
    # Boundary values already supported by the numerical APIs stay legal.
    profile["data"]["alphas"] = [0.1]  # Background alpha need not belong to the template grid.
    profile["image"].update(final_shape=2, final_fov_rg=46.2262, gaussian_sigma_rg=0, rotation_degrees=-720)
    profile["aperture"].update(weighting="flux", window_degrees=0.5)
    profile["eht"].update(angular_offset=9, intensity_fraction=0, inherit_2d_mask_in_1d=True)
    profile["background"].update(alpha=0.44, method="legacy_hist")
    profile["inference"].update(prior="uniform", exposure_hours=0, pair_separation_hours=0, principal_n=2)
    assert read_config(save_config(tmp_path, profile)) == profile
