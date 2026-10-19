"""Derived scales use independent formulae, not historical geometry constants.

Synthetic/unit tests only: no scientific inputs are read and no products run.
"""

from copy import deepcopy
import math
from pathlib import Path

import pytest
import yaml

from eht_axion.pipeline.config import _validate, read_config
from eht_axion.pipeline.constants import PhysicalScale, analysis_geometry
from eht_axion.producer.config import read_production_config, validate_production_config
from eht_axion.producer.delivery import _identity

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def profile():
    return read_config(ROOT / "configs/analysis.yaml")


def independent_scales(physical):
    # IAU nominal solar GM, exact c, parsec in metres; angle in microarcseconds.
    gm = 1.3271244e20 * physical["mass_solar"]
    theta = gm / 299792458.0**2 / (physical["distance_mpc"] * 1e6 * 3.085677581491367e16)
    return theta * (180 / math.pi) * 3600 * 1e6, gm / 299792458.0**3


def test_angular_geometry_and_time_from_independent_formula(profile):
    before = deepcopy(profile)
    result = analysis_geometry(profile)
    theta, seconds = independent_scales(profile["physical"])
    image, aperture, metadata = (result[k] for k in ("image_preprocess", "aperturefit", "metadata"))
    assert image == pytest.approx(dict(
        sigma=profile["image"]["gaussian_fwhm_uas"] / theta / math.sqrt(8 * math.log(2)), init_fov=46.2261,
        final_fov=128 / theta, final_shape=64, rotation_angle=108 * math.pi / 180,
    ))
    assert aperture == dict(
        fov=image["final_fov"], radial_max=image["final_fov"] / 2,
        weighting="intensity", window_degrees=10, angular_bins=360,
    )
    assert metadata["angular_rg_uas"] == pytest.approx(theta)
    assert metadata["tg_seconds"] == pytest.approx(seconds)
    assert metadata["tg_hours"] == pytest.approx(seconds / 3600)
    assert metadata["label_step_hours"] == pytest.approx(0.1 * seconds / 3600)
    assert metadata["omega_conversion_per_hour"] == pytest.approx(3600 / seconds)
    assert metadata["axion_mass_per_alpha_ev"] == pytest.approx(6.582119569e-16 / seconds)
    assert metadata["final_fov_uas"] == pytest.approx(128)
    assert metadata["gaussian_fwhm_uas"] == pytest.approx(profile["image"]["gaussian_fwhm_uas"])
    assert metadata["pixel_size_uas"] == pytest.approx(128 / 64)
    assert metadata["radial_max_uas"] == pytest.approx(64)
    assert metadata["image_mode"] == "angular"
    assert metadata["aperture_mode"] == "fraction"
    assert profile == before


@pytest.mark.parametrize("mass_factor,distance_factor", [(1.25, 1), (1, 0.8), (1.2, 1.1)])
def test_mass_distance_changes_propagate(profile, mass_factor, distance_factor):
    baseline = analysis_geometry(profile)
    profile["physical"]["mass_solar"] *= mass_factor
    profile["physical"]["distance_mpc"] *= distance_factor
    result = analysis_geometry(profile)
    theta, seconds = independent_scales(profile["physical"])
    assert result["metadata"]["angular_rg_uas"] == pytest.approx(theta)
    assert result["metadata"]["tg_seconds"] == pytest.approx(seconds)
    assert result["metadata"]["tg_hours"] == pytest.approx(baseline["metadata"]["tg_hours"] * mass_factor)
    assert result["metadata"]["angular_rg_uas"] == pytest.approx(
        baseline["metadata"]["angular_rg_uas"] * mass_factor / distance_factor)
    for field in ("sigma", "final_fov"):
        assert result["image_preprocess"][field] == pytest.approx(
            baseline["image_preprocess"][field] * distance_factor / mass_factor)
    assert result["image_preprocess"]["init_fov"] == 46.2261
    assert result["aperturefit"]["radial_max"] == pytest.approx(64 / theta)


@pytest.mark.parametrize("filename", ["research.yaml", "production.yaml"])
def test_legacy_rg_geometry_is_identical(filename):
    config = read_config(ROOT / "configs" / filename)
    before = deepcopy(config)
    image, aperture = config["image"], config["aperture"]
    result = analysis_geometry(config)
    assert result["image_preprocess"] == dict(
        sigma=image["gaussian_sigma_rg"], init_fov=image["initial_fov_rg"],
        final_fov=image["final_fov_rg"], final_shape=image["final_shape"],
        rotation_angle=image["rotation_degrees"] * math.pi / 180,
    )
    assert result["aperturefit"] == dict(
        fov=image["final_fov_rg"], weighting=aperture["weighting"],
        radial_max=aperture["radial_max_rg"], window_degrees=aperture["window_degrees"],
        angular_bins=aperture["angular_bins"],
    )
    assert config == before


def test_profile_changes_only_geometry_and_preserves_generated_identity(profile):
    production = read_production_config(ROOT / "configs/production.yaml")
    assert profile.keys() - production.keys() == {"eht_input", "cloud"}
    for section in production.keys() - {"image", "aperture"}:
        assert profile[section] == production[section]
    theta, _ = independent_scales(profile["physical"])
    preserved_fwhm = production["image"]["gaussian_sigma_rg"] * theta * math.sqrt(8 * math.log(2))
    assert profile["image"] == dict(
        initial_fov_rg=46.2261, final_fov_uas=128.0, final_shape=64,
        rotation_degrees=108.0, gaussian_fwhm_uas=pytest.approx(preserved_fwhm),
    )
    assert profile["image"] != production["image"]  # Analysis grid differs, generation identity does not.
    assert analysis_geometry(profile)["image_preprocess"]["sigma"] == pytest.approx(2.16871)
    assert profile["eht_input"]["shape"] == profile["image"]["final_shape"] == 64
    assert profile["eht_input"]["fov_uas"] == profile["image"]["final_fov_uas"] == 128
    assert profile["eht_input"]["pixel_size_uas"] == 2
    assert profile["cloud"] == dict(
        coefficient=0.03, pivot_alpha=0.05, exponent=1.75, reference_norm=0.1,
    )
    assert {k: v for k, v in profile["aperture"].items() if k != "radial_max_fraction"} == {
        k: v for k, v in production["aperture"].items() if k != "radial_max_rg"}
    params = dict(alpha=0.3, resolution=100, labels=[1600, 1610], times=[160, 161], state_count=17)
    assert _identity(profile, **params) == _identity(production, **params)
    assert read_production_config(ROOT / "configs/analysis.yaml") == profile


def test_producer_validation_without_recursive_read(profile, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("in-memory validation must not re-read a configuration")

    monkeypatch.setattr("eht_axion.producer.config.read_config", forbidden)
    validate_production_config(profile)
    del profile["producer"]
    validate_production_config(profile)
    assert analysis_geometry(profile)["metadata"]["image_mode"] == "angular"


@pytest.mark.parametrize("fraction", [0.25, 0.6, 1.0])
def test_fraction_is_independent_aperture_choice(profile, fraction):
    profile["aperture"]["radial_max_fraction"] = fraction
    original = analysis_geometry(profile)
    expected = fraction * original["image_preprocess"]["final_fov"] / 2
    assert original["aperturefit"]["radial_max"] == expected
    profile["image"]["gaussian_fwhm_uas"] = 0
    assert analysis_geometry(profile)["aperturefit"]["radial_max"] == expected
    profile["aperture"].pop("radial_max_fraction")
    profile["aperture"]["radial_max_rg"] = expected
    assert analysis_geometry(profile)["aperturefit"] == original["aperturefit"]


def test_fraction_with_legacy_image_and_zero_kernel():
    config = read_config(ROOT / "configs/research.yaml")
    config["aperture"].pop("radial_max_rg")
    config["aperture"]["radial_max_fraction"] = 0.5
    config["image"]["gaussian_sigma_rg"] = 0
    result = analysis_geometry(config)
    assert result["aperturefit"]["radial_max"] == config["image"]["final_fov_rg"] / 4
    assert result["image_preprocess"]["sigma"] == 0


@pytest.mark.parametrize("section,field,value", [
    ("image", "final_fov_rg", 32.6868), ("image", "gaussian_sigma_rg", 2.16871),
    ("aperture", "radial_max_rg", 16.3434),
])
def test_reject_mixed_modes(profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match="mutually exclusive"):
        analysis_geometry(profile)


@pytest.mark.parametrize("section,field", [
    ("image", "final_fov_uas"), ("image", "gaussian_fwhm_uas"),
    ("aperture", "radial_max_fraction"),
])
@pytest.mark.parametrize("value", [True, False, "1", None, float("nan"), float("inf"), -float("inf")])
def test_new_numeric_fields_strictly_finite(profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError, match=rf"{section}\.{field}.*finite"):
        _validate(profile)


@pytest.mark.parametrize("section,field,value", [
    ("image", "final_fov_uas", 0), ("image", "final_fov_uas", -1),
    ("image", "final_fov_uas", 200), ("image", "gaussian_fwhm_uas", -1),
    ("aperture", "radial_max_fraction", 0), ("aperture", "radial_max_fraction", -0.1),
    ("aperture", "radial_max_fraction", 1.001), ("physical", "distance_mpc", 30),
    ("eht", "min_radius_rg", 20),
])
def test_angular_and_aperture_ranges(profile, section, field, value):
    profile[section][field] = value
    with pytest.raises(ValueError):
        _validate(profile)


@pytest.mark.parametrize("section", ["image", "aperture", "physical"])
def test_alternative_schemas_have_no_optional_or_unknown_fields(profile, section):
    for field in profile[section]:
        changed = deepcopy(profile)
        del changed[section][field]
        with pytest.raises(ValueError, match="Missing or invalid"):
            analysis_geometry(changed)
    profile[section]["typo"] = 1
    with pytest.raises(ValueError, match="Unknown"):
        analysis_geometry(profile)


def test_roundtrip_validation_does_not_add_legacy_keys(profile, tmp_path):
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(profile))
    assert read_config(path) == profile
    assert "gaussian_sigma_rg" not in profile["image"]
    assert "radial_max_rg" not in profile["aperture"]


def test_label_step_is_explicit_and_producer_acceptance_is_preserved(profile):
    # The analysis-only helper uses the chosen label unit; producer validation
    # still enforces its existing label_step_tg=0.1 acceptance contract.
    profile["physical"]["label_step_tg"] = 0.2
    with pytest.raises(ValueError, match="label_step_tg=0.1"):
        validate_production_config(profile)
    del profile["producer"]
    result = analysis_geometry(profile)
    _, seconds = independent_scales(profile["physical"])
    assert result["metadata"]["label_step_tg"] == 0.2
    assert result["metadata"]["label_step_hours"] == pytest.approx(0.2 * seconds / 3600)


def test_explicit_scale_helpers_and_legacy_defaults(profile):
    scale = PhysicalScale(**profile["physical"])
    theta, _ = independent_scales(profile["physical"])
    assert scale.angular_to_rg(128) == pytest.approx(128 / theta)
    fwhm = profile["image"]["gaussian_fwhm_uas"]
    assert scale.gaussian_sigma_rg(fwhm) == pytest.approx(fwhm / theta / math.sqrt(8 * math.log(2)))
    assert scale.gaussian_sigma_rg(0) == 0
    assert PhysicalScale() == scale  # Public legacy defaults are intentionally retained.
    # Analysis does not use these defaults for missing YAML physical values.
    del profile["physical"]["label_step_tg"]
    with pytest.raises(ValueError, match="Missing or invalid"):
        analysis_geometry(profile)


@pytest.mark.parametrize("field,value", [("final_fov_uas", 96), ("final_shape", 128)])
def test_analysis_target_must_match_declared_eht_grid(profile, field, value):
    # Intentionally retain the original input contract to exercise mismatch rejection.
    profile["image"][field] = value
    with pytest.raises(ValueError, match="target grid differs from declared EHT input sky grid"):
        analysis_geometry(profile)


def test_changed_analysis_grid_with_synchronized_eht_contract(profile):
    profile["image"].update(final_fov_uas=96, final_shape=128)
    profile["eht_input"].update(fov_uas=96, shape=128, pixel_size_uas=0.75)
    result = analysis_geometry(profile)
    assert result["image_preprocess"]["final_shape"] == 128
    assert result["metadata"]["final_fov_uas"] == pytest.approx(96)
    assert result["metadata"]["pixel_size_uas"] == pytest.approx(0.75)


def test_analysis_only_grid_without_eht_input_contract(profile):
    # No EHT input is claimed: alternate analysis-only geometry remains supported.
    del profile["eht_input"]
    profile["image"].update(final_fov_uas=96, final_shape=128)
    result = analysis_geometry(profile)
    assert result["metadata"]["pixel_size_uas"] == pytest.approx(0.75)


@pytest.mark.parametrize("field", ["coefficient", "pivot_alpha", "exponent", "reference_norm"])
@pytest.mark.parametrize("value", [True, False, "0.1", None, float("nan"),
                                   float("inf"), -float("inf")])
def test_read_config_rejects_nonfinite_or_coerced_cloud_values(profile, tmp_path, field, value):
    profile["cloud"][field] = value
    path = tmp_path / "invalid-cloud.yaml"
    path.write_text(yaml.safe_dump(profile))
    with pytest.raises(ValueError, match=rf"cloud\.{field}.*finite"):
        read_config(path)


@pytest.mark.parametrize("field", ["coefficient", "pivot_alpha", "reference_norm"])
@pytest.mark.parametrize("value", [0, -0.1])
def test_read_config_requires_positive_cloud_scales(profile, tmp_path, field, value):
    profile["cloud"][field] = value
    path = tmp_path / "invalid-cloud.yaml"
    path.write_text(yaml.safe_dump(profile))
    with pytest.raises(ValueError, match=rf"cloud\.{field}.*> 0"):
        read_config(path)


@pytest.mark.parametrize("exponent", [0, -1.25])
def test_read_config_accepts_finite_nonpositive_cloud_exponent(profile, tmp_path, exponent):
    profile["cloud"]["exponent"] = exponent
    path = tmp_path / "cloud.yaml"
    path.write_text(yaml.safe_dump(profile))
    assert read_config(path) == profile


@pytest.mark.parametrize("section", ["cloud", "eht_input"])
def test_optional_sections_are_complete_strict_contracts(profile, tmp_path, section):
    path = tmp_path / "invalid.yaml"
    for field in profile[section]:
        changed = deepcopy(profile)
        del changed[section][field]
        path.write_text(yaml.safe_dump(changed))
        with pytest.raises(ValueError, match=rf"Missing or invalid.*{section}"):
            read_config(path)
    profile[section]["typo"] = 1
    path.write_text(yaml.safe_dump(profile))
    with pytest.raises(ValueError, match=rf"Unknown.*{section}"):
        read_config(path)


@pytest.mark.parametrize("field", ["cloud", "producer"])
def test_read_config_rejects_cloud_producer_norm_mismatch(profile, tmp_path, field):
    if field == "cloud":
        profile["cloud"]["reference_norm"] = 0.2
    else:
        profile["producer"]["response"]["reference_norm"] = 0.2
    path = tmp_path / "mismatched-norm.yaml"
    path.write_text(yaml.safe_dump(profile))
    with pytest.raises(
        ValueError, match="cloud.reference_norm differs from generated response normalization",
    ):
        read_config(path)


def test_read_config_allows_omitted_optional_contracts(profile, tmp_path):
    del profile["cloud"], profile["eht_input"]
    path = tmp_path / "without-optional-contracts.yaml"
    path.write_text(yaml.safe_dump(profile))
    assert read_config(path) == profile
