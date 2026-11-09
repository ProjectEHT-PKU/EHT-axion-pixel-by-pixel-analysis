"""experiment.yaml: structure, validation, derived geometry and process views."""

from copy import deepcopy
import math

import pytest
import yaml

from eht_axion import experiment
from eht_axion.experiment import (
    analysis_geometry, eht_imaging_parameters, physical_scale, ray_tracing_config, read_experiment,
    validate, validate_ray_tracing,
)

LOCATIONS = dict(input_root="/data/uvfits", input_receipt="/data/receipt.json",
                 feed_rotation_file="/data/feed.json", official_script="/data/pipeline.py")


@pytest.fixture
def config():
    return read_experiment()


def write(tmp_path, value):
    path = tmp_path / "experiment.yaml"
    path.write_text(yaml.safe_dump(value))
    return path


def test_repository_file_has_shared_two_inputs_and_analysis(config):
    assert set(config) == {"shared", "raptor", "eht_imaging", "analysis"}
    assert config["shared"]["image"] == {"shape": 64, "fov_uas": 120.0, "blur_fwhm_uas": 19.0}
    assert "centering" not in config["analysis"]["features"]  # Paper VII: zero shift only


def independent_scales(shared):
    gm = 1.3271244e20 * shared["mass_solar"]
    theta = gm / 299792458.0**2 / (shared["distance_mpc"] * 1e6 * 3.085677581491367e16)
    return theta * (180 / math.pi) * 3600 * 1e6, gm / 299792458.0**3


def test_analysis_geometry_from_independent_formula(config):
    before = deepcopy(config)
    result = analysis_geometry(config)
    theta, seconds = independent_scales(config["shared"])
    image, aperture, metadata = (result[k] for k in ("image_preprocess", "aperturefit", "metadata"))
    assert image == pytest.approx(dict(
        sigma=19 / theta / math.sqrt(8 * math.log(2)), init_fov=46.2261, final_fov=120 / theta,
        final_shape=64, rotation_angle=108 * math.pi / 180))
    assert aperture == dict(fov=image["final_fov"], radial_max=image["final_fov"] / 2,
                            window_degrees=10, angular_bins=360)
    assert metadata["tg_seconds"] == pytest.approx(seconds)
    assert metadata["final_fov_uas"] == pytest.approx(120)
    assert metadata["radial_max_uas"] == 60.0
    assert metadata["pixel_size_uas"] == pytest.approx(1.875)
    assert config == before


def test_mass_distance_propagate(config):
    baseline = analysis_geometry(config)
    config["shared"]["mass_solar"] *= 1.25
    result = analysis_geometry(config)
    assert result["metadata"]["tg_hours"] == pytest.approx(baseline["metadata"]["tg_hours"] * 1.25)
    assert result["image_preprocess"]["final_fov"] == pytest.approx(
        baseline["image_preprocess"]["final_fov"] / 1.25)
    assert physical_scale(config).mass_solar == config["shared"]["mass_solar"]


def test_ray_tracing_subset_excludes_sky_grid_imaging_and_analysis(config):
    subset = ray_tracing_config(config)
    assert set(subset) == {"shared", "raptor"}
    assert set(subset["shared"]) == {"mass_solar", "distance_mpc", "label_step_tg", "alphas", "principal_n"}
    assert subset["raptor"] == config["raptor"]
    assert ray_tracing_config(subset) == subset
    validate_ray_tracing(subset)
    with pytest.raises(ValueError):
        validate_ray_tracing(dict(subset, analysis=config["analysis"]))


def test_eht_imaging_parameters_flat_set_with_shared_grid(config):
    parameters = eht_imaging_parameters(config, LOCATIONS)
    assert list(parameters)[:5] == ["input_root", "input_receipt", "feed_rotation_file", "official_script",
                                    "inputs"]
    assert (parameters["npix"], parameters["fov_uas"], parameters["display_blur_uas"]) == (64, 120.0, 19.0)
    assert parameters["uncertainty_draws"] == 1000
    imaging = {k: v for k, v in config["eht_imaging"].items() if k not in ("refine", "monte_carlo")}
    for key, value in imaging.items():
        assert parameters[key] == value
    # The solver settings of eht_refine/eht_sample are not part of the imaging parameter set.
    assert len(parameters) == len(imaging) + 3 + len(LOCATIONS)
    assert "refine" not in parameters and "monte_carlo" not in parameters
    with pytest.raises(ValueError, match="locations"):
        eht_imaging_parameters(config, {})


def test_duplicate_keys_rejected(tmp_path):
    path = tmp_path / "experiment.yaml"
    path.write_text(experiment.DEFAULT_PATH.read_text() + "\nshared: {}\n")
    with pytest.raises(ValueError, match="Duplicate YAML key"):
        read_experiment(path)


@pytest.mark.parametrize("section,key", [
    ("shared", "alphas"), ("raptor", "plasma"), ("eht_imaging", "fsys"), ("analysis", "features"),
])
def test_missing_field_rejected(config, tmp_path, section, key):
    del config[section][key]
    with pytest.raises(ValueError, match="missing"):
        read_experiment(write(tmp_path, config))


def test_unknown_section_and_field_rejected(config):
    with pytest.raises(ValueError, match="exactly"):
        validate(dict(config, configs={}))
    config["analysis"]["inference"]["observation"] = "measured"
    with pytest.raises(ValueError, match="unknown"):
        validate(config)


@pytest.mark.parametrize("path,value,match", [
    (("shared", "mass_solar"), "6.5e9", "finite number"),
    (("shared", "alphas"), [0.2, 0.1], "increasing"),
    (("shared", "alphas"), [0.123], "two-decimal"),
    (("shared", "principal_n"), 2.0, "integer"),
    (("raptor", "response", "resolution"), 105, "multiple of 10"),
    (("raptor", "snapshots", "stop"), 1605, "step grid"),
    (("raptor", "alpha_table"), {"0.10": "low"}, "alpha_table"),
    (("raptor", "plasma", "inclination_degrees"), 181.0, "<= 180"),
    (("eht_imaging", "restart_rounds"), [1], "restart_rounds"),
    (("eht_imaging", "fixed_sites"), ["AZ"], "overlap"),
    (("eht_imaging", "norm_reg"), 1, "true or false"),
    (("analysis", "features", "polarization_fraction"), 1.0, "< 1"),
    (("analysis", "features", "polarization_fraction"), -0.1, ">= 0"),
    (("analysis", "features", "sector_offset_degrees"), 10, "<= 9"),
    (("analysis", "features", "min_radius_rg"), 100.0, "inside the aperture"),
    (("analysis", "inference", "coupling_max"), 0.001, "coupling_max"),
])
def test_invalid_values_rejected(config, path, value, match):
    target = config
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError, match=match):
        validate(config)


def test_feature_thresholds_always_both_apply(config):
    assert config["analysis"]["features"] == dict(
        intensity_fraction=0.1, polarization_fraction=0.2, min_radius_rg=2.0, pixel_stride=5,
        sector_spacing_degrees=10, sector_offset_degrees=4, resultant_min=0.7)
    config["analysis"]["features"]["polarization_fraction"] = 0  # brightness-only selection
    validate(config)
    config["analysis"]["features"]["mask_mode"] = "intensity"
    with pytest.raises(ValueError, match="unknown"):
        validate(config)


def test_sky_grid_must_fit_in_camera(config):
    config["shared"]["image"]["fov_uas"] = 1000.0
    with pytest.raises(ValueError, match="camera"):
        validate(config)


def test_library_signatures_have_no_numeric_default_values():
    """Scientific values come from experiment.yaml, never from library defaults."""
    import ast
    from pathlib import Path

    allowed = {
        ("raptor/jobs.py", "run_jobs", "threads"), ("raptor/jobs.py", "run_jobs", "workers"),  # execution
        ("analysis/templates.py", "fit_phase_dft", "phase_axis"),  # array layout
        ("raptor/plan.py", "append_task", "phase"),  # baseline tasks have phase 0 by definition
        ("eht/calibration.py", "tlist", "t_gather"),  # mirrors the ehtim method signature
    }
    found = []
    for path in sorted((Path(experiment.__file__).parent).rglob("*.py")):
        relative = str(path.relative_to(Path(experiment.__file__).parent))
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.FunctionDef) or node.name.startswith("_"):
                continue
            args = node.args
            pairs = list(zip(args.args[len(args.args) - len(args.defaults):], args.defaults))
            pairs += [(a, d) for a, d in zip(args.kwonlyargs, args.kw_defaults) if d is not None]
            for arg, default in pairs:
                value = default.value if isinstance(default, ast.Constant) else None
                if (isinstance(value, (int, float)) and not isinstance(value, bool)
                        and (relative, node.name, arg.arg) not in allowed):
                    found.append(f"{relative}:{node.name}({arg.arg}={value})")
    assert not found
