"""Strict scientific YAML contract, without filesystem or runtime options.

Validation preserves supplied values: no defaults, coercion, sorting or new
physical cuts. Small synthetic grids obey the same contract as research.yaml.
"""

import math
from pathlib import Path

import yaml

from .layout import alpha_label


class ScientificConfigLoader(yaml.SafeLoader):
    """Safe YAML with duplicate (including merged) keys rejected."""


def _unique_mapping(loader, node, deep=False):
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str):
            raise ValueError("Scientific YAML configuration keys must be strings")
        if key in result:
            raise ValueError(f"Duplicate YAML configuration key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


ScientificConfigLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)

# Every listed field is required. Paths and execution switches belong outside YAML.
_FIELDS = {
    "physical": {"mass_solar", "distance_mpc", "label_step_tg"},
    "data": {"alphas"},
    "image": {
        "initial_fov_rg", "final_fov_rg", "final_shape", "rotation_degrees", "gaussian_sigma_rg",
    },
    "aperture": {"weighting", "window_degrees", "radial_max_rg", "angular_bins"},
    "eht": {
        "topset_size", "min_radius_rg", "intensity_fraction", "historical_mask",
        "inherit_2d_mask_in_1d", "angular_offset", "covariance_ridge",
    },
    "background": {"alpha", "source_resolution", "lag_frames", "method"},
    "inference": {
        "phase_count", "coupling_points", "coupling_min", "coupling_max", "prior",
        "principal_n", "exposure_hours", "day_gap_hours", "pair_separation_hours",
    },
}


def _keys(value, expected, name):
    if not isinstance(value, dict):
        raise ValueError(f"Missing or invalid scientific configuration section: {name}")
    missing, unknown = expected - value.keys(), value.keys() - expected
    if missing:
        raise ValueError(f"Missing or invalid scientific configuration fields in {name}: {sorted(missing)}")
    if unknown:
        raise ValueError(f"Unknown scientific configuration fields in {name}: {sorted(unknown)}")


def _number(config, section, field, *, minimum=None, maximum=None, open_min=False, open_max=False):
    value = config[section][field]
    name = f"{section}.{field}"
    # bool is an int subclass; strings, dates and nulls must not be coerced.
    if type(value) not in (int, float):
        raise ValueError(f"{name} must be a finite number (not bool or string)")
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    if not finite:
        raise ValueError(f"{name} must be a finite number")
    if minimum is not None and (value < minimum or (open_min and value == minimum)):
        raise ValueError(f"{name} must be {'>' if open_min else '>='} {minimum}")
    if maximum is not None and (value > maximum or (open_max and value == maximum)):
        raise ValueError(f"{name} must be {'<' if open_max else '<='} {maximum}")
    return value


def _integer(config, section, field, minimum=1, maximum=None):
    value = config[section][field]
    if type(value) is not int:
        raise ValueError(f"{section}.{field} must be an integer (not bool or float)")
    return _number(config, section, field, minimum=minimum, maximum=maximum)


def _choice(config, section, field, choices):
    value = config[section][field]
    if not isinstance(value, str) or value not in choices:
        raise ValueError(f"{section}.{field} must be one of {choices}")


def _alpha(value, name):
    if type(value) not in (int, float):
        raise ValueError(f"{name} must be a positive finite numeric alpha with two-decimal precision")
    try:
        return alpha_label(value)
    except (ValueError, ArithmeticError) as exc:
        raise ValueError(
            f"{name} requires a positive finite alpha with canonical two-decimal precision"
        ) from exc


def _validate(config):
    expected = set(_FIELDS) | ({"producer"} if "producer" in config else set())
    _keys(config, expected, "root")
    for section, fields in _FIELDS.items():
        if section == "inference" and isinstance(config[section], dict):
            if {"observation", "observed"} & config[section].keys():
                raise ValueError(
                    "Only zero EVPA day differences are supported; remove observation/observed from YAML"
                )
        _keys(config[section], fields, section)

    for field in _FIELDS["physical"]:
        _number(config, "physical", field, minimum=0, open_min=True)
    alphas = config["data"]["alphas"]
    if not isinstance(alphas, list) or not alphas:
        raise ValueError("data.alphas must be a nonempty list")
    labels = [_alpha(a, "data.alphas") for a in alphas]
    if len(set(labels)) != len(labels) or any(a >= b for a, b in zip(alphas, alphas[1:])):
        raise ValueError("data.alphas must be unique and strictly increasing in canonical identity")

    initial = _number(config, "image", "initial_fov_rg", minimum=0, open_min=True)
    final = _number(config, "image", "final_fov_rg", minimum=0, maximum=initial, open_min=True)
    _integer(config, "image", "final_shape", minimum=2)
    _number(config, "image", "rotation_degrees")
    _number(config, "image", "gaussian_sigma_rg", minimum=0)
    _choice(config, "aperture", "weighting", ("intensity", "flux"))
    _number(config, "aperture", "window_degrees", minimum=0, maximum=360, open_min=True)
    radius = _number(config, "aperture", "radial_max_rg", minimum=0, maximum=final / 2, open_min=True)
    bins = _integer(config, "aperture", "angular_bins")
    if bins != 360:
        raise ValueError("Research profile requires 360 angular bins (one-degree grid)")

    _integer(config, "eht", "topset_size")
    _number(config, "eht", "min_radius_rg", minimum=0, maximum=radius, open_max=True)
    for field in ("intensity_fraction", "covariance_ridge"):
        _number(config, "eht", field, minimum=0, maximum=1, open_max=True)
    for field in ("historical_mask", "inherit_2d_mask_in_1d"):
        if type(config["eht"][field]) is not bool:
            raise ValueError(f"eht.{field} must be a boolean")
    # build_observations uses a fixed ten-bin spacing, independent of window width.
    _integer(config, "eht", "angular_offset", minimum=0, maximum=9)

    _alpha(config["background"]["alpha"], "background.alpha")
    resolution = _integer(config, "background", "source_resolution")
    if resolution not in (100, 400):
        raise ValueError("background.source_resolution must be 100 or 400")
    _integer(config, "background", "lag_frames")
    _choice(config, "background", "method", ("mle", "legacy_hist"))
    _integer(config, "inference", "phase_count", minimum=2)
    _integer(config, "inference", "coupling_points", minimum=2)
    lower = _number(config, "inference", "coupling_min", minimum=0, open_min=True)
    _number(config, "inference", "coupling_max", minimum=lower, open_min=True)
    _choice(config, "inference", "prior", ("log-uniform", "uniform"))
    _integer(config, "inference", "principal_n")
    _number(config, "inference", "day_gap_hours", minimum=0, open_min=True)
    for field in ("exposure_hours", "pair_separation_hours"):
        _number(config, "inference", field, minimum=0)


def read_config(path):
    """Read safe, unambiguous YAML and reject invalid science before stage IO."""
    path = Path(path)
    if path.suffix.lower() not in {".yaml", ".yml"}:
        raise ValueError("Scientific configuration must be YAML (.yaml or .yml)")
    try:
        config = yaml.load(path.read_text(encoding="utf-8"), Loader=ScientificConfigLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid scientific YAML configuration: {path}: {exc}") from exc
    if not isinstance(config, dict):
        raise ValueError("Scientific YAML configuration must contain a mapping")
    _validate(config)
    if "producer" in config:
        from ..producer.config import validate_production_config

        validate_production_config(config)
    return config
