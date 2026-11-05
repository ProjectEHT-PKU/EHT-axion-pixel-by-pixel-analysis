"""Strict production science contract; paths and execution policy stay in callers.

The analysis portion uses the existing validator unchanged. ``producer`` is
optional for validation of analysis-only profiles, but required by the reader.
"""

from pathlib import Path
import math

from ..pipeline.config import _validate as validate_analysis, read_config
from ..pipeline.layout import alpha_label

GRID = {
    "nxlone": [96, 48, 48],
    "xprobmin": [0.17395330712343798, 0.0, 0.0],
    "xprobmax": [7.824046010856292, 0.5, 1.0],
    "hslope": 0.0,
}
BUILD_CHOICES = {
    "model": "bhac", "metric": "mks", "integrator": "rk4", "grid": "SMR",
    "polarized": True, "sfc": True, "magnetic_polarity": "plus", "uniform_output": False,
}
FIELDS = {
    "build": set(BUILD_CHOICES) | {"sigma_cut", "thetae_min", "thetae_max"},
    "plasma": {"mass_unit_g", "r_low", "r_high", "inclination_degrees", "frequency_hz", "step_size"},
    "snapshots": {"start", "stop", "step"},
    "response": {"resolution", "phase_count", "reference_norm"},
    "background": {"resolution"},
    "radial_tables": {"low", "high"},
    "alpha_table": None,
    "grid": set(GRID),
}
ALPHA_TABLE = {"0.10": "low", "0.15": "low", "0.20": "low", "0.25": "low", "0.30": "low",
               "0.35": "high", "0.40": "high", "0.44": "high"}


def canonical_alpha(alpha):
    """Producer map keys omit the filesystem layout's ``alpha-`` prefix."""
    return alpha_label(alpha).removeprefix("alpha-")


def _keys(value, expected, name):
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    if set(value) != set(expected):
        raise ValueError(f"{name}: missing {sorted(set(expected) - set(value))}; "
                         f"unknown {sorted(set(value) - set(expected), key=str)}")


def _number(value, name, *, positive=True):
    if type(value) not in (int, float):
        raise ValueError(f"{name} must be a finite number (not bool or string)")
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    if not finite or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))


def _integer(value, name, minimum=1, maximum=2**31 - 1):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")


def validate_production_config(config):
    """Validate an in-memory profile without coercing, sorting or adding defaults."""
    if not isinstance(config, dict):
        raise ValueError("Production YAML configuration must contain a mapping")
    validate_analysis({key: value for key, value in config.items() if key != "producer"})
    if "producer" not in config:
        return
    p = config["producer"]
    _keys(p, FIELDS, "producer")
    for section, fields in FIELDS.items():
        if fields is not None:
            _keys(p[section], fields, f"producer.{section}")
    build = p["build"]
    for key, expected in BUILD_CHOICES.items():
        actual = build[key]
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(f"producer.build.{key} supports only {expected!r}")
    for key in ("sigma_cut", "thetae_min", "thetae_max"):
        _number(build[key], f"producer.build.{key}")
    if build["thetae_min"] >= build["thetae_max"]:
        raise ValueError("producer.build thetae_min must be less than thetae_max")
    for key, value in p["plasma"].items():
        _number(value, f"producer.plasma.{key}", positive=key != "inclination_degrees")
    if not 0 <= p["plasma"]["inclination_degrees"] <= 180:
        raise ValueError("producer.plasma.inclination_degrees must be in [0, 180]")
    for key in ("start", "stop", "step"):
        _integer(p["snapshots"][key], f"producer.snapshots.{key}", minimum=0 if key != "step" else 1,
                 maximum=999999)
    start, stop, step = (p["snapshots"][key] for key in ("start", "stop", "step"))
    if stop < start or (stop - start) % step:
        raise ValueError("producer.snapshots requires stop >= start on the inclusive step grid")
    if step != 10 or start % 10 or stop % 10:
        raise ValueError("Producer first scope requires snapshots.step=10 and start/stop multiples of 10")
    if config["physical"]["label_step_tg"] != 0.1:
        raise ValueError("Producer requires physical.label_step_tg=0.1")
    if config["inference"]["principal_n"] != 2:
        raise ValueError("Producer first scope requires inference.principal_n=2")
    for section in ("response", "background"):
        value = p[section]["resolution"]
        _integer(value, f"producer.{section}.resolution")
        if value % 10 or (value // 10)**2 > 2**31 - 1:
            raise ValueError(f"producer.{section}.resolution must be a core-compatible multiple of 10")
    _integer(p["response"]["phase_count"], "producer.response.phase_count")
    if p["response"]["phase_count"] != 16:
        raise ValueError("producer.response.phase_count must be 16 (core pi/8 convention)")
    _number(p["response"]["reference_norm"], "producer.response.reference_norm")
    # Research response pilots may explicitly select any finite positive norm.
    # The default profile remains .1; plan/job identities bind the actual norm.
    # This generalizes the producer contract, not the linear-response inference.
    # Explicit research grid generalization, not a core bug fix or inference change.
    # The unchanged production profile remains100; background stays400.
    if p["response"]["resolution"] not in (100, 200, 400) or p["background"]["resolution"] != 400:
        raise ValueError("Producer requires response.resolution in {100,200,400} "
                         "and background.resolution=400")
    if p["background"]["resolution"] != config["background"]["source_resolution"]:
        raise ValueError("producer.background.resolution must match analysis background.source_resolution")
    for table_id, spin in (("low", 0.9375), ("high", 0.99)):
        _keys(p["radial_tables"][table_id], {"spin"}, f"producer.radial_tables.{table_id}")
        value = p["radial_tables"][table_id]["spin"]
        _number(value, f"producer.radial_tables.{table_id}.spin")
        if value != spin:
            raise ValueError(f"producer.radial_tables.{table_id}.spin must be {spin}")
    expected = {canonical_alpha(alpha): ALPHA_TABLE.get(canonical_alpha(alpha))
                for alpha in config["data"]["alphas"]}
    _keys(p["alpha_table"], expected, "producer.alpha_table (canonical string keys)")
    if any(table_id is None or p["alpha_table"][key] != table_id for key, table_id in expected.items()):
        raise ValueError("producer.alpha_table requires supported low/high alpha assignments")
    for key, expected_value in GRID.items():
        actual = p["grid"][key]
        if isinstance(expected_value, list):
            if not isinstance(actual, list) or len(actual) != 3:
                raise ValueError(f"producer.grid.{key} must contain three numbers")
            for value in actual:
                if key == "nxlone":
                    _integer(value, f"producer.grid.{key}")
                else:
                    _number(value, f"producer.grid.{key}", positive=False)
        else:
            _number(actual, f"producer.grid.{key}", positive=False)
        if actual != expected_value:
            raise ValueError(f"producer.grid.{key} differs from the BHAC core template")


def read_production_config(path) -> dict:
    """Read duplicate-safe YAML with full analysis and required producer sections."""
    config = read_config(Path(path))
    validate_production_config(config)
    if "producer" not in config:
        raise ValueError("Production configuration requires producer section")
    return config
