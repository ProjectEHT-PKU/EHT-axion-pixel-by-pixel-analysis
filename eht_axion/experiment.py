"""experiment.yaml: the single parameter file, its validation and derived settings.

Sections (see experiment.yaml): ``shared`` (used by more than one part), the two
independent inputs ``raptor`` (simulation) and ``eht_imaging`` (observation), and
``analysis``, which combines them. Each script reads ``shared`` plus its own
section; the analysis also reads the raptor/eht_imaging values it needs. Validation
checks structure, types and ranges only; it fills no defaults and never compares
a value with a copy of itself. Derived settings (unit conversions, the subset of
values that determines RAPTOR rays, the eht-imaging parameter set) are computed
here, so every script reads the same numbers from the same place.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import yaml

from .physics import PhysicalScale
from .raptor.layout import alpha_label

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATH = ROOT / "experiment.yaml"

PHYSICS_KEYS = ("mass_solar", "distance_mpc", "label_step_tg", "alphas", "principal_n")
SCHEMA = {
    "shared": set(PHYSICS_KEYS) | {"image"},
    "shared.image": {"shape", "fov_uas", "blur_fwhm_uas"},
    "raptor": {"source", "camera_fov_rg", "fluid_spin", "build", "plasma", "snapshots", "response",
               "background", "radial_tables", "alpha_table", "grid"},
    "raptor.source": {"repository", "commit"},
    "raptor.build": {"model", "metric", "integrator", "grid", "polarized", "sfc", "magnetic_polarity",
                     "sigma_cut", "thetae_min", "thetae_max", "uniform_output"},
    "raptor.plasma": {"mass_unit_g", "r_low", "r_high", "inclination_degrees", "frequency_hz",
                      "step_size"},
    "raptor.snapshots": {"start", "stop", "step"},
    "raptor.response": {"resolution", "reference_norm"},
    "raptor.background": {"resolution"},
    "raptor.radial_tables": {"low", "high"},
    "raptor.grid": {"nxlone", "xprobmin", "xprobmax", "hslope"},
    "eht_imaging": {
        "inputs", "gain_tol", "fsys", "solve_sites", "fixed_sites", "p_data", "p_reg", "norm_reg",
        "initial_m", "initial_evpa_deg", "rounds", "restart_rounds", "restart_blur_uas",
        "restart_domain_epsilon", "maxit", "stop", "shift_objective_constant", "uv_min",
        "uncertainty_draws", "seed", "covariance_ddof", "residual_d_sigma",
        "refine", "monte_carlo",
    },
    "eht_imaging.refine": {"outer_rounds", "gap_atol", "maxiter", "check_every", "time_limit",
                           "d_optimizer"},
    "eht_imaging.refine.d_optimizer": {"ftol", "gtol", "maxiter", "maxcor", "maxls", "maxfun"},
    "eht_imaging.monte_carlo": {"initial_p_fraction", "gap_atol", "maxiter", "check_every", "time_limit",
                                "retry_maxiter", "retry_time_limit"},
    "analysis": {"rotation_degrees", "aperture", "features", "cloud", "background", "inference"},
    "analysis.aperture": {"window_degrees", "radial_max_fraction"},
    "analysis.features": {"intensity_fraction", "polarization_fraction",
                          "min_radius_rg", "pixel_stride", "sector_spacing_degrees",
                          "sector_offset_degrees", "resultant_min"},
    "analysis.cloud": {"coefficient", "pivot_alpha", "exponent"},
    "analysis.background": {"lag_frames"},
    "analysis.inference": {"phase_count", "coupling_points", "coupling_min", "coupling_max", "prior",
                           "credibility", "exposure_hours", "day_gap_hours", "pair_separation_hours"},
}
EHT_DAYS = (95, 96, 100, 101)
# The 1D observable uses one-degree aperture bins (an angle grid, not a parameter);
# aperture EVPAs are always intensity weighted (EHT VII).
ANGULAR_BINS = 360
# Solver settings of eht_refine.py / eht_sample.py: recorded in their own plans,
# not in the flat parameter set of the imaging stages.
EHT_SOLVER_SECTIONS = ("refine", "monte_carlo")
# Machine-local eht-imaging locations; recorded in the imaging plan, never in YAML.
EHT_LOCATIONS = ("input_root", "input_receipt", "feed_rotation_file", "official_script")


class _UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML with duplicate (including merged) keys rejected."""


def _unique_mapping(loader, node, deep=False):
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f"Duplicate YAML key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


# --------------------------------------------------------------------------- checks
def _section(config, name, expected=None):
    value = config
    for part in name.split("."):
        if not isinstance(value, dict) or part not in value:
            raise ValueError(f"Missing section: {name}")
        value = value[part]
    expected = SCHEMA[name] if expected is None else expected
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    missing, unknown = expected - value.keys(), value.keys() - expected
    if missing or unknown:
        raise ValueError(f"{name}: missing {sorted(missing)}, unknown {sorted(map(str, unknown))}")
    return value


def _number(value, name, *, low=None, high=None, open_low=False, open_high=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number (not bool or string)")
    if low is not None and (value < low or (open_low and value == low)):
        raise ValueError(f"{name} must be {'>' if open_low else '>='} {low}")
    if high is not None and (value > high or (open_high and value == high)):
        raise ValueError(f"{name} must be {'<' if open_high else '<='} {high}")
    return value


def _positive(value, name):
    return _number(value, name, low=0, open_low=True)


def _integer(value, name, low=1, high=None):
    if type(value) is not int:
        raise ValueError(f"{name} must be an integer")
    return _number(value, name, low=low, high=high)


def _choice(value, name, choices):
    if not isinstance(value, str) or value not in choices:
        raise ValueError(f"{name} must be one of {choices}")
    return value


def _boolean(value, name):
    if type(value) is not bool:
        raise ValueError(f"{name} must be true or false")
    return value


def _validate_shared(config):
    shared = _section(config, "shared")
    for key in ("mass_solar", "distance_mpc", "label_step_tg"):
        _positive(shared[key], f"shared.{key}")
    alphas = shared["alphas"]
    if not isinstance(alphas, list) or not alphas:
        raise ValueError("shared.alphas must be a nonempty list")
    for alpha in alphas:
        if type(alpha) not in (int, float):
            raise ValueError("shared.alphas must contain numbers")
        alpha_label(alpha)  # positive, two-decimal identity
    if any(a >= b for a, b in zip(alphas, alphas[1:])):
        raise ValueError("shared.alphas must be strictly increasing")
    _integer(shared["principal_n"], "shared.principal_n")
    image = _section(config, "shared.image")
    _integer(image["shape"], "shared.image.shape", low=2)
    _positive(image["fov_uas"], "shared.image.fov_uas")
    _number(image["blur_fwhm_uas"], "shared.image.blur_fwhm_uas", low=0)


def normalize_repository(url):
    """'host/owner/name' of a git URL: no scheme, user, port or trailing .git; lower case.

    https://github.com/A/b.git, git@github.com:A/b and ssh://git@github.com/A/b all give
    'github.com/a/b'.
    """
    text = str(url).strip()
    match = re.fullmatch(r"[a-z][a-z0-9+.-]*://(?:[^@/]+@)?([^/:]+)(?::\d+)?/(.*)", text, re.IGNORECASE)
    if match is None:
        match = re.fullmatch(r"(?:[^@/]+@)?([^/:]+):(.*)", text)  # scp-like user@host:path
    host, path = (match[1], match[2]) if match else ("", text)
    path = path.strip("/")
    if path.endswith(".git"):
        path = path[:-4]
    return f"{host}/{path}".strip("/").lower()


def _validate_raptor_source(config):
    source = _section(config, "raptor.source")
    repository, commit = source["repository"], source["commit"]
    if (not isinstance(repository, str) or normalize_repository(repository) != repository
            or not re.fullmatch(r"[a-z0-9.-]+/[a-z0-9_.-]+/[a-z0-9_.-]+", repository)):
        raise ValueError("raptor.source.repository must be host/owner/name in lower case, "
                         "for example github.com/owner/raptor")
    if not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("raptor.source.commit must be a full 40-character git commit")


def _validate_raptor(config, *, with_source=True):
    """The raptor section; the task-plan subset (ray_tracing_config) carries no source."""
    raptor = _section(config, "raptor", None if with_source else SCHEMA["raptor"] - {"source"})
    if with_source:
        _validate_raptor_source(config)
    _positive(raptor["camera_fov_rg"], "raptor.camera_fov_rg")
    _number(raptor["fluid_spin"], "raptor.fluid_spin", low=0, high=1)
    build = _section(config, "raptor.build")
    for key in ("model", "metric", "integrator", "grid", "magnetic_polarity"):
        if not isinstance(build[key], str):
            raise ValueError(f"raptor.build.{key} must be a string")
    for key in ("polarized", "sfc", "uniform_output"):
        _boolean(build[key], f"raptor.build.{key}")
    for key in ("sigma_cut", "thetae_min", "thetae_max"):
        _positive(build[key], f"raptor.build.{key}")
    if build["thetae_min"] >= build["thetae_max"]:
        raise ValueError("raptor.build.thetae_min must be less than thetae_max")
    plasma = _section(config, "raptor.plasma")
    for key, value in plasma.items():
        if key == "inclination_degrees":
            _number(value, "raptor.plasma.inclination_degrees", low=0, high=180)
        else:
            _positive(value, f"raptor.plasma.{key}")
    snapshots = _section(config, "raptor.snapshots")
    start = _integer(snapshots["start"], "raptor.snapshots.start", low=0, high=999999)
    stop = _integer(snapshots["stop"], "raptor.snapshots.stop", low=0, high=999999)
    step = _integer(snapshots["step"], "raptor.snapshots.step")
    if stop < start or (stop - start) % step:
        raise ValueError("raptor.snapshots requires stop >= start on the inclusive step grid")
    response = _section(config, "raptor.response")
    background = _section(config, "raptor.background")
    for name, section in (("response", response), ("background", background)):
        value = _integer(section["resolution"], f"raptor.{name}.resolution")
        if value % 10:  # the core writes 10x10-pixel blocks
            raise ValueError(f"raptor.{name}.resolution must be a multiple of 10")
    _positive(response["reference_norm"], "raptor.response.reference_norm")
    tables = _section(config, "raptor.radial_tables")
    for table_id, entry in tables.items():
        if not isinstance(entry, dict) or set(entry) != {"spin"}:
            raise ValueError(f"raptor.radial_tables.{table_id} must contain exactly spin")
        _number(entry["spin"], f"raptor.radial_tables.{table_id}.spin", low=0, high=1)
    expected = {alpha_label(a).removeprefix("alpha-") for a in config["shared"]["alphas"]}
    table = raptor["alpha_table"]
    if not isinstance(table, dict) or set(table) != expected:
        raise ValueError(f"raptor.alpha_table needs exactly the keys {sorted(expected)}")
    if any(value not in tables for value in table.values()):
        raise ValueError("raptor.alpha_table values must name raptor.radial_tables entries")
    grid = _section(config, "raptor.grid")
    for key in ("nxlone", "xprobmin", "xprobmax"):
        if not isinstance(grid[key], list) or len(grid[key]) != 3:
            raise ValueError(f"raptor.grid.{key} must contain three numbers")
        for value in grid[key]:
            (_integer if key == "nxlone" else _number)(value, f"raptor.grid.{key}")
    _number(grid["hslope"], "raptor.grid.hslope")


def _validate_eht_imaging(config):
    c = _section(config, "eht_imaging")
    if not isinstance(c["inputs"], dict) or set(c["inputs"]) != set(EHT_DAYS):
        raise ValueError(f"eht_imaging.inputs requires exactly the days {EHT_DAYS}")
    for day, entry in c["inputs"].items():
        if (not isinstance(entry, list) or len(entry) != 2 or not all(isinstance(v, str) for v in entry)
                or Path(entry[0]).name != entry[0] or len(entry[1]) != 64):
            raise ValueError(f"eht_imaging.inputs.{day} must be [file name, sha256]")
    for key in ("rounds", "maxit", "uncertainty_draws", "covariance_ddof"):
        _integer(c[key], f"eht_imaging.{key}")
    _integer(c["seed"], "eht_imaging.seed", low=0)
    gain_tol = c["gain_tol"]
    if not isinstance(gain_tol, list) or len(gain_tol) != 2:
        raise ValueError("eht_imaging.gain_tol must be [lower, upper] fractional gain tolerances")
    for value in gain_tol:
        _positive(value, "eht_imaging.gain_tol")
    for key in ("fsys", "restart_blur_uas", "stop", "uv_min", "residual_d_sigma"):
        _positive(c[key], f"eht_imaging.{key}")
    _number(c["initial_m"], "eht_imaging.initial_m", low=0, high=1, open_low=True, open_high=True)
    _number(c["restart_domain_epsilon"], "eht_imaging.restart_domain_epsilon",
            low=0, high=1, open_low=True, open_high=True)
    _number(c["initial_evpa_deg"], "eht_imaging.initial_evpa_deg")
    for key in ("norm_reg", "shift_objective_constant"):
        _boolean(c[key], f"eht_imaging.{key}")
    for key in ("solve_sites", "fixed_sites"):
        if not isinstance(c[key], list) or not c[key] or not all(isinstance(s, str) for s in c[key]):
            raise ValueError(f"eht_imaging.{key} must be a list of station codes")
    if set(c["solve_sites"]) & set(c["fixed_sites"]):
        raise ValueError("eht_imaging solve_sites and fixed_sites overlap")
    for key, names in (("p_data", {"pvis", "m"}), ("p_reg", {"hw", "ptv"})):
        if not isinstance(c[key], dict) or set(c[key]) != names:
            raise ValueError(f"eht_imaging.{key} requires exactly {sorted(names)}")
        for name, value in c[key].items():
            _positive(value, f"eht_imaging.{key}.{name}")
    restarts = c["restart_rounds"]
    if (not isinstance(restarts, list) or len(set(restarts)) != len(restarts)
            or any(type(r) is not int or not 2 <= r <= c["rounds"] for r in restarts)):
        raise ValueError("eht_imaging.restart_rounds must be distinct rounds in [2, rounds]")
    refine = _section(config, "eht_imaging.refine")
    for key in ("outer_rounds", "maxiter", "check_every"):
        _integer(refine[key], f"eht_imaging.refine.{key}")
    for key in ("gap_atol", "time_limit"):
        _positive(refine[key], f"eht_imaging.refine.{key}")
    d_optimizer = _section(config, "eht_imaging.refine.d_optimizer")
    for key in ("maxiter", "maxcor", "maxls", "maxfun"):
        _integer(d_optimizer[key], f"eht_imaging.refine.d_optimizer.{key}")
    for key in ("ftol", "gtol"):
        _positive(d_optimizer[key], f"eht_imaging.refine.d_optimizer.{key}")
    mc = _section(config, "eht_imaging.monte_carlo")
    for key in ("maxiter", "check_every", "retry_maxiter"):
        _integer(mc[key], f"eht_imaging.monte_carlo.{key}")
    for key in ("gap_atol", "time_limit", "retry_time_limit"):
        _positive(mc[key], f"eht_imaging.monte_carlo.{key}")
    _number(mc["initial_p_fraction"], "eht_imaging.monte_carlo.initial_p_fraction", low=0, high=1,
            open_low=True, open_high=True)


def _validate_analysis(config):
    analysis = _section(config, "analysis")
    _number(analysis["rotation_degrees"], "analysis.rotation_degrees")
    aperture = _section(config, "analysis.aperture")
    _number(aperture["window_degrees"], "analysis.aperture.window_degrees", low=0, high=360, open_low=True)
    _number(aperture["radial_max_fraction"], "analysis.aperture.radial_max_fraction",
            low=0, high=1, open_low=True)
    features = _section(config, "analysis.features")
    for key in ("intensity_fraction", "polarization_fraction"):
        _number(features[key], f"analysis.features.{key}", low=0, high=1, open_high=True)
    _number(features["min_radius_rg"], "analysis.features.min_radius_rg", low=0)
    _integer(features["pixel_stride"], "analysis.features.pixel_stride")
    spacing = _integer(features["sector_spacing_degrees"], "analysis.features.sector_spacing_degrees",
                       high=360)
    _integer(features["sector_offset_degrees"], "analysis.features.sector_offset_degrees",
             low=0, high=spacing - 1)
    _number(features["resultant_min"], "analysis.features.resultant_min", low=0, high=1, open_low=True)
    cloud = _section(config, "analysis.cloud")
    for key in ("coefficient", "pivot_alpha"):
        _positive(cloud[key], f"analysis.cloud.{key}")
    _number(cloud["exponent"], "analysis.cloud.exponent")
    background = _section(config, "analysis.background")
    _integer(background["lag_frames"], "analysis.background.lag_frames")
    inference = _section(config, "analysis.inference")
    _integer(inference["phase_count"], "analysis.inference.phase_count", low=2)
    _integer(inference["coupling_points"], "analysis.inference.coupling_points", low=2)
    lower = _positive(inference["coupling_min"], "analysis.inference.coupling_min")
    _number(inference["coupling_max"], "analysis.inference.coupling_max", low=lower, open_low=True)
    _choice(inference["prior"], "analysis.inference.prior", ("log-uniform", "uniform"))
    _number(inference["credibility"], "analysis.inference.credibility", low=0, high=1,
            open_low=True, open_high=True)
    _positive(inference["day_gap_hours"], "analysis.inference.day_gap_hours")
    for key in ("exposure_hours", "pair_separation_hours"):
        _number(inference[key], f"analysis.inference.{key}", low=0)
    geometry = analysis_geometry(config)
    if not geometry["image_preprocess"]["final_fov"] <= config["raptor"]["camera_fov_rg"]:
        raise ValueError("shared.image.fov_uas exceeds the RAPTOR camera field")
    if not features["min_radius_rg"] < geometry["aperturefit"]["radial_max"]:
        raise ValueError("analysis.features.min_radius_rg must lie inside the aperture")


def validate(config):
    """Validate a complete experiment dict (all four sections) in place-free fashion."""
    if not isinstance(config, dict) or set(config) != {"shared", "raptor", "eht_imaging", "analysis"}:
        raise ValueError("experiment.yaml requires exactly: shared, raptor, eht_imaging, analysis")
    _validate_shared(config)
    _validate_raptor(config)
    _validate_eht_imaging(config)
    _validate_analysis(config)
    return config


def validate_ray_tracing(config):
    """Validate the ray-tracing subset stored in task plans."""
    if not isinstance(config, dict) or set(config) != {"shared", "raptor"}:
        raise ValueError("Ray-tracing configuration requires exactly shared and raptor")
    if set(config["shared"]) != set(PHYSICS_KEYS):
        raise ValueError(f"Ray-tracing shared values must be exactly {PHYSICS_KEYS}")
    _validate_shared({"shared": dict(config["shared"], image={"shape": 2, "fov_uas": 1.,
                                                              "blur_fwhm_uas": 0.})})
    _validate_raptor(config, with_source=False)
    return config


def read_experiment(path=DEFAULT_PATH):
    """Read duplicate-key-safe YAML and validate every section before any IO."""
    path = Path(path)
    try:
        config = yaml.load(path.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML: {path}: {exc}") from exc
    return validate(config)


# --------------------------------------------------------------------- derived views
def physical_scale(config):
    shared = config["shared"]
    return PhysicalScale(mass_solar=shared["mass_solar"], distance_mpc=shared["distance_mpc"],
                         label_step_tg=shared["label_step_tg"])


def ray_tracing_config(config):
    """Physics scalars + raptor section without raptor.source: the ray-tracing settings.

    This subset is what task plans and generated-data identities record; the sky grid,
    eht-imaging and analysis settings are excluded. The RAPTOR source (repository and
    commit) is recorded by the plan as raptor_source and checked against every build.
    """
    if not isinstance(config, dict) or not {"shared", "raptor"} <= config.keys():
        raise ValueError("Ray-tracing configuration requires shared and raptor sections")
    return {"shared": {key: config["shared"][key] for key in PHYSICS_KEYS},
            "raptor": {k: v for k, v in config["raptor"].items() if k != "source"}}


def analysis_geometry(config) -> dict:
    """Preprocessing/aperture keyword arguments in gravitational radii.

    The RAPTOR images are blurred with shared.image.blur_fwhm_uas, rotated by
    analysis.rotation_degrees and resampled onto the shared sky grid.
    """
    scale = physical_scale(config)
    image, aperture = config["shared"]["image"], config["analysis"]["aperture"]
    final_fov = scale.angular_to_rg(image["fov_uas"])
    sigma = scale.gaussian_sigma_rg(image["blur_fwhm_uas"])
    radial_max = aperture["radial_max_fraction"] * (final_fov / 2)
    return {
        "image_preprocess": {
            "sigma": sigma,
            "init_fov": config["raptor"]["camera_fov_rg"],
            "final_fov": final_fov,
            "final_shape": image["shape"],
            "rotation_angle": np.deg2rad(config["analysis"]["rotation_degrees"]),
        },
        "aperturefit": {
            "fov": final_fov,
            "radial_max": radial_max,
            "window_degrees": aperture["window_degrees"],
            "angular_bins": ANGULAR_BINS,
        },
        "metadata": dict(
            scale.metadata(),
            initial_fov_uas=config["raptor"]["camera_fov_rg"] * scale.angular_rg_uas,
            final_fov_uas=final_fov * scale.angular_rg_uas,
            gaussian_fwhm_uas=sigma * np.sqrt(8 * np.log(2)) * scale.angular_rg_uas,
            pixel_size_uas=final_fov * scale.angular_rg_uas / image["shape"],
            radial_max_rg=radial_max,
            radial_max_uas=radial_max * scale.angular_rg_uas,
        ),
    }


def eht_imaging_parameters(config, locations):
    """Flat parameter set used by the eht-imaging scripts and recorded in their plans.

    ``locations`` holds the machine-local paths named in EHT_LOCATIONS.
    """
    if set(locations) != set(EHT_LOCATIONS):
        raise ValueError(f"eht-imaging locations require exactly {EHT_LOCATIONS}")
    c, image = config["eht_imaging"], config["shared"]["image"]
    parameters = {key: str(locations[key]) for key in EHT_LOCATIONS}
    parameters["inputs"] = c["inputs"]
    parameters.update(npix=image["shape"], fov_uas=image["fov_uas"])
    for key in SCHEMA["eht_imaging"] - {"inputs", *EHT_SOLVER_SECTIONS}:
        parameters[key] = c[key]
    parameters["blur_fwhm_uas"] = image["blur_fwhm_uas"]
    return {key: parameters[key] for key in _EHT_PARAMETER_ORDER}


_EHT_PARAMETER_ORDER = (
    "input_root", "input_receipt", "feed_rotation_file", "official_script", "inputs", "npix", "fov_uas",
    "gain_tol", "fsys", "solve_sites", "fixed_sites", "p_data", "p_reg", "norm_reg", "initial_m",
    "initial_evpa_deg", "rounds", "restart_rounds", "restart_blur_uas", "restart_domain_epsilon",
    "maxit", "stop", "shift_objective_constant", "uv_min", "blur_fwhm_uas", "uncertainty_draws",
    "seed", "covariance_ddof", "residual_d_sigma",
)
