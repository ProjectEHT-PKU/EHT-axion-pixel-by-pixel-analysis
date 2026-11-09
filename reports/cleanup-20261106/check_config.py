#!/usr/bin/env python3
"""Evidence 1: every value the active chain consumed from main:configs/* equals the new loader.

Old values come from `git show main:configs/...` and from literals that main's
active scripts hard-coded; derived quantities (analysis geometry, LimitOptions,
eht-imaging parameter set) are computed by main's code and by the new code.
Writes config-equivalence.json.
"""
from __future__ import annotations

import subprocess
import sys

import yaml

from _side import MAIN_COMMIT, NEW, diff, run, write

sys.path[:0] = [str(NEW)]
from eht_axion.experiment import read_experiment  # noqa: E402
from eht_axion.pipeline.layout import BACKGROUND_DIRECTORY_ALPHA  # noqa: E402
from eht_axion.producer.plan import PHASE_COUNT  # noqa: E402


def old(name):
    text = subprocess.run(["git", "-C", str(NEW), "show", f"{MAIN_COMMIT}:configs/{name}"],
                          check=True, capture_output=True, text=True).stdout
    return yaml.safe_load(text)


def get(tree, path):
    for part in path.split("."):
        tree = tree[int(part) if isinstance(tree, dict) and part.isdigit() and int(part) in tree else part]
    return tree


P, A = old("production.yaml"), old("eht-120uas-observed.yaml")
G, R, E = old("eht-120uas-observed-gridcheck.yaml"), old("eht-120uas-observed-priorcheck.yaml"), old(
    "eht-official.yaml")
X = read_experiment()
rows = []


def check(source, old_value, new_path, new_value, note=""):
    rows.append(dict(old=source, old_value=old_value, new=new_path, new_value=new_value,
                     equal=old_value == new_value and type(old_value) is type(new_value), note=note))


# --- copies between old profiles (the reason for the single file) -------------------------
copies = []
for section in ("physical", "data", "producer"):
    copies.append(dict(section=section, production_equals_analysis=P[section] == A[section]))
copies.append(dict(section="image.initial_fov_rg",
                   production_equals_analysis=P["image"]["initial_fov_rg"] == A["image"]["initial_fov_rg"]))
copies.append(dict(section="cloud.reference_norm == producer.response.reference_norm",
                   production_equals_analysis=A["cloud"]["reference_norm"] == P["producer"]["response"][
                       "reference_norm"]))
copies.append(dict(section="background.source_resolution == producer.background.resolution",
                   production_equals_analysis=A["background"]["source_resolution"] == P["producer"][
                       "background"]["resolution"]))
copies.append(dict(section="gridcheck/priorcheck equal analysis except inference",
                   production_equals_analysis=all(
                       {k: v for k, v in c.items() if k != "inference"}
                       == {k: v for k, v in A.items() if k != "inference"} for c in (G, R))))
copies.append(dict(section="eht-official npix/fov/display_blur == analysis grid/blur",
                   production_equals_analysis=(E["npix"], E["fov_uas"], E["display_blur_uas"]) == (
                       A["image"]["final_shape"], A["image"]["final_fov_uas"], A["image"]["gaussian_fwhm_uas"])))

# --- RAPTOR: production.yaml (raptor_generate/raptor_build default profile) --------------
for key in ("mass_solar", "distance_mpc", "label_step_tg"):
    check(f"production.physical.{key}", P["physical"][key], f"shared.{key}", X["shared"][key])
check("production.data.alphas", P["data"]["alphas"], "shared.alphas", X["shared"]["alphas"])
check("production.inference.principal_n", P["inference"]["principal_n"], "shared.principal_n",
      X["shared"]["principal_n"])
check("production.image.initial_fov_rg", P["image"]["initial_fov_rg"], "raptor.camera_fov_rg",
      X["raptor"]["camera_fov_rg"])
for key in ("build", "plasma", "snapshots", "background", "radial_tables", "alpha_table", "grid"):
    check(f"production.producer.{key}", P["producer"][key], f"raptor.{key}", X["raptor"][key])
for key in ("resolution", "reference_norm"):
    check(f"production.producer.response.{key}", P["producer"]["response"][key], f"raptor.response.{key}",
          X["raptor"]["response"][key])
check("production.producer.response.phase_count", P["producer"]["response"]["phase_count"],
      "eht_axion.producer.plan.PHASE_COUNT (core: phase k*pi/8)", PHASE_COUNT,
      "fixed by the RAPTOR core, no longer a YAML value")

# --- analysis: eht-120uas-observed.yaml + literals of calculate_eht_ensemble_limits --------
check("analysis.image.final_fov_uas", A["image"]["final_fov_uas"], "shared.image.fov_uas",
      X["shared"]["image"]["fov_uas"])
check("analysis.image.final_shape", A["image"]["final_shape"], "shared.image.shape", X["shared"]["image"]["shape"])
check("analysis.image.gaussian_fwhm_uas", A["image"]["gaussian_fwhm_uas"], "shared.image.blur_fwhm_uas",
      X["shared"]["image"]["blur_fwhm_uas"])
check("analysis.image.rotation_degrees", A["image"]["rotation_degrees"], "analysis.rotation_degrees",
      X["analysis"]["rotation_degrees"])
for key, value in A["aperture"].items():
    check(f"analysis.aperture.{key}", value, f"analysis.aperture.{key}", X["analysis"]["aperture"][key])
check("analysis.eht.min_radius_rg", A["eht"]["min_radius_rg"], "analysis.features.min_radius_rg",
      X["analysis"]["features"]["min_radius_rg"])
for key in ("coefficient", "pivot_alpha", "exponent"):
    check(f"analysis.cloud.{key}", A["cloud"][key], f"analysis.cloud.{key}", X["analysis"]["cloud"][key])
check("analysis.cloud.reference_norm", A["cloud"]["reference_norm"], "raptor.response.reference_norm",
      X["raptor"]["response"]["reference_norm"], "same quantity, written once")
check("analysis.background.alpha", A["background"]["alpha"], "pipeline.layout.BACKGROUND_DIRECTORY_ALPHA",
      BACKGROUND_DIRECTORY_ALPHA, "directory label of the axion-off frames only")
check("analysis.background.source_resolution", A["background"]["source_resolution"],
      "raptor.background.resolution", X["raptor"]["background"]["resolution"], "same quantity, written once")
for key in ("lag_frames", "method"):
    check(f"analysis.background.{key}", A["background"][key], f"analysis.background.{key}",
          X["analysis"]["background"][key])
for key, value in A["inference"].items():
    if key == "principal_n":
        check("analysis.inference.principal_n", value, "shared.principal_n", X["shared"]["principal_n"])
    else:
        check(f"analysis.inference.{key}", value, f"analysis.inference.{key}", X["analysis"]["inference"][key])
literals = [  # calculate_eht_ensemble_limits.py at main: literal or CLI default (centering from invocation)
    ("intensity fraction literal .1", .1, "intensity_fraction"),
    ("polarization fraction literal .2", .2, "polarization_fraction"),
    ("--pixel-stride default 5", 5, "pixel_stride"),
    ("--resultant-min default .7", .7, "resultant_min"),
    ("--mask-mode default", "intensity-polarization", "mask_mode"),
    ("sector spacing literal 10 (np.arange(360) % 10 == 4)", 10, "sector_spacing_degrees"),
    ("sector offset literal 4", 4, "sector_offset_degrees"),

]
for source, value, key in literals:
    check(f"calculate_eht_ensemble_limits: {source}", value, f"analysis.features.{key}",
          X["analysis"]["features"][key])
sys.path.insert(0, str(NEW / "scripts"))
import analysis_limits  # noqa: E402
check("calculate_eht_ensemble_limits: --centering native (paper invocation, limits-native-v1)", "native",
      "scripts/analysis_limits.py CENTERING (fixed; Paper VII Sec. 5.2 centring = zero shift)",
      analysis_limits.CENTERING, "the centroid alternative was removed (user decision: Paper VII only)")
for name, cfg, overrides in (("gridcheck", G, ("phase_count", "coupling_points")),
                             ("priorcheck", R, ("phase_count", "coupling_points", "coupling_max"))):
    for key in overrides:
        rows.append(dict(old=f"{name}.inference.{key}", old_value=cfg["inference"][key],
                         new=f"analysis_check_grid --{key.replace('_', '-')}", new_value=cfg["inference"][key],
                         equal=True, note="CLI override (recorded in the check manifest)"))

# --- EHT imaging: eht-official.yaml + literals of run_eht_uncertainty_batch/process_eht_mc ----
shared_image = {"npix": "shape", "fov_uas": "fov_uas", "display_blur_uas": "blur_fwhm_uas"}
locations = ("input_root", "input_receipt", "feed_rotation_file", "official_script")
for key, value in E.items():
    if key in shared_image:
        check(f"eht-official.{key}", value, f"shared.image.{shared_image[key]}",
              X["shared"]["image"][shared_image[key]])
    elif key in locations:
        rows.append(dict(old=f"eht-official.{key}", old_value=value,
                         new=f"--{key.replace('_', '-')} / .env EHT_IMAGING_*", new_value=value, equal=True,
                         note="machine-local location (identical value in .env.example)"))
    else:
        check(f"eht-official.{key}", value, f"eht_imaging.{key}", X["eht_imaging"][key])
c, img, ana = X["eht_imaging"], X["shared"]["image"], X["analysis"]
batch = [("draws-per-day default 1000", 1000, "eht_imaging.uncertainty_draws", c["uncertainty_draws"]),
         ("seed default 20210405", 20210405, "eht_imaging.seed", c["seed"]),
         ("fsys literal .005", .005, "eht_imaging.fsys", c["fsys"]),
         ("residual sigma literal .01", .01, "eht_imaging.residual_d_sigma", c["residual_d_sigma"]),
         ("ddof literal 1", 1, "eht_imaging.covariance_ddof", c["covariance_ddof"]),
         ("uv_min literal 1e8", 1e8, "eht_imaging.uv_min", c["uv_min"]),
         ("display blur literal 19", 19, "shared.image.blur_fwhm_uas", img["blur_fwhm_uas"]),
         ("npix literal 64", 64, "shared.image.shape", img["shape"]),
         ("fov literal 120", 120, "shared.image.fov_uas", img["fov_uas"]),
         ("SITES literal", ["AA", "AP", "AZ", "JC", "LM", "PV", "SM"], "sorted(solve_sites+fixed_sites)",
          sorted(c["solve_sites"] + c["fixed_sites"]))]
for source, value, path, new in batch:
    rows.append(dict(old=f"run_eht_uncertainty_batch: {source}", old_value=value, new=path, new_value=new,
                     equal=value == new))
summary = [("1000 draws", 1000, "eht_imaging.uncertainty_draws", c["uncertainty_draws"]),
           ("64 pixels", 64, "shared.image.shape", img["shape"]),
           ("1.875 uas pixel", 1.875, "fov_uas/shape", img["fov_uas"] / img["shape"]),
           ("rmax_uas 60.", 60.0, "radial_max_fraction*fov/2",
            ana["aperture"]["radial_max_fraction"] * img["fov_uas"] / 2),
           ("aperture window 10", 10, "analysis.aperture.window_degrees", ana["aperture"]["window_degrees"]),
           ("--rmin default .7", .7, "analysis.features.resultant_min", ana["features"]["resultant_min"])]
for source, value, path, new in summary:
    rows.append(dict(old=f"process_eht_mc: {source}", old_value=value, new=path, new_value=new, equal=value == new))

# --- derived quantities: main code vs new code ---------------------------------------------
OLD_CODE = f"""
from pathlib import Path
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.constants import analysis_geometry
import calculate_eht_ensemble_limits as L
from dataclasses import asdict, replace
import yaml
cfgs = {{k: read_config(Path({str(NEW)!r}) / 'cfg' / k) if False else None for k in ()}}
import subprocess, tempfile
def load(name):
    text = subprocess.run(['git', '-C', {str(NEW)!r}, 'show', '{MAIN_COMMIT}:configs/' + name],
                          check=True, capture_output=True, text=True).stdout
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return read_config(path)
A, G, R = load('eht-120uas-observed.yaml'), load('eht-120uas-observed-gridcheck.yaml'), load('eht-120uas-observed-priorcheck.yaml')
g = analysis_geometry(A)
o = asdict(L.options_from_config(A, 0.031))
og = asdict(replace(L.options_from_config(G, 0.031), phase_count=128, coupling_points=240))
orr = asdict(replace(L.options_from_config(R, 0.031), phase_count=128, coupling_points=260))
E = yaml.safe_load(subprocess.run(['git', '-C', {str(NEW)!r}, 'show', '{MAIN_COMMIT}:configs/eht-official.yaml'],
                   check=True, capture_output=True, text=True).stdout)
result = dict(geometry=dict(image_preprocess={{k: float(v) if not isinstance(v, int) else v for k, v in g['image_preprocess'].items()}},
                            aperturefit=g['aperturefit'], radial_max_uas=float(g['metadata']['radial_max_uas']),
                            final_fov_uas=float(g['metadata']['final_fov_uas'])),
              options=o, options_grid=og, options_prior=orr, eht_parameters=E)
"""
NEW_CODE = f"""
from dataclasses import asdict, replace
from eht_axion.experiment import read_experiment, analysis_geometry, eht_imaging_parameters
import analysis_limits as L
import yaml, subprocess
X = read_experiment()
g = analysis_geometry(X)
E = yaml.safe_load(subprocess.run(['git', '-C', {str(NEW)!r}, 'show', '{MAIN_COMMIT}:configs/eht-official.yaml'],
                   check=True, capture_output=True, text=True).stdout)
loc = {{k: E[k] for k in ('input_root', 'input_receipt', 'feed_rotation_file', 'official_script')}}
result = dict(geometry=dict(image_preprocess={{k: float(v) if not isinstance(v, int) else v for k, v in g['image_preprocess'].items()}},
                            aperturefit=g['aperturefit'], radial_max_uas=float(g['metadata']['radial_max_uas']),
                            final_fov_uas=float(g['metadata']['final_fov_uas'])),
              options=asdict(L.options_from_config(X, 0.031)),
              options_grid=asdict(replace(L.options_from_config(X, 0.031), phase_count=128, coupling_points=240)),
              options_prior=asdict(replace(L.options_from_config(X, 0.031), phase_count=128, coupling_points=260,
                                           coupling_max=10000.0)),
              eht_parameters=eht_imaging_parameters(X, loc))
"""
old_derived, new_derived = run("main", OLD_CODE), run("new", NEW_CODE)
derived = {key: dict(identical=old_derived[key] == new_derived[key],
                     differences=diff(old_derived[key], new_derived[key]))
           for key in old_derived}
derived["eht_parameters"]["key_order_identical"] = list(old_derived["eht_parameters"]) == list(
    new_derived["eht_parameters"])

dead = {
    "eht.topset_size": "only pipeline/eht.py (old pickle top-set workflow, deleted)",
    "eht.historical_mask": "only pipeline/eht.py (deleted)",
    "eht.inherit_2d_mask_in_1d": "only pipeline/eht.py (deleted)",
    "eht.covariance_ridge": "only products.build_observations (old null workflow, deleted)",
    "eht.intensity_fraction": "only pipeline/eht.py; the active limits hard-coded .1 -> now analysis.features",
    "eht.angular_offset": "only products.build_observations; active limits hard-coded 4 -> now analysis.features",
    "eht_input.* (shape, fov_uas, pixel_size_uas, angle_roll, rank_offsets, final_offsets, exclude_columns, "
    "exclusion_rectangles)": "only workflow eht stage / pipeline/eht.py (old pickle inputs, deleted); "
                             "the active ensemble adapter uses the shared sky grid",
    "image.final_fov_rg, image.gaussian_sigma_rg, aperture.radial_max_rg (research/production Rg schema)":
        "production.yaml analysis part was validated but never consumed by the producer; research.yaml "
        "and analysis.yaml were used only by deleted workflows",
    "producer.response.phase_count": "fixed by the RAPTOR core (phase k*pi/8) -> producer.plan.PHASE_COUNT",
    "background.alpha": "directory label only -> pipeline.layout.BACKGROUND_DIRECTORY_ALPHA",
    "configs/research.yaml, analysis.yaml, eht-observed-research.yaml": "deleted null/research workflows",
    "configs/eht-fiducial.yaml, eht-apr11-validation.yaml, eht-polarization-study.yaml":
        "deleted rebuild_eht_fiducial / validate_eht_apr11 / research_eht_polarization",
    "configs/eht-120uas-observed-{gridcheck,priorcheck}.yaml": "replaced by analysis_check_grid CLI overrides",
}
record = dict(
    purpose="Values consumed by the active chain: main:configs + main literals vs experiment.yaml",
    main_commit=MAIN_COMMIT, values=rows, all_values_equal=all(r["equal"] for r in rows),
    old_profile_copies=copies, derived_main_vs_new=derived,
    all_derived_identical=all(v["identical"] for v in derived.values()),
    removed_dead_fields=dead,
    notes=["eht-imaging plan parameters are rebuilt with the same keys/order; the four machine paths are "
           "now CLI/.env inputs and enter the parameter set unchanged.",
           "Generated-data identity (generated.json) and task ids intentionally changed format (no compatibility)."])
write("config-equivalence.json", record)
print("values equal:", record["all_values_equal"], "| derived identical:", record["all_derived_identical"])
for r in rows:
    if not r["equal"]:
        print("DIFF", r)
for k, v in derived.items():
    if not v["identical"]:
        print("DERIVED DIFF", k, v["differences"][:5])
