"""Analysis: the expected axion signal and the GRMHD background from the RAPTOR image cubes.

Reads raptor/images (scripts/raptor/generate.py) and publishes ONE stage (analysis:expected-signal)
in a new output root. For each alpha, snapshot by snapshot and without storing anything in
between, the 16 axion-phase images are blurred, rotated and resampled onto the shared sky grid
and their EVPA is fitted with a first harmonic over the phase (per pixel and per 10-degree sector);
the mean and scatter over the snapshots are the expected signal amplitude A(x) and phase phi(x)
per RAPTOR normalisation (raptor.response.reference_norm).

  expected.npz     1d/... (sectors, 360 one-degree bins) and 2d/... (pixels): alphas,
                   amplitude_mean, amplitude_std, phase_mean, phase_std, amplitude_count, phase_count
  background.npy   axion-off background images on the sky grid (T, 1, N, N, 4)

The GRMHD background sigma is computed by scripts/analysis/constraints.py on exactly the features
used in the likelihood, from background.npy. No ray tracing, no EHT data.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse

import numpy as np

from eht_axion.analysis.expected_signal import (
    _operator, alpha_responses, expected_signal, processed_background,
)
from eht_axion.environment import add_data_root_option, fill_defaults
from eht_axion.experiment import DEFAULT_PATH, analysis_geometry, read_experiment
from eht_axion.provenance import reject_directory_symlinks, stage_directory
from eht_axion.raptor.images import generated_inputs
from eht_axion.raptor.layout import image_cube


def stage_config(config, labels):
    """Recorded configuration of analysis:expected-signal (geometry of the current experiment)."""
    geometry = analysis_geometry(config)
    return {"image_preprocess": geometry["image_preprocess"], "aperturefit": geometry["aperturefit"],
            "alphas": list(config["shared"]["alphas"]), "snapshot_labels": labels,
            "state0_used": False, "phase_count": 16, "minimum_samples": 2,
            "method": "first-harmonic EVPA fit over the 16 phases per snapshot; Gaussian MLE (ddof=0) "
                      "over the snapshots, phases lifted about their circular mean",
            "units": "rad per raptor.response.reference_norm"}


def image_inputs(config, images_root):
    """Verified RAPTOR image cubes of every alpha and of the background, their stage inputs and labels."""
    response_resolution = config["raptor"]["response"]["resolution"]
    background_resolution = config["raptor"]["background"]["resolution"]
    cubes, inputs = {}, []
    for alpha in config["shared"]["alphas"]:
        cubes[alpha] = image_cube(images_root, alpha, response_resolution)
        found, labels = generated_inputs(cubes[alpha], config, alpha=alpha, resolution=response_resolution)
        inputs += found
    cubes[None] = image_cube(images_root, None, background_resolution)
    inputs += generated_inputs(cubes[None], config, alpha=None, resolution=background_resolution)[0]
    return cubes, inputs, labels


def run(config_path, images_root, output_root):
    config = read_experiment(config_path)
    geometry = analysis_geometry(config)
    preprocess, aperture = geometry["image_preprocess"], geometry["aperturefit"]
    images_root = Path(images_root).resolve(strict=True)
    reject_directory_symlinks(output_root)
    output_root = Path(output_root).resolve()
    if output_root.exists():
        raise FileExistsError(f"New output root required: {output_root}")
    if output_root.is_relative_to(images_root) or images_root.is_relative_to(output_root):
        raise ValueError("Output root must not overlap the RAPTOR images")
    cubes, inputs, labels = image_inputs(config, images_root)
    operator = _operator(preprocess["final_shape"], aperture["fov"], aperture["radial_max"],
                         aperture["window_degrees"], aperture["angular_bins"])
    stage = stage_config(config, labels)
    with stage_directory(output_root, "analysis:expected-signal", stage, inputs) as out:
        pixels, sectors = {}, {}
        for alpha in config["shared"]["alphas"]:
            cube = np.load(cubes[alpha], mmap_mode="r", allow_pickle=False)
            pixels[alpha], sectors[alpha] = alpha_responses(cube, preprocess, operator)
            print(f"EXPECTED SIGNAL alpha={alpha:.2f}", flush=True)
        arrays = {}
        for dim, responses in (("1d", sectors), ("2d", pixels)):
            _, products = expected_signal(responses)
            arrays.update({f"{dim}/{name}": value for name, value in products.items()})
        np.savez(out / "expected.npz", **arrays)
        background = np.load(cubes[None], mmap_mode="r", allow_pickle=False)
        np.save(out / "background.npy", processed_background(background, preprocess))
    return output_root


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=DEFAULT_PATH, help="default: configs/experiment.yaml")
    parser.add_argument("--raptor-images", type=Path, help="RAPTOR image cubes (default: raptor/images)")
    parser.add_argument("--output-root", type=Path, help="new directory (default: analysis/expected-signal)")
    add_data_root_option(parser)
    args = parser.parse_args(argv)
    fill_defaults(parser, args, {"raptor_images": ("raptor.images", {}),
                                 "output_root": ("analysis.expected_signal", {})})
    run(args.config, args.raptor_images, args.output_root)


if __name__ == "__main__":
    main()
