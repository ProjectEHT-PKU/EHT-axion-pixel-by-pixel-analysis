"""Analysis: response templates and the GRMHD background from delivered RAPTOR data.

Reads the delivered Stokes cubes (scripts/raptor_merge.py or raptor_generate.py)
and writes, into a NEW output root, one verified stage per step:

  processed/alpha-*      analysis:preprocess   blur, rotate, resample to the sky grid
  response/alpha-*       analysis:fit          per-snapshot first-harmonic EVPA response
  templates/{1d,2d}      analysis:templates    template mean/std over the snapshots
  background             analysis:background   axion-off frames on the sky grid

The GRMHD background sigma is computed by analysis_limits.py on exactly the
features used in the likelihood, from background/. No ray tracing, no EHT data.
"""
from pathlib import Path
import argparse
import json

import _bootstrap  # noqa: F401
from eht_axion.analysis.templates import RunLayout, fit_processed, preprocess_cube, summarize_fits
from eht_axion.experiment import DEFAULT_PATH, analysis_geometry, read_experiment
from eht_axion.provenance import (
    manifest_path, reject_directory_symlinks, stage_directory, verification_session, verify_stage,
)
from eht_axion.raptor.delivery import generated_inputs
from eht_axion.raptor.layout import input_raptor_files


def stage_inputs(directory):
    """A completed stage as an input: its manifest plus every declared output."""
    manifest = verify_stage(directory)
    return [manifest_path(directory)] + [Path(directory) / p["path"] for p in manifest["outputs"]]


def run(config_path, raptor_root, output_root):
    config = read_experiment(config_path)
    geometry = analysis_geometry(config)
    preprocess, aperture = geometry["image_preprocess"], geometry["aperturefit"]
    raptor_root = Path(raptor_root).resolve(strict=True)
    reject_directory_symlinks(output_root)
    output_root = Path(output_root).resolve()
    if output_root.exists():
        raise FileExistsError(f"New output root required: {output_root}")
    if output_root.is_relative_to(raptor_root) or raptor_root.is_relative_to(output_root):
        raise ValueError("Output root must not overlap the RAPTOR input root")
    layout = RunLayout(output_root)
    alphas = config["shared"]["alphas"]
    response_resolution = config["raptor"]["response"]["resolution"]
    background_resolution = config["raptor"]["background"]["resolution"]
    with verification_session():
        verified = set()
        for alpha in alphas:
            cubes = input_raptor_files(raptor_root, alpha, response_resolution)
            inputs, labels = generated_inputs(cubes, config, alpha=alpha, resolution=response_resolution,
                                              verified=verified)
            with stage_directory(layout.processed(alpha), "analysis:preprocess", preprocess, inputs) as out:
                if out is not None:
                    preprocess_cube(cubes, out, **preprocess)
            fit = dict(aperture, snapshot_labels=labels)
            with stage_directory(layout.response(alpha), "analysis:fit", fit,
                                 stage_inputs(layout.processed(alpha))) as out:
                if out is not None:
                    fit_processed(layout.processed(alpha), out, **fit)
        for dim in ("1d", "2d"):
            fits = {a: layout.response(a) / f"response-{dim}.npz" for a in alphas}
            inputs = [p for a in alphas for p in (manifest_path(layout.response(a)), fits[a])]
            with stage_directory(layout.templates(dim), "analysis:templates",
                                 {"dimension": dim, "method": "mle"}, inputs) as out:
                if out is not None:
                    summarize_fits(fits, out)

        cubes = input_raptor_files(raptor_root, None, background_resolution)
        inputs, _ = generated_inputs(cubes, config, alpha=None,
                                     resolution=background_resolution, verified=verified)
        params = dict(preprocess, state_indices=[0])
        with stage_directory(layout.background_processed, "analysis:background",
                             params, inputs) as out:
            if out is not None:
                preprocess_cube(cubes, out, **params)
    summary = dict(experiment=str(Path(config_path).resolve()), raptor_input_root=str(raptor_root),
                   output_root=str(output_root), geometry=geometry)
    print(json.dumps(summary, indent=2), flush=True)
    return output_root


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=DEFAULT_PATH, help="default: experiment.yaml")
    parser.add_argument("--raptor-input-root", type=Path, required=True,
                        help="delivered data root containing raptor/res-*/")
    parser.add_argument("--output-root", type=Path, required=True, help="new directory")
    args = parser.parse_args(argv)
    run(args.config, args.raptor_input_root, args.output_root)


if __name__ == "__main__":
    main()
