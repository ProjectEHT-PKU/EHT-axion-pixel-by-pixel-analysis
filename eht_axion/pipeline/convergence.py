"""Recompute constraints with denser numerical grids, reusing physical products."""

from dataclasses import replace
from pathlib import Path
import argparse
import json
import time
import numpy as np
from ..limits import LimitOptions, calculate_from_products
from .products import load_observations, load_templates
from .common import verify_stage


def check_convergence(
    run_root,
    output_path,
    *,
    dimensions=("1d", "2d"),
    modes=("std", "cor", "cor_grmhd"),
    phase_multiplier=2,
    grid_multiplier=2,
):
    root = Path(run_root)
    verify_stage(root / "07_inference")
    saved = json.loads((root / "07_inference" / "run.json").read_text())
    base = LimitOptions(**saved["options"])
    report = {"reference": str((root / "07_inference").resolve()), "checks": {}}
    for dim in dimensions:
        alphas, amplitude, phase, _ = load_templates(root / "03_templates" / dim / "templates.npz", dim)
        obs = load_observations(root / "06_observations" / dim / "observations.npz")
        for mode in modes:
            stem = f"{mode}_{dim}"
            old = np.array(saved["results"][stem]["cuts"])
            entry = {}
            for name, options in [
                ("coupling", replace(base, coupling_points=base.coupling_points * grid_multiplier)),
                ("phase", replace(base, phase_count=base.phase_count * phase_multiplier)),
            ]:
                started = time.perf_counter()
                _, cuts, _ = calculate_from_products(obs, amplitude, phase, alphas, dim, mode, options)
                entry[name] = {
                    "cuts": cuts.tolist(),
                    "max_delta_dex": float(np.max(abs(cuts - old))),
                    "elapsed_seconds": time.perf_counter() - started,
                    "coupling_points": options.coupling_points,
                    "phase_count": options.phase_count,
                }
                print(stem, name, entry[name]["max_delta_dex"], flush=True)
            report["checks"][stem] = entry
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modes", nargs="+", choices=["std", "cor", "cor_grmhd"], default=["std", "cor"])
    parser.add_argument("--dimensions", nargs="+", choices=["1d", "2d"], default=["1d", "2d"])
    args = parser.parse_args()
    check_convergence(args.run_root, args.output, dimensions=args.dimensions, modes=args.modes)


if __name__ == "__main__":
    main()
