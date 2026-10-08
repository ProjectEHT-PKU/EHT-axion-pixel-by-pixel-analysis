"""Recompute constraints with denser numerical grids, reusing physical products."""

from dataclasses import replace
from pathlib import Path
import json
import time
import numpy as np
from ..limits import LimitOptions, calculate_from_products
from .products import load_observations, load_templates
from .common import verify_stage
from .layout import RunLayout
from .figures import MODES, load_limits


def check_convergence(
    run_root,
    output_path,
    *,
    dimensions=("1d", "2d"),
    modes=("std", "cor", "cor_grmhd"),
    phase_multiplier=2,
    grid_multiplier=2,
):
    layout = RunLayout(run_root)
    verify_stage(layout.inference)
    saved = json.loads((layout.inference / "run.json").read_text())
    base = LimitOptions(**saved["options"])
    report = {"reference": str(layout.inference.resolve()), "checks": {}}
    reference, identity = None, None
    for dim in dimensions:
        alphas, amplitude, phase, _ = load_templates(layout.templates(dim) / "templates.npz", dim)
        obs = load_observations(layout.observations(dim) / "observations.npz")
        old_alphas, old_cuts, metadata = load_limits(layout.inference, dim)
        if not np.array_equal(alphas, old_alphas) or (
            reference is not None and not np.array_equal(reference, old_alphas)
        ):
            raise ValueError("Alpha grid mismatch between inference and template products")
        if (
            metadata["observation"] != ("observed" if base.observed else "null")
            or metadata["prior"] != base.prior
            or (identity is not None and identity != metadata)
        ):
            raise ValueError("Mixed observations, priors or algorithm variants")
        reference, identity = old_alphas, metadata
        for mode in modes:
            stem = f"{mode}_{dim}"
            old = old_cuts[MODES.index(mode)]
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
