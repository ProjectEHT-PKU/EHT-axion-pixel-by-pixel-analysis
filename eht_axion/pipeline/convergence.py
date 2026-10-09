"""Recompute constraints with denser numerical grids, reusing physical products."""

from dataclasses import replace
from pathlib import Path
import json
import time
import numpy as np
from ..limits import LimitOptions, calculate_from_products
from .products import load_observations, load_templates
from .common import verify_stage, manifest_path, require_new_output_file, require_declared_outputs
from .inference import _check_lineage, validate_run_lineage
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
    reuse_imported=False,
):
    for value in (phase_multiplier, grid_multiplier):
        if isinstance(value, bool) or not isinstance(value, int) or value < 2:
            raise ValueError("Grid multipliers must be integers >= 2")
    if (not dimensions or any(dim not in ("1d", "2d") for dim in dimensions)
            or not modes or any(mode not in MODES for mode in modes)):
        raise ValueError("Invalid convergence dimensions/modes")
    require_new_output_file(output_path)
    layout = RunLayout(run_root)
    seen = set()
    manifest = verify_stage(layout.inference, recursive=True, allow_imported=reuse_imported, _seen=seen)
    require_declared_outputs(layout.inference, [layout.inference / "run.json"])
    saved = json.loads((layout.inference / "run.json").read_text())
    recorded = dict(saved["options"])
    legacy_flag = recorded.pop("observed", False)
    if legacy_flag is not False:
        raise ValueError("Measured-day-difference runs are not supported")
    base = LimitOptions(**recorded)
    report = {"reference": str(layout.inference.resolve()), "checks": {}}
    reference, identity = None, None
    required_inputs = [manifest_path(layout.background), layout.background / "summary.json"]
    for dim in ("1d", "2d"):
        for stage in (layout.templates(dim), layout.observations(dim)):
            verify_stage(stage, recursive=True, allow_imported=reuse_imported, _seen=seen)
        required_inputs += [
            manifest_path(layout.observations(dim)), manifest_path(layout.templates(dim)),
            layout.observations(dim) / "observations.npz", layout.templates(dim) / "templates.npz",
        ]
    _check_lineage(manifest, required_inputs)
    recorded_alphas, _, _ = load_limits(layout.inference, "1d")
    validate_run_lineage(run_root, recorded_alphas, allow_imported=reuse_imported)
    for dim in dimensions:
        alphas, amplitude, phase, _ = load_templates(layout.templates(dim) / "templates.npz", dim)
        obs = load_observations(layout.observations(dim) / "observations.npz")
        old_alphas, old_cuts, metadata = load_limits(layout.inference, dim)
        if not np.array_equal(alphas, old_alphas) or (
            reference is not None and not np.array_equal(reference, old_alphas)
        ):
            raise ValueError("Alpha grid mismatch between inference and template products")
        if (
            metadata["observation"] not in {"null", "expected"}
            or metadata["prior"] != base.prior
            or (identity is not None and identity != metadata)
        ):
            raise ValueError("Mixed observations, priors or algorithm variants")
        reference, identity = old_alphas, metadata
        for mode in modes:
            stem = f"{mode}_{dim}"
            old = old_cuts[MODES.index(mode)]
            _, baseline, _ = calculate_from_products(obs, amplitude, phase, alphas, dim, mode, base)
            baseline_delta = float(np.max(abs(baseline - old)))
            if not np.allclose(baseline, old, rtol=0, atol=1e-10):
                raise ValueError(
                    f"Saved baseline differs from current model/algorithm ({stem}, "
                    f"{baseline_delta:.6g} dex); recompute inference before checking grids"
                )
            entry = {"baseline_max_delta_dex": baseline_delta}
            resolved_phases = base.phase_count or (16 if base.legacy else 64)
            for name, options in [
                ("coupling", replace(base, coupling_points=base.coupling_points * grid_multiplier)),
                ("phase", replace(base, phase_count=resolved_phases * phase_multiplier)),
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
