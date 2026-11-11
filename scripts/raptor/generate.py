"""Plan/run/resume RAPTOR jobs and write their Stokes image cubes into raptor/images; never build."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import json
import shutil

from eht_axion.raptor.options import options
from eht_axion.raptor.images import write_images
from eht_axion.raptor.jobs import run_jobs
from eht_axion.raptor.plan import plan_generation


def frozen_plan(args):
    """raptor/plan.json: written once (snapshots hashed then); later runs must plan the same tasks."""
    inputs = (args.config_data, args.fluid_root, {"low": args.radial_low, "high": args.radial_high},
              args.run_root / "build")
    path = args.run_root / "plan.json"
    plan = plan_generation(*inputs)
    roots = (args.fluid_root, args.radial_low.parent, args.radial_high.parent)
    extra = dict(input_roots=[str(p) for p in roots],
                 source_roots=[str(args.source_root)],
                 raptor_source=dict(args.config_data["raptor"]["source"]))
    plan.update(extra)
    if path.exists():
        stored = json.loads(path.read_text())
        light = json.loads(json.dumps(stored))
        light["hash_inputs"] = False
        for record in light["snapshots"] + light["tasks"]:
            record["source_sha256"] = None
        if light != json.loads(json.dumps(plan)):
            raise ValueError(f"{path} plans different tasks; use a new data root")
        return stored, False
    if args.plan:
        return plan, True
    plan = dict(plan_generation(*inputs, hash_inputs=True), **extra)
    return plan, True


def main(argv=None):
    args = options("generate", argv)
    plan, new = frozen_plan(args)
    summary = {"counts": plan["counts"], "requirements": plan["requirements"],
               "warnings": plan["warnings"], "run_root": str(args.run_root),
               "scheduling": {"workers": args.workers, "threads_per_worker": args.threads,
                              "affinity": args.cpu_affinity, "cpu_sets": args.cpu_sets}}
    print(json.dumps(summary, indent=2, allow_nan=False))
    if args.plan:
        return 0
    # Numerical arrays only are a lower bound. Require additional disk headroom;
    # job RSS/elapsed time are not estimated here.
    existing = args.run_root
    while not existing.exists():
        existing = existing.parent
    estimate = plan["requirements"]["output_stokes_lower_bound_bytes"]
    required = 4 * estimate + 2 * 1024**3
    if shutil.disk_usage(existing).free < required:
        raise ValueError(f"Insufficient free space: require at least {required} bytes "
                         "for conservative headroom")
    if new:
        args.run_root.mkdir(parents=True, exist_ok=True)
        with (args.run_root / "plan.json").open("x") as stream:
            json.dump(plan, stream, indent=2, allow_nan=False)
            stream.write("\n")
    if not args.convert_only:
        from eht_axion.raptor.affinity import check_memory

        per_worker = max(int(2.5 * 1024**3),
                         2 * plan["requirements"]["per_process_memory_lower_bound_bytes"])
        check_memory(args.workers, per_worker)
        result = run_jobs(plan, args.run_root / "build", args.run_root, max_jobs=args.max_jobs,
                          threads=args.threads, workers=args.workers,
                          cpu_affinity=args.cpu_affinity, cpu_list=args.cpu_list,
                          timeout_seconds=args.timeout_seconds)
        print(json.dumps(result["counts"], indent=2))
        if not result["complete"]:
            return 1
        if args.max_jobs is not None and args.max_jobs < len(plan["tasks"]):
            print("Selected task prefix completed; the image cubes wait for the remaining jobs.")
            return 0
    products = write_images(plan, args.run_root)
    print(f"Stokes cubes: {args.run_root / 'images'} ({len(products)} products)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
