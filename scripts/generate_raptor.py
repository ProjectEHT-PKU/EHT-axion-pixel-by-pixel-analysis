"""Plan/run/resume RAPTOR jobs or deliver verified HDF; never implicitly build."""
import json
import shutil

import _bootstrap  # noqa: F401
from _producer import options, save_request
from eht_axion.producer.delivery import deliver_stokes
from eht_axion.producer.jobs import run_jobs
from eht_axion.producer.plan import plan_generation


def main(argv=None):
    args = options("generate", argv)
    plan = plan_generation(args.config_data, args.fluid_root,
                           {"low": args.radial_low, "high": args.radial_high},
                           args.run_root / "build")
    plan["input_roots"] = [str(args.fluid_root), str(args.radial_low.parent), str(args.radial_high.parent)]
    plan["source_roots"] = [str(args.source_root)]
    plan["expected_core_commit"] = args.expected_commit
    summary = {"counts": plan["counts"], "requirements": plan["requirements"],
               "warnings": plan["warnings"], "run_root": str(args.run_root),
               "scheduling": {"workers": args.workers, "threads_per_worker": args.threads,
                              "affinity": args.cpu_affinity, "cpu_sets": args.cpu_sets}}
    print(json.dumps(summary, indent=2, allow_nan=False))
    if args.plan:
        return 0
    # Numerical arrays only are a lower bound. Require additional disk headroom;
    # job RSS/elapsed estimates remain unmeasured until an authorized pilot.
    existing = args.run_root
    while not existing.exists():
        existing = existing.parent
    estimate = plan["requirements"]["output_stokes_lower_bound_bytes"]
    required = 4 * estimate + 2 * 1024**3
    if shutil.disk_usage(existing).free < required:
        raise ValueError(f"Insufficient free space: require at least {required} bytes "
                         "for conservative headroom")
    if not args.convert_only:
        from eht_axion.producer.affinity import check_memory

        per_worker = max(int(2.5 * 1024**3),
                         2 * plan["requirements"]["per_process_memory_lower_bound_bytes"])
        check_memory(args.workers, per_worker)
    save_request(args, plan)
    if not args.convert_only:
        result = run_jobs(plan, args.run_root / "build", args.run_root, max_jobs=args.max_jobs,
                          threads=args.threads, workers=args.workers,
                          cpu_affinity=args.cpu_affinity, cpu_list=args.cpu_list,
                          timeout_seconds=args.timeout_seconds)
        print(json.dumps(result["counts"], indent=2))
        if not result["complete"]:
            return 1
        if args.max_jobs is not None and args.max_jobs < len(plan["tasks"]):
            print("Selected task prefix completed; full delivery waits for the remaining jobs.")
            return 0
    products = deliver_stokes(plan, args.run_root)
    print(f"Verified Stokes data: {args.run_root / 'data'} ({len(products)} products)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
