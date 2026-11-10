"""Deliver Stokes cubes from the completed job roots of several hosts into a new root; runs no RAPTOR."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import _bootstrap  # noqa: F401
from _data_root import DATA_ROOT_KEY, add_data_root_option, data_root, fill_defaults, layout_path
from eht_axion.provenance import manifest_path, stage_directory, verification_session
from eht_axion.experiment import validate_ray_tracing
from eht_axion.raptor.delivery import deliver_stokes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, help="Frozen task-plan JSON (default: raptor/plan.json)")
    parser.add_argument("--job-root", type=Path, action="append",
                        help="Source run containing jobs/ and its real build/; repeat for independent runs "
                             "(default: every directory in raptor/runs)")
    parser.add_argument("--output-root", type=Path,
                        help="New delivery container, disjoint from the source job roots "
                             "(default: raptor/delivered)")
    add_data_root_option(parser)
    args = parser.parse_args(argv)
    fill_defaults(parser, args, {"plan": ("raptor_plan", {}), "output_root": ("raptor_delivered", {})})
    if args.job_root is None:
        root = data_root(args)
        if root is None:
            parser.error(f"Provide --job-root, or --data-root / ${DATA_ROOT_KEY}")
        args.job_root = sorted(p for p in layout_path(root, "raptor_runs").iterdir() if p.is_dir())
        if not args.job_root:
            parser.error(f"No run roots in {layout_path(root, 'raptor_runs')}")
    plan = json.loads(args.plan.read_text())
    validate_ray_tracing(plan["config"])
    with verification_session():
        products = deliver_stokes(plan, args.job_root[0], job_roots=args.job_root,
                                  output_root=args.output_root)
        inputs = [args.plan.resolve()]
        for path in products:
            inputs += [manifest_path(path.parent), path, path.parent / "generated.json",
                       path.parent / "sources.json"]
        summary = {"layout": "raptor-merge-v2", "tasks": len(plan["tasks"]),
                   "source_roots": [str(p.resolve()) for p in args.job_root],
                   "science": plan["config"],
                   "products": [str(p.resolve()) for p in products]}
        with stage_directory(args.output_root / "inventory", "raptor:merge",
                             summary, inputs) as temporary:
            if temporary is not None:
                (temporary / "index.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"MERGE_COMPLETE: {len(products)} delivered products in {args.output_root / 'data'}")
    print("Job directories, input snapshots and manifests are only read; no RAPTOR run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
