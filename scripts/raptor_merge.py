"""Merge completed multi-host RAPTOR jobs into canonical Stokes NPY; never launch C."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import _bootstrap  # noqa: F401
from eht_axion.provenance import manifest_path, stage_directory, verification_session
from eht_axion.experiment import validate_ray_tracing
from eht_axion.raptor.delivery import deliver_stokes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True, help="Frozen production-plan JSON")
    parser.add_argument("--job-root", type=Path, action="append", required=True,
                        help="Source run containing jobs/ and its real build/; repeat for independent runs")
    parser.add_argument("--output-root", type=Path, required=True,
                        help="New delivery container, disjoint from original source runs")
    args = parser.parse_args(argv)
    plan = json.loads(args.plan.read_text())
    validate_ray_tracing(plan["config"])
    with verification_session():
        products = deliver_stokes(plan, args.job_root[0], job_roots=args.job_root,
                                  output_root=args.output_root)
        inputs = [args.plan.resolve()]
        for path in products:
            inputs += [manifest_path(path.parent), path, path.parent / "generated.json",
                       path.parent / "sources.json"]
        summary = {"layout": "raptor-multi-source-delivery-v1", "tasks": len(plan["tasks"]),
                   "source_roots": [str(p.resolve()) for p in args.job_root],
                   "science": plan["config"],
                   "products": [str(p.resolve()) for p in products]}
        with stage_directory(args.output_root / "inventory", "producer:merged-inventory",
                             summary, inputs) as temporary:
            if temporary is not None:
                (temporary / "index.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"MERGE_COMPLETE: {len(products)} canonical products in {args.output_root / 'data'}")
    print("Original HDF5, input snapshots and manifests were preserved; no science rerun.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
