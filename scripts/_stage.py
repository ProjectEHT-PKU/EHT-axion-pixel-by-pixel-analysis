"""Shared path/config options for the small research-stage scripts."""

from pathlib import Path
import argparse

import _bootstrap  # noqa: F401
from eht_axion.pipeline.workflow import run_stage


def execute(stage, argv=None):
    parser = argparse.ArgumentParser(description=f"Run the {stage} research stage (no package installation).")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--trust-pickle", action="store_true")
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--reuse-imported", action="store_true")
    args = parser.parse_args(argv)
    output = run_stage(
        args.config,
        args.output_root,
        input_root=args.input_root,
        stage=stage,
        trust_pickle=args.trust_pickle,
        rebuild=args.rebuild,
        reuse_imported=args.reuse_imported,
    )
    print(f"Results: {output}")
    return output
