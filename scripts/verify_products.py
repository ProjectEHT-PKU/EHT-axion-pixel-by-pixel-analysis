"""Verify one completed physical stage; this script never recomputes data."""

import argparse
import json
from pathlib import Path
import _bootstrap  # noqa: F401
from eht_axion.pipeline.common import verify_stage


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--recursive", action="store_true", help="Also validate upstream lineage")
    args = parser.parse_args(argv)
    print(json.dumps(verify_stage(args.directory, recursive=args.recursive), indent=2))


if __name__ == "__main__":
    main()
