"""Verify one completed stage (and with --recursive its whole lineage); never recomputes data.

All checks run in one verification session: every file is hashed once, however many
stages declare it as an input, and must be unchanged when the check ends.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse
import json
from eht_axion.provenance import verification_session, verify_stage


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--recursive", action="store_true", help="Also validate upstream lineage")
    args = parser.parse_args(argv)
    with verification_session():
        manifest = verify_stage(args.directory, recursive=args.recursive)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
