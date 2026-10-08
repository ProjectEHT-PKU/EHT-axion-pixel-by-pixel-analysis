"""Direct research script; numerical functions remain in reusable modules."""

import argparse
from pathlib import Path
import _bootstrap  # noqa: F401
from eht_axion.pipeline.migrate import migrate_inputs, migrate_run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["inputs", "run"])
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args(argv)
    action = migrate_inputs if args.kind == "inputs" else migrate_run
    print(action(args.source, args.destination))


if __name__ == "__main__":
    main()
