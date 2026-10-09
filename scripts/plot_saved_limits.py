"""Direct research script; numerical functions remain in reusable modules."""

import argparse
from pathlib import Path
import _bootstrap  # noqa: F401
from eht_axion.pipeline.figures import plot_results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reuse-imported", action="store_true")
    args = parser.parse_args(argv)
    plot_results(args.input_dir, args.output_dir, reuse_imported=args.reuse_imported)


if __name__ == "__main__":
    main()
