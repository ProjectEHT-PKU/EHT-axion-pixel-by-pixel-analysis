"""Direct research script; numerical functions remain in reusable modules."""

import argparse
from pathlib import Path
import _bootstrap  # noqa: F401
from eht_axion.pipeline.diagnostics import plot_diagnostics


def main(argv=None):

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    plot_diagnostics(args.run_root, args.output_dir)


if __name__ == "__main__":
    main()
