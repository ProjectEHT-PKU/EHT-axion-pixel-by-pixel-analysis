"""Direct research script; numerical functions remain in reusable modules."""

import argparse
from pathlib import Path
import _bootstrap  # noqa: F401
from eht_axion.pipeline.convergence import check_convergence


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modes", nargs="+", choices=["std", "cor", "cor_grmhd"], default=["std", "cor"])
    parser.add_argument("--dimensions", nargs="+", choices=["1d", "2d"], default=["1d", "2d"])
    args = parser.parse_args(argv)
    check_convergence(args.run_root, args.output, dimensions=args.dimensions, modes=args.modes)


if __name__ == "__main__":
    main()
