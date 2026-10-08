"""Generate Figure-4 transport diagnostics from materialized NPY inputs.

No notebook or legacy ImageS is used; see docs/transport-stage.md for conventions.
"""
from __future__ import annotations

import _bootstrap  # noqa: F401

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from eht_axion.pipeline.transport import ORDERS, run_transport_diagnostics


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="directory containing the seven materialized .npy files")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--radii", nargs="+", type=float, default=[4, 5, 6, 7, 8, 9, 10])
    parser.add_argument("--r-width", type=float, default=0.5)
    parser.add_argument("--fov", type=float, default=46.2262)
    parser.add_argument("--sigma", type=float, default=2.16871, help="Gaussian stddev in native FOV rg units")
    parser.add_argument("--rotation-angle", type=float, default=0, help="radians; same for all curves")
    parser.add_argument("--angular-bins", type=int, default=360)
    parser.add_argument("--mask-order", choices=ORDERS, default=ORDERS[0])
    parser.add_argument("--amplitude-method", choices=("legacy-rms", "dft"), default="legacy-rms")
    parser.add_argument("--angular-interpolation", choices=("spin2-periodic", "legacy-angle"),
                        default="spin2-periodic")
    parser.add_argument("--crossing", choices=("first-outward", "legacy-bracket"), default="first-outward")
    parser.add_argument("--spin", type=float, default=0.9375)
    parser.add_argument("--denominator-floor", type=float, default=1e-12)
    try:
        products = run_transport_diagnostics(**vars(parser.parse_args(argv)))
    except (ValueError, OSError, TypeError) as error:
        parser.error(str(error))
    print(json.dumps({name: str(path) for name, path in products.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
