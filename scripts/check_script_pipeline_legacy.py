#!/usr/bin/env python3
"""Validate the preserved legacy numerical branch against saved historical curves."""

import _bootstrap  # noqa: F401
import argparse
import json
from pathlib import Path
import numpy as np
from eht_axion.limits import LimitOptions, calculate_limits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--historical-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {}
    for dim in ("1d", "2d"):
        for mode in ("std", "cor", "cor_grmhd"):
            alpha, cuts = calculate_limits(
                args.historical_root, dimension=dim, mode=mode,
                options=LimitOptions(legacy=True, n_principal=4),
            )
            suffix = "_1d" if dim == "1d" else ""
            with np.load(args.historical_root / "Results" / "data" / f"{mode}{suffix}.npz") as z:
                np.testing.assert_allclose(alpha, z["alphas"], rtol=0, atol=1e-15)
                difference = float(np.max(abs(cuts - z["cuts"])))
            if difference >= 1e-12:
                raise ValueError(f"Legacy regression failed: {mode} {dim}: {difference}")
            report[f"{mode}_{dim}"] = difference
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
