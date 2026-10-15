"""Read-only compare principal-N=2 leading frequencies with existing source CSVs.

Run from repository root:
  pixi run python reports/script-pipeline/compare_211_frequencies.py \
      --radial-root ../radial-function --output-dir /tmp/211-comparison
No eigenvalue/radial/transfer calculation is run. Output directories must be new.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ALPHAS = ((0.1, 0.9375), (0.15, 0.9375), (0.2, 0.9375), (0.25, 0.9375),
          (0.3, 0.9375), (0.35, 0.99), (0.4, 0.99), (0.44, 0.99))


def compare(radial_root):
    sources, tables = {}, {}
    for spin in (0.9375, 0.99):
        path = Path(radial_root) / "source" / f"a={spin}_l1m1_real.txt"
        raw = path.read_bytes()
        declared = int(raw.splitlines()[0])
        table = np.loadtxt(path, delimiter=",", skiprows=1)
        if (table.shape != (declared, 2) or not np.isfinite(table).all()
                or not np.all(np.diff(table[:, 0]) > 0)):
            raise ValueError(f"Invalid source frequency table: {path}")
        tables[spin] = table
        sources[str(spin)] = {"path": str(path.resolve()), "rows": declared,
                              "sha256": hashlib.sha256(raw).hexdigest()}
    rows = []
    for alpha, spin in ALPHAS:
        table = tables[spin]
        if not table[0, 0] <= alpha <= table[-1, 0]:
            raise ValueError("Frequency table interpolation cannot extrapolate")
        ratio = float(np.interp(alpha, table[:, 0], table[:, 1]))
        tab = alpha * ratio
        leading = alpha * (1 - alpha**2 / 8)
        rows.append({"alpha": alpha, "cloud_spin": spin, "Momega_N2": leading,
                     "Momega_table": tab, "delta_Momega_N2_minus_table": leading - tab,
                     "frequency_relative_difference_percent": (leading / tab - 1) * 100,
                     "period_relative_difference_percent": (tab / leading - 1) * 100})
    return {"convention": "dimensionless angular frequency M*omega_R",
            "formula": "alpha*(1-alpha**2/8), principal N=2 (211 leading)",
            "table_method": "linear interpolation of source omega_R/mu, times alpha",
            "relative_difference_denominator": "external interpolated frequency, not binding energy",
            "caveat": "External tables are interpolated source data, not certified exact eigenvalues.",
            "sources": sources, "rows": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radial-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.radial_root)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    (args.output_dir / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output_dir / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result["rows"][0]))
        writer.writeheader()
        writer.writerows(result["rows"])
    for row in result["rows"]:
        print(f"{row['alpha']:.2f} spin={row['cloud_spin']} "
              f"N2={row['Momega_N2']:.12f} table={row['Momega_table']:.12f} "
              f"difference={row['frequency_relative_difference_percent']:+.6f}%")


if __name__ == "__main__":
    main()
