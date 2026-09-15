#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
export PYTHONPATH="$ROOT"
export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-4}

python -m compileall -q eht_axion scripts tests utils
python -m pytest -q
ruff check eht_axion scripts tests
python scripts/verify_legacy_results.py
python scripts/validate_historical_selection.py

for dimension in 1d 2d; do
  observation_dir=audit/rebuilt_2d
  fit_args=()
  if [[ $dimension == 1d ]]; then
    observation_dir=audit/rebuilt_1d
    fit_args=(--fit-parameter-dir audit/rebuilt_1d)
  fi
  for mode in std cor cor_grmhd; do
    python scripts/run_final_limits.py \
      --dimension "$dimension" \
      --mode "$mode" \
      --observation-dir "$observation_dir" \
      "${fit_args[@]}" \
      --output-dir audit/results/rebuilt_corrected_expected
  done
done

python scripts/plot_final_limits.py

python - <<'PY'
import ast
import json
from pathlib import Path

for path in map(Path, ("result.ipynb", "GRMHD/main.ipynb", "EHT-Imaging/main.ipynb")):
    notebook = json.loads(path.read_text())
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if any(line.lstrip().startswith(("%", "!")) for line in source.splitlines()):
            continue
        ast.parse(source, filename=f"{path}:cell-{index}")
print("Notebook code-cell syntax: OK")
PY

echo "All verification gates passed."
