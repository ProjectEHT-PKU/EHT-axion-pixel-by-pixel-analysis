#!/usr/bin/env bash
# End-to-end rerun of the NEW code on the migrated job root (no ray tracing, no imaging),
# then numerical comparison with the original products. Run after migration/migrate_jobs.py.
# Usage: run_verification.sh [first_step]   (steps run in the listed order; <= 8 cores via taskset)
# Record of the 2026-11-09 run: merge/templates/ensemble/limits-native ran with commits
# 64ca9a1..afb7118, checks/comparisons with afb7118; the summary step with the final code.
# The centroid variant was dropped (Paper VII only), see CLEANUP_ZH.md.
set -euo pipefail
W=/home/canoziia/MyProject/PHYS/eht-script-pipeline-cleanup
V=/home/canoziia/GRMHD_DATA/analysis/cleanup-verify-opus-20261106
N=/home/canoziia/GRMHD_DATA/production/m65-211-newids-20261106
A=/home/canoziia/GRMHD_DATA/analysis
OLD_MERGED=/home/canoziia/GRMHD_DATA/production/m65-211-merged
PY=/home/canoziia/MyProject/PHYS/eht-script-pipeline/.pixi/envs/default/bin/python
PYE=/home/canoziia/MyProject/PHYS/eht-script-pipeline/.pixi/envs/eht-rebuild/bin/python
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
T="taskset -c 24-31"
mkdir -p "$V/logs"
FIRST=${1:-merge}
started=0

step() {
  local name=$1; shift
  if [ "$started" = 0 ] && [ "$name" != "$FIRST" ]; then return; fi
  started=1
  local t0; t0=$(date +%s)
  echo "START $name $(date -Is)" >> "$V/logs/timings.txt"
  "$@" > "$V/logs/$name.log" 2>&1
  echo "DONE  $name $(( $(date +%s) - t0 ))s" >> "$V/logs/timings.txt"
}

cd "$W"
# RAPTOR: delivered Stokes cubes from the migrated jobs
step merge $T $PY scripts/raptor_merge.py --plan "$N/plan.json" --job-root "$N/m65-211" \
  --job-root "$N/pilot-m65-t3000-res400-t16" --output-root "$V/raptor-merged"
# Analysis: templates and background
step templates $T $PY scripts/analysis_templates.py --raptor-input-root "$V/raptor-merged/data" \
  --output-root "$V/templates"
# EHT imaging: ensemble summary of the 4000 recorded draws
step ensemble $T $PYE scripts/eht_summarize.py --dataset-root "$A/eht-uncertainty-remote-v1" \
  --output-root "$V/eht-ensemble"
# Analysis: limits (Paper VII: images used as reconstructed, zero shift)
step limits-native $T $PY scripts/analysis_limits.py --templates "$V/templates" --ensemble "$V/eht-ensemble" \
  --output-root "$V/limits-native"
# Analysis: checks
step gridcheck-native $T $PY scripts/analysis_check_grid.py --source "$V/limits-native" --templates "$V/templates" \
  --output-root "$V/native-gridcheck" --phase-count 128 --coupling-points 240
step priorcheck-native $T $PY scripts/analysis_check_grid.py --source "$V/limits-native" --templates "$V/templates" \
  --output-root "$V/native-priorcheck" --phase-count 128 --coupling-points 260 --coupling-max 10000
step covcheck-native $T $PY scripts/analysis_check_covariance.py --source "$V/limits-native" --templates "$V/templates" \
  --output-root "$V/native-covcheck"
# Analysis: summary (native limits are the primary result)
step summary $T $PY scripts/analysis_summarize.py --limits "$V/limits-native" --gridcheck "$V/native-gridcheck" \
  --priorcheck "$V/native-priorcheck" --covcheck "$V/native-covcheck" --output-root "$V/results-summary"

# Comparisons with the original products
C="$T $PY reports/cleanup-20261106/compare_products.py"
step compare bash -c "
cd reports/cleanup-20261106
$C $OLD_MERGED/data $V/raptor-merged/data merged-raptor
$C $A/eht-120uas-observed-templates-v1 $V/templates templates --skip invocations config.yaml manifest.json template-run
$C $A/eht-mc-summary-v1 $V/eht-ensemble eht-ensemble
$C $A/eht-120uas-observed-limits-native-v1 $V/limits-native limits-native
$C $A/eht-120uas-native-gridcheck-v1 $V/native-gridcheck native-gridcheck
$C $A/eht-120uas-native-priorcheck-v1 $V/native-priorcheck native-priorcheck
$C $A/eht-120uas-native-covcheck-v1 $V/native-covcheck native-covcheck
"
echo ALL_DONE >> "$V/logs/timings.txt"
