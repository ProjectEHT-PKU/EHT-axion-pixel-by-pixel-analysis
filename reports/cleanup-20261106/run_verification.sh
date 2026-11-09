#!/usr/bin/env bash
# End-to-end rerun of the NEW code on the migrated job root (no ray tracing, no imaging),
# then numerical comparison with the original products. Run after migration/migrate_jobs.py.
# Usage: run_verification.sh [first_step]   (steps run in order; <= 8 cores via taskset)
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
# Process 1 -> delivered Stokes cubes from the migrated jobs
step merge $T $PY scripts/raptor_merge.py --plan "$N/plan.json" --job-root "$N/m65-211" \
  --job-root "$N/pilot-m65-t3000-res400-t16" --output-root "$V/raptor-merged"
# Process 3.1 templates and background
step templates $T $PY scripts/analysis_templates.py --raptor-input-root "$V/raptor-merged/data" \
  --output-root "$V/templates"
# Process 2 -> 3 boundary: EHT ensemble summary of the 4000 recorded draws
step ensemble $T $PYE scripts/eht_summarize.py --dataset-root "$A/eht-uncertainty-remote-v1" \
  --output-root "$V/eht-ensemble"
# Process 3.2 limits (native = paper; centroid = comparison)
step limits-native $T $PY scripts/analysis_limits.py --templates "$V/templates" --ensemble "$V/eht-ensemble" \
  --output-root "$V/limits-native"
step limits-centroid $T $PY scripts/analysis_limits.py --templates "$V/templates" --ensemble "$V/eht-ensemble" \
  --centering centroid --output-root "$V/limits-centroid"
# Process 3.3/3.4 checks, for both centerings
for c in native centroid; do
  step gridcheck-$c $T $PY scripts/analysis_check_grid.py --source "$V/limits-$c" --templates "$V/templates" \
    --output-root "$V/$c-gridcheck" --phase-count 128 --coupling-points 240
  step priorcheck-$c $T $PY scripts/analysis_check_grid.py --source "$V/limits-$c" --templates "$V/templates" \
    --output-root "$V/$c-priorcheck" --phase-count 128 --coupling-points 260 --coupling-max 10000
  step covcheck-$c $T $PY scripts/analysis_check_covariance.py --source "$V/limits-$c" --templates "$V/templates" \
    --output-root "$V/$c-covcheck"
done
# Process 3.5 summary (same source roles as eht-120uas-results-summary-v1)
step summary $T $PY scripts/analysis_summarize.py --primary "$V/limits-centroid" --gridcheck "$V/centroid-gridcheck" \
  --priorcheck "$V/centroid-priorcheck" --native "$V/limits-native" --covcheck "$V/centroid-covcheck" \
  --output-root "$V/results-summary"

# Comparisons with the original products
C="$T $PY reports/cleanup-20261106/compare_products.py"
step compare bash -c "
cd reports/cleanup-20261106
$C $OLD_MERGED/data $V/raptor-merged/data merged-raptor
$C $A/eht-120uas-observed-templates-v1 $V/templates templates --skip invocations config.yaml manifest.json template-run
$C $A/eht-mc-summary-v1 $V/eht-ensemble eht-ensemble
$C $A/eht-120uas-observed-limits-native-v1 $V/limits-native limits-native
$C $A/eht-120uas-observed-limits-v1 $V/limits-centroid limits-centroid
$C $A/eht-120uas-native-gridcheck-v1 $V/native-gridcheck native-gridcheck
$C $A/eht-120uas-native-priorcheck-v1 $V/native-priorcheck native-priorcheck
$C $A/eht-120uas-native-covcheck-v1 $V/native-covcheck native-covcheck
$C $A/eht-120uas-observed-limits-gridcheck-v1 $V/centroid-gridcheck centroid-gridcheck
$C $A/eht-120uas-observed-limits-priorcheck-v1 $V/centroid-priorcheck centroid-priorcheck
$C $A/eht-120uas-observed-covcheck-v1 $V/centroid-covcheck centroid-covcheck
$C $A/eht-120uas-results-summary-v1 $V/results-summary results-summary
"
echo ALL_DONE >> "$V/logs/timings.txt"
