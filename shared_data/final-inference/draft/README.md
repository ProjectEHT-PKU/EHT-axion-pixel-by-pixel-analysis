# Final-inference dataset — development draft

> **Status: mutable development data, not the paper release.**

This directory contains the minimum validated intermediates currently used to
recompute all six 1-D/2-D expected or observed limit curves. It intentionally
skips raw RAPTOR runs, GRMHD image processing, EHT parameter scans, and topset
reconstruction.

## Get and verify data

```bash
git lfs pull --include='shared_data/final-inference/draft/**'
python scripts/verify_shared_data.py
```

The verifier checks size, SHA-256, and NumPy readability for all 21 files in
`manifest.json`. A Git LFS pointer without downloaded content will fail.

## Run all uncertainty modes

```bash
for dimension in 1d 2d; do
  for mode in std cor cor_grmhd; do
    python scripts/run_final_limits.py \
      --dimension "$dimension" \
      --mode "$mode" \
      --observation-dir "shared_data/final-inference/draft/$dimension" \
      --fit-parameter-dir "shared_data/final-inference/draft/$dimension" \
      --output-dir "results/shared-$dimension"
  done
done
```

Add `--observed` to use measured day-to-day EVPA differences. Without it, the
command calculates the expected/null limit used by the main six-curve figure.

## Contents

- `2d/*_all_data.npy`: corrected 1572-image Stokes topsets for four EHT days;
- `2d/{day}.npz`: daily EVPA mean/std and mean P/I;
- `2d/condition.npy`: 757-pixel final mask;
- `2d/{A_params,phi0_params}_dict.npy`: 2-D GRMHD signal fit summaries;
- `1d/*_EVPA_1d.npy`: corrected ten-degree integrated EVPA ensembles;
- `1d/{day}_1d.npz`: daily 1-D EVPA mean/std;
- `1d/{A_params,phi0_params}_1d_dict.npy`: 1-D GRMHD signal summaries.

## Updating this draft

Only update it after rebuilding and validating local intermediates:

```bash
python scripts/publish_shared_data.py draft --force-draft
python scripts/verify_shared_data.py
```

Commit changed data and its manifest together. Git LFS stores a complete object
for every changed large array, so do not regenerate arrays merely to change file
metadata. The frozen paper dataset will be published separately as `v1`.
