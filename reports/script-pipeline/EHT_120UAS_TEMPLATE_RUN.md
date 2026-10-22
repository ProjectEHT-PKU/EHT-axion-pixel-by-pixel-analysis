# EHT120 template preparation (not science certification)

## Scope and validated preflight
Only new files and a new output root; no old source/model/limits edits, RT,
EHT processing, ensemble duplication, or limits calculation.

Profile: `configs/eht-120uas-observed.yaml`. Launcher:
`scripts/start_eht_120_templates.py`.
Raw input is exactly `/home/canoziia/GRMHD_DATA/production/m65-211-merged/data`.
Output ROOT: `/home/canoziia/GRMHD_DATA/analysis/eht-120uas-observed-templates-v1`.
Old observed template root was inspected for layout only; no old fits reused.

Preflight before launch:
- `pixi run python scripts/start_eht_120_templates.py --validate-only`: passed
  strict config and unchanged scientific-section assertions.
- `pixi run python scripts/process_raptor.py --help`: CLI accepted, stage enum
  checked against `pipeline/workflow.py`: `raptor` executes preprocess/fit/all-alpha
  1d+2d summaries, no EHT stage.
- `pixi run test tests/test_pipeline_raptor.py`: **50 passed**.

Geometry: 120 uas / 64 = 1.875 uas/pixel, FOV31.42189369570729Rg,
aperture radius60uas =15.710946847853645Rg (circle, not square corners).
Initial46.2261Rg, rotation108deg, mass6.5e9, distance16.8Mpc,
referenceNorm0.1, principalN2, cloud coefficient0.03 and all eight alphas retained.
Gaussian FWHM19uas => sigma2.1127445296660023Rg: **additional matching-kernel
approximation** to EHT already-19uas display, not an effective-PSF certification.
Cadence must come from raw provenance: expected8.893247544352288h, lag3 =>
26.679742633056865h, not exact24h; no sqrt2 correction.

EHT contract declares120uas/64/1.875, with all four rank/final offsets zero.
These are not old pixel offsets. angle_roll=-90 explicitly retains existing
theory aperture convention; parent owns native-coordinate/centroid adapter.
No trusted EHT pickle or geometry evidence is fabricated or consumed.

## Execution
Single `verification_session()` wraps canonical `_stage.execute('raptor', ...)`
and background preparation, so shared raw/HDF lineage content hashes are memoized
within this invocation with identity checks, not repeatedly read per alpha.
No persistent hash bypass, imported stage trust, or old output reuse.

Launch command (tmux `eht120-templates`, log `/tmp/eht120-templates/full.log`):
```
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  pixi run python scripts/start_eht_120_templates.py
```
Sequential frame-streaming avoids duplicate workers/large active copies. No RT.

## Parent consumer paths and explicit boundary
- `ROOT/raptor/templates/1d/templates.npz`
- `ROOT/raptor/templates/2d/templates.npz`
  (`response_summary_roots`: corresponding containing directories)
- `ROOT/raptor/response/alpha-0.XX/response-{1d,2d}.npz` (all eight alpha labels)
- `ROOT/raptor/background/stokes.npy` (canonical `RunLayout.background_processed`)
- `ROOT/background-unmasked/{1d,2d}/{summary.json,evpa-difference.npz}`
- `ROOT/template-run/{manifest.json,report.json}` completion/consumer record.

**Background boundary:** canonical workflow `background` requires verified
`ROOT/eht/manifest.json` plus `ROOT/eht/mask.npy`. Neither exists until parent
completes the new EHT adapter. Launcher therefore prepares canonical res400
state0 Stokes and explicitly **unmasked** 1d/2d lag-difference diagnostics via
existing `estimate_background`, without inventing an EHT mask. Those diagnostics
are NOT selected-feature Cbg. Parent must construct selected-feature background
statistics/Cbg from the prepared Stokes using the new ensemble selection, then
run observed std/cor/cor_fixed_background. This worker does not certify or alter
that API. Existing research observed CLI is hard-coded to128uas/64uas radius;
parent must use its matching120 adapter, not that CLI unchanged.

Status at creation: preflight passed; launch/completion evidence appended below.
