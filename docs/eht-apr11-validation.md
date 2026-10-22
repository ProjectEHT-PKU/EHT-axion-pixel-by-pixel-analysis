# Apr11 objective / initialization validation

This is a separate optional working-model validation, NOT four-day imaging or
an axion inference run. Historical executed scripts and products remain intact.

## Single entrypoint

```bash
pixi run -e eht-rebuild python scripts/validate_eht_apr11.py \
  --config configs/eht-apr11-validation.yaml --case baseline \
  --rounds 100 --output-root /NEW/apr11-baseline --plan > /tmp/apr11-plan.json
```

`--plan` validates source identities and writes nothing to the output directory.
It prints `audit_template` with exact source hashes and selected objective hash.
The responsible reviewer must complete the objective audit and set
`complete: true, working_model_approved: true` in a separate JSON evidence file.
This approval means permission for this explicit working-model comparison,
NOT certification that the objective equals the historical PaperVII pipeline.

```bash
pixi run -e eht-rebuild python scripts/validate_eht_apr11.py \
  --config configs/eht-apr11-validation.yaml --case baseline --rounds 100 \
  --audit-file /PATH/to/apr11-objective-audit.json \
  --output-root /NEW/apr11-baseline
```

Run `evpa45` and `restart` with different new output roots and the same audit if
all objective/source bytes are unchanged. Each case uses exactly100rounds by
default, not an undocumented convergence cutoff. A shorter explicit budget is
recorded as partial. No case is selected because it resembles a published image.

- baseline: initial5%fraction,EVPA0°, no blur restart;
- evpa45: onlyinitialEVPA45° changes;
- restart: baselineinitialization,20uasblurrestarts2/12/22/32/42.

The noise, fixed I, normalization, objective weights, station map and D prior
are identical between cases. `--norm-reg/--no-norm-reg` changes the objective
and requires a NEW audit; it is not a harmless implementation switch.

## Independent diagnostics

```bash
pixi run -e eht-rebuild python scripts/diagnose_eht_apr11.py \
  --runs /NEW/apr11-baseline /NEW/apr11-evpa45 /NEW/apr11-restart \
  --base-input /home/canoziia/GRMHD_DATA/analysis/eht-official-I-cleanpol-apr11-v132-v2/official-I-gaincal-clean-pol.uvfits \
  --feed-file /home/canoziia/GRMHD_DATA/EHT_DATA/official-2023-D01-01/reference/feed-rotation.json \
  --output-root /NEW/apr11-diagnostics --norm-reg-audit-pending
```

This performs no imaging or new calibration fit. It validates declared stages,
reapplies final station-indexed D values once to the same base, samples official
direct model visibilities and independently checks residual reductions against
upstream statistics. Reports group by baseline, station and time; overlapping
station/conjugate groups are NOT independent observations or telescope fault
attribution. All uncertainties are the stated diagonal approximation.

Map comparisons share one fixed intersection mask across the same19uasdisplay
kernel. No hidden registration orEVPArotation is fitted. Latecheckpoints80/90
are compared with100, and finalcases withbaseline. The19uas kernel is not a
measuredPSF; comparisonpixelcount is not a numberofindependentresolutionelements.

Completed audit/result report:
`reports/script-pipeline/APR11_CONTROLLED_VALIDATION.md`.

Known mapping gaps: PaperVII HW/PTV numerical normalization, exact systematic
error component convention and D prior were not recovered. Current legacy
D API uses inverse-correcteddata / D-dependent propagated diagonal sigmas,
whereas PaperVII describes forward-Jones model fitting. Keeping those assumptions
fixed lets us test initialization/restart sensitivity; it does not certify them.
No trusted observedaxionlimit should be produced from thisvalidation alone.
