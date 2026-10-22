# Outer calibration endpoint stability: Apr10 and follow-up queue

## Completed Apr10 comparison

Source: `/home/canoziia/GRMHD_DATA/analysis/eht-convex-tightD-100-100round-v1/`.
This starts from the original50-round `eht-paper-PD-100-v1` and performs100
ADDITIONAL rounds with certified convex P solves and tighter numerical D options.
It is not a homogeneous150-round run of a single algorithm. Every saved round-N
P was solved for D(N-1), then D(N) was updated and saved. The final-P call is
aligned to D100. Comparisons below are saved iteration endpoint drift, not
independently refitted matched-D experiments or uncertainty samples.

`check_eht_outer_stability.py` verifies declared source closure, loads native
IQUV FITS, checks identical grid/fixedI, and compares saved additional50/80/100
and final images. No optimization, map translation, EVPA offset fit or image
rescaling. Common19uas IQU additional kernel is NOT a measured effective PSF.

Same fixed553-pixel intersection mask across all four displayed maps:
I>=.1peak AND P>=.2peak. This is a diagnostic support, not a claimed official
mask. Pixel EVPA RMS uses symmetric fractional-P weights and modulo-pi angles.
Whole-image integrated EVPA and netm use coherent sums; not mean pixel angles.

|Saved endpoints|bright EVPA RMS deg|display QU relativeL2|integrated EVPA difference deg|
|---|---:|---:|---:|
|additional50 to additional80|.24730|.0122373|-.84989|
|additional80 to additional100|.05434|.0026754|-.19281|
|additional100 to final-P|.001620|.00007887|-.006068|
|additional50 to final-P|.30316|.0150031|-1.04877|
|additional80 to final-P|.05594|.0027538|-.19888|

Whole-map netm:1.14241% at additional50,1.15898% at80,1.16234% at100,
1.16244% atfinal. Resolvedm:5.55874%,5.58563%,5.59139%,5.59155% respectively.
Diagnostic beta2 centered on the fixedI centroid: abs .0305650->.0306162;
phase-156.4509->-156.5694deg (50->final). These are not certified official beta2
centering conventions.10deg all-radius intensity-weighted apertures at1deg steps
also saved; overlapping bins are not independent samples.

Interpretation: strong decline in late endpoint drift; local and integrated
EVPA stability differ because coherent sums can cancel. Small iterated changes
are not calibration/measurement uncertainties, official acceptance thresholds,
or evidence that poor Apr10 uv coverage is cured.

Outputs:
`/home/canoziia/GRMHD_DATA/analysis/eht-outer-stability-100-v1/`
report.json, apertures.csv, common-mask.npy, stability.png, manifest.json.
Newpure-helper tests4passed (identity, modulo-pi seam, coherent cancellation,
whole-image-vs-mask distinction); new sources Ruff pass.

## Next three days: actual execution, not completed results

Already launched sequential101,95,96, each100additionalrounds+finalP, using
unchanged same solver/gap target and D numerical settings as Apr10. OriginalI
and50-round products reused. Every day automatically runs the same endpoint
comparison after its numerical stage completes. No full4000MC or newlimits.

Newstage roots:
`.../eht-convex-tightD-<day>-100round-v1/`
`.../eht-outer-stability-<day>-v1/`.

Queue tmux:eht-outer-followup; read-only60secondmonitor:eht-outer-monitor.
Logs/status under `/tmp/eht-outer-followup/`.
Independent monitoring subagent waits with no token/turn budget or childGoal.
It checks P certificates, exceptions, exits and published comparisons, and will
report completion or failure without auto-restarting or changing science.

Budget perday is deliberately explicit, not a conclusion that100 is sufficient.
Use actual map/whole-polarization drift and D trajectory for each day; if still
changing, report remaining drift rather than relabel a budget as convergence.
