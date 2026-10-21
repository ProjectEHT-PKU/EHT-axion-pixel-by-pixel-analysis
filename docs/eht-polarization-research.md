# EHT polarimetric image / observed-data research

This optional workflow does **not** replace the NULL-only canonical analysis.
It processes the backed-up original reconstructed FITS; it does not reconstruct
UVFITS or claim to reproduce the official five-pipeline result.

```bash
pixi install -e eht-research
pixi run -e eht-research research-polarization \
  --config configs/eht-polarization-study.yaml \
  --output-root /home/canoziia/GRMHD_DATA/analysis/eht-polarization-study-v1 \
  --trust-pickle --study-id full4800-source-fits
```

It validates candidate identities/file completeness, reads each native IQUV FITS,
adds a19microarcsecond FWHM kernel, and compares unregistered versus I-only
registered maps on one fixed local Apr11 best-I reference. It preserves celestial
WCS and Q/U signs, translating all Stokes together. It saves a float64 mmap and
per-image coherent10degree aperture EVPA, using the same quadrature as the
RAPTOR template fits. Absolute plotting ticks use10uas=30microJy/uas² polarized
brightness and colors represent fractional polarization;120uas is displaycrop,
128uas remains the native data grid. Display mask requires both I>=10%Ipeak and
P>=20%Ppeak; papercaption's AND wording is ambiguous, this is a declared stricter
research display choice, not proof of originalplottingmask. Inference selections
are separately declared and must not be confused with rendering masks.

Variants:
- all4800: full local hyperparameter ensemble;
- legacy_evpa: exact existing1572 selected identities;
- datafit_perday: all five stored metrics<2 independently perday;
- datafit_shared: same identity must pass all fourdays;
- bestI: one locally best mean-I fit image/day, no covariance estimation.

None of these is an official polarimetric Top Set. PaperIV's1572 was the result
of a37500-combination StokesI survey with datafit AND syntheticrecovery gates;
PaperVII chose polarizedfiducials by a six-synthetic-source image/Dterm score.
Its1000-sample uncertainty ensembles came from Dtermperturbations (or posterior
samples), not our4800hyperparametergrid. Stored localχ² describes optimizer
stages before final5uas restoration, not a fresh finalFITS/UVFITS consistency test.
The original restoration remains;19uas additionalkernel is NOT a certifiedPSF.

## Observed research

Use the separate theorypostprocessing profile, leaving oldoutputs intact:

```bash
pixi run raptor --config configs/eht-observed-research.yaml \
  --raptor-input-root /home/canoziia/GRMHD_DATA/production/m65-211-merged/data \
  --output-root /home/canoziia/GRMHD_DATA/analysis/eht-observed-research-templates-v1
pixi run -e eht-research research-observed \
  --config configs/eht-observed-research.yaml \
  --study /home/canoziia/GRMHD_DATA/analysis/eht-polarization-study-v1 \
  --templates /home/canoziia/GRMHD_DATA/analysis/eht-observed-research-templates-v1 \
  --variant datafit_shared \
  --output-root /home/canoziia/GRMHD_DATA/analysis/eht-observed-datafit-shared-v1
```

Observed entrypoint verifies sourcebuild/physical/templategeometry anddeclared
files. It constructs actual later-minus-earlier data forApr6-Apr5 andApr11-Apr10,
thenjointconditionalGaussianlikelihoodwithonecommonphase. It also computes an
exactly matched zero-center comparison. Modesstd/cor only: no coupling-dependent
GRMHDresponsevariance is smuggled into this study. OriginalNULL APIs unchanged.
Populationcovariance of localbranch hyperparameterangles is used only as a
workingcovariance, withzero crossday blocks assumed. Resultant>=0.7,1D10degree
spacing offset4,2D5pixelstride and displaybrightnessintersection are declared
research selections. No automaticridge/eigenvalueclipping. These areconditional
crediblelimits, NOT calibrated observationalexclusions ormedianexpectedlimits.

If branch fails, covariance is singular, residualfit is poor or bounds hit prior
edges, report the failure rather than silently filteringout signal or modifying
uncertainty until a desired limit emerges. Do not interpret daydifference as
axiondetections without sourceevolution/calibration/coverage validation.

Literature/sourceevidence: PHYS/reports/eht-polarization-research/
`OFFICIAL_METHOD.md`, `LOCAL_PROVENANCE.md`. Longrunsusetmux; no defaultanalysis
entrypoint launches this optionalstudy oranyradiative/reconstructionproduction.
