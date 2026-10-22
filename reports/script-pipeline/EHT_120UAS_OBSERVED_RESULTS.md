# EHT 120 μas: observed conditional credible limits

## Main result and status

The primary **actual day-difference** analysis is complete in
`/home/canoziia/GRMHD_DATA/analysis/eht-120uas-observed-limits-v1`.
All six dimension/mode cases and their matched zero-center NULL controls were run
at eight α values. The numbers below are **95% conditional Bayesian upper credible
bounds on the model coupling c**, not a conversion to GeV⁻¹, a detection, a
frequentist exclusion, or a calibrated complete-systematics constraint.
The controls use the same selection and working covariance with zero day differences;
they are not median expected limits from repeated simulated experiments.

Grid, expanded-prior, finite-MC covariance and native-centering checks are complete.
The aggregate stage also completed EXIT0 after all four source manifests committed.
Copy-back verification and final execution receipt are documented separately in
`EHT_120UAS_COMPLETION.md`; this report states the conditional scientific result.

### Actual observed c95

Values are rounded for readability; the precision of saved floating-point outputs
is not physical precision. `std` uses diagonal D-draw scatter, `cor` the full
D-draw covariance, and `cor_fixed_background` adds a selected-feature fixed scalar
pair-difference background variance to the correlated case.

| α | 1D std | 1D cor | 1D cor + fixed background | 2D std | 2D cor | 2D cor + fixed background |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 189.5 | 192.2 | 176.6 | 48.51 | 162.4 | 127.7 |
| 0.15 | 30.48 | 25.39 | 26.58 | 14.11 | 39.41 | 25.93 |
| 0.20 | 8.018 | 7.845 | 7.598 | 6.987 | 16.08 | 8.356 |
| 0.25 | 2.853 | 4.032 | 3.629 | 2.952 | 5.386 | 3.130 |
| 0.30 | 1.408 | 1.116 | 1.798 | 0.9010 | 1.492 | 1.383 |
| 0.35 | 0.8269 | 0.1952 | 0.7953 | 0.3965 | 0.6873 | 0.7771 |
| 0.40 | 0.5043 | 0.1288 | 0.4245 | 0.2552 | 0.5300 | 0.5244 |
| 0.44 | 0.2734 | 0.1203 | 0.2653 | 0.1810 | 0.3737 | 0.3832 |

### Matched NULL-center c95 (same C, not new data)

| α | 1D std | 1D cor | 1D cor + fixed background | 2D std | 2D cor | 2D cor + fixed background |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 38.01 | 31.16 | 56.40 | 38.34 | 21.86 | 45.90 |
| 0.15 | 6.777 | 5.683 | 10.11 | 6.814 | 3.988 | 8.168 |
| 0.20 | 2.180 | 1.942 | 3.374 | 2.155 | 1.345 | 2.670 |
| 0.25 | 0.9278 | 0.9067 | 1.438 | 0.9093 | 0.5840 | 1.130 |
| 0.30 | 0.4400 | 0.4154 | 0.6337 | 0.4353 | 0.2646 | 0.5150 |
| 0.35 | 0.2296 | 0.1920 | 0.3177 | 0.2333 | 0.1361 | 0.2686 |
| 0.40 | 0.1408 | 0.1127 | 0.1950 | 0.1448 | 0.07903 | 0.1667 |
| 0.44 | 0.1109 | 0.08900 | 0.1585 | 0.1127 | 0.05657 | 0.1314 |

The actual centers generally weaken these bounds substantially. Their α dependence
also changes with covariance treatment; diagonal and correlated limits are not
interchangeable. Canonical NULL products and canonical inference were not modified.

## Goodness of fit is essential

These are saved joint coupling/phase grid minima at **α=0.30**, not marginal-posterior
maxima, not a refit by the summary helper. Both epochs share one phase.
The component counts are vector sizes, **not independent-bin counts or a claim of
calibrated χ² degrees of freedom**.

| Dimension | Mode | Components | Zero-signal χ² | Best joint-grid χ² |
|---|---|---:|---:|---:|
| 1D | std | 46 | 39.294 | 30.647 |
| 1D | cor | 46 | 158.817 | 152.953 |
| 1D | cor_fixed_background | 46 | 25.044 | 19.687 |
| 2D | std | 40 | 37.757 | 34.105 |
| 2D | cor | 40 | 196.013 | 171.431 |
| 2D | cor_fixed_background | 40 | 44.895 | 38.059 |

**The D-only correlated conditional model fits poorly, including at its best template
fit. It does not support a formal new-physics detection or exclusion.** The fixed
background cases have much less strained quadratic residuals, but that benchmark
noise prescription is not independently calibrated by an acceptable χ² value.
No likelihood-ratio significance or coverage claim is made.

The background pair σ is **0.0375864 rad = 2.1535° (1D)** and
**0.0367089 rad = 2.1033° (2D)**. It is pooled on the final selected features,
not taken from an unmasked summary, and has no extra √2 correction. It is a fixed
scalar approximation, not the canonical `cor_grmhd` signal-response uncertainty.
The observed engine excludes coupling-squared response/template-fit uncertainty.
The selected background benchmark has **26.679 h** lag versus the assumed **24 h**
day-pair separation; this mismatch is explicitly accepted as a benchmark approximation,
not silently treated as a matched physical background prediction.

## Numerical and finite-MC checks

Sources under `/home/canoziia/GRMHD_DATA/analysis/`:

- `eht-120uas-observed-limits-gridcheck-v1`: phase count 128 and coupling points 240.
  Maximum change from primary is **0.004435 dex**, approximately **1.03%** in c95.
- `eht-120uas-observed-limits-priorcheck-v1`: phase count 128, coupling points 260,
  upper prior bound increased from **1000 to 10000**. Its maximum change from primary
  is 0.004211 dex. This is **not** evidence for a 0.004 dex prior effect: the primary
  grid also differs. Compare directly against the dense grid (the helper saves
  `prior_vs_dense_grid.csv`). The coupling-point count and log spacing still differ,
  so residual differences cannot be attributed exclusively to the prior cutoff.
  The available comparison is consistent with little cutoff dependence at the
  demonstrated grid accuracy, not an exact proof of prior independence. The direct
  prior-versus-dense-grid maximum difference is **0.00031665 dex (~0.073%)**.
- `eht-120uas-observed-covcheck-v1`: with centers, selection and background σ fixed,
  replace the full 1000-draw covariance by first/last 500-draw estimates. Reported
  maximum changes are **0.1061 dex (27.7%) for correlated 1D**, **0.0358 dex (8.6%)
  for fixed-background 1D**, **0.0168 dex (3.94%) for correlated 2D**, and
  **0.0408 dex (9.85%) for fixed-background 2D**; diagonal cases reach 0.0155 dex.
  This is finite-MC covariance sensitivity, **not a standard error of the mean**,
  coverage calibration or a complete covariance-estimation uncertainty treatment.

The primary 1D correlated pair condition numbers are **6.62×10⁷ and 9.63×10⁵**;
2D values are about **9825 and 9639**. D-only scatter describes a geometrically
correlated calibration ensemble, not full measurement noise. Small eigen-directions
strongly weight some residual shapes. Neither this report nor its helper clips
eigenvalues, adds a default diagonal ridge, or changes the likelihood to hide the
poor fit. Finite-MC sensitivity is more consequential than the ~1% grid changes.

Native centering can change the final selection N as well as observations and C.
When complete, its comparison must retain both feature counts and must **not** be
called a pure attribution to centering at fixed data vector.

## Observation and branch diagnostics

The primary common mask retains **N=23 1D features and N=20 2D features**. Pair
vectors have 46 and 40 total components. The likelihood uses matched theory CSR
aperture quadrature on transposed per-draw Stokes, not the diagnostic native
pixel-sector EVPA curves. No daily EVPA alignment, offset removal or extra roll is
applied. Independent days are a working assumption: C_pair=C_later+C_earlier;
MC draw indices across days are not treated as paired realizations.

From the primary observed correlated diagnostics:

| Dimension | Pair | Max absolute day difference | Max marginal pair σ |
|---|---|---:|---:|
| 1D | 96−95 | 18.855° | 18.593° |
| 1D | 101−100 | 17.207° | 14.818° |
| 2D | 96−95 | 16.409° | 9.353° |
| 2D | 101−100 | 11.340° | 12.745° |

No retained observed center is within three marginal σ of its principal-branch seam
for these cases. **This is a local observed-center check only.** The engine uses
unwrapped residual d−c·template; predicted means at large c, low α or some phases
can approach/cross branch seams. A c95 below the upper bound, small last-bin mass,
or a moderate best-fit phase does not prove global validity of a nonwrapped
Gaussian likelihood throughout the prior. A full wrapped-model check is not supplied.

### Context from the completed MC diagnostics

See [MC_DIAGNOSTICS.md](../eht-120uas-observed/MC_DIAGNOSTICS.md) for the full
independent diagnostic report (not a substitute for the selected likelihood data).
Key facts retained here:

- Four days each have **1000 draws (4000 total)**, fixed I and conditional D-term
  calibration scatter; no thermal-noise bootstrap or full imaging uncertainty.
- Early SW diagnostic-sector differences are small; late 101−100 differences in
  parts of the same bright SW region are approximately **5–9°**, with larger
  conditional calibration scatter. Overlapping sectors are not independent tests.
- Day 100 has net polarization **0.995±0.479%** and integrated EVPA local std
  **26.55°**, while β₂ phase std is only **3.61°**. Vector cancellation can broaden
  the global EVPA without erasing spatial polarization information.
- Data frequency is **227.070703125 GHz**, not exactly the nominal theory 230 GHz.
  `display19` means an additional 19 μas kernel product, not a certified final
  19 μas instrument resolution. Small centroid shifts do not establish zero
  registration uncertainty.

## Final aggregate helper (executed successfully)

Repository helper: `eht-script-pipeline/scripts/summarize_eht_ensemble_results.py`.
It reads completed small outputs only, verifies their declared output hashes with
`verify_stage(recursive=False)`, then writes a new stage with byte-identical receipt
snapshots named `*-receipt.json`. These are **receipts**, not rebound canonical
manifests or a fresh raw-lineage calibration certificate. The source manifest SHA
and original path remain explicit. No 97 GB raw replay or static RAPTOR/EHT source
closure hashing is requested.

Run from `eht-script-pipeline` with its Python environment:

```sh
ROOT=/home/canoziia/GRMHD_DATA/analysis
.pixi/envs/default/bin/python scripts/summarize_eht_ensemble_results.py \
  --primary "$ROOT/eht-120uas-observed-limits-v1" \
  --gridcheck "$ROOT/eht-120uas-observed-limits-gridcheck-v1" \
  --priorcheck "$ROOT/eht-120uas-observed-limits-priorcheck-v1" \
  --native "$ROOT/eht-120uas-observed-limits-native-v1" \
  --covcheck "$ROOT/eht-120uas-observed-covcheck-v1" \
  --output-root "$ROOT/eht-120uas-results-summary-v1"
```

The native directory and all other sources were verified before execution. The helper
fails fast if any source lacks a complete manifest and does not poll or wait.
Outputs: `limits.csv`, comparison CSVs, `alpha030-chi2.csv`, `report.json`,
`RESULTS.md`, **`limits-observed-null.png`** (six observed solid curves with matched
NULL dashed controls in 1D/2D panels; log c ordinate), and **`alpha030-fits.png`**
(three mode rows, two dimension columns). The latter reads saved shared-phase
joint-grid predictions, plots data ± effective marginal pair scatter and reports
Cholesky-whitened residual RMS; it performs no refit. Feature order in the 2D plot
is not spatial distance. Existing MC plots and metadata are not overwritten.

Small tests: `PYTHONPATH=scripts .pixi/envs/default/bin/python -m pytest -q tests/test_eht_ensemble_results.py`.
