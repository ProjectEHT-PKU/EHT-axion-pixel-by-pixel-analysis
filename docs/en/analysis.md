# Analysis

## Purpose

This page describes the analysis, which combines the two inputs: it turns the delivered
RAPTOR Stokes cubes into EVPA response templates, selects likelihood features on the EHT
ensembles, builds the day-pair EVPA differences and their covariance, and computes
phase-marginalized posterior upper limits $c_{95}(\alpha)$ on the dimensionless coupling
$c$, together with a zero-data control and numerical checks.

## Inputs

- The delivered RAPTOR data root (`data/` of a delivery or merge), containing
  `raptor/res-100/alpha-<α>/` and `raptor/res-400/background/`
  ([RAPTOR simulation](raptor-simulation.md)).
- The EHT ensemble stage of `scripts/eht/summarize.py` ([EHT imaging](eht-imaging.md)).
- `configs/experiment.yaml`: `shared`, `analysis`, and from the other sections
  `raptor.camera_fov_rg`, `raptor.snapshots`, `raptor.response`, `raptor.background` and
  `eht_imaging.uncertainty_draws` ([Experiment configuration](experiment-config.md)).

## Outputs

- Templates root: preprocessed cubes, per-snapshot response fits, templates (1D and 2D),
  preprocessed axion-off frames.
- Limits root: selected observations, GRMHD background, limits for 2 dimensions × 3
  covariance models × {observed, zero-data control}, `limits.csv`, a figure and
  `report.json`.
- Check roots (grid, prior, covariance), a summary root and the paper root.

File details are on [Data layout](data-layout.md).

## How to run

Default environment. With a data root (`M87_AXION_DATA_ROOT` or `--data-root`), every root option below has a
default from the [default layout](data-layout.md#default-layout); the commands give the roots
explicitly.

```bash
pixi run analysis-templates --raptor-input-root MERGED/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root GRID \
    --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root PRIOR \
    --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance --source LIMITS --templates TEMPLATES --output-root COV
pixi run analysis-summarize --limits LIMITS --gridcheck GRID --priorcheck PRIOR \
    --covcheck COV --output-root SUMMARY
pixi run analysis-paper --limits LIMITS --gridcheck GRID --priorcheck PRIOR --covcheck COV \
    --templates TEMPLATES --ensemble ENSEMBLE --raptor-input-root MERGED/data \
    --raptor-plan plan.json --raptor-run-root RUN_ROOT --output-root PAPER
```

The override values above are examples; the check scripts accept any values that differ
from `analysis.inference`. `--covcheck` is optional for `scripts/analysis/summarize.py`.

## Details

### Image preprocessing

`scripts/analysis/templates.py` verifies each delivered cube against the current experiment
(`generated_inputs`: identity record, shape, labels equal to `raptor.snapshots`) and runs
`preprocess_cube` (`eht_axion/analysis/templates.py`) frame by frame. For each Stokes
channel, `preprocess_frame`:

- convolves with a Gaussian of standard deviation
  $\sigma = \mathrm{FWHM}/\sqrt{8\ln 2}$, with FWHM = `shared.image.blur_fwhm_uas` (19 µas)
  converted to $r_g$ (2.11 $r_g$), using nearest-value boundary extension;
- maps the camera field `raptor.camera_fov_rg` onto the sky grid of `shared.image.shape`
  pixels and `shared.image.fov_uas` (31.42 $r_g$) by bilinear interpolation, rotated by
  `analysis.rotation_degrees` $\theta$ = 108°;
- rotates the polarization: $Q + iU \to e^{2i\theta}(Q + iU)$.

No flux renormalization is applied. Axis 0 of the simulated images is identified with the
west direction (decreasing east offset) and axis 1 with north, matching the transposed
EHT maps (see [Feature selection](#feature-selection)). Response cubes keep all 17 states;
the background cube keeps state 0 only (`analysis:background`).

### Response fit

`fit_processed` drops state 0 and, for each snapshot, fits the 16 EVPA values
$\chi_k$ of the phase states (axion phase $\psi_k = 2\pi k/16$) per pixel (2D) and per
aperture sector (1D, next section) with `fit_phase_dft`:

$$
\chi_k = \bar\chi + A\cos(\psi_k + \phi).
$$

The values are first lifted onto the branch of their circular mean (period $\pi$), then
$A e^{i\phi} = \tfrac{1}{8}\sum_k (\chi_k - \bar\chi)\, e^{-2\pi i k/16}$. A location with a
missing phase is invalid; zero amplitude has an undefined phase. The output per $\alpha$
is `response-2d.npz` and `response-1d.npz` with arrays `amplitude`, `phase`, `baseline`
(radians) of shape (141, …) and `snapshot_labels`.

### Templates

`summarize_fits` combines the 141 snapshots for each $\alpha$ and location: the mean and
population standard deviation (ddof = 0) of $A$, and the circular mean of $\phi$ with the
standard deviation of $\phi$ about it. A location needs at least two finite snapshots.
`templates.npz` holds `alphas` and arrays of shape (8, 360) for 1D and (8, 64, 64) for 2D.
`load_templates` marks a location valid when all four statistics are finite and the
standard deviations are nonnegative for every $\alpha$. The signal model uses only the
means; the standard deviations enter only this validity mask.

### Aperture observable

The 1D observable is the EVPA of the intensity-weighted polarization in a sector
("pizza-slice") aperture:

$$
\chi_\mathrm{ap}(\varphi_k) = \tfrac12 \arg \int_0^{r_\mathrm{max}}\!\!\int_{\varphi_k-w/2}^{\varphi_k+w/2}
 I\,(Q + iU)\; r\,\mathrm{d}\varphi\,\mathrm{d}r ,
$$

with $w$ = `analysis.aperture.window_degrees` (10°),
$r_\mathrm{max}$ = `radial_max_fraction` × half field (60 µas) and 360 sector centres one
degree apart. `aperture_operator` (`eht_axion/analysis/aperture.py`) is a sparse matrix of
bilinear quadrature weights (radial step 0.02 in field units, angular step 0.2°,
trapezoidal end weights). The same operator is applied to the simulated images and, after
transposing to the same axis order, to every EHT draw, so model and data use an identical
observable. After the operator, the bins are rolled by 90 so that bin $k$ is centred
$k$ degrees east of north.

With the quadrature weights $W_{kj}$ of sector $k$ and pixel $j$, the code computes:

- the sector EVPA $\chi_k = \tfrac12\arg\sum_j W_{kj}\, I_j (Q_j + iU_j)$
  (`integrate_stokes`; EHT Paper VII eq. 11 with this quadrature);
- the sector intensity $B_k = \sum_j W_{kj}\, I_j$;
- the intensity-weighted mean polarized brightness
  $\langle\mathcal{P}\rangle_k = \sum_j W_{kj}\, I_j \lvert P_j\rvert / \sum_j W_{kj}\, I_j$
  with $\lvert P\rvert = \sqrt{Q^2+U^2}$ (`sector_brightness`; Paper VII eq. 10 with this
  quadrature).

$B_k$ and $\langle\mathcal{P}\rangle_k$ use the same weights and the same intensity
weighting as the observable; they are used only to select sectors.

### Signal model

For one day pair $p$ the model EVPA difference at feature $x$ is linear in $c$:

$$
s_p(x;\alpha,\phi_0) = b(\alpha)\,\frac{2\sin(\omega\tau/2)\,\mathrm{sinc}(\omega T/2)}{N_\mathrm{ref}}\,
 A(x,\alpha)\,\sin\!\big(\phi_0 + \phi(x,\alpha) + \omega t_p\big),
$$

where

- $\omega = \Omega(\alpha)/t_g$ in rad h$^{-1}$, $\Omega(\alpha) = \alpha(1-\alpha^2/(2n^2))$;
- $\tau$ = `day_gap_hours` (24 h), $T$ = `exposure_hours` (6 h),
  $\mathrm{sinc}(u) = \sin u / u$ (rectangular exposure average);
- $t_1 = 0$ and $t_2$ = `pair_separation_hours` (120 h): the first pair midpoint is the
  phase origin;
- $N_\mathrm{ref}$ = `raptor.response.reference_norm`;
- $b(\alpha) = 0.03\,(\alpha/0.05)^{1.75}$ (`analysis.cloud`);
- $A$ and $\phi$ are the template means, interpolated linearly in $\alpha$ between template
  nodes ($\phi$ on the unit circle) and held constant outside the node range;
- $\phi_0$ is the unknown axion phase on a uniform grid of `phase_count` (64) values in
  $[0, 2\pi)$.

This is `evpa_templates` and `cloud_amplitude_scaling` in `eht_axion/analysis/model.py`.

### Feature selection

`scripts/analysis/limits.py` reads the blurred-image ensemble of each day and selects features without
looking at the size of the day differences. The EHT map $[y, x]$ is transposed to the
template order $[x, y]$; the sky axes must equal the shared grid.

2D features are template pixels that satisfy all of:

- the template is valid;
- on every day, the ensemble-mean Stokes image has $I \ge$ `intensity_fraction` (0.1) ×
  its maximum and $\lvert P\rvert = \sqrt{Q^2+U^2} \ge$ `polarization_fraction` (0.2) × its
  maximum, and the EVPA is finite in every draw;
- the distance from the grid centre is between `min_radius_rg` (2 $r_g$ = 7.6 µas) and
  the aperture radius (60 µas);
- both pixel indices are multiples of `pixel_stride` (5).

1D features are sectors $k$ with $k \bmod$ `sector_spacing_degrees` $=$
`sector_offset_degrees` (4°, 14°, …, 354°) such that the template is valid, the sector
EVPA $\chi_k$ is finite in every draw of every day, and, on every day, the draw means of
$B_k$ and $\langle\mathcal{P}\rangle_k$ (see [Aperture observable](#aperture-observable))
satisfy $\bar B_k \ge$ `intensity_fraction` $\times \max_k \bar B_k$ and
$\overline{\langle\mathcal{P}\rangle}_k \ge$ `polarization_fraction` $\times
\max_k \overline{\langle\mathcal{P}\rangle}_k$. All three statistics are computed from
each EHT draw with the aperture operator.

Then `build_observations` keeps only features whose circular resultant length
$R = \lvert\langle e^{2i\chi}\rangle\rvert$ over the draws is at least `resultant_min` (0.7)
on all four days.

### Observations and covariance

For day $d$ and the kept features, `moments` gives the circular mean
$\bar\chi_d = \tfrac12\arg\langle e^{2i\chi}\rangle$ and the population covariance
(ddof = 0) of the draws about it,

$$
C_d = \frac{1}{N}\sum_{k=1}^{N} r_k r_k^\mathsf{T},\qquad
r_k = \mathrm{wrap}_\pi(\chi_k - \bar\chi_d) - \langle\mathrm{wrap}_\pi(\chi - \bar\chi_d)\rangle .
$$

$C_d$ is the scatter of single draws (the D-term calibration uncertainty), not the
covariance of the mean. The data are the two day-pair differences
$d_1 = \mathrm{wrap}_\pi(\bar\chi_{96} - \bar\chi_{95})$ and
$d_2 = \mathrm{wrap}_\pi(\bar\chi_{101} - \bar\chi_{100})$, with covariances
$C_1 = C_{95} + C_{96}$ and $C_2 = C_{100} + C_{101}$ (days independent, no cross-pair
covariance). Each $C_p$ must be strictly positive definite; there is no regularization.
Branch proximity ($\pi/2 - \lvert d\rvert$ in units of $\sigma$) and eigenvalue
diagnostics are reported but do not change the selection.

### GRMHD background

`selected_background` estimates the intrinsic EVPA variability of the flow on exactly the
selected features. From the preprocessed axion-off frames (141 snapshots, state 0, 400-pixel
cubes resampled to the sky grid) it takes the EVPA $\theta_t$ of each selected feature
(2D: pixel EVPA; 1D: the same aperture operator), the lagged differences
$\Delta_t = \mathrm{wrap}_\pi(\theta_{t+L} - \theta_t)$ with $L$ =
`analysis.background.lag_frames` (3), and the pooled population standard deviation over
all features and frame pairs:

$$
\sigma_\mathrm{bg} = \sqrt{\langle(\Delta - \langle\Delta\rangle)^2\rangle}.
$$

The frame spacing comes from the GRMHD footer times recorded in the delivered cubes'
`generated.json` (1 $t_g$ = 8.89 h, so $L$ corresponds to 26.7 h, close to the one-day
separation of a pair). $\sigma_\mathrm{bg}$ describes a difference directly, so no
$\sqrt2$ factor is applied.

### Covariance models

| Mode | Pair covariance used in the likelihood |
|---|---|
| `std` | $\mathrm{diag}(C_p)$ |
| `cor` | $C_p$ |
| `cor_fixed_background` | $C_p + \sigma_\mathrm{bg}^2\,\mathbb{1}$ |

No mode includes a response-model variance proportional to $c^2$.

### Likelihood and phase marginalization

For each $\alpha$, `calculate_observed` (`eht_axion/analysis/likelihood.py`) evaluates on
the grid of $c$ and $\phi_0$

$$
\chi^2(c,\phi_0) = \sum_{p=1}^{2} \big(d_p - c\,s_p(\phi_0)\big)^\mathsf{T} C_p^{-1} \big(d_p - c\,s_p(\phi_0)\big),
$$

with one shared phase $\phi_0$ for both pairs, and marginalizes the phase with a uniform
discrete prior:

$$
\ln L(c) = \ln\frac{1}{N_\phi}\sum_{j=1}^{N_\phi} e^{-\chi^2(c,\phi_{0,j})/2}.
$$

Constants independent of $c$ and $\phi_0$ (the Gaussian log-determinant and normalization)
are omitted. The residual $d_p - c\,s_p$ is not wrapped again: the likelihood is a
Gaussian on the principal EVPA branch of the data, not a wrapped normal.

### Prior and posterior

The grid of $c$ has `coupling_points` (120) log-spaced values from `coupling_min` (1/137)
to `coupling_max` (1000). With `prior: log-uniform` the prior density is
$\pi(c) \propto 1/c$ on that interval (`uniform`: constant). `normalized_posterior_cdf`
(`eht_axion/analysis/posterior.py`) integrates $L(c)\pi(c)$ with the trapezoidal rule in
$c$; `posterior_upper_limit` inverts the cumulative distribution of the same piecewise-linear
density at `credibility` (0.95), giving $c_{95}$. Results are conditional on the finite grid
and prior; boundary diagnostics (probability mass in the first and last grid interval, peak
at a boundary, limit in the last interval) are reported.

### Zero-data control

For every case the same calculation is repeated with $d_1 = d_2 = 0$ and the same
covariances and features (`zero-data` in file names). It shows the limit that the
covariance and templates alone would give, for comparison with the observed limit.

### Checks

- `scripts/analysis/check_grid.py` re-evaluates the saved observations of a limits root with
  overridden `phase_count`, `coupling_points` and/or `coupling_max`; everything else must
  equal the experiment that produced the source. It reports $\Delta\log_{10}c_{95}$
  against the source. Run once with a denser grid ("grid check") and once with a wider
  range of $c$ ("prior check").
- `scripts/analysis/check_covariance.py` splits the draws of each day into the first and the last
  half (500 each), recomputes $C_d$ about the fixed full-ensemble means, and repeats the
  observed limits for all six cases; it reports $\Delta\log_{10}c_{95}$ against the
  primary result. Features and data stay fixed.

### Summary

`scripts/analysis/summarize.py` reads the completed limits, grid-check, prior-check and optional
covariance-check stages and writes: `limits.csv`; comparisons grid check vs. primary and
prior check vs. grid check (both must keep the same feature counts); the saved joint-grid
best fit of every observed case and every $\alpha$, with $\chi^2$ at $c = 0$ and at the best
grid point (`fit-chi2.csv`, one figure `fits-alpha-<α>.png` per $\alpha$); a limits figure;
`report.json`; and `RESULTS.md`. `RESULTS.md` is generated from the configured
`analysis.inference.credibility` and the numbers of the run (limits table, $\sigma_\mathrm{bg}$,
largest changes in the checks, fit table); it contains no fixed conclusions.

### Paper output

`scripts/analysis/paper.py` (stage `analysis:paper`) collects everything the paper takes from the
pipeline in one directory. It reads the limits, the three check stages, the templates (the
processed cube, response fit and templates at the paper's $\alpha$, and the processed
axion-off background), the EHT ensemble, the delivered axion-off 400-pixel cubes, the RAPTOR
plan and the build's copy of `axion.c`. The radial tables are the files recorded in the plan
(`--radial-low`, `--radial-high` or `EHT_RADIAL_LOW`, `EHT_RADIAL_HIGH` override the path; the
SHA-256 must equal the plan record), and the field window $W(r)=\exp[-(r/r_w)^n]$ is read from
`axion.c`. The limits must have been produced from the current `configs/experiment.yaml`. Nothing is
refitted. The output directory contains:

- `figures/`: eight vector PDFs, `radial_function`, `two_photon`, `signal_maps`,
  `template_maps`, `background`, `observational_uncertainty`, `data_vs_model` and `limits`;
- `numbers.tex`: one `\newcommand` per number quoted in the paper, named `\Num` followed by
  words without digits (for example `\NumLimitObsSectorsCorBgAlphaThirty`), each preceded by a
  comment with its description and source;
- `numbers.json` (schema `analysis-paper-numbers-v1`): the TeX text, the unrounded value and
  the source of every macro;
- `table-limits.tex` and `table-chi2.tex`: complete `tabular` environments with the observed
  $c_{95}$ (three significant digits) and with $\chi^2$ at $c=0$ and at the best grid point.

Derived quantities: the simulated flux is $\sum I\,\Delta\Omega$ of the delivered 400-pixel
axion-off images, with $\Delta\Omega$ from `raptor.camera_fov_rg` and the angular size of
$r_g$; the reconstructed flux, polarization fraction and $\lvert\beta_2\rvert$ are the
`global.npz` values of the ensemble (per-day means); the simulated polarization fraction and
$\lvert\beta_2\rvert$ apply `global_diagnostics` to the processed axion-off frames, transposed
to the EHT orientation; the grey histograms of the background figure pool the lagged
differences of every sector and of every pixel with a defined EVPA. The presentation choices
are constants of the script and are recorded in the manifest: `PAPER_ALPHA` = 0.30 and
`PAPER_SNAPSHOT_LABEL` = 2300 (response maps, template figure, $\chi^2$ table and fit figure),
the EVPA tick settings, the radii at which the radial profiles, the window and the orbital
period are quoted, the 15 µas separation of close pixel pairs and the two-emitter check
points.

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Experiment configuration](experiment-config.md)
- [EHT imaging](eht-imaging.md)
- [RAPTOR simulation](raptor-simulation.md)
- [Data layout](data-layout.md)
- [Glossary](glossary.md)
