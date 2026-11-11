# Analysis

[← Documentation index](index.md)

The analysis pipeline combines delivered `RAPTOR` Stokes cubes with the summarized EHT polarimetric ensembles to place phase-marginalized Bayesian upper limits $c_{95}(\alpha)$ on the dimensionless coupling $c$. All analysis scripts run in the `default` Pixi environment.

## Running the analysis

With `EHT_DATA_ROOT` set in `.env` (or `--data-root`), all input and output paths default to the [standard data layout](data-layout.md#default-layout); `analysis-paper` additionally reads `src/axion.c` from the pinned `RAPTOR` checkout (`EHT_RAPTOR_SOURCE_ROOT`) and the radial tables from the [input layout](data-layout.md#input-layout) (`EHT_INPUT_ROOT`):

```bash
pixi run analysis-templates
pixi run analysis-limits
pixi run analysis-check-grid --check grid --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --check prior --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance
pixi run analysis-summarize
pixi run analysis-paper
```

Explicit directory flags override any layout default:

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
    --raptor-plan plan.json --raptor-source-root RAPTOR --output-root PAPER
```

When `--data-root` / `EHT_DATA_ROOT` is not set, `--covcheck` is optional in `scripts/analysis/summarize.py`; with a data root, `summarize.py` defaults `--covcheck` to `$EHT_DATA_ROOT/analysis/check-covariance`.

## Simulation templates (`scripts/analysis/templates.py`)

### Image preprocessing

`scripts/analysis/templates.py` verifies every delivered cube against `configs/experiment.yaml` (`generated_inputs`) and streams frames through `preprocess_cube` (`eht_axion/analysis/templates.py`). For each frame, `preprocess_frame`:
1. Convolves each Stokes channel on the native camera grid with a circular Gaussian of standard deviation $\sigma = \mathrm{FWHM}/\sqrt{8\ln 2}$ (`shared.image.blur_fwhm_uas` $= 19\,\mu\mathrm{as} \to 2.113\,r_g$) using nearest-value boundary extension (`mode="nearest"`).
2. Resamples from `raptor.camera_fov_rg` ($46.2261\,r_g \approx 176.5\,\mu\mathrm{as}$) onto the $64 \times 64$, $120\,\mu\mathrm{as}$ ($31.42\,r_g$) sky grid by bilinear interpolation (`order=1`, `prefilter=False`), actively rotating the image by $\theta =$ `analysis.rotation_degrees` ($108^\circ$).
3. Rotates linear polarization by the spin-2 factor $Q + iU \mapsto e^{2i\theta}(Q + iU)$.

No flux renormalization is applied. Processed maps have spatial axes `[x, y]` with axis `0` pointing west (`-East`) and axis `1` pointing north (`+North`). Response cubes (`analysis:preprocess`) retain all 17 states; the 400-pixel axion-off background (`analysis:background`) retains state `0` only.

### Sector aperture operator (1D observable)

The 1D observable is the coherent intensity-weighted EVPA integrated over $360$ overlapping "pizza-slice" sectors centred $1^\circ$ apart:

$$
\chi_{\mathrm{ap}}(\varphi_k) = \frac{1}{2}\arg \int_0^{r_{\max}} \int_{\varphi_k - w/2}^{\varphi_k + w/2} I\,(Q + iU)\,r\,\mathrm{d}\varphi\,\mathrm{d}r,
$$

with sector width $w =$ `analysis.aperture.window_degrees` ($10^\circ$) and outer radius $r_{\max} =$ `radial_max_fraction` $\times$ half-FOV ($60\,\mu\mathrm{as} \approx 15.71\,r_g$). `aperture_operator` (`eht_axion/analysis/aperture.py`) constructs a $(360, 64^2)$ CSR sparse matrix of bilinear pixel quadrature weights $W_{kj}$ (radial step `0.02` in `fov` units, i.e. $r_g$; angular step $0.2^\circ$; trapezoidal half-weights at radial and angular endpoints). After applying the operator, bins are rolled by $-90$ so bin $k \in \{0, \dots, 359\}$ is centred $k^\circ$ east of north:

- **Sector EVPA** (`integrate_stokes`, Paper VII eq. 11): $\chi_k = \frac{1}{2}\arg\sum_j W_{kj}\,I_j(Q_j + iU_j)$.
- **Sector intensity** (`sector_brightness`): $B_k = \sum_j W_{kj}\,I_j$.
- **Intensity-weighted polarized brightness** (`sector_brightness`, Paper VII eq. 10): $\langle\mathcal{P}\rangle_k = \sum_j W_{kj}\,I_j\lvert P_j\rvert / B_k$ with $\lvert P_j\rvert = \sqrt{Q_j^2 + U_j^2}$.

The identical sparse operator $W_{kj}$ is applied to the simulated frames and (after transposing `[y, x] -> [x, y]`) to every EHT draw.

### Harmonic response fit and snapshot summary

For each snapshot and $\alpha$, `fit_processed` (stage `analysis:fit`) drops state `0` and fits the 16 phase states ($\psi_k = 2\pi k/16$, $k = 0, \dots, 15$) at every pixel (2D) and sector (1D) via `fit_phase_dft`:

$$
\chi_k = \bar{\chi} + A\cos(\psi_k + \phi).
$$

Values are first lifted onto the $\pi$-branch of their circular mean $\chi_0 = \frac{1}{2}\arg\langle e^{2i\chi_k}\rangle$, giving $\bar{\chi} = \langle\tilde{\chi}_k\rangle$ and $A e^{i\phi} = \frac{1}{8}\sum_{k=0}^{15} (\tilde{\chi}_k - \bar{\chi})\,e^{-2\pi i k/16}$. Any missing phase invalidates that location; $A \le 64\varepsilon_{\mathrm{mach}}$ yields $A = 0$ and $\phi = \mathrm{NaN}$.

`summarize_fits` (stage `analysis:templates`) aggregates across the 141 snapshots for each $\alpha$ (requiring $\ge 2$ finite snapshots): Gaussian MLE (`ddof = 0`) `amplitude_mean` and `amplitude_std`, plus circularly lifted `phase_mean` and `phase_std`. `load_templates` marks a location valid only if all four statistics are finite and both standard deviations are non-negative at every $\alpha$; only the means enter the signal model.

## Observational features and covariance (`scripts/analysis/limits.py`)

### Feature selection

`prepare_limit_inputs` transposes each EHT ensemble map `[y, x]` to template order `[x, y]` (verifying exact grid match) and selects features blind to the observed day-to-day EVPA differences:

- **2D pixel features** must have a valid template, radial distance from the grid centre in $[\texttt{min\_radius\_rg}, r_{\max}] = [7.64\,\mu\mathrm{as}, 60\,\mu\mathrm{as}]$, both spatial indices divisible by `pixel_stride` (`5`), finite EVPA in all 1000 draws of every day, and daily draw-mean Stokes maps satisfying $I \ge 0.1\,I_{\max}$ (`intensity_fraction`) and $\lvert P\rvert \ge 0.2\,\lvert P\rvert_{\max}$ (`polarization_fraction`) on all four days.
- **1D sector features** are candidate bins $k \equiv 4 \pmod{10}$ ($4^\circ, 14^\circ, \dots, 354^\circ$; `sector_offset_degrees`, `sector_spacing_degrees`) that have a valid template, finite sector EVPA $\chi_k$ across all draws of every day, and daily draw-means satisfying $\bar{B}_k \ge 0.1\max_j \bar{B}_j$ and $\overline{\langle\mathcal{P}\rangle}_k \ge 0.2\max_j \overline{\langle\mathcal{P}\rangle}_j$ on all four days.

`build_observations` (`eht_axion/analysis/observations.py`) then retains only features whose circular resultant length $R = \lvert\langle e^{2i\chi}\rangle\rvert$ across draws is at least `resultant_min` ($0.7$) on all four days.

### Daily moments, day-pair differences, and GRMHD background

For each observing day $d \in \{95, 96, 100, 101\}$, `moments` computes the $\pi$-periodic circular mean $\bar{\chi}_d = \frac{1}{2}\arg\langle e^{2i\chi}\rangle$ and the centered population D-draw covariance (`ddof = 0`):

$$
C_d = \frac{1}{N}\sum_{k=1}^{N} r_k r_k^{\mathsf{T}}, \qquad r_k = \operatorname{wrap}_\pi(\chi_k - \bar{\chi}_d) - \bigl\langle\operatorname{wrap}_\pi(\chi - \bar{\chi}_d)\bigr\rangle.
$$

$C_d$ represents single-draw D-term calibration uncertainty (not the standard error of the Monte Carlo mean). Treating days as independent yields two data vectors and block-diagonal pair covariances:

$$
d_1 = \operatorname{wrap}_\pi(\bar{\chi}_{96} - \bar{\chi}_{95}), \quad C_1 = C_{95} + C_{96}; \qquad
d_2 = \operatorname{wrap}_\pi(\bar{\chi}_{101} - \bar{\chi}_{100}), \quad C_2 = C_{100} + C_{101}.
$$

Both $C_1$ and $C_2$ must pass an unregularized Cholesky factorization. Proximity to the $\pm\pi/2$ branch cut ($\pi/2 - \lvert d\rvert$ in units of marginal $\sigma$) and eigenvalue condition numbers are recorded as diagnostics without altering feature selection.

**GRMHD background (`selected_background`):** on the exact retained 1D or 2D features, the 141 preprocessed 400-pixel axion-off frames yield EVPA time series $\theta_t$ and lagged differences $\Delta_t = \operatorname{wrap}_\pi(\theta_{t+L} - \theta_t)$ at lag $L =$ `analysis.background.lag_frames` ($3$ frames $= 3\,t_g = 26.68\,\mathrm{h}$, derived from the BHAC footer times). The pooled centered population standard deviation across all retained features and frame pairs,

$$
\sigma_{\mathrm{bg}} = \sqrt{\bigl\langle(\Delta - \langle\Delta\rangle)^2\bigr\rangle},
$$

gives $\sigma_{\mathrm{bg}}$ for a ~1-day difference directly, so no $\sqrt{2}$ factor is applied.

## Statistical inference

### Two-epoch signal model

In `eht_axion/analysis/model.py` (`evpa_templates`, `cloud_amplitude_scaling`), the model EVPA difference in pair $p \in \{1, 2\}$ at feature $x$ scales linearly with $c$:

$$
s_p(x; \alpha, \phi_0) = b(\alpha)\,\frac{2\sin(\omega\tau / 2)\,\operatorname{sinc}(\omega T / 2)}{N_{\mathrm{ref}}}\,A(x, \alpha)\,\sin\!\bigl(\phi_0 + \phi(x, \alpha) + \omega t_p\bigr),
$$

where $\omega = \Omega(\alpha)/t_g$ ($\mathrm{rad\,h^{-1}}$), $\tau = 24\,\mathrm{h}$ (`day_gap_hours`), $T = 6\,\mathrm{h}$ (`exposure_hours`), $\operatorname{sinc}(u) = \sin u / u$, $t_1 = 0$ and $t_2 = 120\,\mathrm{h}$ (`pair_separation_hours`, placing the phase origin at the first pair midpoint), $N_{\mathrm{ref}} = 0.1$ (`raptor.response.reference_norm`), and $b(\alpha) = 0.03\,(\alpha / 0.05)^{1.75}$. Template means $A$ and $\phi$ are linearly interpolated in $\alpha$ between nodes ($\phi$ on the unit circle) and clamped outside the node range.

### Covariance modes, likelihood, and upper limits

`calculate_observed` (`eht_axion/analysis/likelihood.py`) evaluates three pair-covariance models $\tilde{C}_p$ (none includes a $c^2$-dependent response-model variance):

| Mode | Effective pair covariance $\tilde{C}_p$ |
|---|---|
| `std` | $\operatorname{diag}(C_p)$ |
| `cor` | $C_p$ |
| `cor_fixed_background` | $C_p + \sigma_{\mathrm{bg}}^2\,\mathbb{1}$ |

On the joint grid of $c$ (`120` log-spaced points on $[1/137, 1000]$) and shared phase $\phi_{0,j} = 2\pi j / N_\phi$ ($N_\phi = 64$, $j = 0, \dots, N_\phi - 1$), the Gaussian quadratic form (residuals $d_p - c\,s_p$ are not re-wrapped; fixed log-determinant constants are omitted) is marginalized over $\phi_0$ with a uniform discrete prior:

$$
\chi^2(c, \phi_0) = \sum_{p=1}^{2} \bigl(d_p - c\,s_p(\phi_0)\bigr)^{\mathsf{T}} \tilde{C}_p^{-1} \bigl(d_p - c\,s_p(\phi_0)\bigr), \qquad
\ln L(c) = \ln \frac{1}{N_\phi} \sum_{j=0}^{N_\phi - 1} e^{-\chi^2(c, \phi_{0,j}) / 2}.
$$

`normalized_posterior_cdf` (`eht_axion/analysis/posterior.py`) multiplies $L(c)$ by $\pi(c) \propto 1/c$ (`log-uniform`, or constant for `uniform`), integrates with the trapezoidal rule in $c$, and `posterior_upper_limit` inverts the cumulative distribution of that piecewise-linear density analytically at `credibility = 0.95` to obtain $c_{95}$. Every dimension/mode pair is evaluated both on the observed differences (`observed`) and on a **zero-data control** ($d_1 = d_2 = 0$ with identical features and covariances, `zero-data`).

## Sensitivity checks, summary, and paper outputs

- **Grid and prior checks (`scripts/analysis/check_grid.py`)**: re-evaluates the saved `observations-{1d,2d}.npz` of `analysis:limits` with overridden `--phase-count`, `--coupling-points`, and/or `--coupling-max` (stage `analysis:check-grid`, reporting $\Delta\log_{10} c_{95}$ against the source stage without re-selecting features or re-hashing raw upstream lineage).
- **Split-ensemble covariance check (`scripts/analysis/check_covariance.py`)**: splits each day's 1000 draws into `first500` and `last500`, recomputes $C_d$ around the fixed full-ensemble circular means $\bar{\chi}_d$, and reports $\Delta\log_{10} c_{95}$ across all six observed cases (`analysis:check-covariance`).
- **Summary (`scripts/analysis/summarize.py`)**: verifies `limits`, `gridcheck`, `priorcheck`, and optional `covcheck`, asserts unchanged feature counts across checks, and writes `limits.csv`, `grid_vs_primary.csv`, `prior_vs_dense_grid.csv`, `fit-chi2.csv` ($\chi^2$ at $c=0$ and at the joint-grid best fit), `limits.png`, `fits-alpha-<α>.png`, `report.json`, and a purely data-driven `RESULTS.md` (`analysis:summary`).
- **Paper outputs (`scripts/analysis/paper.py`)**: collects all manuscript assets into stage `analysis:paper` without refitting anything:
  - `figures/`: eight deterministic vector PDFs (`radial_function`, `two_photon`, `signal_maps`, `template_maps`, `background`, `observational_uncertainty`, `data_vs_model`, `limits`);
  - `numbers.tex` and `numbers.json` (`analysis-paper-numbers-v1`): digit-free `\newcommand{\Num...}{...}` macros paired with unrounded values and exact data sources;
  - `table-limits.tex` and `table-chi2.tex`: LaTeX `tabular` environments for observed $c_{95}$ and $\chi^2$ at `PAPER_ALPHA = 0.30`.
  Presentation constants (`PAPER_ALPHA = 0.30`, `PAPER_SNAPSHOT_LABEL = 2300`, EVPA tick settings, radial sample radii `5, 10, 20` $r_g$, close-pixel threshold `15` $\mu\mathrm{as}$, and two-emitter check parameters) are recorded in the manifest, and the radial window $W(r) = \exp[-(r/r_w)^n]$ is parsed directly from `src/axion.c` of the `RAPTOR` checkout (`--raptor-source-root` / `EHT_RAPTOR_SOURCE_ROOT`), which must be clean at `raptor.source.commit` (the commit identifies the file, so it is not a stage input). The radial tables come from `--radial-low` / `--radial-high`, else the [input layout](data-layout.md#input-layout), else the paths in `raptor/plan.json`, and must match the SHA-256 recorded in the plan.
