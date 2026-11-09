# Experiment configuration

## Purpose

`experiment.yaml` is the only parameter file of the experiment. Every scientific value is
written there exactly once, and every script reads it through `eht_axion/experiment.py`.
This page lists every key with its meaning, unit, current value and the code that reads it,
and describes validation and the derived settings.

## Inputs

- `experiment.yaml` in the repository root (the default of every `--config` option).
- An alternative file can be given with `--config PATH`; it must contain exactly the same
  keys.

## Outputs

Nothing is written by the loader. The validated dictionary, or the parts of it that a stage
used, is recorded in stage manifests and plans ([Provenance](provenance.md)):

- `raptor_generate.py` and delivery record the ray-tracing subset (`ray_tracing_config`).
- The EHT imaging stages record the flat parameter set (`eht_imaging_parameters`).
- `analysis_limits.py` records the whole file as `scientific_profile`.

## How to run

The file is read by the scripts; there is no separate command. A quick validation:

```bash
pixi run python -c "from eht_axion.experiment import read_experiment; read_experiment()"
```

## Details

### Structure and reading rule

The file has exactly four sections:

- `shared`: values used by more than one part.
- `raptor`: what determines the ray-traced Stokes images.
- `eht_imaging`: the EHT imaging method, its data identity and its solver settings.
- `analysis`: from templates and EHT ensembles to limits on $c$.

Each script reads `shared` plus its own section; the analysis also reads the `raptor` and
`eht_imaging` values it needs. Machine paths and execution switches are not in the file
(see [Development](development.md) for `.env` and `.env.raptor`).

### Validation

`read_experiment` parses the file with a YAML loader that rejects duplicate keys, then
`validate` checks it. Validation checks structure (exact key sets per section), types
(booleans are not numbers, integers are integers) and ranges. It fills no defaults. Two
cross-checks depend on derived values: the analysis sky field (`shared.image.fov_uas`
converted to $r_g$) must fit inside `raptor.camera_fov_rg`, and
`analysis.features.min_radius_rg` must lie inside the aperture radius.

### Derived settings

| Function in `eht_axion/experiment.py` | Result | Used by |
|---|---|---|
| `physical_scale` | `PhysicalScale`: $t_g = GM/c^3$ (8.893 h for the values below), angular size of $r_g$ (3.819 µas), conversions | analysis scripts, delivery cadence |
| `ray_tracing_config` | `shared` physics keys (`mass_solar`, `distance_mpc`, `label_step_tg`, `alphas`, `principal_n`) plus the whole `raptor` section | RAPTOR plan, delivery identity (`generated.json`) |
| `analysis_geometry` | preprocessing arguments (`sigma`, `init_fov`, `final_fov`, `final_shape`, `rotation_angle`, in $r_g$ and radians) and aperture arguments (`fov`, `radial_max`, `window_degrees`, `angular_bins`) | `analysis_templates.py`, `analysis_limits.py` |
| `eht_imaging_parameters` | flat dictionary of the `eht_imaging` values (without `refine` and `monte_carlo`), the sky grid (`npix`, `fov_uas`, `blur_fwhm_uas`) and four machine locations | EHT imaging stages, `eht_refine.py`, `eht_export_bundle.py` |

`ANGULAR_BINS = 360` (one-degree aperture bins) and `EHT_DAYS = (95, 96, 100, 101)` are
code constants in the same module, not YAML values. `EHT_DAYS` is the only list of the
observing days in the code (the EHT imaging scripts, the ensemble summary and the
analysis import it); validation requires `eht_imaging.inputs` to have exactly these days.

### Section shared

| Key | Value | Unit | Meaning | Read by |
|---|---|---|---|---|
| `mass_solar` | 6.5e9 | $M_\odot$ | M87\* mass | `physical_scale`; RAPTOR `MBH` |
| `distance_mpc` | 16.8 | Mpc | distance to M87\* | `physical_scale`; RAPTOR `DISTANCE` (written in kpc) |
| `label_step_tg` | 0.1 | $t_g$ per label | time per unit of GRMHD snapshot label | plan check of the footer time cadence |
| `alphas` | 0.10 … 0.44 (8 values) | – | gravitational fine-structure constants $\alpha$ | RAPTOR tasks and delivery, templates, likelihood grid |
| `principal_n` | 2 | – | principal quantum number $n$ of the 211 state in $\Omega(\alpha)=\alpha\,(1-\alpha^2/(2n^2))$ | RAPTOR `AXION_OMEGA`; signal model |
| `image.shape` | 64 | pixels | side of the common sky grid | EHT imaging `npix`; analysis resampling; ensemble summary |
| `image.fov_uas` | 120.0 | µas | field of view of the common sky grid (1.875 µas pixels) | EHT imaging (replaces 128 µas in the official script); analysis resampling and aperture |
| `image.blur_fwhm_uas` | 19.0 | µas | FWHM of a circular Gaussian applied to the EHT images (the "blurred" images) and to the simulated images; an additional matching kernel, not an effective PSF | EHT `save_image`; analysis preprocessing $\sigma$ |

### Section raptor

| Key | Value | Unit | Meaning | Read by |
|---|---|---|---|---|
| `camera_fov_rg` | 46.2261 | $r_g$ | camera width (176.5 µas) | RAPTOR `CAM_SIZE_X/Y`; preprocessing `init_fov` |
| `build.model` | `bhac` | – | `setup.sh --code` | `raptor/build.py` |
| `build.metric` | `mks` | – | `--metric` | `raptor/build.py` |
| `build.integrator` | `rk4` | – | `--int` | `raptor/build.py` |
| `build.grid` | `SMR` | – | `--grid` | `raptor/build.py` |
| `build.polarized` | true | – | `--rad=pol` | `raptor/build.py` |
| `build.sfc` | true | – | `--sf=sfc` | `raptor/build.py` |
| `build.magnetic_polarity` | `plus` | – | `--bflip` | `raptor/build.py` |
| `build.sigma_cut` | 1.0 | – | magnetization cut, `--sigma-cut` | `raptor/build.py` |
| `build.thetae_min` | 0.001 | – | `--thetae-min` | `raptor/build.py` |
| `build.thetae_max` | 100.0 | – | `--thetae-max` | `raptor/build.py` |
| `build.uniform_output` | false | – | `--uniform` | `raptor/build.py` |
| `plasma.mass_unit_g` | 3.0e25 | g | mass unit of the GRMHD flow, `M_UNIT` | RAPTOR `model.in` |
| `plasma.r_low` | 1.0 | – | electron temperature parameter `R_LOW` | RAPTOR `model.in` |
| `plasma.r_high` | 1.0 | – | electron temperature parameter `R_HIGH` | RAPTOR `model.in` |
| `plasma.inclination_degrees` | 163.0 | degree | observer inclination | RAPTOR `model.in` |
| `plasma.frequency_hz` | 2.3e11 | Hz | observing frequency (one frequency) | RAPTOR `model.in`; delivery dataset names |
| `plasma.step_size` | 0.005 | – | geodesic integration step, `STEPSIZE` | RAPTOR `model.in` |
| `snapshots.start`, `stop`, `step` | 1600, 3000, 10 | label | inclusive snapshot labels, 141 snapshots | plan, delivery, label checks in the analysis |
| `response.resolution` | 100 | pixels | side of axion-on images and their axion-off baseline | plan, delivery, `analysis_templates.py` |
| `response.reference_norm` | 0.1 | – | axion field normalization `AXION_NORM`; templates are per this norm | RAPTOR `model.in`; signal model |
| `background.resolution` | 400 | pixels | side of the axion-off images for the GRMHD background | plan, delivery, `analysis_templates.py` |
| `radial_tables.low.spin` | 0.9375 | – | declared spin of the low radial table | plan, task identity |
| `radial_tables.high.spin` | 0.99 | – | declared spin of the high table (intentionally differs from the fluid spin 0.9375) | plan, task identity |
| `alpha_table` | 0.10–0.30 `low`, 0.35–0.44 `high` | – | radial table used for each $\alpha$ | plan |
| `grid.nxlone`, `xprobmin`, `xprobmax`, `hslope` | see file | – | BHAC `grid_mks.in` coordinates; compared with the build's `grid_mks.in` and with snapshot footers | plan, job runner |

### Section eht_imaging

| Key | Value | Unit | Meaning | Read by |
|---|---|---|---|---|
| `inputs` | four `[file name, sha256]` | – | identity of the HOPS uvfits per day | `eht_reconstruct.py` (checked against the input receipt and the file bytes) |
| `gain_tol` | [0.02, 0.2] | – | lower and upper fractional gain tolerance of the shared-gain self-calibration (the value of the official script; checked against the recorded I stages) | `eht_reconstruct.py` |
| `fsys` | 0.005 | – | fractional systematic noise added in quadrature to RL/LR errors: $\sigma^2 \to \sigma^2 + (f_\mathrm{sys}\,\lvert(RR+LL)/2\rvert)^2$ | `add_fsys` in fiducial, refine and Monte Carlo; bundle contract |
| `solve_sites` | AZ, LM, PV | – | stations whose D terms are solved | fiducial, refine, Monte Carlo D distribution |
| `fixed_sites` | AA, AP, SM, JC | – | stations already leakage-corrected in the HOPS products (the only list of them in the code) | `load_native` metadata; `sample_dterms` zero-mean residual draws |
| `p_data` | pvis 1.0, m 1.0 | – | ehtim data-term weights of P imaging | fiducial `Imager`; exact P solver of `eht_refine.py` and `eht_sample.py`; bundle contract |
| `p_reg` | hw 100.0, ptv 1.0 | – | ehtim regularizer weights (entropy, total variation) | as `p_data` |
| `norm_reg` | true | – | ehtim regularizer normalization | as `p_data` |
| `initial_m` | 0.1 | – | initial fractional polarization of the fiducial P image | `eht_reconstruct.py` |
| `initial_evpa_deg` | 0.0 | degree | initial EVPA of the fiducial P image | `eht_reconstruct.py` |
| `rounds` | 50 | – | fiducial P/D alternation rounds | `eht_reconstruct.py` |
| `restart_rounds` | 2, 12, 22, 32, 42 | – | rounds that restart P from a blurred fractional polarization (numerical, not from the paper) | `eht_reconstruct.py` |
| `restart_blur_uas` | 20.0 | µas | FWHM of that restart blur | `eht_reconstruct.py` |
| `restart_domain_epsilon` | 1e-8 | – | restart P is scaled by $1-\epsilon$ so that $\lvert P/I\rvert<1$ | `eht_reconstruct.py` |
| `maxit` | 500 | – | ehtim optimizer iteration cap per call (numerical) | `eht_reconstruct.py` |
| `stop` | 1e-8 | – | ehtim convergence tolerance | `eht_reconstruct.py` |
| `shift_objective_constant` | true | – | add $\sum$`p_data` $+\,h_w\ln 2$ to the ehtim objective (constant; minimizer unchanged) | `eht_reconstruct.py` |
| `uv_min` | 1.0e8 | wavelengths | baselines shorter than 100 Mλ (intra-site "zero" baselines) are removed before the shared gain solve and in all P imaging | I stage, fiducial, refine, Monte Carlo |
| `uncertainty_draws` | 1000 | – | D-term draws per day | `eht_sample.py`, `eht_summarize.py`, `analysis_limits.py`, `analysis_check_covariance.py` |
| `seed` | 20210405 | – | Monte Carlo seed; per day `SeedSequence([seed, day])` | `eht_sample.py` |
| `covariance_ddof` | 1 | – | delta degrees of freedom of the four-day D covariance | `eht_sample.py`; bundle contract |
| `residual_d_sigma` | 0.01 | – | per real/imaginary component SD of zero-mean D draws for `fixed_sites` | `eht_sample.py`; bundle contract |
| `residual_d_convention` | `per_component_sd` | – | the only allowed convention for `residual_d_sigma`; kept because it is part of the recorded imaging parameter set | validation, recorded only |
| `refine.outer_rounds` | 100 | – | P solve / D refit alternations before the final P solve | `eht_refine.py` |
| `refine.gap_atol` | 1e-6 | objective units | absolute duality-gap tolerance of each exact P solve | `eht_refine.py` |
| `refine.maxiter` | 20000 | – | P solver iterations per call | `eht_refine.py` |
| `refine.check_every` | 500 | – | iterations between gap certificates | `eht_refine.py` |
| `refine.time_limit` | 300.0 | s | P solver time per call | `eht_refine.py` |
| `refine.d_optimizer` | ftol 1e-13, gtol 1e-9, maxiter 10000, maxcor 50, maxls 50, maxfun 100000 | – | scipy L-BFGS-B options of every D refit | `eht_refine.py` → `fit_forward_d` |
| `monte_carlo.initial_p_fraction` | 0.05 | – | initial P $=0.05\,I$ (real, EVPA 0), no warm start | `eht_sample.py` |
| `monte_carlo.gap_atol` | 1e-6 | objective units | gap tolerance per draw | `eht_sample.py` |
| `monte_carlo.maxiter` | 50000 | – | P solver iterations, first attempt | `eht_sample.py` |
| `monte_carlo.check_every` | 500 | – | iterations between certificates | `eht_sample.py` |
| `monte_carlo.time_limit` | 600.0 | s | time, first attempt | `eht_sample.py` |
| `monte_carlo.retry_maxiter` | 200000 | – | iterations of the second attempt of an uncertified draw | `eht_sample.py` |
| `monte_carlo.retry_time_limit` | 1800.0 | s | time of the second attempt | `eht_sample.py` |

The compact flux density of the official I workflow (0.6 Jy, `zbl = 0.60`) is fixed inside
the official EHT script, which runs unchanged apart from two substitutions
([EHT imaging](eht-imaging.md#i-stage)); P imaging uses the total flux of the fixed I image.

`refine` and `monte_carlo` are recorded in the plans of `eht_refine.py` and
`eht_sample.py`, not in the flat parameter set of the imaging stages.

### Section analysis

| Key | Value | Unit | Meaning | Read by |
|---|---|---|---|---|
| `rotation_degrees` | 108.0 | degree | rotation of the simulated images onto the sky (image and $Q+iU$) | `analysis_geometry` |
| `aperture.window_degrees` | 10 | degree | full width of each aperture sector | theory aperture operator; `eht_summarize.py` |
| `aperture.radial_max_fraction` | 1.0 | – | aperture radius as a fraction of the half field (60 µas) | as above |
| `features.intensity_fraction` | 0.1 | – | keep where mean $I \ge 0.1 \times$ daily maximum | `analysis_limits.py` |
| `features.polarization_fraction` | 0.2 | – | and mean $\lvert P\rvert \ge 0.2 \times$ daily maximum (0 selects by brightness only) | `analysis_limits.py` |
| `features.min_radius_rg` | 2.0 | $r_g$ | minimum distance of a 2D pixel from the image centre (7.6 µas) | `analysis_limits.py` |
| `features.pixel_stride` | 5 | pixels | 2D: every 5th pixel in both axes | `analysis_limits.py` |
| `features.sector_spacing_degrees` | 10 | degree | 1D: one sector every 10° | `analysis_limits.py` |
| `features.sector_offset_degrees` | 4 | degree | 1D: starting at 4° east of north | `analysis_limits.py` |
| `features.resultant_min` | 0.7 | – | minimum circular resultant length on every day | `build_observations` |
| `cloud.coefficient` | 0.03 | – | $b(\alpha)=0.03\,(\alpha/0.05)^{1.75}$ | `cloud_amplitude_scaling` |
| `cloud.pivot_alpha` | 0.05 | – | pivot of $b(\alpha)$ | as above |
| `cloud.exponent` | 1.75 | – | exponent of $b(\alpha)$ | as above |
| `background.lag_frames` | 3 | frames | lag of the axion-off frame differences (3 frames $\approx$ 26.7 h) | `selected_background` |
| `inference.phase_count` | 64 | – | uniform grid of the marginalized axion phase | `evpa_templates` |
| `inference.coupling_points` | 120 | – | log-spaced grid of $c$ | `calculate_observed` |
| `inference.coupling_min` | 1/137 | – | lower bound of the prior on $c$ | as above |
| `inference.coupling_max` | 1000.0 | – | upper bound of the prior on $c$ | as above |
| `inference.prior` | `log-uniform` | – | prior on $c$ (`uniform` also allowed) | posterior |
| `inference.credibility` | 0.95 | – | posterior quantile reported as the limit | posterior |
| `inference.exposure_hours` | 6.0 | h | rectangular time average of one observation | `evpa_templates` |
| `inference.day_gap_hours` | 24.0 | h | separation of the two days of a pair | `evpa_templates` |
| `inference.pair_separation_hours` | 120.0 | h | separation of the two pair midpoints | `evpa_templates` |

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Analysis](analysis.md)
- [EHT imaging](eht-imaging.md)
- [RAPTOR production](raptor-production.md)
- [Development](development.md)
