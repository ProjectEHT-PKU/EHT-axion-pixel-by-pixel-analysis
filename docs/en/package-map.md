# Package map

## Purpose

This page lists every module of the `eht_axion` package and every script in `scripts/`
with a one-line description and its main public functions, so that a reader can go from a
scientific step to the code that implements it. The package mirrors the scripts:
`eht_axion/raptor` serves the `raptor_*` scripts, `eht_axion/eht` the `eht_*` scripts and
`eht_axion/analysis` the `analysis_*` scripts; modules used by more than one part are at
the top level.

## Details

### Top-level modules

| Module | Purpose | Main functions and classes |
|---|---|---|
| `eht_axion/experiment.py` | Load and validate `experiment.yaml`; derived settings | `read_experiment`, `validate`, `validate_ray_tracing`, `physical_scale`, `ray_tracing_config`, `analysis_geometry`, `eht_imaging_parameters`, `ANGULAR_BINS`, `EHT_DAYS` |
| `eht_axion/provenance.py` | Stage manifests, SHA-256 of inputs and outputs, atomic stage directories, verification | `stage_directory`, `verify_stage`, `require_declared_outputs`, `verification_session`, `sha256`, `manifest_path`, `source_closure`, `contained_path` |
| `eht_axion/physics.py` | Physical constants, M87\* unit scales, axion 211 frequency | `PhysicalScale`, `omega` |
| `eht_axion/evpa.py` | EVPA of linear polarization and $\pi$-periodic differences, shared by `eht/` and `analysis/` | `finite_evpa`, `circular_difference` |

### RAPTOR production: eht_axion/raptor

| Module | Purpose | Main functions |
|---|---|---|
| `plan.py` | Read-only snapshot-footer and radial-table checks; task list and task identifiers | `plan_generation`, `task_identity`, `read_snapshot_footer`, `validate_radial_table` |
| `build.py` | Build of the pinned RAPTOR core with its official `setup.sh` | `build_raptor` |
| `jobs.py` | Concurrent, per-task atomic job runner with resume | `run_jobs`, `serialize_parameters`, `parameter_values` |
| `affinity.py` | Disjoint per-worker CPU sets and memory headroom check | `allocate_cpu_sets`, `check_memory`, `launch_prefix` |
| `outputs.py` | Validation of one job's HDF5 and spectrum | `validate_job_output`, `read_core_constants`, `hdf_mapping` |
| `convert.py` | Job HDF5 blocks to a float64 Stokes cube; small array-file helpers | `convert_hdf5` |
| `delivery.py` | Delivery and merge of cubes; checks of delivered cubes against the experiment | `deliver_stokes`, `generated_inputs`, `input_cadence_hours`, `validate_processed_source` |
| `layout.py` | Directory layout of the delivered cubes (`alpha-<α>/` and `background/`) | `alpha_label`, `product_directory`, `input_raptor_files` |

### EHT imaging: eht_axion/eht

| Module | Purpose | Main functions |
|---|---|---|
| `calibration.py` | Paper VII calibration primitives on ehtim 1.3.2: native data, systematic noise, fixed D-term correction, zero-prior forward D fit, four-day D distribution and sampling | `load_native`, `add_fsys`, `make_imaging_obs`, `apply_fixed_d`, `fit_forward_d`, `d_distribution`, `sample_dterms` |
| `convex.py` | Exact convex fixed-I, fixed-D P solver (PDHG) with a duality-gap certificate | `exact_cartesian_problem`, `solve_problem`, `CartesianProblem` |
| `ensemble.py` | Native FITS reader, Paper VII image-averaged quantities, per-day ensemble summary | `load_stokes_fits`, `global_diagnostics`, `summarize` |

### Analysis: eht_axion/analysis

| Module | Purpose | Main functions |
|---|---|---|
| `templates.py` | Preprocessing of simulated images, response fit, template summary and reader, template-root layout | `preprocess_frame`, `preprocess_cube`, `fit_phase_dft`, `fit_processed`, `summarize_fits`, `load_templates`, `RunLayout` |
| `aperture.py` | Intensity-weighted aperture EVPA (1D observable) and sector selection statistics, with one quadrature operator | `aperture_operator`, `integrate_stokes`, `sector_brightness` |
| `model.py` | Axion EVPA signal model | `cloud_amplitude_scaling`, `interpolate_fit_parameters`, `evpa_templates` |
| `observations.py` | Selected EHT features to daily moments and day-pair differences | `ObservationData`, `moments`, `build_observations`, `branch_risk`, `covariance_diagnostics` |
| `likelihood.py` | Inference options, Gaussian likelihood, phase marginalization, limits | `LimitOptions`, `MODES`, `calculate_observed` |
| `posterior.py` | Posterior density, CDF and upper limit on a grid | `normalized_posterior_cdf`, `posterior_upper_limit` |

### Scripts

| Script | Purpose |
|---|---|
| `raptor_build.py` | Build the pinned RAPTOR core into `RUN_ROOT/build` |
| `raptor_generate.py` | Plan tasks, run jobs and deliver Stokes cubes of one run root |
| `raptor_merge.py` | Deliver cubes from the job roots of several hosts into a new root |
| `eht_reconstruct.py` | EHT imaging stages `I` and `fiducial`, one day per output root |
| `eht_refine.py` | Exact convex P refinement and D refits of a fiducial stage |
| `eht_export_bundle.py` | Bundle the four I and refined stages as the Monte Carlo input |
| `eht_sample.py` | D-term Monte Carlo: one exact P image per D draw |
| `eht_summarize.py` | Per-day ensemble summary of the Monte Carlo images |
| `analysis_templates.py` | Preprocessing, response fits and templates from delivered RAPTOR data |
| `analysis_limits.py` | Feature selection, observations, GRMHD background, likelihood and limits |
| `analysis_check_grid.py` | Re-evaluate saved observations on another phase or coupling grid |
| `analysis_check_covariance.py` | Covariance sensitivity from two halves of the draws |
| `analysis_summarize.py` | Summary tables, comparisons and figures |
| `analysis_paper.py` | Figures, quoted numbers (`numbers.tex`, `numbers.json`) and tables of the paper |
| `_paper_figures.py` | Style and plotting functions of the paper figures |
| `verify_products.py` | Verify a completed stage (`pixi run verify`) |
| `_bootstrap.py` | Put the repository root on `sys.path` for direct script execution |
| `_environment.py` | Read `.env` and `.env.raptor` as data with a fixed set of allowed keys |
| `_raptor_options.py` | Command-line options shared by `raptor_build.py` and `raptor_generate.py` |
| `_data_root.py` | `--data-root` / `M87_AXION_DATA_ROOT` and the default roots of the data layout |

### Tests

`tests/` contains synthetic unit and integration tests, named after the module or script
they test (for example `test_raptor_plan.py`, `test_eht_convex.py`, `test_likelihood.py`,
`test_analysis_limits.py`). `test_docs_links.py` checks the links of this documentation.
See [Development](development.md).

## Related pages

- [Documentation index](index.md)
- [Overview](overview.md)
- [Data layout](data-layout.md)
- [Development](development.md)
