# M87* EHT axion analysis — script pipeline

Scripts for two independent inputs — RAPTOR ray-traced images of a GRMHD flow with
an axion field, and EHT 2017 M87* polarimetric reconstructions — and the analysis that
combines them into phase-marginalized Bayesian limits on the axion–photon coupling c.
Not an installable package: run the scripts directly (or through `pixi run`). Chinese summary of the structure and
of the 2026-11 cleanup: [CLEANUP_ZH.md](CLEANUP_ZH.md).

## One parameter file

[`experiment.yaml`](experiment.yaml) holds every scientific value, each exactly once,
in four sections:

| section | used by | content |
|---|---|---|
| `shared` | more than one part | M87* mass/distance, label step, alphas, principal n, common sky grid (64², 120 µas) and the 19 µas matching blur |
| `raptor` | RAPTOR simulation (independent input) | camera width, build options, plasma, snapshots, response/background resolution, radial tables, BHAC grid |
| `eht_imaging` | EHT imaging (independent input) | Paper VII imaging: uvfits identity (name + sha256), weights, rounds, D-term statistics, draws; solver settings of the D refinement (`refine`) and of the D-term Monte Carlo (`monte_carlo`) |
| `analysis` | analysis (combines both inputs) | rotation, aperture, feature-selection rules, cloud model, background lag, inference grid/prior/credibility/times |

Each script reads `shared` plus its own section; the analysis also reads the raptor
and eht_imaging values it needs (e.g. camera width, response norm, draw count).
`eht_axion/experiment.py` validates the file and derives everything else (unit
conversions, the ray-tracing subset that defines RAPTOR job ids, the flat
eht-imaging parameter set). Library functions have no scientific default values.
Machine paths and execution options (workers, threads, memory) are not in the YAML:
they are CLI options or `.env` / `.env.raptor` entries (see `.env.example`,
`.env.raptor.example`).

## Run

RAPTOR simulation and EHT imaging do not depend on each other and can run in any
order or in parallel; the analysis needs the products of both.

```bash
pixi install --locked
pixi run check                                   # tests + lint

# RAPTOR simulation (environment raptor-build; paths from .env.raptor)
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json
pixi run -e raptor-build raptor-generate         # run jobs, deliver data/raptor/res-*/
pixi run -e raptor-build raptor-merge --plan plan.json --job-root RUN_A --job-root RUN_B --output-root MERGED

# EHT imaging (environment eht-rebuild; uvfits locations from .env)
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-refine --day 95 --source-root PD95 --i-root I95 --output-root D95
#   (each of the four days 95, 96, 100, 101; then)
pixi run -e eht-rebuild eht-export-bundle --i-root I95 --i-root I96 --i-root I100 --i-root I101 \
    --d-root D95 --d-root D96 --d-root D100 --d-root D101 --output-root BUNDLE
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC    # 4 days x 1000 D draws
pixi run -e eht-research eht-summarize --dataset-root MC_COPY --output-root ENSEMBLE

# Analysis of both (default environment)
pixi run analysis-templates --raptor-input-root MERGED/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root GRID \
    --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root PRIOR \
    --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance --source LIMITS --templates TEMPLATES --output-root COV
pixi run analysis-summarize --limits LIMITS --gridcheck GRID --priorcheck PRIOR --covcheck COV \
    --output-root SUMMARY
```

The analysis follows EHT Paper VII only: the EHT images are used as reconstructed
(Paper VII §5.2 centring = zero shift in this reproduction); there is no re-centring option.
Every output root must be new. Each step writes stage directories with a
`manifest.json` recording parameters, input/output SHA-256 and code hashes;
`pixi run verify DIR --recursive` re-checks a stage and its whole lineage.

## Code

The package mirrors the scripts: `eht_axion/raptor` serves `raptor_*`, `eht_axion/eht`
serves `eht_*`, `eht_axion/analysis` serves `analysis_*`.

```text
experiment.yaml        the parameter file
eht_axion/
  experiment.py        load/validate experiment.yaml, derived settings
  provenance.py        stage manifests, SHA-256 of inputs/outputs, atomic stage directories
  physics.py           constants, M87* unit scales, axion 211 frequency omega(alpha)
  evpa.py              EVPA validity and pi-periodic differences (shared by eht/ and analysis/)
  raptor/              RAPTOR production
    plan.py              task plan and job ids
    build.py             pinned RAPTOR build
    jobs.py affinity.py  job runner, CPU allocation
    outputs.py           validation of one job's HDF5/spectrum
    convert.py           job HDF5 -> Stokes cube
    delivery.py          delivery/merge of cubes; checks of delivered cubes
    layout.py            delivered-cube directory layout
  eht/                 EHT imaging
    calibration.py       Paper VII primitives: native data, fsys, D terms, D-term sampling
    convex.py            exact convex P solver with duality-gap certificate
    ensemble.py          FITS reader, sector apertures, per-day ensemble summaries
  analysis/            analysis of both inputs
    templates.py         simulated images -> response templates (and their reader)
    aperture.py          intensity-weighted aperture EVPA (1D observable)
    model.py             axion EVPA signal model
    observations.py      EHT features -> daily moments and day-pair differences
    likelihood.py        Gaussian likelihood, limit options, limits on c
    posterior.py         posterior CDF and upper limit on a grid
scripts/               raptor_*, eht_*, analysis_* entry points, verify_products
tests/                 synthetic unit/integration tests
docs/                  historical design notes (written for the earlier configs/ layout)
```

## Original project safety

The original project (notebooks, simulations, archives) lives separately in
`/home/canoziia/MyProject/PHYS/EHT_ZC`; this worktree contains none of its data.
