# M87* EHT axion analysis — script pipeline

Scripts that turn (1) RAPTOR ray-traced images of a GRMHD flow with an axion field
and (2) EHT 2017 M87* polarimetric reconstructions into (3) phase-marginalized
Bayesian limits on the axion–photon coupling c. Not an installable package: run
the scripts directly (or through `pixi run`). Chinese summary of the structure and
of the 2026-11 cleanup: [CLEANUP_ZH.md](CLEANUP_ZH.md).

## One parameter file

[`experiment.yaml`](experiment.yaml) holds every scientific value, each exactly once,
in four sections:

| section | used by | content |
|---|---|---|
| `shared` | all processes | M87* mass/distance, label step, alphas, principal n, common sky grid (64², 120 µas) and the 19 µas matching blur |
| `raptor` | process 1 | camera width, build options, plasma, snapshots, response/background resolution, radial tables, BHAC grid |
| `eht_imaging` | process 2 | Paper VII imaging: uvfits identity (name + sha256), weights, rounds, D-term statistics, draws |
| `analysis` | process 3 | rotation, aperture, feature-selection rules, cloud model, background lag, inference grid/prior/times |

A process reads `shared` and the sections of the processes before it.
`eht_axion/experiment.py` validates the file and derives everything else (unit
conversions, the ray-tracing subset that defines RAPTOR job ids, the flat
eht-imaging parameter set). Machine paths are not in the YAML: they are CLI
options or `.env` / `.env.raptor` entries (see `.env.example`, `.env.raptor.example`).

## Run, in process order

```bash
pixi install --locked
pixi run check                                   # tests + lint

# 1. RAPTOR (environment raptor-build; paths from .env.raptor)
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json
pixi run -e raptor-build raptor-generate         # run jobs, deliver data/raptor/res-*/
pixi run -e raptor-build raptor-merge --plan plan.json --job-root RUN_A --job-root RUN_B --output-root MERGED

# 2. EHT imaging (environment eht-rebuild; uvfits locations from .env)
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC    # 4 days x 1000 D draws
pixi run -e eht-research eht-summarize --dataset-root MC_COPY --output-root ENSEMBLE

# 3. Analysis (default environment)
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

```text
experiment.yaml        the parameter file
eht_axion/
  experiment.py        load/validate experiment.yaml, derived settings
  producer/            RAPTOR build, task plan (job ids), job runner, delivery of Stokes cubes
  pipeline/            stage manifests (common), layout, unit constants, RAPTOR image
                       preprocessing / response fit / templates / background (raptor),
                       template reader (products)
  research/            eht-imaging helpers (eht_official, eht_optimizer), EHT ensemble
                       summaries (ensemble, polarimetry), observed likelihood (observed,
                       ensemble_limits)
  model.py limits.py statistics.py observations.py polarization.py angles.py stokes.py
scripts/               raptor_*, eht_*, analysis_* entry points (process order), verify_products
tests/                 synthetic unit/integration tests
reports/               versioned evidence (reports/cleanup-20261106: equivalence of this layout)
docs/                  historical design notes (written for the earlier configs/ layout)
```

## Original project safety

The original project (notebooks, simulations, archives) lives separately in
`/home/canoziia/MyProject/PHYS/EHT_ZC`; this worktree contains none of its data.
