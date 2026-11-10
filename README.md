# M87* EHT axion analysis — script pipeline

Scripts for two independent inputs — RAPTOR ray-traced images of a GRMHD flow with an
axion field, and EHT 2017 M87* polarimetric images — and the analysis that combines them
into phase-marginalized Bayesian upper limits on a dimensionless coupling c. This is not
an installable package: run the scripts directly or through `pixi run`.

**Documentation: [docs/en/index.md](docs/en/index.md)** (overview, every key of
`experiment.yaml`, RAPTOR simulation, EHT imaging, analysis, provenance, data layout,
package map, development). Chinese translation: [docs/zh/index.md](docs/zh/index.md).

## Parameters

[`experiment.yaml`](experiment.yaml) holds every scientific value, each exactly once, in
the sections `shared`, `raptor`, `eht_imaging` and `analysis`. Each script reads `shared`
plus its own section; the analysis also reads the `raptor` and `eht_imaging` values it
needs. Machine paths and execution options are command-line options or `.env` /
`.env.raptor` entries (see [`.env.example`](.env.example) and
[`.env.raptor.example`](.env.raptor.example)). Details:
[Experiment configuration](docs/en/experiment-config.md).

## Run

RAPTOR simulation and EHT imaging do not depend on each other; the analysis needs the
products of both. With a data root (`M87_AXION_DATA_ROOT` or `--data-root`), every script
reads and writes the default layout of [Data layout](docs/en/data-layout.md#default-layout);
explicit `--...-root` options override it. Every output root must be new.

```bash
pixi install --locked
pixi run check                                   # tests + lint
export M87_AXION_DATA_ROOT=/path/to/m87-axion    # or set it in .env

# RAPTOR simulation (environment raptor-build)
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate --plan --plan-output $M87_AXION_DATA_ROOT/raptor/plan.json
pixi run -e raptor-build raptor-generate
pixi run -e raptor-build raptor-merge

# EHT imaging (environment eht-imaging; each of the days 95, 96, 100, 101)
pixi run -e eht-imaging eht-reconstruct --stage I --day 95 --review-file REVIEW
pixi run -e eht-imaging eht-reconstruct --stage fiducial --day 95 --review-file REVIEW
pixi run -e eht-imaging eht-refine --day 95
pixi run -e eht-imaging eht-export-bundle
pixi run -e eht-imaging eht-sample
pixi run -e eht-summary eht-summarize

# Analysis (default environment)
pixi run analysis-templates
pixi run analysis-limits
pixi run analysis-check-grid --check grid --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --check prior --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance
pixi run analysis-summarize
pixi run analysis-paper                          # figures, numbers and tables of the paper
```

Full command lines, including merging RAPTOR runs and the analysis checks, are in
[RAPTOR simulation](docs/en/raptor-simulation.md), [EHT imaging](docs/en/eht-imaging.md)
and [Analysis](docs/en/analysis.md). Each step writes stage directories with a
`manifest.json`; `pixi run verify DIR --recursive` re-checks a stage and its lineage
([Provenance](docs/en/provenance.md)).

## Layout

```text
experiment.yaml   the parameter file
eht_axion/        library: experiment, provenance, physics, evpa; raptor/, eht/, analysis/
scripts/          raptor_*, eht_*, analysis_* entry points, verify_products
tests/            synthetic unit and integration tests
docs/en/          documentation
```

A module-by-module list is in [Package map](docs/en/package-map.md), and the default data
directory layout is in [Data layout](docs/en/data-layout.md).
