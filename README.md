# M87* EHT axion analysis — script pipeline

Scripts for two independent inputs — RAPTOR ray-traced images of a GRMHD flow with an
axion field, and EHT 2017 M87* polarimetric images — and the analysis that combines them
into phase-marginalized Bayesian upper limits on a dimensionless coupling c. This is not
an installable package: run the scripts directly or through `pixi run`.

**Documentation: [docs/en/index.md](docs/en/index.md)** (overview, every key of
`experiment.yaml`, RAPTOR production, EHT imaging, analysis, provenance, data layout,
package map, development).

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
products of both. Every output root must be new.

```bash
pixi install --locked
pixi run check                                   # tests + lint

# RAPTOR simulation (environment raptor-build)
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate

# EHT imaging (environment eht-rebuild; each of the days 95, 96, 100, 101)
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-refine --day 95 --source-root PD95 --i-root I95 --output-root D95
pixi run -e eht-rebuild eht-export-bundle --i-root I95 ... --d-root D95 ... --output-root BUNDLE
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC
pixi run -e eht-research eht-summarize --sample-root MC --output-root ENSEMBLE

# Analysis (default environment)
pixi run analysis-templates --raptor-input-root RUN_ROOT/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
```

Full command lines, including merging RAPTOR runs and the analysis checks, are in
[RAPTOR production](docs/en/raptor-production.md), [EHT imaging](docs/en/eht-imaging.md)
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
