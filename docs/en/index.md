# Documentation index

## Purpose

This is the entry page of the English documentation of the M87\* EHT axion analysis. The
repository uses two independent inputs: RAPTOR ray-traced images of a GRMHD flow with an
axion field, and EHT 2017 polarimetric images of M87\*. The analysis combines them into
Bayesian upper limits on a dimensionless coupling $c$. Every page below describes the
code as it is in this repository, and every page links back here. A Chinese translation is
in [中文文档](../zh/index.md).

## Pages

| Page | Content |
|---|---|
| [Overview](overview.md) | The three parts, how their products connect, and a data-flow diagram |
| [Experiment configuration](experiment-config.md) | Every key of `experiment.yaml`: meaning, unit, value and the code that reads it |
| [RAPTOR production](raptor-production.md) | Task plan and task identifiers, build, jobs, job outputs, conversion, delivery and merge, CPU affinity, `.env.raptor` |
| [EHT imaging](eht-imaging.md) | Paper VII method as implemented: I stage, shared gains, fiducial P/D alternation, refinement, bundle export, D-term Monte Carlo, ensemble summary, known deviations |
| [Analysis](analysis.md) | Image preprocessing, response fit, templates, aperture observable, feature selection, GRMHD background, covariance models, likelihood, posterior, checks and summary |
| [Provenance](provenance.md) | Stage manifests, hashes, atomic stage directories, verification |
| [Data layout](data-layout.md) | Directories and files that each script writes, with array shapes and units |
| [Package map](package-map.md) | Every module of `eht_axion` and every script, with a one-line purpose |
| [Development](development.md) | Pixi environments and tasks, tests, lint, `.env` files |
| [Glossary](glossary.md) | Symbols and terms used on these pages |

## Where to start

- To understand the experiment: [Overview](overview.md), then [Experiment configuration](experiment-config.md).
- To audit a scientific step: [EHT imaging](eht-imaging.md) or [Analysis](analysis.md),
  together with [Provenance](provenance.md).
- To run the scripts: the "How to run" sections of [RAPTOR production](raptor-production.md),
  [EHT imaging](eht-imaging.md) and [Analysis](analysis.md); environments are described in
  [Development](development.md).
- To find a function or file: [Package map](package-map.md) and [Data layout](data-layout.md).

## Conventions on these pages

- File, key, function and stage names are written as code spans, for example `analysis.features`.
- YAML keys are given with their section path, for example `eht_imaging.refine.gap_atol`.
- Angles of electric-vector position angle (EVPA) are in radians in all arrays; EVPA is
  defined modulo $\pi$. Sky angles in EHT products are measured east of north.
- Lengths on the sky are in microarcseconds (µas) or in gravitational radii $r_g = GM/c^2$;
  times are in hours or in $t_g = GM/c^3$. The light speed $c$ in these definitions is not
  the coupling $c$; see the [Glossary](glossary.md).

## Related pages

- [Repository README](../../README.md)
- [中文文档](../zh/index.md)
- [Overview](overview.md)
- [Glossary](glossary.md)
