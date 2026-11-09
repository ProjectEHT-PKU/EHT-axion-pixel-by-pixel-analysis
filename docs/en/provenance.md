# Provenance

## Purpose

This page explains how every scientific product records where it came from: stage
directories with a `manifest.json` that binds parameters, input files, output files, code
and software environment by SHA-256, the atomic way these directories are created, and how
to verify a stage and its whole lineage. It also describes the two places that use their
own receipts instead (the EHT bundle and the D-term Monte Carlo batch).

## Inputs

- Any completed stage directory, for verification.
- For writing: the stage name, a JSON-serializable parameter record (`config`) and the list
  of input files, supplied by the calling script.

## Outputs

- `manifest.json` in every stage directory.
- During creation only: `<stage>.lock` next to the destination and a temporary directory
  `<stage>.building-*`, both removed when the call ends.

## How to run

```bash
pixi run verify PATH/TO/STAGE              # check this stage's outputs
pixi run verify PATH/TO/STAGE --recursive  # also re-hash every input and upstream stage
```

`scripts/verify_products.py` prints the manifest after a successful check; it never
recomputes data.

## Details

### Stage directories

`stage_directory(destination, stage, config, input_paths, *, verify_lineage=True)`
(`eht_axion/provenance.py`) is a context manager used by every script that writes
scientific products:

- It refuses destinations reached through a directory symlink and destinations that
  contain one of their own inputs.
- It creates `<destination>.lock` exclusively; a second writer fails instead of waiting.
- If the destination already exists, it is verified recursively and its signature is
  compared with the new one: equal means a cache hit (the body is skipped, the context
  yields `None`), different is an error. An existing stage is never replaced.
- Otherwise the body writes into a temporary directory in the same parent. When the body
  finishes, every file is hashed, `manifest.json` is written last, and the temporary
  directory is renamed to the destination. On failure the temporary directory is removed.

### Manifest

| Field | Content |
|---|---|
| `stage` | stage name, for example `raptor:fit` |
| `schema` | manifest schema version (1) |
| `config` | the parameter record of the stage |
| `inputs` | for each input file: absolute `path`, `size`, `sha256` |
| `implementation` | `git_commit` of the repository and `python_source_sha256`, a hash over the relative paths and bytes of every `.py` file under `eht_axion/` |
| `runtime` | Python version and implementation, platform, machine, versions of numpy, scipy, h5py, matplotlib, PyYAML, and the SHA-256 of `pixi.lock` |
| `signature` | SHA-256 of `stage`, `schema`, `config`, `inputs`, `python_source_sha256` and `runtime` |
| `complete` | `true` |
| `outputs` | for each output file: relative `path`, `size`, `sha256`; shape and dtype for `.npy`, shapes and dtypes of every array for `.npz` |
| `elapsed_seconds` | wall time of the body |

The git commit is recorded but not part of the signature; uncommitted edits are covered by
`python_source_sha256`. Scripts add their own source, `experiment.yaml`, `pixi.lock` and,
for EHT stages, the ehtim package files to the inputs where they matter.

### Verification

`verify_stage(directory, recursive=False)` checks that the manifest is complete and that
every declared output exists inside the directory with the recorded size and SHA-256. With
`recursive=True` it also re-hashes every input and, for inputs named `manifest.json` that
belong to a completed stage, verifies that stage recursively; each upstream stage is
visited once.

`require_declared_outputs(directory, paths)` is used whenever a script reads a file of
another stage: the file must be a declared output of that stage, not an arbitrary extra
file.

`verification_session()` hashes each input file once per invocation and checks at the end
that no hashed file changed (device, inode, size, modification and change times).

By default `stage_directory` verifies the full lineage of every input stage before it
starts. The small result steps `analysis_check_grid.py`, `analysis_check_covariance.py`
and `analysis_summarize.py` pass `verify_lineage=False`: they verify the outputs of the
stages they read, and record those manifests as hashed inputs, but do not re-hash the
whole raw lineage again.

### Stage names

Stage names are fixed identifiers recorded in manifests and checked by the consumers:

| Stage | Written by |
|---|---|
| `producer-build` | `raptor_build.py` |
| `producer-job` | `raptor_generate.py` (one per task) |
| `producer:delivery` | `raptor_generate.py`, `raptor_merge.py` (one per product) |
| `producer:merged-inventory` | `raptor_merge.py` |
| `research:eht-official:I`, `research:eht-official:fiducial` | `eht_reconstruct.py` |
| `research:eht-convex:fiducial` | `eht_refine.py` |
| `research:eht-fixed-D-ensemble` | `eht_summarize.py` |
| `raptor:preprocess`, `raptor:fit`, `raptor:summarize`, `raptor:background-preprocess` | `analysis_templates.py` |
| `research:eht120-observed-limits` | `analysis_limits.py` |
| `research:eht120-limit-numerics` | `analysis_check_grid.py` |
| `research:eht120-covariance-check` | `analysis_check_covariance.py` |
| `research:eht120-results-summary` | `analysis_summarize.py` |

### Receipts outside stage directories

- `eht_export_bundle.py` writes `bundle.json` instead of a manifest: every copied file with
  SHA-256, size and original path, and the science contract. The source stages are
  verified (full lineage) before copying, and copies are made read-only.
- `eht_sample.py` writes `batch.json` whose `identity` is the SHA-256 of the batch plan
  (bundle snapshot hashes, science values, hashes of the `eht_axion` sources, of
  `eht_sample.py`, `eht_reconstruct.py`, `pixi.toml`, `pixi.lock` and of the installed ehtim,
  numpy, scipy, astropy and threadpoolctl files, thread count). Every attempt writes a `taskreceipt.json`
  with the hashes of its outputs; `finish.json` lists the receipt of the certified attempt
  of every draw. `eht_summarize.py` re-checks these receipts before it builds the ensemble
  stage.
- `raptor_generate.py` and `raptor_build.py` append a request record to
  `requests/<uuid>.json` in the run root and keep a copy of the experiment file as
  `config.yaml`.

## Related pages

- [Documentation index](index.md)
- [Data layout](data-layout.md)
- [RAPTOR production](raptor-production.md)
- [EHT imaging](eht-imaging.md)
- [Analysis](analysis.md)
