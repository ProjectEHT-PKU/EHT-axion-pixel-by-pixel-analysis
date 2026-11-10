# Provenance

## Purpose

This page explains how every scientific product records where it came from: stage
directories with a `manifest.json` that records the parameters, the git commit of the code
and the SHA-256 of every data file read and written, the atomic way these directories are created, and how
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

`scripts/verify_products.py` runs the check in one verification session (each file is
hashed once, however many stages declare it), prints the manifest after a successful
check, and never recomputes data.

## Details

### Stage directories

`stage_directory(destination, stage, config, input_paths, *, verify_lineage=True)`
(`eht_axion/provenance.py`) is a context manager used by every script that writes
scientific products:

- It refuses destinations reached through a directory symlink and destinations that
  contain one of their own inputs.
- It creates `<destination>.lock` exclusively; a second writer fails instead of waiting.
- If the destination already exists, it is verified (see [Output roots](#output-roots))
  and its recorded stage name, `config` and `inputs` are compared with the new ones: equal
  means the stage is reused (the body is skipped, the context yields `None`), different is
  an error. The code commit is not compared, so a resumable RAPTOR run can continue after
  a commit that does not change its tasks. An existing stage is never replaced.
- A new stage is written only if the tracked files of the repository are unchanged
  (`code_identity`; untracked files do not matter); otherwise it is refused.
- Otherwise the body writes into a temporary directory in the same parent. When the body
  finishes, every file is hashed, `manifest.json` is written last, and the temporary
  directory is renamed to the destination. On failure the temporary directory is removed.

### Manifest

| Field | Content |
|---|---|
| `stage` | stage name, for example `analysis:fit` |
| `schema` | manifest schema version (1) |
| `config` | the parameter record of the stage |
| `inputs` | for each input file: absolute `path`, `size`, `sha256` |
| `implementation` | `git_commit`: the full commit of the repository |
| `runtime` | Python version and implementation, platform, machine, and the version strings of numpy, scipy, h5py, matplotlib and PyYAML (not hashed) |
| `complete` | `true` |
| `outputs` | for each output file: relative `path`, `size`, `sha256`; shape and dtype for `.npy`, shapes and dtypes of every array for `.npz` |
| `elapsed_seconds` | wall time of the body |

The code is identified by the git commit only. Because a stage can be written only from a
clean tracked tree, the commit pins every script, every `eht_axion` module,
`configs/experiment.yaml`, `pixi.toml` and `pixi.lock`; none of these files is hashed or listed
as an input. Third-party packages are pinned by `pixi.lock`; stages that depend on a
specific package record its version (for example `ehtim_version` in the EHT imaging plans).
Files that git does not hold are hashed: data inputs read from outside the repository
(uvfits, GRMHD snapshots, radial tables, the official EHT script, upstream stage files, the
built RAPTOR binary) and every stage output.

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
starts. The small result steps `analysis_check_grid.py`, `analysis_check_covariance.py`,
`analysis_summarize.py` and `analysis_paper.py` pass `verify_lineage=False`: they verify the outputs of the
stages they read, and record those manifests as hashed inputs, but do not re-hash the
whole raw lineage again.

### Output roots

One rule applies to every script: a script that produces a result refuses an output root
that already exists. The exceptions are the resumable parts of the pipeline:
`raptor_build.py`, `raptor_generate.py` and `raptor_merge.py` reuse completed stages in
their run root (a stage with the same stage name, config and input data is reused; a
different stage at the same place is an error), and `eht_sample.py` resumes a batch in the same output root if the
batch identity is unchanged. `analysis_templates.py` writes several stages below one new
root.

When `stage_directory` finds an existing destination it verifies it with the same lineage
setting as for new stages: recursively by default, outputs only with
`verify_lineage=False`.

### Stage names

Stage names are fixed identifiers recorded in manifests and checked by the consumers. They
have the form `<part>:<step>` with the part `raptor`, `eht` or `analysis`. JSON records
outside manifests carry schema strings of the form `<part>-<record>-v<N>`
(`raptor-generated-v2`, `raptor-merge-v2`, `eht-bundle-v2`, `eht-sample-batch-v2`,
`eht-sample-task-v2`, `eht-sample-finish-v2`, `eht-ensemble-v2`, `analysis-limits-v2`,
`analysis-summary-v2`, `analysis-paper-v1`, `analysis-paper-numbers-v1`).

| Stage | Written by |
|---|---|
| `raptor:build` | `raptor_build.py` |
| `raptor:job` | `raptor_generate.py` (one per task) |
| `raptor:delivery` | `raptor_generate.py`, `raptor_merge.py` (one per product) |
| `raptor:merge` | `raptor_merge.py` |
| `eht:i`, `eht:fiducial` | `eht_reconstruct.py` |
| `eht:refine` | `eht_refine.py` |
| `eht:ensemble` | `eht_summarize.py` |
| `analysis:preprocess`, `analysis:fit`, `analysis:templates`, `analysis:background` | `analysis_templates.py` |
| `analysis:limits` | `analysis_limits.py` |
| `analysis:check-grid` | `analysis_check_grid.py` |
| `analysis:check-covariance` | `analysis_check_covariance.py` |
| `analysis:summary` | `analysis_summarize.py` |
| `analysis:paper` | `analysis_paper.py` |

### Receipts outside stage directories

- `eht_export_bundle.py` writes `bundle.json` instead of a manifest: every copied file with
  SHA-256, size and original path, the science contract and the git commit (`code`). The source stages are
  verified (full lineage) before copying, and copies are made read-only.
- `eht_sample.py` writes `batch.json` whose `identity` is the SHA-256 of the batch plan
  (bundle snapshot hashes, science values, the git commit with the versions of ehtim,
  numpy, scipy, astropy and threadpoolctl, thread count). Every attempt writes a `taskreceipt.json`
  with the hashes of its outputs; `finish.json` lists the receipt of the certified attempt
  of every draw. `eht_summarize.py` re-checks these receipts before it builds the ensemble
  stage.
- `raptor_generate.py` and `raptor_build.py` append a request record to
  `requests/<uuid>.json` in the run root and keep a copy of the experiment file as
  `config.yaml`.

## Related pages

- [Documentation index](index.md)
- [Data layout](data-layout.md)
- [RAPTOR simulation](raptor-simulation.md)
- [EHT imaging](eht-imaging.md)
- [Analysis](analysis.md)
