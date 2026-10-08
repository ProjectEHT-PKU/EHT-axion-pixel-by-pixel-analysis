# Canonical naming migration

## Scope

This change reorganizes files and numeric schemas only. It does not change physical
formulas, selection, angular statistics, priors, or numerical limits. Original shared
LFS data and the previous Pixi/Nix run directories are unchanged.

## Actual locations

- Code: `/home/canoziia/MyProject/PHYS/eht-script-pipeline`
- Canonical input view: `data/inputs/` within the code worktree.
- Canonical migrated run: `/home/canoziia/MyProject/PHYS/eht-pipeline-runs/research-001`
- Source run: `eht-pipeline-runs/pixi-research` (previously validated Pixi output).
- Original input archive: `EHT_ZC/shared_data/analysis-inputs/draft` (unchanged).

`draft` remains only in the original publication archive, not the canonical runtime
input path. The new view links to original bytes; it is not a standalone backup.

## Naming rules now in code

- Two-decimal alpha identity: `alpha-0.30`.
- ISO EHT date: `2017-04-05`, not day-of-year95.
- Physical directory names instead of numbered audit stages.
- Large images/ensembles remain `.npy`, preserving mmap behavior.
- Response NPZ contains `amplitude`, `phase`, `baseline`, `snapshot_labels`.
- EHT daily summary NPZ contains named1D/2D mean/std/circular branch/concentration/count,
  intensity mean and1D covariance.
- Final `limits-1d.npz`/`limits-2d.npz` contain `alphas`, `modes`, `cuts` with(mode,alpha)
  axes, plus observation/prior/variant metadata.
- Stage completion marker is `manifest.json`; old `stage.json` remains readable.
- New stage manifests describe array shape/dtype and preserve input/config/code hashes.
- A run-root `manifest.json` and `config.toml` describe configuration/environment, but
  are explicitly not a completion marker. Inspect per-stage manifests.

Full schema: [DATA_LAYOUT_PLAN.md](../../docs/DATA_LAYOUT_PLAN.md).

## Compatibility

- Input readers support the old shared input layout and the new canonical view.
- Response summarization supports old three-column NPY and new named NPZ.
- Observation builders can read old/new daily layouts without changing circular lifts.
- Plot readers support the old six-file layout and the new two-file/mode-axis layout;
  they reject mixed or inconsistent results.
- Existing valid runs are migrated explicitly, never silently relabeled as recomputed.
- Imported products retain original implementation/signature and checksums.
  `--reuse-imported` is required. Changed config/lineage still causes rejection.

## Real-data validation

Evidence: [naming-migration.json](naming-migration.json).

-25 stage checksum manifests verified after migration.
- All16 response NPZs compared field-by-field with original `[A,phase,baseline]` data:
  exact equality, including NaNs.
- All four days'1D/2D mean/std/center/unwrapped_mean/concentration/circular_std/count
  fields compared exactly with original combined summaries.
- All48 final upper limits exactly unchanged (max difference0dex).
- Canonical input view contains all34 source files, all linked after source hash checks.
- Large processed Stokes arrays are symlinks; original datasets were not copied or modified.
- Imported whole-chain replay succeeds with25 explicit imported-reuse messages.
- New-schema synthetic end-to-end computation succeeds, including cache replay,
  diagnostics, covariance, plots, config mismatch checks and two named limit files.
- Complete legacy environment test suite:239 passed. Review fixes cover changed-source
  lineage, explicit imported-upstream permission and directory-symlink rebuild protection.
- Original source run's25 stage input/output lineage chains were independently verified
  before the final migration, not merely rebound to current input hashes.
- Figure creation validates inference's upstream lineage without recomputing inference.
- Root scientific identity is stable; later command/config/environment records append
  to `invocations/` instead of overwriting the run's original identity.

The migrated run adds≈149MiB, mostly converting moderate response arrays into named NPZ.
The canonical input view itself uses≈236KiB. All linked source archives must remain in place.
No expensive raw scientific computation was performed just to rename data.

## Commands

```bash
pixi run migrate inputs --source /path/to/old/analysis-inputs/draft --destination data/inputs
pixi run migrate run --source /path/to/old-run --destination /path/to/research-001

# Inspect/reuse imported results explicitly:
pixi run pipeline --config configs/research.toml --input-root data/inputs \
  --output-root /path/to/research-001 --trust-pickle --reuse-imported

# New actual computation uses new names natively:
pixi run pipeline --config configs/research.toml --input-root data/inputs \
  --output-root /path/to/research-002 --trust-pickle
```

Do not point a migration destination inside its source tree. A partial multi-stage run
migration is not a completed run: inspect the migration log and choose a new destination
on retry; individual completed stages remain checksummed, original data remain untouched.

The transport variant identifiers `jonly/eq_jonly/fake` have deliberately not been
expanded to new physical names: their exact original switch settings are not verified.
