# Data naming migration contract

Implemented on the script-pipeline branch. Validation and actual paths are in
[the naming-migration report](../reports/script-pipeline/NAMING_MIGRATION.md).

This change reorganizes data identity/schema, not physical formulas. Old data and
validated runs remain read-only; importing existing arrays must not copy large NPYs.

## Canonical inputs

```text
data/inputs/
  raptor/res-100/alpha-0.10/stokes.npy
  raptor/res-400/alpha-0.30/stokes-t1600-t1690.npy
  eht/2017-04-05/image-ensemble.pkl
  transport/snapshot-2110-res-400/*.npy
  manifest.json
```

`draft` is only a legacy publication namespace, not an analysis requirement.
Alpha identifiers use two decimals; dates are ISO8601. The adapter records symlink
sources, hashes, units/axis uncertainties without pretending unknown Stokes units
are known. No physical data are renamed inside EHT_ZC.

## Canonical run

```text
config.yaml
manifest.json                       stable run identity (not stage completion)
invocations/*.json                  append-only requested config/environment records
raptor/processed/alpha-0.30/stokes.npy
raptor/response/alpha-0.30/response-2d.npz   amplitude, phase, baseline, snapshot_labels
raptor/response/alpha-0.30/response-1d.npz
raptor/templates/1d/templates.npz     existing named numeric keys
raptor/templates/2d/templates.npz
raptor/background/stokes.npy          state0 only
background/evpa-difference.npz        histogram+sigma, all named
background/summary.json

eht/selection.npz
eht/mask.npy
eht/2017-04-05/stokes.npy
eht/2017-04-05/evpa-1d.npy
eht/2017-04-05/summary.npz
  evpa_{mean,std,center,unwrapped_mean,concentration,circular_std,count}_{1d,2d}
  intensity_mean, covariance_1d

observations/1d/observations.npz
observations/2d/observations.npz
inference/limits-1d.npz              alphas, modes, cuts (mode,alpha), metadata
inference/limits-2d.npz
inference/posterior-{mode}-{dim}.json
inference/run.json
figures/limits.{pdf,png,svg}, limits.csv
diagnostics/...
```

Every completed stage uses `manifest.json` (formerly stage.json). `raptor/templates`
and `observations` retain a dimension subdirectory because each dimension is an
independent validated cache unit. Physical axes and units belong in manifests, not
filenames. `response-*.npz` are moderate products (~15MB combined per alpha); large
Stokes and EVPA ensembles remain mmap-capable NPY.

## Compatibility and migration

**Safety update:** see [engineering safety and cache contract](ENGINEERING_SAFETY.md).
Migration now preflights and publishes an entire run atomically. `origin_inputs`
are independently checksummed live dependencies on original manifests/products,
separate from new-layout `inputs`. Old hash receipts alone are not validated origins.
Imported reuse requires opt-in and preserved source archives; new manifests bind
necessary upstream stage manifests and the runtime/lock identity. Old caches are
not silently re-signed to meet this stronger contract.

- Input resolution accepts canonical layout or original shared analysis-inputs.
- An explicit migration CLI creates a new canonical view from an existing complete
  run; NPY bytes are shared by read-only symlink by default. Moderate fit arrays become
  named NPZ with exact element checks.
- Migration records old manifest + input/output hashes, actual old source identity and
  a new migration code identity. Never relabel old science as newly recomputed.
- New stages verify imported products but do not call them a current-source cache hit.
  `--reuse-imported` is explicitly required; config/input changes still fail. Their original
  scientific source identity is retained separately from the migration code identity.
- Original run directories, original LFS manifests and binary objects stay untouched.
