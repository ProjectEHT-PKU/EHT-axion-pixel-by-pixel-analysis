# Analysis-input datasets

These datasets start after the two expensive upstream applications have run:

- RAPTOR has produced Stokes I/Q/U/V images;
- EHT-Imaging has produced 4800 image candidates for each of four days.

They support reconstruction of image processing, topset selection, EVPA signal
fits, covariance matrices, and final limits. They do **not** contain raw GRMHD
fluid snapshots, RAPTOR source/configuration, EHT visibility/UVFITS data, or the
EHT-Imaging reconstruction setup, and therefore are not raw end-to-end inputs.

## Lifecycle

- `draft/` is mutable while analysis changes are under review.
- `vN/` is an immutable paper snapshot and must never be overwritten.
- Every update must include a matching `manifest.json` with SHA-256, shapes,
  dtypes, physical index meanings, and provenance.

## Download and verify

```bash
git lfs install
git lfs pull --include='shared_data/analysis-inputs/draft/**'
python scripts/verify_shared_data.py shared_data/analysis-inputs/draft
```

The full draft is approximately 20.9 GB. Download subsets when possible:

```bash
# Only EHT four-day 4800-image ensembles (~2.52 GB)
git lfs pull --include='shared_data/analysis-inputs/draft/eht-image-ensemble/**'

# Only RAPTOR res=100 alpha scan (~6.14 GB)
git lfs pull --include='shared_data/analysis-inputs/draft/raptor-stokes/res-100/**'

# Only RAPTOR alpha=.3 res=400 shards (~12.27 GB)
git lfs pull --include='shared_data/analysis-inputs/draft/raptor-stokes/res-400/**'
```

## Update draft

```bash
python scripts/publish_analysis_inputs.py draft --force-draft
python scripts/verify_shared_data.py shared_data/analysis-inputs/draft
git add shared_data/analysis-inputs/draft
```

The res=400 source is sharded by ten time slices. Every shard is independently
loadable with NumPy and below 1 GB. Unchanged LFS content is reused, but every
changed shard creates a new complete LFS object.
