# RAPTOR + EHT-Imaging analysis inputs — development draft

> Mutable development dataset; not a frozen paper release.

## Contents

### RAPTOR Stokes cubes

```text
raptor-stokes/res-100/alpha-*.npy
```

Eight arrays with shape `(141, 17, 100, 100, 4)`:

- 141 GRMHD times from 1600 through 3000 in steps of 10;
- state 0 without axion, states 1–16 with 16 initial phases;
- Stokes axis `[I, Q, U, V]`;
- float64 values.

The high-resolution `(141, 17, 400, 400, 4)` alpha=0.3 source is stored as
independently loadable time-axis shards:

```text
raptor-stokes/res-400/alpha-0.3/stokes-t1600-t1690.npy
...
raptor-stokes/res-400/alpha-0.3/stokes-t3000-t3000.npy
```

Concatenate shards in filename/time order along axis 0 to reconstruct the source.

### EHT image ensembles

```text
eht-image-ensemble/230_{95,96,100,101}.pkl
```

Each contains 4800 image/hyperparameter combinations with `I_data`, `Q_data`,
`U_data`, historical `theta_data`, and `meta_data`. EVPA should be recomputed
from Q/U; the stored historical theta array is zero-filled.

## Verify

```bash
python scripts/verify_shared_data.py shared_data/analysis-inputs/draft
```

`manifest.json` is authoritative for checksums, exact shapes, dtypes, source
indices, time values, metadata fields, and source provenance.
