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

### RAPTOR transport-term diagnostics (Figure 4)

```text
raptor-transport-diagnostics/snapshot-2110-res-400/
├── common.npy
├── common_noaxion.npy
├── jonly.npy
├── eq_jonly.npy
├── fake.npy
├── r_coord_400.npy
└── axion_field_coord_400_16.npy
```

These single-snapshot runs compare full/normal transfer, emissivity-only
transport, equatorial emissivity, and an artificial emissivity reference. They
are required for the plasma-on/off and source-contribution decomposition shown
in Figure 4. The exact upstream RAPTOR parameter files are not present, so the
variant interpretations in the manifest are inferred from `GRMHD/main.ipynb`
and must be confirmed before paper release.

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
