# Correct principal N=4 to N=2 for the article's scalar 211 mode

## Authorized scientific change

The article explicitly defines n as the principal quantum number and selects
(n,l,m)=(2,1,1) throughout (main.tex:151,162,179). The old Notebook used
`alpha*(1-alpha**2/2/(2+1+1)**2)`. Its leading frequency corresponds to N=4,
not the article's 211. In the radial-overtone convention N=nr+l+1, 211 means
nr=0 and l=1, hence N=2. m does not enter the leading denominator.

The current script pipeline now uses N=2 in configs/research.yaml,
model.omega(), model.evpa_templates() and LimitOptions defaults. Explicit
other principal N remain available. Historical Notebook regression tests and
scripts explicitly request N=4 so they remain historical, not new defaults.
The old EHT_ZC repository, notebooks, data, templates, caches and saved results
are untouched. No science rerun or cache re-signing was performed. Existing
inference products have their old identities; do not relabel them N=2.

The updated formula returns dimensionless angular frequency:

    M omega_R = alpha*(1-alpha**2/8)

This is the leading small-alpha 211 approximation, not an exact Kerr spectrum.
The original RAPTOR cubes used time*alpha; this change does not retroactively
correct their integrated ray phases or the BL/KS inconsistency. The new source
accepts independent alpha and omega; its future runner must supply a consistent
frequency. No radial reference table is changed by this correction.

## Comparison with the radial generator's external real-frequency tables

`source/a=..._l1m1_real.txt` column 2 is omega_R/mu. Use the same linear
interpolation rule as get_omega(), then multiply by alpha. Low five alpha nodes
use cloud spin .9375; high three use .99, as in the existing selected profiles.
These are interpolated external source tables, not certified exact eigenvalues.
README explicitly describes spin interpolation for the .9375 source.

Relative difference below is `(omega_N2/omega_table-1)*100`, not binding-energy
error. All eight leading N2 frequencies are higher, so their predicted periods
are slightly shorter.

| alpha | cloud spin | M omega_N2 | M omega_table | frequency difference |
|---:|---:|---:|---:|---:|
| .10 | .9375 | .099875000000 | .099854704950 | +.020325% |
| .15 | .9375 | .149578125000 | .149552296783 | +.017270% |
| .20 | .9375 | .199000000000 | .198970961663 | +.014594% |
| .25 | .9375 | .248046875000 | .247938794413 | +.043592% |
| .30 | .9375 | .296625000000 | .296319073264 | +.103242% |
| .35 | .99 | .344640625000 | .343855826447 | +.228235% |
| .40 | .99 | .392000000000 | .390415500000 | +.405850% |
| .44 | .99 | .429352000000 | .426552036870 | +.656418% |

At alpha=.30 the absolute difference in M omega is .000305926736;
at alpha=.44 it is .002799963130. No limit-shift estimate is inferred from
these frequency differences. Finite exposure and separation phase effects must
be evaluated by the inference code if a rerun is requested.

Reproducible read-only comparison:

```bash
pixi run python reports/script-pipeline/compare_211_frequencies.py \
  --radial-root ../radial-function --output-dir /tmp/new-211-comparison
```

Machine-readable numbers and input hashes: `frequency-211/comparison.json` and
`comparison.csv`. One-source convention research: `frequency-211/CONVENTIONS.md`.
The research report includes Dolan Eq.(32), Detweiler Eq.(26), Baumann Eq.(2.28),
principal vs radial n definitions and higher-order limitations.

## Validation

New regression coverage checks all eight default frequencies, consistent YAML,
LimitOptions and function defaults, explicit N=4 historical formula, and
implicit-vs-explicit N=2 template equality. The independent historical
Notebook posterior oracle remains explicitly N=4. Full test/lint results are
recorded alongside this report, without touching real outputs:

- `pixi run test`: **752 passed in 82.09 seconds**.
- `pixi run lint`: passed.
- Standalone comparison-script Ruff check and `git diff --check`: passed.
- The pre-N2 whole-chain numerical fixture explicitly selects historical N=4;
  it remains an engineering regression, not an oracle for the corrected physics.
- Logs: `frequency-211/tests.log`, `frequency-211/lint.log`.
