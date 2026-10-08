# Pixi portability validation

## Changes

- Added `pixi.toml` and committed multi-platform `pixi.lock`.
- Python3.12 and numeric/HDF5/plot/test dependencies managed by Pixi.
- Editable package installation, no Nix/system Python/manual PYTHONPATH required.
-17 tasks for scientific stages, integrity checks, comparison, tests and lint.
- Default environment excludes statsmodels; optional `legacy` adds historical utilities.
- Agg rendering, two BLAS/OpenMP threads enabled in activation.
- `.pixi/` remains ignored; manifests and lockfile explicitly allowed by Git ignore rules.
- Main README and guide use Pixi; detailed Chinese guide in `docs/PIXI_QUICKSTART.md`.

## Execution environment

Actually tested on Linux x86-64, Pixi0.81.0. Package versions and lock checksum:
[pixi-environment.json](pixi-environment.json).

Lock resolves Linux x86-64/ARM64,macOS Intel/Apple Silicon,Windows x86-64.
Other platforms were not execution-tested. Current lock baseline: glibc2.28,
macOS13,Windows10. Recommend a current Pixi version supporting lockfile format7
(the validation used0.81.0).

## Checks completed

- `pixi install --no-config` created the lock and default environment.
- `pixi run --locked help`, task argument forwarding and console/module loading work.
- Invocation from `/tmp` using `--manifest-path` works.
- Default:204 tests passed,1 optional legacy module skipped; Ruff passed.
- `pixi install --locked -e legacy` succeeds; legacy:209 tests passed with no skips.
- `pixi run --locked verify` successfully validates the earlier Nix-generated figure stage.
- No absolute `/home` or `/nix` paths embedded in pixi.lock; editable project uses `./`.

Recorded console output: [pixi-tests.txt](pixi-tests.txt).

## Real-data run

A separate real-data run uses the Pixi default environment and unchanged scientific
configuration. It starts at all8 raw RAPTOR res100 cubes, res400 time shards for
state0 background, and all four4800-image EHT ensembles. The previous Nix-generated
results remain untouched.

```bash
pixi run --locked pipeline --config configs/research.toml \
  --input-root /home/canoziia/MyProject/PHYS/EHT_ZC/shared_data/analysis-inputs/draft \
  --output-root /home/canoziia/MyProject/PHYS/eht-pipeline-runs/pixi-research \
  --trust-pickle
```

**The real-data run completed successfully.** All25 stage manifests/output checksums
were verified. A full cache replay reported25 cache hits. Six curves agree with the
previous Nix run to a maximum absolute difference of **3.34e-15dex** over all48 points.
Evidence: [pixi-result-comparison.json](pixi-result-comparison.json).

Exact reproducibility across OS/BLAS/library versions is not assumed; numerical
differences were measured. Scientific caveats from the research profile remain
unchanged by switching environment managers. Matplotlib logged a nonfatal semibold
font-weight fallback to700; all PNG/PDF/SVG files were generated.

## Disk use

Default environment≈1.5GiB; default+legacy≈3GiB on this machine. `.pixi/` is local,
not tracked. The separate validation output adds≈3.53GB of physical cache data;
ordinary collaborators can choose one environment and one scientific output directory.
