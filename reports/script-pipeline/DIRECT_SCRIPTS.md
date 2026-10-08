# Direct research scripts, reusable modules, no installable project

## What changed

- Removed `[build-system]`, `[project]`, package discovery and console registration
  from `pyproject.toml`; it now contains pytest/Ruff settings only.
- Removed editable project installation and explicit pip/setuptools installation
  requirements from Pixi; lock regenerated without a local PyPI distribution.
- Removed `eht_axion/pipeline/cli.py` and package `__main__.py` entrypoint.
- Moved research orchestration into `workflow.py::run_stage`, an explicit-path Python
  function without command-line parsing or installed entrypoints.
- Pixi tasks now execute individual `scripts/*.py` files directly.
- Moved argparse/main code out of RAPTOR/EHT/transport/plot/migration/convergence
  modules; reusable computation modules contain numerical functions, not CLI code.
- A small `scripts/_bootstrap.py` locates adjacent source modules from the script's
  own path. No global installation or manual PYTHONPATH required.
- Old editable-install metadata and the generated console command were removed.

The local `eht_axion` namespace still organizes reusable angle, HDF5/image, plotting,
statistics and cache functions. It is intentionally importable from the source tree
and scripts, but is not installed globally or presented as a universal M87 pipeline.

## Validation

- Default Pixi test/lint gate: **255 passed**, Ruff passed.
- All17 direct scripts launch `--help` from outside the repository with PYTHONPATH
  removed, using the Pixi Python interpreter.
- Synthetic end-to-end scientific test invokes actual `scripts/run_all.py` from
  an external working directory, produces six curves/figures and checks caches.
- Low-level HDF5/EHT/transport script parameter/default/error tests preserved.
- Outside source directory, `importlib.metadata` finds no `eht-axion-analysis`
  distribution; no `eht-pipeline` executable is installed; source namespace is not
  globally importable. Direct absolute script paths work as intended.
- `pixi install --locked` passes with updated lockfile.
- AST comparison of all59 non-main functions/classes in the7 affected numerical
  modules confirms function bodies are unchanged. Evidence:
  [direct-script-kernel-check.json](direct-script-kernel-check.json).

No expensive science calculation was rerun solely for this change. Existing arrays,
limits and figures are unchanged; the boundary between scripts and reusable modules
was restructured. Code-signature changes may require a new output root or explicit
imported-product reuse under the existing cache policy.

## How to run now

```bash
pixi install --locked
pixi run check
pixi run raptor --config configs/research.toml --input-root data/inputs --output-root ../eht-pipeline-runs/my-run
pixi run pipeline --config configs/research.toml --input-root data/inputs --output-root ../eht-pipeline-runs/my-run --trust-pickle
```

`pixi run raptor` executes `python scripts/process_raptor.py`.
`pixi run pipeline` executes `python scripts/run_all.py`; it is just a convenience
alias, not a registered CLI. See [REUSABLE_MODULES.md](../../docs/REUSABLE_MODULES.md)
for source-function examples and the full script mapping.
