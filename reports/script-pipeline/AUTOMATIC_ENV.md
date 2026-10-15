# Automatic per-stage dotenv loading

User superseded the previous manual export/source requirement. Direct analysis
entrypoints now parse repository `.env` automatically. Producer build/generate
entrypoints parse `.env` then `.env.raptor`. CLI overrides exported environment,
which overrides file values; producer overlay overrides shared file values.

Implementation uses python-dotenv (Conda dependency in pixi.toml/lock), not a
shell evaluator or handwritten dotenv syntax. `interpolate=False` prevents
variable expansion. Only declared EHT keys are accepted, malformed/bare/duplicate
keys rejected. No arbitrary cwd search, no shell execution, and no mutation of
os.environ. File content is read once then parsed/validated from the same bytes.
File-derived relative paths anchor to the checkout root; CLI/externally exported
relative paths retain cwd semantics. Explicit empty exported values still override
file defaults. This does not silently grant pickle trust or enable rebuild.

`--env-file PATH` replaces default file selection, `--no-env-file` disables file
loading without clearing existing shell variables. Help is processed before any
file read. Low-level expert tools with explicitly required arguments remain
explicit; the shared stage/producer scripts handle the env-backed workflows.

Tracked `.env.example` contains analysis/shared defaults only. New tracked
`.env.raptor.example` contains producer operational paths/run ID/threads. A local
ignored `.env.raptor` was created from that reviewed example; existing local `.env`
values are unchanged (comments updated). Both local files remain untracked.
New deployment copies/reviews templates once, no export/source step afterwards.

Verified with every EHT_* export removed: direct `pixi run -e raptor-build
generate-raptor --plan` finds the copied fluid/radial data and expands 18,330
tasks without executing C, building, or creating a new production run. Analysis
parse_options also resolves the existing .env automatically, without running
any scientific stage.

Validation: **1015 tests passed in131.06s**, including root-only discovery,
producer overlay isolation, CLI/export precedence, quoting/spaces, no variable or
command evaluation, malformed files, explicit file replacement/disable switches,
reload without process contamination, and producer entrypoint parsing. Ruff and
git diff --check passed. python-dotenv is available in default and raptor-build
environments. Logs in `automatic-env-validation/`.

No scientific parameters, datasets, caches or old results changed. No real
snapshot, build or inference was run. Adding the dependency changes the lock
identity, so old caches can be rejected; they were not rewritten or recomputed.
