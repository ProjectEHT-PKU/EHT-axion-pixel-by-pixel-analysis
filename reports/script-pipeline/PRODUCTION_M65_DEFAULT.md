# New production default: MBH=6.5e9 Msun

User authorized changing the new producer's default mass after reviewing recent
EHT literature. `configs/production.yaml` now uses **6.5e9 Msun** and retains
**D=16.8 Mpc**. The shared profile supplies both C model.in MBH and downstream
PhysicalScale, so there is no independent generation-only mass constant.

Primary basis: EHT Collaboration, *The persistent shadow of the supermassive
black hole of M87. II. Model comparisons and theoretical interpretations*,
A&A693,A265(2025), Sect.3 p.6, DOI10.1051/0004-6361/202451296. Its main model
library uses6.5e9/16.8; this is an adopted simulation scale, not a new claim that
6.5 is measured more accurately than6.2. The same paper's combined mass inference
is6.16(+.71,-.57)e9 at fixed distance, and includes both values.

Unchanged: existing research.yaml already uses6.5; image FOV remains explicitly
46.2261Rg for new generation; MUNIT3e25, R1/1, cut1, reference norm.1, phase16,
Nprincipal2. No flux calibration, beam/crop/WCS change or old Stokes rescaling.
At fixed D, changing6.2→6.5 changes Rg/tg/angularRg by4.8387%. Fixed FOV/Rg
therefore does not preserve angular camera width. This is deliberate and recorded,
not hidden relabeling. MUNIT remains a benchmark awaiting flux validation.

Old engineering-validation config snapshots/builds, inputs, outputs and .env
are untouched. Example new run ID isbenchmark-211-m65; trying to use an old
run ID with different config is rejected rather than rewriting its config.

Updated tests assert new profile/model.in values and check that generated data
reject a changed mass6.2. Low-level serialization tests with explicit6.2
synthetic inputs remain valid arbitrary-parameter tests, not production defaults.

Validation: complete synthetic suite and Ruff logs archived in
`production-m65-validation/`; real141-footer/two-table --plan passed and still
expands18,330 tasks. No new build or scientific snapshot was run for this
change. Existing compiled binary accepts mass at runtime; actual future
production will record the new19-fieldmodel.in for every task.
