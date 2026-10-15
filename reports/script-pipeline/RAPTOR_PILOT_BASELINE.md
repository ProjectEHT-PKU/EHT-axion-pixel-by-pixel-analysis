# First real RAPTOR pilot: one100x100 baseline,16 OpenMP threads

User authorized one image to measure elapsed time, memory and CPU. Exactly one
RAPTOR child was executed; no axion phase, no400 image and no full batch.

## Configuration

- New default: MBH6.5e9 Msun,D16.8Mpc, camera46.2261Rg,100x100.
- Snapshot1600, actual footer time15010tg, fluid spin.9375.
- norm0 baseline, outputslot0; thermal, MUNIT3e25,R1/1,cut1,230GHz,step.005.
- Core20056006820b6536c96de790dd5fdaad5ada4eae; official setup/make via Python/Pixi.
- Machine AMD EPYC7R13,48 physical cores/96 logical CPUs,31GiB RAM.
- OMP_NUM_THREADS=16, one process, no competing test suite during the timing.

Successful run:

`../eht-pipeline-runs/production/pilot-m65-t1600-res100-t16-v2/`.

## Measured results

| Metric | Result |
|---|---:|
| RAPTOR process wall time |38.990958742s|
| Pipeline invocation including preflight/hash/output verification |43.071108584s|
| Child user CPU time |488.876754s|
| Child system CPU time |1.855132s|
| CPU time / wall time |12.5858 cores on average|
| CPU percentage with100%=one core |1258.58% average|
| Utilization relative to16 configured threads |78.6612%|
| /proc sampled peak RSS/HWM |2129320KiB =2.03068GiB|
| Linux wait4/rusage peak RSS |2122564KiB =2.02423GiB|
| Observed max threads |16|
| /proc samples |77, at approximately0.5s intervals|
| Process/pipeline exit |0 /0|

The two RSS measurements differ slightly due to measurement/accounting methods;
report approximately **2.03GiB**, not false exactness. Rusage differences cover
only the RAPTOR child execution interval, not the build. /proc samples inspect
actual RAPTOR descendants. CPU percentage above100 is normal for multithreading.
The process timing includes fluid initialization, ray tracing and C output;
full invocation additionally includes Python validation/hashing. There was no
build in the invocation timed at43.07s. No assertion of cold filesystem cache.

## Output validation

Completed job manifest, success receipt and HDF schema/coordinates/finite checks
passed. HDF area integral I=0.7336417058011014Jy matches spectrum
0.7336417058010996Jy. This is one snapshot, not time-average0.6Jy calibration.
Diagnostic counts: negative I pixels0; polarization fraction above1 pixels1.
No new numerical clipping or arbitrary physical acceptance threshold was applied.
Output is a validated serialization, not a full physics accuracy certificate.

## Directory guard corrected before C launch

The first attempt executed no RAPTOR: the runner incorrectly rejected the
canonical-run-v1 descriptive manifest at the shared output container as though
it were an immutable completed stage. The intended production/<run-id> child
is now allowed only beneath an explicit non-stage container manifest; completed
stage ancestors, symlinks, malformed/unknown stage markers remain rejected.

Regression tests cover allowed container and forbidden completed stage ancestor.
`tests/test_producer_jobs.py`:38 passed in12.80s; Ruff andgit diff --check pass.
The pre-fix build/request directory is preserved as evidence; the successful
pilot used a fresh-v2run, not a cache re-sign or overwrite. Official core and
scientific settings did not change for the guard fix.

## Evidence and limits

`raptor-pilot-baseline/`: process/invocation JSON, /proc samples, measurement
harness, subprocess and validation log, regression/lint logs. The harness wraps
the actual runner's subprocess boundary without modifying the binary or science.

This measures a no-axion100x100 image at one snapshot and one thread count only.
Axion field lookup/screen rotation adds work;400x400 and other snapshots may cost
differently. Multi-process throughput and memory contention are unmeasured.
No additional scientific tasks were launched after this pilot.
