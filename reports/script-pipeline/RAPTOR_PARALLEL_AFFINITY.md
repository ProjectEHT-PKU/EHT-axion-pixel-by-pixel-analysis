# Bounded RAPTOR concurrency, CPU pinning and eight-response real pilot

## Scope and configuration

User authorized adding Python concurrency/affinity and a small real pilot on
the14TB disk. C science/source pin2005600 and production science are unchanged.

New operational options:
- `--workers` / EHT_RAPTOR_WORKERS: simultaneous RAPTOR processes.
- `--threads` / EHT_RAPTOR_THREADS: OpenMP threads per process.
- `--cpu-affinity none|physical|logical` / EHT_RAPTOR_CPU_AFFINITY.
- `--cpu-list` / EHT_RAPTOR_CPU_LIST: optional allowed CPU subset.

Physical mode selects one allowed hardware thread per(socket,core), not both
SMT siblings, partitions disjoint worker CPU sets, and rejects insufficient CPUs.
Linux taskset applies the set before exec; OMP_PROC_BIND=close,OMP_PLACES=cores
bind the OpenMP threads. No shell=True or multithreaded preexec_fn. Unbound mode
clears inherited OMP/GOMP binding policy instead of silently changing cache state.

Bounded ThreadPoolExecutor submits at most workers tasks, with a private cwd and
context copy per task. A shared locked invocation hash cache avoids repeat raw
reads; temporary-root ContextVars remain per worker. Production checks host and
cgroup-v2 memory headroom, also for direct API calls. No ProcessPool copies of
arrays/plan/fluid. PID/group and per-task wall time are recorded.

Cancellation stops new launches, kills all active process groups, joins workers
and preserves attempts. Tests include an actual main-thread SIGINT with eight
fake sleeping process groups. Individual ordinary failures stay independent.

Tasks deterministically use slot i%workers, then each slot advances by workers;
this prevents fast cache hits changing pinned CPU identity during resume. CPU
sets/thread binding are conservatively included in job cache identity; use the
same allocation on resume, otherwise choose a new run. Workers alone is not
science identity. max_jobs remains selected task-prefix count, not worker count.

## Disk placement

Local ignored.env.raptor now sets:

`EHT_PRODUCTION_ROOT=/home/canoziia/GRMHD_DATA/production/benchmark-211-m65`.

This was verified to be/dev/sda, the12.7TiB filesystem. Existing project outputs
and baseline pilot are untouched. The default remains workers1/threads2 for
safety; explicitly request8x6 for the tested policy. The reviewed example has
these operational keys. Normal analysis.env output location is unchanged.

## Real pilot performed

Precisely eight selected response tasks, not the first eight of the full plan:
- snapshot1600(time15010tg), alpha.3, phase0..7;
-100x100,norm.1,MBH6.5e9,D16.8Mpc,thermal,MUNIT3e25,R1/1,cut1;
-8processes×6OpenMP threads, physical CPUs0..47;
-sets0–5,6–11,12–17,18–23,24–29,30–35,36–41,42–47;
-machine AMD EPYC7R13,48physical/96logical,31GiB RAM.

Run directory:

`/home/canoziia/GRMHD_DATA/production/pilot-m65-a030-t1600-8x6-pinned/`.

The measurement harness used the real production plan and runner; it explicitly
filtered/asserted the eight tasks, saved the request, measured/probed descendants
and every OpenMP TID mask. No binary or physics modification, no full-batch
launch, no400image, nootheralpha/snapshot. Build and synthetic tests completed
before launching science; no test suite competed with the pilot.

| Metric | Measurement |
|---|---:|
| Total batch wall, including runner preflight/verification |81.9811s|
| Individual RAPTOR wall range |74.8211–79.9080s|
|8images/hour normalized throughput |351.3003/h|
| Aggregate peak sampled RSS |16,977,988KiB =16.19147GiB|
| Per-process peak RSS |about2.024GiB|
| Average aggregate CPU cores from sampled/procCPU |40.9759|
| Observed simultaneous RAPTOR count |8|
| Observed threads per process |6|
| All TID CPU masks confined to assigned set |yes, all8|
| Completed/failed jobs |8/0|

/proc sampled CPU totals are lower estimates missing the final fractional
sampling interval. Memory is simultaneous sumRSS, not sum of noncoincident
per-process peaks. Sample interval.25s. CPU average is about85.4% of48physical
cores; do not confuse with whole-machine logicalCPU percentage.

All eight HDF/schema/finite/coordinates/spectrum tests passed. Iflux range
roughly.7336409–.7336503Jy; negativeI pixels0 for each, polarization>I diagnostic
pixels1 each (also seen in baseline). No arbitrary physical cut was introduced.
This validates serialization, not all physics or time-averageflux calibration.

## What the timing supports

18,048 res100 response images at351.3/h gives about51.37h. If baseline costs are
similar and400images cost16res100 equivalents, whole-plan rough estimate is
58.20h (about2.43days), excluding failures/other-system contention and changes
acrossalpha/snapshot.400cost remains unmeasured. This is not a wall-time promise.

Do not claim a speedup relative to the previous39s baseline: that was noaxion,
oneprocess/16threads, while this pilot contains8axisresponseprocesses/6threads.
Neither thread-scaling nor multi-process speedup was isolated by these two runs.
No additional scientific tasks were launched after this batch.

## Validation and evidence

Synthetic tests cover bounded overlap, private jobs, per-task failures/retry,
sharedsourcehash once, context isolation, malformed/excess CPU selection, SMT
selection, actual taskset probe, bindingreceipt, stable pinnedresume, timeout,
actualSIGINT/groupreaping and source mutation. **1059 tests passed in140.05s**;
Ruff andgit diff --check passed. Full logs in `raptor-parallel-validation/`.

Pilot JSON/results/proc samples/harness and validation logs are archived there.
The real binary/HDF/jobmanifests remain in the14TB run, not committed as code.
Independent read-only review requested scope filtering, clearing inherited
unbound affinity and cgroupheadroom; these were implemented before the real
pilot. All sources/old outputs were preserved; no Git push in this turn.
