"""Measure one explicitly selected pipeline job. No extra simulation or code changes."""
import json
import os
from pathlib import Path
import resource
import sys
import threading
import time

ROOT = Path('/home/canoziia/MyProject/PHYS/eht-script-pipeline')
RUN = ROOT.parent / 'eht-pipeline-runs/production/pilot-m65-t1600-res100-t16-v2'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from eht_axion.producer import jobs
import generate_raptor

metrics = RUN / 'benchmark'
metrics.mkdir(exist_ok=True)
original = jobs._execute
measurements = []
HZ = os.sysconf('SC_CLK_TCK')


def descendants(pid):
    pending, found = [pid], set()
    while pending:
        parent = pending.pop()
        try:
            children = Path(f'/proc/{parent}/task/{parent}/children').read_text().split()
        except (FileNotFoundError, ProcessLookupError):
            continue
        for value in children:
            child = int(value)
            if child not in found:
                found.add(child)
                pending.append(child)
    return found


def process_sample(pid):
    try:
        stat = Path(f'/proc/{pid}/stat').read_text()
        command = stat[stat.index('(')+1:stat.rindex(')')]
        if command != 'RAPTOR':
            return None
        fields = stat[stat.rindex(')')+2:].split()
        status = {}
        for line in Path(f'/proc/{pid}/status').read_text().splitlines():
            key, _, value = line.partition(':')
            status[key] = value.strip()
        return {'pid': pid, 'cpu_seconds': (int(fields[11])+int(fields[12]))/HZ,
                'rss_kib': int(status.get('VmRSS', '0 kB').split()[0]),
                'hwm_kib': int(status.get('VmHWM', '0 kB').split()[0]),
                'threads': int(status.get('Threads', '0')),
                'state': status.get('State')}
    except (OSError, ProcessLookupError, ValueError):
        return None


def measured(command, cwd, env, stdout, stderr, timeout):
    if measurements or Path(command[0]).name != 'RAPTOR':
        raise RuntimeError('Pilot authorizes exactly one RAPTOR process')
    if command[-1] != '0':
        raise RuntimeError('Pilot authorizes baseline slot0 only')
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    stop = threading.Event()
    samples = []
    start = time.perf_counter()

    def monitor():
        previous = {}
        with (metrics / 'samples.jsonl').open('w') as stream:
            while not stop.is_set():
                now = time.perf_counter()
                for pid in descendants(os.getpid()):
                    sample = process_sample(pid)
                    if sample is None:
                        continue
                    sample['elapsed_seconds'] = now-start
                    if pid in previous:
                        last_time, last_cpu = previous[pid]
                        sample['cpu_percent_one_core'] = 100*(sample['cpu_seconds']-last_cpu)/(now-last_time)
                    previous[pid] = (now, sample['cpu_seconds'])
                    samples.append(sample)
                    stream.write(json.dumps(sample)+'\n')
                    stream.flush()
                stop.wait(.5)

    worker = threading.Thread(target=monitor, daemon=True)
    worker.start()
    try:
        result = original(command, cwd, env, stdout, stderr, timeout)
        return result
    finally:
        elapsed = time.perf_counter()-start
        stop.set()
        worker.join(timeout=3)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = (after.ru_utime-before.ru_utime)+(after.ru_stime-before.ru_stime)
        summary = {'command': command, 'cwd_at_execution': str(cwd),
                   'omp_threads': env['OMP_NUM_THREADS'], 'wall_seconds': elapsed,
                   'user_cpu_seconds': after.ru_utime-before.ru_utime,
                   'system_cpu_seconds': after.ru_stime-before.ru_stime,
                   'cpu_seconds': cpu, 'average_cpu_cores': cpu/elapsed,
                   'average_cpu_percent_one_core': 100*cpu/elapsed,
                   'average_utilization_of_16_threads_percent': 100*cpu/(elapsed*16),
                   'wait4_peak_rss_kib': after.ru_maxrss,
                   'sampled_peak_rss_kib': max((s['rss_kib'] for s in samples), default=0),
                   'sampled_peak_hwm_kib': max((s['hwm_kib'] for s in samples), default=0),
                   'sampled_max_threads': max((s['threads'] for s in samples), default=0),
                   'samples': len(samples), 'timing_scope': 'RAPTOR child only; excludes Python preflight/hash/validation',
                   'linux_wait4_rss_unit': 'KiB'}
        measurements.append(summary)
        (metrics / 'raptor-process.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('MEASURED_RAPTOR '+json.dumps(summary), flush=True)


jobs._execute = measured
start = time.perf_counter()
code = generate_raptor.main(['--run-id','pilot-m65-t1600-res100-t16-v2', '--max-jobs','1',
                             '--threads','16','--timeout-seconds','1800'])
if len(measurements) != 1:
    raise RuntimeError('Pilot did not execute exactly one process')
(metrics / 'invocation.json').write_text(json.dumps({
    'exit_code': code, 'wall_seconds_including_python': time.perf_counter()-start,
    'jobs_executed': len(measurements), 'real_scientific_pilot': True
},indent=2)+'\n')
raise SystemExit(code)
