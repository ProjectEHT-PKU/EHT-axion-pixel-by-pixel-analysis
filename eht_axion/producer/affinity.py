"""Linux CPU allocation: explicit disjoint sets, physical-core-first policy."""
from __future__ import annotations

import os
from pathlib import Path
import re
import shutil

TOPOLOGY = Path('/sys/devices/system/cpu')


def parse_cpu_list(value):
    if not isinstance(value, str) or not value or not re.fullmatch(r"\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*", value):
        raise ValueError("CPU list must be integers/ranges such as 0-5,8,10-12")
    result = []
    for token in value.split(','):
        ends = [int(x) for x in token.split('-')]
        start, stop = ends[0], ends[-1]
        if stop < start or stop-start > 100000:
            raise ValueError("Invalid/reversed/oversized CPU range")
        result += list(range(start, stop+1))
    if len(set(result)) != len(result):
        raise ValueError("Duplicate CPU IDs in CPU list")
    return sorted(result)


def available_cpus():
    if not hasattr(os, 'sched_getaffinity'):
        raise ValueError("CPU affinity requires Linux sched_getaffinity")
    return sorted(os.sched_getaffinity(0))


def cpu_topology(cpus):
    rows = []
    for cpu in sorted(cpus):
        folder = TOPOLOGY / f'cpu{cpu}' / 'topology'
        try:
            socket = int((folder / 'physical_package_id').read_text().strip())
            core = int((folder / 'core_id').read_text().strip())
        except (ValueError, OSError) as exc:
            raise ValueError(f"Cannot determine physical core topology for CPU {cpu}") from exc
        rows.append({'cpu': cpu, 'socket': socket, 'core': core})
    return rows


def allocate_cpu_sets(workers, threads, mode='none', cpu_list=None):
    if type(workers) is not int or workers < 1 or type(threads) is not int or threads < 1:
        raise ValueError("workers and threads must be positive integers")
    if mode not in {'none', 'physical', 'logical'}:
        raise ValueError("cpu_affinity must be none, physical or logical")
    if mode == 'none':
        if cpu_list is not None:
            raise ValueError("--cpu-list requires physical or logical CPU affinity")
        if hasattr(os, 'sched_getaffinity') and workers*threads > len(os.sched_getaffinity(0)):
            raise ValueError("workers*threads exceeds available logical CPUs")
        return [None] * workers
    allowed = available_cpus()
    selected = allowed if cpu_list is None else parse_cpu_list(cpu_list)
    if not set(selected).issubset(allowed):
        raise ValueError("CPU list includes CPUs outside this process's allowed affinity")
    if mode == 'physical':
        cores = {}
        for row in cpu_topology(selected):
            cores.setdefault((row['socket'], row['core']), row['cpu'])
        selected = list(cores.values())
    if workers*threads > len(selected):
        raise ValueError(f"Need {workers*threads} {mode} CPUs, only {len(selected)} available")
    if shutil.which('taskset') is None:
        raise FileNotFoundError("CPU pinning requires taskset (util-linux)")
    return [selected[i*threads:(i+1)*threads] for i in range(workers)]


def check_memory(workers, minimum_per_worker_bytes=int(2.5 * 1024**3)):
    """Conservative Linux headroom check, including effective cgroup-v2 limits."""
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be positive")
    memory_file = Path('/proc/meminfo')
    if not memory_file.is_file():
        return
    memory = dict(line.split(':', 1) for line in memory_file.read_text().splitlines())
    available = int(memory['MemAvailable'].split()[0]) * 1024
    # Container root is normally the effective cgroup. Inspect parents too;
    # deployments with a visible unified hierarchy can constrain any ancestor.
    membership = Path('/proc/self/cgroup')
    group = Path('/sys/fs/cgroup')
    if membership.is_file():
        unified = [line[3:] for line in membership.read_text().splitlines() if line.startswith('0::')]
        if unified and unified[0] != '/':
            candidate = group / unified[0].lstrip('/')
            if candidate.is_dir():
                group = candidate
    for folder in [group, *group.parents]:
        if not folder.is_relative_to('/sys/fs/cgroup'):
            break
        limit_file, current_file = folder / 'memory.max', folder / 'memory.current'
        if limit_file.is_file() and current_file.is_file():
            limit = limit_file.read_text().strip()
            if limit != 'max':
                remaining = max(0, int(limit)-int(current_file.read_text().strip()))
                available = min(available, remaining)
    required = workers * minimum_per_worker_bytes + 2 * 1024**3
    if available < required:
        raise ValueError(f"Insufficient available RAM: workers={workers} requires approximately "
                         f"{required} bytes with headroom, available {available}")


def launch_prefix(cpu_set):
    if cpu_set is None:
        return []
    taskset = shutil.which('taskset')
    if taskset is None:
        raise FileNotFoundError("CPU pinning requires taskset")
    if not cpu_set or any(type(cpu) is not int or cpu < 0 for cpu in cpu_set):
        raise ValueError("Invalid CPU allocation")
    return [taskset, '--cpu-list', ','.join(str(cpu) for cpu in cpu_set)]
