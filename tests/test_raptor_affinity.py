"""Synthetic topology and real taskset probe; never run RAPTOR."""
import json
import os
import subprocess
import sys

import pytest

from eht_axion.raptor import affinity


@pytest.fixture
def topology(tmp_path, monkeypatch):
    for cpu in range(8):
        directory = tmp_path / f'cpu{cpu}' / 'topology'
        directory.mkdir(parents=True)
        (directory / 'physical_package_id').write_text('0')
        (directory / 'core_id').write_text(str(cpu % 4))
    monkeypatch.setattr(affinity, 'TOPOLOGY', tmp_path)
    monkeypatch.setattr(affinity, 'available_cpus', lambda: list(range(8)))
    monkeypatch.setattr(affinity.shutil, 'which', lambda name: '/usr/bin/taskset')


def test_physical_does_not_assign_smt_siblings_twice(topology):
    assert affinity.allocate_cpu_sets(2, 2, 'physical') == [[0, 1], [2, 3]]
    with pytest.raises(ValueError, match='only 4'):
        affinity.allocate_cpu_sets(3, 2, 'physical')


def test_logical_explicit_sets_and_available_subset(topology):
    assert affinity.allocate_cpu_sets(2, 2, 'logical', '2-3,6-7') == [[2, 3], [6, 7]]
    assert affinity.allocate_cpu_sets(2, 1, 'physical', '4-7') == [[4], [5]]
    with pytest.raises(ValueError, match='outside'):
        affinity.allocate_cpu_sets(1, 1, 'logical', '9')


@pytest.mark.parametrize('value', ['1,1', '3-1', '-1', '1,,2', '', 'abc', '1-99999999'])
def test_bad_cpu_lists(value):
    with pytest.raises(ValueError):
        affinity.parse_cpu_list(value)


def test_disabled_pinning_rejects_silent_cpu_list(topology):
    assert affinity.allocate_cpu_sets(1, 1) == [None]
    with pytest.raises(ValueError, match='requires'):
        affinity.allocate_cpu_sets(1, 1, cpu_list='0')


@pytest.mark.parametrize('workers,threads', [(0, 1), (True, 1), (1, 0), (1, 1.5)])
def test_invalid_counts(workers, threads):
    with pytest.raises(ValueError):
        affinity.allocate_cpu_sets(workers, threads)


def test_taskset_applies_affinity_before_python_exec():
    if not hasattr(os, 'sched_getaffinity') or affinity.shutil.which('taskset') is None:
        pytest.skip('Linux taskset required')
    cpus = affinity.available_cpus()[:2]
    command = affinity.launch_prefix(cpus) + [sys.executable, '-c',
                'import os,json; print(json.dumps(sorted(os.sched_getaffinity(0))))']
    result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=15)
    assert json.loads(result.stdout) == cpus
    assert set(cpus).issubset(os.sched_getaffinity(0))
