"""Isolated fake toolchain tests: no C compiler or scientific data required."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from eht_axion import provenance
from eht_axion.provenance import verify_stage
from eht_axion.experiment import read_experiment
from eht_axion.raptor.build import build_raptor

BUILD = read_experiment()["raptor"]["build"]


SETUP = r'''#!/bin/bash
set -eu
cp "$RAPTOR"/run/* .
cp "$RAPTOR"/model/bhac/* .
for arg in "$@"; do
    case "$arg" in
        --code=bhac|--metric=mks|--int=rk4|--grid=SMR|--rad=pol|--sf=sfc|--bflip=plus) ;;
        --sigma-cut=*) macro=SIGMA_CUT; value="${arg#*=}"
            if [ "${BAD_MACRO:-0}" != 1 ]; then
                sed -i "s/^#define $macro .*/#define $macro ($value)/" definitions.h
            fi ;;
        --uniform=*) sed -i "s/^#define UNIF .*/#define UNIF (${arg#*=})/" definitions.h ;;
        --thetae-min=*) sed -i "s/^#define THETAE_MIN .*/#define THETAE_MIN (${arg#*=})/" definitions.h ;;
        --thetae-max=*) sed -i "s/^#define THETAE_MAX .*/#define THETAE_MAX (${arg#*=})/" definitions.h ;;
        *) echo "unsupported fixture argument $arg" >&2; exit 3 ;;
    esac
done
make clean
make all
'''

HEADERS = '''#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#define metric (MKSBHAC)
#define int_method (RK4)
#define AMR 0
#define SMR 1
#define POL (1)
#define UNIF (1)
#define SIGMA_CUT (2)
#define THETAE_MIN (1e-3)
#define THETAE_MAX (100)
#endif
'''

FAKE_TOOL = '''import json, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ['TOOL_LOG'], 'a') as f:
    f.write(json.dumps({'tool': name, 'args': args, 'cwd': os.getcwd(),
                        'RAPTOR': os.environ.get('RAPTOR'),
                        'contaminated': [key for key in ('VPATH', 'CPATH', 'MAKEFLAGS', 'MAKEFILES')
                                         if key in os.environ]}) + '\\n')
if '--version' in args:
    print('fake-toolchain 1.0')
elif name == 'h5cc' and args == ['-show']:
    print('fake-cc -I/fake/hdf/include -L/fake/hdf/lib -lhdf5')
elif name == 'h5cc' and '-dM' in args:
    print(pathlib.Path('definitions.h').read_text())
    print(pathlib.Path('model_definitions.h').read_text())
elif name == 'make' and '-nB' in args:
    print('h5cc -c -I' + os.environ['RAPTOR'] + '/src -I' + os.getcwd() + ' -O2 main.c -o build/main.o')
elif name == 'make' and args == ['clean']:
    pathlib.Path('RAPTOR').unlink(missing_ok=True)
elif name == 'make' and args == ['all']:
    if os.environ.get('FAIL_BUILD') == '1':
        sys.exit(9)
    binary = pathlib.Path('RAPTOR')
    binary.write_text('#!/bin/sh\\nexit 0\\n')
    binary.chmod(0o755)
else:
    sys.exit(8)
'''


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    # Isolate cache identity from concurrent edits elsewhere in the checkout.
    # provenance.signature deliberately hashes *all* Python for unknown stages;
    # another agent changing jobs.py between two calls must invalidate a real
    # cache, but must not make this fake-toolchain unit test race with that agent.
    monkeypatch.setattr(provenance, 'revision', lambda stage=None: {
        'git_commit': 'fixture', 'python_source_sha256': 'fixture-implementation',
    })
    monkeypatch.setattr(provenance, 'runtime_identity', lambda: {'toolchain': 'fixture'})
    source = tmp_path / 'core'
    source.mkdir()
    files = {
        'setup.sh': SETUP, 'run/makefile': 'SOURCES=main.c\n', 'run/model.in': 'official runtime input',
        'src/main.c': 'int main(void) {return 0;}\n', 'model/bhac/definitions.h': HEADERS,
        'model/bhac/model_definitions.h': '#define MODEL_DEFINITIONS_H\n#define SFC 1\n#define BPOL (PLUS)\n',
        'model/bhac/model_functions.h': '/* fixture */\n',
        'model/bhac/model_global_vars.h': '/* fixture */\n',
        'model/bhac/model.c': '/* fixture model */\n', 'model/bhac/grid_mks.in': 'grid fixture\n',
        '.gitignore': 'raw.dat\n',
    }
    for name, content in files.items():
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    git(source, 'init', '-q')
    git(source, 'config', 'user.name', 'test')
    git(source, 'config', 'user.email', 'test@example.invalid')
    git(source, 'add', '.')
    git(source, 'commit', '-qm', 'fixture')
    commit = git(source, 'rev-parse', 'HEAD')
    (source / 'raw.dat').write_bytes(b'ignored simulation data must not be copied')
    tools = tmp_path / 'tools'
    tools.mkdir()
    for name in ('h5cc', 'make'):
        path = tools / name
        path.write_text('#!' + sys.executable + '\n' + FAKE_TOOL)
        path.chmod(0o755)
    # rsync is not executed by the fake official setup, but required by the contract.
    (tools / 'rsync').write_text('#!/bin/sh\nexit 0\n')
    (tools / 'rsync').chmod(0o755)
    monkeypatch.setenv('PATH', str(tools) + os.pathsep + os.environ['PATH'])
    monkeypatch.delenv('HDF5_CC', raising=False)
    monkeypatch.setenv('TOOL_LOG', str(tmp_path / 'tool.log'))
    return source, tmp_path / 'build', commit


def test_build_contract_arguments_snapshot_and_stable_receipt(fixture, monkeypatch):
    source, destination, commit = fixture
    for name in ('VPATH', 'CPATH', 'MAKEFLAGS', 'MAKEFILES'):
        monkeypatch.setenv(name, '/external/poison')
    assert build_raptor(source, destination, BUILD, commit) == destination
    record = verify_stage(destination, recursive=True)
    receipt = json.loads((destination / 'receipt.json').read_text())
    assert receipt['source_commit'] == commit
    assert receipt['runtime_cwd'] == str(destination)
    args = receipt['commands'][0][2:]
    assert args == ['--code=bhac', '--metric=mks', '--int=rk4', '--grid=SMR', '--rad=pol',
                    '--sf=sfc', '--bflip=plus', '--sigma-cut=1.0', '--uniform=0',
                    '--thetae-min=0.001', '--thetae-max=100.0']
    assert receipt['actual_macros']['UNIF'] == '0'
    assert receipt['actual_macros']['SIGMA_CUT'] == '1.0'
    assert receipt['tables_are_build_dependencies'] is False
    assert '-lhdf5' in receipt['toolchain']['hdf5_expansion']
    assert not (destination / 'source/.git').exists()
    assert not (destination / 'source/raw.dat').exists()
    assert (destination / 'grid_mks.in').read_text() == 'grid fixture\n'
    assert all(Path(r['path']).is_file() for r in receipt['source_files'])
    assert all('.building-' not in entry['path'] for entry in record['inputs'])
    assert not any('consumed' in key for key in receipt)
    output_paths = {entry['path'] for entry in record['outputs']}
    assert {'RAPTOR', 'definitions.h', 'model_definitions.h', 'grid_mks.in', 'compile.log',
            'macros.txt', 'receipt.json', 'compiler-flags.txt'} <= output_paths
    calls = [json.loads(line) for line in (destination.parent / 'tool.log').read_text().splitlines()]
    assert all(not call['contaminated'] for call in calls)
    macro_call = next(call for call in calls if '-dM' in call['args'])
    assert '-I' + receipt['command_cwd_history'] in macro_call['args']
    assert '-I' + receipt['command_cwd_history'] + '/source/src' in macro_call['args']
    assert '-O2' in macro_call['args']
    assert 'main.c' not in macro_call['args']
    assert '-c' not in macro_call['args']
    assert '-o' not in macro_call['args']
    assert any(call['args'] == ['all'] for call in calls)
    assert 'SIGMA_CUT (2)' in (source / 'model/bhac/definitions.h').read_text()


def test_cache_and_corrupt_snapshot(fixture):
    source, destination, commit = fixture
    build_raptor(source, destination, BUILD, commit)
    before = (destination / 'RAPTOR').stat().st_mtime_ns
    build_raptor(source, destination, BUILD, commit)
    assert (destination / 'RAPTOR').stat().st_mtime_ns == before
    (destination / 'source/src/main.c').write_text('tampered snapshot')
    with pytest.raises(ValueError, match='checksum mismatch'):
        build_raptor(source, destination, BUILD, commit)


@pytest.mark.parametrize('kind', ['tracked', 'untracked', 'commit'])
def test_refuses_invalid_source_before_build(fixture, kind):
    source, destination, commit = fixture
    if kind == 'tracked':
        (source / 'src/main.c').write_text('dirty')
    elif kind == 'untracked':
        (source / 'extra').write_text('untracked')
    else:
        commit = '0' * 40
    with pytest.raises(ValueError, match='dirty|Wrong RAPTOR commit'):
        build_raptor(source, destination, BUILD, commit)
    assert not destination.exists()
    assert not (destination.parent / 'tool.log').exists()


@pytest.mark.parametrize('settings', [
    {'unknown': True}, {}, {'sigma_cut': float('nan')}, {'sigma_cut': 0}, {'thetae_min': 200},
    {'uniform_output': 0}, {'grid': 'smr'}, {'metric': 'unknown'}, {'model': 'harm3d', 'metric': 'cks'},
])
def test_refuses_invalid_settings(fixture, settings):
    source, destination, commit = fixture
    with pytest.raises(ValueError):
        complete = dict(BUILD, **settings) if settings.keys() - {'unknown'} else settings
        build_raptor(source, destination, complete,
                     commit)
    assert not destination.exists()


def test_cache_rejects_actual_implementation_change(fixture, monkeypatch):
    source, destination, commit = fixture
    build_raptor(source, destination, BUILD, commit)
    monkeypatch.setattr(provenance, 'revision', lambda stage=None: {
        'git_commit': 'fixture', 'python_source_sha256': 'changed-implementation',
    })
    with pytest.raises(ValueError, match='Stale cache'):
        build_raptor(source, destination, BUILD, commit)


def test_failure_is_nonzero_and_transactional(fixture, monkeypatch):
    source, destination, commit = fixture
    monkeypatch.setenv('FAIL_BUILD', '1')
    with pytest.raises(subprocess.CalledProcessError) as failed:
        build_raptor(source, destination, BUILD, commit)
    assert 'compile.log tail' in '\n'.join(failed.value.__notes__)
    assert not destination.exists()
    assert not list(destination.parent.glob('build.building-*'))
    monkeypatch.delenv('FAIL_BUILD')
    build_raptor(source, destination, BUILD, commit)
    verify_stage(destination)


def test_mismatch_macro_does_not_publish(fixture, monkeypatch):
    source, destination, commit = fixture
    monkeypatch.setenv('BAD_MACRO', '1')
    with pytest.raises(ValueError, match='macro mismatch: SIGMA_CUT'):
        build_raptor(source, destination, BUILD, commit)
    assert not destination.exists()


def test_rejects_source_symlink_and_unsafe_destination(fixture):
    source, destination, commit = fixture
    with pytest.raises(ValueError, match='outside'):
        build_raptor(source, source / 'build', BUILD, commit)
    (source / 'escape').symlink_to('/etc/passwd')
    git(source, 'add', 'escape')
    git(source, 'commit', '-qm', 'unsafe symlink')
    with pytest.raises(ValueError, match='symlink/submodule'):
        build_raptor(source, destination, BUILD, git(source, 'rev-parse', 'HEAD'))


def test_output_symlink_never_written(fixture):
    source, destination, commit = fixture
    external = destination.parent / 'external'
    destination.symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match='directory symlink'):
        build_raptor(source, destination, BUILD, commit)
    assert not external.exists()
