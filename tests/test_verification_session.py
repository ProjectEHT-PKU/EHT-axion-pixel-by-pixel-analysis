"""Content hashing is invocation-local and rejects mutation, never mtime-trusting."""
import os
import time
from pathlib import Path

import pytest

from eht_axion.pipeline.common import sha256, verification_session


def test_hash_once_per_session_and_again_next_session(tmp_path, monkeypatch):
    source = tmp_path / "fluid"
    source.write_bytes(b"immutable")
    original, reads = Path.open, []

    def opened(path, *args, **kwargs):
        if path == source and args and args[0] == "rb":
            reads.append(path)
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", opened)
    with verification_session():
        expected = sha256(source)
        for _ in range(16):
            assert sha256(source) == expected
    assert len(reads) == 1
    with verification_session():
        assert sha256(source) == expected
    assert len(reads) == 2


def test_parallel_contexts_share_once_only_hash(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from contextvars import copy_context
    import threading

    source = tmp_path / "parallel-fluid"
    source.write_bytes(b"immutable")
    original, lock, calls = Path.open, threading.Lock(), []

    def opened(path, *args, **kwargs):
        if path == source and args and args[0] == "rb":
            with lock:
                calls.append(path)
            time.sleep(.02)
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", opened)
    with verification_session():
        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = [pool.submit(copy_context().run, sha256, source) for _ in range(16)]
            assert len({f.result() for f in futures}) == 1
    assert len(calls) == 1


def test_same_size_mtime_mutation_caught_on_reuse(tmp_path):
    source = tmp_path / "fluid"
    source.write_bytes(b"aaaa")
    before = source.stat()
    with pytest.raises(ValueError, match="changed"):
        with verification_session():
            sha256(source)
            # Container /tmp may update timestamps only at scheduler ticks.
            time.sleep(.05)
            source.write_bytes(b"bbbb")
            os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
            sha256(source)


def test_mutation_after_last_hash_caught_at_exit(tmp_path):
    directory = tmp_path / "fluid.building-case"
    directory.mkdir()
    source = directory / "fluid"
    source.write_bytes(b"aaaa")
    with pytest.raises(ValueError, match="changed"):
        with verification_session():
            sha256(source)
            time.sleep(.05)
            source.write_bytes(b"bbbb")
