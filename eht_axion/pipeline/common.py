"""Content-addressed stage manifests and transactional output directories."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

SCHEMA_VERSION = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_inputs(paths):
    """Hash bytes, never trust mtime as scientific input identity."""
    result = []
    for value in paths:
        path = Path(value).resolve(strict=True)
        if not path.is_file():
            raise ValueError(f"Not a regular input file: {path}")
        result.append({"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)})
    return result


def revision(stage=None):
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    # Include implementation bytes even in an uncommitted development tree.
    h = hashlib.sha256()
    dependencies = {
        "raptor": ["pipeline/raptor.py", "pipeline/constants.py", "polarization.py", "angles.py"],
        "eht": ["pipeline/eht.py", "polarization.py", "angles.py"],
        "observations": ["pipeline/products.py", "observations.py", "angles.py"],
        "inference": [
            "pipeline/inference.py",
            "pipeline/products.py",
            "pipeline/constants.py",
            "observations.py",
            "angles.py",
            "limits.py",
            "model.py",
            "statistics.py",
        ],
        "figures": ["pipeline/figures.py"],
    }
    group = str(stage).split(":")[0]
    files = dependencies.get(group)
    paths = [root / "eht_axion" / p for p in files] if files else list(root.rglob("*.py"))
    paths.append(Path(__file__).resolve())
    for path in sorted(set(paths)):
        if path.is_file():
            h.update(str(path.relative_to(root)).encode())
            h.update(path.read_bytes())
    return {"git_commit": commit, "python_source_sha256": h.hexdigest()}


def signature(stage, config, inputs):
    payload = {
        "stage": stage,
        "schema": SCHEMA_VERSION,
        "config": config,
        "inputs": inputs,
        "implementation": revision(stage),
    }
    identity = dict(payload, implementation=payload["implementation"]["python_source_sha256"])
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return payload, hashlib.sha256(encoded.encode()).hexdigest()


def complete_manifest(directory):
    path = Path(directory) / "stage.json"
    if not path.is_file():
        raise ValueError(f"Stage incomplete (no stage.json): {directory}")
    manifest = json.loads(path.read_text())
    if manifest.get("schema") != SCHEMA_VERSION or not manifest.get("complete"):
        raise ValueError(f"Invalid stage manifest: {path}")
    return manifest


def verify_stage(directory):
    directory = Path(directory).resolve(strict=True)
    manifest = complete_manifest(directory)
    if not manifest.get("outputs"):
        raise ValueError(f"Empty stage: {directory}")
    for entry in manifest["outputs"]:
        path = (directory / entry["path"]).resolve()
        if not path.is_relative_to(directory) or not path.is_file():
            raise ValueError(f"Unsafe/missing stage output: {path}")
        if path.stat().st_size != entry["size"] or sha256(path) != entry["sha256"]:
            raise ValueError(f"Stage output checksum mismatch: {path}")
    return manifest


@contextlib.contextmanager
def stage_directory(destination, stage, config, input_paths=(), *, rebuild=False):
    """Yield temporary destination, or None for a verified cache hit.

    A manifest is committed last, and a stage is never consumed by mere directory
    existence. Replacement preserves old data until the new directory is complete.
    No hardlinks to inputs. Failed new stages are cleaned, old stages are untouched.
    """
    destination = Path(destination).resolve()
    # Replacing a directory must never replace its raw inputs, including via symlinks.
    for value in input_paths:
        path = Path(value).resolve(strict=True)
        if path == destination or path.is_relative_to(destination):
            raise ValueError(f"Output stage contains an input source: {path}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    inputs = describe_inputs(input_paths)
    payload, key = signature(stage, config, inputs)
    lock = destination.with_name(destination.name + ".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"Stage locked: {lock}; inspect owner before clearing stale lock") from exc
    with os.fdopen(fd, "w") as stream:
        stream.write(json.dumps({"pid": os.getpid(), "started": time.time()}))
    temporary = None
    try:
        if destination.exists() and not rebuild:
            manifest = verify_stage(destination)
            if manifest["signature"] != key:
                raise ValueError(f"Stale cache at {destination}; use a new output root or --rebuild")
            print(f"CACHE HIT {stage}: {destination}", flush=True)
            yield None
            return
        temporary = Path(tempfile.mkdtemp(prefix=destination.name + ".building-", dir=destination.parent))
        started = time.perf_counter()
        yield temporary
        outputs = [
            {"path": str(p.relative_to(temporary)), "size": p.stat().st_size, "sha256": sha256(p)}
            for p in sorted(temporary.rglob("*"))
            if p.is_file() and p.name != "stage.json"
        ]
        if not outputs:
            raise ValueError(f"Stage {stage} produced no files")
        manifest = dict(
            payload,
            signature=key,
            complete=True,
            outputs=outputs,
            elapsed_seconds=time.perf_counter() - started,
        )
        (temporary / "stage.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        backup = destination.with_name(destination.name + ".previous")
        if backup.exists():
            raise ValueError(f"Backup already exists: {backup}; inspect before rebuilding")
        if destination.exists():
            destination.rename(backup)
        try:
            temporary.rename(destination)
            temporary = None
        except BaseException:
            if backup.exists():
                backup.rename(destination)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        print(f"COMPLETE {stage}: {destination}", flush=True)
    finally:
        if temporary is not None:
            shutil.rmtree(temporary)
        lock.unlink(missing_ok=True)
