"""Stage manifests and transactional output directories.

Code is identified by the git commit of this repository: a stage can be written only
when the tracked files are unchanged, so the commit pins every script, module,
experiment.yaml and pixi.lock. A manifest records the stage name, its configuration,
the paths it read, the commit and the files it wrote. External input files carry their
SHA-256 where they are declared (experiment.yaml, raptor/plan.json), not in manifests.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
from importlib.metadata import PackageNotFoundError, version
import os
import platform
import shutil
import subprocess
import tempfile
import threading
import time
from pathlib import Path

_HASHES = {}
_HASH_LOCK = threading.Lock()


def sha256(path: Path) -> str:
    """SHA-256 of an input file; repeated calls on an unchanged file read it once."""
    path = Path(path).resolve(strict=True)
    info = path.stat()
    key = (path, info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)
    with _HASH_LOCK:
        if key in _HASHES:
            return _HASHES[key]
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    with _HASH_LOCK:
        _HASHES[key] = digest.hexdigest()
    return _HASHES[key]


def describe_inputs(paths):
    """Resolved paths of the input files (regular files)."""
    result = []
    for value in paths:
        path = Path(value).resolve(strict=True)
        if not path.is_file():
            raise ValueError(f"Not a regular input file: {path}")
        result.append(str(path))
    return result


REPOSITORY = Path(__file__).resolve().parents[1]


def code_identity(root=REPOSITORY):
    """{"git_commit": full commit}; refuse when tracked files differ from that commit.

    Untracked files (for example git-ignored .env or outputs) do not matter; any
    uncommitted change to a tracked file does, because the commit must identify the code.
    """
    def git(*args):
        return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True,
                              text=True).stdout
    try:
        commit = git("rev-parse", "HEAD").strip()
        changes = git("status", "--porcelain", "--untracked-files=no")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"The code must be a git checkout to write results: {root}") from exc
    if changes.strip():
        raise ValueError("Uncommitted changes to tracked files; commit them before writing results "
                         "(the git commit is the only code identity):\n" + changes.rstrip())
    return {"git_commit": commit}


def product_schema(path):
    """Shapes/dtypes of NPY/NPZ products; large NPY is inspected through mmap."""
    import numpy as np

    path = Path(path)
    if path.suffix == ".npy":
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        schema = {"shape": list(a.shape), "dtype": str(a.dtype)}
        if path.name == "stokes.npy":
            schema.update(
                channels=["I", "Q", "U", "V"],
                units="intensity units of the producing stage; not independently verified",
            )
        return schema
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as z:
            fields = {key: {"shape": list(z[key].shape), "dtype": str(z[key].dtype)} for key in z.files}
        return {"arrays": fields, "serialization": "named arrays; no pickle"}
    return {}


def runtime():
    """Plain version strings of the interpreter and key packages (pixi.lock pins them)."""
    packages = {}
    for name in ("numpy", "scipy", "h5py", "matplotlib", "PyYAML"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = "not-installed"
    return {"python": platform.python_version(), "implementation": platform.python_implementation(),
            "platform": platform.system(), "machine": platform.machine(), "packages": packages}


def manifest_path(directory):
    path = Path(directory) / "manifest.json"
    if not path.is_file():
        raise ValueError(f"Stage incomplete (no manifest.json): {directory}")
    return path


def is_stage_manifest(record):
    """A stage manifest names its stage and is marked complete (written last by stage_directory)."""
    return (isinstance(record, dict) and isinstance(record.get("stage"), str) and bool(record["stage"])
            and record.get("complete") is True)


def complete_manifest(directory):
    path = manifest_path(directory)
    manifest = json.loads(path.read_text())
    if not is_stage_manifest(manifest):
        raise ValueError(f"Invalid stage manifest: {path}")
    return manifest


def verify_stage(directory):
    """A completed stage: complete manifest and every declared output present inside it."""
    directory = Path(directory).resolve(strict=True)
    manifest = complete_manifest(directory)
    if not manifest.get("outputs"):
        raise ValueError(f"Empty stage: {directory}")
    for entry in manifest["outputs"]:
        relative = Path(entry["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe relative product path")
        path = (directory / relative).resolve()
        if not path.is_relative_to(directory) or not path.is_file():
            raise ValueError(f"Missing stage output: {path}")
    return manifest


def require_declared_outputs(directory, paths):
    """A stage covers its declared outputs, never arbitrary extra files."""
    directory = Path(directory).absolute()
    record = complete_manifest(directory)
    declared = {entry["path"] for entry in record["outputs"]}
    for value in paths:
        path = Path(value).absolute()
        try:
            relative = str(path.relative_to(directory))
        except ValueError as exc:
            raise ValueError(f"Consumed product outside its stage: {path}") from exc
        if relative not in declared:
            raise ValueError(f"Consumed product is not a declared stage output: {path}")


def source_closure(paths):
    """Every file the given inputs were made from (through upstream manifests)."""
    pending, seen = list(paths), set()
    while pending:
        path = Path(pending.pop()).resolve(strict=True)
        if path in seen:
            continue
        seen.add(path)
        yield path
        if path.name == "manifest.json":
            record = json.loads(path.read_text())
            if is_stage_manifest(record):
                pending.extend(record["inputs"])
                pending.extend(path.parent / entry["path"] for entry in record["outputs"])


def contained_path(root, relative):
    """root/relative for a relative path that must stay inside root (resolved, must exist)."""
    root = Path(root).resolve(strict=True)
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    target = (root / relative).resolve(strict=True)
    if not target.is_relative_to(root):
        raise ValueError(f"Unsafe relative path: {relative}")
    return target


def reject_directory_symlinks(path):
    lexical = Path(path).absolute()
    if any(p.is_symlink() for p in [lexical, *lexical.parents]):
        raise ValueError("Refusing to write/rebuild through a shared directory symlink")


def require_new_output_file(path):
    """Fail before any write through an existing or dangling output link."""
    path = Path(path)
    reject_directory_symlinks(path.parent)
    if path.is_symlink():
        raise ValueError(f"Refusing output file symlink: {path}")
    if path.exists():
        raise FileExistsError(f"refusing to overwrite stage output: {path}")
    return path


def verify_input_stages(paths):
    """Inputs named manifest.json must be completed stages."""
    for value in paths:
        path = Path(value)
        if path.name == "manifest.json":
            record = json.loads(path.read_text())
            if is_stage_manifest(record):
                verify_stage(path.parent)


@contextlib.contextmanager
def stage_directory(destination, stage, config, input_paths=()):
    """Yield a temporary destination, or None when an equal stage already exists.

    The manifest is written last and the directory is published by one rename; a
    failed new stage is removed and an existing stage is never replaced. An existing
    stage is reused only when it records the same stage name, config and input paths;
    this lets the resumable RAPTOR steps skip finished stages. The result scripts
    refuse existing output roots themselves.
    """
    reject_directory_symlinks(destination)
    destination = Path(destination).resolve()
    input_paths = list(input_paths)
    verify_input_stages(input_paths)
    # Replacing a directory must never replace its inputs, including via symlinks.
    for path in source_closure(input_paths):
        if path == destination or path.is_relative_to(destination):
            raise ValueError(f"Output stage contains an input source: {path}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    inputs = describe_inputs(input_paths)
    lock = destination.with_name(destination.name + ".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"Stage locked: {lock}; inspect owner before clearing stale lock") from exc
    with os.fdopen(fd, "w") as stream:
        stream.write(json.dumps({"pid": os.getpid(), "started": time.time()}))
    temporary = None
    try:
        if destination.exists():
            manifest = verify_stage(destination)
            recorded = {k: manifest.get(k) for k in ("stage", "config", "inputs")}
            if recorded != json.loads(json.dumps({"stage": stage, "config": config, "inputs": inputs})):
                raise ValueError(f"A different stage exists at {destination}; use a new output root")
            print(f"EXISTING {stage}: {destination}", flush=True)
            yield None
            return
        code = code_identity()
        temporary = Path(tempfile.mkdtemp(prefix=destination.name + ".building-", dir=destination.parent))
        started = time.perf_counter()
        yield temporary
        files = [p for p in sorted(temporary.rglob("*")) if p.is_file() and p.name != "manifest.json"]
        if any(p.is_symlink() for p in files):
            raise ValueError("Stage products must be regular files, not symlinks")
        if not files:
            raise ValueError(f"Stage {stage} produced no files")
        outputs = [{"path": str(p.relative_to(temporary)), **product_schema(p)} for p in files]
        manifest = dict(
            stage=stage,
            config=config,
            inputs=inputs,
            implementation=code,
            runtime=runtime(),
            complete=True,
            outputs=outputs,
            elapsed_seconds=time.perf_counter() - started,
        )
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(f"Stage appeared while building: {destination}")
        temporary.rename(destination)
        temporary = None
        print(f"COMPLETE {stage}: {destination}", flush=True)
    finally:
        if temporary is not None:
            shutil.rmtree(temporary)
        lock.unlink(missing_ok=True)
