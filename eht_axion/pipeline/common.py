"""Content-addressed stage manifests and transactional output directories."""

from __future__ import annotations

import contextlib
from contextvars import ContextVar
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

SCHEMA_VERSION = 1
_VERIFICATION = ContextVar("content_verification_session", default=None)
_STAGE_TEMPORARIES = ContextVar("active_stage_temporaries", default=())


class _VerificationRecords(dict):
    """Shared across copied worker contexts; protect once-only content reads."""

    def __init__(self):
        super().__init__()
        self.lock = threading.RLock()


def _file_identity(path):
    info = path.stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


@contextlib.contextmanager
def verification_session():
    """Hash shared immutable inputs once per invocation, never across invocations.

    Every reuse checks file identity, including ctime, and successful exit checks
    all hashed files again. This is not a persistent mtime-based content cache.
    Scope this around a batch of jobs; temporary renamed outputs are not memoized.
    """
    current = _VERIFICATION.get()
    if current is not None:
        yield
        return
    records = _VerificationRecords()
    token = _VERIFICATION.set(records)
    try:
        yield
        for path, (identity, _) in records.items():
            if not path.is_file() or _file_identity(path) != identity:
                raise ValueError(f"Input changed during verification session: {path}")
    finally:
        _VERIFICATION.reset(token)


def sha256(path: Path) -> str:
    records = _VERIFICATION.get()
    guard = records.lock if records is not None else contextlib.nullcontext()
    with guard:
        return _sha256_verified(path, records)


def _sha256_verified(path, records):
    path = Path(path).resolve(strict=True)
    # Atomic stage temporaries move on publish; never bind them to session paths.
    memoize = records is not None and not any(
        path.is_relative_to(root) for root in _STAGE_TEMPORARIES.get()
    )
    before = _file_identity(path)
    if memoize and path in records:
        identity, result = records[path]
        if before != identity:
            raise ValueError(f"Input changed during verification session: {path}")
        return result
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if _file_identity(path) != before:
        raise ValueError(f"Input changed while hashing: {path}")
    result = digest.hexdigest()
    if memoize:
        records[path] = (before, result)
    return result


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
    paths = [root / "eht_axion" / p for p in files] if files else list((root / "eht_axion").rglob("*.py"))
    paths.extend([
        Path(__file__).resolve(), root / "eht_axion/pipeline/layout.py",
        root / "eht_axion/pipeline/workflow.py", root / "eht_axion/pipeline/config.py",
        root / "eht_axion/stokes.py",
    ])
    for path in sorted(set(paths)):
        if path.is_file():
            h.update(str(path.relative_to(root)).encode())
            h.update(path.read_bytes())
    return {"git_commit": commit, "python_source_sha256": h.hexdigest()}


def product_schema(path):
    """Shapes/dtypes are metadata; inspect large NPY through mmap, not a RAM copy."""
    import numpy as np

    path = Path(path)
    if path.suffix == ".npy":
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        schema = {"shape": list(a.shape), "dtype": str(a.dtype)}
        if path.name == "stokes.npy":
            schema.update(
                channels=["I", "Q", "U", "V"],
                units="Original intensity units preserved; not independently verified",
            )
        elif path.name == "evpa-1d.npy":
            schema.update(axes=["ensemble", "angle"], units="rad", angle_period="pi")
        return schema
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as z:
            fields = {key: {"shape": list(z[key].shape), "dtype": str(z[key].dtype)} for key in z.files}
        return {"arrays": fields, "serialization": "named arrays; no pickle"}
    return {}


def runtime_identity():
    """Bind stage reuse to the resolved software environment, not only source."""
    root = Path(__file__).resolve().parents[2]
    packages = {}
    for name in ("numpy", "scipy", "h5py", "matplotlib", "PyYAML"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = "not-installed"
    lock = root / "pixi.lock"
    return {
        "python": platform.python_version(), "implementation": platform.python_implementation(),
        "platform": platform.system(), "machine": platform.machine(), "packages": packages,
        "pixi_lock_sha256": sha256(lock) if lock.is_file() else None,
    }


def signature(stage, config, inputs):
    payload = {
        "stage": stage,
        "schema": SCHEMA_VERSION,
        "config": config,
        "inputs": inputs,
        "implementation": revision(stage),
        "runtime": runtime_identity(),
    }
    identity = dict(payload, implementation=payload["implementation"]["python_source_sha256"])
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return payload, hashlib.sha256(encoded.encode()).hexdigest()


def manifest_path(directory):
    directory = Path(directory)
    current, legacy = directory / "manifest.json", directory / "stage.json"
    if current.is_file() and legacy.is_file():
        raise ValueError(f"Ambiguous stage manifests: {directory}")
    if current.is_file():
        return current
    if legacy.is_file():
        return legacy
    raise ValueError(f"Stage incomplete (no manifest.json): {directory}")


def complete_manifest(directory):
    path = manifest_path(directory)
    manifest = json.loads(path.read_text())
    if manifest.get("schema") != SCHEMA_VERSION or not manifest.get("complete"):
        raise ValueError(f"Invalid stage manifest: {path}")
    return manifest


def verify_stage(directory, *, recursive=False, allow_imported=True, _seen=None):
    directory = Path(directory).resolve(strict=True)
    manifest = complete_manifest(directory)
    origin = manifest.get("imported_from") or {}
    if "manifest_sha256" in origin:
        # Earlier migrations kept only a receipt. Do not silently promote that
        # receipt to the new live-verified provenance contract.
        origin_directory = Path(origin["path"]).resolve(strict=True)
        expected_manifest = manifest_path(origin_directory).resolve()
        records = [entry for entry in manifest.get("origin_inputs", [])
                   if Path(entry["path"]).resolve() == expected_manifest]
        if len(records) != 1 or records[0]["sha256"] != origin["manifest_sha256"]:
            raise ValueError("Unverifiable legacy migration origin; migrate again to a new destination")
    if (manifest.get("imported_from") or manifest.get("imported_ancestors")) and not allow_imported:
        raise ValueError("Consuming imported products requires explicit --reuse-imported")
    if recursive:
        seen = set() if _seen is None else _seen
        identity = manifest_path(directory).resolve()
        if identity in seen:
            return manifest
        seen.add(identity)
        for entry in manifest["inputs"] + manifest.get("origin_inputs", []):
            path = Path(entry["path"])
            if not path.is_file() or path.stat().st_size != entry["size"] or sha256(path) != entry["sha256"]:
                raise ValueError(f"Stale upstream lineage: {path}")
            if path.name in {"manifest.json", "stage.json"}:
                try:
                    upstream = json.loads(path.read_text())
                except (ValueError, OSError):
                    upstream = {}
                if upstream.get("complete") and upstream.get("schema") == SCHEMA_VERSION:
                    verify_stage(path.parent, recursive=True, allow_imported=allow_imported, _seen=seen)
    if not manifest.get("outputs"):
        raise ValueError(f"Empty stage: {directory}")
    for entry in manifest["outputs"]:
        path = (directory / entry["path"]).resolve()
        lexical = directory / entry["path"]
        if Path(entry["path"]).is_absolute() or ".." in Path(entry["path"]).parts:
            raise ValueError("Unsafe relative product path")
        external = entry.get("external_symlink")
        if external is not None:
            if not lexical.is_symlink() or path != Path(external).resolve():
                raise ValueError(f"Changed imported product link: {lexical}")
        elif not path.is_relative_to(directory):
            raise ValueError(f"Unsafe external stage output: {path}")
        if not path.is_file():
            raise ValueError(f"Unsafe/missing stage output: {path}")
        if path.stat().st_size != entry["size"] or sha256(path) != entry["sha256"]:
            raise ValueError(f"Stage output checksum mismatch: {path}")
    return manifest


def require_declared_outputs(directory, paths):
    """A verified stage covers declared outputs, never arbitrary extra files."""
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
    """Enumerate all declared source bytes before replacing any output tree."""
    pending, seen = list(paths), set()
    while pending:
        path = Path(pending.pop()).resolve(strict=True)
        if path in seen:
            continue
        seen.add(path)
        yield path
        if path.name in {"manifest.json", "stage.json"}:
            record = json.loads(path.read_text())
            if record.get("complete") and record.get("schema") == SCHEMA_VERSION:
                pending.extend(entry["path"] for entry in record["inputs"]
                               + record.get("origin_inputs", []))
                pending.extend(path.parent / entry["path"] for entry in record["outputs"])


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


def verify_input_stages(paths, *, allow_imported=False):
    """Verify declared stage manifests, including their full scientific ancestry."""
    seen = set()
    for value in paths:
        path = Path(value)
        if path.name in {"manifest.json", "stage.json"}:
            record = json.loads(path.read_text())
            if record.get("complete") and record.get("schema") == SCHEMA_VERSION:
                verify_stage(path.parent, recursive=True, allow_imported=allow_imported, _seen=seen)


@contextlib.contextmanager
def stage_directory(
    destination, stage, config, input_paths=(), *, rebuild=False, reuse_imported=False, imported_from=None,
    origin_paths=(),
):
    """Yield temporary destination, or None for a verified cache hit.

    A manifest is committed last, and a stage is never consumed by mere directory
    existence. Replacement preserves old data until the new directory is complete.
    No hardlinks to inputs. Failed new stages are cleaned, old stages are untouched.
    """
    reject_directory_symlinks(destination)
    destination = Path(destination).resolve()
    input_paths, origin_paths = list(input_paths), list(origin_paths)
    verify_input_stages(
        input_paths + origin_paths, allow_imported=reuse_imported or imported_from is not None,
    )
    # Replacing a directory must never replace its raw inputs, including via symlinks.
    for path in source_closure(input_paths + origin_paths):
        if path == destination or path.is_relative_to(destination):
            raise ValueError(f"Output stage contains an input source: {path}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    inputs = describe_inputs(input_paths)
    payload, key = signature(stage, config, inputs)
    origins = describe_inputs(origin_paths)
    if origins:
        payload["origin_inputs"] = origins
    lock = destination.with_name(destination.name + ".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"Stage locked: {lock}; inspect owner before clearing stale lock") from exc
    with os.fdopen(fd, "w") as stream:
        stream.write(json.dumps({"pid": os.getpid(), "started": time.time()}))
    temporary = None
    temporary_token = None
    try:
        if destination.exists() and not rebuild:
            manifest = verify_stage(destination, recursive=True, allow_imported=reuse_imported)
            if manifest.get("imported_from"):
                if not reuse_imported:
                    raise ValueError("Imported products require explicit --reuse-imported or a new run")
                if (manifest["config"] != json.loads(json.dumps(config)) or manifest["inputs"] != inputs
                        or (origins and manifest.get("origin_inputs", []) != origins)):
                    raise ValueError("Imported product settings/lineage differ; cannot reuse")
                print(f"IMPORTED REUSE {stage}: {destination} (not recomputed)", flush=True)
                yield None
                return
            if manifest["signature"] != key:
                raise ValueError(f"Stale cache at {destination}; use a new output root or --rebuild")
            print(f"CACHE HIT {stage}: {destination}", flush=True)
            yield None
            return
        temporary = Path(tempfile.mkdtemp(prefix=destination.name + ".building-", dir=destination.parent))
        temporary_token = _STAGE_TEMPORARIES.set(_STAGE_TEMPORARIES.get() + (temporary,))
        started = time.perf_counter()
        yield temporary
        outputs = [
            {
                "path": str(p.relative_to(temporary)),
                "size": p.stat().st_size,
                "sha256": sha256(p),
                **product_schema(p),
            }
            for p in sorted(temporary.rglob("*"))
            if p.is_file() and p.name not in {"stage.json", "manifest.json"}
        ]
        for entry in outputs:
            product = temporary / entry["path"]
            if product.is_symlink():
                if imported_from is None:
                    raise ValueError("External product symlinks only allowed in explicit migration")
                entry["external_symlink"] = str(product.resolve(strict=True))
        if not outputs:
            raise ValueError(f"Stage {stage} produced no files")
        manifest = dict(
            payload,
            signature=key,
            complete=True,
            outputs=outputs,
            elapsed_seconds=time.perf_counter() - started,
        )
        ancestors = []
        for entry in inputs + origins:
            path = Path(entry["path"])
            if path.name in {"manifest.json", "stage.json"}:
                upstream = json.loads(path.read_text())
                if upstream.get("imported_from") or upstream.get("imported_ancestors"):
                    ancestors.append(
                        {
                            "manifest": str(path),
                            "sha256": entry["sha256"],
                            "original": upstream.get("imported_from"),
                            "ancestors": upstream.get("imported_ancestors", []),
                        }
                    )
        if ancestors:
            manifest["imported_ancestors"] = ancestors
        if imported_from is not None:
            manifest["imported_from"] = imported_from
            manifest["scientific_recomputation"] = False
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        backup = destination.with_name(destination.name + ".previous")
        if backup.exists() or backup.is_symlink():
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
        if temporary_token is not None:
            _STAGE_TEMPORARIES.reset(temporary_token)
        if temporary is not None:
            shutil.rmtree(temporary)
        lock.unlink(missing_ok=True)
