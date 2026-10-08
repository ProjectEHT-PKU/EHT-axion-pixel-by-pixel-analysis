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
    paths = [root / "eht_axion" / p for p in files] if files else list((root / "eht_axion").rglob("*.py"))
    paths.extend([Path(__file__).resolve(), root / "eht_axion/pipeline/layout.py"])
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
    if (manifest.get("imported_from") or manifest.get("imported_ancestors")) and not allow_imported:
        raise ValueError("Consuming imported products requires explicit --reuse-imported")
    if recursive:
        seen = set() if _seen is None else _seen
        identity = manifest_path(directory).resolve()
        if identity in seen:
            return manifest
        seen.add(identity)
        for entry in manifest["inputs"]:
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


def reject_directory_symlinks(path):
    lexical = Path(path).absolute()
    if any(p.is_symlink() for p in [lexical, *lexical.parents]):
        raise ValueError("Refusing to rebuild through a shared directory symlink")


@contextlib.contextmanager
def stage_directory(
    destination, stage, config, input_paths=(), *, rebuild=False, reuse_imported=False, imported_from=None
):
    """Yield temporary destination, or None for a verified cache hit.

    A manifest is committed last, and a stage is never consumed by mere directory
    existence. Replacement preserves old data until the new directory is complete.
    No hardlinks to inputs. Failed new stages are cleaned, old stages are untouched.
    """
    if rebuild:
        reject_directory_symlinks(destination)
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
            manifest = verify_stage(destination, allow_imported=reuse_imported)
            if manifest.get("imported_from"):
                if not reuse_imported:
                    raise ValueError("Imported products require explicit --reuse-imported or a new run")
                if manifest["config"] != json.loads(json.dumps(config)) or manifest["inputs"] != inputs:
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
        for entry in inputs:
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
