#!/usr/bin/env python3
"""Verify every file in a shared-data manifest after Git LFS checkout."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASET = ROOT / "shared_data" / "final-inference" / "draft"


def sha256(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", nargs="?", type=Path, default=DEFAULT_DATASET)
    return parser.parse_args()


def main() -> None:
    dataset = parse_args().dataset.resolve()
    manifest_path = dataset / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    failures: list[str] = []
    for entry in manifest["files"]:
        path = dataset / entry["path"]
        if not path.is_file():
            failures.append(f"missing: {entry['path']}")
            continue
        if path.stat().st_size != entry["bytes"]:
            failures.append(
                f"size: {entry['path']} expected {entry['bytes']}, got {path.stat().st_size}"
            )
            continue
        actual_hash = sha256(path)
        if actual_hash != entry["sha256"]:
            failures.append(
                f"sha256: {entry['path']} expected {entry['sha256']}, got {actual_hash}"
            )
            continue
        # Loading catches a common failure mode where LFS pointers were checked
        # out without their binary objects.
        try:
            if entry.get("format") == "numpy-pickled-dictionary":
                np.load(path, allow_pickle=True).item()
            else:
                np.load(path, mmap_mode="r", allow_pickle=False)
        except Exception as error:  # pragma: no cover - diagnostic path
            failures.append(f"load: {entry['path']}: {error}")
            continue
        print(f"OK {entry['path']}")
    if failures:
        raise SystemExit("Shared-data verification failed:\n" + "\n".join(failures))
    total = sum(entry["bytes"] for entry in manifest["files"])
    print(
        f"Verified {len(manifest['files'])} files, {total / 2**20:.1f} MiB, "
        f"revision={manifest['revision']}, status={manifest['status']}"
    )


if __name__ == "__main__":
    main()
