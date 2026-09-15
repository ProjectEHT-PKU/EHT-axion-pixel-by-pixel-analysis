#!/usr/bin/env python3
"""Publish validated final-inference inputs into the Git-LFS namespace.

The default ``draft`` dataset may evolve with the analysis. Frozen ``vN``
datasets are immutable: this command refuses to overwrite an existing version.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATASET_ROOT = ROOT / "shared_data" / "final-inference"
DAYS = (95, 96, 100, 101)


def sha256(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def describe_array(path: Path) -> dict[str, object]:
    try:
        if path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as archive:
                arrays = {
                    name: {"shape": list(archive[name].shape), "dtype": str(archive[name].dtype)}
                    for name in archive.files
                }
            return {"arrays": arrays}
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        return {"shape": list(array.shape), "dtype": str(array.dtype)}
    except ValueError:
        # Parameter dictionaries are intentionally small object arrays. Their
        # payload is described by alpha keys and per-alpha value shapes.
        dictionary = np.load(path, allow_pickle=True).item()
        return {
            "format": "numpy-pickled-dictionary",
            "keys": [str(key) for key in sorted(dictionary)],
            "value_shapes": {str(key): list(dictionary[key].shape) for key in sorted(dictionary)},
            "value_dtypes": {str(key): str(dictionary[key].dtype) for key in sorted(dictionary)},
        }


def copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "version",
        nargs="?",
        default="draft",
        help="'draft' or immutable release name such as 'v1'",
    )
    parser.add_argument(
        "--force-draft",
        action="store_true",
        help="confirm replacing the mutable draft dataset",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.version != "draft" and not args.version.startswith("v"):
        raise SystemExit("version must be 'draft' or start with 'v'")
    destination = DATASET_ROOT / args.version
    if destination.exists() and args.version != "draft":
        raise SystemExit(f"refusing to overwrite frozen dataset {destination}")
    if destination.exists() and not args.force_draft:
        raise SystemExit("draft exists; pass --force-draft to update it intentionally")

    sources: dict[Path, Path] = {}
    corrected_2d = ROOT / "audit" / "rebuilt_2d"
    corrected_1d = ROOT / "audit" / "rebuilt_1d"
    grmhd = ROOT / "GRMHD" / "data_out"

    sources[corrected_2d / "condition.npy"] = Path("2d/condition.npy")
    for day in DAYS:
        sources[corrected_2d / f"{day}.npz"] = Path(f"2d/{day}.npz")
        sources[corrected_2d / f"{day}_all_data.npy"] = Path(f"2d/{day}_all_data.npy")
        sources[corrected_1d / f"{day}_1d.npz"] = Path(f"1d/{day}_1d.npz")
        sources[corrected_1d / f"{day}_EVPA_1d.npy"] = Path(f"1d/{day}_EVPA_1d.npy")
    sources[grmhd / "A_params_dict.npy"] = Path("2d/A_params_dict.npy")
    sources[grmhd / "phi0_params_dict.npy"] = Path("2d/phi0_params_dict.npy")
    sources[corrected_1d / "A_params_1d_dict.npy"] = Path("1d/A_params_1d_dict.npy")
    sources[corrected_1d / "phi0_params_1d_dict.npy"] = Path("1d/phi0_params_1d_dict.npy")

    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise SystemExit("missing validated inputs:\n" + "\n".join(missing))
    destination.mkdir(parents=True, exist_ok=True)
    for source, relative in sources.items():
        copy(source, destination / relative)

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    files = []
    for path in sorted(destination.rglob("*.np[yz]")):
        files.append(
            {
                "path": str(path.relative_to(destination)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                **describe_array(path),
            }
        )
    manifest = {
        "dataset": "EHT M87* final-inference inputs",
        "schema_version": 1,
        "revision": args.version,
        "status": "development" if args.version == "draft" else "frozen",
        "generated_from_commit": commit,
        "topset_size": 1572,
        "selected_2d_pixels": 757,
        "angular_window_degrees": 10,
        "angular_spacing_degrees": 10,
        "angular_offset": 4,
        "files": files,
        "total_bytes": sum(item["bytes"] for item in files),
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Published {len(files)} files ({manifest['total_bytes'] / 2**20:.1f} MiB) to {destination}"
    )


if __name__ == "__main__":
    main()
