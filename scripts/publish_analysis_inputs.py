#!/usr/bin/env python3
"""Publish RAPTOR Stokes cubes and EHT image ensembles into Git LFS.

The 12.27 GB res=400 cube is sharded along its time axis so every output is a
standalone NumPy file below common GitHub LFS per-object limits. Unsharded files
are hard-linked when possible to avoid duplicate local working-tree storage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import shutil
import subprocess
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap

ROOT = Path(__file__).resolve().parents[1]
DATASET_ROOT = ROOT / "shared_data" / "analysis-inputs"
ALPHAS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44)
DAYS = (95, 96, 100, 101)
TIME_START = 1600
TIME_STEP = 10


def sha256(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def hardlink_or_copy(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        return "copy"


def array_entry(path: Path, root: Path, **metadata) -> dict[str, object]:
    array = np.load(path, mmap_mode="r", allow_pickle=False)
    return {
        "path": str(path.relative_to(root)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        **metadata,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("revision", nargs="?", default="draft")
    parser.add_argument("--force-draft", action="store_true")
    parser.add_argument("--time-shard-size", type=int, default=10)
    return parser.parse_args()


def validate_revision(args: argparse.Namespace) -> Path:
    if args.revision != "draft" and not args.revision.startswith("v"):
        raise SystemExit("revision must be 'draft' or start with 'v'")
    destination = DATASET_ROOT / args.revision
    if destination.exists() and args.revision != "draft":
        raise SystemExit(f"refusing to overwrite frozen dataset {destination}")
    if destination.exists() and not args.force_draft:
        raise SystemExit("draft exists; pass --force-draft to update it")
    if args.time_shard_size <= 0:
        raise SystemExit("time shard size must be positive")
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def publish_raptor(destination: Path, shard_size: int) -> tuple[list[dict], dict]:
    source_dir = ROOT / "GRMHD" / "data_in"
    output_dir = destination / "raptor-stokes"
    files: list[dict] = []
    sources = []

    for alpha in ALPHAS:
        source = source_dir / f"GRMHD_Data_1600-3000_alpha={alpha}_res=100.npy"
        target = output_dir / "res-100" / f"alpha-{alpha}.npy"
        method = hardlink_or_copy(source, target)
        files.append(
            array_entry(
                target,
                destination,
                alpha=alpha,
                resolution=100,
                time_start=TIME_START,
                time_step=TIME_STEP,
                state_axis="0=no axion; 1..16=axion phases 0..15",
                stokes_axis=["I", "Q", "U", "V"],
            )
        )
        sources.append(
            {
                "source": str(source.relative_to(ROOT)),
                "published": str(target.relative_to(destination)),
                "local_materialization": method,
            }
        )
        print(f"published res=100 alpha={alpha} via {method}", flush=True)

    source = source_dir / "GRMHD_Data_1600-3000_alpha=0.3_res=400.npy"
    cube = np.load(source, mmap_mode="r")
    shard_dir = output_dir / "res-400" / "alpha-0.3"
    shard_dir.mkdir(parents=True, exist_ok=True)
    expected_names = set()
    for start in range(0, cube.shape[0], shard_size):
        stop = min(start + shard_size, cube.shape[0])
        first_time = TIME_START + start * TIME_STEP
        last_time = TIME_START + (stop - 1) * TIME_STEP
        name = f"stokes-t{first_time}-t{last_time}.npy"
        expected_names.add(name)
        target = shard_dir / name
        output = open_memmap(
            target, mode="w+", dtype=cube.dtype, shape=(stop - start,) + cube.shape[1:]
        )
        output[:] = cube[start:stop]
        output.flush()
        del output
        if not np.array_equal(np.load(target, mmap_mode="r"), cube[start:stop]):
            raise RuntimeError(f"shard verification failed: {target}")
        files.append(
            array_entry(
                target,
                destination,
                alpha=0.3,
                resolution=400,
                source_indices=[start, stop],
                time_values=list(range(first_time, last_time + TIME_STEP, TIME_STEP)),
                state_axis="0=no axion; 1..16=axion phases 0..15",
                stokes_axis=["I", "Q", "U", "V"],
            )
        )
        print(f"published and verified {name}", flush=True)
    for obsolete in shard_dir.glob("*.npy"):
        if obsolete.name not in expected_names:
            obsolete.unlink()
    source_summary = {
        "source": str(source.relative_to(ROOT)),
        "bytes": source.stat().st_size,
        "sha256": sha256(source),
        "shape": list(cube.shape),
        "dtype": str(cube.dtype),
        "shard_size": shard_size,
        "shards": len(expected_names),
    }
    del cube
    return files, {"unsharded": sources, "res_400_source": source_summary}


def publish_transport_diagnostics(destination: Path) -> tuple[list[dict], list[dict]]:
    """Publish the single-snapshot transfer-term decomposition used by Fig. 4.

    Exact upstream RAPTOR switch files are not present in this repository; the
    semantic labels below are inferred from the analysis notebook and retained
    as explicit provenance caveats in the manifest.
    """
    source_dir = ROOT / "GRMHD" / "data_in"
    output_dir = destination / "raptor-transport-diagnostics" / "snapshot-2110-res-400"
    definitions = {
        "common": "normal/full transfer configuration",
        "common_noaxion": "no-axion background for the common configuration",
        "jonly": "emissivity-only run; absorption/Faraday propagation disabled",
        "eq_jonly": "equatorial-plane emissivity-only run",
        "fake": "artificial emissivity test run for analytical comparison",
    }
    files = []
    sources = []
    for variant, interpretation in definitions.items():
        source = source_dir / f"grmhdData2110_400_{variant}.npy"
        target = output_dir / f"{variant}.npy"
        method = hardlink_or_copy(source, target)
        files.append(
            array_entry(
                target,
                destination,
                snapshot=2110,
                resolution=400,
                transport_variant=variant,
                interpretation=interpretation,
                interpretation_source="inferred from GRMHD/main.ipynb; upstream RAPTOR configs unavailable",
                state_axis=(
                    "single no-axion image"
                    if variant == "common_noaxion"
                    else "0=no axion; 1..16=axion phases 0..15"
                ),
                stokes_axis=["I", "Q", "U", "V"],
            )
        )
        sources.append(
            {
                "source": str(source.relative_to(ROOT)),
                "published": str(target.relative_to(destination)),
                "local_materialization": method,
            }
        )
        print(f"published transport diagnostic {variant} via {method}", flush=True)

    auxiliary = {
        "r_coord_400.npy": "equatorial-plane crossing radius per screen pixel",
        "axion_field_coord_400_16.npy": "analytical axion field for 16 phases",
    }
    for name, interpretation in auxiliary.items():
        source = source_dir / name
        target = output_dir / name
        method = hardlink_or_copy(source, target)
        files.append(
            array_entry(
                target,
                destination,
                snapshot=2110,
                resolution=400,
                interpretation=interpretation,
            )
        )
        sources.append(
            {
                "source": str(source.relative_to(ROOT)),
                "published": str(target.relative_to(destination)),
                "local_materialization": method,
            }
        )
        print(f"published transport auxiliary {name} via {method}", flush=True)
    return files, sources


def publish_eht(destination: Path) -> tuple[list[dict], list[dict]]:
    source_dir = ROOT / "EHT-Imaging" / "data_in"
    output_dir = destination / "eht-image-ensemble"
    files = []
    sources = []
    for day in DAYS:
        source = source_dir / f"230_{day}.pkl"
        target = output_dir / source.name
        method = hardlink_or_copy(source, target)
        with target.open("rb") as stream:
            data = pickle.load(stream)
        arrays = {
            key: {"shape": list(value.shape), "dtype": str(value.dtype)}
            for key, value in data.items()
            if key != "meta_data"
        }
        files.append(
            {
                "path": str(target.relative_to(destination)),
                "bytes": target.stat().st_size,
                "sha256": sha256(target),
                "format": "pickle",
                "day": day,
                "image_count": len(data["meta_data"]),
                "arrays": arrays,
                "metadata_fields": sorted(data["meta_data"][0]),
            }
        )
        sources.append(
            {
                "source": str(source.relative_to(ROOT)),
                "published": str(target.relative_to(destination)),
                "local_materialization": method,
            }
        )
        print(f"published EHT day={day} via {method}", flush=True)
    return files, sources


def main() -> None:
    args = parse_args()
    destination = validate_revision(args)
    raptor_files, raptor_sources = publish_raptor(destination, args.time_shard_size)
    transport_files, transport_sources = publish_transport_diagnostics(destination)
    eht_files, eht_sources = publish_eht(destination)
    files = sorted(
        raptor_files + transport_files + eht_files, key=lambda item: item["path"]
    )
    commit = subprocess.check_output(
        ["git", "rev-parse", "origin/audit/reproduce-fix-final-results"],
        cwd=ROOT,
        text=True,
    ).strip()
    manifest = {
        "dataset": "EHT axion analysis inputs after RAPTOR and EHT-Imaging",
        "schema_version": 1,
        "revision": args.revision,
        "status": "development" if args.revision == "draft" else "frozen",
        "analysis_commit": commit,
        "scope": {
            "raptor": "Stokes IQUV cubes, 141 times, no-axion + 16 phases",
            "eht": "four days, 4800 image hyperparameter combinations per day",
            "transport_diagnostics": "snapshot 2110 transfer-term decomposition used by Fig. 4",
        },
        "files": files,
        "file_count": len(files),
        "total_bytes": sum(item["bytes"] for item in files),
        "provenance": {
            "raptor": raptor_sources,
            "transport_diagnostics": transport_sources,
            "eht": eht_sources,
        },
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Published {len(files)} files ({manifest['total_bytes'] / 1e9:.3f} GB) to {destination}",
        flush=True,
    )


if __name__ == "__main__":
    main()
