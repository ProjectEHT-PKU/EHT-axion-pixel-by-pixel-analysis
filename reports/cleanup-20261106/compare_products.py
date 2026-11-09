#!/usr/bin/env python3
"""Compare two product roots file by file: arrays bitwise (else max |diff|), JSON/CSV numerics.

Usage: compare_products.py OLD_ROOT NEW_ROOT NAME [--tol 1e-12]
Writes products-NAME.json (reports/cleanup-20261106 and the verify root).
- .npy/.npz: shape, dtype and bitwise equality (NaN == NaN); otherwise max abs diff.
- .json (incl. manifest 'config' and output digests): every numeric leaf compared
  (|diff| <= tol); string leaves that differ are listed separately (paths, hashes).
- .csv: exact text, else numeric cells; .png: byte identity only (informational).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from _side import write

IGNORED_MANIFEST = ("inputs", "implementation", "runtime", "signature", "elapsed_seconds")


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def arrays_equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return dict(equal=False, reason=f"shape/dtype {a.shape}/{a.dtype} vs {b.shape}/{b.dtype}")
    if a.dtype.kind not in "fciub":
        return dict(equal=bool(np.array_equal(a, b)), max_abs_diff=None)
    equal, worst = True, 0.0
    step = max(1, (64 << 20) // max(1, a[0:1].nbytes)) if a.ndim else 1
    for start in range(0, len(a) if a.ndim else 1, step):
        x = np.asarray(a[start:start + step] if a.ndim else a)
        y = np.asarray(b[start:start + step] if b.ndim else b)
        same = np.array_equal(x, y, equal_nan=a.dtype.kind in "fc")
        if not same:
            equal = False
            with np.errstate(invalid="ignore"):
                d = np.abs(x.astype(np.complex128 if a.dtype.kind == "c" else np.float64) - y)
            finite = d[np.isfinite(d)]
            nan_mismatch = bool(np.any(np.isnan(x) != np.isnan(y))) if a.dtype.kind in "fc" else False
            worst = max(worst, float(finite.max()) if finite.size else 0.0, np.inf if nan_mismatch else 0.0)
    return dict(equal=equal, max_abs_diff=worst)


def json_leaves(a, b, path, tol, numeric, strings):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                strings.append(dict(path=f"{path}/{k}", old=a.get(k, "<missing>"), new=b.get(k, "<missing>")))
            else:
                json_leaves(a[k], b[k], f"{path}/{k}", tol, numeric, strings)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            numeric.append(dict(path=path, problem=f"length {len(a)} vs {len(b)}"))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            json_leaves(x, y, f"{path}[{i}]", tol, numeric, strings)
    elif (isinstance(a, (int, float)) and isinstance(b, (int, float))
          and not isinstance(a, bool) and not isinstance(b, bool)):
        if a != b and not (np.isnan(a) and np.isnan(b)) and abs(a - b) > tol:
            numeric.append(dict(path=path, old=a, new=b, abs_diff=abs(a - b)))
        elif a != b:
            numeric.append(dict(path=path, old=a, new=b, abs_diff=abs(a - b), within_tolerance=True))
    elif a != b:
        (numeric if isinstance(a, bool) or isinstance(b, bool) else strings).append(
            dict(path=path, old=a, new=b))


def compare_file(old, new, tol):
    suffix = old.suffix
    if suffix == ".npy":
        a, b = np.load(old, mmap_mode="r", allow_pickle=False), np.load(new, mmap_mode="r", allow_pickle=False)
        return dict(kind="npy", **arrays_equal(a, b))
    if suffix == ".npz":
        with np.load(old, allow_pickle=False) as za, np.load(new, allow_pickle=False) as zb:
            if set(za.files) != set(zb.files):
                return dict(kind="npz", equal=False, reason=f"keys {sorted(za.files)} vs {sorted(zb.files)}")
            fields = {k: arrays_equal(za[k], zb[k]) for k in za.files}
        return dict(kind="npz", equal=all(f["equal"] for f in fields.values()),
                    max_abs_diff=max((f.get("max_abs_diff") or 0.0) for f in fields.values()),
                    unequal_fields={k: v for k, v in fields.items() if not v["equal"]})
    if suffix == ".json":
        a, b = json.loads(old.read_text()), json.loads(new.read_text())
        if (isinstance(a, list) and isinstance(b, list) and all(isinstance(e, dict) and "sha256" in e for e in a + b)):
            # Provenance list of consumed files: compare by content hash, report path-level changes.
            sa, sb = {e["sha256"]: e["path"] for e in a}, {e["sha256"]: e["path"] for e in b}
            return dict(kind="provenance", equal=a == b, numeric_equal=True, same_bytes=len(set(sa) & set(sb)),
                        only_old=sorted(sa[h] for h in set(sa) - set(sb)),
                        only_new=sorted(sb[h] for h in set(sb) - set(sa)))
        if old.name == "manifest.json":
            outputs_a = {e["path"]: e["sha256"] for e in a.get("outputs", [])}
            outputs_b = {e["path"]: e["sha256"] for e in b.get("outputs", [])}
            a = dict(stage=a.get("stage"), config=a.get("config"))
            b = dict(stage=b.get("stage"), config=b.get("config"))
            extra = dict(same_output_names=set(outputs_a) == set(outputs_b),
                         outputs_with_identical_sha=sorted(p for p in outputs_a if outputs_b.get(p) == outputs_a[p]),
                         outputs_with_different_sha=sorted(p for p in outputs_a
                                                           if p in outputs_b and outputs_b[p] != outputs_a[p]))
        else:
            extra = {}
        numeric, strings = [], []
        json_leaves(a, b, "", tol, numeric, strings)
        bad = [n for n in numeric if not n.get("within_tolerance")]
        return dict(kind="manifest" if old.name == "manifest.json" else "json", equal=not numeric and not strings,
                    numeric_equal=not bad, numeric_mismatches=bad[:20], numeric_mismatch_count=len(bad),
                    numeric_within_tol=len(numeric) - len(bad), string_differences=strings[:20],
                    string_difference_count=len(strings), **extra)
    if suffix == ".csv":
        if old.read_bytes() == new.read_bytes():
            return dict(kind="csv", equal=True, numeric_equal=True)
        ra, rb = list(csv.reader(old.open())), list(csv.reader(new.open()))
        numeric = []
        if len(ra) != len(rb):
            return dict(kind="csv", equal=False, numeric_equal=False, reason="row count")
        for i, (x, y) in enumerate(zip(ra, rb)):
            for j, (u, v) in enumerate(zip(x, y)):
                if u != v:
                    try:
                        if abs(float(u) - float(v)) > tol:
                            numeric.append(dict(row=i, col=j, old=u, new=v))
                    except ValueError:
                        numeric.append(dict(row=i, col=j, old=u, new=v))
        return dict(kind="csv", equal=False, numeric_equal=not numeric, mismatches=numeric[:20])
    return dict(kind=suffix.lstrip(".") or "file", equal=sha(old) == sha(new), informational=True)


def compare(old_root, new_root, tol=1e-12, skip=()):
    old_root, new_root = Path(old_root), Path(new_root)
    files_old = {p.relative_to(old_root) for p in old_root.rglob("*") if p.is_file()}
    files_new = {p.relative_to(new_root) for p in new_root.rglob("*") if p.is_file()}
    results = {}
    for rel in sorted(files_old & files_new):
        if any(str(rel).startswith(s) for s in skip):
            continue
        results[str(rel)] = compare_file(old_root / rel, new_root / rel, tol)
    numeric_kinds = ("npy", "npz")
    arrays = {k: v for k, v in results.items() if v["kind"] in numeric_kinds}
    texts = {k: v for k, v in results.items() if v["kind"] in ("json", "manifest", "csv", "provenance")}
    summary = dict(
        old_root=str(old_root), new_root=str(new_root), tolerance=tol,
        only_in_old=sorted(map(str, files_old - files_new)), only_in_new=sorted(map(str, files_new - files_old)),
        compared_files=len(results), arrays=len(arrays),
        arrays_bitwise_identical=sum(v["equal"] for v in arrays.values()),
        arrays_max_abs_diff=max([v.get("max_abs_diff") or 0.0 for v in arrays.values()] or [0.0]),
        json_csv_files=len(texts), json_csv_numeric_equal=sum(v.get("numeric_equal", v["equal"]) for v in texts.values()),
        json_csv_fully_identical=sum(v["equal"] for v in texts.values()),
        other_files={k: v["equal"] for k, v in results.items() if v.get("informational")},
    )
    summary["all_numerics_equal_within_tol"] = (
        summary["arrays_max_abs_diff"] <= tol and summary["json_csv_numeric_equal"] == len(texts))
    summary["all_arrays_bitwise"] = summary["arrays_bitwise_identical"] == len(arrays)
    return dict(summary=summary, files=results)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("old_root", type=Path)
    p.add_argument("new_root", type=Path)
    p.add_argument("name")
    p.add_argument("--tol", type=float, default=1e-12)
    p.add_argument("--skip", nargs="*", default=())
    args = p.parse_args()
    record = compare(args.old_root, args.new_root, args.tol, args.skip)
    write(f"products-{args.name}.json", record)
    print(json.dumps(record["summary"], indent=1))


if __name__ == "__main__":
    main()
