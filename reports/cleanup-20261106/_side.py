"""Run a snippet against the OLD code (main, exported read-only) or the NEW worktree.

Both trees define the package ``eht_axion``, so each side runs in its own
interpreter. The snippet must assign a JSON-serialisable ``result``.
Usage from evidence scripts: ``run('main', code)`` / ``run('new', code)``.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

NEW = Path(__file__).resolve().parents[2]
MAIN = Path(os.environ.get("CLEANUP_MAIN_SRC", "/tmp/cleanupopus/main-src"))
MAIN_COMMIT = "5b1807bc802903630a1704ecef395a58d0a1bdfd"
VERIFY = Path("/home/canoziia/GRMHD_DATA/analysis/cleanup-verify-opus-20261106")


def run(side, code, python=None):
    root = {"main": MAIN, "new": NEW}[side]
    if not (root / "eht_axion").is_dir():
        raise FileNotFoundError(f"{side} source missing: {root} (git archive {MAIN_COMMIT})")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "result.json"
        wrapper = (
            "import sys, json\n"
            "sys.dont_write_bytecode = True\n"
            f"sys.path[:0] = [{str(root)!r}, {str(root / 'scripts')!r}]\n"
            f"{code}\n"
            f"open({str(out)!r}, 'w').write(json.dumps(result, default=repr))\n"
        )
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2")
        subprocess.run([python or sys.executable, "-c", wrapper], check=True, cwd=tmp, env=env)
        return json.loads(out.read_text())


def diff(a, b, path="", tol=0.0, out=None):
    """List leaf differences between two JSON-like values (exact unless tol > 0)."""
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b), key=str):
            if key not in a or key not in b:
                out.append(dict(path=f"{path}/{key}", old=a.get(key, "<missing>"), new=b.get(key, "<missing>")))
            else:
                diff(a[key], b[key], f"{path}/{key}", tol, out)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for k, (x, y) in enumerate(zip(a, b)):
            diff(x, y, f"{path}[{k}]", tol, out)
    elif (isinstance(a, (int, float)) and isinstance(b, (int, float))
          and not isinstance(a, bool) and not isinstance(b, bool)):
        if not (a == b or (tol and abs(a - b) <= tol)):
            out.append(dict(path=path, old=a, new=b))
    elif a != b:
        out.append(dict(path=path, old=a, new=b))
    return out


def write(name, record):
    """Write evidence JSON to reports/cleanup-20261106 and to the verify root."""
    text = json.dumps(record, indent=2, default=repr) + "\n"
    (NEW / "reports/cleanup-20261106" / name).write_text(text)
    VERIFY.mkdir(parents=True, exist_ok=True)
    (VERIFY / "evidence").mkdir(exist_ok=True)
    (VERIFY / "evidence" / name).write_text(text)
    return text
