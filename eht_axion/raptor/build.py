"""Transactional builds using RAPTOR's official setup and make source list.

The returned directory contains RAPTOR, the generated cwd headers/grid inputs,
source/ (private, tracked Git bytes only), compile.log, macros.txt and receipt.json.
No simulation data or axion tables are build inputs. Run with the raptor-build
Pixi environment; this module never installs into the default environment.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

from ..experiment import normalize_repository
from ..provenance import sha256, stage_directory

_FIELDS = {"model", "metric", "integrator", "grid", "polarized", "sfc", "magnetic_polarity",
           "sigma_cut", "thetae_min", "thetae_max", "uniform_output"}


def _settings(settings):
    """Validate the complete raptor.build section (no defaults) against setup.sh options."""
    if not isinstance(settings, dict) or set(settings) != _FIELDS:
        raise ValueError(f"raptor.build requires exactly {sorted(_FIELDS)}")
    result = dict(settings)
    choices = {
        "model": {"bhac", "harm3d", "hamr"}, "metric": {"mks", "cks"},
        "integrator": {"rk4", "rk2", "rk45", "ver"}, "grid": {"SMR", "AMR"},
        "magnetic_polarity": {"plus", "minus"},
    }
    for key, allowed in choices.items():
        if not isinstance(result[key], str) or result[key] not in allowed:
            raise ValueError(f"Invalid raptor.build.{key}: {result[key]!r}")
    if result["model"] != "bhac" and result["metric"] == "cks":
        raise ValueError("CKS requires the BHAC model")
    for key in ("polarized", "sfc", "uniform_output"):
        if type(result[key]) is not bool:
            raise ValueError(f"raptor.build.{key} must be boolean")
    for key in ("sigma_cut", "thetae_min", "thetae_max"):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"raptor.build.{key} must be a finite positive number")
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"raptor.build.{key} must be a finite positive number")
        result[key] = float(value)
    if result["thetae_min"] >= result["thetae_max"]:
        raise ValueError("thetae_min must be less than thetae_max")
    return result


def _run(command, *, cwd=None, env=None, **kwargs):
    return subprocess.run(command, cwd=cwd, env=env, check=True, **kwargs)


def _git(root, *args):
    return _run(["git", "-C", str(root), *args], stdout=subprocess.PIPE).stdout


def _clean_source(root, source):
    """The checkout is clean, HEAD is source['commit'] and a remote points to source['repository']."""
    if _git(root, "rev-parse", "--show-toplevel").decode().strip() != str(root):
        raise ValueError("source_root must be the RAPTOR Git repository root")
    commit = _git(root, "rev-parse", "HEAD").decode().strip()
    if commit != source["commit"]:
        raise ValueError(f"Wrong RAPTOR commit: expected {source['commit']}, found {commit}")
    if _git(root, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Refusing dirty RAPTOR source tree (including untracked files)")
    remotes = subprocess.run(["git", "-C", str(root), "config", "--get-regexp", r"^remote\..*\.url$"],
                             stdout=subprocess.PIPE, text=True, check=False).stdout.split()[1::2]
    if source["repository"] not in {normalize_repository(url) for url in remotes}:
        raise ValueError(f"No git remote of {root} points to {source['repository']} (remotes: {remotes})")
    return commit


def _tracked_source(root, commit):
    records = []
    for entry in _git(root, "ls-tree", "-rz", "--full-tree", commit).split(b"\0"):
        if not entry:
            continue
        metadata, name = entry.split(b"\t", 1)
        mode, kind, oid = metadata.decode().split()
        path = Path(os.fsdecode(name))
        if mode not in {"100644", "100755"} or kind != "blob":
            raise ValueError(f"Unsupported tracked source (symlink/submodule): {path}")
        if path.is_absolute() or ".." in path.parts or ".git" in path.parts:
            raise ValueError(f"Unsafe tracked source path: {path}")
        data = _git(root, "cat-file", "blob", oid)
        records.append((path, mode, data))
    required = {"setup.sh", "run/makefile", "src/main.c", "model/bhac/definitions.h"}
    if not required <= {str(path) for path, _, _ in records}:
        raise ValueError("Incomplete RAPTOR source tree")
    return records


def _toolchain():
    env = os.environ.copy()
    # Do not let a caller's make search paths, injected makefiles, or header paths
    # redirect the official source list outside the private snapshot.
    for name in ("VPATH", "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH", "OBJC_INCLUDE_PATH",
                 "MAKEFLAGS", "MFLAGS", "GNUMAKEFLAGS", "MAKEFILES", "RAPTOR", "CDPATH"):
        env.pop(name, None)
    compiler = shlex.split(env.get("HDF5_CC", "h5cc"))
    if len(compiler) != 1:
        raise ValueError("HDF5_CC must name one compiler executable, not a shell command")
    executable = shutil.which(compiler[0], path=env.get("PATH"))
    if not executable:
        raise FileNotFoundError(f"HDF5 compiler not found: {compiler[0]}; use pixi run -e raptor-build")
    compiler = [executable]
    tools = {}
    for name in ("bash", "make", "rsync"):
        tools[name] = shutil.which(name, path=env.get("PATH"))
        if not tools[name]:
            raise FileNotFoundError(f"Missing build tool: {name}")
    # Explicit CC takes precedence over ambient CC and make's built-in cc.
    env["CC"] = shlex.quote(executable)
    env["HDF5_CC"] = executable
    env["LC_ALL"] = "C"
    identity = {
        "compiler": executable,
        "compiler_version": _run(compiler + ["--version"], env=env, capture_output=True, text=True).stdout,
        "hdf5_expansion": _run(compiler + ["-show"], env=env, capture_output=True, text=True).stdout,
        "make_version": _run([tools["make"], "--version"], env=env, capture_output=True, text=True).stdout,
        "tools": tools,
        "environment_flags": {key: env.get(key, "") for key in ("CPPFLAGS", "CFLAGS", "LDFLAGS")},
    }
    return env, compiler, identity


def _arguments(settings):
    return [
        f"--code={settings['model']}", f"--metric={settings['metric']}",
        f"--int={settings['integrator']}", f"--grid={settings['grid']}",
        "--rad=" + ("pol" if settings["polarized"] else "unpol"),
        "--sf=" + ("sfc" if settings["sfc"] else "nosfc"),
        f"--bflip={settings['magnetic_polarity']}", f"--sigma-cut={settings['sigma_cut']}",
        f"--uniform={int(settings['uniform_output'])}",
        f"--thetae-min={settings['thetae_min']}", f"--thetae-max={settings['thetae_max']}",
    ]


def _verify_macros(text, settings):
    # Empty include guards must not consume the following macro's line.
    macros = dict(re.findall(r"^#define[ \t]+(\w+)[ \t]+([^\n]*)$", text, re.MULTILINE))
    metric = "CKS" if settings["metric"] == "cks" else (
        "MKSBHAC" if settings["model"] == "bhac" else "MKSHARM"
    )
    expected = {
        "metric": metric, "int_method": settings["integrator"].upper(),
        "AMR": int(settings["grid"] == "AMR"), "SMR": int(settings["grid"] == "SMR"),
        "POL": int(settings["polarized"]), "SFC": int(settings["sfc"]),
        "BPOL": settings["magnetic_polarity"].upper(), "SIGMA_CUT": settings["sigma_cut"],
        "UNIF": int(settings["uniform_output"]), "THETAE_MIN": settings["thetae_min"],
        "THETAE_MAX": settings["thetae_max"],
    }
    actual = {}
    for name, value in expected.items():
        raw = macros.get(name, "").strip().strip("()").strip()
        try:
            valid = raw == value if isinstance(value, str) else float(raw) == value
        except ValueError:
            valid = False
        if not valid:
            raise ValueError(f"Official setup macro mismatch: {name} expected {value!r}, got {raw!r}")
        actual[name] = raw
    return actual


def build_raptor(source_root: Path, destination: Path, settings: dict, source: dict) -> Path:
    """Build clean pinned core with official setup; return the verified stage directory.

    Inputs in the manifest never refer to temporary build paths. The RAPTOR source is
    identified by raptor.source (repository and commit; clean checkout at that commit with
    a remote pointing to that repository);
    the private source snapshot copied from that commit is a declared stage output. The
    built binary is hashed: jobs consume it as data. The commands run in the temporary
    stage directory.
    """
    settings = _settings(settings)
    source_root = Path(source_root).resolve(strict=True)
    destination = Path(destination).absolute()
    if destination.resolve() == source_root or destination.resolve().is_relative_to(source_root):
        raise ValueError("Build destination must be outside the core source repository")
    if source_root.is_relative_to(destination.resolve()):
        raise ValueError("Build destination must not contain the core source repository")
    if not isinstance(source, dict) or set(source) != {"repository", "commit"}:
        raise ValueError("raptor.source requires exactly repository and commit")
    commit = _clean_source(source_root, source)
    records = _tracked_source(source_root, commit)
    env, compiler, toolchain = _toolchain()
    config = {"settings": settings, "source_repository": source["repository"], "source_commit": commit,
              "toolchain": toolchain}
    with stage_directory(destination, "raptor:build", config) as work:
        if work is None:
            return destination
        snapshot = work / "source"
        for path, mode, data in records:
            target = snapshot / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            target.chmod(0o755 if mode == "100755" else 0o644)
        env["RAPTOR"] = str(snapshot)
        env["PWD"] = str(work)
        command = [toolchain["tools"]["bash"], str(snapshot / "setup.sh"), *_arguments(settings)]
        commands = [command]
        with (work / "compile.log").open("w") as log:
            log.write("COMMAND " + shlex.join(command) + "\n")
            log.flush()
            try:
                _run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as exc:
                log.flush()
                # stage_directory removes failed temporary trees. Keep useful
                # compiler diagnostics in the exception without publishing a stage.
                exc.add_note("Official setup compile.log tail:\n" + "\n".join(
                    (work / "compile.log").read_text(errors="replace").splitlines()[-40:]
                ))
                raise
        if not (work / "RAPTOR").is_file() or not os.access(work / "RAPTOR", os.X_OK):
            raise ValueError("Official setup did not produce an executable RAPTOR")
        for header in ("definitions.h", "model_definitions.h", "model_functions.h", "model_global_vars.h"):
            if not (work / header).is_file():
                raise ValueError(f"Official setup did not produce {header}")
        flags_command = [toolchain["tools"]["make"], "--no-print-directory", "-nB", "all"]
        commands.append(flags_command)
        flags = _run(flags_command, cwd=work, env=env, capture_output=True, text=True).stdout
        (work / "compiler-flags.txt").write_text(flags)
        # Preserve the official makefile's actual CPPFLAGS and CFLAGS, including
        # any -D options and language standard, when asking for active macros.
        compile_args = None
        for line in flags.splitlines():
            args = shlex.split(line)
            if args and Path(args[0]).name == Path(compiler[0]).name and "-c" in args:
                compile_args = args[1:]
                break
        if compile_args is None:
            raise ValueError("Official make did not report a compiler command")
        preprocess_args = []
        index = 0
        while index < len(compile_args):
            option = compile_args[index]
            if option in {"-o", "-MF", "-MT", "-MQ"}:
                index += 2
                continue
            if option not in {"-c", "-MMD", "-MD", "-MP"} and not option.endswith(".c"):
                preprocess_args.append(option)
            index += 1
        macro_command = compiler + preprocess_args + [
            "-dM", "-E", "-x", "c", "-include", "definitions.h", "-include", "model_definitions.h", "-",
        ]
        commands.append(macro_command)
        macros = _run(macro_command, cwd=work, env=env, input="", capture_output=True, text=True).stdout
        (work / "macros.txt").write_text(macros)
        actual = _verify_macros(macros, settings)
        _clean_source(source_root, source)
        receipt = {
            "source_repository": source["repository"], "source_commit": commit,
            "source_root": str(source_root),
            "source_snapshot": str(destination / "source"),
            "settings": settings, "actual_macros": actual, "toolchain": toolchain,
            "binary": {"path": str(destination / "RAPTOR"), "sha256": sha256(work / "RAPTOR")},
            "compiler_flags_file": str(destination / "compiler-flags.txt"),
            "commands": commands, "command_cwd": str(work),
            "runtime_cwd": str(destination), "tables_are_build_dependencies": False,
        }
        (work / "receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return destination
