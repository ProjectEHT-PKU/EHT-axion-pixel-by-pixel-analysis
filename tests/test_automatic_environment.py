"""Automatic dotenv is deterministic data parsing, not shell evaluation."""
from argparse import Namespace
import os
from pathlib import Path
import subprocess
import sys

import pytest

import _environment
import _raptor_options

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setattr(_environment, "ROOT", tmp_path)
    for name in list(os.environ):
        if name.startswith("EHT_"):
            monkeypatch.delenv(name)
    return tmp_path


def request(**kwargs):
    return Namespace(env_file=kwargs.get("env_file"), no_env_file=kwargs.get("no_env_file", False))


def test_raptor_overlay_only_loaded_for_raptor_scripts(project):
    (project / ".env").write_text("EHT_RAPTOR_RUN_ROOT=../shared\nEHT_RAPTOR_THREADS=1\n")
    (project / ".env.raptor").write_text("EHT_RAPTOR_THREADS=4\nEHT_GRMHD_ROOT=fluid\n")
    assert _environment.resolve_environment(request())["EHT_RAPTOR_THREADS"] == "1"
    values = _environment.resolve_environment(request(), raptor=True)
    assert values["EHT_RAPTOR_THREADS"] == "4"
    assert values["EHT_GRMHD_ROOT"] == str(project / "fluid")
    assert values["EHT_RAPTOR_RUN_ROOT"] == str(project / "../shared")


def test_cli_then_exports_then_files(project, monkeypatch):
    (project / ".env").write_text("EHT_RAPTOR_CPU_LIST=file\nEHT_GRMHD_ROOT=files\n")
    assert _environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "file"
    monkeypatch.setenv("EHT_RAPTOR_CPU_LIST", "export")
    assert _environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "export"


def test_file_changes_reloaded_without_mutating_process(project):
    path = project / ".env"
    path.write_text("EHT_RAPTOR_CPU_LIST=first\n")
    assert _environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "first"
    path.write_text("EHT_RAPTOR_CPU_LIST=second\n")
    assert _environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "second"
    assert "EHT_RAPTOR_CPU_LIST" not in os.environ


def test_custom_file_replaces_defaults_and_can_be_disabled(project):
    (project / ".env").write_text("EHT_RAPTOR_CPU_LIST=default\n")
    (project / ".env.raptor").write_text("EHT_RAPTOR_CPU_LIST=overlay\n")
    custom = project / "custom.env"
    custom.write_text("EHT_RAPTOR_CPU_LIST=custom\n")
    key = "EHT_RAPTOR_CPU_LIST"
    assert _environment.resolve_environment(request(env_file=custom), raptor=True)[key] == "custom"
    assert key not in _environment.resolve_environment(request(no_env_file=True), raptor=True)
    with pytest.raises(ValueError, match="does not exist"):
        _environment.resolve_environment(request(env_file=project / "missing.env"))


def test_no_interpolation_or_command_execution(project):
    path = project / ".env"
    marker = project / "side-effect"
    path.write_text(f'EHT_RAPTOR_CPU_LIST="$(touch {marker})"\nEHT_GRMHD_ROOT=${{HOME}}/fluid\n')
    values = _environment.resolve_environment(request())
    assert values["EHT_RAPTOR_CPU_LIST"] == f"$(touch {marker})"
    assert values["EHT_GRMHD_ROOT"] == str(project / "${HOME}/fluid")
    assert not marker.exists()


@pytest.mark.parametrize("content", [
    "touch side-effect\n", "EHT_RAPTOR_CPU_LIST\n", "EHT_RAPTOR_CPU_LIST=one\nEHT_RAPTOR_CPU_LIST=two\n",
    "PATH=/unexpected\n", 'EHT_RAPTOR_CPU_LIST="unfinished\n', "EHT_RAPTOR_CPU_LIST=x\nrm -rf /\n",
])
def test_malformed_shell_and_unknown_keys_rejected(project, content):
    (project / ".env").write_text(content)
    with pytest.raises(ValueError):
        _environment.resolve_environment(request())


def test_raptor_entrypoint_uses_automatic_files(project, monkeypatch):
    core = project / "core"
    core.mkdir()
    (project / ".env").write_text("EHT_RAPTOR_RUN_ROOT=outputs/run\n")
    (project / ".env.raptor").write_text("EHT_RAPTOR_SOURCE_ROOT=core\n")
    monkeypatch.setattr(_raptor_options, "ROOT", ROOT)
    args = _raptor_options.options("build", [])
    assert args.run_root == project / "outputs/run"
    assert args.source_root == core
    assert args.config == ROOT / "configs" / "experiment.yaml"
    assert args.config_data["shared"]["mass_solar"] == 6.5e9


def test_scripts_help_accepts_env_switches_without_reading_bad_file(project):
    bad = project / "bad.env"
    bad.write_text("touch side-effect\n")
    for script in ("raptor_build.py", "raptor_generate.py", "eht_reconstruct.py"):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / script),
                                 "--env-file", str(bad), "--help"], cwd=project,
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stderr
        assert "--no-env-file" in result.stdout
