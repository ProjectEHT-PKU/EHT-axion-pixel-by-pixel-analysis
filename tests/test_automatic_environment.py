"""Automatic dotenv is deterministic data parsing, not shell evaluation."""
from argparse import Namespace
import os
from pathlib import Path
import subprocess
import sys

import pytest

import _environment
import _producer

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


def test_producer_overlay_only_loaded_for_producer(project):
    (project / ".env").write_text("EHT_OUTPUT_ROOT=../shared\nEHT_RUN_ID=shared\n")
    (project / ".env.raptor").write_text("EHT_RUN_ID=producer\nEHT_GRMHD_ROOT=fluid\n")
    assert _environment.resolve_environment(request())["EHT_RUN_ID"] == "shared"
    prod = _environment.resolve_environment(request(), producer=True)
    assert prod["EHT_RUN_ID"] == "producer"
    assert prod["EHT_GRMHD_ROOT"] == str(project / "fluid")
    assert prod["EHT_OUTPUT_ROOT"] == str(project / "../shared")


def test_cli_then_exports_then_files(project, monkeypatch):
    (project / ".env").write_text("EHT_RUN_ID=file\nEHT_GRMHD_ROOT=files\n")
    assert _environment.resolve_environment(request())["EHT_RUN_ID"] == "file"
    monkeypatch.setenv("EHT_RUN_ID", "export")
    assert _environment.resolve_environment(request())["EHT_RUN_ID"] == "export"


def test_scientific_or_removed_keys_are_rejected(project):
    for key in ("EHT_CONFIG", "EHT_PRODUCTION_CONFIG", "EHT_TRUST_PICKLE"):
        (project / ".env").write_text(f"{key}=value\n")
        with pytest.raises(ValueError, match="Unknown dotenv key"):
            _environment.resolve_environment(request())


def test_file_changes_reloaded_without_mutating_process(project):
    path = project / ".env"
    path.write_text("EHT_RUN_ID=first\n")
    assert _environment.resolve_environment(request())["EHT_RUN_ID"] == "first"
    path.write_text("EHT_RUN_ID=second\n")
    assert _environment.resolve_environment(request())["EHT_RUN_ID"] == "second"
    assert "EHT_RUN_ID" not in os.environ


def test_custom_file_replaces_defaults_and_can_be_disabled(project):
    (project / ".env").write_text("EHT_RUN_ID=default\n")
    (project / ".env.raptor").write_text("EHT_RUN_ID=overlay\n")
    custom = project / "custom.env"
    custom.write_text("EHT_RUN_ID=custom\n")
    assert _environment.resolve_environment(request(env_file=custom), producer=True)["EHT_RUN_ID"] == "custom"
    assert "EHT_RUN_ID" not in _environment.resolve_environment(request(no_env_file=True), producer=True)
    with pytest.raises(ValueError, match="does not exist"):
        _environment.resolve_environment(request(env_file=project / "missing.env"))


def test_no_interpolation_or_command_execution(project):
    path = project / ".env"
    marker = project / "side-effect"
    path.write_text(f'EHT_RUN_ID="$(touch {marker})"\nEHT_GRMHD_ROOT=${{HOME}}/fluid\n')
    values = _environment.resolve_environment(request())
    assert values["EHT_RUN_ID"] == f"$(touch {marker})"
    assert values["EHT_GRMHD_ROOT"] == str(project / "${HOME}/fluid")
    assert not marker.exists()


@pytest.mark.parametrize("content", [
    "touch side-effect\n", "EHT_RUN_ID\n", "EHT_RUN_ID=one\nEHT_RUN_ID=two\n",
    "PATH=/unexpected\n", 'EHT_RUN_ID="unfinished\n', "EHT_RUN_ID=x\nrm -rf /\n",
])
def test_malformed_shell_and_unknown_keys_rejected(project, content):
    (project / ".env").write_text(content)
    with pytest.raises(ValueError):
        _environment.resolve_environment(request())


def test_producer_entrypoint_uses_automatic_files(project, monkeypatch):
    core = project / "core"
    core.mkdir()
    (project / ".env").write_text("EHT_OUTPUT_ROOT=outputs\n")
    (project / ".env.raptor").write_text("EHT_RAPTOR_SOURCE_ROOT=core\nEHT_RUN_ID=automatic\n")
    monkeypatch.setattr(_producer, "ROOT", ROOT)
    args = _producer.options("build", [])
    assert args.run_root == project / "outputs/production/automatic"
    assert args.source_root == core
    assert args.config == ROOT / "experiment.yaml"
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
