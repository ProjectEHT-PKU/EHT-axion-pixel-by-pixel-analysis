"""Automatic dotenv is deterministic data parsing, not shell evaluation."""
from argparse import Namespace
import os
from pathlib import Path
import subprocess
import sys

import pytest

from eht_axion import environment
from eht_axion.raptor import options as raptor_options

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setattr(environment, "ROOT", tmp_path)
    for name in list(os.environ):
        if name.startswith("EHT_"):
            monkeypatch.delenv(name)
    return tmp_path


def request(**kwargs):
    return Namespace(env_file=kwargs.get("env_file"), no_env_file=kwargs.get("no_env_file", False))


def test_single_dotenv_with_three_roots_anchored_to_the_project(project):
    (project / ".env").write_text("EHT_RAPTOR_SOURCE_ROOT=../raptor\nEHT_INPUT_ROOT=inputs\n"
                                  "EHT_DATA_ROOT=/data\nEHT_RAPTOR_THREADS=1\n")
    (project / ".env.raptor").write_text("EHT_RAPTOR_THREADS=4\n")  # no longer read
    values = environment.resolve_environment(request())
    assert values["EHT_RAPTOR_THREADS"] == "1"
    assert values["EHT_RAPTOR_SOURCE_ROOT"] == str(project / "../raptor")
    assert values["EHT_INPUT_ROOT"] == str(project / "inputs")
    assert values["EHT_DATA_ROOT"] == "/data"
    assert environment.PATH_KEYS == {"EHT_RAPTOR_SOURCE_ROOT", "EHT_INPUT_ROOT", "EHT_DATA_ROOT"}


@pytest.mark.parametrize("key", ["EHT_RAPTOR_RUN_ROOT", "EHT_GRMHD_ROOT", "EHT_RADIAL_LOW", "EHT_RADIAL_HIGH",
                                 "EHT_IMAGING_INPUT_ROOT", "EHT_IMAGING_INPUT_RECEIPT",
                                 "EHT_IMAGING_FEED_ROTATION", "EHT_IMAGING_OFFICIAL_SCRIPT"])
def test_obsolete_location_keys_are_rejected(project, key):
    (project / ".env").write_text(f"{key}=/somewhere\n")
    with pytest.raises(ValueError, match="Unknown dotenv key"):
        environment.resolve_environment(request())


def test_cli_then_exports_then_files(project, monkeypatch):
    (project / ".env").write_text("EHT_RAPTOR_CPU_LIST=file\nEHT_INPUT_ROOT=files\n")
    assert environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "file"
    monkeypatch.setenv("EHT_RAPTOR_CPU_LIST", "export")
    assert environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "export"


def test_file_changes_reloaded_without_mutating_process(project):
    path = project / ".env"
    path.write_text("EHT_RAPTOR_CPU_LIST=first\n")
    assert environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "first"
    path.write_text("EHT_RAPTOR_CPU_LIST=second\n")
    assert environment.resolve_environment(request())["EHT_RAPTOR_CPU_LIST"] == "second"
    assert "EHT_RAPTOR_CPU_LIST" not in os.environ


def test_custom_file_replaces_defaults_and_can_be_disabled(project):
    (project / ".env").write_text("EHT_RAPTOR_CPU_LIST=default\n")
    custom = project / "custom.env"
    custom.write_text("EHT_RAPTOR_CPU_LIST=custom\n")
    key = "EHT_RAPTOR_CPU_LIST"
    assert environment.resolve_environment(request(env_file=custom))[key] == "custom"
    assert key not in environment.resolve_environment(request(no_env_file=True))
    with pytest.raises(ValueError, match="does not exist"):
        environment.resolve_environment(request(env_file=project / "missing.env"))


def test_no_interpolation_or_command_execution(project):
    path = project / ".env"
    marker = project / "side-effect"
    path.write_text(f'EHT_RAPTOR_CPU_LIST="$(touch {marker})"\nEHT_INPUT_ROOT=${{HOME}}/inputs\n')
    values = environment.resolve_environment(request())
    assert values["EHT_RAPTOR_CPU_LIST"] == f"$(touch {marker})"
    assert values["EHT_INPUT_ROOT"] == str(project / "${HOME}/inputs")
    assert not marker.exists()


@pytest.mark.parametrize("content", [
    "touch side-effect\n", "EHT_RAPTOR_CPU_LIST\n", "EHT_RAPTOR_CPU_LIST=one\nEHT_RAPTOR_CPU_LIST=two\n",
    "PATH=/unexpected\n", 'EHT_RAPTOR_CPU_LIST="unfinished\n', "EHT_RAPTOR_CPU_LIST=x\nrm -rf /\n",
])
def test_malformed_shell_and_unknown_keys_rejected(project, content):
    (project / ".env").write_text(content)
    with pytest.raises(ValueError):
        environment.resolve_environment(request())


def test_raptor_entrypoint_uses_automatic_files(project, monkeypatch):
    core = project / "core"
    core.mkdir()
    (project / ".env").write_text("EHT_DATA_ROOT=outputs\nEHT_RAPTOR_SOURCE_ROOT=core\n")
    monkeypatch.setattr(raptor_options, "ROOT", ROOT)
    args = raptor_options.options("build", [])
    assert args.run_root == project / "outputs/raptor/runs/main"
    assert args.source_root == core
    assert args.config == ROOT / "configs" / "experiment.yaml"
    assert args.config_data["shared"]["mass_solar"] == 6.5e9


def test_build_never_needs_or_resolves_generation_inputs(project, monkeypatch):
    core = project / "core"
    core.mkdir()
    # The input root does not exist: build must neither require nor resolve it.
    (project / ".env").write_text("EHT_RAPTOR_SOURCE_ROOT=core\nEHT_INPUT_ROOT=missing\nEHT_DATA_ROOT=d\n")
    monkeypatch.setattr(raptor_options, "ROOT", ROOT)
    args = raptor_options.options("build", [])
    assert not hasattr(args, "fluid_root") and not hasattr(args, "radial_low")
    with pytest.raises(SystemExit):
        raptor_options.options("build", ["--fluid-root", str(project)])


def test_generate_inputs_from_input_root_and_explicit_overrides(project, monkeypatch):
    core = project / "core"
    core.mkdir()
    inputs = project / "inputs"
    for key, relative in environment.INPUT_LAYOUT.items():
        if key in ("grmhd", "eht.uvfits"):
            (inputs / relative).mkdir(parents=True, exist_ok=True)
        else:
            (inputs / relative).parent.mkdir(parents=True, exist_ok=True)
            (inputs / relative).write_text("x")
    other = project / "elsewhere" / "other.csv"
    other.parent.mkdir()
    other.write_text("x")
    (project / ".env").write_text("EHT_RAPTOR_SOURCE_ROOT=core\nEHT_INPUT_ROOT=inputs\nEHT_DATA_ROOT=d\n")
    monkeypatch.setattr(raptor_options, "ROOT", ROOT)
    args = raptor_options.options("generate", ["--plan", "--radial-high", str(other)])
    assert args.fluid_root == inputs / "grmhd"
    assert args.radial_low == inputs / "radial/result_a=0.9375_l1m1.txt"
    assert args.radial_high == other
    assert args.run_root == project / "d/raptor/runs/main"
    with pytest.raises(SystemExit):  # neither explicit inputs nor an input root
        raptor_options.options("generate", ["--plan", "--no-env-file", "--source-root", str(core),
                                            "--data-root", str(project / "d")])
    with pytest.raises(SystemExit):  # the input root lacks the snapshots
        raptor_options.options("generate", ["--plan", "--input-root", str(project / "missing")])


def test_scripts_help_accepts_env_switches_without_reading_bad_file(project):
    bad = project / "bad.env"
    bad.write_text("touch side-effect\n")
    for script in ("scripts/raptor/build.py", "scripts/raptor/generate.py", "scripts/eht/reconstruct.py"):
        result = subprocess.run([sys.executable, str(ROOT / script),
                                 "--env-file", str(bad), "--help"], cwd=project,
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stderr
        assert "--no-env-file" in result.stdout
