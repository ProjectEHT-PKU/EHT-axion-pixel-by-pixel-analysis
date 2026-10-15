"""CLI and exports retain priority above automatically parsed dotenv defaults."""

from pathlib import Path
import subprocess
import sys
import os

import pytest
import _stage
import _environment


@pytest.fixture(autouse=True)
def clean_environment(monkeypatch, tmp_path):
    # Never let tests consume the user's actual local dotenv or science outputs.
    monkeypatch.setattr(_environment, "ROOT", tmp_path)
    for name in _stage.ENVIRONMENT_OPTIONS.values():
        monkeypatch.delenv(name, raising=False)


def configured(monkeypatch):
    monkeypatch.setenv("EHT_CONFIG", "configs/research.yaml")
    monkeypatch.setenv("EHT_INPUT_ROOT", "data/inputs")
    monkeypatch.setenv("EHT_OUTPUT_ROOT", "../runs/example")


def test_exported_values_are_used(monkeypatch):
    configured(monkeypatch)
    monkeypatch.setenv("EHT_TRUST_PICKLE", "TRUE")
    args = _stage.parse_options("eht", [])
    assert args.config == Path("configs/research.yaml")
    assert args.input_root == Path("data/inputs")
    assert args.output_root == Path("../runs/example")
    assert args.trust_pickle is True
    assert args.rebuild is False
    assert args.reuse_imported is False


def test_cli_overrides_paths_and_bools_even_when_environment_invalid(monkeypatch):
    configured(monkeypatch)
    monkeypatch.setenv("EHT_TRUST_PICKLE", "not-a-boolean")
    monkeypatch.setenv("EHT_REBUILD", "true")
    monkeypatch.setenv("EHT_REUSE_IMPORTED", "yes")
    args = _stage.parse_options(
        "eht",
        [
            "--config",
            "other.yaml",
            "--input-root",
            "other data",
            "--output-root",
            "other run",
            "--trust-pickle",
            "--no-rebuild",
            "--no-reuse-imported",
        ],
    )
    assert args.config == Path("other.yaml")
    assert args.input_root == Path("other data")
    assert args.output_root == Path("other run")
    assert args.trust_pickle and not args.rebuild and not args.reuse_imported


@pytest.mark.parametrize(
    "value,expected",
    [
        ("true", True),
        ("1", True),
        ("yes", True),
        ("on", True),
        (" FALSE ", False),
        ("0", False),
        ("no", False),
        ("off", False),
        ("", False),
    ],
)
def test_boolean_environment_spellings(monkeypatch, value, expected):
    monkeypatch.setenv("EHT_TRUST_PICKLE", value)
    assert _stage.environment_boolean("EHT_TRUST_PICKLE") is expected


def test_invalid_boolean_reports_environment_name(monkeypatch, capsys):
    configured(monkeypatch)
    monkeypatch.setenv("EHT_TRUST_PICKLE", "maybe")
    with pytest.raises(SystemExit):
        _stage.parse_options("eht", [])
    assert "EHT_TRUST_PICKLE" in capsys.readouterr().err


def test_required_environment_and_no_pickle_trust_by_default(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        _stage.parse_options("all", [])
    assert "EHT_CONFIG" in capsys.readouterr().err
    configured(monkeypatch)
    assert not _stage.parse_options("all", []).trust_pickle
    monkeypatch.delenv("EHT_INPUT_ROOT")
    assert _stage.parse_options("inference", []).input_root is None
    with pytest.raises(SystemExit):
        _stage.parse_options("raptor", [])
    assert "EHT_INPUT_ROOT" in capsys.readouterr().err


def test_shell_exports_reviewed_dotenv(tmp_path):
    import json
    import shutil

    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("Bash shell loading example is tested only where Bash exists")
    env_file = tmp_path / "settings with spaces.env"
    env_file.write_text(
        'EHT_CONFIG="science config.yaml"\nEHT_INPUT_ROOT="input data"\n'
        'EHT_OUTPUT_ROOT="run output"\nEHT_TRUST_PICKLE=false\n'
    )
    script_dir = Path(__file__).resolve().parents[1] / "scripts"
    program = (
        "import sys,json;sys.path.insert(0,sys.argv[1]);from _stage import parse_options;"
        'a=parse_options("raptor",[]);'
        "print(json.dumps(dict(config=str(a.config),input=str(a.input_root),"
        "output=str(a.output_root),trust=a.trust_pickle)))"
    )
    result = subprocess.run(
        [
            bash,
            "-c",
            'set -a; source "$1"; set +a; exec "$2" -c "$3" "$4"',
            "env-test",
            str(env_file),
            sys.executable,
            program,
            str(script_dir),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == dict(
        config="science config.yaml", input="input data", output="run output", trust=False
    )


def test_local_dotenv_is_not_automatically_executed(tmp_path):
    (tmp_path / ".env").write_text("touch unexpected-side-effect\nEHT_CONFIG=do-not-load.yaml\n")
    root = Path(__file__).resolve().parents[1]
    environment = dict(os.environ)
    for name in _stage.ENVIRONMENT_OPTIONS.values():
        environment.pop(name, None)
    environment.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, str(root / "scripts/run_all.py"), "--no-env-file"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode != 0
    assert "EHT_CONFIG" in result.stderr
    assert not (tmp_path / "unexpected-side-effect").exists()
