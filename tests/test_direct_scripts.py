"""Direct scripts work without an editable distribution or source PYTHONPATH."""

import os
from pathlib import Path
import subprocess
import sys
import tomllib

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "run_all",
    "process_raptor",
    "process_eht",
    "estimate_background",
    "build_observations",
    "calculate_limits",
    "plot_limits",
    "plot_saved_limits",
    "verify_products",
    "migrate_data",
    "plot_diagnostics",
    "check_convergence",
    "raptor_tools",
    "transport_diagnostics",
    "process_eht_ensemble",
    "compare_original_limits",
    "check_script_pipeline_legacy",
]


@pytest.mark.parametrize("name", SCRIPTS)
def test_script_help_outside_repository_without_installation(name, tmp_path):
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / f"{name}.py"), "--help"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


def test_workspace_configuration_is_only_in_pixi():
    for name in ("pyproject.toml", "pytest.ini", "ruff.toml", ".ruff.toml"):
        assert not (ROOT / name).exists()
    pixi = tomllib.loads((ROOT / "pixi.toml").read_text())
    assert "pythonpath=. scripts" in pixi["tasks"]["test"]["cmd"]
    assert "testpaths=tests" in pixi["tasks"]["test"]["cmd"]
    assert "--isolated" in pixi["tasks"]["lint"]["cmd"]
    assert "--select E,F,W" in pixi["tasks"]["lint"]["cmd"]
    assert not pixi.get("pypi-dependencies")
    assert pixi["tasks"]["pipeline"]["cmd"] == "python scripts/run_all.py"
    assert not (ROOT / "eht_axion/pipeline/cli.py").exists()
    assert not (ROOT / "eht_axion/pipeline/__main__.py").exists()
