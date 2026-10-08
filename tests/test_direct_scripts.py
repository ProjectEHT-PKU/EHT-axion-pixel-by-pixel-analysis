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


def test_workspace_is_tool_config_not_package():
    tools = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert "build-system" not in tools
    assert "project" not in tools
    pixi = tomllib.loads((ROOT / "pixi.toml").read_text())
    assert not pixi.get("pypi-dependencies")
    assert pixi["tasks"]["pipeline"]["cmd"] == "python scripts/run_all.py"
    assert not (ROOT / "eht_axion/pipeline/cli.py").exists()
    assert not (ROOT / "eht_axion/pipeline/__main__.py").exists()
