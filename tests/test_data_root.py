"""Default roots below the data root: documented layout, precedence, and every script's defaults."""
from argparse import ArgumentParser
import os
from pathlib import Path
import re

import pytest

from eht_axion import environment as env

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def project(tmp_path, monkeypatch):
    """An isolated project root for dotenv files, without exported settings."""
    monkeypatch.setattr(env, "ROOT", tmp_path)
    for name in list(os.environ):
        if name.startswith("EHT_") or name == env.DATA_ROOT_KEY:
            monkeypatch.delenv(name)
    return tmp_path


def documented_layout(path):
    """Paths of the text tree in the 'Default layout' section, relative to <data root>."""
    text = path.read_text(encoding="utf-8")
    block = re.search(r"```text\n<data root>/\n(.*?)```", text, re.S)[1]
    stack, paths = [], set()
    for line in block.splitlines():
        name = line.strip().split()[0]
        depth = (len(line) - len(line.lstrip())) // 2 - 1
        stack[depth:] = [name.rstrip("/")]
        path = "/".join(stack)
        paths |= {path.rsplit("/", k)[0] for k in range(path.count("/") + 1)}
    return paths


@pytest.mark.parametrize("language", ["en", "zh"])
def test_layout_is_the_documented_default_layout(language):
    documented = documented_layout(ROOT / "docs" / language / "data-layout.md")
    for key, relative in env.LAYOUT.items():
        expected = relative.format(day="<day>")
        assert expected in documented, key


def documented_input_layout(path):
    """Paths of the text tree in the 'Input layout' section, relative to <input root>."""
    text = path.read_text(encoding="utf-8")
    block = re.search(r"```text\n<input root>/\n(.*?)```", text, re.S)[1]
    stack, paths = [], set()
    for line in block.splitlines():
        name = line.strip().split()[0]
        depth = (len(line) - len(line.lstrip())) // 2 - 1
        stack[depth:] = [name.rstrip("/")]
        paths.add("/".join(stack))
    return paths


@pytest.mark.parametrize("language", ["en", "zh"])
def test_input_layout_is_the_documented_input_layout(language):
    documented = documented_input_layout(ROOT / "docs" / language / "data-layout.md")
    for key, relative in env.INPUT_LAYOUT.items():
        assert relative in documented, key


def test_input_root_precedence(project, monkeypatch):
    parser = ArgumentParser()
    parser.add_argument("--fluid-root", type=Path)
    env.add_input_root_option(parser)
    env.add_data_root_option(parser)

    def fill(argv):
        args = parser.parse_args(argv)
        env.fill_input_defaults(parser, args, {"fluid_root": "grmhd"})
        return args.fluid_root
    assert fill(["--fluid-root", "x", "--input-root", "/i"]) == Path("x")
    assert fill(["--input-root", "/i"]) == Path("/i/grmhd")
    (project / ".env").write_text("EHT_INPUT_ROOT=inputs\n")
    assert fill([]) == project / "inputs/grmhd"
    monkeypatch.setenv("EHT_INPUT_ROOT", "/exported")
    assert fill([]) == Path("/exported/grmhd")
    monkeypatch.delenv("EHT_INPUT_ROOT")
    with pytest.raises(SystemExit):
        fill(["--no-env-file"])


def parse(argv, defaults, environment=None):
    parser = ArgumentParser()
    parser.add_argument("--limits", type=Path)
    env.add_data_root_option(parser)
    args = parser.parse_args(argv)
    env.fill_defaults(parser, args, defaults, environment)
    return args


def test_precedence_cli_root_then_data_root_from_cli_env_dotenv(project, monkeypatch):
    defaults = {"limits": ("analysis.limits", {})}
    assert parse(["--limits", "x", "--data-root", "/d"], defaults).limits == Path("x")
    assert parse(["--data-root", "/d"], defaults).limits == Path("/d/analysis/limits")
    (project / ".env").write_text("EHT_DATA_ROOT=data\n")
    assert parse([], defaults).limits == project / "data/analysis/limits"
    monkeypatch.setenv("EHT_DATA_ROOT", "/exported")
    assert parse([], defaults).limits == Path("/exported/analysis/limits")
    assert parse(["--data-root", "/cli"], defaults).limits == Path("/cli/analysis/limits")
    assert parse(["--no-env-file"], defaults).limits == Path("/exported/analysis/limits")


def test_missing_root_is_an_error(project):
    with pytest.raises(SystemExit):
        parse([], {"limits": ("analysis.limits", {})})


def capture(monkeypatch, module, name="run"):
    calls = []
    monkeypatch.setattr(module, name, lambda *a, **k: calls.append((a, k)))
    return calls


def test_analysis_scripts_default_roots(project, monkeypatch):
    from scripts.analysis import expected_signal, limits, summarize

    calls = capture(monkeypatch, expected_signal)
    expected_signal.main(["--data-root", "/d"])
    assert calls[0][0][1:] == (Path("/d/raptor/images"), Path("/d/analysis/expected-signal"))

    calls = capture(monkeypatch, limits)
    limits.main(["--data-root", "/d", "--output-root", "/new"])
    args = calls[0][0][0]
    assert (args.ensemble, args.expected_signal, args.output_root) == (
        Path("/d/eht/ensemble"), Path("/d/analysis/expected-signal"), Path("/new"))

    calls = capture(monkeypatch, summarize)
    summarize.main(["--data-root", "/d"])
    args = calls[0][0][0]
    assert (args.limits, args.output_root) == (Path("/d/analysis/limits"), Path("/d/analysis/summary"))


def test_eht_scripts_default_roots(project, monkeypatch):
    from scripts.eht import reconstruct, refine, sample

    args = reconstruct.parse(["--stage", "fiducial", "--day", "96", "--data-root", "/d",
                              "--input-root", "/i"])
    assert (args.i_root, args.output_root) == (Path("/d/eht/intensity/day-96"),
                                               Path("/d/eht/fiducial/day-96"))
    assert (args.uvfits_dir, args.official_script) == (
        Path("/i/eht"), Path("/i/eht/eht-imaging_pipeline.py"))
    args = reconstruct.parse(["--stage", "I", "--day", "95", "--data-root", "/d", "--input-root", "/i",
                              "--official-script", "official.py"])
    assert (args.i_root, args.output_root) == (None, Path("/d/eht/intensity/day-95"))
    assert reconstruct.locations(args)["official_script"] == Path("official.py").resolve()
    with pytest.raises(SystemExit):  # no input root and no explicit locations
        reconstruct.parse(["--stage", "I", "--day", "95", "--data-root", "/d", "--no-env-file"])
    args = refine.parse(["--day", "100", "--data-root", "/d"])
    assert (args.source_root, args.i_root, args.output_root) == (
        Path("/d/eht/fiducial/day-100"), Path("/d/eht/intensity/day-100"), Path("/d/eht/refine/day-100"))

    seen = {}
    monkeypatch.setattr(sample, "controller",
                        lambda args, out: seen.update(i=args.i_root, d=args.d_root, out=out))
    monkeypatch.setattr(sample.fcntl, "flock", lambda *a: None)
    monkeypatch.setattr(sample.os, "open", lambda *a: 0)
    monkeypatch.setattr(sample.os, "close", lambda *a: None)
    monkeypatch.setattr(Path, "mkdir", lambda *a, **k: None)
    sample.main(["--data-root", "/d"])
    assert seen["i"] == [Path(f"/d/eht/intensity/day-{day}") for day in (95, 96, 100, 101)]
    assert seen["d"] == [Path(f"/d/eht/refine/day-{day}") for day in (95, 96, 100, 101)]
    assert seen["out"] == Path("/d/eht/sample")


def test_eht_sample_and_summary_default_roots(project, monkeypatch, tmp_path):
    from scripts.eht import sample, summarize

    data = tmp_path / "d"
    with pytest.raises(SystemExit):  # stops at the worker check, after the roots are filled in
        sample.main(["--data-root", str(data), "--workers", "0"])
    (data / "eht" / "sample").mkdir(parents=True)
    (data / "eht" / "ensemble").mkdir(parents=True)
    with pytest.raises(SystemExit):  # the default output root exists: refused
        summarize.main(["--data-root", str(data)])


def test_raptor_scripts_default_roots(project, monkeypatch, tmp_path):
    from eht_axion.raptor import options as raptor_options

    core = project / "core"
    core.mkdir()
    monkeypatch.setattr(raptor_options, "ROOT", ROOT)
    args = raptor_options.options("build", ["--source-root", str(core), "--data-root", str(tmp_path / "d")])
    assert args.run_root == (tmp_path / "d" / "raptor").resolve()
    (project / ".env").write_text(f"EHT_DATA_ROOT={tmp_path / 'dotenv'}\n")
    args = raptor_options.options("build", ["--raptor-source-root", str(core)])
    assert args.run_root == (tmp_path / "dotenv" / "raptor").resolve()
    with pytest.raises(SystemExit):
        raptor_options.options("build", ["--source-root", str(core), "--run-root", str(tmp_path)])
