"""Default roots below the data root: documented layout, precedence, and every script's defaults."""
from argparse import ArgumentParser
import os
from pathlib import Path
import re

import pytest

import _data_root
import _environment

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def project(tmp_path, monkeypatch):
    """An isolated project root for dotenv files, without exported settings."""
    monkeypatch.setattr(_environment, "ROOT", tmp_path)
    for name in list(os.environ):
        if name.startswith("EHT_") or name == _data_root.DATA_ROOT_KEY:
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
    for key, relative in _data_root.LAYOUT.items():
        expected = relative.format(day="<day>").replace("runs/main", "runs/<name>")
        if key == "raptor_data":
            expected = "raptor/delivered"
        assert expected in documented, key


def parse(argv, defaults, environment=None):
    parser = ArgumentParser()
    parser.add_argument("--limits", type=Path)
    _data_root.add_data_root_option(parser)
    args = parser.parse_args(argv)
    _data_root.fill_defaults(parser, args, defaults, environment)
    return args


def test_precedence_cli_root_then_data_root_from_cli_env_dotenv(project, monkeypatch):
    defaults = {"limits": ("analysis_limits", {})}
    assert parse(["--limits", "x", "--data-root", "/d"], defaults).limits == Path("x")
    assert parse(["--data-root", "/d"], defaults).limits == Path("/d/analysis/limits")
    (project / ".env").write_text("M87_AXION_DATA_ROOT=data\n")
    assert parse([], defaults).limits == project / "data/analysis/limits"
    monkeypatch.setenv("M87_AXION_DATA_ROOT", "/exported")
    assert parse([], defaults).limits == Path("/exported/analysis/limits")
    assert parse(["--data-root", "/cli"], defaults).limits == Path("/cli/analysis/limits")
    assert parse(["--no-env-file"], defaults).limits == Path("/exported/analysis/limits")


def test_missing_root_is_an_error(project):
    with pytest.raises(SystemExit):
        parse([], {"limits": ("analysis_limits", {})})


def capture(monkeypatch, module, name="run"):
    calls = []
    monkeypatch.setattr(module, name, lambda *a, **k: calls.append((a, k)))
    return calls


def test_analysis_scripts_default_roots(project, monkeypatch):
    import analysis_check_covariance
    import analysis_check_grid
    import analysis_limits
    import analysis_summarize
    import analysis_templates

    calls = capture(monkeypatch, analysis_templates)
    analysis_templates.main(["--data-root", "/d"])
    assert calls[0][0][1:] == (Path("/d/raptor/delivered/data"), Path("/d/analysis/templates"))

    calls = capture(monkeypatch, analysis_limits)
    analysis_limits.main(["--data-root", "/d", "--output-root", "/new"])
    args = calls[0][0][0]
    assert (args.ensemble, args.templates, args.output_root) == (
        Path("/d/eht/ensemble"), Path("/d/analysis/templates"), Path("/new"))

    calls = capture(monkeypatch, analysis_check_grid)
    analysis_check_grid.main(["--data-root", "/d", "--check", "prior", "--coupling-max", "1e4"])
    args = calls[0][0][0]
    assert (args.source, args.output_root) == (Path("/d/analysis/limits"), Path("/d/analysis/check-prior"))
    with pytest.raises(SystemExit):
        analysis_check_grid.main(["--data-root", "/d"])

    monkeypatch.setattr(analysis_check_covariance.provenance, "verification_session",
                        lambda: (_ for _ in ()).throw(RuntimeError("roots filled in")))
    with pytest.raises(RuntimeError, match="roots filled in"):
        analysis_check_covariance.main(["--data-root", str(project / "d")])
    (project / "d" / "analysis" / "check-covariance").mkdir(parents=True)
    with pytest.raises(FileExistsError):
        analysis_check_covariance.main(["--data-root", str(project / "d")])

    calls = capture(monkeypatch, analysis_summarize)
    analysis_summarize.main(["--data-root", "/d"])
    args = calls[0][0][0]
    assert (args.limits, args.gridcheck, args.priorcheck, args.covcheck, args.output_root) == (
        Path("/d/analysis/limits"), Path("/d/analysis/check-grid"), Path("/d/analysis/check-prior"),
        Path("/d/analysis/check-covariance"), Path("/d/analysis/summary"))
    analysis_summarize.main(["--limits", "a", "--gridcheck", "b", "--priorcheck", "c", "--output-root", "e"])
    assert calls[1][0][0].covcheck is None


def test_eht_scripts_default_roots(project, monkeypatch):
    import eht_export_bundle
    import eht_reconstruct
    import eht_refine

    args = eht_reconstruct.parse(["--stage", "fiducial", "--day", "96", "--data-root", "/d"])
    assert (args.i_root, args.output_root) == (Path("/d/eht/intensity/day-96"),
                                               Path("/d/eht/fiducial/day-96"))
    args = eht_reconstruct.parse(["--stage", "I", "--day", "95", "--data-root", "/d"])
    assert (args.i_root, args.output_root) == (None, Path("/d/eht/intensity/day-95"))
    args = eht_refine.parse(["--day", "100", "--data-root", "/d"])
    assert (args.source_root, args.i_root, args.output_root) == (
        Path("/d/eht/fiducial/day-100"), Path("/d/eht/intensity/day-100"), Path("/d/eht/refine/day-100"))

    seen = {}
    monkeypatch.setattr(eht_export_bundle, "export",
                        lambda config, i, d, out: seen.update(i=i, d=d, out=out) or {"files": []})
    eht_export_bundle.main(["--data-root", "/d"])
    assert seen["i"] == [Path(f"/d/eht/intensity/day-{day}") for day in (95, 96, 100, 101)]
    assert seen["d"] == [Path(f"/d/eht/refine/day-{day}") for day in (95, 96, 100, 101)]
    assert seen["out"] == Path("/d/eht/bundle")


def test_eht_sample_and_summary_default_roots(project, monkeypatch, tmp_path):
    import eht_sample
    import eht_summarize

    data = tmp_path / "d"
    with pytest.raises(SystemExit):  # stops at the worker check, after the roots are filled in
        eht_sample.main(["--data-root", str(data), "--workers", "0"])
    (data / "eht" / "sample").mkdir(parents=True)
    (data / "eht" / "ensemble").mkdir(parents=True)
    with pytest.raises(SystemExit):  # the default output root exists: refused
        eht_summarize.main(["--data-root", str(data)])


def test_raptor_scripts_default_roots(project, monkeypatch, tmp_path):
    import _raptor_options
    import raptor_merge

    core = project / "core"
    core.mkdir()
    monkeypatch.setattr(_raptor_options, "ROOT", ROOT)
    args = _raptor_options.options("build", ["--source-root", str(core), "--data-root", str(tmp_path / "d")])
    assert args.run_root == (tmp_path / "d" / "raptor" / "runs" / "main").resolve()
    (project / ".env.raptor").write_text(f"EHT_RAPTOR_RUN_ROOT={tmp_path / 'explicit'}\n")
    args = _raptor_options.options("build", ["--source-root", str(core), "--data-root", str(tmp_path / "d")])
    assert args.run_root == (tmp_path / "explicit").resolve()

    data = tmp_path / "m"
    for name in ("main", "background"):
        (data / "raptor" / "runs" / name).mkdir(parents=True)
    (data / "raptor" / "plan.json").write_text('{"config": {}}')
    monkeypatch.setattr(raptor_merge, "validate_ray_tracing",
                        lambda config: (_ for _ in ()).throw(RuntimeError("stop")))
    parser_args = {}
    real_parse = ArgumentParser.parse_args

    def spy(self, *a, **k):
        result = real_parse(self, *a, **k)
        parser_args["args"] = result
        return result
    monkeypatch.setattr(ArgumentParser, "parse_args", spy)
    with pytest.raises(RuntimeError, match="stop"):
        raptor_merge.main(["--data-root", str(data)])
    args = parser_args["args"]
    assert args.job_root == [data / "raptor/runs/background", data / "raptor/runs/main"]
    assert (args.plan, args.output_root) == (data / "raptor/plan.json", data / "raptor/delivered")
