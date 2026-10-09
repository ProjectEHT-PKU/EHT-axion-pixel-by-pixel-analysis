"""Direct numerical tools must validate products and avoid external writes."""

import json
import numpy as np
import pytest

from eht_axion.pipeline.common import stage_directory
from eht_axion.pipeline.figures import plot_results
from eht_axion.pipeline.diagnostics import plot_diagnostics
from eht_axion.pipeline.convergence import check_convergence
from compare_original_limits import compare


MODES = ("std", "cor", "cor_grmhd")


def make_curves(path, *, nan=False):
    for dim in ("1d", "2d"):
        cuts = np.ones((3, 2))
        if nan:
            cuts[0, 0] = np.nan
        np.savez(path / f"limits-{dim}.npz", alphas=[.1, .3], cuts=cuts,
                 modes=MODES, prior="log-uniform", variant="test", observation="null")


def test_plot_refuses_stale_canonical_stage(tmp_path):
    source = tmp_path / "source"
    with stage_directory(source, "test", {}) as t:
        make_curves(t)
    np.savez(source / "limits-1d.npz", alphas=[.1, .3], cuts=np.ones((3, 2))*2,
             modes=MODES, prior="log-uniform", variant="test", observation="null")
    out = tmp_path / "figure"
    with pytest.raises(ValueError, match="checksum"):
        plot_results(source, out)
    assert not out.exists()


def test_independent_canonical_stage_can_be_plotted(tmp_path):
    source = tmp_path / "standalone-curves"
    with stage_directory(source, "curves:validated-bundle", {}) as t:
        make_curves(t)
    out = tmp_path / "figures"
    plot_results(source, out)
    assert (out / "limits.png").is_file()


def test_undeclared_extra_curves_cannot_shadow_verified_stage(tmp_path):
    source = tmp_path / "source"
    with stage_directory(source, "curves:validated-bundle", {}) as t:
        (t / "unrelated.txt").write_text("the only declared output")
    make_curves(source)
    out = tmp_path / "figures"
    with pytest.raises(ValueError, match="not a declared stage output"):
        plot_results(source, out)
    assert not out.exists()


def test_compare_refuses_nonfinite_canonical_curve(tmp_path):
    source = tmp_path / "source"
    with stage_directory(source, "test", {}) as t:
        make_curves(t, nan=True)
    with pytest.raises(ValueError, match="Invalid curves"):
        compare(tmp_path / "old", source, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_compare_refuses_mixed_legacy_prior(tmp_path):
    source, old = tmp_path / "source", tmp_path / "old"
    source.mkdir()
    old.mkdir()
    for dim in ("1d", "2d"):
        for mode in MODES:
            np.savez(source / f"{mode}_{dim}.npz", alphas=[.1, .3], cuts=[0., 1.],
                     prior="uniform" if mode == "cor" else "log-uniform",
                     observation="null", variant="legacy")
            np.savez(old / f"{mode}{'_1d' if dim == '1d' else ''}.npz",
                     alphas=[.1, .3], cuts=[0., 1.])
    with pytest.raises(ValueError, match="Incompatible legacy"):
        compare(old, source, tmp_path / "out")


def test_diagnostics_refuses_stale_stage_before_output(tmp_path):
    stage = tmp_path / "run/eht"
    with stage_directory(stage, "test", {}) as t:
        (t / "value.txt").write_text("valid")
    (stage / "value.txt").write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        plot_diagnostics(tmp_path / "run", tmp_path / "diagnostics")
    assert not (tmp_path / "diagnostics").exists()


@pytest.mark.parametrize("tool", ["plot", "diagnostics", "compare", "convergence"])
def test_auxiliary_tools_reject_output_link_first(tmp_path, tool):
    alias = tmp_path / "alias"
    target = tmp_path / "external"
    alias.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="shared directory symlink"):
        if tool == "plot":
            plot_results(tmp_path / "input", alias)
        elif tool == "diagnostics":
            plot_diagnostics(tmp_path / "run", alias)
        elif tool == "compare":
            compare(tmp_path / "old", tmp_path / "new", alias)
        else:
            check_convergence(tmp_path / "run", alias / "convergence.json")
    assert not target.exists()


def test_direct_legacy_plot_labels_limited_validation(tmp_path):
    source = tmp_path / "legacy"
    source.mkdir()
    for dim in ("1d", "2d"):
        for mode in MODES:
            np.savez(source / f"{mode}_{dim}.npz", alphas=[.1, .3], cuts=[0., 1.],
                     prior="log-uniform", observation="null", variant="legacy")
    out = tmp_path / "figures"
    plot_results(source, out)
    assert "lineage unavailable" in json.loads((out / "figure.json").read_text())["source_validation"]
    with pytest.raises(FileExistsError, match="overwrite"):
        plot_results(source, out)
