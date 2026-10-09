"""Migration regressions use only tiny synthetic legacy trees under tmp_path."""

import json
from pathlib import Path

import numpy as np
import pytest

from eht_axion.pipeline import migrate
from eht_axion.pipeline.common import (
    describe_inputs,
    manifest_path,
    sha256,
    signature,
    stage_directory,
    verify_stage,
)
from eht_axion.pipeline.layout import RunLayout, date_label

DIMS = ("1d", "2d")
MODES = ("std", "cor", "cor_grmhd")


def write_json(path, data):
    path.write_text(json.dumps(data) + "\n")


def read_json(path):
    return json.loads(path.read_text())


def make_run(
    root,
    labels=None,
    *,
    bad_identity=False,
    bad_dimensions=False,
    missing_figures=False,
    alphas=(0.1,),
):
    source = root / "legacy"
    source.mkdir()
    raw = root / "raw.txt"
    raw.write_text("immutable synthetic raw")

    def create(relative, stage, inputs, writer):
        directory = source / relative
        with stage_directory(directory, stage, {}, inputs) as temporary:
            writer(temporary)
        manifest_path(directory).rename(directory / "stage.json")
        return directory

    def processed(t):
        np.save(t / "stokes.npy", np.ones((2, 17, 2, 2, 4)))
        write_json(t / "preprocess.json", {})

    def response(t):
        np.save(t / "fits_1d.npy", np.arange(24.0).reshape(2, 4, 3))
        np.save(
            t / "fits_2d.npy",
            np.arange(12.0 if bad_dimensions else 24.0).reshape(-1, 2, 2, 3),
        )
        write_json(t / "fit.json", labels or {})

    responses = []
    for alpha in alphas:
        p = create(f"01_raptor/alpha-{alpha:g}", "raptor:preprocess", [raw], processed)
        responses.append(
            create(
                f"02_response/alpha-{alpha:g}",
                "raptor:fit",
                migrate.stage_files(p),
                response,
            )
        )
    for dim in DIMS:

        def templates(t, dim=dim):
            np.savez(t / "templates.npz", amplitude_mean=np.ones((len(alphas), 4)))
            write_json(t / "summarize.json", {"dimension": dim})

        create(
            f"03_templates/{dim}",
            "raptor:summarize",
            [f / f"fits_{dim}.npy" for f in responses],
            templates,
        )

    def eht(t):
        np.save(t / "mask.npy", np.ones((2, 2), dtype=bool))
        np.savez(t / "selection.npz", indices=[0, 1])
        write_json(t / "report.json", {})
        for dim, shape in (("1d", (4, 4)), ("2d", (4, 2, 2))):
            fields = {
                key: np.arange(np.prod(shape)).reshape(shape)
                for key in (
                    "mean",
                    "std",
                    "center",
                    "unwrapped_mean",
                    "concentration",
                    "circular_std",
                    "count",
                )
            }
            if dim == "1d":
                fields["covariance"] = np.ones((4, 4, 4))
            else:
                fields["intensity_mean"] = np.ones(shape)
            np.savez(t / f"summary_{dim}.npz", **fields)
        for day in (95, 96, 100, 101):
            np.save(t / f"stokes_{day}.npy", np.ones((2, 2, 2, 4)))
            np.save(t / f"evpa_1d_{day}.npy", np.ones((2, 4)))

    e = create("04_eht", "eht:ensemble", [raw], eht)
    bp = create(
        "01_raptor_background", "raptor:background-preprocess", [raw], processed
    )

    def background(t):
        write_json(t / "background.json", {"sigma": 0.1})
        np.savez(t / "background_diagnostics.npz", values=[1, 2])

    b = create(
        "05_background",
        "raptor:background",
        migrate.stage_files(bp) + [e / "mask.npy"],
        background,
    )
    for dim in DIMS:

        def observations(t):
            np.savez(t / "observations.npz", observation=np.zeros(4))
            write_json(t / "observations.json", {})

        create(
            f"06_observations/{dim}",
            "observations:covariance",
            migrate.stage_files(e) + [source / "03_templates" / dim / "templates.npz"],
            observations,
        )
    inputs = [b / "stage.json", b / "background.json"]
    for dim in DIMS:
        inputs += [
            source / "06_observations" / dim / "stage.json",
            source / "03_templates" / dim / "stage.json",
            source / "06_observations" / dim / "observations.npz",
            source / "03_templates" / dim / "templates.npz",
        ]

    def inference(t):
        for dim in DIMS:
            for mode in MODES:
                np.savez(
                    t / f"{mode}_{dim}.npz",
                    alphas=alphas,
                    cuts=np.full(len(alphas), 0.2),
                    variant="legacy",
                    observation="measured"
                    if bad_identity and mode == "cor" and dim == "2d"
                    else "null",
                    prior="log-uniform",
                )
                write_json(t / f"{mode}_{dim}.json", {"source": "synthetic"})
        write_json(
            t / "run.json",
            {"options": {}, "results": {"std_1d": {"alphas": list(alphas)}}},
        )

    i = create("07_inference", "inference:six-curves", inputs, inference)
    if not missing_figures:
        create(
            "08_figures",
            "figures:limits",
            migrate.stage_files(i),
            lambda t: (t / "limits_by_dimension.png").write_bytes(b"synthetic figure"),
        )
    return source, raw


def snapshot(source):
    return {
        str(p.relative_to(source)): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in source.rglob("*")
        if p.is_file()
    }


def assert_clean(destination):
    assert not destination.exists()
    assert not destination.is_symlink()
    assert not destination.with_name(destination.name + ".lock").exists()
    assert not list(destination.parent.glob(destination.name + ".building-*"))


def test_migrate_run_complete_chain_and_rebound_lineage(tmp_path):
    source, _raw = make_run(tmp_path)
    before = snapshot(source)
    out = tmp_path / "canonical"
    assert migrate.migrate_run(source, out) == out
    assert snapshot(source) == before
    record = read_json(out / "migration.json")
    assert record["scientific_recomputation"] is False
    assert record["implementation"]["kind"] == "migration"
    assert record["implementation"]["sha256"] == sha256(Path(migrate.__file__))
    assert len(record["completed_stages"]) == 11
    layout = RunLayout(out)
    for relative in record["completed_stages"]:
        directory = out / relative
        manifest = verify_stage(directory, recursive=True)
        assert manifest["scientific_recomputation"] is False
        assert manifest["imported_from"]["implementation"]["kind"] == "migration"
        assert "origin_inputs" in manifest
        assert all(
            Path(e["path"]).is_file()
            for e in manifest["inputs"] + manifest["origin_inputs"]
        )
        assert ".building-" not in manifest_path(directory).read_text()
        assert (
            manifest["signature"]
            == signature(manifest["stage"], manifest["config"], manifest["inputs"])[1]
        )
        for field in ("inputs", "origin_inputs"):
            assert manifest[field] == describe_inputs(
                [e["path"] for e in manifest[field]]
            )
        with pytest.raises(ValueError, match="imported"):
            verify_stage(directory, recursive=True, allow_imported=False)
    for dim in DIMS:
        template = verify_stage(layout.templates(dim))
        assert [e["path"] for e in template["inputs"]] == [
            str(manifest_path(layout.response(0.1))),
            str(layout.response(0.1) / f"response-{dim}.npz"),
        ]
        observed = verify_stage(layout.observations(dim))
        expected = migrate.stage_files(layout.eht) + [
            manifest_path(layout.templates(dim)),
            layout.templates(dim) / "templates.npz",
        ]
        assert observed["inputs"] == describe_inputs(expected)
    assert verify_stage(layout.background)["inputs"] == describe_inputs(
        migrate.stage_files(layout.background_processed)
        + [manifest_path(layout.eht), layout.eht / "mask.npy"]
    )
    for dim in DIMS:
        original = np.load(
            source / "02_response/alpha-0.1" / f"fits_{dim}.npy", allow_pickle=False
        )
        with np.load(
            layout.response(0.1) / f"response-{dim}.npz", allow_pickle=False
        ) as z:
            for column, key in enumerate(("amplitude", "phase", "baseline")):
                np.testing.assert_array_equal(z[key], original[..., column])
        with np.load(layout.inference / f"limits-{dim}.npz", allow_pickle=False) as z:
            np.testing.assert_array_equal(z["cuts"], np.full((3, 1), 0.2))
    assert (layout.processed(0.1) / "stokes.npy").is_symlink()
    assert (
        layout.figures / "limits-by-dimension.png"
    ).read_bytes() == b"synthetic figure"
    with np.load(layout.eht / date_label(96) / "summary.npz", allow_pickle=False) as z:
        np.testing.assert_array_equal(z["evpa_mean_1d"], np.arange(4, 8))
    with pytest.raises(FileExistsError):
        migrate.migrate_run(source, out)
    assert snapshot(source) == before


@pytest.mark.parametrize(
    "metadata, expected, label_source, verified",
    [
        ({}, [1600, 1610], "default_convention", False),
        (
            {
                "snapshot_labels": [42, 47],
                "snapshot_labels_source": "producer",
                "snapshot_labels_verified": True,
            },
            [42, 47],
            "original_fit_metadata",
            True,
        ),
        (
            {
                "snapshot_labels": [42, 47],
                "snapshot_labels_source": "caller",
                "snapshot_labels_verified": False,
            },
            [42, 47],
            "original_fit_metadata_unverified",
            False,
        ),
    ],
)
def test_snapshot_labels_provenance(
    tmp_path, metadata, expected, label_source, verified
):
    source, _ = make_run(tmp_path, metadata)
    out = migrate.migrate_run(source, tmp_path / "new")
    response = RunLayout(out).response(0.1)
    report = read_json(response / "response.json")
    assert report["snapshot_labels"] == expected
    assert report["snapshot_labels_source"] == label_source
    assert report["snapshot_labels_verified"] is verified
    for dim in DIMS:
        with np.load(response / f"response-{dim}.npz", allow_pickle=False) as z:
            assert z["snapshot_labels"].tolist() == expected


@pytest.mark.parametrize(
    "kwargs",
    [
        {"bad_identity": True},
        {"bad_dimensions": True},
        {"missing_figures": True},
        {
            "labels": {
                "snapshot_labels": [20, 20],
                "snapshot_labels_source": "producer",
                "snapshot_labels_verified": True,
            }
        },
        {
            "labels": {
                "snapshot_labels": [30, 20],
                "snapshot_labels_source": "producer",
                "snapshot_labels_verified": True,
            }
        },
        {
            "labels": {
                "snapshot_labels": [10],
                "snapshot_labels_source": "producer",
                "snapshot_labels_verified": True,
            }
        },
    ],
)
def test_preflight_rejects_invalid_source_without_building(
    tmp_path, monkeypatch, kwargs
):
    source, _ = make_run(tmp_path, **kwargs)
    before = snapshot(source)
    out = tmp_path / "new"

    def forbidden(*args, **kwargs):
        pytest.fail("preflight must finish before creating a temporary run")

    monkeypatch.setattr(migrate.tempfile, "mkdtemp", forbidden)
    with pytest.raises((ValueError, FileNotFoundError)):
        migrate.migrate_run(source, out)
    assert_clean(out)
    assert snapshot(source) == before


def test_source_stage_identity_and_foreign_lineage_rejected(tmp_path):
    source, _ = make_run(tmp_path)
    p = source / "08_figures/stage.json"
    manifest = read_json(p)
    manifest["stage"] = "not:figures"
    write_json(p, manifest)
    with pytest.raises(ValueError, match="source stage identity"):
        migrate.migrate_run(source, tmp_path / "new")
    manifest["stage"] = "figures:limits"
    foreign = tmp_path / "foreign.txt"
    foreign.write_text("not this run")
    manifest["inputs"] = describe_inputs([foreign])
    write_json(p, manifest)
    with pytest.raises(ValueError, match="different run"):
        migrate.migrate_run(source, tmp_path / "new")
    assert_clean(tmp_path / "new")


@pytest.mark.parametrize("failure", ["copy", "rename", "rebind"])
def test_late_failure_cleans_run_and_allows_retry(tmp_path, monkeypatch, failure):
    source, _ = make_run(tmp_path)
    before = snapshot(source)
    out = tmp_path / "new"
    with monkeypatch.context() as patch:
        if failure == "copy":
            original = migrate.shutil.copy2

            def copy(src, dst, *args, **kwargs):
                if Path(src).name == "limits_by_dimension.png":
                    assert not out.exists()
                    raise OSError("injected disk full")
                return original(src, dst, *args, **kwargs)

            patch.setattr(migrate.shutil, "copy2", copy)
        elif failure == "rename":
            original = Path.rename

            def rename(self, target):
                if Path(target) == out:
                    assert (self / "migration.json").is_file()
                    raise OSError("injected rename failure")
                return original(self, target)

            patch.setattr(Path, "rename", rename)
        else:

            def rebind(*args):
                raise OSError("injected rebind failure")

            patch.setattr(migrate, "_rebind", rebind)
        with pytest.raises(OSError, match="injected"):
            migrate.migrate_run(source, out)
    assert_clean(out)
    assert snapshot(source) == before
    migrate.migrate_run(source, out)
    verify_stage(out / "figures", recursive=True)
    assert snapshot(source) == before


@pytest.mark.parametrize("mutation", ["raw", "manifest", "response"])
def test_original_lineage_remains_verifiable_after_migration(tmp_path, mutation):
    source, raw = make_run(tmp_path)
    out = migrate.migrate_run(source, tmp_path / "new")
    path = {
        "raw": raw,
        "manifest": source / "02_response/alpha-0.1/stage.json",
        "response": source / "02_response/alpha-0.1/fits_1d.npy",
    }[mutation]
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="lineage|checksum"):
        verify_stage(out / "figures", recursive=True)


def test_source_mutation_during_construction_aborts(tmp_path, monkeypatch):
    source, raw = make_run(tmp_path)
    original = migrate.shutil.copy2

    def copy(src, dst, *args, **kwargs):
        result = original(src, dst, *args, **kwargs)
        if Path(src).name == "limits_by_dimension.png":
            raw.write_text("external concurrent edit")
        return result

    monkeypatch.setattr(migrate.shutil, "copy2", copy)
    with pytest.raises(ValueError, match="lineage|changed"):
        migrate.migrate_run(source, tmp_path / "new")
    assert_clean(tmp_path / "new")


def make_inputs(root):
    root.mkdir()
    for alpha in (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44):
        p = root / "raptor-stokes/res-100" / f"alpha-{alpha:g}.npy"
        p.parent.mkdir(parents=True, exist_ok=True)
        np.save(p, np.ones((2, 17, 2, 2, 4)))
    p = root / "raptor-stokes/res-400/alpha-0.3/stokes-t1600-t1610.npy"
    p.parent.mkdir(parents=True)
    np.save(p, np.ones((2, 1, 2, 2, 4)))
    for day in (95, 96, 100, 101):
        p = root / "eht-image-ensemble" / f"230_{day}.pkl"
        p.parent.mkdir(exist_ok=True)
        p.write_bytes(b"synthetic pickle placeholder; never unpickled")
    p = root / "raptor-transport-diagnostics/snapshot-2110-res-400/common.npy"
    p.parent.mkdir(parents=True)
    np.save(p, np.ones((2, 2)))
    files = describe_inputs(sorted(p for p in root.rglob("*") if p.is_file()))
    write_json(
        root / "manifest.json",
        {
            "files": [
                {
                    "path": str(Path(e["path"]).relative_to(root)),
                    "sha256": e["sha256"],
                    "bytes": e["size"],
                }
                for e in files
            ]
        },
    )


def test_migrate_inputs_complete_chain(tmp_path):
    source = tmp_path / "input"
    make_inputs(source)
    before = snapshot(source)
    out = tmp_path / "view"
    migrate.migrate_inputs(source, out)
    manifest = verify_stage(out, recursive=True)
    assert snapshot(source) == before
    dataset = read_json(out / "dataset.json")
    assert len(dataset["files"]) == 14
    assert len(manifest["origin_inputs"]) == 14
    for entry in dataset["files"]:
        p = out / entry["path"]
        assert p.is_symlink() and p.is_file()
        assert sha256(p) == entry["sha256"]
    assert ".building-" not in (out / "manifest.json").read_text()
    with pytest.raises(FileExistsError):
        migrate.migrate_inputs(source, out)
    (source / "raptor-stokes/res-100/alpha-0.1.npy").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="lineage|checksum"):
        verify_stage(out, recursive=True)


@pytest.mark.parametrize("failure", ["checksum", "copy", "rename"])
def test_migrate_inputs_failure_and_retry(tmp_path, monkeypatch, failure):
    source = tmp_path / "input"
    make_inputs(source)
    out = tmp_path / "view"
    before = snapshot(source)
    with monkeypatch.context() as patch:
        if failure == "checksum":
            manifest = read_json(source / "manifest.json")
            manifest["files"][-1]["sha256"] = "0" * 64
            write_json(source / "manifest.json", manifest)
        elif failure == "copy":
            original = migrate.link

            def link(src, dst):
                if Path(src).suffix == ".pkl":
                    raise OSError("injected link failure")
                return original(src, dst)

            patch.setattr(migrate, "link", link)
        else:
            original = Path.rename

            def rename(self, target):
                if Path(target) == out:
                    raise OSError("injected rename failure")
                return original(self, target)

            patch.setattr(Path, "rename", rename)
        with pytest.raises((ValueError, OSError)):
            migrate.migrate_inputs(source, out)
    assert_clean(out)
    if failure == "checksum":
        (source / "manifest.json").write_bytes(before["manifest.json"][0])
    else:
        assert snapshot(source) == before
    migrate.migrate_inputs(source, out)
    verify_stage(out, recursive=True)


@pytest.mark.parametrize("operation", ["run", "inputs"])
@pytest.mark.parametrize("kind", ["destination", "dangling", "ancestor", "nested"])
def test_migration_link_and_tree_safety(tmp_path, operation, kind):
    if operation == "run":
        source, _ = make_run(tmp_path)
        migrate_fn = migrate.migrate_run
    else:
        source = tmp_path / "inputs"
        make_inputs(source)
        migrate_fn = migrate.migrate_inputs
    before = snapshot(source)
    out = tmp_path / "new"
    outside = tmp_path / "outside"
    outside.mkdir()
    if kind == "destination":
        out.symlink_to(outside, target_is_directory=True)
    elif kind == "dangling":
        out.symlink_to(outside / "missing", target_is_directory=True)
    elif kind == "ancestor":
        alias = tmp_path / "alias"
        alias.symlink_to(outside, target_is_directory=True)
        out = alias / "new"
    else:
        out = source / "nested"
    with pytest.raises((ValueError, FileExistsError)):
        migrate_fn(source, out)
    assert snapshot(source) == before
    assert not list(outside.iterdir())


def test_run_rejects_symlinked_stage_directory(tmp_path):
    source, _ = make_run(tmp_path)
    original = source / "08_figures"
    moved = tmp_path / "external-figures"
    original.rename(moved)
    original.symlink_to(moved, target_is_directory=True)
    before = snapshot(moved)
    with pytest.raises(ValueError, match="symlink"):
        migrate.migrate_run(source, tmp_path / "new")
    assert_clean(tmp_path / "new")
    assert snapshot(moved) == before


@pytest.mark.parametrize("operation", ["run", "inputs"])
def test_read_only_source_alias_is_allowed(tmp_path, operation):
    if operation == "run":
        source, _ = make_run(tmp_path)
        migrate_fn = migrate.migrate_run
        target = "figures"
    else:
        source = tmp_path / "inputs"
        make_inputs(source)
        migrate_fn = migrate.migrate_inputs
        target = "."
    alias = tmp_path / "source-alias"
    alias.symlink_to(source, target_is_directory=True)
    before = snapshot(source)
    out = migrate_fn(alias, tmp_path / "new")
    verify_stage(out / target, recursive=True)
    assert snapshot(source) == before


def test_imported_cache_reuse_preserves_origin_checks(tmp_path):
    source, _ = make_run(tmp_path)
    out = migrate.migrate_run(source, tmp_path / "new")
    response = RunLayout(out).response(0.1)
    manifest = verify_stage(response, recursive=True)
    paths = [e["path"] for e in manifest["inputs"]]
    before = snapshot(response)
    with stage_directory(
        response, manifest["stage"], manifest["config"], paths, reuse_imported=True
    ) as temporary:
        assert temporary is None
    assert snapshot(response) == before
    original = source / "02_response/alpha-0.1/fits_1d.npy"
    original.write_bytes(original.read_bytes() + b"changed")
    with (
        pytest.raises(ValueError, match="lineage|checksum"),
        stage_directory(
            response, manifest["stage"], manifest["config"], paths, reuse_imported=True
        ),
    ):
        pytest.fail("changed origin must not be reused")
    assert snapshot(response) == before


def test_multi_alpha_run_rebinding_and_imported_ancestry(tmp_path):
    source, _ = make_run(tmp_path, alphas=(0.1, 0.3))
    out = migrate.migrate_run(source, tmp_path / "new")
    layout = RunLayout(out)
    for dim in DIMS:
        template = verify_stage(layout.templates(dim), recursive=True)
        expected = [
            p
            for a in (0.1, 0.3)
            for p in (
                manifest_path(layout.response(a)),
                layout.response(a) / f"response-{dim}.npz",
            )
        ]
        assert template["inputs"] == describe_inputs(expected)
        assert len(template["imported_ancestors"]) == 2
        for ancestor in template["imported_ancestors"]:
            assert ancestor["sha256"] == sha256(Path(ancestor["manifest"]))
            assert str(out) in ancestor["manifest"]
            for upstream in ancestor["ancestors"]:
                assert upstream["sha256"] == sha256(Path(upstream["manifest"]))
    verify_stage(layout.figures, recursive=True)


def _legacy_from_canonical(canonical, source, alphas):
    """Write real tiny-chain products in the exact original legacy contract."""
    layout = RunLayout(canonical)

    def create(relative, original, inputs, writer):
        manifest = verify_stage(original)
        directory = source / relative
        with stage_directory(
            directory, manifest["stage"], manifest["config"], inputs
        ) as temporary:
            writer(temporary)
        manifest_path(directory).rename(directory / "stage.json")
        return directory

    def copy_processed(original, temporary):
        (temporary / "stokes.npy").write_bytes((original / "stokes.npy").read_bytes())
        (temporary / "preprocess.json").write_bytes(
            (original / "processing.json").read_bytes()
        )

    for alpha in alphas:
        original = layout.processed(alpha)
        inputs = [e["path"] for e in verify_stage(original)["inputs"]]
        p = create(
            f"01_raptor/alpha-{alpha:g}",
            original,
            inputs,
            lambda t, original=original: copy_processed(original, t),
        )
        original = layout.response(alpha)

        def response(t, original=original):
            for dim in DIMS:
                with np.load(original / f"response-{dim}.npz", allow_pickle=False) as z:
                    np.save(
                        t / f"fits_{dim}.npy",
                        np.stack(
                            [z[k] for k in ("amplitude", "phase", "baseline")], axis=-1
                        ),
                    )
            (t / "fit.json").write_bytes((original / "response.json").read_bytes())

        create(
            f"02_response/alpha-{alpha:g}", original, migrate.stage_files(p), response
        )
    for dim in DIMS:
        original = layout.templates(dim)

        def template(t, original=original):
            (t / "templates.npz").write_bytes((original / "templates.npz").read_bytes())
            (t / "summarize.json").write_bytes((original / "summary.json").read_bytes())

        create(
            f"03_templates/{dim}",
            original,
            [
                source / "02_response" / f"alpha-{a:g}" / f"fits_{dim}.npy"
                for a in alphas
            ],
            template,
        )

    def eht(t):
        for name in ("mask.npy", "selection.npz", "report.json"):
            (t / name).write_bytes((layout.eht / name).read_bytes())
        summaries = {dim: {} for dim in DIMS}
        for day in (95, 96, 100, 101):
            original = layout.eht / date_label(day)
            for new, old in (
                ("stokes.npy", f"stokes_{day}.npy"),
                ("evpa-1d.npy", f"evpa_1d_{day}.npy"),
            ):
                (t / old).write_bytes((original / new).read_bytes())
            with np.load(original / "summary.npz", allow_pickle=False) as z:
                for dim in DIMS:
                    for key in (
                        "mean",
                        "std",
                        "center",
                        "unwrapped_mean",
                        "concentration",
                        "circular_std",
                        "count",
                    ):
                        summaries[dim].setdefault(key, []).append(
                            z[f"evpa_{key}_{dim}"]
                        )
                summaries["1d"].setdefault("covariance", []).append(z["covariance_1d"])
                summaries["2d"].setdefault("intensity_mean", []).append(
                    z["intensity_mean"]
                )
        for dim in DIMS:
            np.savez(
                t / f"summary_{dim}.npz",
                **{k: np.stack(v) for k, v in summaries[dim].items()},
            )

    e = create(
        "04_eht",
        layout.eht,
        [e["path"] for e in verify_stage(layout.eht)["inputs"]],
        eht,
    )
    original = layout.background_processed
    bp = create(
        "01_raptor_background",
        original,
        [e["path"] for e in verify_stage(original)["inputs"]],
        lambda t: copy_processed(original, t),
    )

    def background(t):
        (t / "background.json").write_bytes(
            (layout.background / "summary.json").read_bytes()
        )
        (t / "background_diagnostics.npz").write_bytes(
            (layout.background / "evpa-difference.npz").read_bytes()
        )

    b = create(
        "05_background",
        layout.background,
        migrate.stage_files(bp) + [e / "mask.npy"],
        background,
    )
    for dim in DIMS:
        original = layout.observations(dim)

        def observations(t, original=original):
            for name in ("observations.npz", "observations.json"):
                (t / name).write_bytes((original / name).read_bytes())

        create(
            f"06_observations/{dim}",
            original,
            migrate.stage_files(e) + [source / "03_templates" / dim / "templates.npz"],
            observations,
        )
    inputs = [b / "stage.json", b / "background.json"]
    for dim in DIMS:
        inputs += [
            source / "06_observations" / dim / "stage.json",
            source / "03_templates" / dim / "stage.json",
            source / "06_observations" / dim / "observations.npz",
            source / "03_templates" / dim / "templates.npz",
        ]

    def inference(t):
        for dim in DIMS:
            with np.load(
                layout.inference / f"limits-{dim}.npz", allow_pickle=False
            ) as z:
                for k, mode in enumerate(MODES):
                    np.savez(
                        t / f"{mode}_{dim}.npz",
                        alphas=z["alphas"],
                        cuts=z["cuts"][k],
                        **{key: z[key] for key in ("variant", "observation", "prior")},
                    )
                    (t / f"{mode}_{dim}.json").write_bytes(
                        (layout.inference / f"posterior-{mode}-{dim}.json").read_bytes()
                    )
        (t / "run.json").write_bytes((layout.inference / "run.json").read_bytes())

    i = create("07_inference", layout.inference, inputs, inference)
    create(
        "08_figures",
        layout.figures,
        migrate.stage_files(i),
        lambda t: (t / "limits_by_dimension.png").write_bytes(b"synthetic"),
    )


@pytest.mark.parametrize("operation", ["run", "inputs"])
def test_migration_lock_conflict_preserves_owner(tmp_path, operation):
    if operation == "run":
        source, _ = make_run(tmp_path)
        migrate_fn = migrate.migrate_run
    else:
        source = tmp_path / "inputs"
        make_inputs(source)
        migrate_fn = migrate.migrate_inputs
    out = tmp_path / "new"
    lock = out.with_name(out.name + ".lock")
    lock.write_text("another owner")
    before = snapshot(source)
    with pytest.raises(RuntimeError, match="Migration locked"):
        migrate_fn(source, out)
    assert lock.read_text() == "another owner"
    assert snapshot(source) == before
    assert not out.exists()
    assert not list(tmp_path.glob("new.building-*"))


def test_lock_record_write_failure_is_cleaned(tmp_path, monkeypatch):
    source, _ = make_run(tmp_path)
    out = tmp_path / "new"

    def fail(*args, **kwargs):
        raise OSError("injected lock disk full")

    monkeypatch.setattr(migrate.json, "dump", fail)
    with pytest.raises(OSError, match="lock disk full"):
        migrate.migrate_run(source, out)
    assert_clean(out)


@pytest.mark.parametrize("mutation", ["raw", "manifest"])
def test_inputs_source_mutation_during_construction_aborts(
    tmp_path, monkeypatch, mutation
):
    source = tmp_path / "inputs"
    make_inputs(source)
    out = tmp_path / "new"
    original = migrate.link

    def link(src, dst):
        result = original(src, dst)
        if Path(src).name == "common.npy":
            target = source / (
                "manifest.json"
                if mutation == "manifest"
                else "raptor-stokes/res-100/alpha-0.1.npy"
            )
            if mutation == "manifest":
                record = read_json(target)
                record["external_edit"] = True
                write_json(target, record)
            else:
                target.write_bytes(target.read_bytes() + b"external edit")
        return result

    monkeypatch.setattr(migrate, "link", link)
    with pytest.raises(ValueError, match="changed|lineage|checksum"):
        migrate.migrate_inputs(source, out)
    assert_clean(out)


def test_rebind_internal_origins_and_exact_path_prefix(tmp_path):
    root = tmp_path / "private"
    root.mkdir()
    destination = tmp_path / "final"
    # A source that shares a string prefix with root must not be remapped.
    external = tmp_path / "private-original"
    external.mkdir()
    raw = external / "raw.txt"
    raw.write_text("source")
    origin = {"path": str(external), "implementation": {"kind": "migration"}}
    with stage_directory(root / "a", "test", {}, [raw], imported_from=origin) as t:
        (t / "value.txt").write_text("converted")
    with stage_directory(
        root / "b",
        "test",
        {},
        [raw],
        imported_from=origin,
        origin_paths=migrate.stage_files(root / "a"),
    ) as t:
        (t / "value.txt").write_text("converted again")
    migrate._rebind(root, destination, ["a", "b"])
    root.rename(destination)
    manifest = verify_stage(destination / "b", recursive=True)
    assert manifest["inputs"] == describe_inputs([raw])
    assert manifest["origin_inputs"] == describe_inputs(
        migrate.stage_files(destination / "a")
    )
    ancestor = manifest["imported_ancestors"][0]
    assert ancestor["manifest"] == str(destination / "a/manifest.json")
    assert ancestor["sha256"] == sha256(destination / "a/manifest.json")
    assert manifest["imported_from"]["path"] == str(external)
    assert ".building-" not in (destination / "b/manifest.json").read_text()


def test_actual_inference_can_reuse_migrated_run(tmp_path_factory):
    # Reuse the parent's actual deterministic scientific tiny-chain generator;
    # no fixture or files in any other test module are modified.
    from test_pipeline_integrity import tiny_chain

    from eht_axion.pipeline.inference import run_inference
    from eht_axion.pipeline.workflow import read_config

    root, _, canonical, config_path = tiny_chain.__wrapped__(tmp_path_factory)
    source, out = root / "legacy", root / "migrated"
    config = read_config(config_path)
    _legacy_from_canonical(canonical, source, config["data"]["alphas"])
    # Genuine pre-null-only metadata: false must be normalized only in the
    # new manifest, with the original receipt and bytes preserved.
    old_inference = source / "07_inference/stage.json"
    value = json.loads(old_inference.read_text())
    value["config"]["options"]["observed"] = False
    write_json(old_inference, value)
    old_figures = source / "08_figures/stage.json"
    value = json.loads(old_figures.read_text())
    value["inputs"] = describe_inputs([x["path"] for x in value["inputs"]])
    write_json(old_figures, value)
    before = snapshot(source)
    migrate.migrate_run(source, out)
    assert "observed" not in verify_stage(out / "inference")["config"]["options"]
    assert snapshot(source) == before
    stamp = (out / "inference/manifest.json").read_bytes()
    with pytest.raises(ValueError, match="imported"):
        run_inference(out, config)
    assert run_inference(out, config, reuse_imported=True) == out / "inference"
    assert (out / "inference/manifest.json").read_bytes() == stamp
    assert snapshot(source) == before
