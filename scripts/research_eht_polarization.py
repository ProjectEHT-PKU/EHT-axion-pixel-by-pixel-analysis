#!/usr/bin/env python3
"""Standalone official-inspired FITS sensitivity study, NOT an official Top Set reproduction.

No inference is performed. The 4800 hyperparameter images are not the paper's
1000 D-term realizations. Pickle fit statistics precede the last 5 uas restoration;
thresholding them does NOT validate the restored FITS. All outputs are new stages.
"""
from __future__ import annotations

import argparse
import gc
import json
import pickle
from pathlib import Path
import sys

FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
METRICS = ("vis", "cphase", "logcamp", "pvis", "m")
DAYS = (95, 96, 100, 101)
VARIANTS = ("all4800", "legacy_evpa", "datafit_perday", "datafit_shared", "bestI")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--config", required=True, type=Path)
    result.add_argument("--output-root", required=True, type=Path)
    result.add_argument("--trust-pickle", action="store_true", help="trusted executable metadata only")
    result.add_argument("--max-candidates", type=int, help="partial/test study, never a full-ensemble claim")
    result.add_argument("--study-id", help="optional label; does not change partial/full status")
    return result


def load_config(path):
    import yaml

    cfg = yaml.safe_load(Path(path).read_text())
    if not isinstance(cfg, dict) or set(map(int, cfg["days"])) != set(DAYS):
        raise ValueError("config must explicitly declare all four days: 95,96,100,101")
    cfg["days"] = {int(day): row for day, row in cfg["days"].items()}
    for row in cfg["days"].values():
        for key in ("fits_dir", "metadata"):
            if not Path(row[key]).is_absolute():
                raise ValueError(f"{key} must be an explicit absolute input path")
    if not Path(cfg["legacy_selection"]).is_absolute():
        raise ValueError("legacy_selection must be an explicit absolute input path")
    if set(cfg["thresholds"]) != set(METRICS):
        raise ValueError("thresholds must declare vis,cphase,logcamp,pvis,m")
    if any(not 0 < float(v) < float("inf") for v in cfg["thresholds"].values()):
        raise ValueError("thresholds must be finite and positive")
    if not 0 <= int(cfg["max_shift_pixels"]) <= 8 or int(cfg["max_shift_pixels"]) != cfg["max_shift_pixels"]:
        raise ValueError("max_shift_pixels must be an integer in [0,8]")
    if float(cfg["additional_blur_fwhm_uas"]) != 19:
        raise ValueError("this study requires an additional 19 uas FWHM kernel")
    return cfg


def identity(row):
    values = tuple(float(row[key]) for key in FIELDS)
    if any(not -float("inf") < value < float("inf") for value in values):
        raise ValueError("nonfinite metadata identity")
    return values


def candidate_path(directory, key):
    return Path(directory) / ("_".join(str(float(value)) for value in key) + ".fits")


def read_metadata(path):
    import numpy as np

    with Path(path).open("rb") as stream:
        raw = pickle.load(stream)  # Only called after explicit trust acknowledgement.
    rows = raw["meta_data"]
    small = {}
    for row in rows:
        key = identity(row)
        metrics = np.asarray([row[name] for name in METRICS], dtype=float)
        chi2 = float(row["chi2"])
        if key in small or not np.isfinite(metrics).all() or not np.isfinite(chi2):
            raise ValueError(f"{path}: duplicate identity or nonfinite fit statistic: {key}")
        if not np.isclose(chi2, metrics[:3].mean(), rtol=1e-6, atol=1e-10):
            raise ValueError(f"{path}: chi2 is not (vis+cphase+logcamp)/3")
        small[key] = (metrics, chi2)
    del rows, raw
    gc.collect()
    return small


def selections(keys, metadata, legacy, thresholds):
    import numpy as np

    fit = {day: np.asarray([all(metadata[day][k][0][j] < thresholds[m]
                                for j, m in enumerate(METRICS)) for k in keys]) for day in DAYS}
    shared = np.logical_and.reduce(list(fit.values()))
    result = {}
    for day in DAYS:
        result[day] = {"all4800": np.arange(len(keys)),
                       "legacy_evpa": np.flatnonzero([k in legacy for k in keys]),
                       "datafit_perday": np.flatnonzero(fit[day]),
                       "datafit_shared": np.flatnonzero(shared),
                       "bestI": np.asarray([min(range(len(keys)), key=lambda i: metadata[day][keys[i]][1])])}
    return result


def translated(image, offset):
    from eht_axion.research.polarimetry import shift_stokes

    return shift_stokes(image, offset)


def translation(reference, image, limit):
    """Owner's bounded mean-subtracted integer I registration; report cosine Icorr."""
    import numpy as np
    from eht_axion.research.polarimetry import estimate_translation

    offset = np.asarray(estimate_translation(image, reference, maxshiftpx=limit))
    moved = translated(np.repeat(image[..., None], 4, axis=-1), offset)[..., 0]
    denom = np.linalg.norm(reference) * np.linalg.norm(moved)
    return offset, float(np.sum(reference * moved) / denom) if denom else 0.0


def sky_curves(image, operator, area):
    """Same quadrature/observable as RAPTOR 1D, expressed as sky PA and sky EVPA."""
    import numpy as np
    from eht_axion.polarization import integrate_stokes

    legacy_basis = image.transpose(1, 0, 2).copy()
    legacy_basis[..., 1:3] *= -1
    angle = np.roll(integrate_stokes(legacy_basis, operator, "intensity"), -90)
    angle = wrap(angle - np.pi / 2)
    brightness = operator.dot(np.hypot(legacy_basis[..., 1], legacy_basis[..., 2]).ravel()) / area
    return angle, np.roll(brightness, -90)


def bright_mask(image):
    import numpy as np

    intensity, pol = image[..., 0], np.hypot(image[..., 1], image[..., 2])
    return ((intensity >= 0.1 * intensity.max()) & (pol >= 0.2 * pol.max())
            & (intensity > 0) & (pol > 0))


def wrap(angle):
    import numpy as np

    return (angle + np.pi / 2) % np.pi - np.pi / 2


def plot_maps(path, maps, east, north, area):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(1, 4, figsize=(15, 4), constrained_layout=True)
    x, y = np.meshgrid(east, north)
    for ax, day, image in zip(axes, DAYS, maps):
        pol = np.hypot(image[..., 1], image[..., 2]) / area * 1e6
        angle = 0.5 * np.arctan2(image[..., 2], image[..., 1])
        mesh = ax.pcolormesh(x, y, image[..., 0] / area * 1e6, cmap="gray_r", shading="nearest")
        mask = bright_mask(image)
        mask[1::2] = False
        mask[:, 1::2] = False
        fraction = np.divide(np.hypot(image[..., 1], image[..., 2]), image[..., 0],
                             out=np.zeros_like(pol), where=image[..., 0] > 0)
        ax.quiver(x[mask], y[mask], pol[mask]/3*np.sin(angle[mask]), pol[mask]/3*np.cos(angle[mask]),
                  fraction[mask], cmap="turbo", clim=(0, .3), angles="xy", scale_units="xy", scale=1,
                  pivot="middle", headwidth=0, headlength=0, headaxislength=0)
        ax.contour(x, y, pol, levels=[5, 10, 20], colors="cyan", linewidths=.5,
                   linestyles=[":", "--", "-"])
        ax.set(xlim=(60, -60), ylim=(-60, 60), aspect="equal",
               title=f"Apr {day-90}: research", xlabel="East (uas)")
    axes[0].set_ylabel("North (uas)")
    fig.colorbar(mesh, ax=list(axes), label="I (microJy/uas²); 10 uas tick = 30 microJy/uas² P; color=m")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_difference(path, delta, mask, east, north):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    mesh = ax.pcolormesh(east, north, np.ma.array(np.rad2deg(delta), mask=~mask),
                         cmap="twilight", vmin=-90, vmax=90, shading="nearest")
    ax.set(xlim=(max(east), min(east)), ylim=(min(north), max(north)), aspect="equal",
           title="Research: later minus earlier EVPA", xlabel="East (uas)", ylabel="North (uas)")
    fig.colorbar(mesh, ax=ax, label="Delta chi (degrees)")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run(args):
    import numpy as np

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from eht_axion.pipeline.common import stage_directory
    from eht_axion.polarization import aperture_operator

    if not args.trust_pickle:
        raise ValueError("--trust-pickle required: pickle can execute arbitrary code")
    if args.output_root.exists() or args.output_root.is_symlink():
        raise FileExistsError(f"refusing existing output root: {args.output_root}")
    if args.max_candidates is not None and args.max_candidates < 1:
        raise ValueError("--max-candidates must be positive")
    from eht_axion.research.polarimetry import blur_stokes_fwhm, load_stokes_fits

    cfg = load_config(args.config)
    metadata = {day: read_metadata(cfg["days"][day]["metadata"]) for day in DAYS}
    keys = sorted(metadata[DAYS[0]])
    if any(set(metadata[d]) != set(keys) for d in DAYS):
        raise ValueError("four-day metadata identity sets differ")
    if len(keys) != cfg["expected_candidates"]:
        raise ValueError(f"expected {cfg['expected_candidates']} candidates, found {len(keys)}")
    with np.load(cfg["legacy_selection"], allow_pickle=False) as legacy_file:
        legacy = {tuple(float(v) for v in json.loads(str(s))) for s in legacy_file["selected_identity"]}
    if not legacy <= set(keys) or any(len(k) != len(FIELDS) for k in legacy):
        raise ValueError("legacy selected_identity must match exact seven-field tuples")
    full_count = len(keys)
    keys = keys[:args.max_candidates] if args.max_candidates is not None else keys
    chosen = selections(keys, metadata, legacy, cfg["thresholds"])
    reference_key = min(metadata[101], key=lambda k: metadata[101][k][1])
    reference_path = candidate_path(cfg["days"][101]["fits_dir"], reference_key)
    paths = {day: [candidate_path(cfg["days"][day]["fits_dir"], k) for k in keys] for day in DAYS}
    for day in DAYS:
        expected = {candidate_path(cfg["days"][day]["fits_dir"], k).name for k in metadata[day]}
        actual = {p.name for p in Path(cfg["days"][day]["fits_dir"]).glob("*.fits")}
        if actual != expected:
            raise ValueError(f"day {day}: FITS filename mismatch; missing={sorted(expected-actual)[:3]}, "
                             f"unexpected={sorted(actual-expected)[:3]}")
    inputs = [Path(__file__).resolve(), args.config.resolve(), Path(cfg["legacy_selection"])]
    inputs += [Path(cfg["days"][d]["metadata"]) for d in DAYS]
    inputs += sorted({reference_path, *(p for daily in paths.values() for p in daily)})
    import astropy

    config = dict(cfg, study_id=args.study_id, max_candidates=args.max_candidates,
                  astropy_version=astropy.__version__)
    with stage_directory(args.output_root, "research:eht-polarization", config, inputs) as out:
        reference_data = load_stokes_fits(reference_path)
        east, north = (np.asarray(reference_data[k]) for k in ("east_uas", "north_uas"))
        native_reference = np.asarray(reference_data["stokes"])
        size = native_reference.shape[0]
        dx, dy = abs(float(np.diff(east)[0])), abs(float(np.diff(north)[0]))
        area = dx * dy
        if str(reference_data["bunit"]).upper().replace(" ", "") != "JY/PIXEL":
            raise ValueError("native FITS brightness must be Jy/pixel; no implicit unit conversion")
        if not (np.allclose(east, (np.arange(64)-31.5)*-2)
                and np.allclose(north, (np.arange(64)-31.5)*2)):
            raise ValueError("centered CRPIX=32.5 required for the fixed aperture quadrature")
        if (native_reference.shape != (64, 64, 4) or not np.allclose(np.diff(east), -dx)
                or not np.allclose(np.diff(north), dy) or not np.isclose(dx, 2) or not np.isclose(dy, 2)):
            raise ValueError("study requires native 64x64, east decreasing/north increasing, 2 uas pixels")
        def blur(image):
            return blur_stokes_fwhm(image, 19, (dy, dx))

        reference = blur(native_reference)[..., 0]
        operator = aperture_operator(shape=size, fov=128, radial_max=64, window_degrees=10, angular_bins=360)
        report = {"research_only": True, "official_reproduction": False, "partial": len(keys) < full_count,
                  "study_id": args.study_id, "full_candidates_per_day": full_count,
                  "processed_per_day": len(keys), "legacy_identity_count": len(legacy),
                  "bestI_scope": "minimum I chi2 among processed rows; reference always full Apr11 bestI",
                  "selection_identity_fields": FIELDS, "variants": VARIANTS,
                  "reference_identity": reference_key, "reference_day": 101,
                  "metadata_stage": "before final 5 uas restoration; NOT post-restoration validation",
                  "blur": "additional 19 uas FWHM kernel; native 5 uas restoration retained; NOT final beam",
                  "array_modes": ["raw_blurred_unregistered", "aligned_blurred"],
                  "basis": "native celestial IQUV Jy/pixel; no Q/U flips in saved maps",
                  "curve_basis": "native.T, Q/U negated, existing intensity aperture, roll -90; EVPA minus90",
                  "curve_units": "chi radians sky E of N; P scalar aperture integral Jy; PA degrees E of N",
                  "threshold_claim": "strict per-metric diagnostic cuts, NOT official Top Set",
                  "ensemble_claim": "hyperparameter ensemble, NOT 1000 D-term realizations",
                  "padding": "zero; integer I-only registration; aligned maps include extra blur",
                  "registration_objective": "maximum mean-subtracted I cross-correlation within bounds",
                  "Icorr_definition": "cosine similarity after zero-padded translation, no mean subtraction",
                  "artifact_schema": {"aligned_stokes": "{day}_aligned_stokes.npy: float32 N,Y,X,IQUV",
                      "aperture": "{day}_aperture.npy: float64 N,360 aligned sky EVPA radians",
                      "indices": "{day}_{variant}_indices.npy: canonical identity row indices",
                      "modes": "axis0 raw/aligned; raw means blurred but unregistered",
                      "summary": "mean_stokes_curve_chi is NOT circular_mean_chi of individual images"},
                  "days": {}}
        means = {}
        np.save(out / "identities.npy", np.asarray(keys))
        for day in DAYS:
            count = len(keys)
            maps = np.lib.format.open_memmap(out / f"{day}_aligned_stokes.npy", mode="w+",
                                            dtype="float64", shape=(count, size, size, 4))
            angles, brightness = np.empty((2, count, 360)), np.empty((2, count, 360))
            offsets, correlations, drift = np.empty((count, 2), int), np.empty(count), np.empty((count, 2))
            blur_drift = np.empty((count, 2))
            sums = np.zeros((len(VARIANTS), 2, size, size, 4))
            memberships = [{int(i) for i in chosen[day][v]} for v in VARIANTS]
            for i, path in enumerate(paths[day]):
                loaded = reference_data if path == reference_path else load_stokes_fits(path)
                if not (np.allclose(loaded["east_uas"], east) and np.allclose(loaded["north_uas"], north)):
                    raise ValueError(f"inconsistent native WCS: {path}")
                if str(loaded["bunit"]).upper().replace(" ", "") != "JY/PIXEL":
                    raise ValueError(f"expected native Jy/pixel: {path}")
                native = np.asarray(loaded["stokes"])
                image = blur(native)
                original_flux = np.array([native[..., 0].sum(),
                                          np.hypot(native[..., 1], native[..., 2]).sum()])
                blurred_flux = np.array([image[..., 0].sum(), np.hypot(image[..., 1], image[..., 2]).sum()])
                blur_drift[i] = np.divide(blurred_flux-original_flux, original_flux,
                                          out=np.zeros(2), where=original_flux != 0) * 100
                offsets[i], correlations[i] = translation(reference, image[..., 0], cfg["max_shift_pixels"])
                aligned = translated(image, offsets[i])
                maps[i] = aligned
                for mode, data in enumerate((image, aligned)):
                    angles[mode, i], brightness[mode, i] = sky_curves(data, operator, area)
                    for j, membership in enumerate(memberships):
                        if i in membership:
                            sums[j, mode] += data
                before = np.array([image[..., 0].sum(), np.hypot(image[..., 1], image[..., 2]).sum()])
                after = np.array([aligned[..., 0].sum(), np.hypot(aligned[..., 1], aligned[..., 2]).sum()])
                drift[i] = np.divide(after-before, before, out=np.zeros(2), where=before != 0) * 100
            maps.flush()
            del maps
            np.save(out / f"{day}_aperture.npy", angles[1])
            np.savez(out / f"{day}_candidates.npz", chi=angles, p=brightness, offsets_yx=offsets,
                     intensity_correlation=correlations, shift_flux_drift_percent=drift,
                     blur_flux_drift_percent=blur_drift,
                     fit_metrics=np.asarray([metadata[day][k][0] for k in keys]),
                     chi2=np.asarray([metadata[day][k][1] for k in keys]), metric_names=METRICS)
            daily = {"counts": {v: len(chosen[day][v]) for v in VARIANTS}, "mean_mask_fractions": {},
                     "offset_bound_fraction": float(np.mean(
                         np.any(abs(offsets) == cfg["max_shift_pixels"], axis=1))),
                     "offset_yx_quantiles": np.quantile(offsets, [0, .5, 1], axis=0).tolist(),
                     "Icorr_quantiles": np.quantile(correlations, [0, .5, 1]).tolist(),
                     "shift_flux_drift_percent_quantiles_I_P":
                         np.quantile(drift, [0, .5, 1], axis=0).tolist(),
                     "blur_flux_drift_percent_quantiles_I_P":
                         np.quantile(blur_drift, [0, .5, 1], axis=0).tolist()}
            for j, variant in enumerate(VARIANTS):
                indices = chosen[day][variant]
                np.save(out / f"{day}_{variant}_indices.npy", indices)
                if not len(indices):
                    continue  # Empty shared intersections are reported, NEVER relaxed.
                mean = sums[j] / len(indices)
                means[day, variant] = mean
                daily["mean_mask_fractions"][variant] = [float(bright_mask(m).mean()) for m in mean]
                resultant = np.nanmean(np.exp(2j * angles[:, indices]), axis=1)
                mean_curves = [sky_curves(image, operator, area) for image in mean]
                np.savez(out / f"{day}_{variant}_summary.npz", mean_stokes=mean,
                         east_uas=east, north_uas=north,
                         mean_stokes_curve_chi=np.asarray([curve[0] for curve in mean_curves]),
                         mean_stokes_curve_p=np.asarray([curve[1] for curve in mean_curves]),
                         circular_mean_chi=wrap(0.5*np.angle(resultant)),
                         circular_std_chi=np.sqrt(-.5*np.log(np.clip(abs(resultant), 1e-300, 1))),
                         mean_p=np.mean(brightness[:, indices], axis=1))
            report["days"][str(day)] = daily
            print(f"day {day}: {daily['counts']}", flush=True)
        for variant in VARIANTS:
            if not all((day, variant) in means for day in DAYS):
                continue
            for mode, label in enumerate(("raw", "aligned")):
                daily_maps = [means[day, variant][mode] for day in DAYS]
                plot_maps(out / f"{variant}_{label}_means.png", daily_maps, east, north, area)
                for earlier, later in ((95, 96), (100, 101)):
                    a, b = (means[d, variant][mode] for d in (earlier, later))
                    delta = wrap(.5*np.angle((b[..., 1]+1j*b[..., 2]) * (a[..., 1]-1j*a[..., 2])))
                    mask = bright_mask(a) & bright_mask(b)
                    ca, pa = sky_curves(a, operator, area)
                    cb, pb = sky_curves(b, operator, area)
                    np.savez(out / f"{variant}_{label}_{earlier}_{later}_daydiff.npz",
                             delta_chi=delta, mask=mask, curve_delta_chi=wrap(cb-ca),
                             earlier_p=pa, later_p=pb, pa_degrees=np.arange(360))
                    plot_difference(out / f"{variant}_{label}_{earlier}_{later}_daydiff.png",
                                    delta, mask, east, north)
                    report[f"{variant}_{label}_{earlier}_{later}"] = {
                        "mask_fraction": float(mask.mean()),
                        "delta_chi_deg_quantiles": np.rad2deg(np.quantile(delta[mask], [0, .5, 1])).tolist()
                        if mask.any() else [], "masked_pixels": int(mask.sum())}
        (out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main(argv=None):
    args = parser().parse_args(argv)  # --help works without bootstrap, numpy, or astropy.
    if not args.trust_pickle:
        parser().error("--trust-pickle is required")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
