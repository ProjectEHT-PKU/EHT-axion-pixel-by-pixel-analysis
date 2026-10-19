"""Notebook-free EHT selection, fixed relocation and ensemble summaries.

Only trusted historical pickles are accepted.  This module deliberately leaves
cache/transaction management to its caller; report.json is written last.
"""

from __future__ import annotations

import json
import math
import pickle
import tempfile
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap

from eht_axion.angles import circular_difference
from eht_axion.stokes import finite_evpa

from .layout import date_label, input_eht_file

DAYS = (95, 96, 100, 101)
FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
RANK_OFFSETS = {95: (2, 0), 96: (1, 1), 100: (-4, 4), 101: (-1, 1)}
FINAL_OFFSETS = {95: (1, 0), 96: (0, 1), 100: (-1, 2), 101: (0, 0)}
HISTORICAL_RECTANGLES = ((48, 55, 10, 19),)


def validate_geometry_contract(
    expected_shape: int | None = None,
    angular_fov_uas: float | None = None,
    pixel_size_uas: float | None = None,
) -> None:
    """Validate declared angular geometry without IO or a mass/distance model.

    All three values must be supplied together, or all omitted for the legacy
    shape-inferred API. Angular sizes are in microarcseconds, NOT Rg. The caller
    separately converts angular FOV to the explicit ``run_eht_stage(fov=...)``
    value in Rg; this helper neither guesses nor recalibrates that conversion.
    """
    supplied = (expected_shape is not None, angular_fov_uas is not None,
                pixel_size_uas is not None)
    if not any(supplied):
        return
    if not all(supplied):
        raise ValueError("expected_shape, angular_fov_uas and pixel_size_uas "
                         "must be all supplied or all None")
    if (isinstance(expected_shape, (bool, np.bool_))
            or not isinstance(expected_shape, (int, np.integer)) or expected_shape < 2):
        raise ValueError("expected_shape must be an integer >= 2 (not bool or float)")
    for name, value in (("angular_fov_uas", angular_fov_uas), ("pixel_size_uas", pixel_size_uas)):
        if (isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, float, np.integer, np.floating))):
            raise ValueError(f"{name} must be finite and positive (not bool or string)")
        try:
            valid = math.isfinite(value) and value > 0
        except OverflowError:
            valid = False
        if not valid:
            raise ValueError(f"{name} must be finite and positive")
    try:
        extent = int(expected_shape) * float(pixel_size_uas)
        consistent = math.isfinite(extent) and math.isclose(
            extent, float(angular_fov_uas), rel_tol=1e-12, abs_tol=0.0)
    except OverflowError:
        consistent = False
    if not consistent:
        raise ValueError("expected_shape * pixel_size_uas must equal angular_fov_uas")


def _check_declared_shape(raw: Mapping, expected_shape: int | None, path: Path) -> None:
    """Check every supplied Stokes channel before any image calculation."""
    if expected_shape is None:
        return
    for field in ("I_data", "Q_data", "U_data", "V_data"):
        if field == "V_data" and field not in raw:
            continue
        shape = np.shape(raw[field])
        if shape[:2] != (expected_shape, expected_shape):
            raise ValueError(f"{path}: {field} raw shape {shape[:2]} != "
                             f"expected_shape ({expected_shape}, {expected_shape})")


def metadata_identities(metadata: Sequence[Mapping]) -> list[tuple]:
    """Return unique, finite scalar FIELDS tuples (never match by row number)."""
    identities = []
    for index, row in enumerate(metadata):
        values = []
        for field in FIELDS:
            if field not in row:
                raise ValueError(f"metadata row {index} missing {field!r}")
            value = row[field]
            if isinstance(value, np.generic):
                value = value.item()
            if not isinstance(value, (str, int, float, bool)):
                raise ValueError(f"metadata {field!r} must be a scalar")
            if isinstance(value, (int, float)) and not np.isfinite(value):
                raise ValueError(f"metadata {field!r} must be finite")
            values.append(value)
        identities.append(tuple(values))
    if not identities or len(set(identities)) != len(identities):
        raise ValueError("metadata identities must be nonempty and unique")
    return identities


def align_metadata(
    metadata: Sequence[Mapping], reference: Sequence[tuple] | None = None
) -> tuple[list[tuple], np.ndarray]:
    """Validate the identity set and map sorted common identity -> raw index."""
    identities = metadata_identities(metadata)
    try:
        common = sorted(identities) if reference is None else list(reference)
    except TypeError as error:
        raise ValueError("metadata FIELDS values are not mutually sortable") from error
    if set(identities) != set(common):
        missing = set(common) - set(identities)
        extra = set(identities) - set(common)
        raise ValueError(f"metadata identity mismatch: {len(missing)} missing, {len(extra)} extra")
    lookup = {identity: index for index, identity in enumerate(identities)}
    return common, np.asarray([lookup[identity] for identity in common], dtype=np.int64)


def shift_spatial(array: np.ndarray, offset: tuple[int, int]) -> np.ndarray:
    """NaN-padded translation on the last two axes, not np.roll.

    output[..., x, y] = input[..., x + offset[0], y + offset[1]].
    This is the historical padded-slicing convention, including its sign.
    """
    array = np.asarray(array)
    if array.ndim < 2 or len(offset) != 2 or any(int(x) != x for x in offset):
        raise ValueError("translation requires two integer spatial offsets")
    result = np.full(array.shape, np.nan, dtype=np.result_type(array.dtype, np.float64))
    source, target = [], []
    for length, delta in zip(array.shape[-2:], offset):
        delta = int(delta)
        if abs(delta) >= length:
            return result
        source.append(slice(max(delta, 0), min(length + delta, length)))
        target.append(slice(max(-delta, 0), min(length - delta, length)))
    result[(..., *target)] = array[(..., *source)]
    return result


def transform_stokes(raw: Mapping, indices: np.ndarray) -> np.ndarray:
    """Convert a small raw (H,W,N) batch to (N,H,W,I/Q/U/V).

    Flip raw axis 0, then rot90(k=-1).  I has NO sign change; Q/U
    change sign.  Historical EHT input has no V: its output channel is zero.
    If V_data is provided, preserve its sign and apply the same spatial map.
    """
    shape = np.asarray(raw["I_data"]).shape[0]
    result = np.zeros((len(indices), shape, shape, 4), dtype=np.float64)
    for channel, field in enumerate(("I_data", "Q_data", "U_data", "V_data")):
        if field not in raw:
            if channel == 3:
                continue
            raise ValueError(f"missing {field}")
        image = np.asarray(raw[field])[::-1, :, indices].transpose(2, 0, 1)
        if not np.isfinite(image).all():
            raise ValueError(f"{field} contains nonfinite input pixels")
        result[..., channel] = np.rot90(image, k=-1, axes=(1, 2))
        if channel in (1, 2):
            result[..., channel] *= -1
    return result


def stokes_evpa(stokes: np.ndarray) -> np.ndarray:
    """Pixel EVPA in radians; absent/zero polarization is undefined, not zero."""
    return finite_evpa(stokes[..., 1], stokes[..., 2])


def summarize_evpa(
    theta: np.ndarray, *, centering: str = "ensemble", pixel_batch_size: int = 128
) -> dict[str, np.ndarray]:
    """Full-ensemble circular branch, then ordinary mean/std (ddof=0).

    The exposed center and unwrapped_mean MUST also be used for residuals and
    covariance; wrapping residuals a second time around mean changes the branch.
    NaNs (translation padding / empty apertures) are excluded, with count saved.
    Incomplete features have mean/std=NaN so downstream finite checks cannot
    silently use a different, sample-dependent ensemble. A normalized resultant
    within eight float64 epsilons of zero has no numerically resolvable branch:
    all branch/spread statistics are NaN, retaining concentration/count only.
    This is a roundoff guard, not a physical concentration selection cut.
    """
    if centering != "ensemble":
        raise ValueError("only centering='ensemble' is implemented; no historical fit fallback")
    if theta.ndim < 2 or theta.shape[0] == 0 or pixel_batch_size <= 0:
        raise ValueError("theta must have a nonempty ensemble and positive pixel batch size")
    flat = theta.reshape(theta.shape[0], -1)
    keys = ("mean", "std", "center", "unwrapped_mean", "concentration", "circular_std")
    result = {key: np.full(flat.shape[1], np.nan) for key in keys}
    count = np.zeros(flat.shape[1], dtype=np.int64)
    for start in range(0, flat.shape[1], pixel_batch_size):
        end = min(start + pixel_batch_size, flat.shape[1])
        values = np.asarray(flat[:, start:end], dtype=np.float64)
        finite = np.isfinite(values)
        n = finite.sum(axis=0)
        denominator = np.where(n > 0, n, 1)
        vector = np.where(finite, np.exp(2j * np.where(finite, values, 0)), 0).sum(axis=0) / denominator
        center = np.angle(vector) / 2
        unwrapped = center + circular_difference(values, center)
        average = np.where(finite, unwrapped, 0).sum(axis=0) / denominator
        variance = np.where(finite, (unwrapped - average) ** 2, 0).sum(axis=0) / denominator
        concentration = np.clip(np.abs(vector), 0, 1)
        resolved = (n > 0) & (concentration > 8 * np.finfo(np.float64).eps)
        with np.errstate(divide="ignore"):
            circular_std = np.sqrt(-2 * np.log(concentration)) / 2
        parts = (circular_difference(average, 0), np.sqrt(variance), center,
                 average, concentration, circular_std)
        for key, part in zip(keys, parts):
            valid = (n > 0) if key == "concentration" else resolved
            result[key][start:end] = np.where(valid, part, np.nan)
        count[start:end] = n
    for key in ("mean", "std"):
        result[key][count != len(theta)] = np.nan
    result["count"] = count
    return {key: value.reshape(theta.shape[1:]) for key, value in result.items()}


def ensemble_residuals(theta: np.ndarray, summary: Mapping[str, np.ndarray]) -> np.ndarray:
    """Residuals on the saved full-ensemble branch, suitable for covariance."""
    return (summary["center"] + circular_difference(theta, summary["center"])
            - summary["unwrapped_mean"])


def ensemble_covariance(
    theta: np.ndarray, summary: Mapping[str, np.ndarray], *, batch_size: int = 64
) -> np.ndarray:
    """Population covariance for complete, resolved features; invalid stay NaN.

    No pairwise/recentered estimator is used.  Intended for the 360-bin summary
    or explicitly preselected 2-D pixels, not all 4096 pixels at once.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    features = int(np.prod(theta.shape[1:]))
    complete = np.asarray(summary["count"]).reshape(-1) == theta.shape[0]
    for key in ("mean", "std", "center", "unwrapped_mean"):
        complete &= np.isfinite(np.asarray(summary[key]).reshape(-1))
    indices = np.flatnonzero(complete)
    accumulated = np.zeros((len(indices), len(indices)))
    for start in range(0, len(theta), batch_size):
        residual = ensemble_residuals(theta[start:start + batch_size], summary)
        residual = residual.reshape(len(residual), features)[:, indices]
        accumulated += residual.T @ residual
    result = np.full((features, features), np.nan)
    result[np.ix_(indices, indices)] = accumulated / len(theta)
    return result


def build_mask(
    intensity: np.ndarray,
    std: np.ndarray,
    *,
    fov: float = 32.6868,
    min_radius: float = 2.0,
    intensity_fraction: float = 0.1,
    historical_mask: bool = True,
    exclude_columns: int = 5,
    exclusion_rectangles: Sequence[tuple[int, int, int, int]] = HISTORICAL_RECTANGLES,
    require_positive_std: bool = False,
) -> np.ndarray:
    """Four-day intersection plus radius and optional historical exclusions.

    Rectangle bounds are half-open and clipped to the image, so small synthetic
    images are supported.  The final mask allows std=0; ranking cannot divide by it.
    """
    intensity, std = np.asarray(intensity), np.asarray(std)
    if intensity.ndim != 3 or intensity.shape != std.shape:
        raise ValueError("intensity/std must be matching (day,height,width) arrays")
    height, width = intensity.shape[1:]
    if height != width or fov <= 0 or min_radius < 0 or not 0 <= intensity_fraction < 1:
        raise ValueError("invalid mask geometry or intensity fraction")
    finite_i = np.where(np.isfinite(intensity), intensity, -np.inf)
    peak = np.max(finite_i, axis=(1, 2))
    valid = np.isfinite(intensity) & np.isfinite(std)
    valid &= intensity > intensity_fraction * peak[:, None, None]
    if require_positive_std:
        valid &= std > 0
    mask = np.all(valid, axis=0)
    x, y = np.indices((height, width))
    radius = np.hypot(x - (height - 1) / 2, y - (width - 1) / 2) * fov / height
    mask &= radius >= min_radius
    if historical_mask:
        if exclude_columns < 0:
            raise ValueError("exclude_columns must be nonnegative")
        mask[:, :exclude_columns] = False
        for x0, x1, y0, y1 in exclusion_rectangles:
            if not 0 <= x0 <= x1 or not 0 <= y0 <= y1:
                raise ValueError("exclusion rectangles require ordered nonnegative bounds")
            mask[x0:min(x1, height), y0:min(y1, width)] = False
    return mask


def _load_day(path: Path) -> dict:
    # The public entry point enforces trust_pickle before reaching this function.
    with path.open("rb") as stream:
        raw = pickle.load(stream)
    if not isinstance(raw, Mapping) or "meta_data" not in raw:
        raise ValueError(f"{path}: missing metadata mapping")
    n = len(raw["meta_data"])
    expected = None
    for field in ("I_data", "Q_data", "U_data", "V_data"):
        if field == "V_data" and field not in raw:
            continue
        if field not in raw:
            raise ValueError(f"{path}: missing {field}")
        array = np.asarray(raw[field])
        if (array.ndim != 3 or array.shape[0] < 2 or array.shape[0] != array.shape[1]
                or array.shape[2] != n):
            raise ValueError(f"{path}: {field} must have shape (H,H,{n})")
        if expected is not None and array.shape != expected:
            raise ValueError(f"{path}: Stokes shapes differ")
        expected = array.shape
    return raw


def _diagnostics(summary: Mapping[str, np.ndarray]) -> dict:
    result = {"finite_features": int(np.isfinite(summary["mean"]).sum())}
    for key in ("std", "concentration", "circular_std"):
        values = summary[key][np.isfinite(summary[key])]
        result[key] = ({"min": float(values.min()), "median": float(np.median(values)),
                        "max": float(values.max())} if values.size else None)
    return result


def run_eht_stage(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    trust_pickle: bool = False,
    topset_size: int = 1572,
    centering: str = "ensemble",
    weighting: str = "intensity",
    inherit_2d_mask: bool = False,
    fov: float = 32.6868,
    expected_shape: int | None = None,
    angular_fov_uas: float | None = None,
    pixel_size_uas: float | None = None,
    radial_max: float | None = None,
    window_degrees: int = 10,
    angular_bins: int = 360,
    angle_roll: int = -90,
    min_radius: float = 2.0,
    intensity_fraction: float = 0.1,
    historical_mask: bool = True,
    exclude_columns: int = 5,
    exclusion_rectangles: Sequence[tuple[int, int, int, int]] = HISTORICAL_RECTANGLES,
    batch_size: int = 64,
    pixel_batch_size: int = 128,
    rank_offsets: Mapping[int, tuple[int, int]] | None = None,
    final_offsets: Mapping[int, tuple[int, int]] | None = None,
    operator=None,
    integrator: Callable | None = None,
) -> dict:
    """Run four-day identity-matched ranking into the supplied empty EHT dir.

    Each ISO-date directory contains stokes.npy, evpa-1d.npy and summary.npz.
    The caller locates the EHT directory; this function does not resolve run roots.

    Injecting an operator/integrator is an explicit API for tests/alternate
    quadrature.  Otherwise use polarization.aperture_operator/integrate_stokes;
    missing parent APIs are an error, never an implicit historical fallback.
    All persisted arrays are float64 (counts/indices int64, mask bool).

    ``expected_shape``, ``angular_fov_uas`` and ``pixel_size_uas`` form an
    optional all-or-none declaration, checked before IO and against each raw
    Stokes image before calculations. ``fov`` remains a separate explicit Rg
    extent for masks/apertures, not the angular extent. No M/D is required here.
    """
    validate_geometry_contract(expected_shape, angular_fov_uas, pixel_size_uas)
    if expected_shape is not None:
        expected_shape = int(expected_shape)
        angular_fov_uas, pixel_size_uas = float(angular_fov_uas), float(pixel_size_uas)
    if not trust_pickle:
        raise ValueError("pickle can execute code: explicitly set trust_pickle=True / --trust-pickle")
    if centering != "ensemble" or weighting not in ("intensity", "flux"):
        raise ValueError("centering must be 'ensemble'; weighting must be 'intensity' or 'flux'")
    for name, value in (("topset_size", topset_size), ("batch_size", batch_size),
                        ("pixel_batch_size", pixel_batch_size), ("angular_bins", angular_bins)):
        if (isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer)) or value <= 0):
            raise ValueError(f"{name}: sizes must be positive integers")
    if not np.isfinite(fov) or fov <= 0:
        raise ValueError("fov must be finite and positive")
    rank_offsets = dict(RANK_OFFSETS if rank_offsets is None else rank_offsets)
    final_offsets = dict(FINAL_OFFSETS if final_offsets is None else final_offsets)
    if set(rank_offsets) != set(DAYS) or set(final_offsets) != set(DAYS):
        raise ValueError("fixed offsets must specify exactly all four days")
    from eht_axion import polarization

    if integrator is None:
        integrator = polarization.integrate_stokes
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    from .common import reject_directory_symlinks

    reject_directory_symlinks(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(f"output directory must be empty: {output_dir}")
    mask_options = dict(fov=fov, min_radius=min_radius, intensity_fraction=intensity_fraction,
                        historical_mask=historical_mask, exclude_columns=exclude_columns,
                        exclusion_rectangles=exclusion_rectangles)
    report = {
        "schema_version": 1, "days": list(DAYS), "fields": list(FIELDS),
        "input_dir": str(input_dir.resolve()), "topset_size": topset_size,
        "centering": centering, "ddof": 0, "weighting": weighting,
        "inherit_2d_mask": inherit_2d_mask, "fov": fov, "radial_max": radial_max,
        "expected_shape": expected_shape, "angular_fov_uas": angular_fov_uas,
        "pixel_size_uas": pixel_size_uas,
        "geometry": {"mode": "legacy" if expected_shape is None else "declared",
                     "expected_shape": expected_shape, "angular_fov_uas": angular_fov_uas,
                     "pixel_size_uas": pixel_size_uas, "fov_rg": fov},
        "window_degrees": window_degrees, "angular_bins": angular_bins, "angle_roll": angle_roll,
        "rank_offsets": rank_offsets, "final_offsets": final_offsets,
        "relocation": "fixed; never refitted", "mask_options": mask_options,
        "stokes_channels": ["I", "Q", "U", "V"],
        "coordinate_transform": "flip raw axis 0, rot90 k=-1; negate Q/U only",
        "covariance": "population, saved full-ensemble branch; incomplete/unresolved features are NaN",
        "per_day": {},
    }
    # Temporary angle maps bound RAM independently of the four-day ensemble size.
    with tempfile.TemporaryDirectory(prefix=".eht-work-", dir=output_dir) as work_name:
        work = Path(work_name)
        common = None
        orders, rank_summaries, rank_intensity = [], [], []
        shape = n = None
        for day in DAYS:
            path = input_eht_file(input_dir, day)
            raw = _load_day(path)
            _check_declared_shape(raw, expected_shape, path)
            common, order = align_metadata(raw["meta_data"], common)
            current_shape = raw["I_data"].shape[0]
            if shape is not None and current_shape != shape:
                raise ValueError("four-day image shapes differ")
            shape, n = current_shape, len(order)
            if topset_size > n:
                raise ValueError(f"topset_size={topset_size} exceeds ensemble size {n}")
            orders.append(order)
            theta = open_memmap(work / f"rank_{day}.npy", mode="w+", dtype="float64",
                                shape=(n, shape, shape))
            intensity_sum = np.zeros((shape, shape))
            for start in range(0, n, batch_size):
                stop = min(start + batch_size, n)
                stokes = transform_stokes(raw, order[start:stop])
                intensity_sum += stokes[..., 0].sum(axis=0)
                angles = stokes_evpa(stokes)
                theta[start:stop] = shift_spatial(angles, rank_offsets[day])
            theta.flush()
            rank_summaries.append(summarize_evpa(theta, centering=centering,
                                                pixel_batch_size=pixel_batch_size))
            rank_intensity.append(shift_spatial(intensity_sum / n, rank_offsets[day]))
            report["per_day"][str(day)] = {"V_source": "input" if "V_data" in raw else "zero",
                                           "ranking": _diagnostics(rank_summaries[-1])}
            del raw, theta, stokes, angles
        rank_mask = build_mask(np.stack(rank_intensity),
                               np.stack([s["std"] for s in rank_summaries]),
                               require_positive_std=True, **mask_options)
        if not rank_mask.any():
            raise ValueError("ranking mask is empty (check geometry, exclusions and nonzero spread)")
        score = np.zeros(n)
        for day, summary in zip(DAYS, rank_summaries):
            theta = np.load(work / f"rank_{day}.npy", mmap_mode="r")
            for start in range(0, n, batch_size):
                residual = ensemble_residuals(theta[start:start + batch_size], summary)[:, rank_mask]
                score[start:start + len(residual)] += np.sum(
                    (residual / summary["std"][rank_mask]) ** 2, axis=1)
            del theta
            (work / f"rank_{day}.npy").unlink()
        if not np.isfinite(score).all():
            raise ValueError("nonfinite ranking scores")
        sorted_indices = np.argsort(score, kind="stable")
        selected = sorted_indices[:topset_size]
        identity = np.asarray([json.dumps(key, ensure_ascii=False, separators=(",", ":"))
                               for key in common])
        np.savez(output_dir / "selection.npz", fields=np.asarray(FIELDS), identity=identity,
                 selected_identity=identity[selected], sorted_indices=sorted_indices,
                 selected_indices=selected, raw_sorted_indices=np.stack(orders),
                 selected_raw_indices=np.stack(orders)[:, selected], score=score,
                 rank_mask=rank_mask, days=np.asarray(DAYS))
        if operator is None:
            operator = polarization.aperture_operator(
                shape=shape, fov=fov, radial_max=radial_max,
                window_degrees=window_degrees, angular_bins=angular_bins)
        if operator.shape != (angular_bins, shape * shape):
            raise ValueError("aperture operator shape must be (angular_bins, H*H)")
        summaries_2d, final_intensity = [], []
        for day, order in zip(DAYS, orders):
            path = input_eht_file(input_dir, day)
            raw = _load_day(path)
            _check_declared_shape(raw, expected_shape, path)
            # Detect metadata changes between ranking and final-product passes.
            _, checked_order = align_metadata(raw["meta_data"], common)
            if not np.array_equal(order, checked_order) or raw["I_data"].shape[:2] != (shape, shape):
                raise ValueError("input changed between ranking and final-product passes")
            day_dir = output_dir / date_label(day)
            day_dir.mkdir()
            stokes_out = open_memmap(day_dir / "stokes.npy", mode="w+",
                                     dtype="float64", shape=(topset_size, shape, shape, 4))
            theta = open_memmap(work / f"final_{day}.npy", mode="w+", dtype="float64",
                                shape=(topset_size, shape, shape))
            for start in range(0, topset_size, batch_size):
                stop = min(start + batch_size, topset_size)
                stokes = transform_stokes(raw, order[selected[start:stop]])
                relocated = shift_spatial(stokes.transpose(0, 3, 1, 2), final_offsets[day])
                relocated = relocated.transpose(0, 2, 3, 1)
                stokes_out[start:stop] = relocated
                theta[start:stop] = stokes_evpa(relocated)
            stokes_out.flush()
            theta.flush()
            summary = summarize_evpa(theta, centering=centering, pixel_batch_size=pixel_batch_size)
            summaries_2d.append(summary)
            # Ordinary mean preserves NaN borders, as did the old relocation.
            intensity_sum = np.zeros((shape, shape))
            for start in range(0, topset_size, batch_size):
                intensity_sum += stokes_out[start:start + batch_size, ..., 0].sum(axis=0)
            final_intensity.append(intensity_sum / topset_size)
            report["per_day"][str(day)]["summary_2d"] = _diagnostics(summary)
            del raw, stokes_out, theta, stokes, relocated
            (work / f"final_{day}.npy").unlink()
        mask = build_mask(np.stack(final_intensity),
                          np.stack([s["std"] for s in summaries_2d]), **mask_options)
        np.save(output_dir / "mask.npy", mask)
        for k, day in enumerate(DAYS):
            day_dir = output_dir / date_label(day)
            stokes_out = np.load(day_dir / "stokes.npy", mmap_mode="r")
            evpa = open_memmap(day_dir / "evpa-1d.npy", mode="w+", dtype="float64",
                               shape=(topset_size, angular_bins))
            for start in range(0, topset_size, batch_size):
                # Keep all-channel padding distinct from unknown individual
                # channels; integrate_stokes validates and zeroes padding only.
                batch = np.array(stokes_out[start:start + batch_size], copy=True)
                if inherit_2d_mask:
                    batch[:, ~mask, :] = 0
                integrated = np.asarray(integrator(batch, operator, weighting=weighting))
                if integrated.shape != (len(batch), angular_bins):
                    raise ValueError("integrate_stokes must return (batch, angular_bins) EVPA radians")
                evpa[start:start + len(batch)] = np.roll(integrated, angle_roll, axis=-1)
            evpa.flush()
            summary = summarize_evpa(evpa, centering=centering, pixel_batch_size=pixel_batch_size)
            covariance = ensemble_covariance(evpa, summary, batch_size=batch_size)
            np.savez(
                day_dir / "summary.npz",
                **{f"evpa_{key}_2d": value for key, value in summaries_2d[k].items()},
                **{f"evpa_{key}_1d": value for key, value in summary.items()},
                intensity_mean=final_intensity[k], covariance_1d=covariance,
            )
            report["per_day"][str(day)]["summary_1d"] = _diagnostics(summary)
            del stokes_out, evpa
        report.update(ensemble_size=n, shape=shape, rank_pixels=int(rank_mask.sum()),
                      mask_pixels=int(mask.sum()), score_min=float(score.min()),
                      score_max=float(score.max()), batch_size=batch_size,
                      pixel_batch_size=pixel_batch_size)
    (output_dir / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report
