"""Script-only Figure-4 transport diagnostics, ported from GRMHD/main.ipynb.

No notebook or legacy ImageS is imported/read at runtime. See docs/transport-stage.md
for deliberate corrections, inferred switch semantics, and array contracts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections.abc import Sequence

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.optimize import brentq

from .raptor import fit_phase_dft

__all__ = [
    "get_360_x_y", "get_analytical_A_1d", "get_raptor_evpa", "get_analytical_evpa",
    "analytical_stokes", "mask_and_blur_frame", "legacy_rms_amplitude",
    "run_transport_diagnostics", "main",
]
VARIANTS = ("common", "jonly", "eq_jonly", "fake")
ORDERS = ("mask-before-blur", "blur-before-mask")


def _grid(n, fov):
    if n < 2 or not np.isfinite(fov) or fov <= 0:
        raise ValueError("square grid needs N >= 2 and finite positive fov")
    return (np.arange(n) - (n - 1) / 2) * fov / n


def _square(array, channels=False):
    a = np.asarray(array)
    spatial = a.shape[:-1] if channels else a.shape
    if (len(spatial) != 2 or spatial[0] != spatial[1] or spatial[0] < 2
            or (channels and a.shape[-1] != 4) or a.dtype.kind not in "iuf"):
        raise ValueError("expected real square (N,N,4) Stokes" if channels
                         else "expected real square (N,N) map")
    return a


def _cube(stokes):
    a = np.asarray(stokes)
    if a.ndim != 4 or a.shape[0] not in (16, 17):
        raise ValueError("expected (16|17,N,N,4) phase Stokes, noaxion is state0 if present")
    _square(a[0], channels=True)
    return a


def _wrap(x):
    return (x + np.pi / 2) % np.pi - np.pi / 2


def get_360_x_y(r_coord, r_s, *, fov=46.2262, angular_bins=360,
                crossing="first-outward", scan_samples=2049, return_diagnostics=False):
    """Intersect rays (-rho*sin(varphi), rho*cos(varphi)) with BL r=r_s.

    Search only within the pixel-center square and rho <= 2*r_s (old bracket).
    first-outward scans for all finite sign changes and chooses the innermost;
    legacy-bracket uses only the two endpoints, matching the old brentq call.
    Old NaN->0 is retained *only for radius geometry*, not polarization data.
    Missing intersections are NaN, never nearest-neighbor/extrapolated points.
    """
    radius = _square(r_coord)
    if not np.isfinite(r_s) or r_s <= 0:
        raise ValueError("r_s must be finite and positive")
    if crossing not in ("first-outward", "legacy-bracket"):
        raise ValueError("crossing must be first-outward|legacy-bracket")
    if isinstance(angular_bins, bool) or int(angular_bins) != angular_bins or angular_bins < 4:
        raise ValueError("angular_bins must be an integer >= 4")
    if int(scan_samples) != scan_samples or scan_samples < 3:
        raise ValueError("scan_samples must be an integer >= 3")
    angular_bins, scan_samples = int(angular_bins), int(scan_samples)
    grid = _grid(len(radius), fov)
    interp = RegularGridInterpolator((grid, grid), np.where(np.isfinite(radius), radius, 0),
                                     bounds_error=False, fill_value=np.nan)
    coords = np.full((angular_bins, 2), np.nan)
    counts = np.zeros(angular_bins, dtype=int)
    residual = np.full(angular_bins, np.nan)
    for k, phi in enumerate(2 * np.pi * np.arange(angular_bins) / angular_bins):
        direction = np.array([-np.sin(phi), np.cos(phi)])
        # Keep round-off at the square boundary from causing out-of-bounds reads.
        limit = min(2 * r_s, grid[-1] / np.max(np.abs(direction))) * (1 - 1e-12)

        def func(rho):
            return float(interp((rho * direction)[None, :])[0] - r_s)

        samples = np.linspace(0, limit, scan_samples if crossing == "first-outward" else 2)
        values = interp(samples[:, None] * direction) - r_s
        exact = np.flatnonzero(values == 0)
        changes = np.flatnonzero(np.isfinite(values[:-1]) & np.isfinite(values[1:])
                                 & (np.signbit(values[:-1]) != np.signbit(values[1:]))
                                 & (values[:-1] != 0) & (values[1:] != 0))
        roots = [float(samples[j]) for j in exact]
        for j in changes:
            roots.append(brentq(func, samples[j], samples[j + 1], xtol=1e-11))
        counts[k] = len(roots)
        if roots:
            rho = min(roots)
            coords[k] = rho * direction
            residual[k] = func(rho)
    diagnostic = {"crossing_count": counts, "residual": residual,
                  "valid": np.all(np.isfinite(coords), axis=1)}
    return (coords, diagnostic) if return_diagnostics else coords


def legacy_rms_amplitude(evpa):
    """Exact cell-31 RMS convention (including its unusual 16-state centering).

    This includes higher harmonics: sqrt(sum(delta_chi**2)/8), NOT a first
    harmonic fit. For 16 states the old mean uses states[1:] but RMS uses all 16.
    DFT amplitudes and phases are separately supplied in every stage product.
    """
    values = np.asarray(evpa, dtype=float)
    if values.shape[0] not in (16, 17):
        raise ValueError("expected 16 or 17 EVPA states")
    lifted = _wrap(values[1:] - values[0]) + values[0]
    mean = np.mean(lifted, axis=0)
    mean = np.where(mean < -np.pi / 2, mean + np.pi, mean)
    mean = np.where(mean > np.pi / 2, mean - np.pi, mean)
    phases = values[1:] if len(values) == 17 else values
    return np.sqrt(np.sum(_wrap(phases - mean)**2, axis=0) / 8)


def _angular_sample(values, rotation_angle, interpolation):
    """Rotate angular labels from +y CCW, with the active spin-2 EVPA shift."""
    phi = np.arange(values.shape[-1]) * 2 * np.pi / values.shape[-1]
    target = (phi - rotation_angle) % (2 * np.pi)
    if interpolation == "legacy-angle":
        # Deliberately nonperiodic and scalar: exact ImageS method-1 behavior.
        result = np.interp(target, phi, values)
    elif interpolation == "spin2-periodic":
        z = np.exp(2j * values)
        vector = (np.interp(target, phi, z.real, period=2 * np.pi)
                  + 1j * np.interp(target, phi, z.imag, period=2 * np.pi))
        result = np.where(np.abs(vector) > 1e-12, np.angle(vector) / 2, np.nan)
    else:
        raise ValueError("angular_interpolation must be spin2-periodic|legacy-angle")
    return _wrap(result + rotation_angle)


def get_analytical_A_1d(axion_field_coord, r_coord, r_s, *, fov=46.2262,
                         rotation_angle=0.6 * np.pi, angular_bins=360,
                         amplitude_method="legacy-rms", crossing="first-outward",
                         periodic=False):
    """Cell-44 amplitude-map interpolation on BL-radius crossings.

    Default rotation matches the *helper* in the notebook; the production stage
    explicitly passes the same rotation to analytical and transport curves.
    """
    field = np.asarray(axion_field_coord)
    if field.ndim != 3 or field.shape[0] != 16 or field.dtype.kind not in "iuf":
        raise ValueError("axion field must be real (16,N,N), in EVPA radians")
    _square(field[0])
    if field.shape[1:] != np.shape(r_coord):
        raise ValueError("axion field and radius map grids must match")
    if not np.isfinite(rotation_angle):
        raise ValueError("rotation_angle must be finite")
    if amplitude_method not in ("legacy-rms", "dft"):
        raise ValueError("amplitude_method must be legacy-rms|dft")
    amplitude = (legacy_rms_amplitude(field) if amplitude_method == "legacy-rms"
                 else fit_phase_dft(field)[..., 0])
    grid = _grid(field.shape[1], fov)
    coords = get_360_x_y(r_coord, r_s, fov=fov, angular_bins=angular_bins, crossing=crossing)
    sampled = RegularGridInterpolator((grid, grid), amplitude, bounds_error=False,
                                      fill_value=np.nan)(coords)
    phi = np.arange(angular_bins) * 2 * np.pi / angular_bins
    # Preserve cell 44 by default; production can explicitly close the seam.
    return np.interp((phi - rotation_angle) % (2 * np.pi), phi, sampled,
                     **({"period": 2 * np.pi} if periodic else {}))


def mask_and_blur_frame(stokes, image_mask, *, sigma, fov, mask_order):
    """Apply the identical boolean source mask to I,Q,U,V, before or after blur."""
    frame = _square(stokes, channels=True)
    mask = np.asarray(image_mask)
    if mask.shape != frame.shape[:2] or mask.dtype.kind != "b":
        raise ValueError("image_mask must be boolean (N,N), True=included")
    if mask_order not in ORDERS:
        raise ValueError(f"mask_order must be one of {ORDERS}")
    _grid(frame.shape[0], fov)
    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("sigma must be finite and nonnegative")
    source = np.asarray(frame, dtype=float)
    if mask_order == "mask-before-blur":
        source = np.where(mask[..., None], source, 0)
    # A zero sigma still copies; no input mutation. Never blur across channels.
    output = gaussian_filter(source, (sigma * len(frame) / fov,) * 2 + (0,), mode="nearest")
    if mask_order == "blur-before-mask":
        output = np.where(mask[..., None], output, 0)
    return output


def get_raptor_evpa(grmhddata, *, image_mask=None, oneD_coord, sigma=0,
                    init_fov=46.2262, rotation_angle=0, mask_order="mask-before-blur",
                    angular_interpolation="spin2-periodic"):
    """Port method-1: bilinear native-grid Q/U sampling, then angular rotation.

    Returns (16|17, angular_bins) EVPA radians. Like old method-1 this bypasses
    final_shape/final_fov regridding, but unlike ImageS the mask affects all IQUV.
    No aperture integration, no extra -90-bin roll (coordinates already start +y).
    """
    data = _cube(grmhddata)
    n = data.shape[1]
    grid = _grid(n, init_fov)
    coords = np.asarray(oneD_coord, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) < 4:
        raise ValueError("oneD_coord must be (angular_bins>=4,2)")
    if not np.isfinite(rotation_angle):
        raise ValueError("rotation_angle must be finite")
    mask = np.ones((n, n), dtype=bool) if image_mask is None else image_mask
    result = np.empty((len(data), len(coords)))
    for state, frame in enumerate(data):
        processed = mask_and_blur_frame(frame, mask, sigma=sigma, fov=init_fov,
                                        mask_order=mask_order)
        sampled = RegularGridInterpolator((grid, grid), processed, bounds_error=False,
                                          fill_value=np.nan)(coords)
        q, u = sampled[:, 1], sampled[:, 2]
        evpa = np.where((q == 0) & (u == 0), np.nan, np.arctan2(u, q) / 2)
        result[state] = _angular_sample(evpa, rotation_angle, angular_interpolation)
    return result


def analytical_stokes(grmhddata, axion_field_coord):
    """Cell-36 screen model: baseline P,chi; add field[k] to chi, copy state I.

    Field values are treated as EVPA rotations (radians), not multiplied by an
    inferred coupling. Field index k corresponds to psi=2*pi*k/16 for DFT.
    Like the old model V=0. A separate baseline frame produces 17 states.
    """
    data = np.asarray(grmhddata)
    if data.ndim == 3:
        baseline = _square(data, channels=True)
    else:
        data = _cube(data)
        if len(data) != 17:
            raise ValueError("analytical screen synthesis needs explicit noaxion state0 or frame")
        baseline = data[0]
    field = np.asarray(axion_field_coord)
    if field.shape != (16,) + baseline.shape[:2] or field.dtype.kind not in "iuf":
        raise ValueError("axion field must be real (16,N,N) matching baseline")
    chi = 0.5 * np.arctan2(baseline[..., 2], baseline[..., 1])
    p = np.hypot(baseline[..., 1], baseline[..., 2])
    output = np.zeros((17,) + baseline.shape, dtype=float)
    output[..., 0] = baseline[..., 0] if data.ndim == 3 else data[..., 0]
    angles = np.concatenate([chi[None], chi[None] + field], axis=0)
    output[..., 1] = p * np.cos(2 * angles)
    output[..., 2] = p * np.sin(2 * angles)
    return output


def get_analytical_evpa(grmhddata, axion_field_coord, **kwargs):
    """Synthesize cell-36 analytical Stokes and use the same corrected sampler."""
    return get_raptor_evpa(analytical_stokes(grmhddata, axion_field_coord), **kwargs)


def _load(path):
    # LFS pointer error is actionable, without invoking git or materializing data.
    with path.open("rb") as stream:
        if stream.read(6) != b"\x93NUMPY":
            raise ValueError(f"not a materialized NPY (possibly Git LFS pointer): {path}")
    return np.load(path, mmap_mode="r", allow_pickle=False)


def _ratio(numerator, denominator, floor):
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (np.abs(denominator) > floor)
    return np.divide(numerator, denominator, out=np.full(np.broadcast_shapes(
        np.shape(numerator), np.shape(denominator)), np.nan), where=valid)


def run_transport_diagnostics(
    input_dir, output_dir, *, radii=(4, 5, 6, 7, 8, 9, 10), r_width=0.5,
    fov=46.2262, sigma=2.16871, rotation_angle=0, angular_bins=360,
    mask_order="mask-before-blur", amplitude_method="legacy-rms",
    angular_interpolation="spin2-periodic", crossing="first-outward",
    spin=0.9375, denominator_floor=1e-12,
) -> dict[str, Path]:
    """Read four snapshot variants, write figures + transport.npz + transport.json.

    Always compare both mask orders and sigma=(0,sigma). mask_order only selects
    the plotted primary; other choices remain in NPZ. r=4 and r=9 are required
    for the Figure-4-style source map. No notebook execution or parsing occurs.
    """
    radii = np.asarray(radii, dtype=float)
    if (radii.ndim != 1 or not len(radii) or not np.all(np.isfinite(radii))
            or np.any(radii <= 0) or len(np.unique(radii)) != len(radii)
            or not all(np.any(radii == r) for r in (4, 9))):
        raise ValueError("radii must be unique positive values including 4 and 9")
    if mask_order not in ORDERS or amplitude_method not in ("legacy-rms", "dft"):
        raise ValueError("invalid mask_order or amplitude_method")
    if angular_interpolation not in ("spin2-periodic", "legacy-angle"):
        raise ValueError("invalid angular_interpolation")
    for value, name in ((r_width, "r_width"), (sigma, "sigma"),
                        (denominator_floor, "denominator_floor")):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not np.isfinite(rotation_angle) or not np.isfinite(spin) or abs(spin) > 1:
        raise ValueError("rotation must be finite; spin must be in [-1,1]")
    root, destination = Path(input_dir), Path(output_dir)
    names = (*VARIANTS, "common_noaxion", "r_coord_400", "axion_field_coord_400_16")
    sources = {name: root / f"{name}.npy" for name in names}
    inputs = {name: _load(path) for name, path in sources.items()}
    radius = _square(inputs["r_coord_400"])
    n = len(radius)
    _grid(n, fov)
    for name in VARIANTS:
        data = _cube(inputs[name])
        if data.shape[1:3] != radius.shape:
            raise ValueError(f"{name} Stokes and radius grids must match")
    if inputs["common"].shape[0] != 17:
        raise ValueError("common must contain noaxion state0 + 16 phases")
    field = inputs["axion_field_coord_400_16"]
    baseline = inputs["common_noaxion"]
    if baseline.ndim == 4 and baseline.shape[0] == 1:
        baseline = baseline[0]
    _square(baseline, channels=True)
    if baseline.shape[:2] != radius.shape:
        raise ValueError("common_noaxion and radius grids must match")
    analytic = analytical_stokes(inputs["common"], field)
    curves = []
    diagnostics = []
    for r_s in radii:
        coords, diagnostic = get_360_x_y(radius, r_s, fov=fov, angular_bins=angular_bins,
                                         crossing=crossing, return_diagnostics=True)
        curves.append(coords)
        diagnostics.append(diagnostic)
    curves = np.asarray(curves)
    photon_r = float(2 * (1 + np.cos(2 / 3 * np.arccos(spin))))
    photon_coords, photon_diagnostic = get_360_x_y(
        radius, photon_r, fov=fov, angular_bins=angular_bins, crossing=crossing,
        return_diagnostics=True)
    paths = {"diagnostics": destination / "transport.npz", "metadata": destination / "transport.json"}
    paths.update({name: destination / "figures" / f"{name}.png"
                  for name in ("amplitudes", "analytic_ratios", "mask_order_comparison", "source_rings")})
    for path in paths.values():
        if path.exists():
            raise FileExistsError(f"refusing to overwrite stage output: {path}")
        if path.resolve() in {p.resolve() for p in sources.values()}:
            raise ValueError(f"output aliases input: {path}")
    # (variant, radius, mask_order, blur, angle[, fit_column]); analytical is last.
    shape = (len(VARIANTS) + 1, len(radii), 2, 2, int(angular_bins))
    rms = np.full(shape, np.nan)
    fits = np.full(shape + (3,), np.nan)
    native_amplitude = np.full((2, len(radii), int(angular_bins)), np.nan)
    for v, data in enumerate([*(inputs[name] for name in VARIANTS), analytic]):
        for r, r_s in enumerate(radii):
            mask = np.isfinite(radius) & (np.abs(radius - r_s) <= r_width / 2)
            for o, order in enumerate(ORDERS):
                for b, blur in enumerate((0, sigma)):
                    evpa = get_raptor_evpa(data, image_mask=mask, oneD_coord=curves[r],
                                           sigma=blur, init_fov=fov, rotation_angle=rotation_angle,
                                           mask_order=order, angular_interpolation=angular_interpolation)
                    rms[v, r, o, b] = legacy_rms_amplitude(evpa)
                    fits[v, r, o, b] = fit_phase_dft(evpa)
                    if v < 2 and b == 0 and order == mask_order:
                        native_evpa = evpa if rotation_angle == 0 else get_raptor_evpa(
                            data, image_mask=mask, oneD_coord=curves[r], sigma=0, init_fov=fov,
                            rotation_angle=0, mask_order=order,
                            angular_interpolation=angular_interpolation)
                        native_amplitude[v, r] = (legacy_rms_amplitude(native_evpa)
                                                  if amplitude_method == "legacy-rms"
                                                  else fit_phase_dft(native_evpa)[..., 0])
    analytical_rms = np.array([get_analytical_A_1d(
        field, radius, r_s, fov=fov, rotation_angle=rotation_angle,
        angular_bins=angular_bins, amplitude_method="legacy-rms", crossing=crossing,
        periodic=angular_interpolation == "spin2-periodic") for r_s in radii])
    analytical_dft = np.array([get_analytical_A_1d(
        field, radius, r_s, fov=fov, rotation_angle=rotation_angle,
        angular_bins=angular_bins, amplitude_method="dft", crossing=crossing,
        periodic=angular_interpolation == "spin2-periodic") for r_s in radii])
    amplitude = rms if amplitude_method == "legacy-rms" else fits[..., 0]
    denominator = analytical_rms if amplitude_method == "legacy-rms" else analytical_dft
    ratios = _ratio(amplitude, denominator[None, :, None, None, :], denominator_floor)
    transport_ratio = _ratio(amplitude[1], amplitude[0], denominator_floor)
    payload = dict(
        variants=np.array([*VARIANTS, "analytical_screen_common_state0"]), radii=radii,
        mask_orders=np.array(ORDERS), blur_sigmas=np.array([0, sigma]),
        varphi_rad=np.arange(angular_bins) * 2 * np.pi / angular_bins,
        ring_xy=curves, crossing_count=np.array([d["crossing_count"] for d in diagnostics]),
        crossing_residual=np.array([d["residual"] for d in diagnostics]),
        ring_valid=np.array([d["valid"] for d in diagnostics]),
        photon_reference_bl_radius=np.array(photon_r), photon_reference_xy=photon_coords,
        photon_crossing_count=photon_diagnostic["crossing_count"],
        photon_crossing_residual=photon_diagnostic["residual"],
        amplitude_rms=rms, fits_dft=fits, amplitude=amplitude, analytic_ratio=ratios,
        analytical_map_rms=analytical_rms, analytical_map_dft=analytical_dft,
        jonly_common_ratio=transport_ratio,
        source_color_ratio=_ratio(native_amplitude[1], native_amplitude[0], denominator_floor),
        mask_order_amplitude_delta=amplitude[:, :, 0] - amplitude[:, :, 1],
        noaxion_intensity=np.asarray(baseline[..., 0]),
        common_state0_minus_noaxion=np.asarray(inputs["common"][0], dtype=float) - baseline,
    )
    destination.mkdir(parents=True, exist_ok=True)
    paths["amplitudes"].parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(paths["diagnostics"], **payload)
    _figures(paths, payload, primary_order=ORDERS.index(mask_order), fov=fov)
    metadata = dict(
        stage="transport-diagnostics", production_runtime="Python/NPY only; no notebook access",
        source_port={"path": "GRMHD/main.ipynb", "cells_zero_based": [7, 15, 31, 36, *range(43, 65)]},
        inputs={name: {"path": str(sources[name]), "shape": list(a.shape), "dtype": str(a.dtype),
                       "bytes": sources[name].stat().st_size} for name, a in inputs.items()},
        variants=list(VARIANTS), variant_semantics_status="inferred, not verified RAPTOR configuration",
        variant_semantics={"common": "reference transport", "jonly": "emissivity-only?",
                           "eq_jonly": "equatorial emissivity-only?", "fake": "synthetic/screen control?"},
        method_choice={"mask_order_primary": mask_order, "compared_mask_orders": list(ORDERS),
                       "mask_channels": list("IQUV"), "mask_definition": "finite BL r within r_width/2",
                       "amplitude_primary": amplitude_method, "angular_interpolation": angular_interpolation,
                       "analytical_amplitude_interpolation": ("periodic scalar" if angular_interpolation
                                                              == "spin2-periodic" else "legacy scalar"),
                       "radius_crossing": crossing, "radius_nan_geometry": "nonfinite->0 as old code",
                       "sampler": "native-grid bilinear IQUV, equivalent to method-1 before angular remap",
                       "final_shape_bypassed": True, "boundary_blur": "nearest", "sampling_outside": "NaN"},
        phase_convention="DFT mean+A*cos(2*pi*k/16+phi); positive amplitude; state0 excluded",
        phase_origin_status="array-index convention verified synthetically; producer phase origin unverified",
        amplitude_units="EVPA radians", field_interpretation="add field[k] directly to EVPA, cell 36",
        rms_caveat="legacy RMS includes harmonics; 16-state helper centers using states[1:]",
        angular_convention="varphi=0 at +y, CCW: (x,y)=(-rho*sin(varphi),rho*cos(varphi)); no -90 roll",
        rotation_angle_rad=float(rotation_angle), analytic_and_transport_rotation_matched=True,
        notebook_alignment_caveat=("cell44 analytical default 0.6*pi but cells47/49 transport use 0; "
                                   "corrected here"),
        fov_rg=float(fov), sigma_rg_stddev=float(sigma), r_width_bl_rg=float(r_width),
        sigma_origin="ImageS None -> 2.16871 native FOV units; not FWHM, not a verified physical beam",
        radii_bl_rg=radii.tolist(), angular_bins=int(angular_bins),
        denominator_floor_rad=float(denominator_floor),
        finite_ratio_counts=np.sum(np.isfinite(ratios), axis=-1).tolist(),
        missing_crossings=np.sum(~payload["ring_valid"], axis=-1).tolist(),
        multiple_crossings=np.sum(payload["crossing_count"] > 1, axis=-1).tolist(),
        spin=float(spin), photon_reference_bl_radius=photon_r,
        photon_reference_caveat=("old Kerr equatorial formula mapped via r_coord; "
                                 "not verified observed photon ring"),
        source_map_colors=("unblurred jonly/common; BL9 on BL9; "
                           "BL4 ratio on photon-reference curve (legacy proxy)"),
        noaxion_usage=("background and saved common_state0-minus-noaxion; "
                       "analytical screen uses common state0"),
        validation_status="synthetic-tested port, not verified producer switches or physical configuration",
        products={name: str(path) for name, path in paths.items()},
    )
    paths["metadata"].write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return paths


def _figures(paths, data, *, primary_order, fov):
    # Local, headless, no TeX/font dependency and no global pyplot rc mutation.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    degrees = np.rad2deg(data["varphi_rad"])
    radii = data["radii"]
    colors = ("tab:blue", "tab:green", "black", "tab:orange", "tab:purple")
    labels = data["variants"]

    def canvas(rows, cols, size):
        figure = Figure(figsize=size, layout="constrained")
        FigureCanvasAgg(figure)
        return figure, figure.subplots(rows, cols, squeeze=False)

    for name, key, ylabel in (("amplitudes", "amplitude", "A (EVPA rad)"),
                              ("analytic_ratios", "analytic_ratio", "A / analytical map A")):
        fig, axes = canvas(1, 2, (12, 4))
        for ax, r_s in zip(axes.flat, (9, 4)):
            r = int(np.flatnonzero(radii == r_s)[0])
            for v, label in enumerate(labels):
                for b in (0, 1):
                    ax.plot(degrees, data[key][v, r, primary_order, b], color=colors[v],
                            alpha=0.3 if b == 0 else 1, linestyle="--" if b == 0 else "-",
                            label=str(label) if b else None)
            if name == "analytic_ratios":
                ax.axhline(1, color="grey", lw=0.8)
            ax.set(title=f"Source BL r={r_s}; dashed: unblurred, solid: blurred",
                   xlabel="varphi (degrees, +y CCW)", ylabel=ylabel, xlim=(0, 360))
            ax.set_xticks([0, 90, 180, 270, 360])
        axes[0, 0].legend(fontsize=8)
        fig.savefig(paths[name], dpi=160)
    fig, axes = canvas(2, 2, (12, 7))
    for v, ax in enumerate(axes.flat):
        for r_s, color in ((9, "tab:blue"), (4, "tab:red")):
            r = int(np.flatnonzero(radii == r_s)[0])
            for o, order in enumerate(ORDERS):
                ax.plot(degrees, data["analytic_ratio"][v, r, o, 1], color=color,
                        linestyle="-" if o == 0 else "--", label=f"r={r_s}: {order}")
        ax.set(title=str(labels[v]), xlabel="varphi (degrees)", ylabel="blurred A / analytical map A")
        ax.legend(fontsize=7)
    fig.savefig(paths["mask_order_comparison"], dpi=160)
    fig, axes = canvas(1, 1, (8, 7))
    ax = axes[0, 0]
    grid = _grid(data["noaxion_intensity"].shape[0], fov)
    intensity = data["noaxion_intensity"]
    positive = np.isfinite(intensity) & (intensity > 0)
    maximum = float(np.max(intensity[positive])) if np.any(positive) else 1.0
    image = np.full(intensity.shape, np.nan)
    image[positive] = np.log10(intensity[positive] / maximum)
    background = ax.pcolormesh(grid, grid, image.T, cmap="gray", vmin=-4, vmax=0, shading="auto")
    fig.colorbar(background, ax=ax, label="log10(I / Imax)", shrink=0.8)
    for r_s, coords, label in (
        (9, data["ring_xy"][int(np.flatnonzero(radii == 9)[0])], "source BL r=9"),
        (4, data["photon_reference_xy"], "photon-orbit reference; color uses BL r=4 proxy"),
    ):
        r = int(np.flatnonzero(radii == r_s)[0])
        points = np.concatenate([coords, coords[:1]])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        # Source geometry is native: use the separately sampled rotation=0 colors.
        values = data["source_color_ratio"][r]
        valid = np.all(np.isfinite(segments), axis=(1, 2)) & np.isfinite(values)
        line = LineCollection(segments[valid], cmap="RdBu_r", norm=Normalize(0.9, 1.1), linewidth=3)
        line.set_array(values[valid])
        ax.add_collection(line)
        if np.any(valid):
            midpoint = coords[len(coords) // 8]
            if np.all(np.isfinite(midpoint)):
                ax.annotate(label, midpoint, xytext=(5, 5), textcoords="offset points", fontsize=7,
                            color="tab:orange", bbox={"facecolor": "black", "alpha": 0.6})
    # A separate scalar mappable also works when no radius crossings are present.
    from matplotlib.cm import ScalarMappable
    fig.colorbar(ScalarMappable(norm=Normalize(0.9, 1.1), cmap="RdBu_r"), ax=ax,
                 orientation="horizontal", label="unblurred A(jonly) / A(common)")
    ax.set(xlabel="x / rg", ylabel="y / rg", title="Native source coordinates; inferred variant labels",
           xlim=(-16.3434, 16.3434), ylim=(-16.3434, 16.3434), aspect="equal")
    fig.savefig(paths["source_rings"], dpi=160)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="directory containing the seven materialized .npy files")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--radii", nargs="+", type=float, default=[4, 5, 6, 7, 8, 9, 10])
    parser.add_argument("--r-width", type=float, default=0.5)
    parser.add_argument("--fov", type=float, default=46.2262)
    parser.add_argument("--sigma", type=float, default=2.16871, help="Gaussian stddev in native FOV rg units")
    parser.add_argument("--rotation-angle", type=float, default=0, help="radians; same for all curves")
    parser.add_argument("--angular-bins", type=int, default=360)
    parser.add_argument("--mask-order", choices=ORDERS, default=ORDERS[0])
    parser.add_argument("--amplitude-method", choices=("legacy-rms", "dft"), default="legacy-rms")
    parser.add_argument("--angular-interpolation", choices=("spin2-periodic", "legacy-angle"),
                        default="spin2-periodic")
    parser.add_argument("--crossing", choices=("first-outward", "legacy-bracket"), default="first-outward")
    parser.add_argument("--spin", type=float, default=0.9375)
    parser.add_argument("--denominator-floor", type=float, default=1e-12)
    try:
        products = run_transport_diagnostics(**vars(parser.parse_args(argv)))
    except (ValueError, OSError, TypeError) as error:
        parser.error(str(error))
    print(json.dumps({name: str(path) for name, path in products.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
