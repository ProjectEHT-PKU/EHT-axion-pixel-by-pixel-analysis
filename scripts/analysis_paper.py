#!/usr/bin/env python3
"""Analysis: everything the paper uses from the pipeline, in one directory (stage analysis:paper).

Reads finished products only (limits, the three checks, templates, the EHT ensemble,
the delivered axion-off background cubes, the RAPTOR plan and build) and writes:

  figures/*.pdf      vector figures
  numbers.tex        one \\newcommand per number quoted in the paper, with its source
  numbers.json       the same values, unrounded numbers and sources
  table-limits.tex   complete tabular environment of the observed limits
  table-chi2.tex     complete tabular environment of the chi^2 values at PAPER_ALPHA

Nothing is refitted: limits, fits and covariances are copied from the saved products.
The constants below are presentation choices of the paper (which alpha and snapshot are
shown, where radial profiles are sampled); they are recorded in the manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re

import _bootstrap  # noqa: F401
import numpy as np

from _data_root import add_data_root_option, fill_defaults
from _environment import resolve_environment
from eht_axion import provenance
from eht_axion.analysis.likelihood import MODES
from eht_axion.analysis.model import cloud_amplitude_scaling
from eht_axion.analysis.templates import RunLayout, _integrate
from eht_axion.eht.ensemble import global_diagnostics
from eht_axion.experiment import DEFAULT_PATH, EHT_DAYS, physical_scale, read_experiment
from eht_axion.physics import HBAR_EV_S, omega
from eht_axion.raptor.layout import input_raptor_files

PAPER_ALPHA = 0.30  # alpha of the response maps, templates figure, chi^2 table and fit figure
PAPER_SNAPSHOT_LABEL = 2300  # GRMHD label of the response-map figure
TICKS = dict(stride=4, intensity_fraction=0.1, polarization_fraction=0.2, length_uas=9.0)  # EVPA ticks
RADIAL_SAMPLE_RG = (5.0, 10.0, 20.0)  # radii of the radial-amplitude panel
WINDOW_SAMPLE_RG = 20.0  # radius at which W(r) is quoted
ORBIT_RADIUS_RG = 5.0  # radius of the quoted Newtonian orbital period
CLOSE_SEPARATION_UAS = 15.0  # pixel pairs closer than this: quoted median |correlation|
GROWTH_RATE_DENOMINATOR = 24.0  # omega_i = a alpha^9 / (24 t_g), leading order (Detweiler 1980)
TWO_EMITTER = dict(
    left_ratio=0.7,
    right_phi_degrees=45.0,
    check_amplitude=1e-5,
    check_phases=256,
    check_points=((0.7, 1.2, 0.8), (0.9, 1.5, 0.785), (0.4, 0.3, 2.0)),
)
DIMENSIONS = {"1d": "Sectors", "2d": "Pixels"}
MODE_WORDS = {"std": "Std", "cor": "Cor", "cor_fixed_background": "CorBg"}
DATA_WORDS = {"observed": "Obs", "zero-data": "Zero"}
YEAR_HOURS = 365.25 * 24
_ONES = (
    "",
    "One",
    "Two",
    "Three",
    "Four",
    "Five",
    "Six",
    "Seven",
    "Eight",
    "Nine",
    "Ten",
    "Eleven",
    "Twelve",
    "Thirteen",
    "Fourteen",
    "Fifteen",
    "Sixteen",
    "Seventeen",
    "Eighteen",
    "Nineteen",
)
_TENS = ("", "", "Twenty", "Thirty", "Forty", "Fifty", "Sixty", "Seventy", "Eighty", "Ninety")
_WINDOW = re.compile(r"window\s*=\s*exp\(\s*-\s*pow\(\s*r\s*/\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\)\s*\)\s*;")


# ----------------------------------------------------------------------------- formatting
def word(number):
    """1..99 as a CamelCase English word (LaTeX macro names cannot contain digits)."""
    if not 0 < number < 100:
        raise ValueError("word() covers 1..99")
    return _ONES[number] if number < 20 else _TENS[number // 10] + _ONES[number % 10]


def alpha_word(alpha):
    return "Alpha" + word(int(round(alpha * 100)))


def sig(value, digits):
    """Significant digits, keeping trailing zeros; integers above the digits are rounded."""
    value = float(value)
    if value == 0:
        return "0"
    places = digits - 1 - int(np.floor(np.log10(abs(value))))
    if places < 0:
        return f"{round(value, places):.0f}"
    return f"{value:.{places}f}"


def approx(value):
    """Two significant digits, but never rounding the integer part of values >= 100."""
    return f"{value:.0f}" if abs(value) >= 100 else sig(value, 2)


def sci(value, digits):
    """LaTeX math 'm\\times10^{e}' with the mantissa to the given significant digits."""
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = round(value / 10**exponent, digits - 1)
    if abs(mantissa) >= 10:
        exponent += 1
        mantissa /= 10
    return rf"{mantissa:.{digits - 1}f}\times10^{{{exponent}}}"


def order(value):
    return rf"10^{{{int(round(np.log10(value)))}}}"


def percent_change(dex):
    return 100 * (10 ** abs(dex) - 1)


def angle_runs(angles, spacing):
    """'$74^\\circ$--$104^\\circ$ and $144^\\circ$--$334^\\circ$' for runs of equally spaced angles."""
    angles = sorted(int(a) for a in angles)
    runs, start = [], angles[0]
    for previous, current in zip(angles, angles[1:] + [None]):
        if current is None or current - previous != spacing:
            runs.append(rf"${start}^\circ$" if start == previous else rf"${start}^\circ$--${previous}^\circ$")
            start = current
    return runs[0] if len(runs) == 1 else ", ".join(runs[:-1]) + " and " + runs[-1]


class Numbers:
    """Named numbers: TeX text, unrounded value and source of each."""

    def __init__(self):
        self.items = {}

    def add(self, name, tex, value, source, description):
        if not re.fullmatch(r"Num[A-Z][A-Za-z]*", name) or name in self.items:
            raise ValueError(f"Invalid or duplicate macro name: {name}")
        if isinstance(value, np.generic):
            value = value.item()
        self.items[name] = dict(tex=tex, value=value, source=source, description=description)

    def tex(self):
        lines = [
            "% Numbers quoted in the paper, written by scripts/analysis_paper.py.",
            "% Each macro: description; source (file in the data root and field).",
            "% Values containing \\times or ^ are for math mode.",
            "",
        ]
        for name, item in self.items.items():
            lines.append(f"% {item['description']}; source: {item['source']}")
            lines.append(f"\\newcommand{{\\{name}}}{{{item['tex']}}}")
        return "\n".join(lines) + "\n"

    def json(self):
        return dict(schema="analysis-paper-numbers-v1", numbers=self.items)


# ----------------------------------------------------------------------------- inputs
class Inputs:
    """Verified access to the files of completed stages; records every consumed file."""

    def __init__(self):
        self.consumed = []

    def stage(self, directory, stage, files):
        directory = Path(directory).resolve()
        record = provenance.verify_stage(directory, recursive=False)
        if record["stage"] != stage:
            raise ValueError(f"{directory}: expected stage {stage}, found {record['stage']}")
        paths = [directory / f for f in files]
        provenance.require_declared_outputs(directory, paths)
        self.consumed += [provenance.manifest_path(directory)] + paths
        return paths

    def file(self, path):
        path = Path(path).resolve()
        self.consumed.append(path)
        return path


def npz(path):
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def read_json(path):
    return json.loads(Path(path).read_text())


def radial_tables(args, plan, inputs):
    """Radial-table paths: --radial-* > EHT_RADIAL_* > raptor/plan.json; SHA-256 must equal the plan."""
    environment = resolve_environment(args)
    tables = {}
    for key, option, env in (
        ("low", args.radial_low, "EHT_RADIAL_LOW"),
        ("high", args.radial_high, "EHT_RADIAL_HIGH"),
    ):
        record = plan["radial_tables"][key]
        path = option or (Path(environment[env]) if environment.get(env) else Path(record["path"]))
        if provenance.sha256(path) != record["sha256"]:
            raise ValueError(f"Radial table {path} differs from the one recorded in the RAPTOR plan")
        tables[key] = (np.loadtxt(inputs.file(path), delimiter=","), record["spin"])
    return tables


def window_constants(axion_source):
    """W(r) = exp(-(r/r_w)^n) of the pinned RAPTOR core: (r_w, n)."""
    matches = _WINDOW.findall(Path(axion_source).read_text())
    if len(matches) != 1:
        raise ValueError("RAPTOR field window not found exactly once in axion.c")
    return float(matches[0][0]), float(matches[0][1])


# ----------------------------------------------------------------------------- numbers
def add_scales(numbers, config, times_tg, lag_hours, pairs):
    scale = physical_scale(config)
    tg, rg = scale.tg_hours, scale.angular_rg_uas
    numbers.add("NumTgHours", sig(tg, 3), tg, "experiment.yaml shared (physical_scale)", "t_g in hours")
    numbers.add(
        "NumRgMicroarcsec",
        sig(rg, 3),
        rg,
        "experiment.yaml shared (physical_scale)",
        "angular size of r_g in microarcseconds",
    )
    camera = config["raptor"]["camera_fov_rg"] * rg
    numbers.add(
        "NumCameraFovMicroarcsec",
        sig(camera, 3),
        camera,
        "experiment.yaml raptor.camera_fov_rg",
        "RAPTOR camera width in microarcseconds",
    )
    labels = config["raptor"]["snapshots"]
    count = len(times_tg)
    if count != (labels["stop"] - labels["start"]) // labels["step"] + 1:
        raise ValueError("Delivered snapshot count differs from raptor.snapshots")
    numbers.add(
        "NumSnapshots", str(count), count, "raptor/delivered generated.json times_tg", "GRMHD snapshots"
    )
    span = times_tg[-1] - times_tg[0]
    numbers.add(
        "NumSpanTg",
        f"{span:g}",
        span,
        "raptor/delivered generated.json times_tg",
        "time spanned by the snapshots in t_g",
    )
    numbers.add(
        "NumSpanDays",
        sig(span * tg / 24, 2),
        span * tg / 24,
        "generated.json times_tg; t_g",
        "time spanned by the snapshots in days",
    )
    numbers.add(
        "NumSimTimeStartTg",
        f"{times_tg[0]:g}",
        times_tg[0],
        "generated.json times_tg[0]",
        "simulation time of the first snapshot (t_g)",
    )
    numbers.add(
        "NumSimTimeEndTg",
        f"{times_tg[-1]:g}",
        times_tg[-1],
        "generated.json times_tg[-1]",
        "simulation time of the last snapshot (t_g)",
    )
    numbers.add(
        "NumLagHours",
        sig(lag_hours, 3),
        lag_hours,
        "analysis/limits background-selected-1d.json lag_hours",
        "background lag in hours",
    )
    numbers.add(
        "NumBackgroundPairs",
        str(pairs),
        pairs,
        "analysis/limits background-selected-1d.npz differences",
        "lagged snapshot pairs of the background",
    )
    period = 2 * np.pi * ORBIT_RADIUS_RG**1.5 * tg / 24
    numbers.add(
        "NumOrbitalPeriod" + word(int(ORBIT_RADIUS_RG)) + "RgDays",
        sig(period, 2),
        period,
        "2 pi r^{3/2} t_g",
        f"Newtonian orbital period at {ORBIT_RADIUS_RG:g} r_g in days",
    )
    gap = config["analysis"]["inference"]["day_gap_hours"]
    n = config["shared"]["principal_n"]
    alphas = config["shared"]["alphas"]
    for alpha in (alphas[0], alphas[-1]):
        frequency = omega(alpha, n) / tg  # rad per hour
        mass = alpha * HBAR_EV_S / scale.tg_seconds
        numbers.add(
            f"NumAxionMass{alpha_word(alpha)}EV",
            sci(mass, 2),
            mass,
            "alpha hbar / t_g",
            f"axion mass at alpha={alpha:.2f} in eV (math mode)",
        )
        days = 2 * np.pi / frequency / 24
        numbers.add(
            f"NumPeriodDays{alpha_word(alpha)}",
            sig(days, 2),
            days,
            "2 pi / omega(alpha)",
            f"oscillation period at alpha={alpha:.2f} in days",
        )
        factor = 2 * np.sin(frequency * gap / 2)
        numbers.add(
            f"NumDayDiffFactor{alpha_word(alpha)}",
            sig(factor, 2),
            factor,
            "2 sin(omega day_gap_hours / 2)",
            f"one-day difference factor at alpha={alpha:.2f}",
        )
    alpha, spin = alphas[0], config["raptor"]["fluid_spin"]
    growth = GROWTH_RATE_DENOMINATOR * tg / (spin * alpha**9) / YEAR_HOURS
    numbers.add(
        f"NumGrowthTimeYears{alpha_word(alpha)}",
        sci(growth, 1),
        growth,
        "24 t_g / (a alpha^9), a = raptor.fluid_spin",
        f"field e-folding time at alpha={alpha:.2f} in years",
    )
    cloud = config["analysis"]["cloud"]
    for alpha in (alphas[0], PAPER_ALPHA, alphas[-1]):
        b = cloud_amplitude_scaling(
            alpha,
            coefficient=cloud["coefficient"],
            pivot_alpha=cloud["pivot_alpha"],
            exponent=cloud["exponent"],
        )
        numbers.add(
            f"NumCloudB{alpha_word(alpha)}",
            f"{b:.2f}",
            b,
            "experiment.yaml analysis.cloud",
            f"cloud amplitude b at alpha={alpha:.2f}",
        )
    numbers.add(
        "NumPaperAlpha",
        f"{PAPER_ALPHA:.2f}",
        PAPER_ALPHA,
        "analysis_paper.py PAPER_ALPHA",
        "alpha of the response maps, chi^2 table and fit figure",
    )
    numbers.add(
        "NumPaperSnapshot",
        str(PAPER_SNAPSHOT_LABEL),
        PAPER_SNAPSHOT_LABEL,
        "analysis_paper.py PAPER_SNAPSHOT_LABEL",
        "GRMHD label of the response maps",
    )


def radial_profiles(config, tables):
    """[(alpha, spin, high table?, r, R/Rmax)] exactly as assigned by raptor.alpha_table."""
    profiles = []
    for alpha in config["shared"]["alphas"]:
        key = config["raptor"]["alpha_table"][f"{alpha:.2f}"]
        table, spin = tables[key]
        curve = table[np.isclose(table[:, 0], alpha, atol=1e-8, rtol=0)]
        if len(curve) < 2:
            raise ValueError(f"No radial profile for alpha={alpha}")
        curve = curve[np.argsort(curve[:, 1])]
        profiles.append((alpha, spin, key == "high", curve[:, 1], curve[:, 2]))
    return profiles


def add_radial(numbers, profiles, window):
    for alpha, _spin, _high, r, response in profiles:
        peak = r[np.argmax(response)]
        numbers.add(
            f"NumRadialPeakRg{alpha_word(alpha)}",
            f"{peak:.0f}",
            peak,
            "radial table of raptor.alpha_table (raptor/plan.json)",
            f"radius of the maximum of R_211 at alpha={alpha:.2f} (r_g)",
        )
    radius, exponent = window
    numbers.add(
        "NumWindowRadiusRg",
        f"{radius:g}",
        radius,
        "raptor/runs/main/build/source/src/axion.c",
        "field window radius r_w of W(r)=exp[-(r/r_w)^n] (r_g)",
    )
    numbers.add(
        "NumWindowExponent",
        f"{exponent:g}",
        exponent,
        "raptor/runs/main/build/source/src/axion.c",
        "field window exponent n",
    )
    value = np.exp(-((WINDOW_SAMPLE_RG / radius) ** exponent))
    numbers.add(
        "NumWindowAt" + word(int(WINDOW_SAMPLE_RG)) + "Rg",
        f"{value:.2f}",
        value,
        "W(r) of axion.c",
        f"W at {WINDOW_SAMPLE_RG:g} r_g",
    )


def add_signal(numbers, amplitude_deg, sector_mean_deg, sector_fraction):
    numbers.add(
        "NumResponseMinDeg",
        sig(amplitude_deg.min(), 2),
        amplitude_deg.min(),
        "analysis/templates response/alpha-*/response-2d.npz amplitude",
        "minimum response amplitude per unit c of the paper snapshot (deg)",
    )
    numbers.add(
        "NumResponseMaxDeg",
        sig(amplitude_deg.max(), 2),
        amplitude_deg.max(),
        "analysis/templates response/alpha-*/response-2d.npz amplitude",
        "maximum response amplitude per unit c of the paper snapshot (deg)",
    )
    for label, f in (("Min", np.min), ("Max", np.max)):
        value = f(sector_mean_deg)
        numbers.add(
            f"NumSectorTemplate{label}Deg",
            sig(value, 2),
            value,
            "analysis/templates templates/1d/templates.npz amplitude_mean",
            f"{label.lower()} sector template amplitude per unit c (deg)",
        )
    for label, f in (("Min", np.min), ("Median", np.median), ("Max", np.max)):
        value = 100 * f(sector_fraction)
        numbers.add(
            f"NumSectorTemplateScatter{label}Percent",
            f"{value:.0f}",
            value,
            "templates/1d/templates.npz amplitude_std / amplitude_mean",
            f"{label.lower()} fractional snapshot scatter of the sector template (percent)",
        )


def add_model_vs_reconstruction(numbers, flux, recon_flux, model_m, model_beta, recon_m, recon_beta):
    for label, value in (("", flux.mean()), ("Min", flux.min()), ("Max", flux.max())):
        numbers.add(
            f"NumModelFlux{label}Jy",
            f"{value:.2f}",
            value,
            "raptor/delivered res-400 background stokes (sum I dOmega, state 0)",
            f"{'mean' if not label else label.lower()} 230 GHz flux of the simulated images (Jy)",
        )
    numbers.add(
        "NumReconFluxJy",
        f"{recon_flux:.2f}",
        recon_flux,
        "eht/ensemble day-*/global.npz total_I",
        "mean total flux of the reconstructions (Jy)",
    )
    excess = 100 * (flux.mean() / recon_flux - 1)
    numbers.add(
        "NumModelFluxExcessPercent",
        f"{excess:.0f}",
        excess,
        "model flux / reconstruction flux",
        "excess of the mean simulated flux over the reconstructions (percent)",
    )
    low, high = np.percentile(model_m, [5, 95])
    for name, value, text in (
        ("NumModelPolPercent", model_m.mean(), "mean"),
        ("NumModelPolLowPercent", low, "5th percentile"),
        ("NumModelPolHighPercent", high, "95th percentile"),
    ):
        numbers.add(
            name,
            f"{100 * value:.1f}",
            100 * value,
            "analysis/templates background/stokes.npy",
            f"{text} resolved linear polarization fraction of the simulated images (percent)",
        )
    numbers.add(
        "NumModelPolApproxPercent",
        f"{100 * model_m.mean():.0f}",
        100 * model_m.mean(),
        "analysis/templates background/stokes.npy",
        "mean simulated polarization fraction, rounded",
    )
    for name, value in (("NumReconPolMinPercent", min(recon_m)), ("NumReconPolMaxPercent", max(recon_m))):
        numbers.add(
            name,
            f"{100 * value:.1f}",
            100 * value,
            "eht/ensemble day-*/global.npz m_resolved",
            "per-day mean resolved polarization fraction of the reconstructions, extreme day (percent)",
        )
    ratio = model_m.mean() / np.mean(recon_m)
    numbers.add(
        "NumModelToReconPolRatio",
        f"{ratio:.0f}",
        ratio,
        "model / mean of the four days",
        "ratio of simulated to reconstructed polarization fraction, rounded",
    )
    numbers.add(
        "NumModelBetaTwo",
        sig(model_beta.mean(), 2),
        model_beta.mean(),
        "analysis/templates background/stokes.npy (global_diagnostics)",
        "mean |beta_2| of the simulated images",
    )
    for name, value in (("NumReconBetaTwoMin", min(recon_beta)), ("NumReconBetaTwoMax", max(recon_beta))):
        numbers.add(
            name,
            sig(value, 2),
            value,
            "eht/ensemble day-*/global.npz beta2_real, beta2_imag",
            "per-day mean |beta_2| of the reconstructions, extreme day",
        )


def dterm_statistics(observations):
    covariance = observations["pair_covariances"]
    sigma = np.sqrt(np.diagonal(covariance, axis1=-2, axis2=-1))
    correlation = covariance / (sigma[:, :, None] * sigma[:, None, :])
    return sigma, correlation


def add_observations(numbers, obs, info, background, config, branch_minimum):
    spacing = config["analysis"]["features"]["sector_spacing_degrees"]
    for dim, word_dim in DIMENSIONS.items():
        n = len(obs[dim]["selected_indices"])
        numbers.add(
            f"Num{word_dim}",
            str(n),
            n,
            f"analysis/limits observations-{dim}.npz selected_indices",
            f"selected {word_dim.lower()}",
        )
        numbers.add(
            f"Num{word_dim[:-1]}Components",
            str(2 * n),
            2 * n,
            f"observations-{dim}.npz",
            f"length of the {dim} data vector (two pairs)",
        )
        sigma, correlation = dterm_statistics(obs[dim])
        degrees = np.rad2deg(sigma)
        for label, value in (("Min", degrees.min()), ("Median", np.median(degrees)), ("Max", degrees.max())):
            numbers.add(
                f"NumDtermSigma{word_dim}{label}Deg",
                f"{value:.1f}",
                value,
                f"observations-{dim}.npz pair_covariances",
                f"{label.lower()} marginal D-term sigma of the {dim} day differences (deg)",
            )
        iu = np.triu_indices(n, 1)
        medians = [np.median(np.abs(c[iu])) for c in correlation]
        for label, value in (("Min", min(medians)), ("Max", max(medians))):
            numbers.add(
                f"NumDtermCorr{word_dim}{label}",
                f"{value:.2f}",
                value,
                f"observations-{dim}.npz pair_covariances",
                f"median |correlation| between {dim} features, {label.lower()} over the two pairs",
            )
        if dim == "2d":
            coordinates = obs[dim]["coordinates"]
            separation = np.linalg.norm(coordinates[:, None] - coordinates[None], axis=-1)[iu]
            if not np.any(separation < CLOSE_SEPARATION_UAS):
                raise ValueError(f"No selected pixel pair closer than {CLOSE_SEPARATION_UAS:g} uas")
            close = [np.median(np.abs(c[iu])[separation < CLOSE_SEPARATION_UAS]) for c in correlation]
            for label, value in (("Min", min(close)), ("Max", max(close))):
                numbers.add(
                    f"NumDtermCorrClosePixels{label}",
                    f"{value:.1f}",
                    value,
                    "observations-2d.npz pair_covariances, coordinates",
                    f"median |correlation| of pixels closer than {CLOSE_SEPARATION_UAS:g} uas, "
                    f"{label.lower()} over pairs",
                )
        conditions = np.asarray(info[dim]["pair_condition_numbers"], dtype=float)
        if dim == "1d":
            for k, value in enumerate(conditions):
                numbers.add(
                    f"NumConditionSectorsPair{word(k + 1)}",
                    sci(value, 2),
                    value,
                    "observations-1d.json pair_condition_numbers",
                    f"condition number of the sector pair-{k + 1} covariance (math mode)",
                )
        else:
            value = float(np.exp(np.mean(np.log(conditions))))
            numbers.add(
                "NumConditionPixelsOrder",
                order(value),
                value,
                "observations-2d.json pair_condition_numbers",
                "order of magnitude of the pixel pair condition numbers (math mode)",
            )
        largest = np.rad2deg(np.abs(obs[dim]["pair_differences"]).max())
        numbers.add(
            f"NumMaxAbsDiff{word_dim}Deg",
            f"{largest:.0f}",
            largest,
            f"observations-{dim}.npz pair_differences",
            f"largest absolute observed {dim} day difference (deg)",
        )
        sigma_bg, mean_bg = (float(np.rad2deg(background[dim][k])) for k in ("sigma_radians", "mean_radians"))
        numbers.add(
            f"NumSigmaBg{word_dim}Deg",
            f"{sigma_bg:.2f}",
            sigma_bg,
            f"background-selected-{dim}.npz sigma_radians",
            f"GRMHD background sigma, {dim} (deg)",
        )
        numbers.add(
            f"NumMeanBg{word_dim}Deg",
            f"{mean_bg:.2f}",
            mean_bg,
            f"background-selected-{dim}.npz mean_radians",
            f"GRMHD background mean, {dim} (deg)",
        )
    numbers.add(
        "NumSectorAngles",
        angle_runs(obs["1d"]["selected_indices"], spacing),
        [int(a) for a in obs["1d"]["selected_indices"]],
        "observations-1d.npz selected_indices",
        "position angles of the selected sectors (text with inline math)",
    )
    numbers.add(
        "NumMinBranchDistanceSigma",
        f"{branch_minimum:.1f}",
        branch_minimum,
        "observed-*.json branch_risk distance_to_seam_in_sigma",
        "smallest distance of an observed difference from the +-90 deg branch cut (marginal sigma)",
    )


def add_limits(numbers, results):
    alphas = results["observed-std-1d"]["alphas"]
    for data, data_word in DATA_WORDS.items():
        for dim, dim_word in DIMENSIONS.items():
            for mode, mode_word in MODE_WORDS.items():
                c95 = results[f"{data}-{mode}-{dim}"]["c95"]
                for alpha, value in zip(alphas, c95):
                    numbers.add(
                        f"NumLimit{data_word}{dim_word}{mode_word}{alpha_word(alpha)}",
                        sig(value, 3),
                        value,
                        f"analysis/limits report.json results {data}-{mode}-{dim}",
                        f"c95, {data}, {mode}, {dim}, alpha={alpha:.2f}",
                    )
    for alpha in (alphas[0], alphas[-1]):
        index = alphas.index(alpha)
        values = []
        for dim, dim_word in DIMENSIONS.items():
            value = results[f"observed-cor_fixed_background-{dim}"]["c95"][index]
            values.append(value)
            numbers.add(
                f"NumBaselineApprox{dim_word}{alpha_word(alpha)}",
                approx(value),
                value,
                f"report.json results observed-cor_fixed_background-{dim}",
                f"baseline c95 at alpha={alpha:.2f}, {dim}, rounded",
            )
        for label, value in (("Min", min(values)), ("Max", max(values))):
            numbers.add(
                f"NumBaselineApprox{alpha_word(alpha)}{label}",
                approx(value),
                value,
                "report.json results observed-cor_fixed_background-*",
                f"{label.lower()} baseline c95 over 1D and 2D at alpha={alpha:.2f}, rounded",
            )
    observed = [results[f"observed-{mode}-{dim}"]["c95"][0] for dim in DIMENSIONS for mode in MODES]
    for label, value in (("Min", min(observed)), ("Max", max(observed))):
        numbers.add(
            f"NumLimitRange{alpha_word(alphas[0])}{label}",
            sig(value, 2),
            value,
            "report.json results observed-*",
            f"{label.lower()} observed c95 over all cases at the smallest alpha",
        )
    for dim, dim_word in DIMENSIONS.items():
        for mode in ("cor_fixed_background", "cor"):
            ratio = np.asarray(results[f"observed-{mode}-{dim}"]["c95"]) / np.asarray(
                results[f"zero-data-{mode}-{dim}"]["c95"]
            )
            for label, value in (("Min", ratio.min()), ("Max", ratio.max())):
                numbers.add(
                    f"NumRatio{MODE_WORDS[mode]}{dim_word}{label}",
                    sig(value, 2),
                    value,
                    f"report.json results observed/zero-data-{mode}-{dim}",
                    f"{label.lower()} over alpha of c95(observed)/c95(zero data), {mode}, {dim}",
                )


def chi2_rows(limits_dir):
    """Saved joint-grid fits at PAPER_ALPHA: {(dim, mode): (components, chi2 at c=0, best chi2)}."""
    rows = {}
    for dim in DIMENSIONS:
        for mode in MODES:
            saved = read_json(limits_dir / f"observed-{mode}-{dim}.json")
            index = int(np.flatnonzero(np.isclose(saved["alphas"], PAPER_ALPHA))[0])
            fit = saved["fits"][index]
            rows[(dim, mode)] = (
                2 * len(saved["selected_indices"]),
                float(fit["data_quadratic"]),
                float(fit["joint_grid_chi2"]),
            )
    return rows


def add_chi2(numbers, rows):
    gains = {}
    for (dim, mode), (_n, zero, best) in rows.items():
        name = f"{DIMENSIONS[dim]}{MODE_WORDS[mode]}"
        source = f"analysis/limits observed-{mode}-{dim}.json fits at alpha={PAPER_ALPHA:.2f}"
        numbers.add(
            f"NumChiZero{name}",
            f"{zero:.1f}",
            zero,
            source + " data_quadratic",
            f"chi^2 at c=0, {mode}, {dim}",
        )
        numbers.add(
            f"NumChiBest{name}",
            f"{best:.1f}",
            best,
            source + " joint_grid_chi2",
            f"chi^2 at the best grid point, {mode}, {dim}",
        )
        numbers.add(
            f"NumChiGain{name}",
            f"{zero - best:.1f}",
            zero - best,
            source,
            f"chi^2 improvement, {mode}, {dim}",
        )
        gains[(dim, mode)] = zero - best
    case = max(gains, key=gains.get)
    numbers.add(
        "NumChiGainMax",
        f"{gains[case]:.1f}",
        dict(value=gains[case], dimension=case[0], mode=case[1]),
        "largest NumChiGain*",
        f"largest chi^2 improvement (case {case[1]}, {case[0]})",
    )


def add_checks(numbers, grid_report, prior_report, covariance_rows):
    grid = grid_report["max_delta_dex"]
    numbers.add(
        "NumGridCheckMaxDex",
        sig(grid, 2),
        grid,
        "analysis/check-grid report.json max_delta_dex",
        "largest |delta log10 c95| of the dense-grid check",
    )
    numbers.add(
        "NumGridCheckMaxPercent",
        f"{percent_change(grid):.1f}",
        percent_change(grid),
        "check-grid report.json max_delta_dex",
        "the same as a change of c95 (percent)",
    )
    prior = max(
        abs(a - b)
        for key in grid_report["results"]
        for a, b in zip(prior_report["results"][key]["log10_c95"], grid_report["results"][key]["log10_c95"])
    )
    numbers.add(
        "NumPriorCheckMaxDex",
        sig(prior, 1),
        prior,
        "check-prior vs check-grid report.json log10_c95",
        "largest |delta log10 c95| of the wider prior relative to the dense grid",
    )
    largest = {}
    for row in covariance_rows:
        key = (row["dimension"], row["mode"])
        largest[key] = max(largest.get(key, 0.0), abs(float(row["delta_dex"])))
    for (dim, mode), dex in sorted(largest.items()):
        name = f"NumCovCheck{DIMENSIONS[dim]}{MODE_WORDS[mode]}"
        numbers.add(
            name + "Percent",
            f"{percent_change(dex):.1f}",
            percent_change(dex),
            "analysis/check-covariance covariance-sensitivity.csv delta_dex",
            f"largest change of c95 with half the draws, {mode}, {dim} (percent)",
        )
        numbers.add(
            name + "PercentApprox",
            f"{percent_change(dex):.0f}",
            percent_change(dex),
            "covariance-sensitivity.csv delta_dex",
            f"the same, rounded, {mode}, {dim}",
        )
    std = max(largest[(dim, "std")] for dim in DIMENSIONS)
    numbers.add(
        "NumCovCheckStdMaxDex",
        sig(std, 2),
        std,
        "covariance-sensitivity.csv delta_dex (std)",
        "largest |delta log10 c95| of the diagonal model with half the draws",
    )
    numbers.add(
        "NumCovCheckStdMaxPercent",
        f"{percent_change(std):.1f}",
        percent_change(std),
        "covariance-sensitivity.csv delta_dex (std)",
        "the same as a change of c95 (percent)",
    )
    baseline = max(largest[(dim, "cor_fixed_background")] for dim in DIMENSIONS)
    numbers.add(
        "NumCovCheckCorBgMaxPercentApprox",
        f"{percent_change(baseline):.0f}",
        percent_change(baseline),
        "covariance-sensitivity.csv delta_dex (cor_fixed_background)",
        "largest change of the baseline c95 with half the draws, rounded (percent)",
    )


def add_two_emitter(numbers, figures):
    p = TWO_EMITTER
    errors = []
    for ratio, chi12, phi12 in p["check_points"]:
        k = figures.kappa(ratio, chi12, phi12)
        direct = figures.kappa_direct(ratio, chi12, phi12, p["check_amplitude"], p["check_phases"])
        errors.append(abs(k - direct) / k)
    largest = max(errors)
    bound = 10.0 ** np.ceil(np.log10(largest))
    numbers.add(
        "NumTwoEmitterCheckAmplitude",
        order(p["check_amplitude"]),
        p["check_amplitude"],
        "analysis_paper.py TWO_EMITTER",
        "amplitude of the direct two-emitter check (math mode)",
    )
    numbers.add(
        "NumTwoEmitterCheckBound",
        order(bound),
        dict(bound=bound, largest_relative_error=largest),
        "kappa formula vs direct evaluation",
        "bound on the relative error of the two-emitter formula (math mode)",
    )


# ----------------------------------------------------------------------------- tables
def table_limits(results, features):
    alphas = results["observed-std-1d"]["alphas"]
    lines = [
        r"\begin{tabular}{c ccc ccc}",
        r"\hline",
        rf" & \multicolumn{{3}}{{c}}{{{features['1d']} sectors}}"
        rf" & \multicolumn{{3}}{{c}}{{{features['2d']} pixels}} \\",
        r"$\alpha$ & std & cor & cor+bg & std & cor & cor+bg \\",
        r"\hline",
    ]
    for index, alpha in enumerate(alphas):
        values = [
            sig(results[f"observed-{mode}-{dim}"]["c95"][index], 3) for dim in DIMENSIONS for mode in MODES
        ]
        lines.append(f"{alpha:.2f} & " + " & ".join(values) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def table_chi2(rows):
    lines = [
        r"\begin{tabular}{l c ccc ccc}",
        r"\hline",
        r" & & \multicolumn{3}{c}{zero signal} & \multicolumn{3}{c}{best fit} \\",
        r"Observable & $N$ & std & cor & cor+bg & std & cor & cor+bg \\",
        r"\hline",
    ]
    for dim, label in DIMENSIONS.items():
        n = rows[(dim, MODES[0])][0]
        zero = [f"{rows[(dim, m)][1]:.1f}" for m in MODES]
        best = [f"{rows[(dim, m)][2]:.1f}" for m in MODES]
        lines.append(f"{label} & {n} & " + " & ".join(zero + best) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------- run
def background_evpa(frames, dimension, operator):
    """EVPA of every feature of every axion-off frame [t, feature]; 1D uses the aperture operator."""
    values = []
    for frame in frames:
        frame = np.asarray(frame[0])
        if dimension == "1d":
            values.append(_integrate(frame, operator))
        else:
            q, u = frame[..., 1].ravel(), frame[..., 2].ravel()
            theta = 0.5 * np.arctan2(u, q)
            theta[(q == 0) & (u == 0)] = np.nan
            values.append(theta)
    return np.asarray(values)


def run(args):
    import _paper_figures as figures
    from analysis_limits import theory_aperture_operator
    from eht_axion.evpa import circular_difference

    destination = args.output_root.resolve()
    if destination.exists():
        raise FileExistsError(f"New output root required: {destination}")
    config = json.loads(json.dumps(read_experiment(args.config)))
    inputs, numbers = Inputs(), Numbers()
    limits = args.limits.resolve()
    with provenance.verification_session():
        files = ["report.json"] + [
            f"{p}-{d}.{e}"
            for p in ("observations", "background-selected")
            for d in DIMENSIONS
            for e in ("npz", "json")
        ]
        files += [f"observed-{m}-{d}.json" for d in DIMENSIONS for m in MODES]
        inputs.stage(limits, "analysis:limits", files)
        report = read_json(limits / "report.json")
        if report["scientific_profile"] != config:
            raise ValueError("The limits were produced from a different experiment.yaml")
        reports = {}
        for name, path, stage in (
            ("grid", args.gridcheck, "analysis:check-grid"),
            ("prior", args.priorcheck, "analysis:check-grid"),
        ):
            inputs.stage(path, stage, ["report.json"])
            reports[name] = read_json(Path(path) / "report.json")
        inputs.stage(
            args.covcheck, "analysis:check-covariance", ["covariance-sensitivity.csv", "report.json"]
        )
        with (Path(args.covcheck) / "covariance-sensitivity.csv").open(newline="") as stream:
            covariance_rows = list(csv.DictReader(stream))
        layout = RunLayout(args.templates.resolve())
        inputs.stage(layout.processed(PAPER_ALPHA), "analysis:preprocess", ["stokes.npy"])
        inputs.stage(layout.response(PAPER_ALPHA), "analysis:fit", ["response-2d.npz"])
        for dim in DIMENSIONS:
            inputs.stage(layout.templates(dim), "analysis:templates", ["templates.npz"])
        inputs.stage(layout.background_processed, "analysis:background", ["stokes.npy"])
        ensemble = args.ensemble.resolve()
        inputs.stage(
            ensemble,
            "eht:ensemble",
            ["report.json", "axes.npz"] + [f"day-{day}/global.npz" for day in EHT_DAYS],
        )
        shards = input_raptor_files(
            args.raptor_input_root.resolve(), None, config["raptor"]["background"]["resolution"]
        )
        for shard in shards:
            inputs.stage(shard.parent, "raptor:delivery", [shard.name, "generated.json"])
        build = args.raptor_run_root.resolve() / "build"
        (axion_source,) = inputs.stage(build, "raptor:build", ["source/src/axion.c"])
        plan = read_json(inputs.file(args.raptor_plan))
        tables = radial_tables(args, plan, inputs)

        obs = {d: npz(limits / f"observations-{d}.npz") for d in DIMENSIONS}
        info = {d: read_json(limits / f"observations-{d}.json") for d in DIMENSIONS}
        background = {d: npz(limits / f"background-selected-{d}.npz") for d in DIMENSIONS}
        background_info = read_json(limits / "background-selected-1d.json")
        observed = {(m, d): read_json(limits / f"observed-{m}-{d}.json") for d in DIMENSIONS for m in MODES}
        features = {d: len(obs[d]["selected_indices"]) for d in DIMENSIONS}
        rows = chi2_rows(limits)

        times = []
        flux = []
        scale = physical_scale(config)
        rg_radians = np.deg2rad(scale.angular_rg_uas / 3.6e9)
        n400 = config["raptor"]["background"]["resolution"]
        solid_angle = (config["raptor"]["camera_fov_rg"] / n400 * rg_radians) ** 2
        for shard in shards:
            times += read_json(shard.parent / "generated.json")["times_tg"]
            cube = np.load(shard, mmap_mode="r", allow_pickle=False)
            flux += list(np.asarray(cube[:, 0, ..., 0]).sum(axis=(1, 2)) * solid_angle)
        flux = np.asarray(flux)

        add_scales(numbers, config, times, background_info["lag_hours"], len(background["1d"]["differences"]))
        profiles = radial_profiles(config, tables)
        window = window_constants(axion_source)
        add_radial(numbers, profiles, window)

        cloud = config["analysis"]["cloud"]
        b = cloud_amplitude_scaling(
            PAPER_ALPHA,
            coefficient=cloud["coefficient"],
            pivot_alpha=cloud["pivot_alpha"],
            exponent=cloud["exponent"],
        )
        per_c = b / config["raptor"]["response"]["reference_norm"]
        response = npz(layout.response(PAPER_ALPHA) / "response-2d.npz")
        index = int(np.flatnonzero(response["snapshot_labels"] == PAPER_SNAPSHOT_LABEL)[0])
        processed = np.load(layout.processed(PAPER_ALPHA) / "stokes.npy", mmap_mode="r", allow_pickle=False)
        frame = np.asarray(processed[index, 0])
        amplitude_deg = np.rad2deg(response["amplitude"][index] * per_c)
        templates = {d: npz(layout.templates(d) / "templates.npz") for d in DIMENSIONS}
        j = {d: int(np.flatnonzero(np.isclose(templates[d]["alphas"], PAPER_ALPHA))[0]) for d in DIMENSIONS}
        sector_mean = np.rad2deg(templates["1d"]["amplitude_mean"][j["1d"]] * per_c)
        sector_fraction = (
            templates["1d"]["amplitude_std"][j["1d"]] / templates["1d"]["amplitude_mean"][j["1d"]]
        )
        add_signal(numbers, amplitude_deg, sector_mean, sector_fraction)

        axes = npz(ensemble / "axes.npz")
        ensemble_report = read_json(ensemble / "report.json")
        if ensemble_report["bunit"].upper() != "JY/PIXEL":
            raise ValueError("Reconstruction flux requires images in Jy/pixel")
        background_frames = np.load(
            layout.background_processed / "stokes.npy", mmap_mode="r", allow_pickle=False
        )
        model_m, model_beta = [], []
        for item in background_frames:
            values = global_diagnostics(
                np.asarray(item[0]).transpose(1, 0, 2), axes["east_uas"], axes["north_uas"]
            )
            model_m.append(values[2])
            model_beta.append(abs(complex(values[4], values[5])))
        recon = {day: npz(ensemble / f"day-{day}" / "global.npz") for day in EHT_DAYS}
        column = {day: list(recon[day]["columns"]) for day in EHT_DAYS}
        recon_m = [recon[d]["values"][:, column[d].index("m_resolved")].mean() for d in EHT_DAYS]
        recon_beta = [
            np.abs(
                recon[d]["values"][:, column[d].index("beta2_real")]
                + 1j * recon[d]["values"][:, column[d].index("beta2_imag")]
            ).mean()
            for d in EHT_DAYS
        ]
        recon_flux = np.mean([recon[d]["values"][:, column[d].index("total_I")].mean() for d in EHT_DAYS])
        add_model_vs_reconstruction(
            numbers, flux, recon_flux, np.asarray(model_m), np.asarray(model_beta), recon_m, recon_beta
        )
        frequency = ensemble_report["frequency_hz"] / 1e9
        numbers.add(
            "NumObsFrequencyGHz",
            f"{frequency:.2f}",
            frequency,
            "eht/ensemble report.json frequency_hz",
            "reference frequency of the EHT data (GHz)",
        )

        branch = min(
            min(entry["distance_to_seam_in_sigma"])
            for saved in observed.values()
            for entry in saved["branch_risk"]
        )
        add_observations(numbers, obs, info, background, config, branch)
        add_limits(numbers, report["results"])
        add_chi2(numbers, rows)
        add_checks(numbers, reports["grid"], reports["prior"], covariance_rows)
        add_two_emitter(numbers, figures)

        # Unmasked background: every sector / pixel with a defined EVPA in all frames.
        operator = theory_aperture_operator(config)
        lag = config["analysis"]["background"]["lag_frames"]
        unmasked = {}
        for dim in DIMENSIONS:
            theta = background_evpa(background_frames, dim, operator)
            theta = theta[:, np.isfinite(theta).all(axis=0)]
            unmasked[dim] = np.rad2deg(circular_difference(theta[lag:], theta[:-lag])).ravel()

        consumed = sorted(set(inputs.consumed))
        identity = dict(
            schema="analysis-paper-v1",
            sources=dict(
                limits=str(limits),
                gridcheck=str(Path(args.gridcheck).resolve()),
                priorcheck=str(Path(args.priorcheck).resolve()),
                covcheck=str(Path(args.covcheck).resolve()),
                templates=str(layout.root),
                ensemble=str(ensemble),
                raptor_input_root=str(args.raptor_input_root.resolve()),
                raptor_plan=str(Path(args.raptor_plan).resolve()),
                raptor_build=str(build),
            ),
            presentation=dict(
                paper_alpha=PAPER_ALPHA,
                paper_snapshot_label=PAPER_SNAPSHOT_LABEL,
                ticks=TICKS,
                radial_sample_rg=RADIAL_SAMPLE_RG,
                window_sample_rg=WINDOW_SAMPLE_RG,
                orbit_radius_rg=ORBIT_RADIUS_RG,
                close_separation_uas=CLOSE_SEPARATION_UAS,
                growth_rate_denominator=GROWTH_RATE_DENOMINATOR,
                two_emitter={k: v for k, v in TWO_EMITTER.items()},
            ),
            verification="completed stage outputs read; raw lineage inherited, not rehashed; "
            "nothing refitted",
        )
        with provenance.stage_directory(
            destination, "analysis:paper", identity, consumed, verify_lineage=False
        ) as out:
            figure_dir = out / "figures"
            figure_dir.mkdir()
            import matplotlib.pyplot as plt

            with plt.rc_context(figures.STYLE):
                high = [alpha for alpha, _s, is_high, _r, _v in profiles if is_high]
                low = [alpha for alpha, _s, is_high, _r, _v in profiles if not is_high]
                switch = (max(low) + min(high)) / 2 if high and low else None
                horizon = 1 + np.sqrt(1 - config["raptor"]["fluid_spin"] ** 2)
                figures.radial_function(figure_dir, profiles, window, horizon, switch, RADIAL_SAMPLE_RG)
                figures.two_photon(figure_dir, TWO_EMITTER["left_ratio"], TWO_EMITTER["right_phi_degrees"])
                fov = config["shared"]["image"]["fov_uas"]
                figures.signal_maps(figure_dir, frame, amplitude_deg, response["phase"][index], fov, TICKS)
                mean = templates["2d"]["amplitude_mean"][j["2d"]]
                fraction = templates["2d"]["amplitude_std"][j["2d"]] / mean
                spread = np.rad2deg(templates["1d"]["amplitude_std"][j["1d"]] * per_c)
                figures.template_maps(
                    figure_dir, np.rad2deg(mean * per_c), fraction, sector_mean, spread, fov
                )
                selected = {
                    d: (
                        np.rad2deg(background[d]["differences"]),
                        float(np.rad2deg(background[d]["mean_radians"])),
                        float(np.rad2deg(background[d]["sigma_radians"])),
                    )
                    for d in DIMENSIONS
                }
                figures.background(figure_dir, selected, unmasked)
                entries = {}
                for d in DIMENSIONS:
                    sigma, correlation = dterm_statistics(obs[d])
                    coordinates = obs[d]["coordinates"]
                    entries[d] = dict(
                        x=coordinates[:, 0] if d == "1d" else np.arange(1, features[d] + 1),
                        sigma_deg=np.rad2deg(sigma),
                        correlation=correlation,
                        separation=np.linalg.norm(coordinates[:, None] - coordinates[None], axis=-1),
                    )
                figures.observational_uncertainty(figure_dir, entries)
                fits = {}
                for d in DIMENSIONS:
                    saved = observed[("cor_fixed_background", d)]
                    k = int(np.flatnonzero(np.isclose(saved["alphas"], PAPER_ALPHA))[0])
                    covariance = np.asarray(saved["effective_pair_covariances"])
                    fits[d] = dict(
                        x=entries[d]["x"],
                        data_deg=np.rad2deg(obs[d]["pair_differences"]),
                        sigma_deg=np.rad2deg(np.sqrt(np.diagonal(covariance, axis1=-2, axis2=-1))),
                        prediction_deg=np.rad2deg(np.asarray(saved["fits"][k]["prediction_radians"])),
                    )
                figures.data_vs_model(
                    figure_dir, fits, PAPER_ALPHA, config["analysis"]["features"]["sector_spacing_degrees"]
                )
                curves = {
                    (data, mode, d): (
                        report["results"][f"{data}-{mode}-{d}"]["alphas"],
                        report["results"][f"{data}-{mode}-{d}"]["c95"],
                    )
                    for data in DATA_WORDS
                    for mode in MODES
                    for d in DIMENSIONS
                }
                figures.limits(figure_dir, curves, features)
            (out / "numbers.tex").write_text(numbers.tex())
            (out / "numbers.json").write_text(json.dumps(numbers.json(), indent=2, allow_nan=False) + "\n")
            (out / "table-limits.tex").write_text(table_limits(report["results"], features))
            (out / "table-chi2.tex").write_text(table_chi2(rows))
    print(f"PAPER OUTPUT {destination}: {len(numbers.items)} numbers, {len(figures.FIGURES)} figures")
    return destination


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", type=Path, default=DEFAULT_PATH, help="default: configs/experiment.yaml")
    for name, default in (
        ("limits", "analysis/limits"),
        ("gridcheck", "analysis/check-grid"),
        ("priorcheck", "analysis/check-prior"),
        ("covcheck", "analysis/check-covariance"),
        ("templates", "analysis/templates"),
        ("ensemble", "eht/ensemble"),
        ("raptor-input-root", "raptor/delivered/data"),
        ("raptor-plan", "raptor/plan.json"),
        ("raptor-run-root", "raptor/runs/main"),
        ("output-root", "analysis/paper"),
    ):
        p.add_argument("--" + name, type=Path, help=f"default: {default}")
    p.add_argument(
        "--radial-low", type=Path, help="overrides $EHT_RADIAL_LOW and the path in the RAPTOR plan"
    )
    p.add_argument(
        "--radial-high", type=Path, help="overrides $EHT_RADIAL_HIGH and the path in the RAPTOR plan"
    )
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(
        p,
        args,
        {
            "limits": ("analysis_limits", {}),
            "gridcheck": ("analysis_check_grid", {}),
            "priorcheck": ("analysis_check_prior", {}),
            "covcheck": ("analysis_check_covariance", {}),
            "templates": ("analysis_templates", {}),
            "ensemble": ("eht_ensemble", {}),
            "raptor_input_root": ("raptor_data", {}),
            "raptor_plan": ("raptor_plan", {}),
            "raptor_run_root": ("raptor_run", {}),
            "output_root": ("analysis_paper", {}),
        },
    )
    os.environ.setdefault("MPLBACKEND", "Agg")
    run(args)


if __name__ == "__main__":
    main()
