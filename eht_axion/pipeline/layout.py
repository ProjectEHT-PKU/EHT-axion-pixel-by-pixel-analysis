"""Canonical dataset identities and paths, independent of publication status."""

from dataclasses import dataclass
from datetime import date, timedelta
from decimal import Decimal
from pathlib import Path
import numpy as np

DAYS = (95, 96, 100, 101)


def date_label(day):
    if int(day) not in DAYS:
        raise ValueError(f"Unsupported EHT observation day: {day}")
    return (date(2017, 1, 1) + timedelta(days=int(day) - 1)).isoformat()


def alpha_label(alpha):
    value = Decimal(str(alpha))
    if not value.is_finite() or value <= 0 or value != value.quantize(Decimal("0.01")):
        raise ValueError("Alpha identities require positive finite two-decimal values")
    return f"alpha-{value:.2f}"


def input_raptor_files(root, alpha, resolution):
    root = Path(root)
    canonical = root / "raptor" / f"res-{resolution}" / alpha_label(alpha)
    legacy = root / "raptor-stokes" / f"res-{resolution}"
    if resolution == 100:
        candidates = [canonical / "stokes.npy", legacy / f"alpha-{float(alpha):g}.npy"]
        for path in candidates:
            if path.is_file():
                return [path]
    elif resolution == 400:
        for directory in (canonical, legacy / f"alpha-{float(alpha):g}"):
            paths = sorted(directory.glob("stokes-t*-t*.npy"))
            if paths:
                return paths
        if (canonical / "stokes.npy").is_file():
            return [canonical / "stokes.npy"]
    raise FileNotFoundError(f"No RAPTOR input for {alpha_label(alpha)}, res={resolution}: {root}")


def input_eht_file(root, day):
    root = Path(root)
    candidates = [
        root / "eht" / date_label(day) / "image-ensemble.pkl",
        root / date_label(day) / "image-ensemble.pkl",
        root / "eht-image-ensemble" / f"230_{day}.pkl",
        root / f"230_{day}.pkl",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"No EHT ensemble for {date_label(day)} in {root}")


@dataclass(frozen=True)
class RunLayout:
    root: Path

    def __post_init__(self):
        object.__setattr__(self, "root", Path(self.root))

    def processed(self, alpha):
        return self.root / "raptor" / "processed" / alpha_label(alpha)

    def response(self, alpha):
        return self.root / "raptor" / "response" / alpha_label(alpha)

    def templates(self, dimension):
        return self.root / "raptor" / "templates" / dimension

    @property
    def background_processed(self):
        return self.root / "raptor" / "background"

    @property
    def eht(self):
        return self.root / "eht"

    @property
    def background(self):
        return self.root / "background"

    def observations(self, dimension):
        return self.root / "observations" / dimension

    @property
    def inference(self):
        return self.root / "inference"

    @property
    def figures(self):
        return self.root / "figures"


def response_arrays(path):
    """Read named response fields or the historical final-axis triplet."""
    path = Path(path)
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as z:
            return tuple(z[k].copy() for k in ("amplitude", "phase", "baseline"))
    array = np.load(path, mmap_mode="r", allow_pickle=False)
    if array.shape[-1] != 3:
        raise ValueError("Historical response requires final [A, phase, baseline] axis")
    return tuple(array[..., k] for k in range(3))


def load_curves(directory, dimension):
    """Return modes/alphas/cuts/identity with a legacy six-file read adapter."""
    directory = Path(directory)
    modes = ("std", "cor", "cor_grmhd")
    combined = directory / f"limits-{dimension}.npz"
    if combined.is_file():
        with np.load(combined, allow_pickle=False) as z:
            if tuple(z["modes"].tolist()) != modes:
                raise ValueError("Unexpected limit-mode order")
            alpha, cuts = z["alphas"].copy(), z["cuts"].copy()
            identity = {k: str(z[k].item()) for k in ("variant", "observation", "prior")}
        if cuts.shape != (3, len(alpha)):
            raise ValueError("Limit array must have shape (mode, alpha)")
        return alpha, cuts, identity
    arrays, reference, identity = [], None, None
    for mode in modes:
        with np.load(directory / f"{mode}_{dimension}.npz", allow_pickle=False) as z:
            alpha, cuts = z["alphas"].copy(), z["cuts"].copy()
            current = {k: str(z[k].item()) for k in ("variant", "observation", "prior")}
        if reference is not None and (not np.array_equal(reference, alpha) or identity != current):
            raise ValueError("Incompatible legacy curve files")
        reference, identity = alpha, current
        arrays.append(cuts)
    return reference, np.stack(arrays), identity
