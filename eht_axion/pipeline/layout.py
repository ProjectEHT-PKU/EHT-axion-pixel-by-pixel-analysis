"""Canonical dataset identities and paths, independent of publication status."""

from dataclasses import dataclass
from datetime import date, timedelta
from decimal import Decimal, DecimalException
from pathlib import Path
import re
import numpy as np

DAYS = (95, 96, 100, 101)


def date_label(day):
    if int(day) not in DAYS:
        raise ValueError(f"Unsupported EHT observation day: {day}")
    return (date(2017, 1, 1) + timedelta(days=int(day) - 1)).isoformat()


def alpha_label(alpha):
    try:
        value = Decimal(str(alpha))
        if not value.is_finite() or value <= 0 or value != value.quantize(Decimal("0.01")):
            raise ValueError("Alpha identities require positive finite two-decimal values")
    except DecimalException as exc:
        raise ValueError("Alpha identities require positive finite two-decimal values") from exc
    return f"alpha-{value:.2f}"


_SHARD_NAME = re.compile(r"stokes-t([0-9]+)-t([0-9]+)\.npy")
_LABEL_STEP = 10  # Snapshot-label convention, not a verified physical cadence.


def _cube_header(path):
    """Inspect only a read-only mmap/header; never scan or materialize pixels."""
    try:
        array = np.load(path, mmap_mode="r", allow_pickle=False)
    except (ValueError, OSError, EOFError) as exc:
        raise ValueError(f"Invalid RAPTOR input header: {path}: {exc}") from exc
    try:
        if not hasattr(array, "shape") or not hasattr(array, "dtype"):
            # np.load can identify a ZIP archive regardless of a .npy suffix.
            array.close()
            raise ValueError(f"RAPTOR input must be a single NPY cube, not an archive: {path}")
        shape, dtype = array.shape, array.dtype
        if (len(shape) != 5 or any(n == 0 for n in shape) or shape[-1] != 4
                or shape[-2] != shape[-3] or shape[-2] < 2):
            raise ValueError(f"RAPTOR cube must be (T,S,N,N,4), T,S > 0, N >= 2: {path}: {shape}")
        if dtype.kind not in "iuf":
            raise ValueError(f"RAPTOR cube requires real numeric dtype: {path}: {dtype}")
        return shape, dtype
    finally:
        # Release each mapping before opening the next (even with many shards).
        del array


def _ordered_shards(paths):
    shards = []
    for path in paths:
        match = _SHARD_NAME.fullmatch(path.name)
        if match is None:
            raise ValueError(f"Invalid RAPTOR shard name; expected stokes-tSTART-tEND.npy: {path}")
        start, end = map(int, match.groups())
        if start > end or start % _LABEL_STEP or end % _LABEL_STEP:
            raise ValueError(f"RAPTOR shard labels require ordered endpoints on step {_LABEL_STEP}: {path}")
        shards.append((start, end, path))
    shards.sort(key=lambda item: (item[0], item[1]))
    previous_end, reference, identities = None, None, set()
    for start, end, path in shards:
        if previous_end is not None and start != previous_end + _LABEL_STEP:
            raise ValueError(f"RAPTOR shards overlap or have a gap; labels must be continuous: {path}")
        identity = path.stat().st_dev, path.stat().st_ino
        if identity in identities:
            raise ValueError(f"Duplicate RAPTOR input shard source: {path}")
        identities.add(identity)
        shape, dtype = _cube_header(path)
        expected = (end - start) // _LABEL_STEP + 1  # Inclusive endpoints.
        if shape[0] != expected:
            raise ValueError(f"RAPTOR shard time length {shape[0]} != label frame count {expected}: {path}")
        current = shape[1:], dtype
        if reference is not None and current != reference:
            raise ValueError(f"RAPTOR shards must have identical non-time shape and dtype: {path}")
        reference, previous_end = current, end
    return [path for _, _, path in shards]


def input_raptor_files(root, alpha, resolution):
    """Resolve automatic sources; res400 rejects ambiguity before choosing data.

    Shard endpoints are inclusive snapshot labels on the current step-10 grid.
    Neither spatial resolution nor a production time length (e.g. 141) is
    imposed on synthetic cubes. Single files carry no filename time labels.
    """
    root = Path(root)
    canonical = root / "raptor" / f"res-{resolution}" / alpha_label(alpha)
    legacy = root / "raptor-stokes" / f"res-{resolution}"
    if resolution == 100:
        # Preserve the existing canonical-before-legacy identity priority.
        candidates = [canonical / "stokes.npy", legacy / f"alpha-{float(alpha):g}.npy"]
        for path in candidates:
            if path.is_file():
                _cube_header(path)
                return [path]
    elif resolution == 400:
        sources = []
        for directory in (canonical, legacy / f"alpha-{float(alpha):g}"):
            whole = directory / "stokes.npy"
            if whole.is_file():
                sources.append(("whole", [whole]))
            # Include malformed time-like names so they cannot disappear silently.
            paths = list(directory.glob("stokes-t*.npy"))
            if paths:
                sources.append(("shards", paths))
        if len(sources) > 1:
            locations = [str(paths[0].parent) + f" ({kind})" for kind, paths in sources]
            raise ValueError(f"Ambiguous res400 RAPTOR input sources: {', '.join(locations)}")
        if sources:
            kind, paths = sources[0]
            if kind == "shards":
                return _ordered_shards(paths)
            _cube_header(paths[0])
            return paths
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
        if identity["observation"] not in {"null", "expected"}:
            raise ValueError("Only zero-day-difference results are supported")
        if cuts.shape != (3, len(alpha)):
            raise ValueError("Limit array must have shape (mode, alpha)")
        return alpha, cuts, identity
    arrays, reference, identity = [], None, None
    for mode in modes:
        with np.load(directory / f"{mode}_{dimension}.npz", allow_pickle=False) as z:
            alpha, cuts = z["alphas"].copy(), z["cuts"].copy()
            current = {k: str(z[k].item()) for k in ("variant", "observation", "prior")}
        if current["observation"] not in {"null", "expected"}:
            raise ValueError("Only zero-day-difference results are supported")
        if reference is not None and (not np.array_equal(reference, alpha) or identity != current):
            raise ValueError("Incompatible legacy curve files")
        reference, identity = alpha, current
        arrays.append(cuts)
    return reference, np.stack(arrays), identity
