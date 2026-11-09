"""Directory layout of the delivered RAPTOR Stokes cubes (raptor/res-*/alpha-*/)."""

from decimal import Decimal, DecimalException
from pathlib import Path
import re

# The axion-off (state 0) background frames are delivered under this alpha
# directory name. It is only a directory label: no axion field enters them.
BACKGROUND_DIRECTORY_ALPHA = 0.3

_SHARD_NAME = re.compile(r"stokes-t([0-9]+)-t([0-9]+)\.npy")


def alpha_label(alpha):
    try:
        value = Decimal(str(alpha))
        if not value.is_finite() or value <= 0 or value != value.quantize(Decimal("0.01")):
            raise ValueError("Alpha identities require positive finite two-decimal values")
    except DecimalException as exc:
        raise ValueError("Alpha identities require positive finite two-decimal values") from exc
    return f"alpha-{value:.2f}"


def input_raptor_files(root, alpha, resolution):
    """Delivered Stokes cube(s) for one alpha/resolution, in time order.

    Response products are one ``stokes.npy``; background products are time
    shards ``tA-tB/stokes-tA-tB.npy``. Identity, labels and shapes are checked
    by raptor.delivery.generated_inputs, not inferred from file names.
    """
    base = Path(root) / "raptor" / f"res-{resolution}" / alpha_label(alpha)
    whole = base / "stokes.npy"
    shards = sorted(base.glob("t*/stokes-t*.npy"))
    if whole.is_file() and shards:
        raise ValueError(f"Ambiguous RAPTOR products (whole cube and shards): {base}")
    if whole.is_file():
        return [whole]
    if not shards:
        raise FileNotFoundError(f"No RAPTOR product for {alpha_label(alpha)}, res={resolution}: {root}")
    named = []
    for path in shards:
        match = _SHARD_NAME.fullmatch(path.name)
        if match is None or path.parent.name != f"t{match[1]}-t{match[2]}":
            raise ValueError(f"Invalid RAPTOR shard name: {path}")
        named.append((int(match[1]), path))
    return [path for _, path in sorted(named)]
