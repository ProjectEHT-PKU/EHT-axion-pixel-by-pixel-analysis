"""Directory layout of the RAPTOR image cubes (res-*/{alpha-*,background}/)."""

from decimal import Decimal, DecimalException
from pathlib import Path



def alpha_label(alpha):
    try:
        value = Decimal(str(alpha))
        if not value.is_finite() or value <= 0 or value != value.quantize(Decimal("0.01")):
            raise ValueError("Alpha identities require positive finite two-decimal values")
    except DecimalException as exc:
        raise ValueError("Alpha identities require positive finite two-decimal values") from exc
    return f"alpha-{value:.2f}"


def product_directory(root, alpha, resolution):
    """root/res-<resolution>/alpha-<alpha>, or .../background for alpha None."""
    name = "background" if alpha is None else alpha_label(alpha)
    return Path(root) / f"res-{resolution}" / name


def image_cube(root, alpha, resolution):
    """The Stokes image cube (T, S, N, N, 4) of one alpha (None: axion-off background)."""
    path = product_directory(root, alpha, resolution) / "stokes.npy"
    if not path.is_file():
        raise FileNotFoundError(f"No RAPTOR product: {path}")
    return path
