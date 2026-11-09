"""Reader for the response-template summaries written by pipeline.raptor.summarize_fits."""

import numpy as np


def load_templates(path, dimension):
    with np.load(path, allow_pickle=False) as z:
        alphas = z["alphas"].copy()
        keys = ["amplitude_mean", "amplitude_std", "phase_mean", "phase_std"]
        values = [z[k].copy() for k in keys]
    if dimension not in {"1d", "2d"} or alphas.ndim != 1 or not len(alphas):
        raise ValueError("Invalid template dimension/alpha grid")
    if not np.isfinite(alphas).all() or np.any(np.diff(alphas) <= 0):
        raise ValueError("Template alpha grid must be finite and increasing")
    expected_ndim = 2 if dimension == "1d" else 3
    if any(v.shape != values[0].shape or v.ndim != expected_ndim for v in values):
        raise ValueError("Template summary shapes differ")
    if values[0].shape[0] != len(alphas):
        raise ValueError("Template alpha axis mismatch")
    amplitude = {float(a): np.stack((values[0][k], values[1][k])) for k, a in enumerate(alphas)}
    phase = {float(a): np.stack((values[2][k], values[3][k])) for k, a in enumerate(alphas)}
    valid = np.ones(values[0].shape[1:], dtype=bool)
    for v in values:
        valid &= np.isfinite(v).all(axis=0)
    valid &= (values[1] >= 0).all(axis=0) & (values[3] >= 0).all(axis=0)
    return alphas, amplitude, phase, valid
