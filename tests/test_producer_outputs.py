"""Synthetic serialization checks only; no raytracing or real fluid inputs."""
from copy import deepcopy

import h5py
import numpy as np
import pytest

from eht_axion.producer.outputs import (
    CORE_CONSTANTS, hdf_mapping, read_core_constants, validate_job_output,
)
from eht_axion.pipeline.raptor import convert_hdf5


def make_task(source=None, *, phase=None, resolution=20):
    enabled = phase is not None
    p = dict(MBH=6.2e9, DISTANCE=16800., M_UNIT=3e25, R_LOW=1., R_HIGH=1., INCLINATION=163.,
             IMG_WIDTH=resolution, IMG_HEIGHT=resolution, CAM_SIZE_X=46.2261, CAM_SIZE_Y=46.2261,
             FREQS_PER_DEC=1., FREQ_MIN=230e9, STEPSIZE=.005, AXION_ALPHA=.2 if enabled else 0.,
             AXION_OMEGA=.197 if enabled else 0., AXION_NORM=.1 if enabled else 0.,
             AXION_PHASE=phase or 0, AXION_TABLE="disabled", MAX_LEVEL=1)
    return dict(id=f"response-{phase}" if enabled else "baseline", kind="response" if enabled else "baseline",
                label=1600, time_tg=1022., source_path=str(source), resolution=resolution,
                alpha=.2 if enabled else None, omega=.197 if enabled else 0., norm=.1 if enabled else 0.,
                phase_index=phase or 0, output_slot=phase + 1 if enabled else 0, parameters=p,
                table_id="low" if enabled else None, table_path=None)


def write_products(directory, p, slot):
    """Self-contained helper also embedded in the fake executable tests."""
    from pathlib import Path
    import h5py
    import numpy as np

    output = Path(directory) / "output"
    output.mkdir()
    n = p["IMG_WIDTH"]
    b = n // 10
    ix = np.array([10 * (block // b) + pixel // 10 for block in range(b*b) for pixel in range(100)])
    iy = np.array([10 * (block % b) + pixel % 10 for block in range(b*b) for pixel in range(100)])
    shape = (b*b, 100)
    # Nonsymmetric image catches transpose, mirror, lost channels and U sign flips.
    intensity = (1. + ix * 2 + iy * 3).reshape(shape) * 1e20
    arrays = [intensity, intensity * .2, intensity * -.3, intensity * .01]
    with h5py.File(output / f"img_data_{slot}.h5", "w") as h:
        for name, a in zip("IQUV", arrays):
            h[f"{name}{p['FREQ_MIN']:e}"] = a
        h[f"tau{p['FREQ_MIN']:e}"] = np.full(shape, .1)
        h[f"tauF{p['FREQ_MIN']:e}"] = np.full(shape, -.1)
        h["alpha"] = (-p["CAM_SIZE_X"] / 2 + (ix + .5) * p["CAM_SIZE_X"] / n).reshape(shape)
        h["beta"] = (-p["CAM_SIZE_Y"] / 2 + (iy + .5) * p["CAM_SIZE_Y"] / n).reshape(shape)
    rg = 6.674e-8 * (p["MBH"] * 1.989e33) / 2.99792458e10 / 2.99792458e10
    area = p["CAM_SIZE_X"] * p["CAM_SIZE_Y"] / n**2 * (rg / (p["DISTANCE"] * 3.086e21))**2
    row = [p["FREQ_MIN"], *[a.sum() * area for a in arrays]]
    np.savetxt(output / f"spectrum_{slot}_{p['INCLINATION']:.2f}.dat", [row], fmt="%+.15e")


@pytest.fixture
def products(tmp_path):
    task = make_task()
    write_products(tmp_path, task["parameters"], 0)
    return tmp_path, task


def test_valid_output_and_explicit_converter_mapping(products, tmp_path):
    directory, task = products
    metadata = validate_job_output(directory, task)
    assert metadata["intensity_units"] == "Jy/sr"
    assert metadata["core_constants"] == CORE_CONSTANTS
    assert metadata["axes"][:2] == ["x/alpha", "y/beta"]
    assert metadata["state_slot"] == 0
    result = convert_hdf5(directory / metadata["hdf5"], tmp_path / "converted", **metadata["mapping"])
    cube = np.load(result["stokes"])
    expected = (1 + np.arange(20)[:, None] * 2 + np.arange(20)[None, :] * 3) * 1e20
    for channel, factor in enumerate([1, .2, -.3, .01]):
        np.testing.assert_allclose(cube[..., channel], expected * factor)


@pytest.mark.parametrize("channel", ["I", "Q", "U", "V", "tau", "tauF", "alpha", "beta"])
def test_missing_or_nonfinite_channel_rejected(products, channel):
    directory, task = products
    name = channel if channel in {"alpha", "beta"} else f"{channel}2.300000e+11"
    with h5py.File(directory / "output/img_data_0.h5", "r+") as h:
        h[name][0, 0] = np.nan
    with pytest.raises(ValueError, match="Nonfinite"):
        validate_job_output(directory, task)
    with h5py.File(directory / "output/img_data_0.h5", "r+") as h:
        del h[name]
    with pytest.raises(ValueError, match="exactly"):
        validate_job_output(directory, task)


@pytest.mark.parametrize("mutation", ["block", "pixel", "transpose", "offset"])
def test_coordinate_order_and_centers_are_not_repaired(products, mutation):
    directory, task = products
    with h5py.File(directory / "output/img_data_0.h5", "r+") as h:
        a = h["alpha"][...]
        if mutation == "block":
            a = a[::-1]
        elif mutation == "pixel":
            a = a[:, ::-1]
        elif mutation == "transpose":
            a = a.reshape(-1, 10, 10).transpose(0, 2, 1).reshape(-1, 100)
        else:
            a += 1e-8
        h["alpha"][...] = a
    with pytest.raises(ValueError, match="centers/order"):
        validate_job_output(directory, task)


def test_flux_area_and_spectrum_precision(products):
    directory, task = products
    metadata = validate_job_output(directory, task)
    spec = directory / metadata["spectrum"]
    rows = np.loadtxt(spec, ndmin=2)
    rows[0, 1:] *= 1.000001
    np.savetxt(spec, rows, fmt="%.15e")
    with pytest.raises(ValueError, match="flux mismatch"):
        validate_job_output(directory, task)


@pytest.mark.parametrize("content", ["1 2 3 4 5\n1 2 3 4 5\n", "230e9 1 2\n", "nan 1 2 3 4\n"])
def test_malformed_spectrum(products, content):
    directory, task = products
    (directory / "output/spectrum_0_163.00.dat").write_text(content)
    with pytest.raises(ValueError):
        validate_job_output(directory, task)


def test_no_arbitrary_physical_cut(products):
    directory, task = products
    with h5py.File(directory / "output/img_data_0.h5", "r+") as h:
        h["I2.300000e+11"][...] *= -1
    spec = directory / "output/spectrum_0_163.00.dat"
    rows = np.loadtxt(spec, ndmin=2)
    rows[0, 1] *= -1
    np.savetxt(spec, rows, fmt="%.15e")
    metadata = validate_job_output(directory, task)
    assert metadata["diagnostics"]["negative_I_pixels"] == 400
    assert metadata["scientific_validation"] is False


def test_external_hdf_link_and_output_symlink_rejected(products, tmp_path):
    directory, task = products
    with h5py.File(directory / "output/img_data_0.h5", "r+") as h:
        del h["alpha"]
        h["alpha"] = h5py.ExternalLink("elsewhere.h5", "alpha")
    with pytest.raises(ValueError, match="link"):
        validate_job_output(directory, task)
    spec = directory / "output/spectrum_0_163.00.dat"
    spec.rename(tmp_path / "external.dat")
    spec.symlink_to(tmp_path / "external.dat")
    with pytest.raises(ValueError, match="unsafe"):
        validate_job_output(directory, task)


def test_constants_read_from_effective_header(tmp_path):
    header = tmp_path / "definitions.h"
    header.write_text("".join(f"#define {k} ({v})\n" for k, v in CORE_CONSTANTS.items()))
    assert read_core_constants(header) == CORE_CONSTANTS
    header.write_text("#define GGRAV complicated_expression\n")
    with pytest.raises(ValueError, match="constant"):
        read_core_constants(header)


@pytest.mark.parametrize("resolution", [0, 11, True, 20.0])
def test_mapping_rejects_unsupported_resolution(resolution):
    with pytest.raises(ValueError):
        hdf_mapping(resolution, 230e9)


def test_amr_rejected(products):
    directory, task = products
    task = deepcopy(task)
    task["parameters"]["MAX_LEVEL"] = 2
    with pytest.raises(ValueError, match="MAX_LEVEL"):
        validate_job_output(directory, task)
