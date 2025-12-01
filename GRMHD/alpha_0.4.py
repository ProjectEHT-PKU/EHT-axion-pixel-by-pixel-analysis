RERUN_ANALYSIS = True

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
from copy import deepcopy
import os

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), "..")))
from utils.libanaly import Plot_Background, plot_quivers, CircleCalc, ImageS
from utils.raytracing import calculate_source_position
from scipy.optimize import curve_fit


def getGRMHDDataAfterImageS(
    grmhdDataOrigin,
    rotation_angle=0.6 * np.pi,
    gaussian_kernel=None,
    init_fov=46.2262,
    final_shape=64,
    image_mask=1,
    oneD_coord=None,
    oneD_method=1,
    final_fov=32.6868,
):
    """
    ImageS blurs the data image, and calculates many useful quantities. This function calculates all the data at once.
    """
    shape = grmhdDataOrigin.shape
    shape = (shape[0], shape[1], final_shape, final_shape, 4)
    grmhdDataAfterImageS = np.zeros(shape)
    EVPA = np.zeros((shape[0], shape[1], final_shape, final_shape))
    EVPA_1d = np.zeros((shape[0], shape[1], 360))
    EVPA_1d_10deg = np.zeros((shape[0], shape[1], 360))
    for t in range(shape[0]):
        print(t)
        for i in range(shape[1]):
            tmp = ImageS(
                grmhdDataOrigin[t, i],
                rotation_angle=rotation_angle,
                gaussian_kernel=gaussian_kernel,
                init_fov=init_fov,
                final_shape=final_shape,
                mask=image_mask,
                final_fov=final_fov,
                calc_1D=(oneD_coord is None),
            )
            grmhdDataAfterImageS[t, i, :, :, 0] = tmp.I
            grmhdDataAfterImageS[t, i, :, :, 1] = tmp.Q
            grmhdDataAfterImageS[t, i, :, :, 2] = tmp.U
            grmhdDataAfterImageS[t, i, :, :, 3] = tmp.V
            EVPA[t, i] = tmp.evpa / 180 * np.pi  # rad
            if oneD_coord is None:
                EVPA_1d[t, i] = np.roll(
                    tmp.evpa_phi / 180 * np.pi, -90
                )  # rad 从pcolormesh显示的顶端（y轴正向）逆时针算起的
                EVPA_1d_10deg[t, i] = np.roll(
                    tmp.evpa_phi_10deg / 180 * np.pi, -90
                )  # rad 从pcolormesh显示的顶端（y轴正向）逆时针算起的
            else:
                if oneD_method == 1:
                    EVPA_1d[t, i] = tmp.EVPA_1d_from_coord(
                        oneD_coord
                    )  # 从pcolormesh显示的顶端（y轴正向）逆时针算起的
                    EVPA_1d_10deg[t, i] = EVPA_1d[t, i]
                elif oneD_method == 2:
                    EVPA_1d[t, i] = tmp.EVPA_1d_from_coord_2(
                        oneD_coord
                    )  # 从pcolormesh显示的顶端（y轴正向）逆时针算起的
                    EVPA_1d_10deg[t, i] = EVPA_1d[t, i]

    return grmhdDataAfterImageS, EVPA, EVPA_1d, EVPA_1d_10deg


if RERUN_ANALYSIS:
    for alpha in [0.4]:
        _grmhd_data_origin = np.load(
            f"data_in/GRMHD_Data_1600-3000_alpha={alpha}_res=100.npy"
        )
        _grmhdDataAfterImageS, _EVPA, _EVPA_1d, _EVPA_1d_10deg = (
            getGRMHDDataAfterImageS(_grmhd_data_origin)
        )
        np.save(
            f"data_out/grmhdDataAfterImageS_alpha={alpha}_res=100.npy",
            _grmhdDataAfterImageS,
        )
        np.save(f"data_out/EVPA_alpha={alpha}_res=100.npy", _EVPA)
        np.save(f"data_out/EVPA_1d_alpha={alpha}_res=100.npy", _EVPA_1d)
        np.save(f"data_out/EVPA_1d_10deg_alpha={alpha}_res=100.npy", _EVPA_1d_10deg)
