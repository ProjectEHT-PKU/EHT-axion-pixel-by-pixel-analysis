import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.integrate
import scipy.stats as stats
from scipy.ndimage import gaussian_filter
from copy import deepcopy
from scipy import integrate, optimize
from statsmodels.stats.diagnostic import lilliefors


class CircleCalc:
    @staticmethod
    def normalize_params(lower_bound, upper_bound):
        # a1,a2 通过c*a+d变成lower~upper的范围
        c = 2 * np.pi / (upper_bound - lower_bound)
        d = -(lower_bound + upper_bound) / (upper_bound - lower_bound) * np.pi
        return c, d

    @staticmethod
    def _normalize(a):
        return (a + np.pi) % (2 * np.pi) - np.pi

    @staticmethod
    def normalize(a, lower_bound=-np.pi, upper_bound=np.pi):
        c, d = CircleCalc.normalize_params(lower_bound, upper_bound)
        return (CircleCalc._normalize(c * a + d) - d) / c

    @staticmethod
    def _circle_sub(a1, a2):
        s1 = CircleCalc.normalize(a1) - CircleCalc.normalize(a2)
        return CircleCalc.choose_small_abs(
            CircleCalc.choose_small_abs(s1, s1 - 2 * np.pi), s1 + 2 * np.pi
        )

    @staticmethod
    def circle_sub(a1, a2, lower_bound=-np.pi, upper_bound=np.pi):
        c, d = CircleCalc.normalize_params(lower_bound, upper_bound)
        return (CircleCalc._circle_sub(c * a1 + d, c * a2 + d) - d) / c

    @staticmethod
    def circle_sub_evpa(a1, a2):
        return CircleCalc.circle_sub(
            a1, a2, lower_bound=-np.pi / 2, upper_bound=np.pi / 2
        )
        # return CircleCalc.circle_sub(2 * a1, 2 * a2) / 2

    @staticmethod
    def choose_small_abs(a1, a2):
        x = deepcopy(a1)
        x[np.abs(x) > np.abs(a2)] = a2[np.abs(x) > np.abs(a2)]
        return x


def gaussian(x, *param):
    """
    input: x, mean, sigma, amplitude
    output: Gaussian function value at x
    """
    mean = param[0]
    sigma = param[1]
    amplitude = param[2]
    return (
        amplitude
        / np.sqrt(2.0 * np.pi * sigma * sigma)
        * np.exp(-np.power(x - mean, 2.0) / 2.0 / sigma / sigma)
    )


def plot_gaussian(data, bins="auto", mean=None, sigma=None):
    """
    可视化一维高斯数据分布。

    该函数会执行以下操作：
    1. 绘制数据的直方图，并进行归一化处理。
    2. 使用scipy.stats.norm.fit()从数据中估算出最佳的高斯分布参数（均值和标准差）。
    3. 在直方图上叠加绘制拟合出的高斯分布概率密度函数(PDF)曲线。
    4. 添加标题、标签和图例，使图形信息完整。

    参数:
    data (array-like): 一个包含数值数据的一维数组或列表。
    """
    # 1. 使用scipy拟合数据，得到最匹配的均值(mean)和标准差(sigma)
    if mean is None or sigma is None:
        mean, sigma = stats.norm.fit(data)

    xmin, xmax = mean - 5 * sigma, mean + 5 * sigma
    # 2. 绘制直方图
    # bins='auto' 会自动选择合适的柱子数量
    # density=True 是关键，它将直方图的高度归一化，使得总面积为1，从而可以与概率密度函数曲线叠加
    plt.hist(
        data, bins=bins, density=True, alpha=0.7, color="skyblue", range=(xmin, xmax)
    )

    # 3. 绘制拟合的高斯分布PDF曲线
    # 创建一个x轴的取值范围，覆盖数据的主要区域
    plt.xlim(xmin, xmax)
    x = np.linspace(xmin, xmax, 100)
    # 使用计算出的mean和sigma生成高斯分布的PDF值
    p = stats.norm.pdf(x, mean, sigma)
    plt.plot(x, p, "k", linewidth=2)

    # 4. 添加图形元素，使其更具可读性
    title = f"μ = {mean:.2f}, σ = {sigma:.2f}"
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)

    # 显示图形
    plt.show()


# def vonmises_mean(d1evpaarr):
#     """
#     These function input is a array of EVPA angles in radians from -0.5pi to 0.5pi.
#     """
#     kappa_fit, mu_fit, _ = stats.vonmises.fit(d1evpaarr * 2)  # 只拟合 kappa 和 mu
#     return mu_fit / 2


def move_angle_to_center(d1anglearr, lower_bound=-np.pi, upper_bound=np.pi):
    c, d = CircleCalc.normalize_params(lower_bound, upper_bound)
    _, mu_fit, _ = stats.vonmises.fit(c * d1anglearr + d)
    mu_fit = (mu_fit - d) / c
    return (
        CircleCalc.circle_sub(
            d1anglearr, mu_fit, lower_bound=lower_bound, upper_bound=upper_bound
        )
        + mu_fit
    )


def move_evpa_to_center(d1evpaarr):
    return move_angle_to_center(
        d1evpaarr, lower_bound=-np.pi / 2, upper_bound=np.pi / 2
    )


def fit_gaussian(
    d1arr,
    rm_nan=True,
    minimum_effective_number=40,
    plot_failed=True,
    angle=False,
    lower_bound=-np.pi,
    upper_bound=np.pi,
    p_value_threshold=0.002,
    **kwargs,
):
    """
    d1arr: 1D array
    evpa_angle: if True, the input array is EVPA angles in radians from -0.5pi to 0.5pi.
    p_value_threshold: threshold for the Lilliefors test to accept the null hypothesis that the data comes from a normal distribution.
    return: mean, sigma of the fitted Gaussian distribution, or np.nan, np.nan if the fit is not good.
    """
    if rm_nan:
        _d1arr = d1arr[np.isfinite(d1arr)]
    else:
        _d1arr = d1arr
    if len(_d1arr) < minimum_effective_number:
        print(f"number left less than minimum effective number: {len(_d1arr)}")
        return np.array([np.nan, np.nan])
    if angle:
        _d1arr = move_angle_to_center(
            _d1arr, lower_bound=lower_bound, upper_bound=upper_bound
        )
    else:
        pass
    _, p_value = lilliefors(_d1arr)
    if p_value > p_value_threshold:
        return np.asarray(stats.norm.fit(_d1arr))
    else:
        if plot_failed:
            plot_gaussian(_d1arr)
        return np.array([np.nan, np.nan])


def fit_gaussian2(
    d1arr,
    rm_nan=True,
    minimum_effective_number=40,
    plot_failed=True,
    angle=False,
    lower_bound=-np.pi,
    upper_bound=np.pi,
    bins=40,
    threshold=1,
    plot=False,
    nomean=False,
    **kwargs,
):
    if rm_nan:
        _d1arr = d1arr[np.isfinite(d1arr)]
    else:
        _d1arr = d1arr
    if len(_d1arr) < minimum_effective_number:
        print(f"number left less than minimum effective number: {len(_d1arr)}")
        return np.array([np.nan, np.nan])
    if angle:
        _d1arr = move_angle_to_center(
            _d1arr, lower_bound=lower_bound, upper_bound=upper_bound
        )
    else:
        pass
    _mean, _sigma = np.asarray(stats.norm.fit(_d1arr))
    counts, bins = np.histogram(
        _d1arr, bins=bins, density=True, range=(_mean - 5 * _sigma, _mean + 5 * _sigma)
    )
    xarray = np.array([(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)])
    try:
        if nomean:

            def _gaussian_noMean(x, sigma):
                return gaussian(x, _mean, sigma, 1)

            popt, pcov = optimize.curve_fit(
                _gaussian_noMean,
                xarray,
                counts,
                p0=[np.std(_d1arr)],
                maxfev=50000,
                # gtol=1e-10,
            )
            pcovsigma = pcov[0][0]
            poptsigma = popt[0]
        else:

            def _gaussian_noAmp(x, mean, sigma):
                return gaussian(x, mean, sigma, 1)

            popt, pcov = optimize.curve_fit(
                _gaussian_noAmp,
                xarray,
                counts,
                p0=[np.mean(_d1arr), np.std(_d1arr)],
                maxfev=50000,
                # gtol=1e-10,
            )
            pcovsigma = pcov[1][1]
            poptsigma = popt[1]
    except Exception as e:
        print(f"badfit, error", e)
        if plot_failed:
            plot_gaussian(_d1arr, bins=bins)
        return np.array([np.nan, np.nan])

    if np.sqrt(abs(pcovsigma)) > abs(poptsigma) * threshold:
        print(f"badfit, {np.sqrt(abs(pcovsigma))}, {abs(poptsigma)} * {threshold}")
        if plot_failed:
            plot_gaussian(
                _d1arr,
                bins=bins,
                mean=_mean if nomean else popt[0],
                sigma=abs(poptsigma),
            )
        return np.array([np.nan, np.nan])
    else:
        if plot:
            plot_gaussian(
                _d1arr,
                bins=bins,
                mean=_mean if nomean else popt[0],
                sigma=abs(poptsigma),
            )
        return np.array([_mean if nomean else popt[0], np.abs(poptsigma)])


def Plot_Background(fig, ax, I, fov=32.6868):
    shape = I.shape[0]
    xgrid = np.linspace(-fov / 2, fov / 2, shape, endpoint=False)
    ygrid = np.linspace(-fov / 2, fov / 2, shape, endpoint=False)
    diff = np.diff(xgrid)[0] / 2.0
    xgrid += diff
    ygrid += diff
    xmesh, ymesh = np.meshgrid(xgrid, ygrid, indexing="ij")

    bkg = ax.pcolormesh(
        xmesh, ymesh, I / np.nanmax(I), cmap="hot", vmin=0, vmax=1, shading="auto"
    )
    ax.set_aspect("equal")
    clb = fig.colorbar(bkg, ax=ax, fraction=0.046, pad=0.05)
    clb.ax.tick_params(labelsize=45)
    clb.set_label("$I/I_{max}$", labelpad=-40, y=1.1, rotation=0, size=100)


def plot_quivers(image, P, EVPA, fov=32.6868):

    plt.ioff()
    figure, ax0 = plt.subplots(figsize=(30, 30))
    # ax0.set_xlim(-15, 15)
    ax0.set_xlabel("x[$r_g$]", size=80)
    # ax0.set_ylim(-15, 15)
    ax0.set_ylabel("y[$r_g$]", size=80)
    plt.tick_params(labelsize=45)
    # plt.xticks([-10, -5, 0, 5, 10])
    # plt.yticks([-10, -5, 0, 5, 10])
    Plot_Background(figure, ax0, image, fov=fov)

    shape = image.shape[0]
    xgrid = np.linspace(-fov / 2, fov / 2, shape, endpoint=False)
    ygrid = np.linspace(-fov / 2, fov / 2, shape, endpoint=False)
    diff = np.diff(xgrid)[0] / 2.0
    xgrid += diff
    ygrid += diff
    xmesh, ymesh = np.meshgrid(xgrid, ygrid, indexing="ij")

    # amp=np.sqrt(image_high_res[:,:,1]**2+image_high_res[:,:,2]**2)
    # max_LP=np.max(amp)
    max_LP = np.max(P)
    vx = P * np.cos(EVPA) / max_LP * 3
    vy = P * np.sin(EVPA) / max_LP * 3
    skip = int(round(shape / 32))
    ax0.quiver(
        xmesh[::skip, ::skip],
        ymesh[::skip, ::skip],
        vx[::skip, ::skip],
        vy[::skip, ::skip],
        color="green",
        headwidth=1,
        headlength=0.1,
        width=0.006,
        units="width",
        scale=16,
    )
    ax0.quiver(
        xmesh[::skip, ::skip],
        ymesh[::skip, ::skip],
        -vx[::skip, ::skip],
        -vy[::skip, ::skip],
        color="green",
        headwidth=1,
        headlength=0.1,
        width=0.006,
        units="width",
        scale=16,
    )
    plt.suptitle("$R_g$", size=80, y=0.94)
    return figure, ax0


class ImageS:

    Intensity_max = 0.00175

    def __init__(
        self,
        data,
        blur=True,
        rotation_angle=0.6 * np.pi,
        gaussian_kernel=None,
        init_fov=46.2262,
        final_fov=32.6868,
        final_shape=64,
        mask=1,
        calc_1D=True,
    ):
        # self.data = np.flip(data.swapaxes(0,1),axis=1)
        # self.data = data.swapaxes(0,1)
        # 把图片逆时针旋转的角度(pcolormesh画出来的)
        self.calc_1D = calc_1D
        self.mask = mask
        self.rotation_angle = rotation_angle
        self.gaussian_kernel = 2.16871 if gaussian_kernel is None else gaussian_kernel
        self.data = deepcopy(data)
        self.shape = data[:, :, 0].shape[0]
        self._final_shape = final_shape
        # self.res = round(self.shape*2/80.)
        self._init_fov = init_fov
        self.res = self.shape / self._init_fov
        self._final_fov = final_fov
        self._final_res = self._final_shape / self._final_fov
        self.Get_Grid()
        self.Filter(blur=blur)
        self.max_LP = np.max(self.amp)

    def Get_Grid(self):
        xgrid = np.linspace(
            -self._init_fov / 2, self._init_fov / 2, self.shape, endpoint=False
        )
        ygrid = np.linspace(
            -self._init_fov / 2, self._init_fov / 2, self.shape, endpoint=False
        )
        diff = np.diff(xgrid)[0] / 2.0
        xgrid += diff
        ygrid += diff
        self.xmesh, self.ymesh = np.meshgrid(xgrid, ygrid, indexing="ij")

    def Filter(self, blur=True):
        I0, Q0, U0, V0 = (
            self.data[:, :, 0],
            self.data[:, :, 1],
            self.data[:, :, 2],
            self.data[:, :, 3],
        )
        if blur == True:
            gausian_kernel = self.gaussian_kernel  # * 1.558722, 20uas = 2.16871
            # gausian_kernel=0 # kernel?
            self.I = gaussian_filter(I0, gausian_kernel * self.res, mode="nearest")
            self.I *= self.mask
            self.Q = gaussian_filter(Q0, gausian_kernel * self.res, mode="nearest")
            self.U = gaussian_filter(U0, gausian_kernel * self.res, mode="nearest")
            self.V = gaussian_filter(V0, gausian_kernel * self.res, mode="nearest")
            self.I0 = I0
            iQ = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.Q, kx=1, ky=1
            )
            self.old_iQ = iQ
            iU = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.U, kx=1, ky=1
            )
            self.old_iU = iU
            iI = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.I, kx=1, ky=1
            )
            iV = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.V, kx=1, ky=1
            )
            iI0 = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.I0, kx=1, ky=1
            )
            rotation_angle = self.rotation_angle
            # rotation_angle = 0
            self.shape = self._final_shape
            self.res = self._final_res
            x, y = np.meshgrid(
                np.arange(self.shape), np.arange(self.shape), indexing="ij"
            )
            _reciprocal_res = 1 / self.res
            xv = (x - (self.shape - 1) / 2.0) * _reciprocal_res * np.cos(
                rotation_angle
            ) + (y - (self.shape - 1) / 2.0) * _reciprocal_res * np.sin(rotation_angle)
            yv = -(x - (self.shape - 1) / 2.0) * _reciprocal_res * np.sin(
                rotation_angle
            ) + (y - (self.shape - 1) / 2.0) * _reciprocal_res * np.cos(rotation_angle)
            self.I = iI.ev(xv, yv)
            self.Q = iQ.ev(xv, yv) * np.cos(2 * rotation_angle) - iU.ev(
                xv, yv
            ) * np.sin(2 * rotation_angle)
            self.U = iQ.ev(xv, yv) * np.sin(2 * rotation_angle) + iU.ev(
                xv, yv
            ) * np.cos(2 * rotation_angle)
            self.V = iV.ev(xv, yv)
            self.I0 = iI0.ev(xv, yv)
            self.data = [
                idata[18:81, 18:81] for idata in self.data
            ]  # ? 有问题，不过没有什么用
            xgrid = np.linspace(
                -self._final_fov / 2, self._final_fov / 2, self.shape, endpoint=False
            )
            ygrid = np.linspace(
                -self._final_fov / 2, self._final_fov / 2, self.shape, endpoint=False
            )
            diff = np.diff(xgrid)[0] / 2.0
            xgrid += diff
            ygrid += diff
            self.xmesh, self.ymesh = np.meshgrid(xgrid, ygrid, indexing="ij")
            self.evpa = (180.0 / np.pi) * 0.5 * np.arctan2(self.U, self.Q)
            self.amp = np.sqrt(self.Q**2 + self.U**2)
            if self.calc_1D:
                self.Create_EVPA_Dataset()
        else:
            gausian_kernel = 0  # * 1.558722, 20uas = 2.16871
            # gausian_kernel=0 # kernel?
            self.I = gaussian_filter(I0, gausian_kernel * self.res)
            self.Q = gaussian_filter(Q0, gausian_kernel * self.res)
            self.U = gaussian_filter(U0, gausian_kernel * self.res)
            self.V = gaussian_filter(V0, gausian_kernel * self.res)
            self.I0 = I0
            iQ = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.Q
            )
            iU = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.U
            )
            iI = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.I
            )
            iV = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.V
            )
            iI0 = scipy.interpolate.RectBivariateSpline(
                self.xmesh[:, 0], self.ymesh[0, :], self.I0
            )
            rotation_angle = 0.6 * np.pi
            # rotation_angle = 0
            x, y = np.meshgrid(
                np.arange(self.shape), np.arange(self.shape), indexing="ij"
            )
            xv = (x - (self.shape - 1.0) / 2.0) * 63 / (
                self.shape - 1
            ) * 0.518838 * np.cos(rotation_angle) + (
                y - (self.shape - 1.0) / 2.0
            ) * 63 / (
                self.shape - 1
            ) * 0.518838 * np.sin(
                rotation_angle
            )
            yv = -(x - (self.shape - 1.0) / 2.0) * 63 / (
                self.shape - 1
            ) * 0.518838 * np.sin(rotation_angle) + (
                y - (self.shape - 1.0) / 2.0
            ) * 63 / (
                self.shape - 1
            ) * 0.518838 * np.cos(
                rotation_angle
            )
            # xv = np.array([[(x-(self.shape-1.)/2.0)*63/(self.shape-1)*0.518838*np.cos(rotation_angle)+(y-(self.shape-1.)/2.0)*63/(self.shape-1)*0.518838*np.sin(rotation_angle) for y in range(self.shape)] for x in range(self.shape)])
            # yv = np.array([[-(x-(self.shape-1.)/2.0)*63/(self.shape-1)*0.518838*np.sin(rotation_angle)+(y-(self.shape-1.)/2.0)*63/(self.shape-1)*0.518838*np.cos(rotation_angle) for y in range(self.shape)] for x in range(self.shape)])
            # yv = np.array([[(x-31.5)*0.518838*np.cos(-0.*np.pi)-(y-31.5)*0.518838*np.sin(-0.*np.pi) for x in range(64)] for y in range(64)])
            # xv = np.array([[(x-31.5)*0.518838*np.sin(-0.*np.pi)+(y-31.5)*0.518838*np.cos(-0.*np.pi) for x in range(64)] for y in range(64)])
            # print(f'x_gap:{np.max(xv)-np.min(xv)}')
            # print(f'y_gap:{np.max(yv)-np.min(yv)}')
            self.I = iI.ev(xv, yv)
            self.Q = iQ.ev(xv, yv) * np.cos(2 * rotation_angle) - iU.ev(
                xv, yv
            ) * np.sin(2 * rotation_angle)
            self.U = iQ.ev(xv, yv) * np.sin(2 * rotation_angle) + iU.ev(
                xv, yv
            ) * np.cos(2 * rotation_angle)
            self.V = iV.ev(xv, yv)
            self.I0 = iI0.ev(xv, yv)
            # self.shape = 64
            # self.data = [idata[18:81,18:81] for idata in self.data]
            xgrid = np.linspace(-16.3434, 16.3434, self.shape, endpoint=False)
            ygrid = np.linspace(-16.3434, 16.3434, self.shape, endpoint=False)
            diff = np.diff(xgrid)[0] / 2.0
            xgrid += diff
            ygrid += diff
            self.xmesh, self.ymesh = np.meshgrid(xgrid, ygrid, indexing="ij")
            self.evpa = (180.0 / np.pi) * 0.5 * np.arctan2(self.U, self.Q)
            self.amp = np.sqrt(self.Q**2 + self.U**2)
            # self.Create_EVPA_Dataset()

    def Plot_EVPA(self, ax, vmax=Intensity_max):

        scal = np.max(self.amp)
        vx = self.amp * np.cos(self.evpa * np.pi / 180.0) / scal
        vy = self.amp * np.sin(self.evpa * np.pi / 180.0) / scal

        skip = int(round(self.shape / 16))

        ax.axis("equal")
        ax.pcolormesh(
            self.xmesh,
            self.ymesh,
            self.I,
            cmap="afmhot",
            vmin=0,
            vmax=vmax,
            shading="auto",
        )
        ax.quiver(
            self.xmesh[::skip, ::skip],
            self.ymesh[::skip, ::skip],
            vx[::skip, ::skip],
            vy[::skip, ::skip],
            headwidth=1,
            headlength=1,
            width=0.005,
            color="green",
            units="width",
            scale=16,
        )

    def Plot_Quiver11(self, ax, color_in, max_LP, i, vmax=Intensity_max):
        mycmap = plt.get_cmap("gist_rainbow")
        colors = [mycmap(j) for j in np.linspace(0, 1, 17)]
        mynorm = matplotlib.colors.BoundaryNorm(
            range(0, 16, 1), ncolors=mycmap.N, clip=True
        )
        scal = np.max(self.amp)
        vx = self.amp * np.cos(self.evpa * np.pi / 180.0) / max_LP * 4
        vy = self.amp * np.sin(self.evpa * np.pi / 180.0) / max_LP * 4
        skip = int(round(self.shape / 8))
        c = i
        if i == 0:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                color=color_in,
                headwidth=1,
                headlength=0.1,
                width=0.0015,
                units="width",
                scale=16,
            )
        else:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                c,
                headwidth=1,
                headlength=0.1,
                width=0.0015,
                cmap=mycmap,
                norm=mynorm,
                units="width",
                scale=16,
            )

    def Plot_Quiver11_reverse(self, ax, color_in, max_LP, i, vmax=Intensity_max):
        mycmap = plt.get_cmap("gist_rainbow")
        colors = [mycmap(j) for j in np.linspace(0, 1, 17)]
        mynorm = matplotlib.colors.BoundaryNorm(
            range(0, 16, 1), ncolors=mycmap.N, clip=True
        )
        scal = np.max(self.amp)
        vx = -self.amp * np.cos(self.evpa * np.pi / 180.0) / max_LP * 4
        vy = -self.amp * np.sin(self.evpa * np.pi / 180.0) / max_LP * 4
        skip = int(round(self.shape / 8))
        c = i
        if i == 0:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                color=color_in,
                headwidth=1,
                headlength=0.1,
                width=0.0015,
                units="width",
                scale=16,
            )
        else:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                c,
                headwidth=1,
                headlength=0.1,
                width=0.0015,
                cmap=mycmap,
                norm=mynorm,
                units="width",
                scale=16,
            )

    def Plot_Quiver(self, ax, color_in, max_LP, i, vmax=Intensity_max):
        scal = np.max(self.amp)
        vx = self.amp * np.cos(self.evpa * np.pi / 180.0) / max_LP
        vy = self.amp * np.sin(self.evpa * np.pi / 180.0) / max_LP
        skip = int(round(self.shape / 16))
        c = i
        if i == 16:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                color="green",
                headwidth=1,
                headlength=0.1,
                width=0.006,
                units="width",
                scale=16,
            )
        else:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                color="green",
                headwidth=1,
                headlength=0.1,
                width=0.006,
                units="width",
                scale=16,
            )

    def Plot_Quiver_reverse(self, ax, color_in, max_LP, i, vmax=Intensity_max):
        mycmap = plt.get_cmap("gist_rainbow")
        colors = [mycmap(j) for j in np.linspace(0, 1, 17)]
        mynorm = matplotlib.colors.BoundaryNorm(
            range(0, 16, 1), ncolors=mycmap.N, clip=True
        )
        scal = np.max(self.amp)
        vx = -self.amp * np.cos(self.evpa * np.pi / 180.0) / max_LP * 4
        vy = -self.amp * np.sin(self.evpa * np.pi / 180.0) / max_LP * 4
        skip = int(round(self.shape / 8))
        c = i
        if i == 16:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                color=color_in,
                headwidth=1,
                headlength=0.1,
                width=0.0025,
                units="width",
                scale=16,
            )
        else:
            ax.quiver(
                self.xmesh[::skip, ::skip],
                self.ymesh[::skip, ::skip],
                vx[::skip, ::skip],
                vy[::skip, ::skip],
                c,
                headwidth=1,
                headlength=0.1,
                width=0.0015,
                cmap=mycmap,
                norm=mynorm,
                units="width",
                scale=16,
            )

    def Plot_Background(self, fig, ax, vmax=Intensity_max):
        scal = np.max(self.amp)
        vx = self.amp * np.cos(self.evpa * np.pi / 180.0) / scal
        vy = self.amp * np.sin(self.evpa * np.pi / 180.0) / scal
        skip = int(round(self.shape / 20))
        bkg = ax.pcolormesh(
            self.xmesh,
            self.ymesh,
            self.I / np.max(self.I),
            cmap="hot",
            vmin=0,
            vmax=1,
            shading="auto",
        )

        clb = fig.colorbar(bkg, ax=ax, fraction=0.046, pad=0.05)
        clb.ax.tick_params(labelsize=45)
        clb.set_label("$I/I_{max}$", labelpad=-40, y=1.1, rotation=0, size=100)

    def EVPA_1d_from_coord(self, giving_coords):
        # iQ = scipy.interpolate.RectBivariateSpline(
        #     self.xmesh[:, 0], self.ymesh[0, :], self.Q, kx=1, ky=1
        # )
        # iU = scipy.interpolate.RectBivariateSpline(
        #     self.xmesh[:, 0], self.ymesh[0, :], self.U, kx=1, ky=1
        # )
        # new_x = giving_coords[:, 0] * np.cos(self.rotation_angle) - giving_coords[
        #     :, 1
        # ] * np.sin(self.rotation_angle)
        # new_y = giving_coords[:, 0] * np.sin(self.rotation_angle) + giving_coords[
        #     :, 1
        # ] * np.cos(self.rotation_angle)
        Qi = self.old_iQ.ev(giving_coords[:, 0], giving_coords[:, 1])
        Ui = self.old_iU.ev(giving_coords[:, 0], giving_coords[:, 1])
        _tmp = (
            np.interp(
                (np.linspace(0, 2 * np.pi, 360, endpoint=False) - self.rotation_angle)
                % (2 * np.pi),
                np.linspace(0, 2 * np.pi, 360, endpoint=False),
                0.5 * np.arctan2(Ui, Qi),
            )
            + self.rotation_angle
        )
        return ((_tmp + np.pi / 2) % np.pi) - np.pi / 2

    def EVPA_1d_from_coord_2(self, giving_coords):
        iQ = scipy.interpolate.RectBivariateSpline(
            self.xmesh[:, 0], self.ymesh[0, :], self.Q, kx=1, ky=1
        )
        iU = scipy.interpolate.RectBivariateSpline(
            self.xmesh[:, 0], self.ymesh[0, :], self.U, kx=1, ky=1
        )
        new_x = giving_coords[:, 0] * np.cos(self.rotation_angle) - giving_coords[
            :, 1
        ] * np.sin(self.rotation_angle)
        new_y = giving_coords[:, 0] * np.sin(self.rotation_angle) + giving_coords[
            :, 1
        ] * np.cos(self.rotation_angle)
        Qi = iQ.ev(new_x, new_y)
        Ui = iU.ev(new_x, new_y)
        _tmp = np.interp(
            (np.linspace(0, 2 * np.pi, 360, endpoint=False) - self.rotation_angle)
            % (2 * np.pi),
            np.linspace(0, 2 * np.pi, 360, endpoint=False),
            0.5 * np.arctan2(Ui, Qi),
        )
        return ((_tmp + np.pi / 2) % np.pi) - np.pi / 2

    def Create_EVPA_Dataset(self):

        iQ = scipy.interpolate.RectBivariateSpline(
            self.xmesh[:, 0], self.ymesh[0, :], self.Q, kx=1, ky=1
        )
        iU = scipy.interpolate.RectBivariateSpline(
            self.xmesh[:, 0], self.ymesh[0, :], self.U, kx=1, ky=1
        )
        iI = scipy.interpolate.RectBivariateSpline(
            self.xmesh[:, 0], self.ymesh[0, :], self.I, kx=1, ky=1
        )

        Rres = 10
        Rmax = 16.343
        Rmin = 0.2

        Pres = 10
        Pmax = 1 * 2 * np.pi
        Pmin = 1.0 / 360.0 * 2 * np.pi

        dr = Rmin / Rres
        dphi = Pmin / Pres

        r = np.arange(np.round(Rmax / Rmin) * Rres + 1) * dr
        phi = (np.arange(np.round(Pmax / Pmin) * Pres + 1) - Pres / 2) * dphi

        window = np.ones((Rres + 1, Pres + 1))
        window[0, :] *= 0.5
        window[-1, :] *= 0.5
        window[:, 0] *= 0.5
        window[:, -1] *= 0.5

        window0 = np.ones(Pres + 1)
        window0[0] *= 0.5
        window0[-1] *= 0.5

        x_mesh = r[:, None] * np.cos(phi[None, :])
        y_mesh = r[:, None] * np.sin(phi[None, :])

        Ii = iI.ev(x_mesh, y_mesh)
        Qi = iQ.ev(x_mesh, y_mesh)
        Ui = iU.ev(x_mesh, y_mesh)

        integrand_Q = r[:, None] * Qi * Ii * dr * dphi
        integrand_U = r[:, None] * Ui * Ii * dr * dphi

        # fast integration, for each point there are 6*6 blocks to be integrated, the border overlaps.
        # fast?
        Rn = round(Rmax / Rmin)
        Pn = round(Pmax / Pmin)
        self.info = ((Rn, Rres, Rmax, Rmin), (Pn, Pres, Pmax, Pmin))

        integrated_Q = np.array(
            [
                [
                    (
                        integrand_Q[
                            Rres * i : Rres * (i + 1) + 1, Pres * j : Pres * (j + 1) + 1
                        ]
                        * window
                    ).sum()
                    for j in range(Pn)
                ]
                for i in range(Rn)
            ]
        )
        integrated_U = np.array(
            [
                [
                    (
                        integrand_U[
                            Rres * i : Rres * (i + 1) + 1, Pres * j : Pres * (j + 1) + 1
                        ]
                        * window
                    ).sum()
                    for j in range(Pn)
                ]
                for i in range(Rn)
            ]
        )

        integrated_Q_phi = np.array(
            [
                (integrand_Q[:, Pres * j : Pres * (j + 1) + 1] * window0).sum()
                for j in range(Pn)
            ]
        )
        integrated_U_phi = np.array(
            [
                (integrand_U[:, Pres * j : Pres * (j + 1) + 1] * window0).sum()
                for j in range(Pn)
            ]
        )

        integrated_Q_phi_10deg = np.array(
            [
                np.sum([integrated_Q_phi[(i + k - 1) % Pn] for k in range(3)])
                for i in range(Pn)
            ]
        )
        integrated_U_phi_10deg = np.array(
            [
                np.sum([integrated_U_phi[(i + k - 1) % Pn] for k in range(3)])
                for i in range(Pn)
            ]
        )
        # 这里的evpa都是从pcolormesh显示的右边（x轴正向）逆时针算起的
        self.evpa_rp = (180.0 / np.pi) * 0.5 * np.arctan2(integrated_U, integrated_Q)
        self.evpa_phi = (
            (180.0 / np.pi) * 0.5 * np.arctan2(integrated_U_phi, integrated_Q_phi)
        )
        self.evpa_phi_10deg = (
            (180.0 / np.pi)
            * 0.5
            * np.arctan2(integrated_U_phi_10deg, integrated_Q_phi_10deg)
        )

        self.r = r[::Rres]

        low = int(np.round(Pres / 2))
        high = int(np.ceil(Pres / 2))
        self.phi = 0.5 * (phi[low::Pres] + phi[high::Pres])

        """

        iIQ_polar = lambda r,phi : r*iQ(r*np.cos(phi), r*np.sin(phi))*iI(r*np.cos(phi), r*np.sin(phi))
        iIU_polar = lambda r,phi : r*iU(r*np.cos(phi), r*np.sin(phi))*iI(r*np.cos(phi), r*np.sin(phi))
       
        
         
        N=100
        dphi = 2*np.pi/N
        phis = dphi * np.arange(0,N+1)
        IQ_smear = np.zeros(N)
        IU_smear = np.zeros(N)
                
        for i in range(N):
            if type(r) in [float,int]:
                IQ_smear[i] = iIQ_polar(r,phis[i])
                IU_smear[i] = iIU_polar(r,phis[i])
            elif type(r) in [list]:
                rmin=r[0]
                rmax=r[1]
                IQ_smear[i] = scipy.integrate.dblquad(iIQ_polar,phis[i],phis[i+1],rmin,rmax)[0]
                IU_smear[i] = scipy.integrate.dblquad(iIU_polar,phis[i],phis[i+1],rmin,rmax)[0]
        self.evpa = (180./np.pi)*0.5*np.arctan2(IU_smear,IQ_smear)
        self.phis = phis[:-1]
        
        """
