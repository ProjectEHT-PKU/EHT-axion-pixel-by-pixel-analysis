import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import warnings
from scipy.integrate import quad, IntegrationWarning


def calculate_source_position(
    obs_rho: float, obs_phi: float, r_o: float, theta_o: float, a: float, M: float = 1.0
):
    """
    根据给定的观测者位置(r_o, theta_o)和守恒量(lambda, eta)，
    计算光线在克尔时空中传播到theta_s = pi/2时的源位置(r_s, phi_s)。
    此版本使用解析/多项式求根方法计算转折点，更加高效和稳健。

    参数:
    lambda_val (float): 守恒量 lambda = L/E。
    eta_val (float): 守恒量 eta = Q/E^2。
    r_o (float): 观测者的径向坐标。
    theta_o (float): 观测者的极向坐标 (弧度)，假定 theta_o < pi/2。
    a (float): 黑洞自旋参数 (0 <= a <= M)。默认为 0.99。
    M (float): 黑洞质量。默认为 1.0。

    返回:
    tuple[float, float] | None: 一个包含(r_s, phi_s)的元组，如果计算失败则返回None。
    """
    _epsilon_ = 1e-10
    # ------------------ 1. 定义势函数和基本参数 ------------------
    theta_s = np.pi / 2.0
    phi_o = 0.0
    t_o = 0.0

    if np.abs(np.sin(obs_phi)) < 0.00013:
        warnings.warn("obs_phi为0或π，可能导致计算不稳定。")
    _alpha = obs_rho * np.sin(obs_phi)
    _beta = obs_rho * np.cos(obs_phi)
    lambda_val = -_alpha * np.sin(theta_o)
    eta_val = (_alpha**2 - a**2) * np.cos(theta_o) ** 2 + _beta**2

    event_horizon_plus = M + np.sqrt(M**2 - a**2)

    def delta(r):
        return r**2 - 2 * M * r + a**2

    def R_func(r):
        term1 = (r**2 + a**2 - a * lambda_val) ** 2
        term2 = delta(r) * (eta_val + (lambda_val - a) ** 2)
        return term1 - term2

    def Theta_func(theta):
        return eta_val + a**2 * np.cos(theta) ** 2 - lambda_val**2 / np.tan(theta) ** 2

    def get_theta_c():
        # ------------------ 2. 解析/多项式方法计算转折点 ------------------

        # --- 极向转折点 theta_c ---
        # Θ(θ)=0 可化为 u=cos²θ 的一元二次方程: A*u² + B*u + C = 0
        theta_c = None
        A_u = a**2
        B_u = lambda_val**2 + eta_val - a**2
        C_u = -eta_val

        discriminant = B_u**2 - 4 * A_u * C_u
        if discriminant >= 0:
            # 解出 u = cos²θ
            u1 = (-B_u + np.sqrt(discriminant)) / (2 * A_u)
            u2 = (-B_u - np.sqrt(discriminant)) / (2 * A_u)
            # 选取物理上允许的解 (0 <= u <= 1)，事实上，一定有一个解满足这个条件?，但可能由于精度影响导致略大于1。需要将其手动修正为1。
            # 如果有解位于theta_o和theta_s之间，则说明光子无法到达theta_s
            # 正常的转折点应该是比theta_o更小的，并且要选其中最接近theta_o的值
            for u_sol in [u1, u2]:
                if 1 < u_sol < 1 + 1e-10:  # 如果解略大于1，则将其修正为1
                    warnings.warn(f"u_sol={u_sol}略大于1，已修正为1")
                    u_sol = 1.0
                if 0 <= u_sol <= 1:
                    _theta_c = np.arccos(np.sqrt(u_sol))
                    if theta_o < _theta_c < theta_s:
                        raise ValueError("光子无法到达theta_s")
                    elif _theta_c <= theta_o:
                        if theta_c is None or _theta_c > theta_c:
                            theta_c = _theta_c

        return theta_c

    def get_r_c():

        # --- 径向转折点 r_c ---
        # R(r)=0 是一个四次多项式: r⁴ + c₃r³ + c₂r² + c₁r + c₀ = 0
        r_c = None
        c4 = 1.0
        c3 = 0.0
        c2 = a**2 - eta_val - lambda_val**2
        c1 = 2 * M * (eta_val + (lambda_val - a) ** 2)
        c0 = -(a**2) * eta_val

        # 使用numpy.roots求解多项式的全部根
        poly_coeffs = [c4, c3, c2, c1, c0]
        all_roots = np.roots(poly_coeffs)

        # 筛选出在 (r_+, r_o) 之间的实数根
        real_roots = [r.real for r in all_roots if np.isreal(r)]
        physical_roots = [r for r in real_roots if event_horizon_plus < r < r_o]

        if physical_roots:
            r_c = max(physical_roots)  # 近心点是最大的物理根

        return r_c

    theta_c = get_theta_c()  # 这个看起来应该一定会有一个值
    if theta_c:
        # warnings.warn("theta_c已略微移动以避免奇点")
        theta_c += _epsilon_
    r_c = get_r_c()  # 这个不一定
    if r_c:
        # warnings.warn("r_c已略微移动以避免奇点")
        r_c += _epsilon_
    # print(f"theta_c: {theta_c}, r_c: {r_c}")

    if theta_c is None:
        print(f"r_c: {r_c}, theta_c: {theta_c}")
        raise ValueError("转折点是None？")

    def sqrt_R(r):
        R_val = R_func(r)
        return np.sqrt(R_val) if R_val >= 0 else np.sqrt(R_val + 0j)

    def sqrt_Theta(theta):
        Theta_val = Theta_func(theta)
        return np.sqrt(Theta_val) if Theta_val >= 0 else np.sqrt(Theta_val + 0j)

    # ------------------ 3. 计算路径长度积分 G_theta ------------------
    def integrand_G_theta(theta, sign_theta):
        return 1.0 / sqrt_Theta(theta) * sign_theta

    def integrand_I_r(r, sign_r):
        return 1.0 / sqrt_R(r) * sign_r

    def integrand_G_phi(theta, sign_theta):
        return 1.0 / np.sin(theta) ** 2 / sqrt_Theta(theta) * sign_theta

    def integrand_I_phi(r, sign_r):
        return a * (2 * M * r - a * lambda_val) / delta(r) / sqrt_R(r) * sign_r

    def integrand_G_t(theta, sign_theta):
        return np.cos(theta) ** 2 / sqrt_Theta(theta) * sign_theta

    def integrand_I_t(r, sign_r):
        return (
            (r**2 * delta(r) + 2 * M * r * (r**2 + a**2 - a * lambda_val))
            / delta(r)
            / sqrt_R(r)
            * sign_r
        )

    if _beta > 0:  # theta会先变小再变大(从观测者到源来看)
        sign_theta = -1.0
        G_theta_1, _ = quad(integrand_G_theta, theta_o, theta_c, sign_theta)
        G_theta_2, _ = quad(integrand_G_theta, theta_c, theta_s, -sign_theta)
        G_theta = G_theta_1 + G_theta_2
        G_phi_1, _ = quad(integrand_G_phi, theta_o, theta_c, sign_theta)
        G_phi_2, _ = quad(integrand_G_phi, theta_c, theta_s, -sign_theta)
        G_phi = G_phi_1 + G_phi_2
        G_t_1, _ = quad(integrand_G_t, theta_o, theta_c, sign_theta)
        G_t_2, _ = quad(integrand_G_t, theta_c, theta_s, -sign_theta)
        G_t = G_t_1 + G_t_2

    else:  # theta会单调变大(从观测者到源来看)
        sign_theta = 1.0
        sign_r = -1.0
        G_theta, _ = quad(integrand_G_theta, theta_o, theta_s, sign_theta)
        G_phi, _ = quad(integrand_G_phi, theta_o, theta_s, sign_theta)
        G_t, _ = quad(integrand_G_t, theta_o, theta_s, sign_theta)

    sign_r = -1.0
    if (
        r_c is not None
        and (_I_r_test := quad(integrand_I_r, r_o, r_c, sign_r)[0]) < G_theta
    ):
        # r方向经过转折点
        _I_r = (
            lambda _r_s: quad(integrand_I_r, r_c, _r_s, -sign_r)[0]
            - G_theta
            + _I_r_test
        )
        r_s = brentq(_I_r, r_c, r_o)
        I_phi_1, _ = quad(integrand_I_phi, r_o, r_c, sign_r)
        I_phi_2, _ = quad(integrand_I_phi, r_c, r_s, -sign_r)
        I_phi = I_phi_1 + I_phi_2
        I_t_1, _ = quad(integrand_I_t, r_o, r_c, sign_r)
        I_t_2, _ = quad(integrand_I_t, r_c, r_s, -sign_r)
        I_t = I_t_1 + I_t_2
    else:
        # r方向没有经过转折点
        _I_r = lambda _r_s: quad(integrand_I_r, r_o, _r_s, sign_r)[0] - G_theta
        _bound = r_c or event_horizon_plus
        if _I_r(r_o) * _I_r(_bound) > 0:
            raise ValueError("光线无法到达theta_s，可能是因为参数设置不当。")
        r_s = brentq(_I_r, r_o, _bound)
        I_phi, _ = quad(integrand_I_phi, r_o, r_s, sign_r)
        I_t, _ = quad(integrand_I_t, r_o, r_s, sign_r)

    phi_s = phi_o + float(phi_o + I_phi + lambda_val * G_phi)
    t_s = t_o + float(I_t + a**2 * G_t)

    return r_s, phi_s, t_s
