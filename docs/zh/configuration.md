# 配置参考

[← 返回文档索引](index.md)

本仓库将配置严格划分为两类：写在 [`configs/experiment.yaml`](../../configs/experiment.yaml) 中的科学参数，以及通过 `.env` 和命令行选项给出的本机路径与并发执行设置（三个路径根目录：`EHT_RAPTOR_SOURCE_ROOT`、`EHT_INPUT_ROOT`、`EHT_DATA_ROOT`）。

## 参数文件结构与校验

`configs/experiment.yaml` 是全仓库唯一的科学参数文件（所有 `--config` 选项的默认值）。每个脚本都通过 `eht_axion/experiment.py` 中的 `read_experiment()` 加载并校验该文件：

```bash
pixi run python -c "from eht_axion.experiment import read_experiment; read_experiment()"
```

若通过 `--config PATH` 指定替代 YAML 文件，该文件必须包含完全相同的四个顶层节及嵌套键集合：

- `shared`：黑洞物理尺度、$\alpha$ 网格、束缚态量子数，以及公共的 $64 \times 64$ 天空网格。
- `raptor`：`RAPTOR` 源码仓库与提交号、相机几何、编译开关、等离子体模型、快照范围、图像分辨率、径向表映射和 BHAC 网格坐标。
- `eht_imaging`：公开发布的 HOPS `uvfits` 文件身份、校准容差、成像权重、基准交替轮次、凸优化精化预算以及仪器偏振泄漏（D 项，D-term）蒙特卡罗（Monte Carlo）设置。
- `analysis`：天空旋转角、扇形孔径几何、特征选择阈值、超辐射（superradiant）轴子云标度、广义相对论磁流体力学（GRMHD）背景滞后帧数以及贝叶斯推断（Bayesian inference）的网格与先验设置。

**校验规则（`eht_axion/experiment.py`）：**
- 使用自定义的 `SafeLoader`（`_UniqueKeyLoader`），拒绝任何重复的 YAML 键（包括合并键）。
- 强制检查每个字典节的精确键集合（`SCHEMA`）、严格的 Python 类型（要求 `int` 或 `float` 时拒绝 `bool`；整数必须为 `int`）、有限性及定义域上下界；不填充任何默认值。
- 检查两项跨节几何约束：`shared.image.fov_uas` 换算为引力半径 $r_g$ 后（$31.42\,r_g$）不得超过 `raptor.camera_fov_rg`（$46.2261\,r_g$）；且 `analysis.features.min_radius_rg`（$2.0\,r_g$）必须严格小于孔径外半径（`radial_max_fraction * final_fov / 2`，即 $15.71\,r_g$）。

### 派生视图与模块常量

`eht_axion/experiment.py` 提供四个写入阶段清单（manifest）和任务计划的派生视图，以及两个固定代码常量：

| 符号 | 返回结果 | 调用方 |
|---|---|---|
| `physical_scale(config)` | `PhysicalScale`（`tg_seconds` $= 32\,015.7\,\mathrm{s}$、`tg_hours` $= 8.893\,\mathrm{h}$、`angular_rg_uas` $= 3.819\,\mu\mathrm{as}$ 及单位换算方法） | 各分析脚本、`RAPTOR` 帧时间间隔（`input_cadence_hours`） |
| `ray_tracing_config(config)` | `shared` 中的物理标量（`mass_solar`、`distance_mpc`、`label_step_tg`、`alphas`、`principal_n`）加上不含 `raptor.source` 的 `raptor` 节 | `RAPTOR` 任务计划、`generated.json` 图像数据立方身份记录 |
| `analysis_geometry(config)` | `image_preprocess`（以 $r_g$ 和 $\mathrm{rad}$ 为单位的 `sigma`、`init_fov`、`final_fov`、`final_shape`、`rotation_angle`）、`aperturefit`（`fov`、`radial_max`、`window_degrees`、`angular_bins`）和 `metadata` | `scripts/analysis/expected_signal.py`、`scripts/analysis/limits.py` |
| `eht_imaging_parameters(config, locations)` | 由 `EHT_LOCATIONS`、`shared.image`（`npix`、`fov_uas`、`blur_fwhm_uas`）以及不含 `refine` 和 `monte_carlo` 的 `eht_imaging` 键组成的扁平字典 | `scripts/eht/reconstruct.py`、`scripts/eht/refine.py`、`scripts/eht/export_bundle.py` |
| `EHT_DAYS = (95, 96, 100, 101)` | 2017 年观测年积日（`eht_imaging.inputs` 必须精确匹配该集合） | 全部 EHT 成像、汇总与分析模块 |
| `ANGULAR_BINS = 360` | 1D 扇形孔径观测量中 $1^\circ$ 间隔的扇区分格数 | `analysis_geometry`、`scripts/analysis/limits.py` |

## 参数详表：`configs/experiment.yaml`

### `shared`

| 键 | 当前取值 | 单位 | 物理含义与校验约束 | 读取位置 |
|---|---|---|---|---|
| `mass_solar` | `6500000000.0` | $M_\odot$ | M87\* 质量（$> 0$） | `physical_scale`；`RAPTOR` `MBH` |
| `distance_mpc` | `16.8` | $\mathrm{Mpc}$ | 到 M87\* 的距离（$> 0$；乘以 $1000$ 后以 $\mathrm{kpc}$ 写入 `RAPTOR` `DISTANCE`） | `physical_scale`；`RAPTOR` `DISTANCE` |
| `label_step_tg` | `0.1` | $t_g$ / 标签 | 每个 GRMHD 快照标签单位对应的预期尾部时间步长（$> 0$） | `plan_generation` 尾部时间间隔检查；`PhysicalScale` |
| `alphas` | `[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]` | – | 严格递增的两位小数正值 $\alpha$ 网格（共 8 个值） | `RAPTOR` 规划与图像数据立方、预期信号、似然、论文脚本 |
| `principal_n` | `2` | – | $\Omega(\alpha) = \alpha(1 - \alpha^2/(2n^2))$ 中的主量子数 $n \ge 1$ | `RAPTOR` `AXION_OMEGA`；`analysis/model.py` |
| `image.shape` | `64` | 像素 | 正方形天空网格的边长（$\ge 2$） | EHT `npix`、系综汇总、分析重采样 |
| `image.fov_uas` | `120.0` | $\mu\mathrm{as}$ | 天空网格的视场（field of view, FOV，每像素 $1.875\,\mu\mathrm{as}$；替换官方 EHT 脚本中的 `128` $\mu\mathrm{as}$） | EHT `fov_uas`、分析重采样与孔径 |
| `image.blur_fwhm_uas` | `19.0` | $\mu\mathrm{as}$ | 同时施加于 EHT 模糊图像与模拟图像的圆高斯匹配核半高全宽（full width at half maximum, FWHM，$\ge 0$；非等效点扩散函数，PSF） | EHT `save_image`；分析预处理的 $\sigma$ |

### `raptor`

| 键 | 当前取值 | 单位 | 物理含义与校验约束 | 读取位置 |
|---|---|---|---|---|
| `source.repository` | `github.com/projecteht-pku/raptor` | – | 规范化小写 `host/owner/name`；本地工作副本至少有一个 Git 远程地址规范化后等于此值 | `build_raptor`、`generate.py` 计划的 `raptor_source`、`run_jobs` |
| `source.commit` | `20056006820b6536c96de790dd5fdaad5ada4eae` | – | 洁净本地 `RAPTOR` 工作副本 `HEAD` 必须精确匹配的 40 位十六进制提交号 | `build_raptor`、`run_jobs`、`write_images`、`paper.py`（`src/axion.c`） |
| `camera_fov_rg` | `46.2261` | $r_g$ | 正方形相机视场宽度（$176.5\,\mu\mathrm{as}$） | `RAPTOR` `CAM_SIZE_X/Y`；`analysis_geometry` `init_fov` |
| `fluid_spin` | `0.9375` | – | 无量纲黑洞自旋 $a \in [0, 1]$，与每个 BHAC 快照尾部核对 | `read_snapshot_footer`、`paper.py` |
| `build.model` | `bhac` | – | `setup.sh --code`（`bhac`、`harm3d`、`hamr`；`run_jobs` 要求 `MKSBHAC`） | `build_raptor` |
| `build.metric` | `mks` | – | `setup.sh --metric`（`mks`、`cks`） | `build_raptor` |
| `build.integrator` | `rk4` | – | `setup.sh --int`（`rk4`、`rk2`、`rk45`、`ver`） | `build_raptor` |
| `build.grid` | `SMR` | – | `setup.sh --grid`（`SMR`、`AMR`） | `build_raptor` |
| `build.polarized` | `true` | – | `setup.sh --rad=pol`（`run_jobs` 要求为 `true`） | `build_raptor` |
| `build.sfc` | `true` | – | `setup.sh --sf=sfc` | `build_raptor` |
| `build.magnetic_polarity` | `plus` | – | `setup.sh --bflip`（`plus`、`minus`） | `build_raptor` |
| `build.sigma_cut` | `1.0` | – | 磁化强度截断 `setup.sh --sigma-cut`（$> 0$） | `build_raptor` |
| `build.thetae_min` | `0.001` | – | 无量纲电子温度下限 `--thetae-min`（$> 0$，且 $< \text{thetae\_max}$） | `build_raptor` |
| `build.thetae_max` | `100.0` | – | 无量纲电子温度上限 `--thetae-max` | `build_raptor` |
| `build.uniform_output` | `false` | – | `setup.sh --uniform=0`（`run_jobs` 要求为 `false`） | `build_raptor` |
| `plasma.mass_unit_g` | `3.0e+25` | $\mathrm{g}$ | 吸积流质量标度 `M_UNIT`（$> 0$） | `RAPTOR` `model.in` |
| `plasma.r_low` | `1.0` | – | 离子–电子温度比参数 `R_LOW`（$> 0$） | `RAPTOR` `model.in` |
| `plasma.r_high` | `1.0` | – | 离子–电子温度比参数 `R_HIGH`（$> 0$） | `RAPTOR` `model.in` |
| `plasma.inclination_degrees` | `163.0` | $\mathrm{deg}$ | 观测者倾角 `INCLINATION`（$\in [0, 180]$） | `RAPTOR` `model.in`、谱文件名 |
| `plasma.frequency_hz` | `230.0e+9` | $\mathrm{Hz}$ | 单一观测频率 `FREQ_MIN`（$> 0$） | `RAPTOR` `model.in`、HDF5 数据集名 |
| `plasma.step_size` | `0.005` | – | 测地线积分步长 `STEPSIZE`（$> 0$） | `RAPTOR` `model.in` |
| `snapshots.{start, stop, step}` | `1600, 3000, 10` | 标签 | 含端点的快照标签等差网格（共 $141$ 个快照，`0 <= start <= stop <= 999999`） | `plan_generation`、`write_images`、`generated_inputs` |
| `response.resolution` | `100` | 像素 | 有轴子图像与 100 像素基线图像的正方形边长（10 的整数倍） | `plan_generation`、`write_images`、`expected_signal.py` |
| `response.reference_norm` | `0.1` | – | `AXION_NORM` 参考场振幅（$> 0$） | `RAPTOR` `model.in`、`limits.py`、`paper.py` |
| `background.resolution` | `400` | 像素 | 无轴子背景图像的正方形边长（10 的整数倍） | `plan_generation`、`write_images`、`expected_signal.py`、`paper.py` |
| `radial_tables.low.spin` | `0.9375` | – | `low` 径向 CSV 表声明的黑洞自旋 $a \in [0, 1]$ | `plan_generation`、`paper.py` |
| `radial_tables.high.spin` | `0.99` | – | `high` 径向 CSV 表声明的黑洞自旋 $a \in [0, 1]$（有意与 `fluid_spin` 不同） | `plan_generation`、`paper.py` |
| `alpha_table` | `0.10`–`0.30`: `low`；`0.35`–`0.44`: `high` | – | 将全部两位小数字符串键 `'0.10'`、…、`'0.44'` 映射到 `low` 或 `high` | `plan_generation`、`run_jobs`、`paper.py` |
| `grid.nxlone` | `[96, 48, 48]` | 网格数 | BHAC `grid_mks.in` 基础计算域网格数（快照块大小 `nx` 必须整除 `nxlone`） | `read_snapshot_footer`、`run_jobs` |
| `grid.xprobmin` | `[0.17395330712343798, 0.0, 0.0]` | MKS 坐标 | `grid_mks.in` 中的 `xprobminx1..3`（乘以 $2\pi$ 之前的原始 BHAC 角度坐标） | `run_jobs` |
| `grid.xprobmax` | `[7.824046010856292, 0.5, 1.0]` | MKS 坐标 | `grid_mks.in` 中的 `xprobmax1..3` | `run_jobs` |
| `grid.hslope` | `0.0` | – | `grid_mks.in` 中的 `hslope` | `run_jobs` |

### `eht_imaging`

注意：总强度重建所用的致密源零基线流量密度（`zbl = 0.60` Jy）固定在官方 EHT 成像脚本内部，不在 `configs/experiment.yaml` 中；偏振阶段使用阶段 `I` 图像的总流量 `fixed.total_flux()`。

| 键 | 当前取值 | 单位 | 物理含义与校验约束 | 读取位置 |
|---|---|---|---|---|
| `inputs` | `{95: [...], 96: [...], 100: [...], 101: [...]}` | – | 每天 HOPS 文件 `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits` 的 `[filename, sha256]` | `reconstruct.py`（核对 `$EHT_INPUT_ROOT/eht/` 中文件字节的 SHA-256） |
| `gain_tol` | `[0.02, 0.2]` | – | 共享增益自校准 `selfcal` 的 `[lower, upper]` 相对增益容差 | `reconstruct.py`（`--stage I`） |
| `fsys` | `0.005` | – | 按平方和加到 `RL`/`LR` 误差上的相对系统噪声（$> 0$） | `reconstruct.py`、`refine.py`、`sample.py` 中的 `add_fsys` |
| `solve_sites` | `['AZ', 'LM', 'PV']` | – | 求解 D 项的台站列表（与 `fixed_sites` 互不相交） | `reconstruct.py`、`refine.py`、`sample.py` |
| `fixed_sites` | `['AA', 'AP', 'SM', 'JC']` | – | 已在公开发布的 HOPS `uvfits` 中预先完成泄漏改正的台站 | `load_native` 元数据、`sample_dterms` |
| `feed_rotation` | `{AA: [1, 0, 0], …, SM: [1, -1, 45]}` | –, –, deg | 每个台站的场旋转系数 `[fr_par, fr_elev, fr_off]`（EHT Paper VII 表 1） | `load_native`、`load_gaincal_base` |
| `p_data` | `{pvis: 1.0, m: 1.0}` | – | $P$ 目标函数中的数据项权重（$> 0$） | `reconstruct.py`、`convex.py` |
| `p_reg` | `{hw: 100.0, ptv: 1.0}` | – | 正则项权重：偏振熵权重 `hw` 与偏振全变差（total variation, TV）权重 `ptv`（$> 0$） | `reconstruct.py`、`convex.py` |
| `norm_reg` | `true` | – | 启用 `ehtim` 正则项归一化（`convex.py` 要求为 `true`） | `reconstruct.py`、`convex.py` |
| `initial_m` | `0.1` | – | 基准阶段第 1 轮的初始线偏振度 $m_0 \in (0, 1)$ | `reconstruct.py` |
| `initial_evpa_deg` | `0.0` | $\mathrm{deg}$ | 基准阶段第 1 轮的初始均匀电矢量位置角（EVPA）$\chi_0$ | `reconstruct.py` |
| `rounds` | `50` | – | 基准 $P$/$D$ 交替轮数（$\ge 1$） | `reconstruct.py` |
| `restart_rounds` | `[2, 12, 22, 32, 42]` | – | `[2, rounds]` 内互不相同的轮次，从模糊后的复偏振度重新初始化 $P$ | `reconstruct.py` |
| `restart_blur_uas` | `20.0` | $\mu\mathrm{as}$ | `restart_rounds` 所用圆高斯核的 FWHM（$> 0$） | `reconstruct.py` |
| `restart_domain_epsilon` | `1.0e-08` | – | 收缩因子 $\epsilon \in (0, 1)$，重启时将 $P$ 乘以 $1 - \epsilon$ 以严格满足 $\lvert P/I\rvert < 1$ | `reconstruct.py` |
| `maxit` | `500` | 迭代步 | `fiducial` 阶段单次调用 `ehtim` 优化器的迭代上限 | `reconstruct.py` |
| `stop` | `1.0e-08` | – | `fiducial` 阶段 `ehtim` 优化器收敛容差 | `reconstruct.py` |
| `shift_objective_constant` | `true` | – | 向 `ehtim` $P$ 目标函数加上常数 $w_{\mathrm{pvis}} + w_m + h_w\ln 2$（极小点不变） | `reconstruct.py` |
| `uv_min` | `100000000.0` | $\lambda$ | 最短基线长度截断（$100\,\mathrm{M}\lambda$），在共享增益校准前及全部 $P$ 成像中剔除同址零基线 | `reconstruct.py`、`refine.py`、`sample.py` |
| `uncertainty_draws` | `1000` | 次/天 | 每天的 D 项蒙特卡罗抽样次数（四天共 $4000$ 次） | `sample.py`、`eht/summarize.py`、`limits.py`、`check_covariance.py` |
| `seed` | `20210405` | – | 基础非负整数随机数种子；第 $d$ 天使用 `SeedSequence([seed, day])` | `sample.py` |
| `covariance_ddof` | `1` | – | 四天 $(\operatorname{Re} D, \operatorname{Im} D)$ 样本协方差的自由度修正（delta degrees of freedom，$\in [0, 3]$） | `d_distribution`、`sample.py` |
| `residual_d_sigma` | `0.01` | – | `fixed_sites` 在每个实部/虚部分量（$\operatorname{Re}$、$\operatorname{Im}$）上的零均值高斯标准差 | `sample_dterms`、`sample.py` |
| `refine.outer_rounds` | `100` | 轮 | 最终 $P$ 求解前精确 $P$ 求解与前向 $D$ 重拟合的外层交替轮数 | `refine.py` |
| `refine.gap_atol` | `1.0e-06` | 目标函数 | 每次原始–对偶混合梯度（PDHG）$P$ 求解的对偶间隙（duality gap）证书绝对容差 | `refine.py` |
| `refine.maxiter` | `20000` | 迭代步 | `refine.py` 中单次 $P$ 求解的 PDHG 最大迭代步数 | `refine.py` |
| `refine.check_every` | `500` | 迭代步 | 计算对偶间隙下界的迭代间隔 | `refine.py` |
| `refine.time_limit` | `300.0` | $\mathrm{s}$ | `refine.py` 中单次 $P$ 求解的墙钟时间上限 | `refine.py` |
| `refine.d_optimizer.{ftol, gtol, maxiter, maxcor, maxls, maxfun}` | `1e-13, 1e-9, 10000, 50, 50, 100000` | – | 每次前向 $D$ 重拟合传给 `scipy.optimize.minimize` 的 `L-BFGS-B` 选项 | `refine.py` → `fit_forward_d` |
| `monte_carlo.initial_p_fraction` | `0.05` | – | 每次蒙特卡罗抽样冷启动的初始值 $P = 0.05\,I$（实数，$\in (0, 1)$） | `sample.py` |
| `monte_carlo.gap_atol` | `1.0e-06` | 目标函数 | 每次抽样获得数值证书的对偶间隙绝对容差 | `sample.py` |
| `monte_carlo.maxiter` | `50000` | 迭代步 | 第一次尝试的 PDHG 迭代步上限 | `sample.py` |
| `monte_carlo.check_every` | `500` | 迭代步 | 对偶间隙检查的迭代间隔 | `sample.py` |
| `monte_carlo.time_limit` | `600.0` | $\mathrm{s}$ | 第一次尝试的墙钟时间上限 | `sample.py` |
| `monte_carlo.retry_maxiter` | `200000` | 迭代步 | 第一次尝试未获证书时第二次尝试的迭代步上限 | `sample.py` |
| `monte_carlo.retry_time_limit` | `1800.0` | $\mathrm{s}$ | 第二次尝试的墙钟时间上限 | `sample.py` |

### `analysis`

| 键 | 当前取值 | 单位 | 物理含义与校验约束 | 读取位置 |
|---|---|---|---|---|
| `rotation_degrees` | `108.0` | $\mathrm{deg}$ | 将 `RAPTOR` 相机坐标系映射到天空的主动旋转角 $\theta$（同时旋转图像并将 $Q + iU$ 乘以 $e^{2i\theta}$） | `analysis_geometry` |
| `aperture.window_degrees` | `10` | $\mathrm{deg}$ | 每个扇形孔径的角向全宽 $w \in (0, 360]$ | `analysis_geometry` → `aperture_operator` |
| `aperture.radial_max_fraction` | `1.0` | – | 孔径外半径 $r_{\max}$ 占半视场（$60\,\mu\mathrm{as}$）的比例（$\in (0, 1]$） | `analysis_geometry` → `aperture_operator`、2D 半径掩码 |
| `features.intensity_fraction` | `0.1` | – | 当日平均亮度相对于当日最大值的最低比例阈值（$\in [0, 1)$） | `limits.py` |
| `features.polarization_fraction` | `0.2` | – | 当日平均偏振亮度相对于当日最大值的最低比例阈值（$\in [0, 1)$；取 `0` 时仅按亮度筛选） | `limits.py` |
| `features.min_radius_rg` | `2.0` | $r_g$ | 2D 像素特征距图像中心的内半径截断（$7.64\,\mu\mathrm{as}$） | `limits.py` |
| `features.pixel_stride` | `5` | 像素 | 2D 像素特征在两个空间轴上的降采样步长 | `limits.py` |
| `features.sector_spacing_degrees` | `10` | $\mathrm{deg}$ | 候选 1D 扇区之间的角度间隔（$\in [1, 360]$） | `limits.py`、`paper.py` |
| `features.sector_offset_degrees` | `4` | $\mathrm{deg}$ | 起始扇区的北偏东角度下标（$4^\circ, 14^\circ, \dots, 354^\circ$；$\in [0, \text{spacing}-1]$） | `limits.py` |
| `features.resultant_min` | `0.7` | – | 全部四天要求的圆合成长度（circular resultant length）$R = \lvert\langle e^{2i\chi}\rangle\rvert \in (0, 1]$ 下限 | `build_observations` |
| `cloud.coefficient` | `0.03` | – | $b(\alpha) = \text{coefficient} \times (\alpha / \text{pivot\_alpha})^{\text{exponent}}$ 的系数（$> 0$） | `cloud_amplitude_scaling` |
| `cloud.pivot_alpha` | `0.05` | – | $b(\alpha)$ 的参考点 $\alpha$（$> 0$） | `cloud_amplitude_scaling` |
| `cloud.exponent` | `1.75` | – | $b(\alpha)$ 的幂律指数 | `cloud_amplitude_scaling` |
| `background.lag_frames` | `3` | 帧 | 无轴子 EVPA 差分的滞后帧数 $L \ge 1$（$3\,t_g = 26.68\,\mathrm{h}$） | `selected_background`、`paper.py` |
| `inference.phase_count` | `64` | 格点 | $[0, 2\pi)$ 上用于对 $\phi_0$ 相位边缘化（phase marginalization）的均匀离散网格数 $N_\phi \ge 2$ | `evpa_templates`（可由 `check_grid.py` `--phase-count` 覆盖） |
| `inference.coupling_points` | `120` | 格点 | `[coupling_min, coupling_max]` 上 $c$ 的对数间隔网格点数（$\ge 2$） | `calculate_observed`（可由 `check_grid.py` `--coupling-points` 覆盖） |
| `inference.coupling_min` | `0.0072992700729927005` | – | $c$ 的先验下界 $c_{\min} = 1/137 > 0$ | `calculate_observed` |
| `inference.coupling_max` | `1000.0` | – | $c$ 的先验上界 $c_{\max} > c_{\min}$ | `calculate_observed`（可由 `check_grid.py` `--coupling-max` 覆盖） |
| `inference.prior` | `log-uniform` | – | $c$ 的先验密度：`log-uniform`（对数均匀先验，$\pi(c) \propto 1/c$）或 `uniform`（均匀先验） | `normalized_posterior_cdf` |
| `inference.credibility` | `0.95` | – | 定义可信上限（credible upper limit）$c_{95}$ 的后验分位数（$\in (0, 1)$） | `posterior_upper_limit`、`RESULTS.md` |
| `inference.exposure_hours` | `6.0` | $\mathrm{h}$ | $\operatorname{sinc}(\omega T / 2)$ 中的矩形曝光时长 $T \ge 0$ | `evpa_templates` |
| `inference.day_gap_hours` | `24.0` | $\mathrm{h}$ | 日对内部两天的间隔 $\tau > 0$（$96 - 95$ 与 $101 - 100$） | `evpa_templates`、`paper.py` |
| `inference.pair_separation_hours` | `120.0` | $\mathrm{h}$ | 两个日对中点之间的时间间隔 $t_2 - t_1 \ge 0$ | `evpa_templates` |

## 本机环境变量（`.env`）

所有本机路径都只由三个路径根目录派生，可写在仓库根目录下一个被 git 忽略的 `.env` 中（由 [`.env.example`](../../.env.example) 复制修改）、在 shell 中导出，或通过命令行传入：

1. `EHT_RAPTOR_SOURCE_ROOT`（`--source-root` / `--raptor-source-root`）：由 `raptor.source` 固定版本的本地 `RAPTOR` Git 工作副本。
2. `EHT_INPUT_ROOT`（`--input-root`）：按固定[输入目录布局](data-layout.md#输入目录布局)（`INPUT_LAYOUT`）存放的只读外部输入：GRMHD 快照、径向表、官方 EHT `uvfits` 文件及成像脚本。
3. `EHT_DATA_ROOT`（`--data-root`）：按[默认目录布局](data-layout.md#默认目录布局)（`LAYOUT`）存放的全部生成产物。

`eht_axion/environment.py` 将 `.env` 严格作为键值数据解析：

- **不执行 Shell 代码、不做变量展开**：不展开 `$VAR` 或 `$(...)`；任何语法错误、未知变量名、缺少 `=value` 或重复键都会抛出 `ValueError`。
- **不在调用者当前工作目录中搜索**：默认 dotenv 文件固定为仓库根目录下的 `ROOT / ".env"`。`--env-file PATH` 用另一个文件替换它；`--no-env-file` 禁用 dotenv 读取。
- **优先级顺序**：`CLI option > exported environment (os.environ) > .env`。显式的单文件选项（`--fluid-root`、`--radial-low`、`--radial-high`、`--uvfits-dir`、`--official-script` 或任一阶段的 `--...-root`）优先于由对应根目录派生的默认路径。
- **相对路径基准**：写在 `.env` 中的相对路径统一以仓库根目录为基准；命令行传入的相对路径（以及 shell 导出的值）以调用者的当前工作目录为基准解析。

| 变量名 | 类型 | 回退值 / 默认值 | 读取脚本 | 说明 |
|---|---|---|---|---|
| `EHT_RAPTOR_SOURCE_ROOT` | 路径 | `build.py`、`generate.py`、`paper.py` 必填 | `scripts/raptor/build.py`、`generate.py`、`scripts/analysis/paper.py` | 处于 `raptor.source.commit` 的 `github.com/projecteht-pku/raptor` 本地工作副本（`--source-root`；`paper.py` 中为 `--raptor-source-root`） |
| `EHT_INPUT_ROOT` | 路径 | 未显式给出全部输入文件时必填 | `scripts/raptor/generate.py`、`scripts/eht/reconstruct.py`、`scripts/analysis/paper.py` | 只读[输入目录布局](data-layout.md#输入目录布局)的根目录：`grmhd/`（`--fluid-root`）、`radial/` 径向表（`--radial-low`、`--radial-high`）与 `eht/`（`--uvfits-dir`、`--official-script`）；可由 `--input-root` 覆盖 |
| `EHT_DATA_ROOT` | 路径 | 未显式传 `--...-root` 时必填 | 所有流水线脚本 | [默认目录布局](data-layout.md#默认目录布局)的数据根目录，包括 `build.py` / `generate.py` 使用的 `raptor/`；可由 `--data-root` 覆盖 |
| `EHT_RAPTOR_THREADS` | 数值（`int`） | `2` | `scripts/raptor/generate.py` | 每个 `RAPTOR` 进程的 OpenMP 线程数（`--threads`） |
| `EHT_RAPTOR_WORKERS` | 数值（`int`） | `1` | `scripts/raptor/generate.py` | 并发运行的 `RAPTOR` 工作进程数（`--workers`） |
| `EHT_RAPTOR_CPU_AFFINITY` | 数值（`str`） | `none`（`.env.example` 示例为 `physical`） | `scripts/raptor/generate.py` | CPU 亲和性（CPU affinity）策略：`none`、`physical` 或 `logical`（`--cpu-affinity`） |
| `EHT_RAPTOR_CPU_LIST` | 数值（`str`） | 未设置（全部允许的 CPU） | `scripts/raptor/generate.py` | 亲和性设为 `physical` 或 `logical` 时允许分配的 CPU 编号/区间，如 `0-47`（`--cpu-list`） |
| `EHT_RAPTOR_TIMEOUT_SECONDS` | 数值（`float`） | 未设置（无单作业超时） | `scripts/raptor/generate.py` | 单个 `RAPTOR` 作业的墙钟超时秒数（`--timeout-seconds`） |
