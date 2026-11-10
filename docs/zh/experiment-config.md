# 实验配置

## 用途

`configs/experiment.yaml` 是本实验唯一的参数文件。每个科学取值在其中只写一次，每个脚本都通过
`eht_axion/experiment.py` 读取它。本页列出每个键的含义、单位、当前取值和读取它的代码，并说明
校验规则和派生设置。

## 输入

- 仓库中的 `configs/experiment.yaml`（所有 `--config` 选项的默认值）。
- 可以用 `--config PATH` 指定另一个文件；它必须包含完全相同的键。

## 输出

加载器不写任何文件。校验后的字典，或某个阶段用到的那部分，记录在阶段清单和计划中
（[来源记录](provenance.md)）：

- `raptor_generate.py` 和交付步骤记录光线追踪子集（`ray_tracing_config`）。
- EHT 成像各阶段记录扁平参数集（`eht_imaging_parameters`）。
- `analysis_limits.py` 把整个文件记录为 `scientific_profile`。

## 运行方法

该文件由各脚本读取，没有单独的命令。快速校验：

```bash
pixi run python -c "from eht_axion.experiment import read_experiment; read_experiment()"
```

## 细节

### 结构与读取规则

文件恰好有四个节：

- `shared`：多个部分共用的取值。
- `raptor`：决定光线追踪 Stokes 图像的取值。
- `eht_imaging`：EHT 成像方法、其数据身份和求解器设置。
- `analysis`：从模板和 EHT 系综到 $c$ 的上限。

每个脚本读取 `shared` 和它自己的节；分析部分还读取它需要的 `raptor` 和 `eht_imaging` 取值。
机器路径和执行开关不在该文件中（`.env` 和 `.env.raptor` 见[开发](development.md)）。

### 校验

`read_experiment` 用一个拒绝重复键的 YAML 加载器解析文件，然后由 `validate` 检查。校验内容包括
结构（每个节的键集合必须完全一致）、类型（布尔值不算数字，整数必须是整数）和取值范围。它不填充
任何默认值。有两项交叉检查依赖派生值：分析的天空视场（`shared.image.fov_uas` 换算成 $r_g$）必须
落在 `raptor.camera_fov_rg` 之内，且 `analysis.features.min_radius_rg` 必须小于孔径半径。

### 派生设置

| `eht_axion/experiment.py` 中的函数 | 结果 | 使用者 |
|---|---|---|
| `physical_scale` | `PhysicalScale`：$t_g = GM/c^3$（按下列取值为 8.893 h）、$r_g$ 的角大小（3.819 µas）及换算 | 分析脚本、交付的帧间隔 |
| `ray_tracing_config` | `shared` 中的物理键（`mass_solar`、`distance_mpc`、`label_step_tg`、`alphas`、`principal_n`）加上不含 `raptor.source` 的 `raptor` 节 | RAPTOR 计划、交付身份（`generated.json`） |
| `analysis_geometry` | 预处理参数（`sigma`、`init_fov`、`final_fov`、`final_shape`、`rotation_angle`，单位为 $r_g$ 和弧度）和孔径参数（`fov`、`radial_max`、`window_degrees`、`angular_bins`） | `analysis_templates.py`、`analysis_limits.py` |
| `eht_imaging_parameters` | `eht_imaging` 取值（不含 `refine` 和 `monte_carlo`）、天空网格（`npix`、`fov_uas`、`blur_fwhm_uas`）和四个机器位置组成的扁平字典 | EHT 成像各阶段、`eht_refine.py`、`eht_export_bundle.py` |

`ANGULAR_BINS = 360`（一度宽的孔径分格）和 `EHT_DAYS = (95, 96, 100, 101)` 是同一模块中的代码
常量，不是 YAML 取值。`EHT_DAYS` 是代码中唯一的观测日列表（EHT 成像脚本、系综汇总和分析都导入
它）；校验要求 `eht_imaging.inputs` 恰好包含这些日期。

### shared 节

| 键 | 取值 | 单位 | 含义 | 读取者 |
|---|---|---|---|---|
| `mass_solar` | 6.5e9 | $M_\odot$ | M87\* 质量 | `physical_scale`；RAPTOR `MBH` |
| `distance_mpc` | 16.8 | Mpc | 到 M87\* 的距离 | `physical_scale`；RAPTOR `DISTANCE`（以 kpc 写入） |
| `label_step_tg` | 0.1 | 每个标签的 $t_g$ | GRMHD 快照标签每个单位对应的时间 | 计划中对快照尾部时间间隔的检查 |
| `alphas` | 0.10 … 0.44（8 个取值） | – | 引力精细结构常数 $\alpha$ | RAPTOR 任务和交付、模板、似然网格 |
| `principal_n` | 2 | – | $\Omega(\alpha)=\alpha\,(1-\alpha^2/(2n^2))$ 中 211 态的主量子数 $n$ | RAPTOR `AXION_OMEGA`；信号模型 |
| `image.shape` | 64 | 像素 | 公共天空网格的边长 | EHT 成像 `npix`；分析重采样；系综汇总 |
| `image.fov_uas` | 120.0 | µas | 公共天空网格的视场（像素 1.875 µas） | EHT 成像（替换官方脚本中的 128 µas）；分析重采样和孔径 |
| `image.blur_fwhm_uas` | 19.0 | µas | 施加于 EHT 图像（即“模糊”图像）和模拟图像的圆高斯核的半高全宽（FWHM）；这是一个额外的匹配核，不是等效点扩散函数（PSF） | EHT `save_image`；分析预处理的 $\sigma$ |

### raptor 节

| 键 | 取值 | 单位 | 含义 | 读取者 |
|---|---|---|---|---|
| `source.repository` | `github.com/projecteht-pku/raptor` | – | RAPTOR 仓库，写成 host/owner/name（小写）；本地工作副本必须有一个 git 远程仓库指向这里 | `raptor_build.py`（远程仓库检查）；编译记录；计划中的 `raptor_source` |
| `source.commit` | `2005600…`（40 位十六进制） | – | RAPTOR 提交；本地工作副本必须没有修改且恰好处于该提交 | `raptor_build.py`；编译记录；计划中的 `raptor_source`，由 `run_jobs` 和交付步骤与每个编译比对 |
| `camera_fov_rg` | 46.2261 | $r_g$ | 相机宽度（176.5 µas） | RAPTOR `CAM_SIZE_X/Y`；预处理 `init_fov` |
| `fluid_spin` | 0.9375 | – | BHAC 吸积流的黑洞自旋；每个快照尾部都必须给出该值 | `read_snapshot_footer`（计划、作业运行器） |
| `build.model` | `bhac` | – | `setup.sh --code` | `raptor/build.py` |
| `build.metric` | `mks` | – | `--metric` | `raptor/build.py` |
| `build.integrator` | `rk4` | – | `--int` | `raptor/build.py` |
| `build.grid` | `SMR` | – | `--grid` | `raptor/build.py` |
| `build.polarized` | true | – | `--rad=pol` | `raptor/build.py` |
| `build.sfc` | true | – | `--sf=sfc` | `raptor/build.py` |
| `build.magnetic_polarity` | `plus` | – | `--bflip` | `raptor/build.py` |
| `build.sigma_cut` | 1.0 | – | 磁化强度截断，`--sigma-cut` | `raptor/build.py` |
| `build.thetae_min` | 0.001 | – | `--thetae-min` | `raptor/build.py` |
| `build.thetae_max` | 100.0 | – | `--thetae-max` | `raptor/build.py` |
| `build.uniform_output` | false | – | `--uniform` | `raptor/build.py` |
| `plasma.mass_unit_g` | 3.0e25 | g | GRMHD 吸积流的质量单位，`M_UNIT` | RAPTOR `model.in` |
| `plasma.r_low` | 1.0 | – | 电子温度参数 `R_LOW` | RAPTOR `model.in` |
| `plasma.r_high` | 1.0 | – | 电子温度参数 `R_HIGH` | RAPTOR `model.in` |
| `plasma.inclination_degrees` | 163.0 | 度 | 观测者倾角 | RAPTOR `model.in` |
| `plasma.frequency_hz` | 2.3e11 | Hz | 观测频率（单一频率） | RAPTOR `model.in`；交付数据集名称 |
| `plasma.step_size` | 0.005 | – | 测地线积分步长，`STEPSIZE` | RAPTOR `model.in` |
| `snapshots.start`、`stop`、`step` | 1600、3000、10 | 标签 | 含端点的快照标签，共 141 个快照 | 计划、交付、分析中的标签检查 |
| `response.resolution` | 100 | 像素 | 有轴子图像及其无轴子基线的边长 | 计划、交付、`analysis_templates.py` |
| `response.reference_norm` | 0.1 | – | 轴子场归一化 `AXION_NORM`；模板以该归一化为单位 | RAPTOR `model.in`；信号模型 |
| `background.resolution` | 400 | 像素 | 用于 GRMHD 背景的无轴子图像边长 | 计划、交付、`analysis_templates.py` |
| `radial_tables.low.spin` | 0.9375 | – | 低自旋径向表声明的自旋 | 计划、任务身份 |
| `radial_tables.high.spin` | 0.99 | – | 高自旋径向表声明的自旋（有意与 `fluid_spin` 不同） | 计划、任务身份 |
| `alpha_table` | 0.10–0.30 用 `low`，0.35–0.44 用 `high` | – | 每个 $\alpha$ 使用的径向表 | 计划 |
| `grid.nxlone`、`xprobmin`、`xprobmax`、`hslope` | 见文件 | – | BHAC `grid_mks.in` 坐标；与编译目录中的 `grid_mks.in` 和快照尾部比对 | 计划、作业运行器 |

### eht_imaging 节

| 键 | 取值 | 单位 | 含义 | 读取者 |
|---|---|---|---|---|
| `inputs` | 四个 `[file name, sha256]` | – | 每天 HOPS uvfits 的身份 | `eht_reconstruct.py`（与输入回执和文件字节比对） |
| `gain_tol` | [0.02, 0.2] | – | 共享增益自校准的增益相对容差下限和上限（官方脚本的取值） | `eht_reconstruct.py` |
| `fsys` | 0.005 | – | 按平方和加到 RL/LR 误差上的相对系统噪声：$\sigma^2 \to \sigma^2 + (f_\mathrm{sys}\,\lvert(RR+LL)/2\rvert)^2$ | 基准、精化和蒙特卡罗中的 `add_fsys`；打包约定 |
| `solve_sites` | AZ、LM、PV | – | 求解 D 项的台站 | 基准、精化、蒙特卡罗 D 分布 |
| `fixed_sites` | AA、AP、SM、JC | – | HOPS 产品中已做过泄漏改正的台站（代码中唯一的此类列表） | `load_native` 元数据；`sample_dterms` 的零均值残差抽样 |
| `p_data` | pvis 1.0，m 1.0 | – | P 成像中 ehtim 数据项的权重 | 基准 `Imager`；`eht_refine.py` 和 `eht_sample.py` 的精确 P 求解器；打包约定 |
| `p_reg` | hw 100.0，ptv 1.0 | – | ehtim 正则项权重（熵、全变差） | 同 `p_data` |
| `norm_reg` | true | – | ehtim 正则项归一化 | 同 `p_data` |
| `initial_m` | 0.1 | – | 基准 P 图像的初始偏振度 | `eht_reconstruct.py` |
| `initial_evpa_deg` | 0.0 | 度 | 基准 P 图像的初始 EVPA | `eht_reconstruct.py` |
| `rounds` | 50 | – | 基准 P/D 交替的轮数 | `eht_reconstruct.py` |
| `restart_rounds` | 2、12、22、32、42 | – | 从模糊后的偏振度重新启动 P 的轮次（数值手段，不出自论文） | `eht_reconstruct.py` |
| `restart_blur_uas` | 20.0 | µas | 该重启模糊的 FWHM | `eht_reconstruct.py` |
| `restart_domain_epsilon` | 1e-8 | – | 重启的 P 乘以 $1-\epsilon$，使 $\lvert P/I\rvert<1$ | `eht_reconstruct.py` |
| `maxit` | 500 | – | 每次调用 ehtim 优化器的迭代上限（数值设置） | `eht_reconstruct.py` |
| `stop` | 1e-8 | – | ehtim 收敛容差 | `eht_reconstruct.py` |
| `shift_objective_constant` | true | – | 给 ehtim 目标函数加上 $\sum$`p_data` $+\,h_w\ln 2$（常数；极小点不变） | `eht_reconstruct.py` |
| `uv_min` | 1.0e8 | 波长 | 短于 100 Mλ 的基线（同址“零”基线）在共享增益求解之前以及所有 P 成像中被剔除 | I 阶段、基准、精化、蒙特卡罗 |
| `uncertainty_draws` | 1000 | – | 每天的 D 项抽样数 | `eht_sample.py`、`eht_summarize.py`、`analysis_limits.py`、`analysis_check_covariance.py` |
| `seed` | 20210405 | – | 蒙特卡罗随机种子；每天用 `SeedSequence([seed, day])` | `eht_sample.py` |
| `covariance_ddof` | 1 | – | 四天 D 协方差的自由度修正（delta degrees of freedom） | `eht_sample.py`；打包约定 |
| `residual_d_sigma` | 0.01 | – | `fixed_sites` 的零均值 D 抽样在每个实部/虚部分量上的标准差 | `eht_sample.py`；打包约定 |
| `refine.outer_rounds` | 100 | – | 最终 P 求解之前 P 求解 / D 重拟合的交替次数 | `eht_refine.py` |
| `refine.gap_atol` | 1e-6 | 目标函数单位 | 每次精确 P 求解的对偶间隙绝对容差 | `eht_refine.py` |
| `refine.maxiter` | 20000 | – | 每次调用 P 求解器的迭代次数 | `eht_refine.py` |
| `refine.check_every` | 500 | – | 两次间隙证书之间的迭代数 | `eht_refine.py` |
| `refine.time_limit` | 300.0 | s | 每次调用 P 求解器的时间 | `eht_refine.py` |
| `refine.d_optimizer` | ftol 1e-13，gtol 1e-9，maxiter 10000，maxcor 50，maxls 50，maxfun 100000 | – | 每次 D 重拟合所用的 scipy L-BFGS-B 选项 | `eht_refine.py` → `fit_forward_d` |
| `monte_carlo.initial_p_fraction` | 0.05 | – | 初始 P $=0.05\,I$（实数，EVPA 为 0），不做热启动 | `eht_sample.py` |
| `monte_carlo.gap_atol` | 1e-6 | 目标函数单位 | 每次抽样的间隙容差 | `eht_sample.py` |
| `monte_carlo.maxiter` | 50000 | – | 第一次尝试的 P 求解器迭代数 | `eht_sample.py` |
| `monte_carlo.check_every` | 500 | – | 两次证书之间的迭代数 | `eht_sample.py` |
| `monte_carlo.time_limit` | 600.0 | s | 第一次尝试的时间 | `eht_sample.py` |
| `monte_carlo.retry_maxiter` | 200000 | – | 未获证书的抽样第二次尝试的迭代数 | `eht_sample.py` |
| `monte_carlo.retry_time_limit` | 1800.0 | s | 第二次尝试的时间 | `eht_sample.py` |

官方 I 流程中的致密源流量密度（0.6 Jy，`zbl = 0.60`）固定在官方 EHT 脚本内部；该脚本除两处替换
外原样运行（[EHT 成像](eht-imaging.md#i-阶段)）；P 成像使用固定 I 图像的总流量。

`refine` 和 `monte_carlo` 记录在 `eht_refine.py` 和 `eht_sample.py` 的计划中，而不在成像各阶段的
扁平参数集中。

### analysis 节

| 键 | 取值 | 单位 | 含义 | 读取者 |
|---|---|---|---|---|
| `rotation_degrees` | 108.0 | 度 | 把模拟图像旋转到天空上的角度（图像和 $Q+iU$） | `analysis_geometry` |
| `aperture.window_degrees` | 10 | 度 | 每个孔径扇区的全宽 | `analysis_geometry` → 孔径算子（模板、扇区 EVPA 和扇区选择） |
| `aperture.radial_max_fraction` | 1.0 | – | 孔径半径占半视场（60 µas）的比例 | 同上 |
| `features.intensity_fraction` | 0.1 | – | 保留平均 $I \ge 0.1 \times$ 当日最大值的位置 | `analysis_limits.py` |
| `features.polarization_fraction` | 0.2 | – | 并且平均 $\lvert P\rvert \ge 0.2 \times$ 当日最大值（取 0 时只按亮度选取） | `analysis_limits.py` |
| `features.min_radius_rg` | 2.0 | $r_g$ | 2D 像素到图像中心的最小距离（7.6 µas） | `analysis_limits.py` |
| `features.pixel_stride` | 5 | 像素 | 2D：两个轴向上每隔 5 个像素取一个 | `analysis_limits.py` |
| `features.sector_spacing_degrees` | 10 | 度 | 1D：每 10° 一个扇区 | `analysis_limits.py` |
| `features.sector_offset_degrees` | 4 | 度 | 1D：从北偏东 4° 开始 | `analysis_limits.py` |
| `features.resultant_min` | 0.7 | – | 每一天的最小圆合成长度（circular resultant length） | `build_observations` |
| `cloud.coefficient` | 0.03 | – | $b(\alpha)=0.03\,(\alpha/0.05)^{1.75}$ | `cloud_amplitude_scaling` |
| `cloud.pivot_alpha` | 0.05 | – | $b(\alpha)$ 的参考点 | 同上 |
| `cloud.exponent` | 1.75 | – | $b(\alpha)$ 的指数 | 同上 |
| `background.lag_frames` | 3 | 帧 | 无轴子帧差分的滞后（3 帧 $\approx$ 26.7 h） | `selected_background` |
| `inference.phase_count` | 64 | – | 被边缘化的轴子相位的均匀网格 | `evpa_templates` |
| `inference.coupling_points` | 120 | – | $c$ 的对数间隔网格 | `calculate_observed` |
| `inference.coupling_min` | 1/137 | – | $c$ 的先验下界 | 同上 |
| `inference.coupling_max` | 1000.0 | – | $c$ 的先验上界 | 同上 |
| `inference.prior` | `log-uniform` | – | $c$ 的先验（也允许 `uniform`） | 后验 |
| `inference.credibility` | 0.95 | – | 作为上限报告的后验分位数 | 后验 |
| `inference.exposure_hours` | 6.0 | h | 一次观测的矩形时间平均 | `evpa_templates` |
| `inference.day_gap_hours` | 24.0 | h | 一个日对中两天的间隔 | `evpa_templates` |
| `inference.pair_separation_hours` | 120.0 | h | 两个日对中点之间的间隔 | `evpa_templates` |

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [分析](analysis.md)
- [EHT 成像](eht-imaging.md)
- [RAPTOR 模拟](raptor-simulation.md)
- [开发](development.md)
