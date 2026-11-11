# EHT 偏振成像

[← 返回文档索引](index.md)

EHT 成像流水线基于 `ehtim == 1.3.2` 复现 EHT Paper VII 的 2017 年事件视界望远镜（Event Horizon Telescope, EHT）M87\* 偏振图像重建方法（其中总强度成像沿用 Paper IV），随后通过精确凸优化偏振求解器精化每日的仪器偏振泄漏（instrumental polarization leakage，即 D 项，D-term）解，并执行每天 1000 次抽样的 D 项蒙特卡罗（Monte Carlo）模拟，汇总得到的图像系综作为[分析](analysis.md)的观测输入。本部分不读取任何 `RAPTOR` 模拟产品。

## 前置条件与输入数据

`scripts/eht/reconstruct.py`、`refine.py`、`export_bundle.py` 和 `sample.py` 在 `eht-imaging` Pixi 环境（含 `ehtim 1.3.2`、`astropy`）中运行；`scripts/eht/summarize.py` 仅需 `eht-summary`（`astropy`）。`scripts/eht/reconstruct.py` 从只读输入根目录（`EHT_INPUT_ROOT` / `--input-root`，见[输入目录布局](data-layout.md#输入目录布局)）的 `eht/` 目录读取三个本机数据路径；每一项都可用显式选项覆盖，相对路径以当前工作目录为基准解析：

- **HOPS `uvfits` 目录**（`--uvfits-dir`，默认 `$EHT_INPUT_ROOT/eht`）：包含第 `95`、`96`、`100`、`101` 天（2017 年 4 月 5、6、10、11 日）公开发布的四个低频段文件 `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits`，其文件名和 SHA-256 摘要必须与 `configs/experiment.yaml` 中的 `eht_imaging.inputs` 一致。
- **馈源旋转表**（`--feed-rotation-file`，默认 `$EHT_INPUT_ROOT/eht/feed-rotation.json`）：含 `"status": "verified"`、`"units": {"fr_off": "...degrees..."}` 以及各台站 `fr_par`、`fr_elev`、`fr_off` 的 JSON。
- **官方 EHT 成像脚本**（`--official-script`，默认 `$EHT_INPUT_ROOT/eht/eht-imaging_pipeline.py`）：官方总强度成像流水线脚本。

在任何成像计算之前，所选观测日 `uvfits` 文件的 SHA-256 会与 `eht_imaging.inputs` 核对；`uvfits` 文件、馈源旋转表和官方脚本均记入阶段输入，其路径作为成像参数 `input_root`、`feed_rotation_file` 与 `official_script`（`EHT_LOCATIONS`）记录下来。

## 运行 EHT 成像流水线

在 `.env` 中设置 `EHT_INPUT_ROOT` 与 `EHT_DATA_ROOT` 后，对全部四个观测日（`95`、`96`、`100`、`101`）依次运行逐日重建与精化，然后打包导出、执行蒙特卡罗抽样并汇总系综：

```bash
# Optional dry-run plan verification (no ehtim import or disk writes)
pixi run -e eht-imaging eht-reconstruct --stage I --day 95 --plan

# 1–3. Per-day fixed-I, fiducial P/D, and convex D refinement
for day in 95 96 100 101; do
  pixi run -e eht-imaging eht-reconstruct --stage I --day $day
  pixi run -e eht-imaging eht-reconstruct --stage fiducial --day $day
  pixi run -e eht-imaging eht-refine --day $day
done

# 4. Bundle the four fixed-I and refined-D stages for Monte Carlo sampling
pixi run -e eht-imaging eht-export-bundle

# 5. Run the 4 x 1000 D-term Monte Carlo (resumable in $EHT_DATA_ROOT/eht/sample)
pixi run -e eht-imaging eht-sample --workers 64 --threads 2

# 6. Summarize certified blurred FITS draws into $EHT_DATA_ROOT/eht/ensemble
pixi run -e eht-summary eht-summarize
```

除 `eht/sample` 外，所有输出目录都必须是全新目录。`scripts/eht/sample.py` 支持 `--workers`（默认 `64`）、`--threads`（默认 `2`）、`--max-active-memory-gib`（默认 `100`）、`--resource-wait`（默认 `300` s；若内存准入持续受阻则以状态码 `3` 退出，可安全续跑）、`--task-limit` 和 `--resume`。只要批次计划（`batch.json`：打包文件、科学参数、代码提交号、依赖包版本、线程数）不变，对同一输出目录重新运行 `eht-sample` 会自动断点续跑。

## 图像重建方法

### 原生圆偏振数据载入（`eht_axion/eht/calibration.py`）

`load_native` 以原生圆偏振表示（`polrep="circ"`，相关乘积 `RR`、`LL`、`RL`、`LR`）直接载入每个 HOPS `uvfits` 文件，不经由斯托克斯（Stokes）参量中转：
- 在扫描平均前，只要四个相关乘积或误差中任一项非有限或满足 $\sigma \le 0$，整条可见度（visibility）记录即被丢弃。
- 从 `feed-rotation.json` 写入台站馈源旋转系数（`fr_par`、`fr_elev`、`fr_off`），将台站 D 项置为 `0j`，并设置校准标志 `ampcal = phasecal = frcal = True`、`dcal = False`（`fixed_sites` 中的 ALMA/APEX/SMA/JCMT 已在 HOPS 产品中预改正，`solve_sites` 中的 `AZ`、`LM`、`PV` 尚未改正）。
- 对可见度做一次相干扫描平均（`avg_coherent(0, scan_avg=True)`）。
- `make_imaging_obs` 封装 `Obsdata`，既修复了 NumPy 2 下等长 `tlist` 分组丢失结构化数据类型的问题，也修正了 `ehtim` 1.3.2 在计算圆偏振比值误差 `msigma` 时 $\operatorname{Var}(I) = (\sigma_{RR}^2 + \sigma_{LL}^2)/4$ 漏掉 $1/4$ 因子的问题。

### 阶段 `I`：固定总强度与共享增益

`scripts/eht/reconstruct.py --stage I` 发布阶段 `eht:i`：
1. 复制官方 EHT 成像脚本，并在源码中**仅做两处精确字符串替换**后运行：`"ttype     = 'nfft'"` $\to$ `"ttype     = 'direct'"`（直接傅里叶变换）以及 `'fov       = 128*eh.RADPERUAS'` $\to$ `'fov       = 120*eh.RADPERUAS'`（视场，field of view, FOV）。其余所有设置——包括致密源零基线流量 `zbl = 0.60` Jy 和生成最终图像 `im_out` 的 $5\,\mu\mathrm{as}$ 恢复模糊——全部保持原样。代码校验 `im_out` 严格为正、有限且位于 $64 \times 64$、$120\,\mu\mathrm{as}$ 网格上，保存为 `final-fixedI.fits`。
2. 在扫描平均后的原生圆偏振数据上剔除短于 `uv_min`（$10^8\,\lambda = 100\,\mathrm{M}\lambda$）的同址基线，以 `final-fixedI.fits` 为模型，通过 `eh.selfcal(..., pol='I', method='both', ttype='direct', gain_tol=[0.02, 0.2])` 对全部在测台站（含 `PV`）求解**左右旋相等（$R = L$）的共享台站增益**。
3. 将共享增益施加（`applycal`）到原生圆偏振可见度上，写出 `clean-gaincal-pol.uvfits`、`shared-gains.npz` 和 `polarization-input.json`。

后续全部阶段将总强度 $I$ 固定为 `final-fixedI.fits`，仅求解线偏振 $P = Q + iU$（固定圆偏振 $V = 0$）。

### 阶段 `fiducial`：基准 $P$ 成像与前向 $D$ 拟合交替

`scripts/eht/reconstruct.py --stage fiducial` 发布阶段 `eht:fiducial`。在施加了 `uv_min` 截断的增益校准数据上，从均匀线偏振 $P = m_0 I e^{2i\chi_0}$（`initial_m = 0.1`、`initial_evpa_deg = 0.0`）和 `dterms = {}` 出发，执行 `rounds`（50）轮交替，$r = 1, \dots, 50$：
1. 当 $r \in \{2, 12, 22, 32, 42\}$（`restart_rounds`）时，由强度加权模糊后的复线偏振度重新初始化 $P$：
   $$P_{\mathrm{init}} = (1 - \epsilon)\,I\,\frac{K_{20} * P}{K_{20} * I},$$
   其中圆高斯核半高全宽（full width at half maximum, FWHM）为 `restart_blur_uas`（$20\,\mu\mathrm{as}$），收缩量 $\epsilon = 10^{-8}$（`restart_domain_epsilon`）保证初值严格落在物理域 $\lvert P/I\rvert < 1$ 内。
2. 将相对系统噪声 `fsys = 0.005` 按平方和加到 `RL` 和 `LR` 方差上（`add_fsys`：$\sigma_{RL,LR}^2 \to \sigma_{RL,LR}^2 + \bigl(f_{\mathrm{sys}}\lvert(RR + LL)/2\rvert\bigr)^2$），并通过逆琼斯（Jones）矩阵一次性施加当前 D 项（`apply_fixed_d`，期间临时将台站表 SEFD 置零以保留输入的对角边缘误差）。
3. 使用 `ehtim` `Imager` 求解 $P$（`pol='P'`、`transform=['mcv']`、`data_term={'pvis': 1.0, 'm': 1.0}`、`reg_term={'hw': 100.0, 'ptv': 1.0}`、`norm_reg=True`、`ttype='direct'`、`maxit=500`、`stop=1e-8`；`shift_objective_constant=True` 为目标函数加上常数 $w_{\mathrm{pvis}} + w_m + h_w\ln 2$，不改变极小点）。
4. 针对**未做 D 改正**但已加 `fsys` 的基础可见度数据和新得到的 $P$ 图像，调用 `fit_forward_d` 重新拟合 `solve_sites`（`AZ`、`LM`、`PV`）的 D 项：该函数是对 `ehtim.calibrating.pol_cal_new.leakage_cal_new`（`pol_fit=['RL', 'LR']`、`ttype='direct'`）的私有两行适配，将 D 项高斯先验及其梯度显式置零。

第 50 轮结束后，执行最后一次 `final-P-only` 求解，使交付的 `final-raw-IQU.fits` 和 `final-blurred-IQU.fits`（$19\,\mu\mathrm{as}$ 模糊）与 `final-dterms.json` 严格对应。

### 凸优化精化（`scripts/eht/refine.py` 与 `eht_axion/eht/convex.py`）

在固定 $I > 0$ 和固定 D 项时，启用正则项归一化的 `ehtim` 偏振目标函数在物理圆盘域 $\lvert P_x\rvert \le I_x$ 上关于 $P \in \mathbb{C}^{64 \times 64}$ 是严格凸函数：

$$
J(P) = \sum_{t \in \{\mathrm{pvis}, m\}} \frac{w_t}{2 N_t} \sum_{j=1}^{N_t} \frac{\lvert (A_t P)_j - d_{t,j}\rvert^2}{\sigma_{t,j}^2} + \frac{h_w}{F} \sum_x I_x\,h\!\left(\frac{\lvert P_x\rvert}{I_x}\right) + p_{\mathrm{tv}}\,\frac{\operatorname{TV}(P)}{\tau_{\mathrm{TV}}},
$$

其中偏振熵核 $h(m) = \frac{1}{2}\bigl[(1+m)\ln(1+m) + (1-m)\ln(1-m)\bigr]$，$F = \sum_x I_x$，$\tau_{\mathrm{TV}}$ 为 `ehtim` 按波束标度的全变差（total variation, TV）归一化常数，数据项三元组 $(d_t, \sigma_t, A_t)$ 直接取自 `ehtim` 的 `Imager._data_tuples`。$J(P)$ 等于 `ehtim` 目标函数加上常数 $w_{\mathrm{pvis}} + w_m + h_w\ln 2$。

`solve_problem` 在约化变量 $x = P / I$ 上采用对角预条件的原始–对偶混合梯度法（primal–dual hybrid gradient, PDHG，步长安全系数 $\theta = 0.99 < 1$，径向熵邻近算子采用带保护区间的牛顿–二分法精确求解），并每隔 `check_every`（`500`）步计算一次可行的 Fenchel–Rockafellar 对偶下界。当最佳原始目标值与最佳对偶下界之差满足 $\text{best\_primal} - \text{best\_lower} \le \text{gap\_atol}$（`1e-6`）时，返回的解即获得**数值对偶间隙（duality gap）证书**（浮点证书，非区间算术证明）。

`scripts/eht/refine.py`（阶段 `eht:refine`）以基准阶段的 `final-raw-IQU.fits` 和 `final-dterms.json` 为起点，执行 `outer_rounds = 100` 轮带证书的 PDHG $P$ 求解（`maxiter = 20000`、`time_limit = 300` s）与前向 `L-BFGS-B` D 项重拟合（参数 `eht_imaging.refine.d_optimizer`，逐次记录于 `D-optimizer-records.json`），最后执行第 101 次 `final-P` 求解以对齐最终 D 项。

### 打包导出、D 项蒙特卡罗与系综汇总

1. **打包导出（`scripts/eht/export_bundle.py`）**：检查四个 `eht:i` 阶段和四个 `eht:refine` 阶段，确认每天均满足 `all_P_numerically_certified == True`，将以只读权限（`0444`）复制的 `<day>/fixedI.fits`、`<day>/gaincal-pol.uvfits`、`<day>/fiducial-dterms.json`、`lineage/<day>/*` 和 `feed-rotation.json` 连同 `bundle.json`（模式 `eht-bundle`）导出到打包目录。
2. **D 项蒙特卡罗（`scripts/eht/sample.py`）**：
   - 对 `solve_sites` 中的每个台站和每种圆偏振旋向（`R`、`L`），由四天精化解计算 $(\operatorname{Re} D, \operatorname{Im} D)$ 的均值与 $2 \times 2$ 样本协方差（`covariance_ddof = 1`，见 `d_distribution`）。
   - 每天以随机数种子 `SeedSequence([20210405, day])` 抽取 `uncertainty_draws = 1000` 组 D 项（`sample_dterms`）：`AZ`、`LM`、`PV` 服从二元正态分布；预改正台站 `AA`、`AP`、`SM`、`JC` 在每个实部和虚部分量上独立抽取标准差为 `residual_d_sigma = 0.01` 的零均值正态残差。不重抽样热噪声，也不对每次抽样重拟合 D 项。
   - 每次抽样从冷启动初值 $P = 0.05\,I$（`initial_p_fraction`）出发求解精确凸优化 $P$ 问题（`gap_atol = 1e-6`、`maxiter = 50000`、`time_limit = 600` s）；若第一次尝试未获证书，则自动以 `retry_maxiter = 200000`、`retry_time_limit = 1800` s 重试一次。每次抽样写出 `best-raw-IQU.fits`、`best-blurred-IQU.fits`（$19\,\mu\mathrm{as}$ 模糊）、`certificate.json` 和 `taskreceipt.json`；仅当全部 $4000$ 次抽样均获得证书时才写出 `finish.json`。
3. **系综汇总（`scripts/eht/summarize.py`）**：逐一核对 `batch.json`、`finish.json`、以及全部 $4000$ 份 `taskreceipt.json` 与 `certificate.json`，发布阶段 `eht:ensemble`，生成 `axes.npz`（`east_uas`、`north_uas`）、`report.json` 及逐日目录 `day-<day>/{stokes.npy, evpa.npy, global.npz}`：
   - `stokes.npy`：float64 `(1000, 64, 64, 4)`（原生 FITS 方向下的 `[draw, y, x, IQUV]`）；
   - `evpa.npy`：float64 `(1000, 64, 64)`，取 $[0, \pi)\,\mathrm{rad}$ 内的电矢量位置角（EVPA；$I \le 0$ 或 $Q = U = 0$ 处为 `NaN`）；
   - `global.npz`：在 $19\,\mu\mathrm{as}$ 模糊图像上按全视场像素求和计算的 Paper VII 第 5.3 节全局诊断量（`total_I`、`m_net`、`m_resolved`、`evpa_rad`、`beta2_real`、`beta2_imag`），用于与 Paper VII 对比（由 `paper.py` 读取，`limits.py` 不使用）。

## 与 EHT Paper VII 的已知差异

- **方法复现而非逐比特图像复现**：官方总强度脚本运行前做了两处适配替换（`ttype = 'direct'` 替换 `'nfft'`，以及 `fov = 120` $\mu\mathrm{as}$ 替换 `128` $\mu\mathrm{as}$）。
- **本地求解的 D 项与 Paper VII 表 5 存在差异**：`reconstruct.py` 与 `refine.py` 在本地求解出的台站 D 项不同于 EHT Paper VII 表 5 所列数值；最大差异出现在 2017 年 4 月 10 日（第 `100` 天）的 `PV` 台站，其 $D_R$ 相差约 19 个百分点。代码不与表 5 比对，也不代入表 5 的数值。
- **图像定中心采用零平移**：EHT Paper VII 第 5.2/5.3 节通过最大化模糊后 Stokes $I$ 与 Paper IV 共识图像（consensus image）的互相关来平移对齐每幅图像；在本复现中这一对齐为零平移：所有 $\beta_2$ 诊断量和全部 1D/2D 分析特征均直接围绕原生网格中心计算，不提供重定中心选项。
- **数值重启计划与求解器预算**：`restart_rounds`、`restart_blur_uas`、`maxit` 以及 `eht_imaging.refine` 和 `eht_imaging.monte_carlo` 中的 PDHG / `L-BFGS-B` 精化与蒙特卡罗预算均为数值实现层面的选择，而非论文规定的物理参数。
