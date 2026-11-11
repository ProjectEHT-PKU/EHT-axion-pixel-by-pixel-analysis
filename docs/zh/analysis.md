# 分析

[← 返回文档索引](index.md)

分析流水线将已交付的 `RAPTOR` 斯托克斯（Stokes）数据立方与汇总后的事件视界望远镜（EHT）偏振系综结合起来，计算无量纲耦合常数 $c$ 的相位边缘化（phase-marginalized）贝叶斯可信上限（credible upper limit）$c_{95}(\alpha)$。所有分析脚本均在 `default` Pixi 环境中运行。

## 运行分析流水线

在 `.env` 中设置 `EHT_DATA_ROOT`（或通过 `--data-root` 传入）后，全部输入与输出路径默认采用[标准数据布局](data-layout.md#默认目录布局)；`analysis-paper` 还会从固定版本的 `RAPTOR` 工作副本（`EHT_RAPTOR_SOURCE_ROOT`）读取 `src/axion.c`，并从[输入目录布局](data-layout.md#输入目录布局)（`EHT_INPUT_ROOT`）读取径向表：

```bash
pixi run analysis-templates
pixi run analysis-limits
pixi run analysis-check-grid --check grid --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --check prior --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance
pixi run analysis-summarize
pixi run analysis-paper
```

显式传入目录选项即可覆盖默认路径：

```bash
pixi run analysis-templates --raptor-input-root MERGED/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root GRID \
    --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root PRIOR \
    --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance --source LIMITS --templates TEMPLATES --output-root COV
pixi run analysis-summarize --limits LIMITS --gridcheck GRID --priorcheck PRIOR \
    --covcheck COV --output-root SUMMARY
pixi run analysis-paper --limits LIMITS --gridcheck GRID --priorcheck PRIOR --covcheck COV \
    --templates TEMPLATES --ensemble ENSEMBLE --raptor-input-root MERGED/data \
    --raptor-plan plan.json --raptor-source-root RAPTOR --output-root PAPER
```

当未设置 `--data-root` / `EHT_DATA_ROOT` 时，`scripts/analysis/summarize.py` 的 `--covcheck` 为可选参数；设置了数据根目录时，`summarize.py` 默认将 `--covcheck` 设为 `$EHT_DATA_ROOT/analysis/check-covariance`。

## 模拟响应模板（`scripts/analysis/templates.py`）

### 图像预处理

`scripts/analysis/templates.py` 按当前 `configs/experiment.yaml` 逐一校验交付的数据立方（`generated_inputs`），并通过 `preprocess_cube`（`eht_axion/analysis/templates.py`）流式逐帧处理。对每一帧，`preprocess_frame`：
1. 在原生相机网格上，用标准差 $\sigma = \mathrm{FWHM}/\sqrt{8\ln 2}$（半高全宽 `shared.image.blur_fwhm_uas` $= 19\,\mu\mathrm{as} \to 2.113\,r_g$）的圆高斯核对每个 Stokes 通道做卷积，边界采用最近值延拓（`mode="nearest"`）。
2. 通过双线性插值（`order=1`、`prefilter=False`）将相机视场 `raptor.camera_fov_rg`（$46.2261\,r_g \approx 176.5\,\mu\mathrm{as}$）重采样到 $64 \times 64$、$120\,\mu\mathrm{as}$（$31.42\,r_g$）的天空网格上，同时将图像主动旋转 $\theta =$ `analysis.rotation_degrees`（$108^\circ$）。
3. 按自旋–2（spin-2）变换因子旋转线偏振：$Q + iU \mapsto e^{2i\theta}(Q + iU)$。

不对流量做重新归一化。预处理后的图像具有空间轴 `[x, y]`，第 `0` 轴指向西（`-East`），第 `1` 轴指向北（`+North`）。有轴子响应数据立方（阶段 `analysis:preprocess`）保留全部 17 个状态；400 像素无轴子背景数据立方（阶段 `analysis:background`）仅保留状态 `0`。

### 扇形孔径算子（1D 观测量）

1D 观测量是在 $360$ 个中心间隔 $1^\circ$ 的重叠“披萨片”扇形孔径上积分得到的强度加权相干电矢量位置角（electric-vector position angle, EVPA）：

$$
\chi_{\mathrm{ap}}(\varphi_k) = \frac{1}{2}\arg \int_0^{r_{\max}} \int_{\varphi_k - w/2}^{\varphi_k + w/2} I\,(Q + iU)\,r\,\mathrm{d}\varphi\,\mathrm{d}r,
$$

其中扇区全宽 $w =$ `analysis.aperture.window_degrees`（$10^\circ$），外半径 $r_{\max} =$ `radial_max_fraction` $\times$ 半视场（$60\,\mu\mathrm{as} \approx 15.71\,r_g$）。`aperture_operator`（`eht_axion/analysis/aperture.py`）构造一个形状为 $(360, 64^2)$ 的压缩稀疏行（CSR）双线性求积（bilinear quadrature）权重矩阵 $W_{kj}$（径向步长取 `fov` 单位即 $r_g$ 下的 `0.02`，角向步长 $0.2^\circ$，径向与角向两端点取梯形半权重）。施加算子后将分格循环移位 $-90$，使第 $k \in \{0, \dots, 359\}$ 格的中心位于北偏东 $k^\circ$：

- **扇区 EVPA**（`integrate_stokes`，Paper VII 式 11）：$\chi_k = \frac{1}{2}\arg\sum_j W_{kj}\,I_j(Q_j + iU_j)$。
- **扇区总强度**（`sector_brightness`）：$B_k = \sum_j W_{kj}\,I_j$。
- **强度加权偏振亮度**（`sector_brightness`，Paper VII 式 10）：$\langle\mathcal{P}\rangle_k = \sum_j W_{kj}\,I_j\lvert P_j\rvert / B_k$，其中 $\lvert P_j\rvert = \sqrt{Q_j^2 + U_j^2}$。

完全相同的稀疏算子 $W_{kj}$ 同时作用于模拟帧以及转置为 `[x, y]` 顺序后的每次 EHT 抽样（`[y, x] -> [x, y]`）。

### 谐波响应拟合与快照汇总

对每个快照和每个 $\alpha$，`fit_processed`（阶段 `analysis:fit`）剔除状态 `0`，通过 `fit_phase_dft` 在每个像素（2D）和每个扇区（1D）上对 16 个相位状态（$\psi_k = 2\pi k/16$，$k = 0, \dots, 15$）做一次谐波拟合：

$$
\chi_k = \bar{\chi} + A\cos(\psi_k + \phi).
$$

先把各相位的 EVPA 提升到其圆均值（circular mean）$\chi_0 = \frac{1}{2}\arg\langle e^{2i\chi_k}\rangle$ 所在的 $\pi$ 周期分支上，记基线均值 $\bar{\chi} = \langle\tilde{\chi}_k\rangle$，再由离散傅里叶变换计算 $A e^{i\phi} = \frac{1}{8}\sum_{k=0}^{15} (\tilde{\chi}_k - \bar{\chi})\,e^{-2\pi i k/16}$。缺少任一相位的位置记为无效；当 $A \le 64\varepsilon_{\mathrm{mach}}$ 时置 $A = 0$、$\phi = \mathrm{NaN}$。

`summarize_fits`（阶段 `analysis:templates`）对每个 $\alpha$ 汇总 141 个快照（要求 $\ge 2$ 个有限值快照）：以高斯极大似然估计（Gaussian maximum-likelihood estimation, MLE，`ddof = 0`）计算振幅均值 `amplitude_mean` 和总体标准差 `amplitude_std`，并在圆分支提升后计算相位圆均值 `phase_mean` 和标准差 `phase_std`。`load_templates` 仅在四个统计量对全部 $\alpha$ 均有限且标准差非负时将该空间位置标记为有效；信号模型只读取均值。

## 观测特征与协方差（`scripts/analysis/limits.py`）

### 特征选择

`prepare_limit_inputs` 将每天的 EHT 系综图 `[y, x]` 转置为模板轴序 `[x, y]`（核对网格坐标完全一致），并在不查看日间 EVPA 差值大小的前提下筛选特征：

- **2D 像素特征**必须满足：模板有效；距网格中心的径向距离落在 $[\texttt{min\_radius\_rg}, r_{\max}] = [7.64\,\mu\mathrm{as}, 60\,\mu\mathrm{as}]$ 之内；两个空间下标均为 `pixel_stride`（`5`）的整数倍；在全部四天的全部 1000 次抽样中 EVPA 均有限；且在全部四天中，系综平均 Stokes 图满足 $I \ge 0.1\,I_{\max}$（`intensity_fraction`）和 $\lvert P\rvert \ge 0.2\,\lvert P\rvert_{\max}$（`polarization_fraction`）。
- **1D 扇区特征**为满足 $k \equiv 4 \pmod{10}$（即 $4^\circ, 14^\circ, \dots, 354^\circ$；`sector_offset_degrees`、`sector_spacing_degrees`）的候选扇区，并要求：模板有效；扇区 EVPA $\chi_k$ 在全部四天的每次抽样中均有限；且在全部四天中，抽样均值满足 $\bar{B}_k \ge 0.1\max_j \bar{B}_j$ 和 $\overline{\langle\mathcal{P}\rangle}_k \ge 0.2\max_j \overline{\langle\mathcal{P}\rangle}_j$。

随后 `build_observations`（`eht_axion/analysis/observations.py`）仅保留在全部四天中抽样圆合成长度（circular resultant length）$R = \lvert\langle e^{2i\chi}\rangle\rvert$ 均不低于 `resultant_min`（$0.7$）的特征。

### 逐日统计矩、日对差分与 GRMHD 背景

对每个观测日 $d \in \{95, 96, 100, 101\}$，`moments` 计算 $\pi$ 周期圆均值 $\bar{\chi}_d = \frac{1}{2}\arg\langle e^{2i\chi}\rangle$ 以及去中心化的总体协方差（population covariance，`ddof = 0`）：

$$
C_d = \frac{1}{N}\sum_{k=1}^{N} r_k r_k^{\mathsf{T}}, \qquad r_k = \operatorname{wrap}_\pi(\chi_k - \bar{\chi}_d) - \bigl\langle\operatorname{wrap}_\pi(\chi - \bar{\chi}_d)\bigr\rangle.
$$

$C_d$ 刻画单次 D 项抽样的校准不确定度（而非蒙特卡罗均值的标准误）。将不同观测日视为相互独立，即得两组观测日对差分向量与分块对角的日对协方差矩阵：

$$
d_1 = \operatorname{wrap}_\pi(\bar{\chi}_{96} - \bar{\chi}_{95}), \quad C_1 = C_{95} + C_{96}; \qquad
d_2 = \operatorname{wrap}_\pi(\bar{\chi}_{101} - \bar{\chi}_{100}), \quad C_2 = C_{100} + C_{101}.
$$

$C_1$ 与 $C_2$ 都必须在不加任何正则化项的条件下通过 Cholesky 正定分解。代码同时报告距离 $\pm\pi/2$ 分支割线的距离（以边缘标准差 $\sigma$ 为单位的 $\pi/2 - \lvert d\rvert$）与特征值条件数作为诊断量，但不据此剔除特征。

**GRMHD 背景散布（`selected_background`）：** 严格在最终选定的 1D 或 2D 特征集合上，从 141 帧预处理后的 400 像素无轴子图像中提取 EVPA 时间序列 $\theta_t$，并在滞后帧数 $L =$ `analysis.background.lag_frames`（$3$ 帧 $= 3\,t_g = 26.68\,\mathrm{h}$，帧间隔由 BHAC 尾部时间给出）处计算差分 $\Delta_t = \operatorname{wrap}_\pi(\theta_{t+L} - \theta_t)$。将全部选定特征和帧对合并，求去中心化总体标准差：

$$
\sigma_{\mathrm{bg}} = \sqrt{\bigl\langle(\Delta - \langle\Delta\rangle)^2\bigr\rangle},
$$

由于 $\sigma_{\mathrm{bg}}$ 直接刻画约一天间隔的帧间差值，因此不再额外乘以 $\sqrt{2}$ 因子。

## 统计推断

### 双历元信号模型

在 `eht_axion/analysis/model.py`（`evpa_templates`、`cloud_amplitude_scaling`）中，第 $p \in \{1, 2\}$ 个日对在特征 $x$ 处的模型 EVPA 差值关于耦合常数 $c$ 呈线性关系：

$$
s_p(x; \alpha, \phi_0) = b(\alpha)\,\frac{2\sin(\omega\tau / 2)\,\operatorname{sinc}(\omega T / 2)}{N_{\mathrm{ref}}}\,A(x, \alpha)\,\sin\!\bigl(\phi_0 + \phi(x, \alpha) + \omega t_p\bigr),
$$

其中角频率 $\omega = \Omega(\alpha)/t_g$（$\mathrm{rad\,h^{-1}}$），日对内两天的间隔 $\tau = 24\,\mathrm{h}$（`day_gap_hours`），单次观测矩形曝光时长 $T = 6\,\mathrm{h}$（`exposure_hours`），$\operatorname{sinc}(u) = \sin u / u$，两日对中点时刻 $t_1 = 0$ 与 $t_2 = 120\,\mathrm{h}$（`pair_separation_hours`，即以第一个日对中点作为相位原点），参考归一化 $N_{\mathrm{ref}} = 0.1$（`raptor.response.reference_norm`），轴子云振幅标度 $b(\alpha) = 0.03\,(\alpha / 0.05)^{1.75}$。模板均值 $A$ 与 $\phi$ 在模板节点之间对 $\alpha$ 做线性插值（$\phi$ 在单位圆上插值），在节点范围之外保持边界常数。

### 协方差模式、似然函数与可信上限

`calculate_observed`（`eht_axion/analysis/likelihood.py`）对三种有效日对协方差模型 $\tilde{C}_p$ 进行评估（三种模式均不含正比于 $c^2$ 的响应模型方差）：

| 模式 | 有效日对协方差 $\tilde{C}_p$ |
|---|---|
| `std` | $\operatorname{diag}(C_p)$ |
| `cor` | $C_p$ |
| `cor_fixed_background` | $C_p + \sigma_{\mathrm{bg}}^2\,\mathbb{1}$ |

在 $c$（区间 $[1/137, 1000]$ 上的 `120` 个对数等距格点）与两组日对共享相位 $\phi_{0,j} = 2\pi j / N_\phi$（$N_\phi = 64$，$j = 0, \dots, N_\phi - 1$）组成的联合网格上，构造主值分支高斯二次型（残差 $d_p - c\,s_p$ 不再重新回绕，且略去与 $c$ 和 $\phi_0$ 无关的对数行列式常数），并以均匀离散先验对 $\phi_0$ 边缘化：

$$
\chi^2(c, \phi_0) = \sum_{p=1}^{2} \bigl(d_p - c\,s_p(\phi_0)\bigr)^{\mathsf{T}} \tilde{C}_p^{-1} \bigl(d_p - c\,s_p(\phi_0)\bigr), \qquad
\ln L(c) = \ln \frac{1}{N_\phi} \sum_{j=0}^{N_\phi - 1} e^{-\chi^2(c, \phi_{0,j}) / 2}.
$$

`normalized_posterior_cdf`（`eht_axion/analysis/posterior.py`）将似然 $L(c)$ 乘以先验密度 $\pi(c) \propto 1/c$（`log-uniform` 对数均匀先验；`uniform` 则为常数），对 $c$ 采用梯形法则数值积分归一化后验密度；`posterior_upper_limit` 在 `credibility = 0.95` 处对该分段线性后验密度的累积分布函数（cumulative distribution function, CDF）求解析逆，得到 95 % 可信上限 $c_{95}$。每种维度与协方差组合都同时在真实观测差值（`observed`）和**零数据对照**（zero-data control，固定相同特征与协方差并令 $d_1 = d_2 = 0$，标记为 `zero-data`）上运行。

## 敏感性检验、结果汇总与论文输出

- **网格与先验检验（`scripts/analysis/check_grid.py`）**：直接读取 `analysis:limits` 保存的 `observations-{1d,2d}.npz`，用命令行覆盖的 `--phase-count`、`--coupling-points` 和/或 `--coupling-max` 重新计算全部上限（阶段 `analysis:check-grid`），报告相对于源阶段的 $\Delta\log_{10} c_{95}$，不重新筛选特征，也不重复遍历上游原始数据哈希。
- **半样本协方差检验（`scripts/analysis/check_covariance.py`）**：将每天的 1000 次抽样拆分为前一半 `first500` 与后一半 `last500`，围绕固定的全系综圆均值 $\bar{\chi}_d$ 重新计算 $C_d$，并报告六种观测情形下的 $\Delta\log_{10} c_{95}$（阶段 `analysis:check-covariance`）。
- **汇总报告（`scripts/analysis/summarize.py`）**：校验 `limits`、`gridcheck`、`priorcheck` 及可选的 `covcheck`，确认各检验未改变特征数量，写出 `limits.csv`、`grid_vs_primary.csv`、`prior_vs_dense_grid.csv`、`fit-chi2.csv`（记录 $c=0$ 与联合网格最佳拟合点处的 $\chi^2$）、`limits.png`、`fits-alpha-<α>.png`、`report.json` 以及完全由数据驱动生成的 `RESULTS.md`（阶段 `analysis:summary`）。
- **论文输出（`scripts/analysis/paper.py`）**：在不做任何重新拟合的前提下，将论文所需的全部流水线资产汇总到阶段 `analysis:paper`：
  - `figures/`：8 幅确定性矢量 PDF 图（`radial_function`、`two_photon`、`signal_maps`、`template_maps`、`background`、`observational_uncertainty`、`data_vs_model`、`limits`）；
  - `numbers.tex` 与 `numbers.json`（模式 `analysis-paper-numbers-v1`）：不含阿拉伯数字的 `\newcommand{\Num...}{...}` LaTeX 宏定义，以及对应的未舍入原值和数据来源；
  - `table-limits.tex` 与 `table-chi2.tex`：观测上限 $c_{95}$ 表以及 `PAPER_ALPHA = 0.30` 处的 $\chi^2$ 完整 `tabular` 环境。
  展示参数（`PAPER_ALPHA = 0.30`、`PAPER_SNAPSHOT_LABEL = 2300`、EVPA 偏振短线参数、径向取样半径 `5, 10, 20` $r_g$、近邻像素对阈值 `15` $\mu\mathrm{as}$ 以及双发射源检验参数）全部记录在阶段清单中，径向窗口函数 $W(r) = \exp[-(r/r_w)^n]$ 则直接从 `RAPTOR` 工作副本（`--raptor-source-root` / `EHT_RAPTOR_SOURCE_ROOT`）的 `src/axion.c` 解析，该工作副本必须洁净地处于 `raptor.source.commit`（提交号即标识该文件，因此它不作为阶段输入）。径向表依次取自 `--radial-low` / `--radial-high`、[输入目录布局](data-layout.md#输入目录布局)或 `raptor/plan.json` 中记录的路径，且必须与计划中记录的 SHA-256 一致。
