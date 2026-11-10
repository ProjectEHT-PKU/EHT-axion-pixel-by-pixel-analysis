# 分析

## 用途

本页说明把两个输入结合起来的分析部分：它把交付的 RAPTOR Stokes 数据立方转换为 EVPA 响应模板，
在 EHT 系综上选取似然特征，构造日对 EVPA 差及其协方差，并计算无量纲耦合常数 $c$ 的、对相位边缘化
的后验上限 $c_{95}(\alpha)$，同时给出零数据对照（zero-data control）和数值检验。

## 输入

- 交付的 RAPTOR 数据根目录（一次交付或合并的 `data/`），含 `raptor/res-100/alpha-<α>/` 和
  `raptor/res-400/background/`（[RAPTOR 生产](raptor-production.md)）。
- `eht_summarize.py` 的 EHT 系综阶段（[EHT 成像](eht-imaging.md)）。
- `experiment.yaml`：`shared`、`analysis`，以及来自其他节的 `raptor.camera_fov_rg`、
  `raptor.snapshots`、`raptor.response`、`raptor.background` 和 `eht_imaging.uncertainty_draws`
  （[实验配置](experiment-config.md)）。

## 输出

- 模板根目录：预处理后的数据立方、逐快照的响应拟合、模板（1D 和 2D）、预处理后的无轴子帧。
- 上限根目录：选定的观测量、GRMHD 背景、2 种维度 × 3 种协方差模型 × {观测数据, 零数据对照} 的
  上限、`limits.csv`、一幅图和 `report.json`。
- 检验根目录（网格、先验、协方差）、一个汇总根目录和论文根目录。

文件细节见[数据布局](data-layout.md)。

## 运行方法

使用默认环境。设置了数据根目录（`M87_AXION_DATA_ROOT` 或 `--data-root`）时，下面每个根目录选项都有取自
[默认布局](data-layout.md#默认布局)的默认值；这里的命令显式给出了这些根目录。

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
    --raptor-plan plan.json --raptor-run-root RUN_ROOT --output-root PAPER
```

上面的覆盖值只是示例；检验脚本接受任何与 `analysis.inference` 不同的取值。对
`analysis_summarize.py` 而言，`--covcheck` 是可选的。

## 细节

### 图像预处理

`analysis_templates.py` 按当前实验校验每个交付的数据立方（`generated_inputs`：身份记录、形状、
标签等于 `raptor.snapshots`），并逐帧运行 `preprocess_cube`（`eht_axion/analysis/templates.py`）。
对每个 Stokes 通道，`preprocess_frame`：

- 与标准差为 $\sigma = \mathrm{FWHM}/\sqrt{8\ln 2}$ 的高斯核卷积，其中 FWHM =
  `shared.image.blur_fwhm_uas`（19 µas）换算成 $r_g$（2.11 $r_g$），边界按最近值延拓；
- 以双线性插值把相机视场 `raptor.camera_fov_rg` 映射到 `shared.image.shape` 像素、
  `shared.image.fov_uas`（31.42 $r_g$）的天空网格上，并旋转 `analysis.rotation_degrees`
  $\theta$ = 108°；
- 旋转偏振：$Q + iU \to e^{2i\theta}(Q + iU)$。

不做流量重新归一化。模拟图像的第 0 轴对应西向（东向偏移减小），第 1 轴对应北向，与转置后的 EHT
图一致（见[特征选择](#特征选择)）。响应数据立方保留全部 17 个状态；背景数据立方只保留状态 0
（`analysis:background`）。

### 响应拟合

`fit_processed` 去掉状态 0，并对每个快照，用 `fit_phase_dft` 对每个像素（2D）和每个孔径扇区
（1D，见下一节）拟合各相位状态（轴子相位 $\psi_k = 2\pi k/16$）的 16 个 EVPA 值 $\chi_k$：

$$
\chi_k = \bar\chi + A\cos(\psi_k + \phi).
$$

这些值先被提升到其圆均值所在的分支（周期 $\pi$），然后
$A e^{i\phi} = \tfrac{1}{8}\sum_k (\chi_k - \bar\chi)\, e^{-2\pi i k/16}$。缺少某个相位的位置无效；
振幅为零时相位无定义。每个 $\alpha$ 的输出为 `response-2d.npz` 和 `response-1d.npz`，含形状为
(141, …) 的数组 `amplitude`、`phase`、`baseline`（弧度）以及 `snapshot_labels`。

### 模板

`summarize_fits` 对每个 $\alpha$ 和每个位置合并 141 个快照：$A$ 的均值和总体标准差（ddof = 0），
以及 $\phi$ 的圆均值和 $\phi$ 围绕该均值的标准差。一个位置至少需要两个有限值的快照。
`templates.npz` 存放 `alphas`，以及 1D 形状为 (8, 360)、2D 形状为 (8, 64, 64) 的数组。当四个统计量
都有限且对每个 $\alpha$ 标准差都非负时，`load_templates` 把该位置标记为有效。信号模型只使用均值；
标准差只进入这个有效性掩码。

### 孔径观测量

1D 观测量是扇形（“披萨片”）孔径内强度加权偏振的 EVPA：

$$
\chi_\mathrm{ap}(\varphi_k) = \tfrac12 \arg \int_0^{r_\mathrm{max}}\!\!\int_{\varphi_k-w/2}^{\varphi_k+w/2}
 I\,(Q + iU)\; r\,\mathrm{d}\varphi\,\mathrm{d}r ,
$$

其中 $w$ = `analysis.aperture.window_degrees`（10°），$r_\mathrm{max}$ = `radial_max_fraction` ×
半视场（60 µas），360 个扇区中心相隔一度。`aperture_operator`（`eht_axion/analysis/aperture.py`）
是一个由双线性求积权重构成的稀疏矩阵（径向步长为视场单位的 0.02，角向步长 0.2°，端点用梯形权重）。
同一个算子作用于模拟图像，并在转置为相同轴顺序后作用于每次 EHT 抽样，因此模型和数据使用完全相同的
观测量。施加算子后，分格循环移位 90，使第 $k$ 格的中心位于北偏东 $k$ 度。

记扇区 $k$ 与像素 $j$ 的求积权重为 $W_{kj}$，代码计算：

- 扇区 EVPA $\chi_k = \tfrac12\arg\sum_j W_{kj}\, I_j (Q_j + iU_j)$（`integrate_stokes`；即采用
  这一求积的 EHT Paper VII 式 11）；
- 扇区强度 $B_k = \sum_j W_{kj}\, I_j$；
- 强度加权的平均偏振亮度
  $\langle\mathcal{P}\rangle_k = \sum_j W_{kj}\, I_j \lvert P_j\rvert / \sum_j W_{kj}\, I_j$，其中
  $\lvert P\rvert = \sqrt{Q^2+U^2}$（`sector_brightness`；即采用这一求积的 Paper VII 式 10）。

$B_k$ 和 $\langle\mathcal{P}\rangle_k$ 使用与观测量相同的权重和相同的强度加权；它们只用于选取扇区。

### 信号模型

对一个日对 $p$，特征 $x$ 处的模型 EVPA 差对 $c$ 是线性的：

$$
s_p(x;\alpha,\phi_0) = b(\alpha)\,\frac{2\sin(\omega\tau/2)\,\mathrm{sinc}(\omega T/2)}{N_\mathrm{ref}}\,
 A(x,\alpha)\,\sin\!\big(\phi_0 + \phi(x,\alpha) + \omega t_p\big),
$$

其中

- $\omega = \Omega(\alpha)/t_g$，单位 rad h$^{-1}$，$\Omega(\alpha) = \alpha(1-\alpha^2/(2n^2))$；
- $\tau$ = `day_gap_hours`（24 h），$T$ = `exposure_hours`（6 h），
  $\mathrm{sinc}(u) = \sin u / u$（矩形曝光平均）；
- $t_1 = 0$，$t_2$ = `pair_separation_hours`（120 h）：第一个日对的中点为相位原点；
- $N_\mathrm{ref}$ = `raptor.response.reference_norm`；
- $b(\alpha) = 0.03\,(\alpha/0.05)^{1.75}$（`analysis.cloud`）；
- $A$ 和 $\phi$ 是模板均值，在模板节点之间对 $\alpha$ 线性插值（$\phi$ 在单位圆上插值），在节点
  范围之外保持常数；
- $\phi_0$ 是未知的轴子相位，取 $[0, 2\pi)$ 上 `phase_count`（64）个值的均匀网格。

这对应 `eht_axion/analysis/model.py` 中的 `evpa_templates` 和 `cloud_amplitude_scaling`。

### 特征选择

`analysis_limits.py` 读取每天的模糊图像系综，在不看日差大小的前提下选取特征。EHT 图 $[y, x]$ 被
转置为模板顺序 $[x, y]$；天空坐标轴必须等于公共网格。

2D 特征是同时满足下列全部条件的模板像素：

- 模板有效；
- 在每一天，系综平均 Stokes 图像满足 $I \ge$ `intensity_fraction`（0.1）× 其最大值，且
  $\lvert P\rvert = \sqrt{Q^2+U^2} \ge$ `polarization_fraction`（0.2）× 其最大值，并且 EVPA 在每次
  抽样中都有限；
- 到网格中心的距离介于 `min_radius_rg`（2 $r_g$ = 7.6 µas）和孔径半径（60 µas）之间；
- 两个像素下标都是 `pixel_stride`（5）的倍数。

1D 特征是满足 $k \bmod$ `sector_spacing_degrees` $=$ `sector_offset_degrees`（4°、14°、…、354°）
的扇区 $k$，并要求：模板有效；扇区 EVPA $\chi_k$ 在每一天的每次抽样中都有限；并且在每一天，$B_k$
和 $\langle\mathcal{P}\rangle_k$ 的抽样均值（见[孔径观测量](#孔径观测量)）满足
$\bar B_k \ge$ `intensity_fraction` $\times \max_k \bar B_k$ 和
$\overline{\langle\mathcal{P}\rangle}_k \ge$ `polarization_fraction` $\times
\max_k \overline{\langle\mathcal{P}\rangle}_k$。这三个统计量都由每次 EHT 抽样经孔径算子计算。

然后 `build_observations` 只保留那些在全部四天中、各抽样的圆合成长度
$R = \lvert\langle e^{2i\chi}\rangle\rvert$ 都不小于 `resultant_min`（0.7）的特征。

### 观测量与协方差

对第 $d$ 天和保留的特征，`moments` 给出圆均值 $\bar\chi_d = \tfrac12\arg\langle e^{2i\chi}\rangle$
以及各抽样围绕它的总体协方差（ddof = 0），

$$
C_d = \frac{1}{N}\sum_{k=1}^{N} r_k r_k^\mathsf{T},\qquad
r_k = \mathrm{wrap}_\pi(\chi_k - \bar\chi_d) - \langle\mathrm{wrap}_\pi(\chi - \bar\chi_d)\rangle .
$$

$C_d$ 是单次抽样的散布（D 项校准不确定度），不是均值的协方差。数据是两个日对差
$d_1 = \mathrm{wrap}_\pi(\bar\chi_{96} - \bar\chi_{95})$ 和
$d_2 = \mathrm{wrap}_\pi(\bar\chi_{101} - \bar\chi_{100})$，协方差为 $C_1 = C_{95} + C_{96}$ 和
$C_2 = C_{100} + C_{101}$（各天独立，两个日对之间无协方差）。每个 $C_p$ 都必须严格正定；不做正则化。
分支接近程度（以 $\sigma$ 为单位的 $\pi/2 - \lvert d\rvert$）和特征值诊断会被报告，但不改变选择。

### GRMHD 背景

`selected_background` 恰好在所选特征上估计吸积流本身的 EVPA 变化。它从预处理后的无轴子帧（141 个
快照、状态 0、重采样到天空网格的 400 像素数据立方）中取出每个所选特征的 EVPA $\theta_t$（2D：像素
EVPA；1D：同一个孔径算子），计算滞后差分 $\Delta_t = \mathrm{wrap}_\pi(\theta_{t+L} - \theta_t)$，
其中 $L$ = `analysis.background.lag_frames`（3），再对所有特征和帧对合并求总体标准差：

$$
\sigma_\mathrm{bg} = \sqrt{\langle(\Delta - \langle\Delta\rangle)^2\rangle}.
$$

帧间隔取自交付数据立方的 `generated.json` 中记录的 GRMHD 尾部时间（1 $t_g$ = 8.89 h，因此 $L$ 对应
26.7 h，接近一个日对中两天的间隔）。$\sigma_\mathrm{bg}$ 直接描述一个差值，因此不乘 $\sqrt2$ 因子。

### 协方差模型

| 模式 | 似然中使用的日对协方差 |
|---|---|
| `std` | $\mathrm{diag}(C_p)$ |
| `cor` | $C_p$ |
| `cor_fixed_background` | $C_p + \sigma_\mathrm{bg}^2\,\mathbb{1}$ |

没有任何模式包含与 $c^2$ 成正比的响应模型方差。

### 似然与相位边缘化

对每个 $\alpha$，`calculate_observed`（`eht_axion/analysis/likelihood.py`）在 $c$ 和 $\phi_0$ 的
网格上计算

$$
\chi^2(c,\phi_0) = \sum_{p=1}^{2} \big(d_p - c\,s_p(\phi_0)\big)^\mathsf{T} C_p^{-1} \big(d_p - c\,s_p(\phi_0)\big),
$$

两个日对共用同一个相位 $\phi_0$，并用均匀离散先验对相位边缘化：

$$
\ln L(c) = \ln\frac{1}{N_\phi}\sum_{j=1}^{N_\phi} e^{-\chi^2(c,\phi_{0,j})/2}.
$$

与 $c$ 和 $\phi_0$ 无关的常数（高斯分布的对数行列式和归一化）被略去。残差 $d_p - c\,s_p$ 不再做
回绕：似然是数据的 EVPA 主值分支上的高斯分布，而不是回绕正态分布。

### 先验与后验

$c$ 的网格有 `coupling_points`（120）个对数间隔的取值，从 `coupling_min`（1/137）到
`coupling_max`（1000）。`prior: log-uniform` 时，先验密度在该区间上为 $\pi(c) \propto 1/c$
（`uniform`：常数）。`normalized_posterior_cdf`（`eht_axion/analysis/posterior.py`）对 $c$ 用梯形
法则积分 $L(c)\pi(c)$；`posterior_upper_limit` 在 `credibility`（0.95）处对同一分段线性密度的累积
分布求逆，得到 $c_{95}$。结果以有限网格和先验为条件（即条件贝叶斯上限）；边界诊断（第一个和最后
一个网格区间内的概率质量、峰值位于边界、上限位于最后一个区间）会被报告。

### 零数据对照

对每种情形，用 $d_1 = d_2 = 0$ 以及相同的协方差和特征重复同样的计算（文件名中为 `zero-data`）。
它给出仅凭协方差和模板会得到的上限，用于与观测上限比较。

### 检验

- `analysis_check_grid.py` 用覆盖后的 `phase_count`、`coupling_points` 和/或 `coupling_max` 对一个
  上限根目录中保存的观测量重新计算；其余一切必须等于产生源目录的那次实验。它报告相对于源目录的
  $\Delta\log_{10}c_{95}$。用更密的网格运行一次（“网格检验”），用更宽的 $c$ 范围再运行一次（“先验
  检验”）。
- `analysis_check_covariance.py` 把每天的抽样分成前一半和后一半（各 500 次），围绕固定的全系综均值
  重新计算 $C_d$，并对全部六种情形重复观测上限的计算；它报告相对于主结果的
  $\Delta\log_{10}c_{95}$。特征和数据保持不变。

### 汇总

`analysis_summarize.py` 读取已完成的上限、网格检验、先验检验阶段以及可选的协方差检验阶段，并写出：
`limits.csv`；网格检验与主结果、先验检验与网格检验的比较（二者都必须保持相同的特征数）；每种观测
情形、每个 $\alpha$ 保存下来的联合网格最佳拟合，以及 $c = 0$ 和最佳网格点处的 $\chi^2$
（`fit-chi2.csv`，每个 $\alpha$ 一幅图 `fits-alpha-<α>.png`）；一幅上限图；`report.json`；以及
`RESULTS.md`。`RESULTS.md` 由配置的 `analysis.inference.credibility` 和本次运行的数字（上限表、
$\sigma_\mathrm{bg}$、各检验中的最大变化、拟合表）生成；其中不含任何固定的结论。

### 论文输出

`analysis_paper.py`（阶段 `analysis:paper`）把论文从流水线取用的全部内容收集到一个目录中。它读取
上限、三个检验阶段、模板（论文所用 $\alpha$ 的预处理数据立方、响应拟合和模板，以及预处理后的无轴子
背景）、EHT 系综、交付的 400 像素无轴子数据立方、RAPTOR 计划以及编译目录中的 `axion.c` 副本。径向表
是计划中记录的文件（`--radial-low`、`--radial-high` 或 `EHT_RADIAL_LOW`、`EHT_RADIAL_HIGH` 可覆盖其
路径；SHA-256 必须等于计划中的记录），场窗口 $W(r)=\exp[-(r/r_w)^n]$ 从 `axion.c` 中读取。上限必须
由当前的 `experiment.yaml` 产生。不做任何重新拟合。输出目录包含：

- `figures/`：八幅矢量 PDF，即 `radial_function`、`two_photon`、`signal_maps`、`template_maps`、
  `background`、`observational_uncertainty`、`data_vs_model` 和 `limits`；
- `numbers.tex`：论文中引用的每个数字对应一个 `\newcommand`，命名为 `\Num` 后接不含数字的英文单词
  （例如 `\NumLimitObsSectorsCorBgAlphaThirty`），每个之前有一行注释给出其说明和来源；
- `numbers.json`（模式 `analysis-paper-numbers-v1`）：每个宏的 TeX 文本、未经舍入的数值和来源；
- `table-limits.tex` 和 `table-chi2.tex`：完整的 `tabular` 环境，分别给出观测上限 $c_{95}$（三位
  有效数字），以及 $c=0$ 和最佳网格点处的 $\chi^2$。

派生量：模拟流量是交付的 400 像素无轴子图像的 $\sum I\,\Delta\Omega$，其中 $\Delta\Omega$ 由
`raptor.camera_fov_rg` 和 $r_g$ 的角大小得到；重建图像的流量、偏振度和 $\lvert\beta_2\rvert$ 是系综
`global.npz` 中的取值（按天平均）；模拟图像的偏振度和 $\lvert\beta_2\rvert$ 对转置到 EHT 方向的
预处理无轴子帧施加 `global_diagnostics` 得到；背景图中的灰色直方图合并了每个扇区以及每个 EVPA 有定义
的像素的滞后差分。展示上的选择是脚本中的常量，并记录在清单中：`PAPER_ALPHA` = 0.30 和
`PAPER_SNAPSHOT_LABEL` = 2300（响应图、模板图、$\chi^2$ 表和拟合图）、EVPA 短线设置、径向轮廓、窗口和
轨道周期的取值半径、近邻像素对的 15 µas 间隔，以及两发射体检验所用的点。

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [实验配置](experiment-config.md)
- [EHT 成像](eht-imaging.md)
- [RAPTOR 生产](raptor-production.md)
- [数据布局](data-layout.md)
- [术语表](glossary.md)
