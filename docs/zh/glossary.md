# 术语表

## 用途

本页定义文档和代码中使用的符号和术语，给出其单位以及（若有）对应的 YAML 键或代码名称，使其他页面
可以直接使用而不必重复定义。

## 细节

### 物理量

| 术语 | 定义 |
|---|---|
| $t_g$ | 引力时间 $GM/c^3$；对 $M = 6.5\times10^9\,M_\odot$ 为 8.893 h（`PhysicalScale.tg_hours`） |
| $r_g$ | 引力半径 $GM/c^2$；在 16.8 Mpc 处为 3.819 µas（`PhysicalScale.angular_rg_uas`） |
| $\alpha$ | 质量为 $\mu$ 的轴子的引力精细结构常数 $GM\mu/(\hbar c)$（`shared.alphas`）；对 M87\* 的质量，$\mu c^2 = \alpha \times 2.06\times10^{-20}$ eV |
| 211 态 | $n = 2$、$l = 1$、$m = 1$ 的轴子云束缚态 |
| $\Omega(\alpha)$ | 无量纲频率 $M\omega = \alpha(1-\alpha^2/(2n^2))$，$n$ = `shared.principal_n`（`physics.omega`）；$\omega = \Omega/t_g$ |
| $c$（耦合常数） | 标定轴子引起的 EVPA 信号大小的无量纲模型耦合常数；在 `likelihood.py` 中称为 `g` 和 `couplings`。代码不把它换算成以 GeV$^{-1}$ 为单位的轴子–光子耦合 |
| $c_{95}$ | $c$ 在 `analysis.inference.credibility`（0.95）处的后验上限 |
| $b(\alpha)$ | 轴子云振幅标度 `analysis.cloud.coefficient` $\times(\alpha/$`pivot_alpha`$)^{\text{exponent}}$ |
| $N_\mathrm{ref}$ | RAPTOR 轴子归一化 `raptor.response.reference_norm`；模板以该归一化为单位 |
| 轴子相位 | RAPTOR 相位序号 $k$，相位为 $k\pi/8$，$k = 0,\dots,15$；在似然中为未知相位 $\phi_0$ |

### 偏振测量

| 术语 | 定义 |
|---|---|
| Stokes I、Q、U、V | 强度和偏振；RAPTOR 数据立方单位为 Jy/sr，EHT FITS 的单位见其 `BUNIT` |
| $P$ | 线偏振 $Q + iU$；$\lvert P\rvert = \sqrt{Q^2+U^2}$ |
| EVPA | 电矢量位置角（electric-vector position angle）$\chi = \tfrac12\arg(Q+iU)$，单位为弧度，以 $\pi$ 为模定义（`evpa.finite_evpa`） |
| $\mathrm{wrap}_\pi$ | 回绕到 $[-\pi/2, \pi/2)$ 的差值（`evpa.circular_difference`） |
| 合成长度（resultant length） | 各抽样上的 $R = \lvert\langle e^{2i\chi}\rangle\rvert$；EVPA 全部相同时为 1 |
| $\lvert m\rvert_\mathrm{net}$、$\langle\lvert m\rvert\rangle$ | 净偏振度和平均偏振度 $\lvert\sum P\rvert/\sum I$ 与 $\sum\lvert P\rvert/\sum I$（Paper VII 式 12、13；列 `m_net`、`m_resolved`） |
| $\beta_2$ | 整个视场上的 $\sum P\,e^{-2i\varphi}/\sum I$，$\varphi$ 为绕网格中心、从北向东量度的位置角（Paper VII 式 14） |

### EHT 成像

| 术语 | 定义 |
|---|---|
| 第 95、96、100、101 天 | 2017 年的年积日：4 月 5、6、10、11 日 |
| 日对 | (95, 96) 和 (100, 101)；数据是差值 96−95 和 101−100 |
| D 项（D terms） | 每个台站的复数仪器偏振泄漏 $D_R$、$D_L$（小数） |
| 求解台站、固定台站（solve sites、fixed sites） | AZ、LM、PV（此处求解其 D 项）和 AA、AP、SM、JC（在 HOPS 数据中已预先改正） |
| fsys | 加到 RL 和 LR 误差上的相对系统误差（`eht_imaging.fsys`） |
| 固定 I（fixed I） | I 阶段的总强度图像；之后的所有步骤只求解 $P$ |
| 基准（fiducial） | `scripts/eht/reconstruct.py --stage fiducial` 的 P/D 交替结果 |
| 精化 D 项 | `scripts/eht/refine.py` 处理后的基准 D 项（精化阶段的 `final-dterms.json`） |
| 模糊图像 | 与 FWHM 为 `shared.image.blur_fwhm_uas`（19 µas）的圆高斯核卷积后的图像；文件名为 `*-blurred-IQU.fits` |
| 抽样（draw） | 一组 D 项样本及其获证书的 P 图像；每天 1000 次 |
| 间隙证书（gap certificate） | 精确 P 求解器的最佳原始目标值减去最佳对偶下界，不超过 `gap_atol` |

### 分析

| 术语 | 定义 |
|---|---|
| 状态（state） | RAPTOR 数据立方第 1 轴上的序号：0 为无轴子，$k+1$ 为轴子相位 $k\pi/8$ |
| 快照（snapshot） | 一个 GRMHD 时刻；共 141 个快照，相隔 1 $t_g$ |
| 响应（response） | 逐快照在 16 个相位上的一次谐波拟合 $\bar\chi + A\cos(\psi+\phi)$ |
| 模板（template） | 响应振幅 $A$ 和相位 $\phi$ 对快照的平均 |
| 1D、2D | 孔径扇区观测量（360 个一度宽的扇区）和像素观测量（64 × 64 网格） |
| 特征（feature） | 用于似然的一个选定扇区（1D）或像素（2D） |
| GRMHD 背景 | $\sigma_\mathrm{bg}$，所选特征上滞后无轴子帧的 EVPA 散布 |
| 零数据对照（zero-data control） | 把日差设为零后的同一似然（文件名中为 `zero-data`） |

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [分析](analysis.md)
- [EHT 成像](eht-imaging.md)
