# EHT 成像

## 用途

本页说明 EHT 成像输入：用 ehtim 1.3.2 复现 EHT Paper VII 的 EHT 2017 M87\* 偏振成像方法（其中
总强度流程沿用 Paper IV），随后对仪器偏振泄漏（D 项）解做蒙特卡罗，每天得到 1000 幅偏振图像。
这些图像按天组成的系综就是[分析](analysis.md)的 EHT 输入。EHT 成像不读取任何 RAPTOR 产品。

## 输入

- 第 95、96、100、101 天（2017 年 4 月 5、6、10、11 日）的四个 HOPS 低频段 uvfits 文件
  `SR2_M87_2017_<day>_lo_hops_zbl-dtcal+selfcal.uvfits`，位于目录 `EHT_IMAGING_INPUT_ROOT`。
  它们的文件名和 SHA-256 即 `eht_imaging.inputs`。
- 输入回执（`EHT_IMAGING_INPUT_RECEIPT`）：一个 JSON，含 `complete: true` 和 `files` 列表
  （`path`、`sha256`、`size`），必须与 `eht_imaging.inputs` 一致。
- 馈源旋转表（`EHT_IMAGING_FEED_ROTATION`）：一个 JSON，含 `status: verified`，并对每个台站给出
  `fr_par`、`fr_elev`、`fr_off`（度）。
- 官方 EHT 成像脚本（`EHT_IMAGING_OFFICIAL_SCRIPT`）。
- `configs/experiment.yaml`：`shared.image` 和 `eht_imaging` 节（[实验配置](experiment-config.md)）。

这四个位置是机器路径，在命令行或 `.env` 中给出。

## 输出

- I 阶段（每天一个）：固定的总强度图像 `final-fixedI.fits` 和经增益校准的原生圆偏振数据
  `clean-gaincal-pol.uvfits`。
- 基准阶段（每天一个）：每一轮的 P 图像和 D 项，以及最终 D 项。
- 精化阶段（每天一个）：精确 P 图像、精化后的 D 项（`final-dterms.json`）和证书。
- 打包（bundle）：四个 I 产品、四个精化 D 项文件和馈源表，以及 `bundle.json`。
- 蒙特卡罗批次：每次抽样一幅获证书的 P 图像（4 × 1000），附回执。
- 系综汇总：每天的 `stokes.npy`、`evpa.npy`、`global.npz`，以及 `axes.npz` 和 `report.json`。

文件细节见[数据布局](data-layout.md)。

## 运行方法

使用 `eht-imaging` 环境（ehtim 1.3.2、astropy）；`scripts/eht/summarize.py` 只需要 `eht-summary`
（astropy）。每条命令写出一个新的根目录。设置了数据根目录（`EHT_DATA_ROOT` 或 `--data-root`）时，下面每个根目录选项都有取自
[默认布局](data-layout.md#默认布局)的默认值；这里的命令显式给出了这些根目录。

```bash
# one root per day and stage; repeat for 95, 96, 100, 101
pixi run -e eht-imaging eht-reconstruct --stage I --day 95 --output-root I95 --plan
pixi run -e eht-imaging eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-imaging eht-reconstruct --stage fiducial --day 95 --i-root I95 \
    --output-root PD95 --review-file REVIEW
pixi run -e eht-imaging eht-refine --day 95 --source-root PD95 --i-root I95 --output-root D95

# after the four days
pixi run -e eht-imaging eht-export-bundle --i-root I95 --i-root I96 --i-root I100 --i-root I101 \
    --d-root D95 --d-root D96 --d-root D100 --d-root D101 --output-root BUNDLE
pixi run -e eht-imaging eht-sample --bundle-root BUNDLE --output-root MC --workers 64 --threads 2
pixi run -e eht-summary eht-summarize --sample-root MC --output-root ENSEMBLE
```

`scripts/eht/reconstruct.py --plan` 只校验输入、计算哈希并打印计划。实际执行需要 `--review-file`，即一个
记录外部审查决定的非空文件；它的字节成为阶段输入，而它的存在并不等于自动批准。

`scripts/eht/sample.py` 的选项：`--workers`（默认 64）、每个工作进程的 `--threads`（默认 2）、
`--max-active-memory-gib`（默认 100）、`--resource-wait`（默认 300 s）、`--task-limit`、
`--resume`。使用同一输出根目录的运行会续跑；它要求批次身份（打包、科学取值和代码）完全相同。

`scripts/eht/summarize.py` 直接读取 `scripts/eht/sample.py` 的输出根目录（见[系综汇总](#系综汇总)）。在另一台主机
上运行蒙特卡罗并把其输出根目录复制回来，属于流水线之外的操作事项。

## 细节

### 数据载入

`load_native`（`eht_axion/eht/calibration.py`）以原生圆偏振（RR、LL、RL、LR）载入 HOPS 文件，不
转换为 Stokes 量。只要四个乘积或其误差中有一个缺失或非正，整行就被丢弃。表中的馈源旋转系数写入
台站表，D 项置零，校准标志为 `ampcal = phasecal = frcal = True`、`dcal = False`：台站 AA、AP、SM、
JC 在 HOPS 产品中已做过泄漏改正，AZ、LM、PV 则没有。数据做一次相干的扫描平均。

### I 阶段

`scripts/eht/reconstruct.py --stage I`：

- 在原始文件上运行官方 EHT 成像脚本，源码中恰好做两处替换：傅里叶变换类型 `nfft` 改为 `direct`，
  视场 128 µas 改为 `shared.image.fov_uas`（120 µas）。脚本的最终图像 `im_out`（经其自身的恢复
  模糊之后）即固定 I 图像；它必须已经为正、有限，并位于 64 × 64、120 µas 的网格上。脚本的其他
  设置全部原样使用；特别是 0.6 Jy 的致密源流量密度（`zbl = 0.60`）固定在官方脚本内部，不是
  `configs/experiment.yaml` 的取值；
- 求解**共享台站增益**：在剔除了短于 `uv_min` 的基线的原生扫描平均数据上，以固定 I 图像为模型，
  对所有台站（包括 PV）做 ehtim 自校准，参数为 `method='both'`、增益容差
  `eht_imaging.gain_tol`（[0.02, 0.2]，官方脚本的取值），左右旋增益相等；
- 把这些增益施加到原生圆偏振数据上，写出 `clean-gaincal-pol.uvfits`、`shared-gains.npz` 和
  `polarization-input.json`。

从这里开始总强度图像固定不变；之后的所有步骤只求解线偏振 $P = Q + iU$（取 $V = 0$）。

### 基准 P/D 交替

`scripts/eht/reconstruct.py --stage fiducial` 读取同一天的 I 阶段。基础数据是施加了 `uv_min` 的增益校准
数据。从 $P = m_0\, I\, e^{2i\chi_0}$（`initial_m`、`initial_evpa_deg`）开始，`rounds`（50）轮中的
每一轮：

- 在 `restart_rounds` 所列轮次，从以 `restart_blur_uas`（20 µas）高斯模糊后的 I 加权偏振度重新
  启动 P，并乘以 $1-\epsilon$（`restart_domain_epsilon`）以保持 $\lvert P/I\rvert < 1$；
- 把系统误差 `fsys` 加到 RL/LR 误差上（`add_fsys`）；
- 用一次逆 Jones 改正施加当前 D 项（`apply_fixed_d`）；
- 用 ehtim `Imager` 对 P 成像（`pol='P'`，变换 `mcv`，数据项 `p_data`：`pvis` 和 `m`；正则项
  `p_reg`：`hw` 熵和 `ptv` 全变差；`norm_reg`；`maxit`；`stop`；直接傅里叶变换；固定 I 图像的
  流量）；启用 `shift_objective_constant` 时给目标函数加上一个常数，不改变极小点；
- 针对同一份未改正的含噪基础数据和新的 P 图像，用 ehtim 官方的前向 RL/LR 泄漏拟合
  （`fit_forward_d`）重新拟合 `solve_sites`（AZ、LM、PV）的 D 项，其中 D 先验及其梯度被替换为零。

最后一次只解 P 的求解使交付的图像与最后的 D 项一致。每一轮的 D 项以及原始图像和模糊图像都被写出；
`termination.json` 记录停止条件是固定的轮数预算。

### 基准 D 项的精化

`scripts/eht/refine.py` 用精确的凸 P 求解器（`eht_axion/eht/convex.py`）精化某一天的基准 D 项（精化阶段）。
在 I 和 D 固定时，正则项归一化的 ehtim P 目标函数（`norm_reg` 为 true；求解器只实现这种情形）是
$P$ 的凸函数：

$$
J(P) = \sum_{t\in\{\mathrm{pvis},m\}} \frac{w_t}{2N_t}\sum_j \frac{\lvert (A_t P)_j - d_{t,j}\rvert^2}{\sigma_{t,j}^2}
 + \frac{h_w}{F}\sum_x I_x\, h\!\left(\frac{\lvert P_x\rvert}{I_x}\right)
 + p_\mathrm{tv}\,\frac{\mathrm{TV}(P)}{\tau_\mathrm{TV}},
$$

其中权重 $w_\mathrm{pvis}$、$w_m$（`eht_imaging.p_data`）和 $h_w$、$p_\mathrm{tv}$
（`eht_imaging.p_reg`；取值 1、1、100、1），$h(m) = \tfrac12[(1+m)\ln(1+m) + (1-m)\ln(1-m)]$，
$F$ 为总流量，$\tau_\mathrm{TV}$ 为 ehtim 的 TV 归一化。$J$ 等于 ehtim 目标函数加上常数
$w_\mathrm{pvis} + w_m + h_w\ln 2$（有单元测试针对多组权重与 ehtim 比对）。数据元组就是官方 ehtim
成像器的数据元组，用相同的权重构造。`solve_problem` 运行一种对角预条件的原始–对偶混合梯度
（PDHG）方法，并每隔 `check_every` 次迭代计算一个可行的对偶下界；当最佳原始目标值减去最佳下界
不超过 `gap_atol` 时，返回的 P 即获得证书。这是浮点证书，不是区间算术。

从基准的最终 P 和 D 出发，`refine.outer_rounds`（100）轮中的每一轮把 P 求解到间隙容差
（`eht_imaging.refine` 的 `gap_atol`、`maxiter`、`check_every`、`time_limit`），并用 L-BFGS-B 选项
`eht_imaging.refine.d_optimizer` 通过 `fit_forward_d` 重新拟合 D；每次优化器调用都记录在
`D-optimizer-records.json` 中。最后一次 P 求解使用最后的 D。`termination.json` 记录
`all_P_numerically_certified`；不声称 P/D 联合最优。

### 打包导出

`scripts/eht/export_bundle.py` 读取四个 I 根目录和四个精化根目录，逐一校验（完整来源链），并检查：它们的
成像参数等于当前 `configs/experiment.yaml` 的 `eht_imaging_parameters`，每个精化阶段都源自同一天的 I 阶段，
每个精化阶段的每次 P 求解都获得证书，四个 I 阶段使用的是同一张馈源表。它以只读方式复制：

- 每天的 `fixedI.fits`、`gaincal-pol.uvfits` 和 `fiducial-dterms.json`；
- 来源链文件（`lineage/<day>/I-manifest.json`、`D-manifest.json`、`D-plan.json`、
  `D-termination.json`、`D-steps.json`）；
- `feed-rotation.json`；`bundle.json` 还记录代码的 git 提交。

`bundle.json` 列出每个文件的 SHA-256 和大小，并在模式 `eht-bundle-v2` 下记录科学约定：`fsys`、
`covariance_ddof`、`residual_component_sigma`（`residual_d_sigma`）、天空网格和 P 目标函数权重。

### D 项蒙特卡罗

`scripts/eht/sample.py` 校验打包（每个文件的哈希；`fsys`、`covariance_ddof`、`residual_component_sigma`
以及 P 目标函数权重与 `configs/experiment.yaml` 比对），然后：

- 构造 D 项分布（`d_distribution`）：对 `solve_sites` 中的每个台站和每种旋向 R、L，在四天的精化
  解上计算 $(\mathrm{Re}\,D, \mathrm{Im}\,D)$ 的均值和 2 × 2 协方差，`covariance_ddof` = 1；
- 每天以种子 `SeedSequence([seed, day])` 抽取 `uncertainty_draws` = 1000 组 D（`sample_dterms`）：
  对 `solve_sites`（AZ、LM、PV）用二元正态分布；对 `fixed_sites`（AA、AP、SM、JC）在每个实部和
  虚部分量上用标准差为 `residual_d_sigma` = 0.01 的零均值正态分布。热噪声不重新抽样；
- 对每次抽样：加上 `fsys`，一次性施加抽到的 D 项，并以 $P = 0.05\,I$（`initial_p_fraction`）为
  起点、用 `eht_imaging.monte_carlo`（`gap_atol`、`maxiter` 50000、`check_every`、`time_limit`
  600 s）求解精确的凸 P 问题（权重 `p_data`、`p_reg`）。未获证书的抽样有一次第二次尝试，参数为
  `retry_maxiter` 200000 和 `retry_time_limit` 1800 s；两次尝试都未获证书或两次出错的任务报告为
  已耗尽；
- 对每个获证书的抽样，写出原始图像和模糊图像（`best-raw-IQU.fits`、`best-blurred-IQU.fits`，
  模糊核为 `shared.image.blur_fwhm_uas`）以及一份证书。

批次身份是计划的 SHA-256（打包快照、科学取值、git 提交和程序包版本、执行线程数）；只有当全部
4000 次抽样都获得证书时才写出 `finish.json`。

### 系综汇总

`scripts/eht/summarize.py` 读取一个已完成的 `scripts/eht/sample.py` 输出根目录（`--sample-root`）。它要求
`batch.json` 带有自身的身份（批次计划的 SHA-256），其记录的科学取值等于 `scripts/eht/sample.py` 从当前
`configs/experiment.yaml` 派生的取值，并且 `finish.json` 具有相同的身份、对每次抽样列出一份获证书的回执。
它检查这些回执、证书以及每个被读取文件的 SHA-256，然后读取每天 1000 幅模糊图像（不再额外模糊），
并在 `day-<day>/` 下按天写出：

- `stokes.npy`：形状 (1000, 64, 64, 4)，轴依次为抽样、y、x、IQUV，采用原生 FITS 方向（东向偏移
  沿 x 减小，北向沿 y 增大）；
- `evpa.npy`：EVPA $\tfrac12\arg(Q+iU)$，以 $\pi$ 为模，单位为弧度；在 $I \le 0$ 或 $Q+iU=0$ 处为
  NaN；
- `global.npz`：每次抽样的 EHT Paper VII 第 5.3 节图像平均量（`global_diagnostics`，列为
  `total_I`、`m_net`、`m_resolved`、`evpa_rad`、`beta2_real`、`beta2_imag`），按整个 120 µas
  视场上的像素求和计算：
  - $\lvert m\rvert_\mathrm{net} = \lvert\sum_i P_i\rvert/\sum_i I_i$（Paper VII 式 12）；
  - $\langle\lvert m\rvert\rangle = \sum_i\lvert P_i\rvert/\sum_i I_i$（式 13；列 `m_resolved`）；
  - EVPA $\tfrac12\arg\sum_i P_i$；
  - $\beta_2 = \sum_i P_i\,e^{-2i\varphi_i}/\sum_i I_i$（式 14），$\varphi$ 为绕网格中心、从北向东
    量度的位置角。论文取 $\rho_\mathrm{min} = 0$，$\rho_\mathrm{max}$ 覆盖整个视场，因此
    $I_\mathrm{ring}$ 即总流量；代码的做法相同。

  所用图像为模糊图像（`shared.image.blur_fwhm_uas` = 19 µas，即 Paper VII 第 5.1 节对 eht-imaging
  使用的核）。与论文唯一的不同是图像中心：论文把每幅图像定中心于使其模糊后的 Stokes I 与 Paper IV
  一致图像（consensus image）的互相关取极大的像素偏移处；这里的中心是网格中心。这些量用于与
  Paper VII 比较，分析部分不使用它们。

`axes.npz` 存放 `east_uas` 和 `north_uas`；`report.json` 记录所用约定。

### 图像中心

EHT Paper VII 第 5.2 节通过互相关对齐图像。在本复现中这一对齐是零平移：图像按重建结果直接用于
公共网格，分析部分没有重新定中心的选项。

### 与 EHT Paper VII 的已知偏差

- 这是方法上的复现，不是对已发表图像的逐比特复现。官方 I 脚本以直接傅里叶变换和 120 µas 视场
  运行。
- 本地求解的 D 项与 Paper VII 表 5 不同。最大的差异出现在 2017 年 4 月 10 日（第 100 天）的 PV，
  其 $D_R$ 相差约 19 个百分点。代码不与表 5 做比较。
- `restart_rounds`、`maxit` 和精化预算是数值上的选择，不是论文中的取值。
- $\beta_2$（以及所有分析特征）的图像中心是网格中心，而不是 Paper VII 第 5.2/5.3 节的互相关中心
  （见[图像中心](#图像中心)）。

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [实验配置](experiment-config.md)
- [分析](analysis.md)
- [来源记录](provenance.md)
- [数据布局](data-layout.md)
