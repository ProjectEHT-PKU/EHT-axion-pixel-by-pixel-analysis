# 四日联合受控研究：已完成产物汇总

## 结论先行

这是**冻结 EHT 条件下的模型比较，不是科学正确性认证或覆盖率校准**。八个 α 全部完成；主表保留 672 行
（每 α：2 维度 × (H2 的 12 个模型 + H3 的 9 个模型) × observed/NULL = 84）。
不按约束强弱挑模型，不把新全模型自动视为更正确。

α=0.3 的 2D joint empirical B,T：H2=1.76224，
H3=0.313781；而 H3 joint Gaussian=1.8054。
从约 1.7622 到 0.31378 的强收紧，面对 Gaussian 约 1.8054 的结果，首先是**抽样分布/模型选择巨大敏感性警报**，
不是应当选用最强上限的理由。下面纳入已完成的有限窗口与峰处 ESS 诊断，但这些检查不等于该联合统计模型已获论文级验证。

## 冻结了什么、没有做什么

四日实测数据、特征选择（1D N=23；2D N=20）、逐日 EHT Monte Carlo 协方差、I/error 与超参数全部冻结。
未重建 EHT、未补充 D/latent-EHT 不确定性、未做 RT、未调 χ²、未缩放响应、未改 canonical NULL。
NULL 是同一条件模型中的 matched_zero 对照，不是新观测，也不是“最强模型”的选择依据。
所有 c95 都是原代码耦合参数 c 的条件 95% posterior 上限；**不另造 GeV 单位**。
先验一致为 [1/137,1000] 上 log-uniform；主计算 coupling/phase=120/64，α=.3 网格对照为 240/128。
遗留 H2 原六项（1D/2D × observed 模型控制）按原网格恢复；α=.3 legacy c95 为 1d: 1.7979244, 2d: 1.3827718。
这是实际受控结果，不为追求更强约束而调参。

H2=[day2−day1,day4−day3]；H3 再加入 day3−day2。
四日最多提供三个独立时间差组合，**不是六个独立两两差**。
观测协方差严格按 H diag(C_daily) Hᵀ 构造；共享日产生的负交叉块保留。
H3/H2 在相同模型假设下的改变也不能全归因于“信息量增加”：新增时间组合同时接受同一时间状态模型约束。
144 个 observed 同模型 H3/H2 比值中 68 个 >1（上限变弱）；最大 5.62394
出现在 1d, α=0.44, independent_gaussian_moments，图中不隐藏这些情况。

## α=0.3 全因子表（observed c95）

| 模型（固定构造顺序） | 1D H2 | 1D H3 | 2D H2 | 2D H3 |
| --- | --- | --- | --- | --- |
| mean+iso | 1.8594 | 1.53315 | 1.41613 | 1.43955 |
| mean+diag | 1.91617 | 1.50703 | 1.62033 | 1.51135 |
| mean+full | 1.52504 | 1.37074 | 1.62271 | 1.68182 |
| mean+full+offset | 1.35735 | 1.20878 | 1.53601 | 1.67853 |
| frozen mixture | 1.62045 | 1.71394 | 1.84858 | 1.79014 |
| evolving mixture | 1.45852 | 1.56274 | 1.8079 | 1.78344 |
| joint empirical B,T | 1.90728 | 3.34672 | 1.76224 | 0.313781 |
| joint Gaussian | 0.935208 | 1.16749 | 1.7349 | 1.8054 |
| indep Gaussian | 1.07436 | 1.76656 | 1.96262 | 1.89366 |

mean+iso/diag/full 在同一 mean response 下只改变背景预测协方差；ISO 方差池固定为 H2，
H3 也使用该池，避免加入第三差时顺便改变噪声定义。
frozen/evolving mixture 使用完整复杂响应状态，不将 A 与 φ 独立抽样：
每个状态为复数 (α,state,day,feature)，同一个全局相位与整状态作用于所有特征和时间差。
冻结响应 141 状态，演化响应 124 重叠窗口，先验均为均匀状态权重；不能把这些状态当独立观测或用 √N 缩小误差。
这里称 response，不将其改名为 π/FDT。

背景默认采用**零均值模型**：从经验背景减去样本均值；这与固定 offset 都是估计模型选择，不是校准改正。
响应均值也由有限状态估计；均值响应和响应混合是不同模型，不是消除真实系统误差的保证。
以下单列固定背景均值的控制，范数有限且非零：

| 维度 | H | 零均值 c95 | 固定 offset c95 | offset 范数 (rad) |
| --- | --- | --- | --- | --- |
| 1d | 2 | 1.52504 | 1.35735 | 0.0308347 |
| 1d | 3 | 1.37074 | 1.20878 | 0.0893169 |
| 2d | 2 | 1.62271 | 1.53601 | 0.0279344 |
| 2d | 3 | 1.68182 | 1.67853 | 0.0833897 |

## 分离 B/T 交叉协方差，而不是混淆分布改变

joint empirical 使用配对的 B+c Re(exp(iψ)T) 混合，仅加 C_EHT，**不再加 C_bg**，避免双计数。
joint Gaussian 在每个 c、ψ 对同一和量矩匹配，保留线性 c 的 B/T 交叉项、c² 项及参数相关 logdet。
它与 empirical 的区别既有分布形状也有依赖结构的表现，不能只称“加入 cross covariance 的收益”。
真正隔离 cross covariance 的对照是 joint Gaussian 对 independent Gaussian：后者边缘矩完全相同，仅令 B/T cross block 为零。

| 维度 | H | joint Gaussian | independent Gaussian | joint/independent |
| --- | --- | --- | --- | --- |
| 1d | 2 | 0.935208 | 1.07436 | 0.87048 |
| 1d | 3 | 1.16749 | 1.76656 | 0.660885 |
| 2d | 2 | 1.7349 | 1.96262 | 0.883969 |
| 2d | 3 | 1.8054 | 1.89366 | 0.953394 |

背景协方差来自 124 重叠窗口的预测矩（ddof=0，不除以窗口数）；H2/H3 向量维数为 2N/3N（1D 46/69，2D 40/60）。
full C_bg 秩范围 40–69；
加入冻结 C_EHT 后固定协方差最小特征值全表最小 1.59457e-09，为 SPD。
没有任意 ridge/shrinkage，也没有因模型拟合好坏而更换协方差。

## 已完成的有限窗口敏感性（α=.3）

前/后各 62 个窗口分别重估本子集的矩；不是独立观测误差条，也不是增加新数据。保留 H2/H3 及 Gaussian 背景控制，不能只检查 empirical 模型。

| 维度 | H | 模型 | 全窗口 | 前62 | 后62 | 最大 \|Δdex\| 对全窗口 |
| --- | --- | --- | --- | --- | --- | --- |
| 1d | 2 | mean_background_full | 1.52504 | 1.31888 | 2.00804 | 0.119493 |
| 1d | 2 | evolving_mixture_full | 1.45852 | 1.3037 | 2.2236 | 0.183144 |
| 1d | 2 | joint_background_evolving_mixture | 1.90728 | 2.10909 | 1.98609 | 0.0436803 |
| 1d | 2 | joint_gaussian_moments | 0.935208 | 0.521187 | 1.60683 | 0.253915 |
| 1d | 3 | mean_background_full | 1.37074 | 0.828654 | 1.80954 | 0.218581 |
| 1d | 3 | evolving_mixture_full | 1.56274 | 0.717028 | 2.1158 | 0.33835 |
| 1d | 3 | joint_background_evolving_mixture | 3.34672 | 2.45226 | 2.09439 | 0.203561 |
| 1d | 3 | joint_gaussian_moments | 1.16749 | 0.702009 | 1.88199 | 0.220912 |
| 2d | 2 | mean_background_full | 1.62271 | 1.29692 | 1.63831 | 0.0973291 |
| 2d | 2 | evolving_mixture_full | 1.8079 | 1.32615 | 1.99938 | 0.134582 |
| 2d | 2 | joint_background_evolving_mixture | 1.76224 | 1.5604 | 1.79427 | 0.0528298 |
| 2d | 2 | joint_gaussian_moments | 1.7349 | 1.04919 | 2.05578 | 0.218421 |
| 2d | 3 | mean_background_full | 1.68182 | 1.25009 | 1.34782 | 0.128838 |
| 2d | 3 | evolving_mixture_full | 1.78344 | 1.36074 | 1.86465 | 0.117486 |
| 2d | 3 | joint_background_evolving_mixture | 0.313781 | 0.50349 | 1.4267 | 0.657707 |
| 2d | 3 | joint_gaussian_moments | 1.8054 | 1.26217 | 1.09285 | 0.218011 |

2D H3 joint empirical 全窗口约 0.31378，前半约 0.5035、后半约 1.4267；后半/全窗口约 4.55 倍（0.6577 dex）。1D H3 从约 3.3467 到约 2.45/2.09。Gaussian 背景加 evolving response 的 1D H3 也变化约 0.338 dex（2.18 倍），所以不能把其他全协方差模型统称稳健。完整 observed/NULL 64 行见 statecheck/sensitivity.csv。

## 已完成的 latent-state 峰值诊断

| 维度 | H | 峰处 state ESS | 最大状态权重 | 最小 component χ² | 向量维数 |
| --- | --- | --- | --- | --- | --- |
| 1d | 2 | 1.702353 | 0.7453375 | 87.52313 | 46 |
| 1d | 3 | 1.000554 | 0.999723 | 686.6479 | 69 |
| 2d | 2 | 1.325108 | 0.8616747 | 124.6847 | 40 |
| 2d | 3 | 1.007831 | 0.9960999 | 1000.873 | 60 |

这是保存的单个 coupling/phase 峰处诊断，**不是 posterior-wide ESS**。H3 的约 1.00055/1.00783 ESS 与 0.999723/0.996100 最大权重表明几乎单状态支配；H2 也集中，但 component χ² 约 87.5/46、124.7/40，仍不应称好拟合；H3 约 686.65/69、1000.87/60 更是强烈失配警报。component 最小 χ² 不是 mixture goodness-of-fit，不能直接转作正式 p 值；这些读数既未证明科学模型成立，也不构成覆盖率认证。

## 原生时间、角度分支与数值离散限制

输入 141 帧，跨度 1245.054656 h，原生 cadence 8.893248 h，
124 个重叠四日窗口。实际相对时间 [0,24,120,144] h，曝光 6 h。
曝光短于原生采样；所谓积分仅对分段线性插值器精确，**曝光未被原生时间分辨率充分解析**。
冻结/演化响应是时间假设包络，不是时间过程已验证。
背景原生分辨率 400、响应 100，末端网格同为 64 并不消除原生分辨率差异；本结果是该不匹配下的 benchmark。
rmax=1、强度加权、10° 窗口等冻结设置见输入快照，不代表应重调这些参数。
局部最短角差提升不是跨 96 h 间隙的物理绕转测定。
输入背景非零 edge winding 为 1D 219、2D 179；背景及 observed closure winding 都为 0，不能把 edge winding 说成闭合失败。

## 网格、先验边界与 posterior 自检

α=.3 coarse→dense 的全 84 项最大 |Δlog10(c95)|=0.012062 dex，
约 2.816%，出现在 1d H3
joint_background_evolving_mixture observed。这是 coupling 与 phase 一起加密的合成变化，不能虚构单轴误差分解。
按模型分解（包含 observed/NULL 与两个维度/H）：

| 模型 | 最大 \|dense/coarse\| dex | 对应变化幅度 |
| --- | --- | --- |
| mean+iso | 0.00398738 | 0.922% |
| mean+diag | 0.0028348 | 0.655% |
| mean+full | 0.0038877 | 0.899% |
| mean+full+offset | 0.00390247 | 0.903% |
| frozen mixture | 0.00493634 | 1.143% |
| evolving mixture | 0.00510314 | 1.182% |
| joint empirical B,T | 0.012062 | 2.816% |
| joint Gaussian | 0.00420127 | 0.972% |
| indep Gaussian | 0.00410861 | 0.951% |

边界表覆盖每 α 的全部 84 行；lower/upper 是首/末一个**网格区间**的 posterior 概率，不是离散点质量。
实际全表 upper mass 仅为 0 到约 1e-15 的浮点尾差，没有上先验截断质量堆积证据，不能将 low-α 直接称为上界先验主导。
likelihood 的边界峰均应结合 posterior-checks.csv 中 lower/upper 方向解读；下界处的似然峰也不自动等于 posterior 全由先验主导。
首区间质量最大随 α 从约 1.36% 到 5.54%，说明不能忽略低耦合先验区间的作用。
未进行新先验实验，不声称先验下界无影响或先验独立性。逐项下/上质量、峰在 lower/interior/upper 均保留在 CSV。

| α | lower bin mass 范围 | upper bin mass 范围 | likelihood 峰在边界（/84） |
| --- | --- | --- | --- |
| 0.1 | 3.3e-47–0.0136 | 0–8.88e-16 | 40 |
| 0.15 | 1.89e-65–0.0171 | 0–1.11e-15 | 41 |
| 0.2 | 1.76e-72–0.021 | 0–9.99e-16 | 40 |
| 0.25 | 1.72e-70–0.0253 | 0–7.77e-16 | 39 |
| 0.3 | 6.37e-43–0.0304 | 0–6.66e-16 | 39 |
| 0.35 | 2.77e-33–0.0369 | 0–7.77e-16 | 41 |
| 0.4 | 2.03e-25–0.0468 | 0–7.77e-16 | 43 |
| 0.44 | 4.03e-23–0.0554 | 0–6.66e-16 | 44 |

逐 case NPZ 检查总数 840（allalpha 672 + coarse α=.3 84 + dense α=.3 84，按源计数，并非主表额外行）。
allalpha 实际另有 16 个 recovery-reference NPZ（8 α × 2 维度），连同其他 declared outputs 均经 SHA 校验；不混入 672 行主表。
以上逐 case NPZ 已检查：coupling/density/CDF/loglike 一维轴一致；phase_loglike 轴正确；
density 非负且按 c 梯形积分归一，CDF 单调且端点 0/1，与标准 log-uniform posterior 计算一致；
NPZ 与 CSV c95 相符，未对错误轴求和。边界质量与 likelihood peak 标记逐项核对。
这只是数值产品内部一致性检查，不是模型真伪或覆盖率证明。

## 图与来源复现

以下图位于 `/home/canoziia/GRMHD_DATA/analysis/eht-joint-days-results-v1`：

1. ![全部 α：实线 observed / 虚线 NULL](/home/canoziia/GRMHD_DATA/analysis/eht-joint-days-results-v1/c95-vs-alpha.png)
2. ![α=.3 全因子](/home/canoziia/GRMHD_DATA/analysis/eht-joint-days-results-v1/alpha030-factorial.png)
3. ![固定 H2 背景协方差几何](/home/canoziia/GRMHD_DATA/analysis/eht-joint-days-results-v1/background-geometry.png)
4. ![H3/H2 比值，>1 为变弱](/home/canoziia/GRMHD_DATA/analysis/eht-joint-days-results-v1/h3-h2-ratios.png)

图 3 在固定 H2、mean response 下对比 scalar/diag/full c95，以及平均/最小方差主轴；
它是方差椭球摘要，不将 trace 误称为椭圆面积。图 4 包含全部九个共享模型、全部八个 α。

消费 roots 为 `eht-joint-days-{inputs,alpha030,allalpha,alpha030-grid}-v1`。
每个源用 verify_stage(recursive=False) 校验所有已声明数值输出的 SHA，再检查消费文件为 declared outputs。
只消费处理后文件，不递归重哈希约 97 GB 原始父链。
源 manifest 字节原样复制为 source-receipt-*.json（非 canonical manifest/stage 名），
记录原始路径与 SHA；新汇总阶段仍可递归核验这些收据和处理后产品，而不引入 raw ancestors。
主表仅采用 allalpha 的 672 行，不重复堆叠 α=.3 的 84 行；单独确认 coarse α=.3 重复产物一致。
报告、CSV、图与 provenance 都由同一批已验证处理后输入生成；此报告副本与输出阶段内文件字节相同。

statecheck：已接入 /home/canoziia/GRMHD_DATA/analysis/eht-joint-days-statecheck-v1 的已完成输出（见 statecheck/ 与 provenance）。可用 `--statecheck ROOT` 显式接入完成产物；默认不等待、不触发额外实验。
其内容不被自动提升为已通过的科学结论。
