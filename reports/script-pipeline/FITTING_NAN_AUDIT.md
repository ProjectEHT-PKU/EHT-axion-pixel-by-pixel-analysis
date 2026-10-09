# 拟合算法与NaN处理专项核查

基线`b57bb64`。本轮仅核查，不修改执行代码、配置或已有科学结果。

## 1. 实际使用的算法

### 16相位响应：第一谐波DFT

`eht_axion/pipeline/raptor.py:346–376`，`fit_phase_dft`：

1. 若输入17态，先丢弃无轴子的state0，保留16个ψk=2πk/16。
2. EVPA周期π；使用mean(exp(2iχ))确定局部branch，展开到该branch。
3. 取baseline平均和第一谐波系数：

```text
z = (2/16) Σ (χk − meanχ) exp(−iψk)
A = |z|, phi = arg(z)
χmodel = baseline + A cos(ψ+phi)
```

对于完整均匀相位且branch正确的数据，它等价于常数+cos+sin的等权线性最小二乘，不是无初值控制的非线性curve_fit。数学归一化、相位符号正确。保存的是第一谐波投影，**没有保存拟合残差、高谐波或拟合质量门槛**。

EVPA跨度大、出现绕圈或低偏振相消时，局部branch可选错；非零resultant不是完整有效性证明。自己的纯cos A1.4rad反例被拟合为A≈.690rad，残差≈.664rad；这不是当前参考注入数据的实际幅度。

### 时间分布：Gaussian MLE矩估计

`raptor.py:466–529`，`summarize_fits`：

- 幅度用Welford累积finite样本的平均和population σ（ddof0）。
- 相位先unit-circle中心、2π局部展开，再同样估平均和σ。
- 这是高斯假设下的MLE，**不是对直方图用curve_fit，也没有证明分布高斯**。
- σ是时间状态散布，不是拟合参数标准误，不除sqrt141。
- 没保存Cov(A,phi)、联合复系数分布、时间相关或正态性质量。

### 背景

`raptor.py:534–665`默认method=mle，池化有限的无轴子EVPA(t+lag)−EVPA(t)差，估自由平均和population σ。histogram只用于诊断，不决定默认σ。可选legacy_hist才用curve_fit并检查拟合误差门槛。

### EHT观测统计

`pipeline/eht.py:119–194`不是正弦拟合：确定periodπ局部branch后估ensemble矩和covariance。当前要求同一特征全部候选样本完整，机器精度退化resultant设NaN；不完整特征不进入后续观测选点。

## 2. NaN逐阶段语义

| 阶段 | 当前处理 | 判断 |
|---|---|---|
| Q/U到EVPA | finite且不同时为零才有角；否则NaN | 正确。不能用Inf的atan2有限结果冒充有效。 |
| 1d积分 | 全四通道NaN padding当零flux；部分未知I/Q/U拒绝；IQU乘积float64 | 正确区分空间边界与未知观测，不把任意NaN补0。 |
| 单快照16态DFT | 任一相位NaN/Inf整点响应无效 | 保守且数学一致；不暗用nanmean破坏完整正交DFT。若想缺态拟合要另做带mask最小二乘。 |
| A数值为0 | A=0，phi=NaN，baseline可有效 | phi无定义正确；不应把这个phiNaN当作模拟失败。 |
| 时间summary | A和phi分别按finite样本累计；不足2则NaN | 有实质样本语义风险，见下一节。 |
| load_templates | 四个mean/std跨全部alpha有限才有效；不看保存counts | 对partial数据保护不足；一个alpha无效会排除所有alpha中的该位置。 |
| inference | 先按selected_indices取模板，再要求所有选择值finite | 正确，未选择的NaN不会污染推断；不是自动把NaN变成零信号。 |
| background | 仅两时刻都有有效角的差进入池化；零差保留 | 数学正确，可能发生依赖亮度/状态的missingness选择，需报告覆盖。 |

## 3. 自己复现的样本边界问题

用全新临时目录、人工141时间片，每位置只需minimum_samples=2：

- **只有2/141个快照有效**，仍得到有限mean/std，load_templates标为valid。
- **139个零幅度+2个非零响应**：A_count141，phi_count2，meanA≈.002837但meanphi只来自两次非零；仍标valid。
- **幅度与相位出现在完全不同两组快照**：分别count2，仍组成一个valid幅相模板（人工/导入异常可触发，正常DFT在A>0时两者同时有效）。
- **全部真实零响应**：A=0且phi无定义，phase summary NaN；最后整位置被排除，而不是表示有效零响应。
- 导入negative幅度可进入summary及load_templates；标准DFT abs系数不会产生负幅，但兼容入口缺语义拒绝。

前两项和全零相位属于内部可发生的边界；不应简单“统一把NaN改0”或“相位NaN一律失败”。需要有效响应/零响应/缺失响应分开表示。

## 4. 当前已有真实中间产物的只读核对

检查路径：`../eht-pipeline-runs/raptor/response/alpha-*/response-{1d,2d}.npz`及templates。

这些是用户已生成、修复前缓存。这里只读数值内容，不将其声明为通过当前源码cache/provenance验证；没运行推断或重新预处理。

- 八个alpha，1d各360方向、2d各4096pixels。
- 每位置A_count和phi_count均为**141**。
- 响应没有zeroA，未见A有限而非零、phi缺失。
- **上节缺样本/零幅度边界目前没有在这批已存产物中触发。**

### 单谐波质量抽查

每alpha只读processed的第一、中间、最后三个快照（共24个），未加载完整大cube：

- 所有4096pixels均有有效拟合；参考注入A最大约.116554rad。
- 16态EVPA resultant最小约.98628，未见严重branch退化。
- RMS残差/A中位约**0.012%–0.119%**；99分位约**0.20%–2.91%**。
- 个别弱响应像素的relative residual可达约35%；不能由全图中位数证明每个像素好。
- 这只抽查2d，未完成141快照/全部1d拟合质量检查，也不是最终EHT mask选点统计。

**判断：当前参考注入下第一谐波算法总体有支持，不是首先应该换成非线性curve_fit。更应补残差/高谐波诊断。**

## 5. 比NaN更值得重视的统计近似

当前nominal响应相当于meanA*exp(i meanphi)。直接平均复响应则是mean(A exp(i phi))，一般不同。

只读完整响应得：

- 1d两种复均值的全方向median相对差约.7%–2.2%，max约2.0%–4.5%。
- 2d全图部分位置差异大得多，含低强度/弱响应区域；不能当最终有效mask误差。
- 按目前offset4选36方向，abs corr(A,phi)的median约**.14–.41**，最大约**.54–.84**。
- 2d全图phase σ最大约1.63–1.86rad；屏幕2≤r≤9Rg的几何圈中也有少量σ>.5rad。这个圈**不是最终EHT选择mask**。

当前`model.py:163–173`只传播两项：

```text
K²[ sin²θ σA² + meanA² cos²θ σphi² ]
```

一般一阶传播还有：

```text
2 K² meanA sinθ cosθ Cov(A,phi)
```

真实响应中的相关性说明不能默认忽略。共同nuisance phase平均可能部分减弱交叉项，但参数logdet、posterior和空间相关不是线性平均，不能据此认为影响必定为0。

同样，宽相位的非线性平均还有E[sinphi]与sin(Ephi)的差，不只补交叉项即可保证完整正确。

## 6. 建议，不在本轮实施

1. 将每快照响应保存为cos/sin系数或复系数z=Aexp(iφ)；A=0对应z=0，无需定义任意phi。
2. 显式区分valid-zero、valid-nonzero、missing，并保存原因、joint_count、总时间片数及coverage。默认生产全141是否严格要求完整，应由分析配置明确，而非以最低2样本代替。
3. 对完整16态保存residual RMS、harmonic power、branch/resultant及弱响应标志；缺态拟合如需要用显式mask线性最小二乘，不改DFT为nanmean。
4. 用复系数的联合均值/covariance或snapshot直接传播，对照现幅相delta近似。仍需分别考虑空间/跨daypair/time相关。
5. 对真实最终mask统计质量和比较六曲线，量化影响后再决定科学估计器。不要把本轮百分比诊断直接当c95偏差。
6. 背景保存每pixel/时间有效计数、missing原因，检验池化/曝光/lag对应，而不是把缺失或真实零差当同一种NaN。

## 结论

DFT和MLE矩估计的代数正确；不存在“所有NaN都被填0进入拟合”的行为。但当前时间summary的partial/zero-response语义和counts gating不够完整，算法近似也缺联合响应及质量诊断。当前已有产物没有缺样本问题；实际A–phi相关性则已发现，值得优先用复响应统计作方法对照。本轮不自动改变算法或选点。
