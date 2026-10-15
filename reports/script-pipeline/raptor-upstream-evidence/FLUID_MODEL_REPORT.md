# RAPTOR 流体、模型及运行参数独立只读审查

## 0. 结论先行：重跑前需要解决什么

**本审查不认为 `SIGMA_CUT=1` 是无条件、唯一正确的物理值。**它是当前主运行的已复核二进制阈值，也是 EHT V/VIII 和本项目引用的 Mizuno 等人模型采用的、排除高磁化不可靠热力学区域的保守基准。`1000` 则保留了大量 `1<σ≤1000` 的流体，包括可能由数值 density/pressure floors 人工填充的漏斗；可用于特定诊断、模型对照或宽松 cut 的压力测试，但不能仅凭“保留更多流体”就称为更真实的 fiducial。

建议按以下优先级处理：

1. **Must-fix／重跑阻断：把运行配置及数据格式锁定为可追溯 manifest。**本地 core 模板为 `SIGMA_CUT=1000`，本地主 run 和云端两处为 `1`。用户提供的父审查已复核主 run binary（SHA 前缀 `7137ba…`）实际比较常量为 `1.0`，`eq_fakevalue` binary（`bd738b…`）为 `1000.0`。不要用 core 模板覆盖主 run 后仍把结果当同一实验。重新构建的实际 include 路径、header 内容、最终阈值和运行目录必须一起记录；本审查不重复二进制反汇编。
2. **Must-fix／模型说明：0.6 Jy 校准不能靠文章文字或 `M_UNIT=3e25` 直接认定完成。**`Run.py` 固定质量归一化，没有求解 0.6 Jy 的过程；当前 `R_LOW=R_HIGH=1`、cut、温度配方及外边界组合的实际 flux 校准仍未数据验证。更改 cut／电子温度配方后原则上需要重新校准该模型的时间平均 flux，而不是只把输出图乘常数。
3. **Must-fix／输入适配：不要直接用仓库自带的 Sgr A* `model.in` 启动项目版 reader。**该模板仍是旧字段布局，缺三个轴子参数；项目版按位置读入三个额外字段，之后才读 `MAX_LEVEL`。主工作流实际用 `Run.py` 临时生成 M87* 输入，本地运行结束后删除该文件。这两类输入不能混为一谈。
4. **Must-fix／如果更换 GRMHD 数据：按 footer 和模拟参数重新确认 `NSPIN`、SFC、网格和坐标。**当前 mounted footer 已独立确认 `neqpar[3]=0.9375`、`neqpar[5]=20`；因此当前数据应保持 `NSPIN=3`。不能因注释含 “Mizuno” 就改为 `5`，否则本数据读到的是 `a=20`。`hslope`、根网格尺寸和 SFC 不是可任意调整的辐射物理参数。
5. **必须显式决定而非偷偷同步：thermal Faraday fit 和电子温度闭合。**当前 thermal `rho_V` 是旧 Dexter 型 fit；最新 main 默认改用 Shcherbakov 型 fit。2023 年 thermal 符号改动随后被上游撤回，不能只摘中间提交“修正”主 thermal 链。另，当前 `theta_e` 闭合并不与 EHT V/VIII 的电子／离子分配公式相同，也不严格等于 Mizuno 2021 文中 R–β 的近似实现。这是继承的模型选择，不是 fork 偷改，但应在论文中准确陈述。
6. **Fiducial 建议：**若目的是延续当前主 run、并以“可信热力学区域的 thermal MAD 全等离子体模型”作为基准，优先保持 cut=1，保留吸收和 Faraday 系数；在 flux、温度配方和解析输入确认后，再以 cut=10 做文献支持的首要 systematic，并把 cut=1000 作为高磁化污染／历史诊断压力测试。若研究目标明确是漏斗辐射或对产生等其他物理，则 cut=1 可能排除真正有辐射的物质，需要额外的电子／正电子与加热模型；此时不能只把 cut 调到1000就宣称解决了物理问题。

本报告只审查模型／参数及与 cut 有关的最小调用链，不重新审查轴子算法或光线积分器。下面给出证据和不确定性。

---

## 1. 范围、基线、证据等级与实施限制

### 1.1 四处源码范围

下文使用以下缩写：

- **LC（local core）**：`/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel/raptor`；模型在 `model/bhac/`。
- **LR（local run）**：同一挂载目录下的 `raptor_run`；其根目录有独立的 `definitions.h`、`model.c`、`model_definitions.h`、`grid_mks.in`。
- **CC（cloud core）**：`/tmp/eht-github-raptor-GckNId/raptor`；模型在 `model/bhac/`。
- **CR（cloud run）**：`/tmp/eht-github-raptor-GckNId/run/raptor_run`。

四处不意味着四份独立物理实现：本次逐文件 SHA-256 比较确认四处 `model.c`、`model_definitions.h`、`grid_mks.in`、`model_functions.h`、`model_global_vars.h` 分别完全相同；`definitions.h` 唯一内容差异是 LC cut=1000，而 LR/CC/CR cut=1。两 core 的 `src/pol_emission.c` 也完全相同。

### 1.2 Git 归属方法

完整上游仓库：`/tmp/eht-raptor-upstream-audit-lVJjW7/upstream`。

- 上游 main：`b91e937`。
- 项目 fork project-main：`b6df33c`。
- 实测 `git merge-base b91e937 b6df33c`：`cdc2547b868ed8d6ae0c9a5428bc1ea86fce901a`。

先比较 `git diff <merge-base> b6df33c` 识别 fork 的实际贡献，再独立比较 `<merge-base> b91e937` 识别上游后来改变。不能拿 fork 与今天 main 的全部差异都归罪为“自定义问题”。本次保留了这两份 diff。

### 1.3 哪些是亲自检查，哪些是沿用已复核证据

- **亲自检查**：四处小型源码／header／网格文件、运行脚本、Git 原始提交与 diff；文章两份指定 TeX 的模型／图4相关段落；EHT V/VIII 缓存论文；另以有超时的网络请求取得并检查 Mizuno 2021 与 RAPTOR II 一手论文；仅做常数公式的代数数值核算。
- **引用独立审查结果，不冒称本次新实验**：mounted footer/尺度结果见 `/tmp/eht-mounted-constants-qu_7_khw/comparison.json`；主 run 与 eq_fakevalue ELF 阈值采用用户转交的父审查复核结果。
- **没有做**：修改源代码／挂载数据／远端；执行旧 RAPTOR binary；构建或运行辐射转移；读取 GB 级 GRMHD 主体；重新计算图像／轴子科学数据／flux／sigma 直方图；重复父 agent 的大反汇编或另一 agent 的积分器审查。
- 因而“源码参数是什么”证据很强；“该参数对当前数据的数值影响多大”除二进制阈值外，本次不作伪造的定量结论。

---

## 2. 实际 fork 变化与上游后续变化

### 2.1 可归属于 fork 的模型变化

`diff merge-base→b6df33c` 显示：

1. `model/bhac/model.c`：只有新增的注释调试输出；**没有改变** `get_fluid_params` 的 sigma cut、温度、单位、primitive 转换或读取物理公式。文件 mode 变化不是物理变化。
2. `model_definitions.h`：`SFC 0→1`；`NSPIN` 的 active 值仍为 `3`，只是新增了“某种 Mizuno 数据要5，本文件夹数据要3”的注释。
3. `grid_mks.in`：根网格 `128×64×64→96×48×48`；log-radius 下界 `0.17→0.17395330712343798`；上界 `8.1117280833→7.824046010856292`；`hslope 0.25→0`。这些是数据格式／坐标适配，而非可以随意“恢复上游”改善物理的参数。
4. `definitions.h`：`USERSPEC 0→1`；Camera 中移除 `pdf/avg`；fork commit 本身 **没有把 SIGMA_CUT 从1改成1000**。
5. `src/metric.c`：本次 diff 范围只有空白，没发现 active metric 公式的 fork 改写。
6. `src/core.c`：增加三个轴子输入字段、对应日志及调用参数；这导致旧 `model.in` 模板布局不再匹配。
7. `src/pol_emission.c`、`src/constants.c`：fork 相对 merge-base 没有变化。

因此 LC 的 cut=1000 是当前 mounted core 相对云端／fork快照的额外状态，不能归为 `b6df33c` 的共性设置；它的实验用途应结合历史／专用 diagnostic run 说明，不能用它覆盖主运行身份。

### 2.2 上游后来改变，不是 fork 自定义

模型部分，`b91e937` 在 `a=neqpar[NSPIN]` 后增加 `fabs(a)>1` 报错并退出。它能抓住当前数据若误用 `NSPIN=5` 导致的 `a=20`，但不改变正确读取时的物理结果。也不检查 `NSPIN` 是否越界，不检查非有限值，更不能识别“读错但碰巧落在[-1,1]内”的字段。

上游的 `grid_mks.in` 与 `model_definitions.h` 在 main 上仍保留原示例网格、`hslope=.25`、`SFC=0`。这不构成当前数据应该改回这些值的证据。

---

## 3. 四处参数对照表

### 3.1 编译参数与流体模型

“上游”列写 merge-base／main；若只写一个值，则二者相同。行号以 CC 的当前文件为准；四处相同文件可直接对应。

| 项目／证据位置 | LC | LR | CC | CR | 上游 merge-base／main | 文章／一手依据及结论 |
|---|---:|---:|---:|---:|---|---|
| `SIGMA_CUT`，`definitions.h:111` | **1000** | **1** | **1** | **1** | 1／1 | EHT V fn115、§7低密度处理；EHT VIII §4.2；Mizuno §2.3均支持1作为可信热力学区域基准，但不是宇宙中σ>1无辐射的定律。主binary=1已由父复核。 |
| `THETAE_MAX/MIN`，`:112–113` | 100／1e-3 | 同 | 同 | 同 | 同 | 无本项目专门校准证据；是整块 plasma 排除窗，不是温度 clamp。需检查边界是否切到主要发射区。 |
| `DF`，`:79` | TH | TH | TH | TH | TH | main.tex:418热各向同性电子；Mizuno §2.3热 Maxwell–Jüttner。KAPPA/POWER不活跃。 |
| `EMISUSER`，`:45` | 0 | 0 | 0 | 0 | 0 | polarized主链直接 `evaluate_coeffs_single`；不能因文件中有混合模型函数就说用了非热电子。 |
| `sigma_min`，`model.c:1174` | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | **不是SIGMA_CUT**；供混合 thermal/kappa 权重，当前主polarized调用不使用该混合权重。 |
| `POL/RAD_TRANS/LIGHT_TRANSPORT` | 1／1／1 | 同 | 同 | 同 | 同 | 主模型确实开启偏振转移；关闭吸收／Faraday的诊断代码是注释，不是主链状态。 |
| `metric/logscale`，`definitions.h:102–107` | MKSBHAC／1 | 同 | 同 | 同 | 同 | `r=exp(x1)`，BHAC θ映射；不能混用MKSHARM的hslope定义。 |
| `SFC`，`model_definitions.h:11` | 1 | 1 | 1 | 1 | **0／0** | fork数据排序适配；必须与模拟输出根block顺序匹配，不能仅以leaf总数一致证明正确。 |
| `NSPIN`，`:33` | 3 | 3 | 3 | 3 | 3／3 | mounted footer独立确认index3=.9375、index5=20；当前3合理，注释5只适用另一格式。 |
| `BPOL`，`:37` | PLUS | PLUS | PLUS | PLUS | PLUS | 极性翻转会改变Faraday/V等有符号观测量，不只是B²；无证据要求翻转。 |
| `NPRIM/NP` | 8／8 | 同 | 同 | 同 | 同 | primitive存8量；文件`nwini`并非固定8，LFAC/XI取末两列，格式须确认。 |
| `rcam`，`definitions.h:122` | 1e4 Rg | 同 | 同 | 同 | 同 | 远场相机，**不是16.8Mpc本身**，也不是RT外边界。适合作为基准，有限距离与相位收敛需另测。 |
| `RT_OUTER_CUTOFF`，`:129` | 50 Rg | 同 | 同 | 同 | 同 | 等离子体转移范围，不是流体域外缘；EHT VIII通常取50–100Rg且讨论外层Faraday screen。应做50→100等systematic。 |
| `horizon_marg`，`:135` | .01 | 同 | 同 | 同 | 同 | 实际`CUTOFF_INNER=r+×1.01`，不是绝对+.01Rg。数值边界，需收敛，不是额外辐射物理。 |
| `num_frequencies/FREQS` | 1／FREQLOG | 同 | 同 | 同 | 同 | 单频；`FREQS_PER_DEC=1`不表示已模拟宽带。 |
| `USERSPEC`，`:55` | 1 | 1 | 1 | 1 | **0／0** | user spectrum函数体被整体注释，当前没有额外谱的物理贡献；不等于EMISUSER开启。 |
| `DEXTER` 宏 | 无 | 无 | 无 | 无 | 无／**0** | 旧源码无开关而直接用Dexter型thermal rhoV；main默认换Shcherbakov fit。详见§7。 |
| spin越界警告 | 无 | 无 | 无 | 无 | 无／有 | 上游b91新增健壮性保护，当前正确NSPIN时不改变数值。 |

### 3.2 运行时网格文件

`grid_mks.in` 是程序从**当前工作目录**以 `fopen("grid_mks.in","r")` 读取的文件，不是编进 binary 后永远不变的常量。

| 项目 | LC | LR | CC | CR | 上游 merge-base／main | 结论／待证据 |
|---|---|---|---|---|---|---|
| `nxlone[1:3]` | 96,48,48 | 同 | 同 | 同 | 128,64,64 | 是根level网格，并非整个AMR数据总分辨率；必须匹配生产模拟。 |
| `xprobmin[1:3]` | .17395330712343798,0,0 | 同 | 同 | 同 | .17,0,0 | 当前`r_min≈1.19Rg`，在a=.9375地平线以内；与地平线外停止光线并不冲突。 |
| `xprobmax[1:3]` | 7.824046010856292,.5,1 | 同 | 同 | 同 | 8.1117280833,.5,1 | 当前`r_max≈2500Rg`；θ/φ输入随后乘2π成为[0,π]/[0,2π]。 |
| `hslope` | **0** | **0** | **0** | **0** | **.25／.25** | BHAC映射`θ=x2+.5h sin(2x2)`，当前为普通角坐标。不是HARM映射中h=0的含义。 |

footer证明了spin和block size；**单靠footer不能证明根网格、hslope和SFC顺序全部正确**。还需要生产模拟的 parameter 文件／坐标说明／既有可信校验。切勿从“上游示例=.25”推出“当前0是bug”。

### 3.3 运行输入：必须区分模板与生成器

LC/CC 的值取 `raptor/run/model.in`；LR/CR 的有效配置取 `raptor_run/Run.py:64–81`（主工作流生成器）。CR还保存了一份旧Sgr A* `model.in`，本地LR没有该文件；其旧内容与core模板相同，**不是M87主生成器参数**。

| 运行项 | LC旧模板 | LR主Run.py | CC旧模板 | CR主Run.py | 上游旧模板（merge-base／main） | 依据／判定 |
|---|---:|---:|---:|---:|---:|---|
| `MBH` [Msun] | 4.14e6 | **6.2e9** | 4.14e6 | **6.2e9** | 4.14e6 | 主run的M87尺度合理；EHT V/VIII用6.2e9，Mizuno引用文用6.5e9，不能悄悄替换后只改图标签。 |
| `DISTANCE` [kpc] | 8.127 | **16800** | 8.127 | **16800** | 8.127 | =16.8Mpc；EHT V/VIII例16.9，Mizuno16.8；都是带约定的尺度，当前一致性应优先。 |
| `M_UNIT` [g] | 1e18 | **3e25** | 1e18 | **3e25** | 1e18 | 是流体密度归一化，不是BH质量；文章main:355声称0.6Jy，尚未实证该数满足目标。 |
| `R_LOW` | 1 | **1** | 1 | **1** | 1 | 有文献基准支持，非唯一电子加热模型。 |
| `R_HIGH` | 1 | **1** | 1 | **1** | 1 | Rlow=Rhigh时trat恒1，β依赖退化为常温度比；不是β依赖效果已经验证。 |
| `INCLINATION` [deg] | 50 | **163** | 50 | **163** | 50 | main:367、old:162与EHT正spin视向支持163°约定。 |
| `IMG_WIDTH/HEIGHT` | 200/200 | **100/100** | 200/200 | **100/100** | 200/200 | 400高分辨率分支定义存在但`RUN_HIGH_RESOLUTION=False`，不能说主Run.py默认400。 |
| `CAM_SIZE_X/Y` [Rg] | 40/40 | **46.2261/46.2261** | 40/40 | **46.2261/46.2261** | 40/40 | 按当前旧常数≈168.414μas全宽；这是FOV，不是rcam。 |
| `FREQ_MIN` [Hz] | 230e9 | **230e9** | 230e9 | **230e9** | 230e9 | 对应EHT1.3mm和文章230GHz目标。 |
| `FREQS_PER_DEC` | 1 | 1 | 1 | 1 | 1 | num_frequencies=1下只有f=0，实际就是230GHz。 |
| `STEPSIZE` | .01 | **.005** | .01 | **.005** | .01 | 当前run更小步长；不是准确性已证明，本审查不重评积分算法。 |
| `MAX_LEVEL` | 1 | 1 | 1 | 1 | 1 | 图像网格参数，不是GRMHD AMR层数。 |
| 三个轴子输入字段 | **缺失** | 生成 | **缺失** | 生成 | 原reader无需 | 项目reader新增字段后模板未更新；不得直接启动旧模板。 |
| 本地临时`model.in` | — | **不存在** | — | 存旧模板 | — | Run.py执行结束`os.remove("model.in")`，不存在不意味着运行没输入。 |

补充：LR实际循环`axion_omega_list=[0.1]`，CR为`[.3,.25,.2,.15,.1]`；参数字典里的`.4`是被循环覆盖的初始值，不是两处当前批处理都在跑.4。此处只记录运行覆盖，不重新审查轴子定义。

---

## 4. SIGMA_CUT 的精确定义、作用与物理论证

### 4.1 代码的σ是什么

`model.c:1063–1182` 先插值ρ、内部能和磁场／速度，再用metric重建四速度和共动磁四矢量 `B_u/B_d`。用于cut的量不是简单欧氏三磁场平方：

```c
Bsq = fabs(B_u[0]*B_d[0] + B_u[1]*B_d[1]
         + B_u[2]*B_d[2] + B_u[3]*B_d[3]) + 1e-6;
B = sqrt(Bsq) * B_unit;
```

理想物理含义是共动 `b²=b^μ b_μ`，代码额外取绝对值并加小量。BHAC采用吸收4π的无量纲磁场；其单位转换是 `B_unit=c sqrt(4πρ_unit)`，因此在忽略regularizer时：

\[
σ=b^2/ρ_{\rm code}=B_{\rm cgs}^2/(4πρ_{\rm cgs}c^2).
\]

**不是**磁压／气压（那是与β有关），也**不是**含焓的 `σ_h=b²/(ρh)`。Mizuno文中的重联加热部分另用σ_h，不能把两者的阈值直接等同。

实际代码另保存：

```c
modvar.sigma = Bsq / (rho + 1e-6);
modvar.sigma_min = 1.0;
// 但cut实际用的是：
Bsq / (rho + 1e-20) > SIGMA_CUT
```

因此低密度时 `modvar.sigma` **不等于cut真正比较的σ**。若未来记录diagnostic只输出前者，可能误以为某些“σ<1”的点被错误删除。应区分 `sigma_cut_ratio` 与混合模型所用的regularized sigma。

`SIGMA_CUT` 使用严格 `>`：恰等阈值时不由此条件删除；文献有时写`σ<1`或`σ≥1`，不要过度解释连续流体中理想边界点，但数值ceiling附近可能有聚集，需留意。

### 4.2 被cut后关闭什么，为什么轴子还能传播

`model.c:1193–1199` 的任一条件成立都会：

```c
modvar.n_e = 0;
return 0;
```

条件是磁化比值过大、`r>50`、`theta_e>100`或`theta_e<1e-3`。它**不删掉磁场几何／四速度，也没有把光线终止**。

项目 `src/pol_rte_integrator.c:973–978` active调用条件为：

```c
if (get_fluid_params(X_u, &modvar) || r_current < RT_OUTER_CUTOFF)
    pol_integration_step(...);
```

所以在`r<50`范围内，即使sigma cut返回0，仍会进入偏振积分。当前thermal链中：

- `jI/jQ/jV` 与`n_e`成正比；`jU=0`。
- 热吸收`aI/aQ/aV`通过对应j除以Planck函数，随n_e归零；`aU=0`。
- 热Faraday conversion／rotation `rhoQ/rhoV`含`n_e`，随之归零；`rhoU=0`。
- 因而在有效、有限的温度／场／频率条件下，**关闭的是全部等离子体发射、吸收和Faraday系数，不只是“spine emission”**。末尾注释“excludes all spine emission”是目的性简写，不是完整行为定义。
- 轴子端点增量不乘`n_e`，继续传入偏振更新。**已经存在的偏振光仍可以经过cut区域积累轴子旋转；cut区域本身不再产生新光。**该事实不意味着从零Stokes状态可以凭轴子生光。
- 这也是将caller从旧返回值门控改成OR半径门控的合理动机之一；本审查只确认这一接口语义，不重做轴子积分算法审查。

两个重要限制：

1. 温度为0、NaN或极端下的Planck／Bessel等计算未必保证`0×非有限=0`，不能把数学上n_e=0当作所有数值异常都被安全旁路的证明。
2. 更早的“域外／找不到grid／r<1”等返回0分支未统一清零或重置modvar；在OR调用条件下需要特别处理vacuum状态。这是接口健壮性问题，交由已完成的pol_rte审查结论整合，不能误认cut末尾这个分支也没有设置n_e。

### 4.3 为什么1有依据，为什么不能宣布唯一正确

**一手依据，不是源码默认值投票：**

- **EHT V**，2019，arXiv:1906.11242，fn115及§7低密度处理：喷流核心密度由数值floors设置，RT对`B²/ρ>1`的区域置零；floor触发区域的温度不再可靠。该文同时承认这种cut可能低估毫米辐射、从而高估吸积率和jet power，尤其高磁化MAD的低／高密度混合边界存在不确定性。
- **EHT VIII**，2021，arXiv:2105.01173，§4.2（PDF p9）：定义`σ=B²/(4πρc²)`，高σ区域内部能难以准确演化，主库cut=1；另做cut=10全模型补充库，其所用model-scoring指标总体分布没有改变。注意：这不是“每像素Q/U、轴子振幅／相位及时间协方差不变”的证明。
- **本项目直接引用的 Mizuno et al. 2021**，arXiv:2106.09272，§2.3（PDF p4）和相关图4讨论（PDF p6）：明确用`σ=b²/ρ`，只考虑`σ<1`发射，因为σ>1的密度／温度可能被simulation floors影响。§2.1（p3）另给模拟density/pressure floors和`σ_max=100` ceiling。这是比“上游header=1”更贴近本项目数据来源的理由。

因此cut=1是**数值可信性选择**，不是等离子体在σ=1处存在真实发射关断。真实漏斗可能有电子／正电子对、磁重联加热或非热粒子，可能辐射；单流体GRMHD加简单thermal闭合未必可信地描述它们。

### 4.4 对1000的合理定位

- cut=1000保留`1<σ≤1000`，更接近“不做磁化排除”的宽松对照，不等于无cut：仍有1000以上、温度窗、r=50等限制。
- 若模拟确实采用如Mizuno文中`σ_max=100`的ceiling，1000超过该ceiling，通常会放行几乎所有因高磁化被数值处理的区域；但当前snapshot的真实σ分布不能从文章或footer的一个100字段直接证明，插值／重建及B²小量也会改变比值。
- **合理用途**：历史eq_fakevalue／人为发射对照；有明确额外热力学模型时探索高σ贡献；测试结论是否受cut强烈控制；作为“包含不可信漏斗”的极端系统误差压力测试。
- **不合理推断**：“σ>1被删除，所以一定应放到1000”；“所有代码默认1，所以1000绝对错误”；“1000结果更亮所以更符合观测”；“调MUNIT把总flux归一后每像素轴子结果必然一致”。这些都没有足够物理／数值证据。
- 1000不是EHT VIII已测试的10；不能把“cut10总体评分稳定”的结论外推三个数量级。

### 4.5 推荐实验矩阵（仅计划，本次未运行）

分成两类比较，避免混淆：

1. **固定MUNIT诊断**：同snapshot、同电子配方、同camera和频率，cut=1/10/1000；看高σ新增发射、吸收、Faraday区、总flux和像素级变化。它回答“同一density模型放行高σ会怎样”。
2. **各模型独立flux校准后的科学systematic**：对各cut模型用同一稳态时间窗重新确定一个时间恒定MUNIT，使平均230GHz compact flux符合预先声明目标（文章当前为0.6Jy）；保留快照间flux涨落，再比较偏振形态、RM／Faraday深度、轴子拟合量、Fig4统计和最终约束。它回答“同样满足flux要求的不同热力学可信域如何改变结论”。

至少cut1和10，若资源允许加1000作为压力测试。不要逐snapshot强行归一0.6Jy后还宣称保留了真实flux变率；不要因EHT评分稳定就省略像素级检验。若1000使图4或约束大变，先判断是不是人工floor污染，而不是自动把三者等权当物理先验。

---

## 5. 热电子闭合、R_LOW/HIGH、单位和floors

### 5.1 实际电子模型：Rβ框架，但当前退化成常比值

`model.c:1164–1182`：

\[
β=\frac{u(\gamma-1)}{(Bsq+10^{-6})/2},\quad β_{\rm trans}=1,
\]
\[
R(β)=\frac{R_{\rm HIGH}β^2+R_{\rm LOW}}{1+β^2},\qquad
Θ_e=\frac{u}{ρ}\frac{m_p/m_e}{3[1+R(β)]}.
\]

注意代码变量名`b2`在此只是`β²`，**不是磁四矢量模平方Bsq**。低β趋向R_LOW，高β趋向R_HIGH；这些参数是离子／电子温度比的模型参数，不是分布函数kappa值，也不是sigma cut。

当前两者都是1，故`R(β)=1`，`Θ_e≈306.0 u/ρ`，温度依赖内部能与密度，但不再通过R(β)依赖磁化／β。这是一种热电子基准，不是“完全没有加热假设”。它并没有从BHAC文件直接读取或演化electron entropy：reader只转成8个primitive量，再由u/ρ和R重建电子温度。

### 5.2 与EHT/Mizuno公式不能一概等同

当前闭合的`1/[3(1+R)]`从merge-base继承，fork未改，main也未改。它与以下一手配方不同：

- **EHT V Eq.(7)、EHT VIII Eq.(16)**：非相对论离子Γp=5/3、相对论电子Γe=4/3时，`Θ_e=(u/ρ)(mp/me)×2/[3(2+R)]`（忽略电子质量微小修正）。当前与此的比值为`(R+2)/[2(R+1)]`，R=1时当前是该EHT配方的**0.75倍**。
- **Mizuno 2021 §2.3 Eq.(8)及fn2**：Rβ后处理中用`p_g`近似代替`p_g-p_e`，等价于`Θ_e≈(p_g/ρ)(mp/me)/R`；Γg=4/3时为`(u/ρ)(mp/me)/(3R)`。当前对R=1为该近似的**0.5倍**。该文自身说明这一近似对应Ti≈Tg、Ti≫Te，所以R=1时的近似前提本来就需要小心。

当前闭合可理解为采用共同Γ=4/3并按两组分温度比分享总内部能的有效配方；不应未经说明就改成EHT公式。但若文章要声称严格复现EHT V/VIII或Mizuno的电子温度，就必须把这些差异写清并验证。**这是需要决定和记录的物理模型选择，不是本审查仅凭文献就下达的唯一代码修复。**改变闭合会同时影响发射频率位置、Faraday深度、温度cut通过率及flux校准。

### 5.3 R_LOW/HIGH 的fiducial与systematic

- EHT V固定Rlow=1并比较Rhigh=1,10,20,40,80,160。
- EHT VIII比较Rlow=1,10，Rhigh同一系列；不同评分方法对Rlow的偏好并不一致。冷电子会提升Faraday深度并改变去偏振／EVPA模式。高spin MAD并不意味着Rhigh=1自动与所有多波段／polarimetric约束相容。
- Mizuno 2021采用Rlow=1、Rhigh=1,5,10,20,40,80,160，比较两种电子加热模型；其结果支持探索一组有效配方，不证明某一参数对每像素时间变率完全等价。
- 本项目若首先复现历史基准，可保留1/1，但须称为选定benchmark；若最终要把约束称为稳健的M87模型推断，应至少加入较冷电子的合理Rhigh分支和/或Rlow=10，并对各模型重新做时间平均flux校准与偏振一致性检查。不能无数据指定“Rhigh必须80”等单一值。

### 5.4 Ne/B单位和BH质量尺度

`model.c:846–855`：

\[
L_u=GM_{BH}/c^2,\quad T_u=L_u/c,\quad
ρ_u=M_{UNIT}/L_u^3,\quad U_u=ρ_uc^2,
\]
\[
B_u=c\sqrt{4πρ_u},\qquad Ne_u=ρ_u/(m_p+m_e).
\]

这里`B_u`单位符号指B_unit而不是磁四矢量分量。对应全电离氢、一质子一电子的设定，不含pair multiplicity；`n_e=ρ_code Ne_unit+1e-6 cm^-3`。若研究高σ漏斗中的pair plasma，不能沿用此Ne关系又宣称已经模拟对产生。

仅用当前header旧物理常数、MBH=6.2e9、MUNIT=3e25代数核算：

- `L_unit≈9.1573806913e14 cm`；
- `T_unit≈8.4849261 h`；
- `RHO_unit≈3.90667648e-20 g/cm³`；
- `Ne_unit≈2.3343889e4 cm⁻³`（每单位code density）；
- `B_unit≈21.0053 G`（每单位共动code磁场）。

这些是**转换因子，不是当前图像的发射加权Ne或B测量**。footer的M=1是几何无量纲质量，不应把运行MBH改成1Msun，也不意味着模拟数据已经有cgs密度尺度。此前独立常数报告已经确认生产尺度6.2/16.8；本次只沿用和核算，不重推整条成像pipeline。

更改MBH或MUNIT不仅改变图像角尺度：`Ne_unit∝MUNIT/MBH³`，`B_unit∝sqrt(MUNIT/MBH³)`；局部u/ρ给出的Θe不直接依赖MUNIT，但synchrotron谱位置与Faraday会改变。整体Stokes输出乘一个0.6Jy比例因子不能替代重新求解这样的非线性转移。

### 5.5 floor/regularizer的区分

四处和上游均有以下继承项：

| 类别 | 代码行为 | 不应误称 |
|---|---|---|
| cgs电子数小量 | `n_e=ρ Ne_unit+1e-6` | 不是BHAC无量纲density floor；也不保证ρ正。 |
| 磁四矢量模小量 | `Bsq=abs(bμbμ)+1e-6` | 不只是NaN保护；给B和cut引入实际最小值。 |
| beta分母 | `(Bsq+1e-6)/2` | Bsq已含小量，实际又加一次。 |
| stored sigma分母 | `ρ+1e-6` | 与cut的`ρ+1e-20`不同。 |
| primitive内部能恢复 | 负UU先用DS熵恢复；仍负则`UU=1e-14` | 不是电子温度下限，也不是可靠物理加热。 |
| BHAC演化floors | 模拟时density/pressure重置、磁化ceiling | 与RAPTOR这些后处理regularizer、SIGMA_CUT、sigma_min均不同。 |

特别地，在上述MUNIT下Bsq小量对应`B≈0.0210G`的数值底值；当真实b²很小时，cut比值会近似`1e-6/ρ`。因此cut1可能把ρ≲1e-6的点排除，即便其真实磁化比未必很高；cut1000则把这一效应移到ρ≲1e-9。这个推导说明为什么未来应核查低密度/弱场区，不是本次已经证明主发射区受污染。不能仅因为小量数字小就认为物理上处处可忽略。

本reader没有统一的正ρ强制clamp，也没有全面finite/positive输入验证。温度比较遇NaN可能全部为false，`n_e=0`也不能保证后续运算永不NaN。对于重跑的健壮版本，应优先做好合法状态与vacuum处理；但本次未改代码，也未把继承问题归为fork自定义。

---

## 6. BHAC读取、坐标和边界是否合理

### 6.1 文件读取与NSPIN

`model.c:529–608` 从文件末尾40字节读nleafs、level、维数、nwini、nws、neqpar数量、iteration和时间，再回退读取各方向block大小与neqpar。`LFAC=nwini-2`、`XI=nwini-1`，随后依赖固定的D/S/B/DS索引做conserved→primitive：

- `ρ=D/LFAC`；
- `u=(XI/LFAC²-ρ)/γ`（原表达式约掉γ−1后的等价式）；
- 从spatial metric、conserved momentum和磁场重建速度；
- 负内部能用DS熵恢复，失败后赋小值。

这是BHAC非均匀AMR conserved数据reader，**不是任意GRMHD primitive dump reader**，也不是BHAC postrad格式通用reader。换数据只改文件名/NSPIN不一定足够。

已确认mounted footer多个快照的block size为16×8×8，γ=4/3、index3=.9375、index5=20。当前NSPIN=3有直接数据依据；index5为20也恰能解释为什么main新增spin警告重要。注释“5用于Mizuno GRMHD”说明曾有另一格式意图，应保留其限定语境，不应把当前数据作者／论文引用等同于那个格式。

现有 reader 对`fopen`错误使用`file_id<0`，而不是`file_id==NULL`，多数`fread/fseek/fscanf`也未检查返回值；坏路径或坏footer可能崩溃／读错，而不是清晰报错。这些从merge-base继承且本模型diff中main未修。重跑前至少应确保输入是预先验证的文件；后续工程修复宜纳入专门任务。

### 6.2 SFC和坐标网格不是可任选“参数”

SFC=1时按Morton顺序遍历level-one blocks，再递归读取AMR forest。若生产数据实际是普通i/j/k次序，用错SFC可能把不同空间位置的流体互换，而leaf总数仍能相等。`nleafs==igrid`只是结构计数检查，不是空间顺序正确性的证明。

当前96×48×48是根网格；block16×8×8意味着根block个数为6×6×6。更高AMR有效分辨率不能仅据根尺寸推定。hslope由运行目录的grid文件读取，不来自已检查的footer。

MKSBHAC用`r=exp(x1)`、`θ=x2+.5h sin(2x2)`；hslope=0使θ=x2。MKSHARM则含π与`1-h`映射，两个参数同名但定义不同。当前φ范围是0–2π，而不是看到文件xprobmax3=1就认为只模拟1弧度。当前metric公式相对fork祖先未有实质改写。

### 6.3 rcam、RT_OUTER、地平线margin分别控制什么

- **rcam=1e4Rg**：远处相机初始半径；观测角尺寸转换仍用SOURCE_DIST=16.8Mpc。不是把相机真正放在地球，也不是流体域外缘。其有限距离误差应通过rcam变化检查，但没有证据要求无条件改500或更大。
- **RT_OUTER_CUTOFF=50Rg**：关闭该半径外的等离子体转移并控制项目轴子积分调用范围。当前GRMHD域可到约2500Rg；“有外层数据”不等于“外层Faraday被使用”。EHT VIII §4.2明确说，在其近面向M87条件下，外边界只要r≳40Rg，其所考察的观测量对边界不强敏感；但同时明确不能排除r≳100Rg处额外Faraday物质，其效应不含在该模型中。这不构成本项目像素协方差／轴子拟合的误差上界。先保留50用于复现，增加100或更远的受控边界systematic，检验的是Faraday/吸收累计而非仅外层发射。
- **horizon_marg=.01**：`constants.c:28–29`实现为`r+=1+sqrt(1-a²)`、`CUTOFF_INNER=r+×1.01`。a=.9375时分别约1.3479853、1.3614651Rg，差约.01348Rg。它是相对margin，不是注释直观暗示的固定+.01；当前不应只凭“BHAC是horizon penetrating”就把停止半径移到地平线内。需要数值收敛与算法审查共同决定。
- **THETAE窗**：范围为`1e-3≤Θe≤100`，即约5.93e6–5.93e11K。EHT VIII一zone常见发射电子Θe约2–20，说明100并非离谱的上界，但实际尾部分布、Faraday贡献和cut交叠仍需检查；不保证每snapshot只有非发射区域受影响。

---

## 7. thermal活跃链与上游2023 pol_emission修复

### 7.1 主链究竟用哪一分布

当前header `DF=TH`；偏振step在`pol_rte_integrator.c:635`直接调用`evaluate_coeffs_single`，该函数调用`j_I/j_Q/j_V/rho_Q/rho_V/a_I/a_Q/a_V`，后者在`src/pol_emission.c`中由`#if(DF==...)`选择thermal分支。

同文件中 `evaluate_coeffs_user` 虽定义了thermal+kappa混合，并用

`eps=(1-exp(-beta^-2))*(1-exp(-(sigma/sigma_min)^2))`

这样的权重，但**当前偏振step没有调用它**。`sigma_min=1`是混合权重平滑开通尺度，不是删除流体的cut；`USERSPEC=1`又是另一独立开关。不能把这三个数字／宏混为一谈。

### 7.2 按最终diff归纳，不能只看中间commit标题

本次核查了提交的原始diff：

| 上游提交 | 实际改变 | 当前thermal主链影响 |
|---|---|---|
| `6e30fa5`（2023-04-13） | 新增DEXTER开关和Shcherbakov型rhoV fit；最初引入错误预处理行`#(DEXTER)` | 更新后的默认rhoV改变，会影响主thermal Faraday；旧fork尚无这段语法错误。 |
| `4d08208`（2023-04-17，"bug in pol_emis"） | 把`#(DEXTER)`修成`#if(DEXTER)` | **修的是6e刚引入的预处理语法，不是证明fork现有thermal有同一数值bug。**不能脱离前置提交单独套用。 |
| `6782822`（2023-05-11） | `get_w(theta_e)`从theta_e改为`(kappa-3)theta_e/kappa` | 只影响kappa模型；当前TH不调用get_w。 |
| `d84589d`（2023-07-24） | thermal `rho_Q`、`j_Q`翻号 | 是中间历史状态，不能当最新修复摘取。 |
| `eeb9ff3`（2023-07-26） | **撤回上述thermal两次翻号**；修kappa/power的jQ/aQ及kappa aV | main最终thermal jQ/rhoQ符号与fork祖先相同。非热修复当前不活跃，未来切DF才必须处理。 |

`git diff merge-base→main -- src/pol_emission.c` 的最终thermal实质变化集中在`rho_V_thermal`的fit选择；没有保留d845对thermal jQ/rhoQ的翻号。这是“上游未同步”的事实，不是本项目写错了thermal符号的证据。

### 7.3 thermal Faraday fit确实变了，不能说后来修复都不影响

旧fork：

\[
\mathrm{fit}_{old}=[K_0(1/Θ_e)-\Delta J_5(X_e)]/K_2(1/Θ_e).
\]

main默认 `DEXTER=0`：

\[
\mathrm{fit}_{new}=[K_0/K_2]\,[1-0.11\ln(1+0.035X_e)],
\]

并在K2不正时为K0/K2提供fallback；`DEXTER=1`才选旧式表达。其余rhoV前因子包含`n_e B cos(theta_B)/ν²`形式，fit变化会进入当前Faraday传播，不能被“DF=TH所以不受非热修复影响”一句掩盖。

建议：

- 历史复现baseline保留旧fit并明确记录。
- 科学重跑如果采用main的fit，需要把它当作可验证模型更新，与旧fit做受控对比，尤其冷电子／大Faraday深度像素。不要悄悄替换后把差异归因sigma。
- 若同步main但希望保留旧fit，应显式定义DEXTER=1；仅拷贝main的pol_emission到没有DEXTER宏的旧header下，未定义宏在`#if`中按0处理，会**静默选择新fit**。
- 改用KAPPA/POWER时，还需同步get_w和非热符号补丁，并检查θB≈π/2附近sign patch的数值处理及完整偏振约定。单摘函数符号不保证一整套代码约定一致；本次不重做polarized基底／积分器审查。

---

## 8. 文章模型／图4应如何解释这些参数

已读取指定原文件：

- `/home/canoziia/MyProject/PHYS/EHTaxion_imaging_pipeline_constraint/main.tex`，重点355、367、405–418行。
- `/home/canoziia/MyProject/PHYS/EHTaxion_imaging_pipeline_constraint/old_images/old.tex_copy`，重点140–162、189–198、345–348行及相关图描述。

### 8.1 当前稿的active叙述与诊断比较

`main.tex:355`明确写：MAD、aJ=.9375、包含emission/absorption/Faraday并归一到0.6Jy。`:367`区分“有／无plasma传播效应”、赤道薄层发射、lensed photons与smearing等比较。这些是受控诊断的科学意图，不能因为源码留着关闭系数的注释，就把正式full-plasma主run认作无吸收或无Faraday。

同样，eq_fakevalue等专用目录若用于人为赤道发射／去传播的教学或因果分解实验，其不同cut未必是正式主模型的错误；但其结果不能无标识混入main的full-MAD物理样本。主run ELF=1、诊断ELF=1000的区别应纳入图表来源manifest。

旧稿`:140/160`的“230GHz optically thin”主要解释toy equatorial disk；旧稿`:189`明确分四种场景：解析toy、加入GR的toy、加入盘厚度、最后才是MAD GRMHD加完整plasma，并写了100个时片／330天；旧稿`:348`又描述选用GRMHD较厚／wash-out显著。两种语境不能拿来推导“正式模型一定无吸收”，也不能反过来断言“光学厚”已经由当前tau图验证。当前完整plasma模型与toy diagnostic需要明确分栏；旧稿的100／330与当前稿占位999／1234也不能冒充同一批数据的已验证统计。

### 8.2 Fig4不是sigma cut稳健性的替代证明

`main.tex:405–416`把Fig4定位为GRMHD快照引起的轴子振幅／相位理论不确定性；文中快照数999和跨度1234天仍有明显待修改标记。旧稿`:193–198`也是强调GRMHD本征EVPA变率与发射位置变化。这里的“std”是所选电子配方、sigma cut、网格、flux尺度等**条件下**的时间样本统计。

因此：

- 不能把固定cut=1样本的时间std叫作已经包含cut1→1000的系统误差。
- Rhigh、temperature closure、Faraday fit和RT外边界变化会改变被采样的背景过程，不能只给同一组Fig4误差条。
- EHT VIII的cut10评分分布测试不等同于本项目Fig4每像素振幅／相位、非高斯尾部及协方差稳健。
- 重跑后应同时更新论文中实际模型配方、平均flux目标及校准办法、快照样本及时间尺度、是否包括上述systematics；本次不重新分析图4数据。

---

## 9. 如何证明旧binary用了哪个header：当前能说与不能说

makefile显示core源通过`RAPTOR/src`取得，模型文件和`definitions.h`通过run目录提供；`setup.sh`还会把模板复制到run目录，并根据选项改metric、SFC、BPOL等。它默认选项包括CKS和nosfc，故不应对当前MKS/SFC=1配置直接无参数重新setup。

重要区别：

- 编译命令/include路径说明**若重新完整编译**会搜索哪里；旧binary的当前行为不由当前header文件自动决定。
- header变化未必使已有object自动重编；本makefile的object规则只显式依赖`.c`，没有完整header依赖生成。这是需交由build审查落实的可追溯风险，本报告不重复构建实验。
- 保存binary中`get_fluid_params`的实际比较常量，配合该函数反汇编定位，是确认编译后cut的直接证据。已由父审查确认主run=1、eq_fakevalue=1000，本次采用该结论。
- **实际常量=1不能唯一证明具体哪份文本header的历史路径**：LR、CC、CR内容相同，甚至历史上其他路径也可能相同。若没有DWARF、编译日志、include trace、object／header归档，最多证明“行为与cut1配置一致”，而不能伪造完整build来源。
- 科学复现真正需要的是重建manifest：core commit和本地diff、run模型/headers/grid/input的SHA、实际宏表、编译器/flags/include依赖、binary SHA以及运行目录。父agent已专门build审查，本次不重跑。

---

## 10. 推荐的重跑配置决策表与验收标准

| 决策项 | 初始fiducial建议 | 首要systematic／何时更改 | 当前证据不足 |
|---|---|---|---|
| sigma cut | 主模型沿用1 | 10为首要文献对照；1000为历史／漏斗污染压力测试，有额外物理时可另立模型 | 当前各σ区间对I/Q/U/V、Faraday和轴子拟合的贡献未量化 |
| 电子分布 | TH，完整吸收和Faraday开启 | 若研究非热/对产生，另立可验证物理模型并处理上游非热修复 | 不能仅切DF或cut就代表PIC／pair physics |
| Rlow/Rhigh | 若先复现历史则1/1且注明benchmark | 增加较冷电子分支；取值按flux、形态、偏振和多波段一致性筛选 | 1/1不是已验证最佳M87解 |
| Θe闭合 | 历史配方可作为复现基准，但必须明确公式 | 与EHT式或数据来源一致的电子加热配方比较并再校准 | 无法宣称现公式严格等于EHT/Mizuno |
| thermal rhoV fit | 历史旧Dexter fit明确记录 | main默认Shcherbakov fit做受控更新对照 | 未知对本样本高Faraday像素的变化量 |
| MBH/D | 复现先保持6.2e9/16.8Mpc | 更改尺度必须同步时频/角尺度/单位和flux归一化 | 不宜只改标签为6.5 |
| MUNIT | 3e25仅作校准初猜 | 每模型、同时间窗求一个常量使平均flux达声明目标 | 0.6Jy未获本次数据验证 |
| NSPIN/SFC/grid/hslope | 当前3/1/96×48×48/0，只对当前数据成立 | 换数据按生产参数／footer重新适配 | footer不足以单独证实SFC/hslope |
| rcam | 1e4Rg用于基准 | 有限距离／时间相位收敛 | 未测变化量 |
| RT_OUTER | 50Rg用于复现 | 100或更远边界与外层Faraday评估 | 不能凭内区发射紧凑忽略传播 |
| Θe窗 | 1e-3–100用于复现 | 分析被删区域的发射／传播贡献、必要时有依据扩大 | 未检查温度尾部与floors交叠 |
| horizon margin | .01相对margin用于复现 | 在算法许可下缩小margin做收敛 | 不指定无条件唯一更优值 |
| 230GHz单频 | 对应本文目标 | 多频数据须真实生成多个频率，不只改每decade参数 | 单频没有检验RM频率行为 |

建议验收分三层：

1. **配置验收**：输入字段完整且reader成功读取；metric、footer spin、网格、单位、cut、DF、温度公式／窗口、Faraday fit均记录；主模型与diagnostic输出物理隔离。
2. **数值验收**：同模型无轴子基准下finite/positive流体，vacuum安全；守恒量转换／场归一性／grid lookup无异常；步长、边界、相机、分辨率收敛。此层可结合其他agent已完成报告，本次未运行。
3. **科学验收**：恒定MUNIT下平均flux与时间变率、I/Q/U/V和偏振形态、Faraday深度等与目标一致；cut／R配方／rhoV fit／外边界systematic传递到Fig4与最终约束。不能用一张“像M87”的图或总flux单指标替代。

---

## 11. 一手来源与本次证据文件

### 11.1 论文

1. **EHT Collaboration 2019, Paper V: Physical Origin of the Asymmetric Ring**，arXiv:[1906.11242](https://arxiv.org/abs/1906.11242)，ApJL 875 L5。重点Eq.(7)/(8)、fn115、§7低密度数值处理。缓存 `/tmp/eht-constant-sources/1906.11242.{pdf,txt}`。请勿误把1906.11241（Paper IV成像）当作Paper V。
2. **EHT Collaboration 2021, Paper VIII: Magnetic Field Structure near the Event Horizon**，arXiv:[2105.01173](https://arxiv.org/abs/2105.01173)，ApJL 910 L13。重点§4.1/4.2、Eq.(15)/(16)、σ=1/10测试、每模型独立密度校准、§6.4。缓存 `/tmp/eht-constant-sources/2105.01173.{pdf,txt}`。
3. **Mizuno et al. 2021, Comparison of the ion-to-electron temperature ratio prescription: GRMHD simulations with electron thermodynamics**，arXiv:[2106.09272](https://arxiv.org/abs/2106.09272)，MNRAS 506,741–758，正是项目references.bib中的`Mizuno:2021esc`。重点§2.1/2.3、Eq.(7)/(8)、fn2、σcut和Fig4。新只读下载缓存 `/tmp/eht-fluid-model-review/papers/2106.09272.{pdf,txt}`。
4. **Bronzwaer et al. 2020, RAPTOR II: Polarized radiative transfer in curved spacetime**，arXiv:[2007.03045](https://arxiv.org/abs/2007.03045)，A&A641 A126。重点§2和Appendix C系数来源／thermal分布；并非对当前axion fork算法自动背书。缓存 `/tmp/eht-fluid-model-review/papers/2007.03045.{pdf,txt}`。

本次两次新下载均使用curl connect timeout=10s、max-time=50s、无重试，外层命令timeout=115s；没有无超时网络搜索。

### 11.2 审查产物

- 本报告：`/tmp/eht-fluid-model-review/REPORT.zh-CN.md`。
- 四处逐文件路径／SHA／宏／输入／Run.py证据：`four-location-model-evidence.json`；报告完成前再次读取逐一校验，所记录四处模型文件SHA全部仍匹配。
- fork实际贡献：`fork-owned.diff`。
- 上游后续模型／pol_emission变化：`upstream-later.diff`。
- merge-base/fork/main相关文件SHA：`git-hashes.txt`。
- 指定文章原路径已列§8；mounted footer／尺度沿用结果：`/tmp/eht-mounted-constants-qu_7_khw/comparison.json`。

**最终判断**：当前主运行cut=1、热电子、完整plasma开启有明确源码、已复核binary以及相关一手文献依据；并无理由把它因LC模板1000而自动判错。真正尚未完成的是这套具体电子温度闭合／R=1／MUNIT3e25／RT50组合的flux与像素级稳健性验证。以cut1为可追溯fiducial、cut10为首要systematic、cut1000为有标签的压力测试，是目前最有依据的重跑起点；任何更强结论都需要额外数据检验，而不是仅凭默认header或注释。
