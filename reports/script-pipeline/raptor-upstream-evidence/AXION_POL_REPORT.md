# RAPTOR 轴子/偏振核心只读审查

## 0. 范围、证据和结论

本报告专注偏振/轴子核心、调用链、沿 ray 积分和相关 run/test provenance，不替代父审查的完整文件 manifest、模型流体或 header 全面审查。没有修改挂载、GitHub 克隆或项目源；没有运行任何旧 RAPTOR binary，没有读取并重算真实 cube。只在 `/tmp/axion-pol-audit` 写报告、复制少量对象作静态检查、编译执行毫秒级合成 C 探针；下载检索文献也只写 `/tmp`。

四处必审路径：

- **L-core** `/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel/raptor`
- **L-run** 同父目录 `/raptor_run`
- **G-core** `/tmp/eht-github-raptor-GckNId/raptor`
- **G-run** `/tmp/eht-github-raptor-GckNId/run/raptor_run`

此外仅核实同父目录 `raptor_run_test`、`raptor_run_test_eq_fakevalue` 的相关源码/对象来源，以及 EHT 目录下其它同名 core 是否可以作为来源；不把它们的物理开关算进主代码。

### 总判断

1. **主 core 不是 j-only，也不是薄赤道假发射模式**：关闭吸收/Faraday、赤道截取、常数假发射、早时发射关闭、pixel 编码全部是注释。主分支实际计算并使用热同步辐射发射、吸收、Faraday rotation/conversion，保留 GR 平行输运。
2. **主 core 也不是图像末端乘一次相位的纯 screen**：沿 ray 每个进入 `pol_integration_step` 的 cell 先等离子体传输再施加局部轴子旋转；新旧辐射一起旋转。因此可含多处发射/不同相位叠加和等离子体反馈。但它是算子分裂近似，并非已经验证的严格联合 full 求解。
3. **有 active 且必须解决的问题**：`&&→||` 把无效 fluid 也送进 plasma step；新端点 `path_counter+1` 首次读未初始化的 ray 记录且与向相机推进不同边；NormR 无并发安全初始化/读入和边界检查；质量 alpha 与时间频率混用；BL 轴子场被直接在 KS 时间/方位上计算；固定 25Rg window 会改变物理场。
4. **独立 test 的确存在 active 物理筛选**：`raptor_run_test` 保存的偏振对象包含高度转整数后条件性清零全部 11 个系数的代码，主对象没有该块。另一个 fakevalue 目录当前偏振对象反而与主对象字节相同。不能靠目录名恢复历史实验开关，更不能把它们混成生产 core。
5. **本轮不认证现存 cube 的历史正确性或污染幅度**。源码、对象、Notebook 的证据等级不同；没有完整构建 provenance，不能由今天源码断言每个旧产品当时用了哪组开关。论文 3%/5% 等叙述也不能替代数值回归。

本报告位置采用 **G-core 当前文件行号**；L-core 在 800 行之前大致一致，之后因注释增加约偏移 13 行。关键对应：G `824/960/973/997` 对 L `837/973/986/1010`。完整 scoped hash/token 证据见 `identity.json`，编号源码见 `*.numbered.txt`。

## 1. 正确的上游基线：两种 diff 不能混为一谈

上游完整克隆 `/tmp/eht-raptor-upstream-audit-lVJjW7/upstream`：

- `main = b91e937b7b746eb24992b6c1dd8831c5dea0cef1`
- `project-main = b6df33c69e650f3617845f4324ba642bdb878582`
- 本轮重新 `git merge-base main project-main` 得到 `cdc2547b868ed8d6ae0c9a5428bc1ea86fce901a`。

### 1.1 自定义改动：merge-base → fork

见 `axion-core-base-fork.diff`。本次关注的 custom 包括：轴子参数读入/传递、两个 plasma solver 的 delta rotation、NormR 查表、场表达式、endpoints 和 `&&→||`，以及大量注释实验块。`gr_integrator.c`、`GRmath.c`、`pol_emission.c`、`tetrad.c` **对 merge-base 没有项目改动**。因此其继承问题不应归咎于本项目 custom。

本地完整 fork 历史还显示：

- `d11f522 feat: add functions for pol` 首次加入这套轴子逻辑，已经是 `||`，关闭 plasma/薄层块当时也为注释。
- `6b1d2b1` 改读入循环 `38×490 → 23×2010`，但当时静态数组仍 `50×500×3`，这个历史中间状态有越界；**不是当前最终源码**。
- `9831508` 改到 `23×2400`，短暂用局部 int 作为 static 数组维数；后续修成宏。
- `1b72ce6 fix conatant`、`b6df33c fix list out of bounds` 完成当前宏和半径搜索上界。
- `b6df33c` 只给半径 j 搜索增加上界，没有修复质量 i 搜索上界，也没有修复 ray endpoint。

不能拿 commit 名“fix list out of bounds”当成全部越界都已解决的证据，也不能把中间历史状态直接说成现存 cube 已使用的 binary。

### 1.2 上游 main → fork 的额外差异

见 `axion-core-main-fork.diff` 和 `upstream-only-pol.diff`。

在本轮三个主要文件里，merge-base→main 只有 `pol_rte_integrator.c` 的后续上游改动：commit `6782822 adding user emiss for pol` 加入

```c
#if (EMISUSER)
    evaluate_coeffs_user(...);
#else
    evaluate_coeffs_single(...);
#endif
```

fork 保留 merge-base 的无条件 `evaluate_coeffs_single(...)`。所以 main→fork 看起来“删除 EMISUSER 分支”，**并非我们的 custom 故意关掉了混合电子模型**，而是 fork 未合入上游后来更新。当前 G-run/L-run `EMISUSER=0, DF=TH`，实际基准物理是 single thermal；若将来启用 EMISUSER，不能只改宏就以为它生效，需合入/实现调用选择并另外验证 user 分支。其它上游新增模型流体/header 差异留父审查，不在本报告重复。

## 2. 四处代码身份与实际调用链

### 2.1 L-core 与 G-core

- `core.c/functions.h/gr_integrator.c/metric.c/camera.c/main.c` 字节相同。
- `pol_rte_integrator.c` 字节不同，但去掉 `//`、`/*...*/` 后词法 tokens 相同；本轮独立重核。
- 差异是格式/注释：hslope/R0 打印文本，赤道注释阈值 `.01→.1`、注释高斯 ring emissivity，以及注释 `GRMHD_cut_400` 输出等。**本地没有额外 active j/a/r 关闭开关**。
- G-core pol SHA256 `b58044a7cbafdfe509df5b1b612c9cf9c5d5101628c00ffb755aa7f6db8a55c4`；L-core `2e401fea06031f922e88432a255edd63996e32d9dac3585a540976aa75fe13d9`。

### 2.2 L-run 与 G-run

本轮涉及的 `makefile/run.sh/definitions.h/NormR.csv` 字节相同。主要 wrapper 差异：G-run `Run.py:96` 扫 `[.3,.25,.2,.15,.1]`；L-run 当前只扫 `[.1]`。这是实际 active 扫描选择，不是新物理公式。两者 `Run.py` 参数字典虽默认 `AXION_OMEGA=.4`，`__main__` 循环覆盖为上述列表；直接 import 调用 `Run` 时是否覆盖则由调用者决定。

G-run `model.in` 仍是 **旧无轴子字段的 Sgr A* 示例**：MBH 4.14e6、i50、200²、230GHz、步长.01，然后直接 MAX_LEVEL；L-run 该文件当前不存在。当前 `core.c` 在 MAX_LEVEL 前要求三个轴子字段，直接使用这个 G-run 示例会错读：AXION_NORM 会把 MAX_LEVEL 的数值1吃掉，后续扫描 EOF；没有检查 `fscanf` 返回值。因此 **不能直接用仓库 model.in 或 run.sh 当生产基准**。run.sh 还指向 `harm3d.txt`，不是本 wrapper 的 BHAC M87 快照。

主 `Run.py:62–82` 的参数：MBH6.2e9 Msun、距离16800kpc、M_UNIT3e25g、R_LOW=R_HIGH=1、i163°、100²（高分辨率默认关闭）、FOV46.2261Rg、230GHz、STEPSIZE .005、norm0 baseline 和 norm.1 的16个相位，max_level1。M_UNIT 的流体后果/总流量验证留父审查；这里仅记录传到 ray/core 的值。草稿用了 M87 6.5e9Msun，不能在频率↔天单位转换时静默混用6.2和6.5。

`f.write` 没有换行不宜凭视觉直接判错。本轮只提取 f.write 表达式构造 `/tmp/generated-model.in`，合成 C 用同样 `%s %s %lf` 扫描，17项全部成功：数值后紧接下一个字母标签仍可作为扫描边界。**因此没有报告“当前 wrapper 因无换行必然错读”的假阳性**；但应改成标准换行文件并验证标签/字段数/范围，避免依赖这种脆弱语法。

### 2.3 核心调用链与哪些参数真正控制物理

```text
main.c: read_model(argv) → init_model() → set_constants()
       → init_camera / prerun_refine
       → calculate_image_block(camera, frequencies, block)
core.c: 每像素 OpenMP parallel for
       → malloc(9*max_steps) lightpath
       → integrate_geodesic(alpha,beta,...,&steps,...)
       → 每频率 radiative_transfer_polarized(...,axion_omega,norm,phase,block,pixel)
pol_rte: p=steps-1 ... 1
       → get_fluid_params
       → pol_integration_step（如当前条件满足）
           pitch/单位变换 → single thermal j,a,r → plasma tetrad
           → NormR/场值 → f→Stokes
           → stiffness → RK4 或 trapezoid → 每步 delta rotation
           → pol_frac guard → Stokes→f → tau/tauF
       → rk4_step_f 正向平行输运
       → 相机 tetrad → IQUV * frequency³
```

轴子参数仅传入 `#if(POL)` 的偏振路径，POL=0 时本 fork 不模拟轴子作用。`block/pixel` 当前主代码只供被注释的诊断输出/编码使用，不改变 active 辐射。`TIME_INIT` 读自 argv[3]；在此链中不是轴子场的时间项，wrapper 用它当输出序号：0为baseline，1..16为相位产品。轴子相位来自 `phase_index`，ray 时间来自 `X_u[0]`。这不是在 core 内把 GRMHD 演化时间和16个轴子相位同时推进的 slow-light 模拟；固定快照上的沿 ray time delay 被包含，快照演化由外层换文件实现。

`DEBUG=1` 是 header 的诊断宏，主 pol 里没有以 DEBUG 自动关闭物理的分支。`LIGHT_TRANSPORT=1/RAD_TRANS=1` 在本次 src 检索没有实际引用，不能把切这个宏当成有效物理切换；父 header 审查可继续核对其它模块。相机新增 printf 只是尺寸日志。

## 3. 每类 custom 改动：旧→新、意图、状态和处置

### C01 参数/原型/调用新增（active，保留但显式定义物理身份）

位置：`core.c:21,66–69,100,111,137`，`functions.h:106–110,335–337`，`main.c:139`；pol `588–595,919–922`。

旧：只有 MBH/M_UNIT/TIME_INIT 等，`radiative_transfer_polarized(...,tau,tauF)`，`calculate_image_block(...,frequencies)`。
新：增加 `axion_norm,axion_omega,phase_index`；传 block/pixel；新增 X_u_next。

意图：在固定流体快照上扫描云强度、质量/频率、初相位。物理扩展应保留，不应作为 debug 删除。必须将混名 `axion_omega` 分成质量查表 alpha 与 Re(Mω)，定义 `norm = 2π g a_max = 2π b g f_a`（当前 Notebook 也标 bc），phase用 radians 或明确16等分索引。读入需检查返回值和参数合法性，输出 metadata 中记录这些值及源表hash。只打印 omega/norm/phase 不足以追溯完整配置。

### C02 RK4/trapezoid 末端增加 delta rotation（active，科学核心，不应直接删除）

位置：`pol_rte_integrator.c:409–505,507–573`。

旧：solver 得 x1..x4 后 `S_A=(x1,x2,x3,x4)`。
新：

```c
delta_a = alphaf - alphai;
Inew = x1; Vnew = x4;
Qnew = x2*cos(2*delta_a) + x3*sin(2*delta_a);
Unew = -x2*sin(2*delta_a) + x3*cos(2*delta_a);
```

即 `Pnew=(x2+i*x3) exp(-2i delta_a)`，**EVPA变化为 −delta_a，不是 +delta_a**。这是 Stokes 的2χ旋转，因子2本身正确。阈值 `if(fabs(alphaf)>0.5e-7)` 和原始 identity else 都是注释，active旋转无该阈值；保持小场连续性是合理的，不要重新启用该不连续阈值。

纯旋转保持 I、V、Q²+U²，probe `delta=.2` 得 Q=.921060994003、U=−.389418342309、χ=−.2；没有直接生成V/吸收。**不能把这个数学旋转矩阵孤立判成符号错**：当前 alphaf 取 ray 的更深端 p+1，alpha_i取p，正向为 p→p−1；向后差分配上负旋转在连续极限给 `dχ=+d(ga)`，与草稿 `rho'_V=rho_V−2d(ga)/ds` 的矩阵一致。但命名误导且离散错配，见C08。**修 endpoint 时必须同时规范 delta定义和旋转符号；只把 +1 改成 −1 而保持全部符号不变会反转物理约定**。

两个 solver 都是先完成 plasma/source 步再一次性旋转新旧辐射：这是 Lie splitting，不是把 axionρV统一放进4阶RK stages或隐式matrix。本轮从源码原样抽出 solver，构造常数 polarized source jQ=1、线性轴子相位、长度1：对联合解析解误差 n=8/16/32/64 为 .105275/.052603/.026297/.013148，约每加倍网格减半。故即便叫 RK4，**含轴子的源项分裂整体也只呈一阶**；没有源项且只有同轴Q/U旋转时局部指数则可精确。若有Faraday conversion/二向色吸收，矩阵一般不对易，也有分裂误差。

建议：保留指数旋转作为精确 axion 子算子；主 full 可用正确端点的 Strang splitting/联合求解，并做步长收敛；旧 Lie 方法留 `legacy_split` 显式开关以复现旧数据，不能无验证宣称 full 四阶。不要同时重新打开注释 `rV+=...`，否则重复计入轴子。一个清晰的forward约定是令 ΔA=A(x_{p−1})−A(x_p)，再用 Pnew=Pplasma exp(+2iΔA)；或保留当前负指数并改用相反定义的delta。局部plasma/frame与真空输运如何分裂仍须整体验证，不能只机械替换索引。

### C03 局部单位变量和 k2（active计算但无 active 物理用途）

位置：`621` 新 `double dlambda=*dl_current;`，`693` 新 k2。

旧：保留 k_old并缩放 k和步长；新：额外保存未缩放 dlambda、计算 `sqrt(k0²-k1²-k2²-k3²)`。

这两个新变量仅用于注释的 rhoV/printf；编译优化通常去掉它们。k2不是曲时空光子可用的Minkowski范数，log-r坐标下更不能这么算；可能产生NaN。主full应移到显式诊断实验里或移除无用计算，**不要启用它来“修正”轴子系数/路径长度**。现有k、dl能量缩放以及C_CONST是继承实现，本项目没有另改；轴子delta已经是无量纲角度，不应再乘CGS频率或L_unit。

### C04 新 NormR 静态读入（active，保留查表能力但必须重构初始化/验证）

位置：`654–683`。旧无轴子表；新静态 `double Axion_Radial[23][2400][3]`、首次 plasma step 从工作目录 `NormR.csv` 顺序 fscanf。

已确认 L-run和G-run表均为 55200行=23个alpha(.09..31)×2400个半径(1..240.9,dr.1)，字节一致，并与 `radial-function/reference/result_a=0.9375_l1m1.txt` 完全一致；SHA256 `0a352558f8b846f2481ab2acf35c2d6f8b30ce1f83f2e6f5196660f23433c353`。

问题是确定源码事实：

- `fopen/fscanf` 均不检查。缺文件会对NULL FILE操作；短表/坏列会留下静态零值并触发零分母/假场。
- `flagg` 与静态表在 `core.c` OpenMP像素线程中无同步首次读写，是C数据竞争。多个线程可能重叠改表；某线程看到flag=1时其它线程仍在写表。即使大家读同一文件、最后数字通常相同，也不构成语言级安全保证。
- **norm=0 baseline仍会读表、做查表及非有限运算**；0乘NaN不是0，所以“关轴子”并不完全隔离表问题。
- 维数/自旋/单位/归一化都不存于输入metadata；只换同名CSV无法保证兼容。高alpha .99成品为15个alpha(.30..44)而不是23个，不能直接覆盖当前NormR让这段代码无修改读入。

重跑前：单线程在ray parallel区之前完整读入不可变表，动态或验证维数、列数、单调网格、finite值、视界/单位/自旋和hash；norm0可直接绕过轴子计算。不要在每条ray首次交互里懒加载。

### C05 质量/半径两维插值（active，域内公式可保留；边界必须修）

位置：`699–730`。旧无插值；新在质量、半径上线性插值后再次插质量。

内部在所有slice共享同一r网格的条件下是正常双线性插值。此NormR满足共享网格，不应随意用另一个近似径向形状替掉。边界问题：

- i从1开始，循环没有 `i<omega_count`。参数>.31时访问i23；独立合成23×2400表、原式搜索，UBSan明确报告 `index 23 out of bounds`。当前扫描 .1..3 在域内；字典默认.4/高alpha扫描不在域内。也不检查NaN/低alpha，低端会默默线性外推。
- 当前 r 循环 j从1，条件 `table[i][j][1]<r && j<=r_count-2`，循环后--j。对合法i，这次上界修复实际上使最大j=2398，因此读取 j+1=2399在界内；**不能误报当前这个r loop仍必然越界**。但r超上限仍用最后两点线性外推；不是明确置零，也未评估外推误差。低端同样外推。
- 两个质量slice复用由slice i得到的j，依赖所有r网格一致；应验证而非默认。

建议full用域外报错/物理明确的外推策略，质量超过上限立即拒绝；支持动态网格时逐slice定位，分母非零，查询值finite。b6的半径guard是应保留的防护，不是应删除的debug。先检查索引范围再访问数组更清晰；这里当前j-loop虽然把数组访问写在&&左侧，但合法i时j最多到2399，尚未超过2400，不据此额外报越界。

### C06 新场坐标、角向因子与时间相位（active，必须统一物理定义）

位置：`684–698,735–736`。新式：

```text
r = exp(X[1]); theta先做BHAC映射后被 theta=X[2]覆盖
phi = X[3]; t = X[0]
A(x)= norm/(2π) sin(theta) R(alpha,r)
      cos[-phi+(1022+t)*axion_omega+phase_index*π/8]
      exp[-(r/25)^15]
```

**alpha独立确认**：本轮读取 `radial-function/provenance.json`，原 `.1.tar` 解出文件的SHA256与provenance一致，再读原producer（不是只看已修库或上一轮报告）。`print_cable:mu=.01*i+Mu_min; Result[i][j][0]=mu`；`Function_R:M=1, omega=get_omega(mu,a)`；`get_omega` 对source第二列插值后乘mu。CSV第一列是质量 `alpha=Mμ`，不是频率。主 RAPTOR直接把同一个 `axion_omega` 用在cos(t×...)；**这属于 ω≈α 近似/实现不一致，不能把全部alpha标签改叫omega消除**。

原producer得实频率：spin.9375 alpha.1→.09985470495，.3→.296319073264；spin.99 alpha.4→.3904155，.44→.426552036870。源表未被本轮认定为绝对oracle，但质量和频率身份已证实。`reports/script-pipeline/RADIAL_ARCHIVES_AUDIT.md` 的同一结论与本轮独立证据一致。原producer取 `double omega=get_omega` 丢弃imag并用 `cabs(result)`，当前CSV仅提供实径向模长，不包含复径向相位；这是受控近似需说明，而不是完整复本征函数。

- **坐标风险**：主run metric=MKSBHAC，`metric.c:113–114`的真实角度是 X2+.5hslope sin2X2；轴子代码 active后赋值覆盖。当前run `grid_mks.in` hslope=0，因此这个覆盖在本run没有数值影响；不能谎称当前一定因θ映射错已污染。换非零hslope则确定错误，换MKSHARM/CKS更不兼容；full应使用metric-aware坐标转换。
- **更重要的时间/方位**：草稿 `main.tex:143–156,309–313` 和旧稿 `155–159` 明确BL分离场，轴子代码直接用BHAC的KS t、phi，无KS→BL变换。现有 `GRmath.c:260–309` 已给变换：tBL=tKS−F(r)，phiBL=phiKS−G(r)。若仍使用实径向包络，正确BL相位与当前相差 `G(r)−ωF(r)`，不是常数；以a=.9375,ω=.3，当前KS相位−BL相位在r3/5/10/20为 .577/.914/1.355/1.784rad。一个1022全局偏移不能消除此径向差。若有意定义KS toy field，必须明确它不等于草稿BL本征场。要实现完整复径向解还需同时处理ArgR，不能只做坐标变换就宣称完整KG解。
- `sinθ, m=1`：与211球谐主导近似相符；不是精确spheroidal角函数。当前保留为显式近似即可，不必因为简单就删。更高alpha/自旋要说明误差。
- **1022偏移**：没有本轮证据证明这个值来自哪一份历史相机距离/零点。它作为所有cell共享常数，只改变同一alpha的整体初相，16等分扫描且相位完全边缘化时通常不改振幅；不应断言它破坏所有约束。但绝对相位、与screen/不同alpha比较时必须同约定。固定rcam=1e4却写1022，尤其需要来源说明，不应硬当真实快照时刻。full应显式 `t_ref`，旧值留legacy配置。
- **phase_index/8*π**：0..15对应16等分整周期，数学正确；保存产品index=phase+1，不等于相位在物理时间中的天数。
- **norm/(2π)**：若norm=bc=2π b g f_a，则是ga_max，符合本草稿的g a FtildeF/2 convention。norm=.1给峰值局部角尺度.0159155rad≈.9119°（还乘R sinθ W）。不要与常见L∝g_standard aFtildeF/4时的g_standard/2公式混淆。main.tex现在引入b(alpha)；core没有另乘b，后处理若再乘/不乘都须知道norm已包含bc。

建议：full分开alpha和omega；明确BL场或KS场；角坐标metric-aware；相位零点参数化并统一所有分析；不要擅自改全局频率再保留旧响应模板声称完全修正。

### C07 场 window：`exp[-(r/25)^15]`（active，不是注释，需科学选择）

位置：`735–736`。旧无轴子；新用径向本征包络之外再乘window。另两个 `radiali=exp[-(r/20)^15]` 和均匀screen-like `cos(phase+.137373737π) exp[-(r/15)^20]` 是注释，不能混淆。

W(9)=.999999779，W(15)=.999529926，W(20)=.965427402，W(22)=.863316549，W(24)=.581533686，W(25)=e^-1，W(26)=.165143001，W(28)=.004196243，W(30)=2.04e-7。

意图可能是让场在有限积分域外很快归零、避免未积分真空云区，但**没有查到作者明确注释证明唯一意图**。这是改变物理场，不是无害性能优化；使低alpha大半径云主体人为消失，并制造额外window梯度。对仅在r<10发射、全程真空传播且积分正确的jonly，端点总导数可令外部形状不影响发射点最终旋转；所以不能只看云峰很远就断言亮环信号必然大变。但外部也发射/有conversion或二向色吸收/轴子传播被gate截断时不能用该抵消保证。

**full建议无人工25Rg window**，完整处理到observer或解析补齐真空段，扩展径向表/定义远场有效性；如果保留window为稳定数值近似，必须独立命名 `axion_window_radius/power` 并做无window、25、50等的信号收敛，而不是藏在表达式。历史复现保留旧window配置，不能直接删除并混用旧cube比较。

### C08 新ray端点和`&&→||`（active，最高优先级整改）

位置：`960,973–978`；L对应973、986。

旧：只取当前位置p，向前步长来自p−1；有效流体且r<cut才做plasma。
新：`X_u_next=lightpath[(p+1)*9+i]`，`if(get_fluid_params(...) || r<RT_OUTER_CUTOFF)`，将next传入求场。

**端点**：`gr_integrator.c:integrate_geodesic` 每轮先存lightpath[steps]位置/速度，存该次负向步长，再steps++。有效存储记录是0..steps−1；终止后的新photon状态没存到lightpath[steps]。transfer首次p=steps−1却读steps：通常malloc容量尚在，但内容未初始化；当steps=max_steps时还读出malloc范围。源代码确定，不能描述成“必然每条ray堆越界”，也不能因为容量未满就说安全。此问题本项目新增。

正向transfer的几何推进是p→p−1，轴子差分却取p→p+1。若所有点都连续做轴子步，当前单发射点p0之后累计χ为 A1−A_(p0+1)，而理想端点是 A0−A_p0（这里A=ga）。因此即便符号在连续极限抵消，也有一格错位及边界项。合成非线性A=[0,.1,.4,.9,1.6]，从p3到相机，当前公式给χ=−1.5而正确端点χ=−.9；linear A例子可偶然抵消掩盖错误。该例只证明离散逻辑，不是现存图像误差幅度。必须规范forward端点+delta符号成套修，ASan/UBSan测试steps=0/1/max_steps和最后record。

**gate**：`A&&B → A||B` 的结果：

- r<50时，即使get_fluid_params失败仍做pitch、frequency、tetrad、coefficients。
- 对域外/查grid失败，model有早return0且不清空全部modvar，可能沿用上一个cell的B/U/n_e；初始modvarU/B全0则tetrad/频率退化。对sigma/T cut，model只n_e=0返回；正常也不是代表“无轴子”的真空。
- r>=50而fluid返回true时也可进入；在当前model后端r>cut通常返回0，所以不能宣称本后端必然计算cut外plasma。
- 可能意图是让轴子在无等离子体的r<cut区仍传播，这个需求合理，但调用完整plasma step不是正确实现。
- 即便保留||，pitch完全对齐时函数前面的继承return也会跳过轴子。若换回上游&&，则所有无fluid区都跳过轴子。故 **不能只回滚&&就称轴子full修好**。

正确分离：valid-fluid∧plasma-domain才计算j/a/r；axion独立根据场域沿ray更新，包括真空、sigma-cut、pitch-aligned段；GR transport照常。无fluid时不能构造来自陈旧B/U的plasma tetrad，可在平行输运/observer-safe screen basis旋转f或用稳定tetrad处理。

### C09 所有实验/注释块的处置（当前主代码均不active）

| 位置G | 旧→新增注释表达式/作用 | 科学用途、full处置 |
|---|---|---|
| 647–652及805附近 | `if(X[0]<−10023) jI=jQ=jU=jV=0` | 按传播时刻截发射、可能隔离次级/远端源；不属于full。可保留具名 emission_time_gate diagnostic，并记录阈值依据；不是合法的泛化“只保留直接光子”方案 |
| 497–504,563–572 | `if(abs(alphaf)>5e−8)`/else不旋转 | 不连续小场阈值，full不要启用。只是历史尝试，不要误报当前有阈值 |
| 731–732 | R替换成20Rg平顶指数窗 | toy场/端点校验；full不用，另设toy profile |
| 738–740 | 均匀角向/径向15Rg窗、固定相位偏移 | 常量screen/校验总导数近似，full不用。若恢复需与真实mode分开命名 |
| 742–756及其它printf | thi_cut_profile/log/GRMHD_cut | 诊断可有科学价值，保留显式输出开关；多线程应分文件/有序写，不直接并行append混行 |
| 756–764 | `rV+=.../dlambda/k2`、`rV+=norm exp(-r)cos(...)` | 早期将轴子塞入Faraday尝试；当前不active。前者重复norm/(2π)、错误k2/单位路径，不能照抄激活；若做联合solver应重新从rhoV=−2dA/ds及invariant单位推导，禁与delta旋转双计 |
| 766–774 | `rQ=rU=rV=aI=aQ=aU=aV=0` | jonly有意义：隔离等离子体传播与发射分布叠加。full必须off；Fig4比较应保留为明确diagnostic mode，而非删掉实验能力 |
| 776–804 | `if(abs(r cosθ)>.01) j,a,r全部0; else jI=1e−30,jQ=jU=3e−31,jV=0,a=r=0` | thin equator/fake source实验。既删去离赤道发射又改变j源，不能把整块称eq_jonly。应拆成 emission mask 和 source model 两个独立开关；`abs`是int函数必须改fabs，否则阈值根本不是.01Rg |
| 811–814 | `S0=block*64+pixel; S1=S2=S0/10; S3=0` | 像素索引/图像布局debug；full关闭。当前tot_pixels=100，*64也不是当前正确像素编号；不应用于物理图 |

本地注释中的 `_ring=5.5,_sigma=1,_scale=exp[-(r-ring)^2/(2σ²)]+.01` 和 `jI=1e−30*_scale,jQ=1e−31*_scale,jU=0`，是径向假发射形状对比，适合fake source，不是“被注释所以必须永久删除”的垃圾。保留版本化diagnostic配置可支持文章，禁止在生产中靠编辑注释切模式。

### C10 其它本范围自定义

- `camera.c:35` 新printf block/pixel尺寸：active但无物理效应，移到日志verbosity即可；full可保留/默认关闭，不需重算。
- `metric.c` custom仅空行/空格：没有活跃metric表达式改动；别把diff数量解读为项目改变Kerr几何。
- `io.c:21` spec_folder[64]→[128]：active缓冲区长度，不改变物理；保留。完整I/O和输出manifest不重复父审查。
- 删除/增加解释注释、空白不改active逻辑；本地/remote token一致证明这里没有隐蔽本地物理关闭。

## 4. 继承实现：不是custom，但重跑需要知道的数值风险

这些在merge-base和上游main中已存在；不得列成项目故意关物理。

1. **stiffness检查** `363–407`：使用sqrt而不是csqrt处理complex z；实sqrt接收负实部会NaN。其稳定多项式按正tau写，阈值.99；还缺实际solver中的C_CONST。对物理aI≥0，正特征值常让它总判STIFF。probe零系数/纯ρV/纯吸收均返回1，零系数也走trapezoid。不能照注释说“低Faraday一定RK4，强Faraday隐式Euler”：实际隐式是trapezoid，且切换不符合理论描述。
2. 上游RAPTOR II, Sec3/AppendixA明确稳定性依赖完整 Δλ×矩阵特征值；现在代码的物理单位/符号/complex处理需验证。主full至少记录solver使用率、用常系数解析测试和step convergence确认；可暂保守固定正确trapezoid并验证精度，勿在同一回归中悄悄大改switch再把差异都归轴子修复。
3. **early return pitch对齐** 是继承，但加入轴子后不该同时跳过纯轴子旋转；需要新链分离。
4. `tau+=aI dl scale`、`tauF+=|rV|dl scale` 在Ipol阈值分支内部；tauF没有统计新delta旋转，只是绝对plasma Faraday depth，不是净EVPA或轴子效应。不能用tauF≈0证明所有旋转都关闭/很小。
5. pol_frac>1时除以(pol_frac+.005)及j偏振renorm都是继承的非物理修正保护，不是我们故意屏蔽轴子；应记录触发次数而非直接删。full在纯旋转应保持I/V/|P|，若仅轴子导致clip频发需查basis/solver。
6. `POLARIZATION_ACTIVE` 以 Ipol>1e−100 而非I决定。最终IQUV也只有flag真才写，纯I而P=V=0的测试可被输出为0；是继承限制，不应把“pure unpolarized source”用作未经修正的total-I回归oracle。
7. `get_fluid_params` 会对传入X做角wrap及X2−1e−6；X_u_next直接取未处理lightpath，而当前X经model处理，两端坐标处理非完全一致。该细节与full的场求值独立化一起处理，避免求值函数改变ray坐标；不在此重审fluid模型。

## 5. 纯screen近似、Fig4与独立test provenance

### 5.1 screen和full不是同一东西

物理上jonly并不等于给最终像素统一旋转：多个发射点的P会按各自轴子场/相位加权相干叠加。即使a=r=0，也有disk thickness、LP、多像素beam averaging导致washout或增强。只有已知单一发射位置或适当dominant-source近似时才有每像素 `P'=P exp(2iΔχ_emit)` screen。

`EHT_ZC/GRMHD/main.ipynb` cell36 `get_analytical_evpa` 正是后处理screen：取baseline EVPA和|P|，每相位加 `axion_field_coord[i−1]`，重构Q/U，不沿ray重新吸收/Faraday。cell43读四种 common/fake/eq_jonly/jonly产品；cell53把eq_jonly r9标(1)、r4标(1)+(2)，jonly r9标(1)+(3)、r4标(1)+(2)+(3)，cell54/56/60/64比较jonly/common。这是明确的科学诊断设计，**不是看到jonly就全删**。

`get_analytical_evpa`本身只从baseline取|P|、设Q/U，V=0；这是linear-only screen构造，不是完整IQUV变换的通用oracle。相关低层EVPA测试还存在 `arctan2(Q,U)` 与主 `arctan2(U,Q)` 的不同写法（cell77/78等），镜像会保振幅但改变相位/handedness，应作为Fig4验证单列，别把诊断图中相位差都归core。

### 5.2 文章主张应如何理解

已读两份完整草稿，不把 `%`注释与作者宏混同。`main.tex` 的 `\yx/\yz/\zc` 定义是彩色输出，**不是LaTeX注释掉的文字**。诸如“Change to the correct number”“Check this number”仍是可见待核实意见。

- `main.tex:309` 写plasma差异≲3%并带作者“改正确数字”；`367` 又写大多时刻至多5%；都未通过本轮数据重新验证，不能当作已证实误差预算。
- `355–367`明确full包含j/a/Faraday，薄赤道只用于隔离第(3)组；平滑用于第(4)组。这支持保留具名诊断开关。
- `old.tex_copy:189`四列顺序为toy解析→GR薄赤道→厚度screening→包含全部plasma的GRMHD，最后一列full并不是建议把全套关掉。
- 图序在两个草稿里不一致：main.tex字面“Fig.4”在408附近是GRMHD uncertainty占位；用户关心common/jonly/eq_jonly/fake实际对应`fig: equator1dplot`及Notebook这组图。报告以标签/产品名而不是易变图号定位。
- 旧稿后面“relatively optically thick”“wash-out使所有信号更小”是过度叙述；新稿说明矢量相消可放大信号。**光学厚度（吸收）和几何厚度（发射位置分散）不是同一概念**。Faraday仅Q/U同轴旋转时与轴子旋转对易，本身不自动破坏单光子端点总导数；conversion/二向色吸收、源项叠加才是关键。不能只写“有Faraday所以topological都失效”。

推荐模式定义：

| 模式 | 发射 | 吸收/Faraday | 轴子 | GR | 科学含义 |
|---|---|---|---|---|---|
| full/common | 实际GRMHD全体有效cell的thermal j | 全开 | 全ray场耦合 | 全开 | 主fiducial |
| jonly | 同一GRMHD j | aI/aQ/aU/aV/rQ/rU/rV全部0 | 全ray | 全开 | 隔离传播plasma影响，仍有厚度/LP |
| eq_jonly | 同一GRMHD j但按明确定义赤道薄层筛选/积分 | 全0 | 全ray，不随发射mask关闭 | 全开 | 厚度诊断；保留所有赤道穿越不等于只保留首次direct |
| fake | 常数/高斯ring等显式toy source，可搭配独立薄层mask | 通常全0，明确记录 | 真实场或toy profile另选 | 全开或显式flat | 总导数/几何/插值回归，不能当fiducial |
| screen | baseline图像；单一发射映射 | 无联合再传播 | 每像素一次后处理角 | mapping里已含GR | 快速近似/对比，不替代full |

如果为了jonly关闭a/r，仍应让轴子在非发射/真空段传播。eq_jonly只屏蔽j即可，不能把轴子一起按薄层gate掉。若要关闭相应区域plasma可同时清a/r，但这也是单独物理选择。

### 5.3 当前本地独立test到底能证明什么

两个run目录 **没有独立 `.c` 偏振源码**。makefile均用 `$(RAPTOR)/src`，不写死core路径；Run.py也不编译、不设置RAPTOR。今天源码/目录名不足以确定历史binary或产品来源。

- `raptor_run_test/Run.py`只选快照切片[100:101]、alpha.3、100²、output_test；没有Python层jonly切换。
- `raptor_run_test_eq_fakevalue/Run.py`只选[151:152]、alpha.3、400²、output_test_eq_fakevalues；同样没有Python层fake emissivity选项。
- 两者当前 `model.c` 与主run字节相同；test header与主相同；fakevalue header唯一相关差异是SIGMA_CUT1→1000。因此fakevalue不能无条件与主fiducial对照：高磁化辐射区开关已不同，虽完整模型影响留父审查。
- 主run偏振对象 SHA256 `bbc254498d649fe728e6f091694181df9226b993beeee2c3fc2fb45d5e038da6`；fakevalue当前偏振对象 **完全相同**；test对象 SHA256 `22c32d9fd468a25a106e3845080497c63a66b95bfeae2da85ff6092bd5a16694`。
- test对象的RK4和trapezoid函数机器码section与主对象相同，差异在pol_integration_step（大小0x101d vs主0xf0b）。只读反汇编 `test-step-asm.txt` 偏移0x854..0x910明确可见 `ri*cos(thi)`、浮点转int(`cvttsd2si`)、非零则11个系数存零。实际筛选相当于int(z)≠0（有效厚度约|z|<1Rg），不是注释期望的 .01/.1Rg；这是调用C整数abs后与小于1阈值比较的典型语义。**可以确认此对象有active高度筛选，但不能仅机器码恢复作者原先到底写了>.01还是>.1，更不能确认四种命名cube分别对应哪个对象**。
- 当前test对象没有常数替换源的证据，清零分支以外保留原coeff；不能因test名字把它叫fake。也没有证据它全域关闭a/r，所以不能直接叫标准eq_jonly。
- fakevalue当前偏振对象相同不表示其整个binary相同：header、model.o及binary不同。更不能证明旧文件夹下历史产物从来没用过fake source。
- EHT同层 `test_raptor_difference/src` 当前是带DM_mass/DM_cross_section/channel/pikachu和暗物质jI插入的另一个模型，HEAD79d7605且工作区有修改；有active GRMHD_cut输出。它不是这些轴子四模式的可信同源源码，不应拿来填补缺失test源码。EHT同层旧`raptor`停在d11f522，沿用早期数组/读法，同样不是当前b6主core来源。

`raptor_run_test_eq_fakevalue/GRMHD_cut.ipynb` cell6用forward存档轨迹中最后θ<π/2的穿越附近插值得赤道映射；cell7使用 **0.29915625** 而不是主core .3生成screen相位：这是alpha.3的下游omega4值，和producer本征频率 .296319073也不同。该事实直接说明screen/full对比尚未统一时间频率，不能把所有比值差归plasma。

同cell用 `last_idx_1/last_idx_2`，需验证“最后”确实是从observer看的第一次穿越，跨φ=0包裹用相位而不是直接φ平均、两点真夹根/分母非0；本轮未重算mapping，不声称具体像素已经错误。

重跑前每个mode需固定一个版本化源/配置+对象/编译器+输入hash，产物metadata写mode。不得通过注释编辑共享core、复制旧.o、不清build来制造common/jonly差异。

## 6. full fiducial 开关表（重跑合同）

| 物理/实现项 | 当前主core/run实际状态 | 建议主full fiducial | 其它模式/备注 |
|---|---|---|---|
| POL | 1 | 1 | POL0没有轴子通道 |
| 发射jI,jQ,jV | evaluate_coeffs_single，DF=TH；jU=0 | 实际GRMHD thermal全开 | jU=0是B-aligned tetrad对称性，不是把天空U关闭 |
| 吸收aI,aQ,aV | 全开；aU=0 | 全开 | jonly/eq_jonly显式全清a |
| Faraday rV,rQ | 全开；rU=0 | 全开 | rQ是conversion，jonly全清r |
| EMISUSER | 无调用选择，当前宏0 | 0并记录 | 若启用须实现选择，非单改宏 |
| fluid质量/cut | ||绕过有效性返回 | 有效fluid才作plasma | 无效区轴子仍独立传播 |
| 轴子norm | baseline0及response.1 | baseline0+小信号.1，并做norm线性回归 | 记录norm究竟是bc还是c；禁止重复乘b |
| 轴子质量alpha | 叫omega，.1..3查询NormR | 独立alpha，拒绝域外 | 高alpha表动态尺寸及spin可追溯 |
| 振荡频率ω | 直接用alpha | 独立Re(Mω)或明确ω≈α的legacy | full/screen/天单位卷积一致 |
| 云形状 | R(alpha,r) sinθ，m1，实模长 | 显式211包络近似，说明误差 | 完整复模式另版本 |
| 云坐标 | BHAC KS t/phi；θ覆盖 | 与BL草稿一致或明确KS toy | hslope0当前θ无影响 |
| 云25Rg窗 | active W25,15 | 默认关；全ray/真空补齐 | 旧窗只legacy或有收敛证据的显式近似 |
| delta rotation | plasma后一次、backward差分 | 修端点/符号、受控split或联合解 | legacy单独记录 |
| 轴子rV注释方案 | off | off（若用split） | 联合solver只能一处计入 |
| early time cut −10023 | off | off | 具名diagnostic |
| 赤道mask | off | off | eq_jonly明确fabs/厚度/收敛 |
| 常数/高斯fake发射 | off | off | fake独立toy source |
| 均匀toy axion profile | off | off | 回归screen测试 |
| 轴子阈值5e−8 | off | off | 不引入小场不连续性 |
| pixel编码输出 | off | off | 只调试布局 |
| GR geodesic/平行输运 | RK4；max_order100 | 保留；步长/高阶阶数收敛 | 不能声称只含first equatorial |
| STEPSIZE | .005 | 经.005/.0025/.00125合成→小样本收敛再定 | 当前不能认证已足够 |
| 相位采样 | phase0..15，π/8 | 保留、metadata radians/index | 与baseline file index区分 |
| t_ref | 硬编码1022 | 显式可追溯参数 | 全局phase nuisance可吸收，不需盲删 |
| screen图像旋转 | 不在主core；Notebook另做 | 不替代full | Fig4比较同ω同坐标同normalization |
| SIGMA_CUT | 主run1 | 固定父审查确定的fiducial | fakevalue1000只独立diagnostic，不混主 |
| build/输入 | VPATH环境、旧对象可能复用 | 新干净build，固定RAPTOR路径和全部hash | 不运行现存oldbinary |

此表是建议的重跑合同，不表示本轮已实施任何修复。

## 7. Must-fix 与验证顺序

### P0：在任何真实full cube重跑前必须解决

1. **ray端点/存储合同**：禁读lightpath[steps]未初始化记录；定义forward段p→p−1、delta与Stokes旋转符号成套；补steps0/1/max_steps边界与纯真空总导数单元测试。
2. **plasma与axion解耦**：移除 `||` 对无效fluid的强行plasma调用；真空、sigma-cut、pitch-aligned区域仍处理轴子；不能只改回&&。
3. **NormR初始化、错误和域边界**：并行区外一次安全加载；验证行列/grid/spin/hash；质量超界明确报错；norm0跳过轴子；高alpha表不能直接换同名文件。
4. **明确物理模型**：alpha≠ω拆参数；BL/KS及R模长/复相位约定定案；full与screen/相位拟合/天单位time averaging用同ω。若选择ω≈α旧模型，须明示，不能一边producer用本征ω一边宣称下游精确一致。
5. **25Rg窗处理**：决定是legacy数值近似还是物理fiducial，测试窗口/真空补齐收敛。未验证不能无说明继续或无说明删除。
6. **run输入和产物完整性**：禁止旧无轴子model.in；读入逐项验证。Run.py当前用os.system不检查退出码、见旧文件即跳过、末尾创建done标记不能证明成功，需新目录+返回码+输出shape/finite/metadata/17相位完整性检查，不能污染/覆盖旧结果。
7. **独立诊断模式可追溯**：不靠目录名/注释/旧.o，输出mode与开关；test中int abs高度筛选不能作为薄层定量比较。生产清理构建放新工作区，不改原证据目录。

### P1：科学误差预算/主full发布前必须验证

- 验证stiffness（继承）/单位/complex特征值问题；统计显隐式调用率及clip次数，常系数j-only/absorption/Faraday rotation/conversion解析回归。先把轴子修复和上游数值solver变动分开对照。
- 证明split收敛，不能用RK4名称或norm很小替代；精确旋转并不保证有source/不对易plasma时整体高阶。
- 单点源、多点源、纯I/纯P、云外真空传播、轴子norm0、倍增norm、交换正反传播、镜像U符号这些合成测试通过后才做小像素/小snapshot回归。
- common/jonly/eq_jonly/fake统一same snapshot、resolution、metric、sigma cut、table、alpha/omega、t_ref、mask、beam；先差分比较再统计文章3%/5%/30%主张，低P/相消像素相位应mask/不确定度加权。
- 物理norm定义与b(alpha)后处理一次计入；M87质量/观测时间换算统一；完整复本征场或实包络近似标清楚。
- 图像mirror/IAU basis和EVPA符号按旧稿末尾“RAPTOR y-flipped、−iU”备注明确；相位差π不等于振幅错，但不能全靠它解释径向KS/BL差或频率差。

### 不应“必须删除”的东西

- 轴子delta指数旋转（科学核心）；合理的r上界防护；本征径向查表；主发射/吸收/Faraday；诊断jonly/eq_jonly/fake能力；GR parallel transport；继承的finite/physical guards。
- 1022偏移不是已经证实的数值错误，需显式来源/统一相位；注释输出可保留开关；窗口需科学评估，不是按debug标签一刀切。

### 本轮证据边界

**确定**：基线关系；四处源码/表身份；主关闭块均注释；active ||、错位端点、表race/未检查；alpha身份与timealpha；25Rg窗；BL/KS实现不一致；test对象active高度筛选；当前fakevalue偏振对象等于主；synthetic UB/rotation/splitting结果。

**没有证明**：每份旧cube构建时的源/开关；真实cube受影响百分比；作者1022/25Rg窗唯一原意；当前源码.4实际被生产执行；所有alpha径向表都坏；文章plasma差异一定超3%或5%；修复后约束必然变强/变弱。上述不得写成已证实结论。

## 8. 文献和本轮probe

先网络检索 `https://export.arxiv.org/api/query?search_query=all:RAPTOR%20AND%20all:polarized&max_results=5`，返回RAPTOR II **arXiv:2007.03045, A&A641 A126, DOI10.1051/0004-6361/202038573**；检索XML已保存。随后下载PDF、pdftotext，查Sec3.1–3.3和AppendixA/B：标准vacuum/plasma分开、向相机forward transfer、B-aligned tetrad jU/aU/rU=0、RK4/implicit trapezoid；没有把“implicit Euler”源码注释当算法事实。

缓存 `2105.04572.txt:118–154`（此前工作）Eq1/2也明确真空端点公式、BL场、bc/(2π)归一化和需联合plasma的限定。草稿二份已读，其作者注释/旧版本主张作为意图证据，不作为当前active开关证据。

主要证据文件：

- `identity.json`、`axion-core-base-fork.diff`、`axion-core-main-fork.diff`、`upstream-only-pol.diff`
- `pol_rte_integrator.c.numbered.txt`、`core.c.numbered.txt`、`functions.h.numbered.txt`、`gr_integrator.c.numbered.txt`、`GRmath.c.numbered.txt`
- `integrator_probe.c/.log`：从当前solver抽取，纯旋转、stiffness和source split；未调用模型/真实ray。
- `lookup_probe.c/.log`：synthetic23×2400表，.31正常，.4界外UBSan报错。
- `generated-model.in`、`parser_probe.c/.log`：只构造wrapper输出文本验证扫描，未执行Run函数/旧binary。
- `main-step-asm.txt`、`test-step-asm.txt`、`test-step-asm.diff`、复制的pol.o及section bins：只读静态对象检查。原件未改。
- `axion-fig4-notebook.txt`、`test-equatorial-screen-cells.txt`、`test-fakevalue-main-cells.txt`。

补充说明：首次objcopy尝试直接dump挂载对象时工具尝试创建临时文件，被只读文件系统拒绝；没有改动原件。后续所有objcopy操作均对/tmp字节副本进行。没有执行任何存量对象/可执行文件。
