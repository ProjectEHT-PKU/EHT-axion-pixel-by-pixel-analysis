# 两份 Radial_function 上传包核查

## 范围

用户上传`/home/canoziia/MyProject/PHYS/Radial_function.1.tar`与`.2.tar`。本轮读取tar inventory，手动复制小型regular source/table到`/tmp/eht-radial-audit`；不跟随tar链接、不覆盖原目录，不执行上传binary。主助手另在临时目录编译了只读C探针，只调用频率/单点径向函数，未运行原main，也未重新生成生产表或RAPTOR图像。

未改pipeline算法、配置、原EHT_ZC或运行结果。另一个GPT-6.1 Sol只读独立审查了两版数学与边界；实验与hash比较由主助手完成。

原包SHA256：

```text
.1: ad261bfab7acebd110179a781323acfaea22a7c94d2f7ad0bd105273bd3795ca
.2: f28aa84fd894540b0c2d4a4430b95aed6bc0b31106789195f0ac1f4d1f511a3c
```

本轮挂载目录`vmdk/`已为空，因此没有重新访问之前挂载的producer文件。不过原EHT_ZC根目录两张表均仍在：后续补查已确认NormR.csv与.9375输出、NormR_0.99.csv与.99输出分别完整SHA256相同。初查遗漏了本地NormR_0.99.csv，不能据挂载为空断言无法做高alpha表的精确核验。

## 1. 哪份更新，差在哪里

- `.1.tar`较新：Axion_Radial.c mtime2025-09-09、main.c2025-09-17。
- `.2.tar`较旧：这两个文件均2025-07-22。
- `.1`大部分体积是打包了`.venv`；科学文件约16.1MB。`.2`约13.3MB，不应以435MBvs13MB判断算法复杂性。
- 两版共享Axion_Radial.h、README、三个Mathematica Notebook及source频率表；核心变化集中在输出网格、归一化和main默认参数。

| 设置 | .2旧版 | .1新版 |
|---|---|---|
| main默认自旋 | .9375 | .99 |
| 默认质量参数区间 | .1–.316 | .30–.44 |
| 质量网格 | 40点linspace | .01步长 |
| 半径网格 | r=1,2,…,1200 | r=1…240.9，dr=.1 |
| 输出归一化 | 注释，未执行 | 按max归一化，但顺序有问题 |
| 径向输出精度 | 6位小数 | 17位小数 |
| 输出.a=.99 | main写result/50.txt | 写result/result_a=.99_l1m1.txt |

旧版`(600-2)/598`整数运算恰等于1，实际半径为2..1200，不是到600；不能按表达式注释猜网格。

## 2. 现有NormR对应哪份

**.1包中 `result/result_a=0.9375_l1m1.txt` 与本机 `EHT_ZC/NormR.csv` 字节完全一致：**

```text
SHA256 0a352558f8b846f2481ab2acf35c2d6f8b30ce1f83f2e6f5196660f23433c353
55200行 = 23个质量参数(.09–.31) × 2400个r采样
```

这证明当前本机表在新版包中有直接成品来源，不只是网格相似。新版main当前设.99/.30–.44，不是生成这份.9375表时的参数；成品存在不能证明当前main配置曾经未变。

两包顶层都保留了相同的2021 `NormR.csv`：17080行、35参数、每参数488个r点，**不是当前RAPTOR的55200行表**。因此不能拿包里顶层同名NormR直接替换生产表。

.2的result.a=.9375是48000行旧粗网格，且与其result_old里的2021同spin CSV字节一致。它不是当前2400点表。

新版.a=.99表：36000行、.30–.44、2400半径点，与本机`EHT_ZC/NormR_0.99.csv`全字节一致，SHA256为`0ff8b6c816c997f32ef4d2029f52fd6cabd22990dab85362a54e35cbc62f3815`。

### 已找到的切换调用，而非要求用户回忆

`EHT_ZC/Plot/fig1.ipynb`（cell零基）：

- cell2读取`../NormR.csv`为R_data；cell3读取`../NormR_0.99.csv`为R_data_99。
- cell4的plot_R_r_curves显式selected_alphas=[.1,.15,.2,.25,.3]，selected_alphas_dashed=[.35,.4,.44]。
- cell5调用plot_R_r_curves(R_data,R_data_99)。**八个研究节点的分配已明确：.30用.9375表，高三点用.99表。**
- cell7的连续径向图调用R_data[1:-1]和R_data_99[:-1]：实际分别画.10–.30、.30–.43，在.30有重叠，不是严格无重叠的一张拼接CSV。不能把这段绘图解释为producer全范围唯一切换阈值。

主限制计算`EHT_ZC/result.ipynb`cell33读取八alpha的EVPA_cos_fit_result；cell36/50分别统计2d/1d，调用链读取的是RAPTOR已生成响应，不直接用NormR再生成高alpha响应。`GRMHD/alpha_0.44.py:78–80`读取同名.44原cube。

主模拟节点最大为.44，未发现.45独立cube。result.ipynb的.45多为轴刻度；cell105有旧五成像方法参考段`linspace(.1,.45,100)`和calculate_ul95_with_data调用，非主八点GRMHD链，且当前Notebook源码未发现该函数定义。make_interpolator(cell23)对超出.44的请求夹到末节点。因此“图到.45”不能证明独立.45RAPTOR计算。

## 3. alpha和omega：现在已明确，不应再含糊

两版生成器：

```c
// Function_R
double M = 1.0;
double omega = get_omega(mu, a);
get_Intermediate_variables(M, a, omega, mu);

// print_cable
mu = ...;
Result[i][j][0] = mu;
```

`get_omega`读取source频率表，将第二列插值后乘mu，区分mu²和omega²构造q。

**所以CSV第一列为M=1单位的质量mu，即alpha=Mmu；不是频率omega。** source第二列表示omega/mu；get_omega返回无量纲Momega。

这纠正父助手上轮“实际扫描是omega”的过度结论：RAPTOR变量虽然叫axion_omega，但它用来查的这张表实际按alpha生成。将同一个值直接放进cos(time×axion_omega)，属于omega≈alpha近似或实现不一致，不能靠把全部标签重新解释为omega来消除。

### 频率例子（由文件表按C同样线性规则计算，避开非法边界）

| spin | alpha | 表给Re(Momega) | 下游omega4 | 211领先阶omega2 |
|---:|---:|---:|---:|---:|
| .9375 | .15 | .1495522968 | .1498945313 | .1495781250 |
| .9375 | .30 | .2963190733 | .2991562500 | .2966250000 |
| .99 | .30 | .2962465998 | .2991562500 | .2966250000 |
| .99 | .40 | .3904155000 | .3980000000 | .3920000000 |
| .99 | .44 | .4265520369 | .4373380000 | .4293520000 |

这些表近似更接近主量子数2的领先式，而不是当前n4。source文件是插值/提取频率，不是本轮精确本征值求解；不能把表本身视为无误oracle。

## 4. 确认的问题及实验

### 4.1 两版共同：频率插值下边界越界

`Axion_Radial.c:87–94,155–162`：找到首个≥mu的点后无条件索引减1。

对于a=.9375、alpha=.1，imag表首个点恰=.1，i变−1。主助手编译只读probe加`-fsanitize=undefined,bounds`，得到：

```text
Axion_Radial.c:94:64: runtime error: index -1 out of bounds for type 'double [500][2]'
```

低于imag表最小.1的.09同样非法；新版现存.9375结果含.09。即使虚部后续被丢弃，越界仍是C未定义行为，不能说无害。正常区间内部插值公式本身正确。

### 4.2 两版共同：复杂频率虚部被丢弃

`Function_R:223`使用`double omega=get_omega(...)`，返回的imag被截掉。编译器也明确报告complex→double discards imaginary component。

这是确定执行行为，可能意图采用微小增长率下的实频率近似；但不能称完整复本征函数。print_cable另算的omega_complex没有参与径向计算。

### 4.3 新版：max先计入视界内值，再清零视界内点

`Axion_Radial.c:248–268`先对r=1..240.9计算Function_R和max，再把r<r+设0，最后除先前max。外部径向表达式也被错误用于归一化分母的视界内点。

主助手在同源临时副本只读计算a=.99的两个2400点slice：

| alpha | 视界内max所在r | 最终外部max | r5当前表 | 若仅按外部max归一化的r5 |
|---:|---:|---:|---:|---:|
| .30 | 1.1 | .68697765795 | .36662188472 | .53367366533 |
| .31 | 1.1 | .84272392620 | .47895819407 | .56834531355 |

前两列最终max和r5与上传.a=.99成品精确吻合。因此不是只凭静态代码猜测，**该错误确实进入这些已存slice**。不能据此声称所有alpha都受同程度影响；.32起输出max1，与上轮摘要一致。主链低alpha采用.a=.9375表时，不能将这里.a=.99低端问题机械套过去。

### 4.4 新版：round生成末节点可超出Mu_max

例如print_all_cable设.1–.316，实际.01网格末点.32，超出imag表.316；get_omega返回0却不终止输出。新版main当前.30–.44没有此问题。

### 4.5 其他实现防护缺口

fopen/fscanf、行数500上限、零/非有限归一化分母没有完整验证。大局部数组可能遇栈容量限制；不宣称当前必然崩溃。旧版被注释get_max只扫200点不属于活跃路径错误。

## 5. 数学近似与质量边界

对照已缓存Dolan arXiv0705.2880 Eq33–45：Leaver径向前因子、q负支和三项递推未发现明显抄错。仍有：

- angular eigenvalue固定lambda=2，是l=m=1球谐最低阶近似，不是完整spheroidal特征值；
- forward recurrence只截断1000项，没有收敛/远场稳定性诊断；
- 没在生成器里解连分式本征条件，而是使用外部source频率表；
- README明确某些spin数据由其它spin插值；l=m=2实频率借用l=m=1，只改flag不够；
- 活跃a=.99 imag表乘alpha后的增长率与Dolan论文TableI最大值不一致，不能直接称精确复现；该imag目前被Function_R丢弃，因此未直接进入径向形状。

这些是近似/来源可信度问题，不应与确定越界/归一化bug混成一个“新版不对、旧版就对”的结论。

## 6. alpha=.1边界修复后的小规模重算

按用户要求，在临时副本只将两个插值括号的无条件`--i/--j`改为`if(index>0)--index`。
保留原实频率计算、1000项递推、2400半径点、原归一化顺序，不实施其它算法修正。
编译使用gcc C99/O2及undefined,bounds sanitizer，计算a=.9375、alpha=.1单切片，约2秒完成且无sanitizer错误。

与原EHT_ZC/NormR.csv的alpha=.1切片对比：

- 最大绝对差：**8.659739592076221e-15**；最大非零点相对差：**8.949642790722247e-15**。
- 两者全部有限、峰值均1、峰值半径均197.6Rg。
- 并非字节相等，2374个值存在浮点末位差；跨编译/数学库的重算不能把这些差别唯一归因于边界修复。
- 原完整CSV SHA256重算前后相同，未写入原文件。

**结论：这项边界修复没有造成有意义的径向形状变化；不需要为它替换当前.1切片或重跑RAPTOR。**
这不证明所有频率/递推近似正确，只解决了先前“0.1越界是否已污染当前径向值”的疑虑。
证据保存在`radial-evidence/alpha-0.10-check/`，生成切片仅保留在/tmp，未作为新生产表发布。

## 7. 结论与下一步

- **哪份新？.1。哪份包含当前本机NormR成品？也是.1的result/，不是顶层同名CSV。**
- **哪份整体正确？不能二选一认定。**新版输出格式匹配当前RAPTOR，但边界和归一化有确定问题；旧版保留同样的核心频率边界问题、未归一化且网格不匹配。
- alpha/omega的物理身份已由生成器确认。应分离alpha查表和omega时间推进，不反解alpha来掩盖producer混用。
- 建议保留两个上传包不动，以新版作来源基线，修边界/视界外归一化/数据有效性并补收敛后，先对照已有径向表，评估哪些RAPTOR产品实际需要重算。不能只改下游频率或替换一个CSV就声称旧模拟已修正。

本轮未实施上述修复。相关C探针、日志、source snippets和JSON证据随报告保存；源包本身保持不变。
