# M87* 常数与分析参数核对表

**图像参数专项更新：**见[视场与Gaussian来源核查](IMAGE_GEOMETRY_PROVENANCE.md)。
现已找到历史HDF5的46.2261Rg坐标、三个参数与128μas/20μas/2.355的精确比例关系，
以及历史raw输出不同M/D角尺度的证据。该专项报告优先于下文首轮“仅待查”的状态；
实际主cube与EHT pickle的统一标定仍未确认，参数未改。

## 范围与原则

本表用于逐项验证此前代码中已算好的常数，不因其是旧值就默认正确，也不因论文出现不同值就直接替换。

代码基线：`74ebae1`。本轮仅查资料、提取与重新计算；**未修改科学代码、未重算数据、未改变共享manifest**。

核对对象分为：

1. **观测量/文献采用值**：质量、距离、喷流方向等，有不确定度和适用数据集。
2. **导出换算常数**：引力时间、角尺度、α到质量、Gaussian FWHM到σ。
3. **模拟配置/分析选择**：自旋、倾角、FOV、blur、曝光近似、mask、先验等。文献有类似设置不证明本次RAPTOR配置相同。
4. **拟合产物**：GRMHD日差sigma等，应回到输入与拟合方法验证，不由EHT论文提供标准答案。

这是**首轮证据与待办清单**，不是“全部常数已经查完”。随上游各阶段整理继续补充。

## 1. 已下载核对的原始资料

论文PDF下载到`/tmp/eht-constant-sources/`，用`pdftotext -layout`提取；PDF大文件不加入仓库。下列页码为实际PDF页码。本MD保存稳定URL、章节/页码及核对依据，/tmp缓存不作为长期依赖。

| 编号 | 一手资料 | 本次核对位置 |
|---|---|---|
| EHT-I | [First M87 EHT Results I: The Shadow of the Supermassive Black Hole](https://arxiv.org/abs/1906.11238)，[PDF](https://arxiv.org/pdf/1906.11238)，[DOI](https://doi.org/10.3847/2041-8213/ab0ec7) | PDF p1摘要质量；p2距离；p8 Table 1和质量推导 |
| EHT-IV | [First M87 EHT Results IV: Imaging the Central Supermassive Black Hole](https://arxiv.org/abs/1906.11241)，[PDF](https://arxiv.org/pdf/1906.11241)，[DOI](https://doi.org/10.3847/2041-8213/ab0e85) | p13–14 §6.2 eht-imaging的64²/128μas；p8–9 beam；p19 compact flux |
| EHT-V | [First M87 EHT Results V: Physical Origin of the Asymmetric Ring](https://arxiv.org/abs/1906.11242)，[PDF](https://arxiv.org/pdf/1906.11242)，[DOI](https://doi.org/10.3847/2041-8213/ab0f43) | p3 §2的PA≈288°与17°；§3模型库质量/距离/角尺度 |
| EHT-VII | [First M87 EHT Results VII: Polarization of the Ring](https://arxiv.org/abs/2105.01169)，[PDF](https://arxiv.org/pdf/2105.01169)，[DOI](https://doi.org/10.3847/2041-8213/abe71d) | p10 §5.1的10°窗口、1°步长、强度加权Eq.(10–11)、19/20/23μas；p22 §D.1的120μas/64²与0.6Jy |
| EHT-VIII | [First M87 EHT Results VIII: Magnetic Field Structure near the Event Horizon](https://arxiv.org/abs/2105.01173)，[PDF](https://arxiv.org/pdf/2105.01173)，[DOI](https://doi.org/10.3847/2041-8213/abe4de) | p8模型质量/距离/flux/倾角；p9 §4.4的20μas FWHM；§5自旋网格 |
| IAU-B3 | [IAU 2015 Resolution B3](https://arxiv.org/abs/1510.07674)，[PDF](https://arxiv.org/pdf/1510.07674) | nominal solar mass parameter `(GM)_sun^N = 1.3271244e20 m³/s²` |
| CODATA | [NIST Fundamental Physical Constants, ASCII table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt) | `c=299792458 m/s`、`hbar=6.582119569...e-16 eV s` |

后续还需查EHT-II/III的具体观测时间/频段、原始RAPTOR图像元数据及实际模拟配置。6小时矩形平均目前是方法近似，不能因有四个观测日就视为精确观测window。

## 2. 文献采用值并不唯一：禁止混用

| 用途 | M | D | 其他 |
|---|---|---|---|
| EHT-I的M87*质量结果/距离采用值 | `(6.5±0.7)×10^9 M_sun` | `16.8±0.8 Mpc` | angular gravitational radius `θ_g=3.8±0.4 μas`；M与D/θ_g并非三个独立测量 |
| EHT-V/VIII模拟库的标定 | `6.2×10^9 M_sun` | `16.9 Mpc` | `θ_g≈3.62 μas`；是库的采用值，不是另一次独立质量测量 |
| 当前本文草稿`main.tex:192–193` | `6.5×10^9 M_sun` | 需查原模拟配置/成像元数据明确 | `a*=0.9375`、`i=163°`为本文基准模拟选择 |

不能从“EHT论文有6.2和16.9”就自动认定当前RAPTOR数组必须采用6.2；也不能把质量不确定度传播时机械当成与距离独立。

## 3. 统一重算公式与结果

本次采用IAU nominal solar GM计算，避免把G与太阳质量的各自近似重复引入；取parsec为`3.085677581491367e16 m`。

```text
r_g = (GM)_sun^N × (M/M_sun) / c²
t_g = r_g/c
θ_g = (r_g/D) × (180/π) × 3600 × 10^6        [μas]
omega_hour = omega_dimensionless × 3600/t_g  [hour^-1]
m_a[eV] = α × hbar[eV s]/t_g[s]
σ_angle = FWHM_angle / sqrt(8 ln 2)
σ_rg = σ_angle / θ_g
FOV_rg = FOV_angle / θ_g
pixel_rg = FOV_rg/Npix
```

| M / D采用值 | t_g (hour) | θ_g (μas) | `3600/t_g` (hour^-1) | `m_a/α` (eV) | σ对应20μas FWHM (r_g) |
|---|---:|---:|---:|---:|---:|
| 6.5e9 / 16.8Mpc | 8.89324754 | 3.81899325 | 0.112444863 | 2.05590426e-20 | 2.22394161 |
| 6.5e9 / 16.9Mpc | 8.89324754 | 3.79639566 | 0.112444863 | 2.05590426e-20 | 2.23717936 |
| 6.2e9 / 16.9Mpc | 8.48278997 | 3.62117740 | 0.117885743 | 2.15538349e-20 | 2.34542997 |

派生量可以为计算保留精度，但物理精度不能高于质量/距离等输入精度。

## 4. 代码常数逐项核对

### 4.1 质量、时间与轴子质量换算

| 代码值 | 位置 | 核对结果/状态 |
|---|---|---|
| `0.112411` | `eht_axion/model.py:13`；旧`result.ipynb`轴子日差模板 | 若M=6.5e9，应为0.112444863 hour^-1；旧值偏低约0.0301%，相当于t_g=8.89593h。很接近，可能是旧常数取整；精确计算来源待确认，暂不改 |
| `2.05591e-20` | `result.ipynb` cell105等，α→轴子质量轴 | 与M=6.5e9的`hbar/t_g=2.05590426e-20eV`吻合，差约0.00028%；它是rest-energy质量单位下的换算，不是角频率或Hz |
| `8h` | `result.ipynb::get_EVPA_diff_xday`注释 | t_g取整近似；准确为8.893h，不能视为精确值 |
| 每标签`0.1r_g`，每10标签保留一帧 | 草稿`old_images/old.tex_copy:741`；`Introduce.md:2` | 应规范为时间`t_g`；相邻NumPy帧间隔1t_g≈8.893h，141帧首末140t_g≈51.877天 |
| `step=30`、`r_step=3` | `result.ipynb` cell10 | 背景日差用了3帧，真实间隔≈26.68h；比24h长11.17%。近似有历史说明，不是无法考证，但需后续决定怎样统一 |
| `6 / 24 / 5` hours/days | `model.py::evpa_templates`；旧模板 | 与草稿6小时平均、一天日差、两组中心五天间隔一致；真实observing-window待EHT-II/III与visibility时刻验证。5天相位推进额外`/2`是公式问题，不是物理常数问题 |

### 4.2 角尺度、FOV与Gaussian blur

| 代码值 | 位置 | 核对结果/状态 |
|---|---|---|
| `final_fov=32.6868` | `ImageS`、alpha脚本、`polarization.py` | 单位r_g。按6.5e9/16.8对应124.83067μas；不恰为120或128μas。若原成像确是128μas，该值隐含θ_g=3.91595384μas；若120μas，则隐含3.67120673μas。**需查本次pickle元数据/原重建设置**，不能仅凭论文替换 |
| `init_fov=46.2262` | 同上 | 与`32.6868×sqrt(2)`几乎相等（约1.8ppm差），支持为旋转裁切的几何选择；不是M87*观测常数。其实际RAPTOR角FOV仍需配置/metadata核验 |
| `0.51073125` | EHT mask半径、`result.ipynb`等 | 精确等于`32.6868/64`；内部一致，但依赖最终FOV的外部标定正确 |
| `gaussian_kernel=2.16871` | `utils/libanaly.py:378,411–417` | 传给`scipy.ndimage.gaussian_filter`前乘pixel/r_g，因此该值是**σ[r_g]**，不是FWHM。FWHM=5.10692178r_g，按6.5e9/16.8对应约19.50μas；与20μas差约2.5%，需核对当时θ_g与期望beam |
| 注释`20uas=2.16871` | `utils/libanaly.py:411` | 必须补足“20μas是FWHM、2.16871是σ及使用的θ_g”；否则单位不同而误读 |
| `0.518838` | `ImageS`的`blur=False`分支 | 与主分支pixel_rg=0.51073125不同；当前默认主路径不走，后续整理非blur对照时核对 |
| `16.343/16.3434` | 一维积分/旧坐标 | 来自最终FOV约一半；积分离散终点还有独立quadrature问题，不能当独立物理常数 |

文献对照：

- EHT-IV eht-imaging是64²、128μas；EHT-VII eht-imaging是64²、120μas，**不同研究/阶段配置不同**。
- EHT-VIII比较模拟图使用**20μas FWHM**圆Gaussian。
- EHT-VII p10给eht-imaging、DMC、THEMIS的匹配Gaussian为**19μas FWHM**，LPCAL20、polsolve23；不能把20μas视为所有算法唯一精确beam。

### 4.3 几何与模拟参数

| 值 | 项目内出处 | 文献核对/状态 |
|---|---|---|
| `a*=0.9375` | 草稿基准模型；光追代码 | EHT模拟网格含±0.9375（文中部分四舍五入成0.94）；这是选择的模拟自旋，**不是EHT精确测得的M87*自旋** |
| `i=163°` | 草稿、旧稿注释、raytracing | EHT-V/VIII正自旋/顺行模型采用163°（也探索158/168），与远离观测者的角动量向量一致；对应approaching-jet与视线约17°，但不能不说明轴的正负方向就互换 |
| `PA≈288°` | 旧稿注释approaching jet | EHT-V p3引用Walker et al.(2018)，喷流天球方向约288°；文中也展示约±10°范围。PA通常east of north，不是数组x轴角 |
| `108°=0.6π` | 旧稿spin方向、ImageS旋转 | 288°的反方向是108°，数值吻合；但**数组旋转角是否应直接等于天球PA，需要原图基底/翻轴/零方向证明**。文献方向本身不能验证整段Q/U变换 |
| `0.6 Jy` | 草稿、EHT-imaging超参 | EHT-IV/VII有0.6Jy采用值；EHT-VIII模拟平均compact flux采用约0.5Jy。需确认当前RAPTOR归一化，不应自动强行统一 |
| `230 GHz` / 1.3mm | 论文与数据命名230 | 与EHT2017一致；真实low/high band中心与带宽待EHT-II/III查证，不能把文件名当精确频率metadata |

### 4.4 分析定义、模型参数与拟合产物（不能用EHT替代确认）

| 参数/定义 | 当前来源 | 核对方向 |
|---|---|---|
| 1D角窗口10°、采样步1° | 旧积分/新polarization | EHT-VII p10明确采用10°angular section、以1°步长滑动；选择36个不重叠窗口是后续分析选择，不自动证明统计独立 |
| `I×Q / I×U`的额外I加权 | 旧ImageS和草稿 | **EHT-VII p10 Eq.(11)确实采用intensity-weighted coherent EVPA**，支持旧方法有明确文献定义；不能只凭Q/U为Stokes就称重复权重bug。Eq.(11)的Q/U比值及EVPA基底须结合坐标约定核验，不可机械替换atan2顺序 |
| `AXION_NORM=0.1` | 旧result/新model | RAPTOR参考轴子强度，不是EHT测量常数；需原模拟配置验证归一化`2πg f`等约定 |
| `omega(alpha)`中的`4²` | 旧result/新model | 束缚态模型参数；草稿描述211态，但n的约定需要Dolan/源径向函数独立核验，不能由EHT论文确认 |
| `b=0.03(α/0.05)^(7/4)` | 草稿引用Baryakhtar等、model.py | 轴子云理论输入，需查对应理论论文和适用条件 |
| `grmhd_sigma_degrees=1.9752611467515768` | 旧GRMHD日差拟合/limits.py | 数据拟合产物；需回到R4输入、mask、时间差和拟合，不存在EHT统一“正确常数” |
| I>0.1peak、r≥2、人工剔除区域 | EHT处理脚本 | 空间筛选/分析选择；EHT-VII展示图有I>10%peak等cut，但不能据此证明本项目mask全体选择等价 |
| 1572topset、4800候选图、141时间片 | 当前项目数据 | 数据规模，不是M87*物理常数；应按数组/metadata/selection验证 |
| coupling prior边界1/137、1000、offset4 | 最终推断参数 | 先验/抽样选择，需统计与理论论证；不是EHT参数 |

### 4.5 旧代码清单补充：需要单独核验的约定与版本差异

以下为代码/草稿直接提取，尚未完成相应理论文献和生成配置的独立核对：

| 数值/约定 | 位置 | 待核对事项 |
|---|---|---|
| `2*(1+cos(2/3*arccos(0.9375)))` | `GRMHD/main.ipynb::photon_sphere`、`Plot/fig1.ipynb::r_ps` | Kerr赤道圆光子轨道的一支；不能混称为唯一photon sphere、屏幕critical curve或亮环半径。还需结合光追中a符号判断分支 |
| `R_photon_ring≈5.14266475` | `GRMHD/rotate.ipynb`人工像素读数；`result.ipynb`注释 | 图像平面读图估计，不是同一BL轨道半径，也不是独立测量常数 |
| `c=2πg f`与`c=2πb g f` | 当前result注释/信号乘b；旧稿`old.tex_copy:310` | **定义冲突**：差一个随α变化的b；限制纵轴与g转换必须确认采用哪个定义 |
| `10*2π` | `Plot/fig4.ipynb`振幅归一化 | 参考轴子幅度的单位响应换算，不是rad→degree |
| `1022.0` | `GRMHD/GRMHD_cut.ipynb`解析场相位 | 固定时间/相位原点，未找到解释；自由全局相位可吸收不代表跨模板逐点比较可忽略 |
| `0.29915625` | 同上 | 与现`omega(0.3)`一致，重复硬编码的模型频率 |
| `π/8`、`6π/8` | 16状态相位注入；`Plot/fig4.ipynb`相位参考 | 前者为16个相位步长，后者为135°图示参考扣除；均非喷流PA |
| `r_o=500, M=1, a=-0.9375` | 主Notebook光追调用 | 无量纲有限观察者与几何单位，不是500Mpc或1太阳质量 |
| 光追失败返回`(200,0,0)` | Notebook `_calculate_source_position()` | 200是失败占位半径，需显式过滤，不是物理源位置 |
| BL圈选`r_width=0.5`、`8.75≤r≤9.25` | `GRMHD/main.ipynb`诊断 | 源平面圈选，区别于屏幕mask；旧实现还有`w-0.5≤r≤w`。ImageS当前mask只作用I，不能假定Q/U已圈选 |
| toy时延项沿x vs旧稿沿y | `GRMHD/main.ipynb`解析phase；旧稿约174–181行 | 可能是轴定义不同，需显式变换，不能仅核对17°数值 |
| `atan2(x,y)` vs `[-ρsinφ,ρcosφ]` | obs_phi与`get_360_x_y()` | 按通常x向右、y向上前者从+y顺时针、后者从+y逆时针；加上天球RA方向后须整体核验 |
| `mean=-0.147, sigma=1.97`、保存拟合`1.95526/1.97526` | `result.ipynb`日差直方图与输出 | 图示近似、不同拟合结果与最终默认sigma不能混称同一精确数据产物 |
| 17×17搜索、径向步0.2、10个子采样、固定offset | EHT定位/积分函数 | 分析算法选择/配准产物，不是M87测量量；两阶段offset不同需要保留来源 |
| `threshold=0.053` | EHT旧1D fit | 拟合sigma相对误差接受阈值，不是5.3%观测系统误差 |
| 786 vs1572，`P>0.1e9` vs相对peak cut | Introduce、旧result copy、当前topset及Plot/fig3 | 不同分析版本，绝对与相对阈值不可替换；旧图不一定代表当前样本 |
| coupling上界100 vs1000，40 vs60点 | 旧稿/旧单元与最终函数 | prior/grid版本不同；需要标记effective配置 |
| `f_a=10^15GeV`、670像素、`I>0.15Imax`、`a=0.99` | 旧稿旧图注 | 旧结果假设，不能直接用于当前0.9375/757像素结果 |

当前代码所称的“常数”实际还包括很多参数约定、拟合产物和历史版本差异；逐层整理时按此清单确认，不能简单统一替换成一个EHT数值。

## 5. 待逐项完成的核对任务

- [x] 一手EHT-I/IV/V/VII/VIII论文下载与关键段落核对。
- [x] 区分6.5e9/16.8与6.2e9/16.9两套文献采用值。
- [x] 独立重算t_g、θ_g、α→m_a及beam的σ/FWHM。
- [x] 找到项目内模拟标签单位说明。
- [ ] 验证本次RAPTOR实际M/D/FOV/frequency/flux/倾角和数据字段单位，不能由通用论文代替。
- [ ] 验证4800图pickle对应的实际FOV、pixel scale、天空方向与Stokes基底；必要时找原成像脚本/UVFITS/FITS。
- [ ] 查EHT-II/III及实际观测时刻，确认6h矩形平均近似的精度。
- [ ] 确认Gaussian目标是19还是20μas及是否已经在观测图中处理，避免重复blur。
- [ ] 用坐标/基底测试确认108°、翻轴、Q/U负号与PA的关系。
- [ ] 核对轴子理论频率态、云幅度、耦合归一化，不将其归给EHT。
- [ ] 核对R4拟合sigma及26.68h近似与24h目标的差异。
- [ ] 逐个整合重复常数为配置/派生量；在方法确认前不修改数值。

## 6. 后续规则

1. 每个修改都要记录旧值、来源、公式、新值和影响阶段。
2. 改质量/距离可能同时改变时间、角尺度、轴子质量轴和blur单位，不能只改一处。
3. 模拟库既定标定可能不能在后处理随意重标，尤其辐射传输归一化/电子温度/Faraday项；需要与生成配置对应。
4. 物理值在误差范围内接近，不等于实现使用了同一套单位；内部一致性与外部测量一致性分别验证。
5. 所有结论保留“确认/近似/待查/冲突”状态，不将猜测写成已验证事实。
