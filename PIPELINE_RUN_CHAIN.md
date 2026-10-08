# 从 RAPTOR / EHT 成像输出到最终结果：当前代码运行链

## 0. 文档目的、起点和状态

本文是**当前代码的真实运行链与输入/输出清单**，用于接下来逐段整理旧代码、精简和修复问题。它不是一键 pipeline，也不是对现有科学结果正确性的认证。

- 代码核对基线：`74ebae1`，分支 `data/analysis-inputs-draft`。
- 分类参见 [`CODE_CLASSIFICATION.md`](CODE_CLASSIFICATION.md)。
- M87*物理常数、单位换算和分析参数的文献核对参见 [`M87_CONSTANTS_AUDIT.md`](M87_CONSTANTS_AUDIT.md)；来源包含EHT-I/IV/V/VII/VIII原论文。
- 本文不移动文件，不修改算法，也不实际执行昂贵计算。
- Notebook cell 编号从 **0** 开始；后续增删单元后编号可能变化，应同时按函数名定位。

### “原始数据”在本文中的含义

**RAPTOR 路线起点：**RAPTOR 已运行完成的 `img_data_*.h5`，或这些 HDF5 已转换得到的 IQUV `.npy`。

**EHT 路线起点：**EHT-Imaging 已完成重建的四天各 4800 张候选图及超参数 `.pkl`。

不包含更上游的 GRMHD 流体模拟快照、RAPTOR 执行本身、EHT visibility/UVFITS 校准与 4800 图成像任务；当前仓库没有这些完整上游运行链。

### 约定

```text
α = 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44
D = 95, 96, 100, 101     # April 5, 6, 10, 11
T = 141                 # 1600到3000，索引标签间隔10
S = 17                  # 无轴子状态 + 16个初相位
Stokes通道 = I, Q, U, V
```

### 模拟时间标签与物理时间：已找到的项目内说明

来源不是文件名推断，而是以下三处相互印证：

1. 草稿旧稿 [`../EHTaxion_imaging_pipeline_constraint/old_images/old.tex_copy`](../EHTaxion_imaging_pipeline_constraint/old_images/old.tex_copy) 第741行：

   ```text
   001500.dat-001501.dat is 0.1rg, 001600.dat-001610.dat is 1rg
   ```

2. [`Introduce.md`](Introduce.md) 第2行：`141个吸积盘，时间间隔1R_g`。
3. `result.ipynb` cell10：`grmhd_step=10`，并用`30 * 0.1 * 8h = 1day`解释标签与近似时间；这是以8小时近似引力时间的旧计算。

这里时间单位写作`rg/R_g`时，应规范为引力时间 \(t_g=GM/c^3=R_g/c\)，而不是长度。于是：

```text
标签差1   → 0.1 t_g
标签差10  → 1 t_g    → 当前NumPy相邻两帧
标签差30  → 3 t_g    → 旧代码的“一天”差分
```

草稿 `main.tex:192–193` 使用黑洞质量 `M=6.5×10^9 M_sun`，因此：

```text
t_g                         ≈ 8.89325小时
一个标签单位                 ≈ 0.889325小时
相邻NumPy帧                  ≈ 8.89325小时
旧代码3个帧间隔              ≈ 26.67974小时，不是精确24小时
1600→3000的首末跨度          = 140 t_g ≈ 51.87728天
```

141帧只有140个首末间隔。可明确写相对时间：

```python
elapsed_hours = (snapshot_label - 1600) * 0.1 * tg_hours
```

此式不需要假定标签0对应哪个绝对时间原点。

`model.py`中的`TIME_CONVERSION=0.112411`对应`1/t_g≈0.1124 hour^-1`（其倒数8.89593小时），与上述质量尺度近似一致，但略有数值差别。**轴子模板里的精确24小时与GRMHD背景里3帧≈26.68小时是两个不同的时间处理，不应混称为同一个精确“一天”。**

`GRMHD/main.ipynb` cell6把数组时间步写成`0.1 R_g`，与当前每隔10标签采样的数组、Introduce和旧稿说明不一致；该注释混淆了单个标签步与数组采样间隔。当前草稿 `main.tex:416,422` 的snapshot数/总天数/配对数仍是`999/1234/888`占位，不能作为实际数据规模证据。

上述证据已足以说明项目原有时间约定；后续若有原始HDF5时间属性或GRMHD输出配置，可进一步作独立核对。

后文路径均相对仓库根目录，除非特别说明。

### 环境与运行目录

```bash
cd /path/to/EHT-axion-pixel-by-pixel-analysis
python -m pip install -e '.[dev]'
export PYTHONPATH="$PWD"
export MPLBACKEND=Agg
```

- 最终推断/新脚本通常从仓库根目录运行。
- `GRMHD/alpha_*.py` 和 `GRMHD/main.ipynb` 使用相对路径，需从 `GRMHD/` 运行。
- `EHT-Imaging/main.ipynb` 需从 `EHT-Imaging/` 运行。
- 上游旧脚本直接 import `h5py`，但 `pyproject.toml` 没有声明它；需要另外安装 `h5py`。
- 运行 Notebook 还需 Jupyter/IPython 环境，不能假定 `pip install -e .` 已包含。
- 有耗时任务时使用后台 tmux，限制 BLAS 线程；不要从 `Run All` 开始盲跑混合 Notebook。

## 1. 总览：两路输入汇合

```text
RAPTOR HDF5                           EHT-Imaging四天4800图pickle
     │ R0                                         │
     ▼                                            │
原始IQUV NumPy cubes                              │
     │ R1                                         │ E0/E1
     ▼                                            ▼
处理后IQUV、2D/旧1D EVPA                 topset排名 → 1572图 → 平移/EVPA统计/mask
     │ R2                                │              │
     ▼                                   │              ├── 2D观测输入
141时间片的2D A/φ                        │              │
     │ R3                                └──────┐       │
     ▼                                          │       │
2D A/φ均值与散布                                 │       │
                                                │       │
处理后RAPTOR IQUV ────────────────┐               │       │
                               │ J1             │       │
                               ▼                ▼       │
                            新1D积分/拟合 ← corrected EHT IQUV
                               │ J2                     │
                               ▼                        │
                            1D观测/模板输入              │
                               │                        │
无轴子EVPA时间序列 ── R4 ── GRMHD日差误差                 │
                               │                        │
                               └─────────┬──────────────┘
                                         │ S1/S2
                                         ▼
                           std/cor/cor_grmhd × 1D/2D
                                         │
                                         ▼
                              每个α的log10(c95)
                                         │ F1
                                         ▼
                               六曲线图 + CSV/JSON
```

图4等离子体/发射贡献拆解是独立诊断旁路，见第8节；不是这条六曲线链的必要步骤。

## 2. 共享输入与本地旧路径的对应

### 2.1 已发布的 analysis-inputs

```text
shared_data/analysis-inputs/draft/
├── raptor-stokes/
│   ├── res-100/alpha-{α}.npy                 # 8个完整数组
│   └── res-400/alpha-0.3/stokes-t*-t*.npy    # 15个时间分片
├── eht-image-ensemble/230_{D}.pkl           # 4天各4800候选图
├── raptor-transport-diagnostics/           # 图4旁路
└── manifest.json
```

下载及完整性检查：

```bash
git lfs install
git lfs pull --include='shared_data/analysis-inputs/draft/**'
python scripts/verify_shared_data.py shared_data/analysis-inputs/draft
```

此校验器验证整个 manifest；只下载子集时，其他文件会报告 missing。不要把部分下载误认为全量已验证。

### 2.2 当前旧代码实际读取的位置

| 共享文件 | 当前上游代码要求的位置 |
|---|---|
| `raptor-stokes/res-100/alpha-{α}.npy` | `GRMHD/data_in/GRMHD_Data_1600-3000_alpha={α}_res=100.npy` |
| `raptor-stokes/res-400/alpha-0.3/stokes-t*-t*.npy` | 先按时间axis拼成 `GRMHD/data_in/GRMHD_Data_1600-3000_alpha=0.3_res=400.npy`，或将来给处理入口增加分片支持 |
| `eht-image-ensemble/230_{D}.pkl` | `EHT-Imaging/data_in/230_{D}.pkl` |

**断点 A：**当前处理脚本没有统一的 `--input-dir` 来直接读取这些共享路径，也没有上游处理入口原生读取 res=400 shards。现在需手工放入旧路径；这是整理时优先解决的接口问题。

res=400 的分片按真实时间标签排序，前14个 shape为 `(10,17,400,400,4)`，最后1个为 `(1,17,400,400,4)`。若必须重建原文件，应用 `open_memmap` 分片写入，不能在有限内存机器上直接 `np.concatenate` 分配整个12.27 GB数组。

## 3. RAPTOR 路线

### R0：HDF5 → 原始 NumPy IQUV cube

| 项目 | 内容 |
|---|---|
| 代码文件 | `GRMHD/rotate.ipynb` |
| 函数 | `readh5()`；旧批量原型 `readh5ALL()` |
| 输入 | 各时间片/相位的 `.../data{time}.dat/img_data_{state}.h5` |
| 单图输出（内存） | `(resolution,resolution,4)`，通道按IQUV解释 |
| Notebook明确保存示例 | `data/grmhdData1600.npy`，shape `(17,100,100,4)` |
| 当前下游要求的汇总输出 | `GRMHD/data_in/GRMHD_Data_1600-3000_alpha={α}_res={res}.npy` |
| 汇总shape | res100：`(141,17,100,100,4)`；res400：`(141,17,400,400,4)`；float64 |

`readh5` 将HDF5 block/pixel重排成完整图像。这个阶段按意图是格式转换/像素重排，不做后面的Gaussian blur。

**断点 B：**不能直接运行旧 `readh5ALL()` 就认为能得到发布的全部cubes：

- 它把实际时间标签 `t=1600...` 直接当数组索引，而不是从0开始的时间索引；
- `range(start,end,gap)` 不含3000，且分配140片，而发布数组有141片；
- 其中硬编码res100及外部绝对路径；
- 按HDF5 key列表位置而非明确Stokes字段名读通道。

因此目前没有确认可靠的完整HDF5批量转换CLI。发布的NumPy数组可以跳过R0；若要从HDF5重新开始，必须先整理和验证这一段，**不能声称当前原型已证明全部cubes的生成过程**。

转换正确、完整且元数据保留后，后续主链不再读取HDF5；HDF5主要作为来源档案和转换核查依据。

### R1：原始 IQUV → 处理后的图像与 EVPA

| 项目 | 内容 |
|---|---|
| 代码入口 | `GRMHD/alpha_0.1.py`、`alpha_0.15.py`、`alpha_0.2.py`、`alpha_0.25.py`、`alpha_0.3.py`、`alpha_0.35.py`、`alpha_0.4.py`、`alpha_0.44.py`；高分辨率为`alpha_0.3_high_res.py` |
| 旧Notebook对应 | `GRMHD/main.ipynb` cell7的`getGRMHDDataAfterImageS()`和cell9 |
| 共用实现 | `utils/libanaly.py::ImageS` |
| 输入 | R0输出的9个原始IQUV cubes |
| 处理 | Gaussian blur；旋转图像及Q/U偏振基底；插值/裁切；计算2D EVPA与旧入口1D EVPA |
| 输出目录 | `GRMHD/data_out/` |

每个输入α/res产生：

```text
grmhdDataAfterImageS_alpha={α}_res={res}.npy   (141,17,64,64,4)
EVPA_alpha={α}_res={res}.npy                  (141,17,64,64)，rad
EVPA_1d_alpha={α}_res={res}.npy               (141,17,360)，rad
EVPA_1d_10deg_alpha={α}_res={res}.npy         (141,17,360)，rad
```

典型命令（需要先确认处理参数与科学定义）：

```bash
cd GRMHD
python alpha_0.3.py
python alpha_0.3_high_res.py
```

其余α逐个运行相应脚本。九个脚本大段重复，将来应合并。

重要说明：

- 文件名中的`res=400`指原始输入分辨率，不代表输出400；默认`final_shape=64`。
- 默认旋转`0.6π`，原始FOV `46.2262`，最终FOV `32.6868`；这些需要有物理/坐标说明。
- 旧入口会额外计算1D EVPA；当前新1D重建主要从processed IQUV重新积分，不能把两套结果混用。
- 原积分的I权重/新积分的直接Q/U、角窗口中心等定义仍待确认。本链只记实际代码，不默认新旧任何一套正确。
- 脚本顶层`RERUN_ANALYSIS=True`，import它也可能触发重算，不能把它直接当安全库模块导入。

### R2：每个时间片的16轴子相位 → 2D A/φ拟合

| 项目 | 内容 |
|---|---|
| 代码文件 | `GRMHD/main.ipynb` |
| 相关单元 | cell10加载中间量；cell15定义；cell18调用 |
| 函数 | `fit_cos_evpa()`、`estimate_cos_params()`、`fit_evpa_time_series()` |
| 输入 | `GRMHD/data_out/EVPA_alpha={α}_res=100.npy`；旧1D拟合还读`EVPA_1d_10deg_*` |
| 方法 | 用无轴子状态作角度展开参考，对16初相位做余弦/DFT拟合 |
| 输出 | 见下表 |

```text
EVPA_cos_fit_result_alpha={α}.npy
shape = (141,64,64,3)
最后轴 = [A, φ, baseline EVPA]

EVPA_cos_fit_result_1d_10deg_alpha={α}.npy
shape = (141,360,3)
```

运行时只选取必要定义和调用单元，不执行整个混合Notebook。`fit_evpa_parameter_scan()`（旧`different_omega/`扫描）与`process_alpha_0_3_res_400()`（输出256²诊断拟合）不是最终默认六曲线的必需步骤。

### R3：141个时间片的2D A/φ → 模板均值与散布

| 项目 | 内容 |
|---|---|
| 代码文件 | `result.ipynb` |
| 函数/单元 | cell33加载余弦拟合；cell36的`get_A_phi0_params()`与保存语句 |
| 输入 | R2的`EVPA_cos_fit_result_alpha={α}.npy`；该旧函数还读取processed IQUV |
| 处理 | 对每个像素141个A和φ样本拟合分布，处理φ的周期性 |
| 输出 | `GRMHD/data_out/A_params_dict.npy`、`phi0_params_dict.npy` |
| 输出schema | dict：α→array `(2,64,64)`；第0行mean，第1行std |

**断点 C：**这一段尚未完全迁出旧Notebook，当前final-inference发布器直接复制本地已有字典，不会从RAPTOR原始cube重新生成它。

此外，当前`fit_gaussian2`对`p_value_threshold`的行为与历史生成保存字典时不同。不能只按旧cell运行就保证与保存字典相等；必须确认要保留的拟合/筛点定义。

一维字典的旧算法是同一Notebook cell50的`get_A_phi0_1d_params()`；当前新一维链改用J1/J2，二者不要同时作为同一结果的生成步骤。

### R4：无轴子时间演化 → GRMHD日差误差

| 项目 | 内容 |
|---|---|
| 旧代码 | `result.ipynb` cell8–22 |
| 函数 | `get_EVPA_diff_xday()`、`get_EVPA_diff_1day()`、`fit_gaussian2()` |
| 输入 | `EVPA_alpha=0.3_res=400.npy`的state0；processed res400图像；EHT空间condition |
| 当前旧mask路径 | `EHT-Imaging/data_out/condition.npy` |
| 处理 | 取指定模拟时间间隔EVPA差，选点，拟合统一散布 |
| 输出 | 内存标量`_cor_grmhd_tmp_`（degree），内存矩阵`cor_grmhd`；没有独立规范的误差参数文件 |

**断点 D：**最终runner没有调用R4，而在`eht_axion/limits.py::LimitOptions`中写死：

```python
grmhd_sigma_degrees = 1.9752611467515768
```

因此“从原始数据重算full模式”还必须把这一估计和它的mask/模拟时间换算显式串入，不能只复用常数就说全部从头重算。项目旧稿明确说明每标签为0.1 t_g、每采样间隔为1 t_g，旧代码以t_g≈8小时把3个采样间隔近似为一天。按草稿质量6.5×10^9 M_sun，实际差分间隔约26.68小时；需决定保留这一离散近似还是估计精确24小时对应的背景误差。

R4依赖EHT mask，应在EHT处理完成之后执行；当前旧mask与新757像素mask也需要明确选哪个。

## 4. EHT 路线

### E0：历史EHT处理中间量（现代码的额外依赖）

| 项目 | 内容 |
|---|---|
| 代码 | `EHT-Imaging/main.ipynb` |
| 输入 | `EHT-Imaging/data_in/230_{D}.pkl` |
| 功能 | 定位、topset、EVPA拟合、保存Stokes和一维观测量 |
| 输出 | `EHT-Imaging/data_out/{D}_all_data.npy`、`{D}_all_data_mean.npy`、`{D}.npz`、`condition.npy`、`{D}_EVPA_1d.npy`、`{D}_1d.npz` |

这一Notebook是旧算法出处，不是当前最终runner依赖。**但是**下面新重建脚本仍强制读取其中的历史mask/1D数据用于对照报告，所以干净环境仍会遇到这些依赖。

注意：当前Notebook的排名对齐和`ImageS`已被改动，重新运行生成的结果不一定是过去的“历史基准”。不能为了满足对照路径而随意重新生成后称其为旧结果。

**整理目标：**把历史对照设为可选，不应要求研究人员先跑旧链才能跑新链。

### E1：4800图 → corrected topset排名

| 项目 | 内容 |
|---|---|
| 代码 | `scripts/validate_topset_ranking.py` |
| 输入 | `EHT-Imaging/data_in/230_{D}.pkl`，四天各4800图与meta_data |
| 处理 | 按共同超参数字段排序；空间/Stokes变换；日EVPA分布拟合；平移对齐后计算偏离分数 |
| 临时输出 | `audit/topset_work/order_{D}.npy`、`theta_{D}.npy`、`fit_{D}.npz` |
| 下游必要输出 | `audit/topset_ranking_comparison.npz`中的`corrected_order` |
| 审计输出 | 同NPZ中的legacy分数/顺序与`topset_ranking_comparison.json` |

```bash
python scripts/validate_topset_ranking.py
```

它虽然名为validate，当前却实际生产了E2需要的排名，不能当纯测试删除。

### E2：排名 → 四天1572图及2D观测输入

| 项目 | 内容 |
|---|---|
| 代码 | `scripts/rebuild_corrected_topset.py` |
| 输入 | E1的`corrected_order`；四天原始4800图pickle |
| 额外强制输入 | `EHT-Imaging/data_out/condition.npy`，仅用于旧新区域比较 |
| 处理 | 每天选择共同超参索引的1572图；旋转；拟合EVPA；固定offset重定位；构建mask |
| 输出目录 | `audit/rebuilt_2d/` |

```text
topset_indices.npy         (1572,)
{D}_all_data.npy           (1572,64,64,4)，当前V通道为0占位
{D}_all_data_mean.npy      (64,64,4)
{D}.npz:
    arr_0                 theta_mean (64,64)，rad
    arr_1                 theta_std  (64,64)，rad
    arr_2                 P_mean     (64,64)
    arr_3                 I_mean     (64,64)
condition.npy             (64,64)，bool；当前757像素
report.json               处理报告
```

```bash
python scripts/rebuild_corrected_topset.py
```

当前第二阶段offset固定为`{95:(1,0),96:(0,1),100:(-1,2),101:(0,0)}`，不是随新topset重新估计的中心。

**断点 E：**对照旧mask是强制步骤；缺少它时，科学输出可能已经写出，但脚本会在最后失败。并且`topset_indices.npy`先于四天输出写出，不能用它存在来证明整批完成。

## 5. 两路联合的一维重建

### J1：处理后的RAPTOR/EHT Stokes → 新1D EVPA与A/φ

| 项目 | 内容 |
|---|---|
| 代码 | `scripts/rebuild_1d_products.py` |
| 共用库 | `eht_axion/polarization.py`、`angles.py` |
| EHT主输入 | `audit/rebuilt_2d/{D}_all_data.npy` |
| RAPTOR主输入 | `GRMHD/data_out/grmhdDataAfterImageS_alpha={α}_res=100.npy`（R1产物） |
| 额外强制对照输入 | `EHT-Imaging/data_out/{D}_EVPA_1d.npy`；`GRMHD/data_out/EVPA_1d_10deg_alpha={α}_res=100.npy` |
| 处理 | 当前使用直接Q/U径向积分、10bin角窗口、角轴roll(-90)；对RAPTOR16相位做余弦拟合 |
| 输出目录 | `audit/rebuilt_1d/` |

```text
{D}_EVPA_1d.npy                       (1572,360)，rad
EVPA_1d_alpha={α}.npy                 (141,17,360)，rad
EVPA_cos_fit_result_1d_alpha={α}.npy   (141,360,3)，[A,φ,baseline]
report.json
```

```bash
python scripts/rebuild_1d_products.py
```

重要限制：

- `audit/rebuilt_2d/topset_indices.npy`不存在时，代码会静默回退到旧`EHT-Imaging/data_out`；这不是可靠的corrected批次门禁。
- 一维积分没有读取2D `condition.npy`；1D和2D不自动使用同一空间区域。
- 新旧1572图成员和顺序可能不同，逐行相减只能作为集合变更粗诊断，不能解释为同一张图的修复前后差。
- 当前Q/U积分定义、10度窗口中心、径向边界有待逐项确认/修复，不应因输出目录叫rebuilt就默认正确。

### J2：1D ensemble → 日均值/散布与模板参数

| 项目 | 内容 |
|---|---|
| 代码 | `scripts/fit_rebuilt_1d_products.py` |
| 输入 | J1的四天EHT EVPA和八组GRMHD余弦拟合 |
| 处理 | 圆周居中/展开后使用普通Gaussian均值和std（unbinned MLE） |
| 输出目录 | `audit/rebuilt_1d/` |

```text
{D}_1d.npz:
    arr_0                  theta_mean (360,)，rad
    arr_1                  theta_std  (360,)，rad
A_params_1d_dict.npy        dict α→(2,360)，[mean,std]
phi0_params_1d_dict.npy     dict α→(2,360)，[mean,std]
fit_report.json
```

```bash
python scripts/fit_rebuilt_1d_products.py
```

当前无正态性/集中度有效性门禁；360个finite bin不等于360个都符合Gaussian假设。它也不等同于此前histogram-fit算法。

## 6. 统计输入汇合与共享快照

### S1：最终统计实际需要哪些文件

| 维度 | 观测输入 | 模板输入 | 当前存放位置 |
|---|---|---|---|
| 2D | `condition.npy`、四天`{D}.npz`、`{D}_all_data.npy` | `A_params_dict.npy`、`phi0_params_dict.npy` | 观测来自`audit/rebuilt_2d/`；模板来自`GRMHD/data_out/` |
| 1D | 四天`{D}_1d.npz`、`{D}_EVPA_1d.npy` | `A_params_1d_dict.npy`、`phi0_params_1d_dict.npy` | `audit/rebuilt_1d/` |
| full模式额外量 | 统一GRMHD日差sigma | 目前来自硬编码常数，非文件 | `eht_axion/limits.py` |

发布共享快照不是科学计算步骤：

```bash
python scripts/publish_shared_data.py draft --force-draft
```

代码：`scripts/publish_shared_data.py`，读取上表本地中间量，输出：

```text
shared_data/final-inference/draft/
├── 1d/   # 10个文件，包含全部1D观测/模板输入
├── 2d/   # 11个文件，包含全部2D观测/模板输入
└── manifest.json
```

```bash
python scripts/verify_shared_data.py
```

它校验发布文件，不代表所有物理参数、批次来源和统计定义已正确。数据发布的事务性/hardlink隔离等问题也待后续整理。

### S2：最终推断入口与内存计算

代码：

```text
scripts/run_final_limits.py
  → eht_axion/limits.py::calculate_limits
      → observations.py
      → model.py
      → statistics.py
```

**输入**：S1中的1D/2D观测与模板文件。

**运行时生成但不单独保存的量**：

- 选点后的日均值/std；1D固定offset4、每10度取一点，当前36点；2D当前757点。
- 每日观测矩阵与两组日差矩阵。
- April5→6、April10→11的观测差（observed模式）。
- 各α、初相位、耦合下的理论模板和模型误差。
- coupling/phase likelihood与posterior（当前只保存最后上限，不保存全部推断诊断）。

**输出文件**：

```text
{mode}_{dimension}.npz
{mode}_{dimension}.json
```

NPZ包含`alphas`、`cuts`、`variant`、`observation`、`prior`；`cuts`是`log10(c95)`。每条默认8个α点。

### 当前可直接运行的六曲线命令

以下从共享final-inference开始，跳过R0–J2。为了与**当前固定绘图路径**对齐，显式指定输出目录：

```bash
for dimension in 1d 2d; do
  for mode in std cor cor_grmhd; do
    python scripts/run_final_limits.py \
      --dimension "$dimension" \
      --mode "$mode" \
      --observation-dir "shared_data/final-inference/draft/$dimension" \
      --fit-parameter-dir "shared_data/final-inference/draft/$dimension" \
      --output-dir audit/results/rebuilt_corrected_expected
  done
done
```

**注意：该命令会覆盖仓库已保存的同名审计结果。本文件只是记录，不授权自动覆盖或重算。**后续应给绘图增加输入目录参数，使新计算写到独立研究输出目录。

- 默认不使用实际观测差，是零日差/null假设；不自动等同于重复实验的median expected。
- 加`--observed`后应写到不同目录，例如`audit/results/rebuilt_corrected_observed`，不能覆盖expected文件。
- `--legacy`用于历史回归，不应把它的新旧标签当成科学正确性证明。
- runner默认输出`audit/results/corrected_expected/`，**不是**绘图固定读取的`rebuilt_corrected_expected/`；漏写`--output-dir`会有画旧数据风险。

**待修复/确认：**5天相位公式、频率态与时间单位、模型/观测跨时相关性、std/cor边际方差定义、posterior分位反演等。当前数据流可执行，不表示这些问题已解决。

## 7. 最终绘图与数值结果

### F1：限制文件 → 六曲线图和比较图

| 项目 | 内容 |
|---|---|
| 代码 | `scripts/plot_final_limits.py` |
| corrected固定输入 | `audit/results/rebuilt_corrected_expected/{std,cor,cor_grmhd}_{1d,2d}.npz` |
| legacy固定输入 | `audit/results/legacy_expected/{std,cor,cor_grmhd}_{1d,2d}.npz` |
| 输出目录 | `audit/figures/` |

```bash
python scripts/plot_final_limits.py
```

输出：

```text
final_limits_corrected.pdf
final_limits_corrected.png
legacy_vs_corrected.pdf
legacy_vs_corrected.png
final_limits.csv
final_limits.json
```

其中legacy输入也被比较图和CSV/JSON表强制读取，缺少它时即便只想要新六曲线图，整个脚本也可能中途失败。后续应拆分“只画新结果”和“与历史比较”。

旧最终绘图来源在`result.ipynb`的`plotting_script()`和1D/2D比较单元。`Plot/fig*.ipynb`是其他信号/误差诊断图，不是此处新六曲线绘图入口。

## 8. 科学诊断旁路，不要混入六曲线主依赖

### 图4：等离子体开关与发射贡献拆解

代码出处：`GRMHD/main.ipynb`，主要cell43–64：

- 加载`grmhd_data_paths`；
- `get_360_x_y()`、`get_analytical_A_1d()`；
- `get_raptor_evpa()`、`get_analytical_evpa()`；
- `raptor_A_1d` / `raptor_A_1d_blurred`与彩色环、四曲线绘图。

输入共享目录：

```text
shared_data/analysis-inputs/draft/raptor-transport-diagnostics/
  snapshot-2110-res-400/
    common.npy
    common_noaxion.npy
    jonly.npy
    eq_jonly.npy
    fake.npy
    r_coord_400.npy
    axion_field_coord_400_16.npy
```

旧Notebook实际读的是`GRMHD/data_in/grmhdData2110_400_*.npy`和两份coord文件，当前也未配置共享输入目录。

输出主要是内存振幅比和诊断图；不少单元只`plt.show()`，不能假定都存在统一规范的保存文件。其他较早对照会保存到`Plot/images/Axion Signal/evpa_ratio_*.png`。

`common/jonly/eq_jonly/fake`的开关精确语义仍需原RAPTOR配置确认。此旁路不需要加入最终likelihood数据依赖，但画相关论文图时必须保留。

### 原始光线轨迹 → 赤道面参考

代码：`GRMHD/GRMHD_cut.ipynb`与`utils/raytracing.py`。

输入：`GRMHD_cut*`轨迹表；输出/用途：`r_coord*.npy`、`axion_field_coord*.npy`等解析对照数据。当前已共享的小型coord数组可跳过大轨迹表再生。这不是普通IQUV格式转换，也不是六曲线默认入口。

## 9. 验证代码应该插在哪些阶段

| 代码 | 检查什么 | 必要输入/输出 |
|---|---|---|
| `scripts/verify_shared_data.py` | 下载文件大小/hash/可读性 | manifest与所列文件；不计算曲线 |
| `scripts/inspect_intermediates.py` | 源中间量shape、有限值、Stokes均值一致性 | 旧GRMHD/data_out及EHT/data_out；输出`audit/intermediate_validation.json` |
| `scripts/validate_fit_products.py` | 抽样重拟合A/φ | 源cosine fit+参数字典；输出`audit/fit_product_validation.json` |
| `scripts/validate_historical_selection.py` | 从保存Stokes签名反推旧topset | 原4800pickle、旧all_data、排名NPZ；输出`audit/historical_selection_validation.json` |
| `scripts/verify_legacy_results.py` | 数值回归 | 重算legacy输出与`Results/data/`；输出`audit/legacy_reproduction.json` |
| `scripts/check_grid_convergence.py` | 60/120 coupling网格 | 本地rebuilt输入与基准结果；输出`audit/grid_convergence.json` |
| `scripts/run_ablation.py` | 分别切换统计定义 | 当前主要用旧中间量；输出`audit/ablation_2d_cor_grmhd.json` |
| `scripts/verify_all.sh` | 串联测试/回归/重算/绘图/Notebook语法 | 依赖预先生成的旧与rebuilt中间量；不是从R0开始的构建命令 |
| `tests/` | 被覆盖函数行为 | 小型合成输入；不能替代物理定义审查 |

## 10. 此前Notebook路线的对应关系（算法出处，不建议直接Run All）

如果以“以前的代码”为整理起点，原链是：

```text
GRMHD/rotate.ipynb             HDF5 → NumPy
          ↓
GRMHD/main.ipynb 或 alpha_*.py 原始Stokes → processed IQUV/EVPA
          ↓
GRMHD/main.ipynb               相位余弦拟合 → 每个时间片A/φ
          ↓
result.ipynb                   参数分布拟合 → A/φ mean/std
          │
          ├── 无轴子时间差 → GRMHD统一日差sigma
          │
EHT-Imaging/main.ipynb         4800图 → topset → mean/std/观测图像
          │
          └─────────────┐
                        ▼
result.ipynb             2D/1D观测矩阵 + 模板 → get_final_cut/get_final_cut_1d
                        ↓
Results/data/            六组限制文件
                        ↓
result.ipynb             限制曲线/1D-2D比较图
```

`result.ipynb`的旧最终阶段对应：

| 原单元/函数 | 输入 | 主要内存/文件产物 |
|---|---|---|
| cell63：2D观测加载与残差矩阵 | `EHT-Imaging/data_out/{D}.npz`和`{D}_all_data.npy` | 四天`theta_mean/std`、`theta_cor`；全4096像素矩阵 |
| cell89：`mean_topset/std_topset` | 上述日统计 | 展平后的日均值/std |
| cell90：`get_final_cut(mode)` | 日差/矩阵、轴子预测和散布、`cor_grmhd`、condition | 2D上限；写`Results/data/{mode}.npz`、`Results/figure/{mode}.pdf`及额外美化图 |
| cell94–96 | mode依次`std/cor/cor_grmhd` | 三组2D结果 |
| cell98：1D观测加载 | `EHT-Imaging/data_out/{D}_1d.npz`和`{D}_EVPA_1d.npy` | 四天1D均值/std/cor与两组日差；`resolution`改成360 |
| cell99：`get_final_cut_1d(mode)` | 上述1D观测、1D模板、统一GRMHDsigma | 1D上限；写`Results/data/{mode}_1d.npz`和对应PDF |
| cell100–102 | mode依次`std/cor/cor_grmhd` | 三组1D结果 |
| cell112–114：读取full结果并比较 | `Results/data/cor_grmhd{,_1d}.npz`等 | `new_static/contract.pdf`；不是三模式六条统一图 |

旧`Results/data`的每个NPZ主要包含`alphas`与`cuts`，2D文件没有`_2d`后缀，新CLI文件则明确带`_2d`：

| 旧文件 | 新CLI对应文件名 |
|---|---|
| `std.npz`、`cor.npz`、`cor_grmhd.npz` | `std_2d.npz`、`cor_2d.npz`、`cor_grmhd_2d.npz` |
| `std_1d.npz`、`cor_1d.npz`、`cor_grmhd_1d.npz` | 同名 |

旧Notebook还含未定义变量、旧786方案、绘图残留和全局变量复用，因此**这是一份依赖图，不是授权直接顺序执行这些cell的操作脚本**。整理时应把所需函数的输入/输出显式化，并用原保存结果做回归；回归一致不等于旧公式正确。

## 11. 接下来按什么顺序整理

每步先确认旧算法意图，再明确输入输出、迁出必要代码、加测试；不要一次性改所有层级。

1. **R0：HDF5转换**——如果已确认只从发布NumPy开始，可先保留为归档工具；否则先验证通道和坐标映射。
2. **R1：RAPTOR预处理**——合并九个alpha脚本，支持共享目录及时间shards，明确blur/旋转/FOV/输出分辨率。
3. **R2/R3：模板拟合**——迁出旧Notebook的余弦和2D分布拟合，明确新旧1D估计器和筛点规则。
4. **E1/E2：EHT处理**——统一共享输入路径，把旧对照设为可选，验证重定位和topset身份对齐。
5. **J1/J2：1D处理**——确认科学积分定义、窗口和区域，再修数值边界，统一两路生成。
6. **R4：GRMHD误差**——生成明确的数据产物而不是硬编码sigma，说明模拟时间/观测时间换算与mask。
7. **S2：统计推断**——首先确认时间相位公式，再处理误差相关性、先验和数值积分。
8. **F1：画图**——输入目录显式参数化；只读数值结果，不混入推断。
9. **发布/验证**——将上述产物按manifest版本化共享，验证不再依赖作者机器上的隐藏路径。

### 本文明确标出的当前断点

- A：共享输入路径/shards与旧上游入口不兼容。
- B：完整HDF5汇总原型有索引和元数据问题，不能直接运行。
- C：2D模板参数生成仍在旧Notebook，新发布器只复制已有参数。
- D：GRMHD日差sigma没有从原始数据自动接入最终runner。
- E：旧对照文件是新重建的强制依赖，可能在输出半成品后失败。
- 绘图固定目录与runner默认目录不一致。

**这些是整理任务的边界，不是本步已经实现的修复。**
