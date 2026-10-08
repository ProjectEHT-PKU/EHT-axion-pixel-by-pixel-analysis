# 无 notebook EHT stage

实现：`eht_axion/pipeline/eht.py`。仅依赖 NumPy 与父侧的 SciPy polarization
operator；不导入 notebook、`utils.libanaly` 或 statsmodels，不调用历史 Gaussian
histogram fit，也不读取历史输出做 mandatory 对照或静默 fallback。

## CLI 与 transaction 集成

```bash
python -m eht_axion.pipeline.eht \
  --input-dir data/inputs/eht \
  --output-dir /path/to/new/run/eht \
  --topset-size 1572 \
  --trust-pickle
```

**pickle 可执行任意代码。** CLI 必须显式传 `--trust-pickle`；函数调用必须传
`trust_pickle=True`。仅对可信来源输入使用此开关；它不是安全沙箱。

直接 CLI 的输出目录必须为空或不存在，拒绝覆盖已有结果。直接运行不提供原子
transaction/cache：失败时可能留下部分产品，但不会写 `report.json`，临时角度文件
会清理。`report.json` 最后写入，但不能替代父侧校验后的 `manifest.json`。
命名总约定见 [DATA_LAYOUT_PLAN.md](DATA_LAYOUT_PLAN.md)。函数接收调用者指定的
EHT 产品目录，不自行查找 run root 或追加 `eht/`；observations 维度目录也由父侧定位。

父侧可用 `common.stage_directory` 包装同一函数，传入其空暂存目录：

```python
from pathlib import Path
from eht_axion.pipeline.common import stage_directory
from eht_axion.pipeline.eht import DAYS, run_eht_stage
from eht_axion.pipeline.layout import input_eht_file

inputs = Path("data/inputs/eht")
config = dict(trust_pickle=True, topset_size=1572,
              centering="ensemble", weighting="intensity")
with stage_directory(
    Path("runs/example/eht"), "eht", config,
    [input_eht_file(inputs, day) for day in DAYS],
) as temporary:
    if temporary is not None:
        report = run_eht_stage(inputs, temporary, **config)
```

缓存 config 必须覆盖调用者的全部科学参数（包括非默认 offsets、mask、operator
设置等），不能仅使用本示例的最小 config。使用 injected operator/integrator 的
调用者负责在 config 中声明其内容/版本。EHT 函数自身不管理 cache，不根据目录存在
复用任何科学输出；父侧负责输入与源码 hash、完整性验证及 transaction。

### Python API

```python
run_eht_stage(
    input_dir, output_dir, *,
    trust_pickle=False,
    topset_size=1572,
    centering="ensemble",
    weighting="intensity",       # 或 "flux"
    inherit_2d_mask=False,
    fov=32.6868,
    radial_max=None,
    window_degrees=10,
    angular_bins=360,
    angle_roll=-90,
    min_radius=2.0,
    intensity_fraction=0.1,
    historical_mask=True,
    exclude_columns=5,
    exclusion_rectangles=((48, 55, 10, 19),),
    batch_size=64,
    pixel_batch_size=128,
    rank_offsets=None,           # None 使用下表 fixed offsets
    final_offsets=None,
    operator=None,
    integrator=None,
) -> dict
```

所有参数都具有上述函数默认值。CLI 暴露 `--centering`、`--weighting`、
`--inherit-2d-mask`、`--fov`、`--radial-max`、`--window-degrees`、`--angular-bins`、
`--angle-roll`、`--min-radius`、`--intensity-fraction`、`--no-historical-mask`、
`--exclude-columns`、`--batch-size`、`--pixel-batch-size`；任意 rectangle 和
自定义 fixed offsets 通过 Python API 配置，不做自动拟合。

`angle_roll` 的单位是 **bin**。默认 `360` bins 时 `-90` 即历史的 -90° 排列调整；
自定义角度网格必须显式选择对应 roll，例如 180 bins 的 -90° 是 -45 bins。

未注入时明确调用父侧 API：

```python
operator = polarization.aperture_operator(
    shape=shape, fov=fov, radial_max=radial_max,
    window_degrees=window_degrees, angular_bins=angular_bins,
)  # CSR, (angular_bins, shape**2)
evpa = polarization.integrate_stokes(stokes, operator, weighting=weighting)
# evpa: (batch, angular_bins), radians, EVPA period pi
```

这些父侧函数缺失或结果形状错误会直接报错，不退回 `integrated_evpa`、旧
`radial_integration_matrix`、逐图 spline 或 flux weighting。

## 输入与 identity 对齐

固定四天顺序：`(95, 96, 100, 101)`，对应 ISO 日期
`(2017-04-05, 2017-04-06, 2017-04-10, 2017-04-11)`。
日期标签统一使用 `layout.date_label(day)`；输入由
`layout.input_eht_file(input_dir, day)` 解析，支持 canonical
`{date}/image-ensemble.pkl` 和旧 `230_{day}.pkl`，也可传父侧支持的输入根目录。
只读取原数据，不重命名、复制或写回历史 pickle。每个 pickle 包含：

- `meta_data`：长度 N 的 mapping 序列；
- `I_data`、`Q_data`、`U_data`：相同形状 `(H, H, N)`；H≥2；
- 可选 `V_data`：与 I 相同形状。没有 V 时输出 V=0，padding 除外。

生产输入为每一天 N=4800、H=64；实现允许小型 synthetic 输入，不硬编码 N/H。
`topset_size` 必须不超过 N，四天 H 与 identity 集合必须一致。变换批次检查原始
Stokes 必须有限；不把坏输入静默转换成零。

identity 是共同的完整 tuple：

```python
FIELDS = ("zbl", "prior_fwhm", "sys_noise", "MEM", "TV", "T2V", "l1")
identity = tuple(metadata[field] for field in FIELDS)
```

每一天先验证字段存在、有限 scalar、tuple 唯一。第一天 lexicographic sort 定义
canonical identity 顺序；其余三天必须与其 **identity set 完全相同**，随后按同一
identity lookup 映射 raw 行号。重复、缺失、新增 identity 均报错，不能只盲目各自排序
后假定相同位置对应同一图。`chi2` 等不在 FIELDS 中的字段不是 identity。最终产品
第二遍读 pickle 时重新验证 metadata/顺序及图像 shape。

输出 identity 是 JSON tuple 编码后的 Unicode 字符串，每个元素可用 `json.loads`
读回；所有 NPZ 均可 `allow_pickle=False` 加载，避免 object arrays 与 metadata
数值/字符串类型混合造成的静默精度损失。

## 空间变换、对齐排名与 fixed relocation

逐通道遵循 `scripts/validate_topset_ranking.py`：

1. 在 raw `(H,H,N)` 上翻转第 0 轴；
2. 转为 `(N,H,H)`；
3. 空间轴 `np.rot90(..., k=-1)`；
4. **I 无负号，Q/U 有负号**。可选 V 只做空间变换、不改符号。

`shift_spatial(array, (dx,dy))` 作用于最后两个空间轴，定义为
`output[...,x,y] = input[...,x+dx,y+dy]`，超出边界填 NaN，不使用周期 roll。
Stokes 先转成 `(N,4,H,H)` 后 shift，再转回 `(N,H,H,4)`；通道轴不平移。

| day | rank fixed `(dx,dy)` | final fixed `(dx,dy)` |
| --- | --- | --- |
| 95 | `(2,0)` | `(1,0)` |
| 96 | `(1,1)` | `(0,1)` |
| 100 | `(-4,4)` | `(-1,2)` |
| 101 | `(-1,1)` | `(0,0)` |

**两组 offsets 均为保留的固定历史参数，不重新拟合**。它们不是一组参数重复使用，
不从 topset 重新估计，也不是把 2D summary 移位、但 per-image residual 留在原坐标。

排名流程：

1. 四天全 ensemble 各自产生变换后的 pixel EVPA、I ordinary mean。
2. 每张 EVPA 与 I mean 使用 rank offsets 移位，summary 与 residual 均在同一对齐
   坐标计算。
3. 构造四天共同 ranking mask，要求每一天 std 有限且 **严格正**；采用下述强度、
   radius、可选人工排除规则。
4. 每个共同 identity 的 score 为四天、共同像素的 `sum((residual/std)**2)`。
   无额外 I 权重；residual 使用下述同一 ensemble branch。
5. `np.argsort(score, kind="stable")` 升序，平分按 canonical identity 顺序。
   取前 `topset_size`，输出四天严格对应的相同 identities。
6. 从 raw 输入重新构造 topset，使用 **final offsets** relocation，产生最终 summary。

ranking mask 为空、score 非有限、metadata 不匹配、topset 超长等都显式失败，不改用
历史 selection 或旧缓存。

## Ensemble mean / spread / covariance

EVPA 单位为弧度，周期为 π。默认且当前唯一支持的 `centering="ensemble"`：

```python
center = angle(mean(exp(2j * theta), axis=0)) / 2
unwrapped = center + circular_difference(theta, center)
unwrapped_mean = mean(unwrapped, axis=0)
mean = circular_difference(unwrapped_mean, 0)  # [-pi/2, pi/2)
std = std(unwrapped, axis=0, ddof=0)
residual = unwrapped - unwrapped_mean
concentration = abs(mean(exp(2j * theta), axis=0))
circular_std = sqrt(-2 * log(concentration)) / 2
```

这里 mean 是 circular center 附近 **ordinary mean**，不是直接 circular mean，也不是
Gaussian histogram 拟合。std 为 ordinary population spread，不是 circular_std。
concentration 理论范围 `[0,1]`；完全抵消时 circular_std 可为 `inf`，诊断不隐藏它。
`report.json` 的聚合统计只对有限值求 min/median/max，无有限值时为 null；原数组
完整保存诊断。精确零 resultant 时 center 使用 NumPy 的确定性 angle 约定，不宣称
该点具有物理上唯一的圆周中心，可据 concentration 在后续增加明确的筛选。

零偏振 `Q==U==0` 的 EVPA 为 NaN，不把 `atan2(0,0)` 的零当有效测量。
summary `count` 记录每个 feature 的有限样本数。为防止后续在不同特征使用不同
ensemble，任何 count 小于 N 的 feature，其 `mean/std` 均标 NaN；辅助诊断仍可描述
可用样本。2D mask 因此排除任一样本缺失/零偏振的像素。最终 mask 本身不要求
std>0；需要正定噪声的 inference 还必须显式拒绝零 variance 或采用有记录的策略。

`ensemble_residuals(theta, summary)` 使用保存的 center/unwrapped_mean。不能对这些
residual 再做 circular wrap，也不能用某个 batch 的 center 替代全 ensemble center。
`ensemble_covariance(theta, summary, batch_size=...)` 为相同分支 residual 的
`R.T @ R / N`，不使用 pairwise deletion；不完整 feature 对应整行整列 NaN。
每个日期的 `summary.npz` 以 `covariance_1d` 保存该 covariance，其 diagonal
与 `evpa_std_1d`² 一致。
不默认生成全部 4096 像素的巨大 2D covariance；父侧可用选中像素和保存的 center/
unwrapped_mean 构造所需 covariance。

**observations 集成：** `products.build_observations` 读取保存的 `center` 与
`unwrapped_mean`，严格按上述分支构造 residual，不围绕输出 mean 重新 unwrap，
不新增 resample。宽分布越过 branch 时，重新居中会改变 covariance，不能视为
舍入误差。函数继续从 mmap ensemble 提取选中的 features 计算 covariance，并
检查 residual mean 接近零；不会改写输入或复制旧 NPY。

读取按 day 支持两种布局：优先使用 `{date}/summary.npz` 和同目录
`stokes.npy` / `evpa-1d.npy`；没有该 summary 时使用旧
`summary_{1d,2d}.npz`（四天 stacked named keys）和
`stokes_{day}.npy` / `evpa_1d_{day}.npy`，允许逐日迁移。canonical summary
存在但 ensemble 缺失或 key 错误时显式失败，不静默换用旧 ensemble。
2D 的 `mask.npy` 在两种布局中都位于 EHT 顶层。
`load_templates` 与 `load_observations` 保持现有 named-key schema；父侧传入
`raptor/templates/{dimension}/templates.npz` 与
`observations/{dimension}/observations.npz` 的具体路径。

## 2D mask 与 1D aperture

2D mask 是四天 final topset 的共同条件：

- 每一天 mean I、EVPA std 有限；mean I **严格大于**各自当日 peak I 的 0.1 倍；
- 像素中心的 radius **≥2 μas**；网格中心 `(H-1)/2`，pixel scale `fov/H`；
- 默认排除所有列 `[0:5]`；
- 默认排除矩形行 `[48:55]`、列 `[10:19]`。

`--no-historical-mask` 关闭两项人工排除，但保留 radius 和 I 阈值。
`--exclude-columns` 配置列数；Python `exclusion_rectangles` 的元素为
`(row_start,row_end,col_start,col_end)`，half-open、非负有序、截到有效 image bounds，
小图不会发生负索引或人为缩放矩形。默认排除前 5 列会使 H≤5 的 mask 为空；tiny
测试应显式禁用人工排除或配置 bounds，不自动改变科学参数。

1D 直接在 **final relocated topset Stokes** 上调用父侧 coherent aperture integral：

- `weighting="intensity"`（默认）：`0.5 * arg ∫ I(Q+iU)`，EHT-VII 的 QI/UI 定义；
- `weighting="flux"`：`0.5 * arg ∫ (Q+iU)`，显式可选的 polarized-flux 定义。

这两种权重是不同统计量，不根据 Q/U 名字自动猜测 fractional/flux。integrator 的
zero integrated polarization 是 NaN。空间 NaN padding 在积分时作为零贡献；保存的
Stokes 保持 NaN padding，不污染其他像素。

默认 **1D 不继承 2D mask**，符合旧流程的 mask 分离方式，包括不继承 radius≥2 和
人工排除。`--inherit-2d-mask` 才会在积分之前将共同 mask 外所有 Stokes 置零。
最终每张 1D EVPA 做 `np.roll(..., -90, axis=-1)`，然后才做 ensemble summary。
该 angle roll 与空间 relocation 完全不同，必须保留。

## 产品契约

生产 H=64、K=1572、A=360；selection 四天轴顺序始终为 `(95,96,100,101)`。
浮点为 float64；indices/count 为 int64；文件均无 Python object arrays。
以下路径均相对于传入的 `output_dir`（通常为 `run/eht`），`{date}` 为上列四个日期。
不再新写顶层 `summary_1d.npz` / `summary_2d.npz` 或旧 flat ensemble 文件。

| 文件 | 内容 |
| --- | --- |
| `{date}/stokes.npy` | `(K,H,H,4)`，通道 I,Q,U,V；`open_memmap` 写入，可 mmap 读取 |
| `{date}/evpa-1d.npy` | `(K,A)`，积分后已 roll 的 EVPA radians；mmap |
| `{date}/summary.npz` | 当日 1D/2D named summaries、mean I、1D covariance，keys 见下表；无 day 轴 |
| `mask.npy` | `(H,H)` bool，最终 2D 共同 mask |
| `selection.npz` | 见下表；可 `allow_pickle=False` 读取 |
| `report.json` | 参数、fixed offsets、输入位置、变换约定、mask 像素数、score 范围、每 day 的 spread/concentration 聚合诊断 |

每个 `{date}/summary.npz` 含完整 16 个 named keys：

| keys | shape / 含义 |
| --- | --- |
| `evpa_mean_2d`, `evpa_std_2d`, `evpa_center_2d`, `evpa_unwrapped_mean_2d`, `evpa_concentration_2d`, `evpa_circular_std_2d`, `evpa_count_2d` | 各 `(H,H)`，以上 ensemble 统计 |
| `evpa_mean_1d`, `evpa_std_1d`, `evpa_center_1d`, `evpa_unwrapped_mean_1d`, `evpa_concentration_1d`, `evpa_circular_std_1d`, `evpa_count_1d` | 各 `(A,)`，同一统计定义 |
| `intensity_mean` | `(H,H)`，final relocated topset 的 ordinary mean I，保留 NaN borders |
| `covariance_1d` | `(A,A)`，保存的 full-ensemble branch population covariance |

四天 2D 小型 summaries 暂存内存，在最终 mask 确定、当日 1D summary 生成后合写
当日 NPZ；Stokes/EVPA 大数组仍逐批写入 mmap-capable NPY。

`selection.npz`：

| key | 含义 |
| --- | --- |
| `fields (7,)` | identity FIELDS 名称 |
| `identity (N,)` | canonical sorted identity，JSON string |
| `selected_identity (K,)` | 选中 identity，顺序与输出 Stokes/EVPA 相同 |
| `sorted_indices (N,)` | 全体 score 排序；index 指向 canonical identity |
| `selected_indices (K,)` | `sorted_indices[:K]` |
| `raw_sorted_indices (4,N)` | 各 day canonical identity → raw index |
| `selected_raw_indices (4,K)` | 各 day 选中 identity → raw index |
| `score (N,)` | canonical identity 顺序的排名 score |
| `rank_mask (H,H)` | 用于 ranking 的共同 bool mask，不等同于 final mask |
| `days (4,)` | 四天轴标签 |

## 内存与验证范围

按 day 读可信 pickle，变换/积分按 `batch_size=64` 分块；EVPA summary 按
`pixel_batch_size=128` feature 分块，避免一次构造 N×H² 复数 phase 与 residual。
仅存储一日 raw Stokes 加小批次数组；四天 ranking angle 写到 output 内临时 mmap，
评分后逐一删除。final angle mmap 每 day summary 后删除。没有四天 raw ensemble
常驻 RAM，也没有全 4800 图的四通道 padded 大副本。

对标准 float64、4800×64²、三通道 raw pickle，其数组约 450 MiB（有 V 则约 600
MiB）；单个 angle mmap 约 150 MiB。该布局以生产输入峰值 <2 GB 为设计目标，但
pickle 内部对象和 OS 文件缓存使精确 RSS 依赖环境；本文不声称已进行真实数据峰值
测量。父侧可后台测量；此任务仅运行 synthetic。临时磁盘最多约四天 angle
mmap（600 MiB），另需最终 Stokes、EVPA、summary 和 transaction 备份空间。

```bash
pixi run python -m pytest tests/test_pipeline_eht.py tests/test_pipeline_products.py -q
```

测试覆盖：translation 方向与 NaN 边界、I/Q/U 变换符号、metadata set 不匹配/重复/
缺字段/非有限值、跨 ±π/2 的 mean branch、宽分布 ordinary mean、population
covariance diagonal、zero polarization、不完整 features、历史及 tiny mask bounds、
四天 mini 全流程对齐 score 与 selection、固定两组 offsets、mmap 输出、-90 roll、
显式 flux/mask inheritance、拒绝未授权 pickle，以及真实父侧 aperture 的 tiny CLI
流程。新增覆盖 canonical/legacy pickle 输入、逐日期目录与全部 named summary keys、
不写旧顶层 summary、两种 observations 维度的新/旧/混合布局、保存 branch 下的
mean/std/covariance/pair 契约、只读 mmap ensemble 与不改写旧输入。
不会访问真实 4800 图数据、执行 notebook 或进行昂贵历史对照。
