# EHT_ZC 代码分类与用途清单

## 说明与范围

本文记录整理代码前的用途分类，依据当前文件内容和实际调用关系，而不是仅按目录名分类。

- 审查基线：`74ebae1`，分支 `data/analysis-inputs-draft`。
- 范围：Git 跟踪的 **41 个 `.py`、1 个 `.sh`、16 个 `.ipynb`**，共 58 个代码文件。
- 本步只分类，不移动、删除代码，不改变算法。
- “当前运行”只说明它被使用，**不代表物理定义、公式和实现已经正确**。
- 现有 `eht_axion/` 和 `scripts/` 是后来整理新增的实现；旧 Notebook 是此前算法的重要出处。计划以此前代码为基础逐步整理时，需要保留二者对应关系。
- 不能把旧文件、单独运行的诊断工具、审计工具一概叫作死代码。

### 状态标记

| 状态 | 含义 |
|---|---|
| 主流程 | 当前从共享 final-inference 输入计算六条曲线/画图的入口及依赖 |
| 上游再生 | 从 RAPTOR/EHT 图像输出生成最终统计输入时需要的代码 |
| 科学诊断 | 用于理解/验证物理效应及生成诊断图，不是六曲线的必要计算步骤 |
| 工具 | 数据管理、校验、测试、回归等辅助入口 |
| 历史实现 | 保存此前算法与结果；不是当前推荐入口，但可能包含尚未完全迁出的功能 |
| 暂不使用 | 当前运行链中不调用的实验、备份或空文件；暂不整理，不等于立即删除 |

## 一、五类用途

### 1. 数据处理

目标：把成像输出转换成统计输入，不在此阶段求最终耦合上限。

**RAPTOR 一侧：**

```text
Stokes IQUV
→ 旋转/模糊/重采样/裁切
→ 二维 EVPA
→ 径向及角窗口积分得到一维 EVPA
→ 对 16 个轴子相位拟合振幅 A、相位 φ
→ 估计 GRMHD 演化及轴子模板的散布/误差
```

**EHT 一侧：**

```text
四天各 4800 个候选图及超参数
→ 候选图排序与 topset 选择
→ 黑洞中心定位/平移
→ EVPA 代表值与散布
→ 空间 mask、像素/角向误差矩阵
```

输出：EVPA ensemble、均值/标准差、模板参数、mask、协方差或其他明确命名的误差矩阵。

### 2. 统计推断

```text
观测量 + 模板 + 误差模型
→ 观测时间与曝光响应
→ 耦合参数下的信号预测
→ likelihood
→ 初相位边缘化
→ posterior
→ 每个 α 的 95% 上限
```

当前组合：`std / cor / cor_grmhd` × `1d / 2d`；另区分零观测日差（null）与实际日差（observed）。

这里必须显式说明频率/时间单位、差分顺序、误差相关性、先验和参数网格。输出数值结果，不做论文排版。

### 3. 画图与科学诊断

- **结果绘图**：读取已算好的限制曲线并绘制，不在画图函数里重新做统计推断。
- **科学诊断**：信号分布、拟合质量、解析/RAPTOR 对照、等离子体开关、首次赤道穿越/透镜像/非赤道发射贡献等。

截图所指的“图 4 传输贡献拆解”位于 `GRMHD/main.ipynb` 的诊断部分，**不等同于** `Plot/fig4.ipynb`（后者主要画误差）。

### 4. 数据管理与复现工具

共享数据导出、分片、manifest、SHA-256、Git LFS、数据版本和来源记录。

本类不应暗中改变科学量的定义。共享目录为 `shared_data/`；其他本地中间数组仍按 `.gitignore` 排除。

### 5. 验证与测试

单元测试、合成信号解析答案、旧结果回归、输入完整性、参数网格收敛和消融检查。

测试通过只证明被覆盖的行为；不自动证明概率模型或物理定义成立。

### 依赖关系

```text
                     数据管理
                        │
                        ▼
输入成像结果 → 数据处理 → 统计推断 → 结果绘图
                  │           │
                  └── 科学诊断图

              验证与测试覆盖各阶段
```

## 二、当前最终六曲线的运行链

```text
shared_data/final-inference/draft/{1d,2d}/
  ↓
scripts/run_final_limits.py
  ↓
eht_axion/limits.py
  ├── observations.py  → angles.py
  ├── model.py
  └── statistics.py
  ↓
每个 mode/dimension 的 NPZ + JSON
  ↓
scripts/plot_final_limits.py
```

**当前链路限制：**runner 的默认输出与绘图脚本固定读取目录不一致；重算后不能直接假定绘图读的是新结果。这里只记分类，暂不修复。

共享 `analysis-inputs` → 全部统计输入的再生链，目前尚未形成统一 CLI：仍有旧路径、历史对照文件和 Notebook 阶段依赖。

## 三、按类别分配现有代码

### 第 1 类：数据处理

| 文件 | 主要职责 | 状态/备注 |
|---|---|---|
| `scripts/validate_topset_ranking.py` | 从 4800 图计算旧/对齐后排名，输出 corrected topset 顺序 | 上游再生 + 第5类审计；虽名为 validate，实际生产了后续重建输入 |
| `scripts/rebuild_corrected_topset.py` | 根据排名选 1572 图，旋转、平移、拟合日均值/散布并生成 mask | 上游再生 |
| `scripts/rebuild_1d_products.py` | EHT/GRMHD Stokes 转一维 EVPA，并做轴子相位余弦拟合 | 上游再生；读取 processed GRMHD，不是原始 RAPTOR cube 的统一入口 |
| `scripts/fit_rebuilt_1d_products.py` | 对一维 EVPA 和振幅/相位 ensemble 拟合统计参数 | 上游再生 |
| `eht_axion/observations.py` | 加载日均值/散布、选点、构造日差和观测误差矩阵 | 主流程；属于数据处理与推断的接口 |
| `eht_axion/polarization.py` | 生成径向积分算子、角窗口及 Stokes 合成 EVPA | 上游再生共用库 |
| `eht_axion/angles.py` | EVPA/相位的周期差、圆均值、展开 | 共用基础库；供第1/2/5类调用 |
| `GRMHD/alpha_0.1.py` | α=0.1、res=100 图像处理与 EVPA 保存 | 上游旧入口；与其余 alpha 脚本大段重复 |
| `GRMHD/alpha_0.15.py` | α=0.15、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.2.py` | α=0.2、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.25.py` | α=0.25、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.3.py` | α=0.3、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.35.py` | α=0.35、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.4.py` | α=0.4、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.44.py` | α=0.44、res=100 图像处理与 EVPA 保存 | 上游旧入口 |
| `GRMHD/alpha_0.3_high_res.py` | α=0.3、原始 res=400 图像处理 | 上游旧入口；原始输入400不等于最终输出400，默认仍重采样到64 |
| `EHT-Imaging/main.ipynb` | EHT 4800→topset、定位、EVPA 拟合、一维数据 | 历史主实现；算法出处，不是最终 runner 的依赖 |
| `GRMHD/main.ipynb` | 图像处理、余弦拟合、解析对照等 | 混合 Notebook，详见下一节；不可整个归为死代码 |
| `utils/libanaly.py` | 图像、圆统计、Gaussian拟合和绘图工具 | 混合共用库，详见下一节 |

### 第 2 类：统计推断

| 文件 | 主要职责 | 状态/备注 |
|---|---|---|
| `scripts/run_final_limits.py` | 选择维度/误差模式/观测假设/输入目录，执行并保存上限 | 主流程 CLI |
| `eht_axion/limits.py` | 联合两组日差 likelihood、耦合扫描、相位边缘化和上限 | 主流程计算 |
| `eht_axion/model.py` | α频率、云场缩放、曝光/日差模板及模型误差传播 | 主流程计算；待逐项复核时间公式 |
| `eht_axion/statistics.py` | Gaussian likelihood、posterior积分和分位点 | 主流程共用库；部分helper未被生产调用 |
| `eht_axion/__init__.py` | package及角度API导出 | 包结构，不是独立算法阶段 |
| `result.ipynb` | 旧数据处理 + 推断 + 绘图/实验 | 历史核心算法出处；部分上游量尚未完全迁出，详见下一节 |

### 第 3 类：画图与科学诊断

| 文件 | 主要职责 | 状态/备注 |
|---|---|---|
| `scripts/plot_final_limits.py` | 六曲线图、legacy对照和CSV/JSON表 | 主流程结果绘图；目前输入目录固定 |
| `Plot/fig1.ipynb` | 轴子径向函数随半径/α的图 | 科学诊断/论文图；使用径向CSV表 |
| `Plot/fig2.ipynb` | 二维轴子振幅、相位、EVPA和强度图 | 科学诊断/论文图；多处仍读旧路径 |
| `Plot/fig3.ipynb` | 一维信号、振幅/相位和观测EVPA分布 | 科学诊断/历史绘图；含786图方案数据 |
| `Plot/fig4.ipynb` | GRMHD、轴子、观测误差及直方图 | 科学诊断/混合绘图；不是截图中的传输拆解图4 |
| `Plot/eht-imaging.ipynb` | EHT图像及covariance-距离诊断 | 科学诊断 |
| `Plot/fig?.ipynb` | 四天图像/EVPA比较的探索性绘图 | 暂不使用；保留核对后再决定归档 |
| `GRMHD/GRMHD_cut.ipynb` | 光线轨迹分组、赤道面交点/场参考生成与对照 | 科学诊断 + 第1类辅助数据处理；最终六曲线不用 |
| `utils/raytracing.py` | Kerr光线追溯至源位置 | 科学诊断共用库；最终六曲线不用 |

图4专项数据目录：

```text
shared_data/analysis-inputs/draft/raptor-transport-diagnostics/
```

该旁路主要在 `GRMHD/main.ipynb`，涉及 `common/jonly/eq_jonly/fake`、赤道穿越坐标和解析场；不是最后六曲线的必要输入。

### 第 4 类：数据管理与复现工具

| 文件 | 主要职责 | 状态/备注 |
|---|---|---|
| `scripts/publish_shared_data.py` | 导出 final-inference 共享数据与manifest | 维护工具 |
| `scripts/publish_analysis_inputs.py` | 导出主RAPTOR数组、res400分片、EHT集合和传输诊断 | 维护工具 |
| `scripts/verify_shared_data.py` | 按manifest验证大小、SHA和文件可读性 | 合作者下载后工具；也属于第5类完整性验证 |

相关配置/文档不是科学计算代码：

- `.gitignore`：仅共享命名空间放行科学二进制，其余本地数据保持忽略。
- `.gitattributes`：共享命名空间的 Git LFS 规则。
- `pyproject.toml`：安装依赖、package发现、pytest/Ruff配置。
- `requirements.txt`：此前环境依赖快照。
- `shared_data/**/README.md`：数据用途、下载和发布流程。
- `shared_data/**/manifest.json`：shape/dtype/hash/来源等元数据。

### 第 5 类：验证与测试

| 文件 | 主要职责 | 状态/备注 |
|---|---|---|
| `scripts/inspect_intermediates.py` | 中间量结构、hash、有限值和Stokes均值一致性 | 审计工具 |
| `scripts/validate_fit_products.py` | 从源振幅/相位样本重拟合抽查参数 | 审计工具 |
| `scripts/validate_historical_selection.py` | 从四天图像签名反推旧topset索引 | 历史选择回归 |
| `scripts/verify_legacy_results.py` | 比较六条legacy重算值和旧保存值 | 数值回归；不证明公式正确 |
| `scripts/check_grid_convergence.py` | 60/120耦合网格对比 | 收敛诊断；不覆盖初相位网格收敛 |
| `scripts/run_ablation.py` | 逐项切换新旧统计实现并比较 | 消融诊断 |
| `scripts/verify_all.sh` | 串联测试、回归、六曲线和绘图 | 审计集成入口；依赖已有本地中间量，不是从共享analysis-inputs开始的完整pipeline |
| `tests/test_angles.py` | 周期差、圆均值、展开 | 单元测试 |
| `tests/test_legacy_utils.py` | CircleCalc、ImageS、Gaussian检验接口 | 单元测试 |
| `tests/test_limits.py` | legacy分位范围、模板shape与两组误差 | 单元测试；未验证完整观测时间解析答案 |
| `tests/test_model.py` | 相位/振幅插值 | 单元测试 |
| `tests/test_observations.py` | 1D固定offset选择 | 单元测试；未覆盖完整观测loader |
| `tests/test_polarization.py` | 窗口与积分接口 | 单元测试；未覆盖完整径向quadrature正确性 |
| `tests/test_statistics.py` | Gaussian、posterior及部分helper | 单元测试；部分helper未被最终算法调用 |

## 四、混合文件需要按块整理，不能整份判死

Notebook cell 编号按当前文件 **从0开始**，后续增删单元后可能变化；应同时使用函数名定位。

### `GRMHD/main.ipynb`

| 功能块 | 类别 | 后续处理建议 |
|---|---|---|
| cell 7–10：`getGRMHDDataAfterImageS`、主数组预处理 | 第1类数据处理 | 上游算法来源，应整理到可配置模块 |
| cell 15–19：DFT余弦拟合、时间序列/参数扫描 | 第1类数据处理 | 保留主时间序列拟合；旧parameter scan另标历史 |
| cell 29–37：径向函数、解析源位置和人工模板 | 第3类科学诊断，兼辅助处理 | 与最终限制分离 |
| cell 42–64：解析/RAPTOR比较、传输项贡献及彩色环图 | 第3类科学诊断 | 对应图4旁路，需要单独整理入口 |
| 之后的圈选、fake/jonly变体、偏振矢量、κ图等 | 第3类诊断/探索 | 分实验确认，暂不进入最终pipeline |
| 注释掉的旧fit/save片段 | 历史残片 | 确认无必要信息后归档，不作为生产代码 |

### `result.ipynb`

| 功能块 | 类别 | 后续处理建议 |
|---|---|---|
| GRMHD一天/8小时EVPA差与统一散布估计 | 第1类数据处理 | 需显式提取；新runner目前直接使用固定sigma数值 |
| cell 36：`get_A_phi0_params()` | 第1类数据处理 | **2D参数字典生成的旧算法来源**，不能整个Notebook判成无用 |
| cell 50：`get_A_phi0_1d_params()` | 第1类数据处理 | 旧1D参数生成，与新MLE实现需要明确选择 |
| `axion_evpa_diffs_*`、`get_final_cut*` | 第2类统计推断 | 旧公式出处，后续逐项对照新模块复核 |
| `plotting_script`、限制对比图 | 第3类结果绘图 | 与计算分离 |
| EXPERIMENT、786-topset残留及未定义变量探索单元 | 历史/暂不使用 | 暂不整理，不作为运行入口 |

### `utils/libanaly.py`

| 功能 | 类别 | 当前用途 |
|---|---|---|
| `CircleCalc`、角度居中 | 第1/2类共用基础 | 旧预处理、Gaussian拟合和测试使用 |
| `gaussian`、`fit_gaussian`、`fit_gaussian2` | 第1类拟合 | 上游EHT脚本实际主要用`fit_gaussian2`；其余函数不能只因在库里就视为主依赖 |
| `ImageS`正常模糊/旋转分支 | 第1类图像处理 | GRMHD alpha脚本与Notebook使用 |
| `Create_EVPA_Dataset` | 第1类积分 | 旧spline实现/诊断；新1D批处理主要用`eht_axion/polarization.py` |
| `EVPA_1d_from_coord*` | 第3类圈选诊断 | 非当前最终重建主分支 |
| `Plot_Background`、`plot_quivers`、`Plot_*` | 第3类绘图 | 诊断Notebook使用，最终runner不需要 |

## 五、历史备份、实验和空文件

以下不是当前六曲线的必要入口；先不改、不删，后续按需要归档：

| 文件 | 原用途 | 分类 |
|---|---|---|
| `result copy.ipynb` | 旧最终处理/推断/绘图备份 | 第1/2/3类历史实现 |
| `EHT-Imaging/main copy.ipynb` | 更早EHT处理备份，含旧import/定位实现 | 第1类历史实现 |
| `EHT-Imaging/main copy 2.ipynb` | EHT处理另一备份 | 第1类历史实现 |
| `GRMHD/Gaussian Blur.ipynb` | 内嵌旧ImageS类、模糊与相位拟合试验 | 第1/3类历史实验 |
| `GRMHD/rotate.ipynb` | HDF5读取、数组汇总、坐标/矢量检查 | 第1/3类历史上游来源；硬编码路径，不是可靠生产入口 |
| `test.ipynb` | 光追、kgeo对照、拟合与绘图实验 | 第3/5类探索，不能当pytest测试集 |
| `Plot/main.py` | 空文件 | 无实现，当前可确认无计算用途 |

### 当前可确认的未用/重复函数或常量

“未用”限于当前仓库生产调用，不代表外部API一定无人调用。

- `eht_axion/model.py::ALPHAS`：未使用，与`limits.py::ALPHAS`重复。
- `eht_axion/statistics.py::combine_independent_log_likelihoods`：仅测试调用；实际主算法直接`sum`。
- `eht_axion/statistics.py::marginalize_log_likelihood`：未见调用；实际主算法直接`logsumexp`。
- `eht_axion/angles.py::unwrap_around`：导出且测试，但当前生产链无调用。
- `select_independent_1d_indices(offset=None)`：备用自动选择分支；当前主入口固定offset=4。
- `_legacy_*` 和 legacy/ablation switches：正常新结果不用，但历史回归明确使用，**不列为死代码**。

## 六、数据、结果、文档与代码的边界

以下目录/文件不是算法代码，不按五类代码硬归属：

| 路径 | 类型 |
|---|---|
| `shared_data/final-inference/draft/` | 最终统计输入，Git LFS共享数据 |
| `shared_data/analysis-inputs/draft/` | RAPTOR输出、EHT4800图、图4专项输入，Git LFS共享数据 |
| `GRMHD/data_in/`、`GRMHD/data_out/` | 本地上游/中间科学数据 |
| `EHT-Imaging/data_in/`、`EHT-Imaging/data_out/`、`EHT-Imaging/result/` | 本地EHT输入/中间量及历史副本 |
| `audit/rebuilt_1d/`、`audit/rebuilt_2d/`、`audit/topset_work/` | 本地审计生成中间量，非正式数据命名空间 |
| `Results/`、`Results_nocondition/`、`new_static/`、`M_all/`、`Fit_result/` | 历史数值结果/图像，非运行入口 |
| `audit/results/`、`audit/figures/` | 本轮实现生成的数值/图像与回归产物 |
| `R.csv`、`R.500.csv`、`NormR.csv`、`NormR_0.99.csv` | 径向波函数科学数据；解析诊断/绘图使用 |
| `README.md`、`Introduce.md`、`EHT-Imaging/table.md` | 使用说明、历史方法说明和超参数分布记录 |
| `audit/REPORT.md` | 上轮审计记录；其中“已修复”判断仍需结合后续复核，不自动作为最终方法规范 |

## 七、下一步整理原则

1. 先确认原算法与观测量定义，再精简和修复，不能因新实现名为corrected就默认正确。
2. 优先整理：数据处理 → 统计推断 → 结果绘图，保证输入输出接口明确。
3. 科学诊断独立入口；不要让图4/fake/jonly等旁路成为六曲线的隐含依赖。
4. 历史Notebook保留为算法出处，把仍有用的功能逐块迁出；不要整份删除。
5. 数据管理/校验与科学计算分离。
6. 测试直接覆盖真正被调用的函数及物理解析答案，而不只测试未调用helper。
7. 当前已发现的时间相位、误差定义、积分定义和路径问题，在后续步骤逐项处理；**本分类文档不宣称这些问题已修复**。
