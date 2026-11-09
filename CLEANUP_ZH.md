# 流水线整理说明（2026-11-06，分支 `cleanup/experiment-config`）

本次整理的目标：让物理学家能逐项审核整个流水线。所有科学参数集中到仓库根目录的
**一个**参数文件 `experiment.yaml`，每个数值只写一次；删除旧的“空假设（null-only）”
工作流和各类旁支研究脚本；删除只用于“检查副本是否一致”的代码；脚本按所属部分（raptor_ / eht_ / analysis_）
重新命名。旧产品的兼容路径全部删除（按用户决定），旧的 RAPTOR 作业通过一个**不属于
流水线**的一次性迁移脚本转换为新格式，然后用新代码从头到尾重算并与旧结果逐数比较。

两个互相独立的输入，以及把二者结合起来的分析（RAPTOR 模拟与 EHT 成像之间没有先后关系，
可以按任意顺序或同时进行；只有分析同时需要二者的产品）：

* **RAPTOR 模拟**（输入之一）：带轴子场的 GRMHD 流的偏振射线追踪图像；
* **EHT 成像**（输入之一）：按 EHT Paper VII 对 2017 年 M87* uvfits 数据做固定 I 的偏振成像、
  D 项求解，以及每天 1000 次 D 项抽样（共 4000 幅图像）；
* **分析**（结合二者）：模拟图像预处理（19 µas 模糊、旋转 108°、重采样到 64²/120 µas）、逐快照
   一阶谐波 EVPA 响应拟合、141 个快照上的模板、由无轴子快照对给出的背景 σ、EHT 系综上的
   特征选择、高斯似然、对相位边缘化的贝叶斯 c 上限。分析只遵循 EHT Paper VII：EHT 图像按重建
   结果直接使用（Paper VII 第 5.2 节的居中在本复现中等于零平移，即原来的 `native`）；
   原有的“质心居中（centroid）”备选方案已按用户决定从代码中整体删除。

---

## 删除了什么

### 文件

| 类别 | 删除 | 说明 |
|---|---|---|
| `configs/` | 全部 11 个 YAML | 由 `experiment.yaml` 取代 |
| `scripts/` | 41 个旧脚本 + `_stage.py`、`start_eht_120_templates.py` | 旧 null 工作流（run_all、process_raptor、process_eht、estimate_background、build_observations、calculate_limits、plot_* 等）、Apr11/fiducial 重建、joint-day 与背景再评估研究、迁移/传输诊断等 |
| `eht_axion/pipeline/` | convergence, diagnostics, migrate, transport, eht, inference, figures, workflow, config | 旧工作流各阶段；`workflow.py` 的 RAPTOR 分支由 `scripts/analysis_templates.py` 取代 |
| `eht_axion/research/` | background_validation, joint_days, state_ensemble, time_sampling | 已删除研究脚本的支撑模块 |
| `eht_axion/producer/config.py` | 整个文件 | 生产配置校验并入 `eht_axion/experiment.py`；硬编码副本（GRID、alpha 表、spin、build 默认值）删除 |
| `tests/` | 40 个测试文件 | 只测试已删除代码的测试；混合测试只删除相关部分 |

保留模块中删除的死代码（静态可达性分析：从全部保留脚本出发，沿函数名引用递归求可达的
`eht_axion` 顶层定义，不可达者删除）：`limits.calculate_from_products / calculate_limits` 及其 legacy 辅助函数
（只被旧 inference 使用；`research/observed.py` 只用到 `LimitOptions` 与 `_validate_request`）、
`observations.load_1d/2d_observations / select_independent_1d_indices`、
`model.load_fit_parameters`（历史 pickle 读入）、`statistics` 中三个未用函数、
`polarization.circular_window / integrated_evpa / radial_integration_matrix`、
`polarimetry` 中四个未用函数、`angles.circular_mean / unwrap_around`、
`eht_official` 中三个别名、`eht_optimizer.solve_fixed_p`、
`products.build_observations / load_observations`、`layout` 的 legacy 读入与 EHT pickle 路径、
`common` 中的 `stage.json`/“imported product” 兼容逻辑，以及质心居中所需的
`ensemble.centroid_offset / translate_ensemble`、`polarimetry.shift_stokes`。
另外删除了库函数签名里**隐藏的科学默认值**（`LimitOptions`、`evpa_templates`、`omega`、
`cloud_amplitude_scaling`、`calculate_observed` 的 alpha 网格、常量 `ALPHAS / AXION_NORM /
TIME_CONVERSION`）——它们是 YAML 数值的副本；现在所有值都由调用者从 `experiment.yaml` 传入。

### 行数（`.py` 与 `.yaml`）

| 目录 | 整理前（main 5b1807b） | 整理后 | 变化 |
|---|---|---|---|
| `eht_axion/` | 42 个文件，10 055 行 | 29 个文件，5 230 行 | −48 % |
| `scripts/` | 60 个文件，9 578 行 | 17 个文件，3 018 行 | −68 % |
| `tests/` | 75 个文件，13 036 行 | 37 个文件，6 027 行 | −54 % |
| `configs/` → `experiment.yaml` | 11 个文件，886 行 | 1 个文件，133 行 | −85 % |

测试：整理前 1715 个（eht-rebuild 环境），整理后 659 个，全部通过；默认环境 648 通过、
11 个因缺少 ehtim/astropy 跳过。`ruff`（E,F,W，行宽 110）无告警（main 上有 49 条）。

---

## `experiment.yaml` 的结构

规则：每个脚本读取 `shared` 和它自己的段；分析脚本另外读取它所需要的 raptor / eht_imaging
中的值（例如相机宽度、轴子归一化、抽样数）。代码中没有段与段之间的先后规则。

| 段 | 一句话说明 |
|---|---|
| `shared` | 不止一个部分用到的值：M87* 质量与距离、快照标签步长、alpha 网格、主量子数 n，以及公共天空网格（64 像素、120 µas）和 19 µas 匹配模糊（EHT 的 display19 产品与模拟图像用的是**同一个**物理量：天空平面上对 I/Q/U/V 施加的圆高斯核 FWHM，已核对 ehtim `blur_circ` 与 `pipeline.raptor.preprocess_frame`） |
| `raptor` | RAPTOR 模拟（独立输入）：相机宽度 46.2261 rg、build 选项、等离子体参数、快照区间、响应/背景分辨率与轴子归一化、径向函数表及 spin、alpha→表的对应、BHAC 网格 |
| `eht_imaging` | EHT 成像（独立输入）：原 `eht-official.yaml` 去掉移入 `shared` 的 npix/fov/display_blur 以及四个机器路径；保留 uvfits 文件名 + sha256（数据身份）、权重、轮数、重启方案、D 项统计、抽样数与种子 |
| `analysis` | 分析（结合两个输入）：旋转角、孔径、特征选择规则（强度/偏振阈值 0.1/0.2、像素步长 5、10° 扇区间隔与 4° 偏移、最小半径 2 rg、mask 模式、resultant ≥ 0.7）、云振幅 b(α) 参数、背景时滞与方法、推断网格/先验/时间 |

设计决定（均为“每个值只写一次”）：

* `cloud.reference_norm` 与 `producer.response.reference_norm` 是同一个量（模板振幅以 RAPTOR
  轴子归一化为单位）→ 只写在 `raptor.response.reference_norm`。
* `background.source_resolution` 与 `producer.background.resolution` 相同 → 只写在 `raptor.background`。
* `background.alpha: 0.3` 只是无轴子背景帧的目录名 → 变成 `pipeline/layout.py` 的常量
  `BACKGROUND_DIRECTORY_ALPHA`（不是科学参数）。
* RAPTOR 的 16 个相位（k·π/8）由 RAPTOR 内核决定 → `producer/plan.py` 的 `PHASE_COUNT`，不在 YAML 中。
* 网格/先验检查（phase 128 & 240 点；phase 128 & 260 点 & coupling_max 1e4）不再是 YAML 副本，
  而是 `analysis_check_grid.py` 的命令行覆盖 `--phase-count/--coupling-points/--coupling-max`，
  覆盖值记录在该检查的 manifest 里。
* 不再有居中开关：Paper VII 第 5.2 节的居中在本复现中就是零平移，`analysis_limits.py` 固定使用
  未平移（native）的 EHT 系综（常量 `CENTERING = 'native'`，仍记录在结果的 report/manifest 中）；
  `eht_summarize.py` 不再生成 icentroid 产品；`analysis_summarize.py` 以 native 上限为主结果
  （参数 `--limits`，不再需要质心输入）。
* 机器路径不进 YAML：RAPTOR 路径在 `.env.raptor`，EHT 成像的 uvfits 目录、输入收据、馈源旋转表、
  官方脚本在 `.env`（`EHT_IMAGING_*`）或命令行；它们仍按原键名记录进成像计划的 `parameters`。
  uvfits 文件名与 sha256 属于数据身份，留在 `eht_imaging.inputs`。
* 原先 `eht_sample`（批量 D 项蒙特卡洛）与 `eht_summarize`（系综汇总）中硬编码的抽样数、种子、
  fsys、D 项统计、站点、网格、显示模糊、孔径半径/宽度、resultant 阈值，现在都从 `experiment.yaml` 读取。
* 名称 `shared` 保留（它同时包含天空网格与模糊核，不只是物理常数）。

---

## 运行

RAPTOR 模拟与 EHT 成像互不依赖，可以任意顺序或同时运行；分析需要二者的产品。

```bash
pixi install --locked && pixi run check
# RAPTOR 模拟（环境 raptor-build；路径来自 .env.raptor）
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json
pixi run -e raptor-build raptor-generate
pixi run -e raptor-build raptor-merge --plan plan.json --job-root RUN_A --job-root RUN_B --output-root MERGED
# EHT 成像（环境 eht-rebuild；uvfits 位置来自 .env）
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC
pixi run -e eht-research eht-summarize --dataset-root MC_COPY --output-root ENSEMBLE
# 分析（默认环境，需要上面两类产品）
pixi run analysis-templates --raptor-input-root MERGED/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root GRID --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root PRIOR --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance --source LIMITS --templates TEMPLATES --output-root COV
pixi run analysis-summarize --limits LIMITS --gridcheck GRID --priorcheck PRIOR --covcheck COV --output-root SUMMARY
```

脚本一览（按所属部分）：`raptor_build / raptor_generate / raptor_merge`，
`eht_reconstruct / eht_refine / eht_refine_d / eht_sample / eht_summarize`，
`analysis_templates / analysis_limits / analysis_check_grid / analysis_check_covariance / analysis_summarize`，
以及 `verify_products`。`analysis_templates.py` 取代了原来的 `start_eht_120_templates.py` + `_stage.py` +
`workflow.run_stage` 的 RAPTOR 分支，产生同名阶段与相同参数（`raptor:preprocess`、`raptor:fit`、
`raptor:summarize`、`raptor:background-preprocess`、`raptor:background-unmasked`）。

溯源（provenance）保持不变：每个阶段目录有 `manifest.json`，记录参数、输入/输出 SHA-256、代码哈希、
运行环境；`pixi run verify DIR --recursive` 重新检查整条血统。只删除了纯“副本一致性”检查。

---

## 不兼容旧产品；作业 id 与一次性迁移

按用户决定，整理后的代码**不含**任何接受旧格式的路径：不接受旧的 `generated.json` 身份、旧作业 id、
旧配置布局、`stage.json` 或 “imported” 产品。

* **作业 id**：仍是 SHA-256，但只由影响 RAPTOR 结果的量计算：快照（文件名 + BHAC 尾部元数据）、
  19 个 model.in 数值（不含表格路径；含分辨率、alpha、omega、norm、相位、等离子体与相机参数）、
  径向表内容 sha256 与 spin、网格、build 设置（`producer/plan.py::task_identity`）。分析、天空网格、
  EHT 成像参数不进入 id。生产计划只保存“射线追踪子集”（`shared` 的物理量 + `raptor` 段），
  `generated.json` 为 schema v2。
* **迁移**：`reports/cleanup-20261106/migration/migrate_jobs.py`（不在 `scripts/` 或 `eht_axion/` 中，
  不属于流水线）。它用新代码生成新计划，把两个旧作业根（`/home/canoziia/PHYS-production/results/m65-211`
  的 18 234 个作业、`/home/canoziia/GRMHD_DATA/production/pilot-m65-t3000-res400-t16` 的 96 个作业）
  逐一对应到新任务，并对每个作业断言：旧 task.json 除 id 外与新代码会记录的完全相同、model.in 字节相同、
  grid 与 build 相同、源快照与径向表 sha256 与旧记录相同、旧输出收据成功。然后在
  `/home/canoziia/GRMHD_DATA/production/m65-211-newids-20261106/` 下按新代码的布局写出
  `<根>/{build, jobs/<新 id>, config.yaml, requests/}` 与顶层 `plan.json`、`MIGRATION.json`：
  build 文件和作业输出（HDF5、spectrum、model.in、grid_mks.in、日志、success.json）都是**硬链接**
  （同一 inode，原文件只读不改），task.json、output.json、manifest.json 由新代码写出。
  两个旧根的 RAPTOR 二进制字节不同（不同构建环境），因此新根下仍分两个子根，各自带自己的 build；
  build 字节与 sha256 记录在 `MIGRATION.json`。attempts/ 尝试收据留在原根中。

---

## 等价性证据

所有证据 JSON 同时保存在 `reports/cleanup-20261106/` 与
`/home/canoziia/GRMHD_DATA/analysis/cleanup-verify-opus-20261106/evidence/`。

### 配置（`config-equivalence.json`）

脚本 `check_config.py`：把 main（5b1807b）`configs/` 中被活跃链路读取的每个值、以及 main 活跃脚本中的
字面量（特征选择阈值、批量抽样常数、系综汇总常数等）与新加载器逐项对比：**101 项全部相等**。
用 main 代码与新代码分别计算的派生量——分析几何（模糊 σ、视场、旋转、孔径半径）、`LimitOptions`
（主结果、网格检查、先验检查三套）、成像参数集（含键顺序）——**全部逐位相同**。文件中还列出旧配置之间
的副本关系（production/analysis/gridcheck/priorcheck/eht-official 中重复的值）以及删除的死字段及原因：
`eht.topset_size / historical_mask / inherit_2d_mask_in_1d / covariance_ridge`、整个 `eht_input` 块
（只服务于已删除的旧 pickle/null 工作流）、Rg 模式字段、`response.phase_count`、`background.alpha` 等。

### RAPTOR 计划（`raptor-plan-equivalence.json`）

`check_raptor_plan.py`：用真实输入（流体 `/home/canoziia/GRMHD_DATA/a094-high`，径向表取自归档计划
`/home/canoziia/PHYS-production/control/full-plan.json`），`hash_inputs=False`，分别用 main 代码和新代码
生成计划：**18 330 个任务一一对应**（141 快照 × (2 个基线 + 8 alpha × 16 相位)），每个任务的
`parameters`（19 个 model.in 值）、源路径、时间、kind、分辨率、alpha、omega、norm、相位、输出槽、
径向表 id/路径/spin **全部相同**；与归档计划逐任务比较也全部相同；任务顺序、快照尾部元数据、径向表元数据、
build 设置、网格、setup 参数集合、资源估计均相同。唯一变化是任务 id（按设计）。

### EHT 成像（`eht-plan-equivalence.json`，`eht-sample-equivalence.json`）

`check_eht_plan.py`（eht-rebuild 环境，只用 `--plan`，不做成像）：对 `eht-paper-I-{95,96,100,101}-v1`
按记录的命令行重放，新计划与 main 计划、与记录的 `plan.json` 相比，除 `cli`（配置路径、四个位置参数）与
代码文件输入外**全部相同**（parameters、raw、draws、label、scope、display，数据输入与 ehtim 库文件的路径
和 sha256）。`eht-paper-PD-*-v1` 和 MC 试点无法由 main 或新代码重放（其 I 阶段 manifest 列出的 main
检出代码在记录后已改动，血统校验失败——两边相同）；对它们比较了参数集：**与记录完全相同**。
`check_eht_sample.py`：从 `experiment.yaml` 读取的批量 D 项蒙特卡洛设置重放 `eht-uncertainty-remote-v1`：
4×1000 次 D 项抽样与任务列表，新代码与 main 代码**逐位相同**；与远程主机记录的结果相比最大差
5.6×10⁻¹⁷（BLAS/libm 末位差异），分布元数据与 bundle sha256 相同。

### 迁移（`migration-receipt.json`；完整映射在新根的 `MIGRATION.json`）

18 330/18 330 个旧作业一一匹配，0 个失败；每个作业通过 6 项检查（链接字节等于旧记录、同 inode、
output.json 除 build 路径外相同、源与径向表 sha256 等于旧记录、输出文件集合相同、task.json 除 id 外相同）。
耗时 26.6 分钟（8 核），新根占用约 1 GB 实际磁盘（13 GB 为硬链接）。

### 新代码端到端重算（产品比较 `products-*.json`）

在 `/home/canoziia/GRMHD_DATA/analysis/cleanup-verify-opus-20261106/` 下用新代码依次运行
（`run_verification.sh`，≤8 核，耗时见括号）：

| 步骤 | 新输出 | 与原产品比较 | 结果 |
|---|---|---|---|
| 合并 raptor_merge（61 min） | `raptor-merged` | `/home/canoziia/GRMHD_DATA/production/m65-211-merged/data` | 23 个 Stokes 立方体**全部逐位相同**；JSON 中只有 generated.json（schema v2）、sources.json（新 id）、路径字符串不同 |
| 模板/背景 analysis_templates（22 min） | `templates` | `eht-120uas-observed-templates-v1` | 29 个数组**逐位相同**；42 个 JSON/manifest 数值全同（只差输入路径字符串）；阶段名与参数相同 |
| EHT 系综 eht_summarize（14 min） | `eht-ensemble` | `eht-mc-summary-v1` | 89 个数组**逐位相同**；report 与 manifest config 相同；consumed-inputs 只差代码文件与新增的 experiment.yaml |
| 上限 analysis_limits（22 min，主要是全链血统校验） | `limits-native` | `eht-120uas-observed-limits-native-v1` | 16 个 NPZ **逐位相同**；limits.csv、全部 observed/null JSON 数值全同；PNG 字节相同；只差 `scientific_profile`（新 schema）与背景源路径 |
| 网格检查（10 s） | `native-gridcheck` | `eht-120uas-native-gridcheck-v1` | 12 个数组逐位相同，JSON/CSV 数值全同，PNG 相同 |
| 先验检查（5 s） | `native-priorcheck` | `eht-120uas-native-priorcheck-v1` | 同上 |
| 协方差检查（1 s） | `native-covcheck` | `eht-120uas-native-covcheck-v1` | CSV 与 report 数值全同（只差源路径） |
| 汇总 analysis_summarize（6 s，最终代码） | `results-summary` | 原 native 上限 | `summary-consistency.json`：96 个上限值与原 native limits.csv 相同，网格差最大值 0.004026 dex 与原检查相同 |

“逐位相同”= 形状、dtype 相同且 `np.array_equal`（NaN 视为相等）；容差 1e-12 实际上未被用到。

### 测试与静态检查

`pytest`（eht-rebuild 环境）659 通过；默认环境 648 通过 / 11 跳过；`ruff check --select E,F,W`
（行宽 110）无告警。

---

## 未验证或遗留的问题

* **按用户决定跳过的重算**：质心（centroid）上限、质心网格/先验/协方差检查以及基于质心的旧汇总
  `eht-120uas-results-summary-v1` 未重算——质心方案已整体删除，分析只遵循 Paper VII。
* **在中途改动代码之后未重跑、靠构造保证等价的部分**：
  - `limits-native` 及三个 native 检查由提交 afb7118 的代码产生；之后删除质心开关只删掉了
    `centering == 'centroid'` 分支和 `--centering` 参数，native 路径（目录 `native/display19`、
    阈值、步长、扇区、似然）代码与数值完全不变；
  - `eht-ensemble` 由删除 icentroid 之前的 `eht_summarize.py` 产生；删除的分支只生成 icentroid 目录与
    `registrations` 报告字段，native 汇总调用的输入与参数不变，因此 native 文件逐位不变；
  - 新的检查脚本会拒绝 schema 不同的源（例如由旧 `experiment.yaml`（含 `centering` 键）产生的
    `limits-native`）——这是“不兼容旧产品”的设计，因此最终代码下需要从 `analysis_limits` 起重跑才能
    继续做新检查（约 25 分钟，主要是全链血统哈希）。
* **未执行的科学计算**：没有运行 RAPTOR 射线追踪，也没有运行 EHT 重建/成像（只做 `--plan` 重放与
  D 项抽样重放）；`eht_refine.py / eht_refine_d.py` 只有单元测试，未在真实数据上验证。
* **原文件只读**：迁移只读取并硬链接原作业文件；内容与 mtime 不变，但硬链接会使 inode 的 ctime 与
  链接计数改变（硬链接的固有行为）。原 `attempts/` 收据未迁移。
* **运行中的两次中断**：第一次 `eht_summarize` 因一个变量名遮蔽错误失败（已修复并加了合成测试）；
  第一次 `limits-native` 运行时 tmux 服务进程意外退出（原因不明，内存充足），重新启动后正常完成。
* **库中保留的方法开关**：`LimitOptions` 仍保留 `legacy`、`std_covariance` 等方法字段（它们写进结果 JSON，
  删除会改变记录格式）；活跃链路固定使用非 legacy 路径。`read_snapshot_footer` 中流体 spin 0.9375 的
  数据检查仍是代码常量。
* **文档**：`docs/` 下的历史文档按要求未改，仍描述旧的 `configs/` 布局；README 已更新并注明。
* **合并耗时**：`raptor_merge` 对 18 330 个作业逐一做递归血统校验，约 1 小时（单进程）；如需加速可并行，
  但本次未改动这部分逻辑。
