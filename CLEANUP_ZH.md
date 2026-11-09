# 流水线整理说明（2026-11-06，分支 `cleanup/experiment-config`）

本次整理的目标：让物理学家能逐项审核整个流水线。所有科学参数集中到仓库根目录的
**一个**参数文件 `experiment.yaml`，每个数值只写一次；删除旧的“空假设（null-only）”
工作流和各类旁支研究脚本；删除只用于“检查副本是否一致”的代码；脚本按三个过程的顺序
重新命名。旧产品的兼容路径全部删除（按用户决定），旧的 RAPTOR 作业通过一个**不属于
流水线**的一次性迁移脚本转换为新格式，然后用新代码从头到尾重算并与旧结果逐数比较。

三个过程：

1. **RAPTOR 模拟**：带轴子场的 GRMHD 流的偏振射线追踪图像；
2. **EHT 成像**：按 EHT Paper VII 对 2017 年 M87* uvfits 数据做固定 I 的偏振成像、
   D 项求解，以及每天 1000 次 D 项抽样（共 4000 幅图像）；
3. **分析**：模拟图像预处理（19 µas 模糊、旋转 108°、重采样到 64²/120 µas）、逐快照
   一阶谐波 EVPA 响应拟合、141 个快照上的模板、由无轴子快照对给出的背景 σ、EHT 系综上的
   特征选择、高斯似然、对相位边缘化的贝叶斯 c 上限。论文结果使用 `centering: native`
   （Paper VII 第 5.2 节的居中在本复现中等于零平移）。

---

## 1. 删除了什么

### 1.1 文件

| 类别 | 删除 | 说明 |
|---|---|---|
| `configs/` | 全部 11 个 YAML | 由 `experiment.yaml` 取代 |
| `scripts/` | 41 个旧脚本 + `_stage.py`、`start_eht_120_templates.py` | 旧 null 工作流（run_all、process_raptor、process_eht、estimate_background、build_observations、calculate_limits、plot_* 等）、Apr11/fiducial 重建、joint-day 与背景再评估研究、迁移/传输诊断等 |
| `eht_axion/pipeline/` | convergence, diagnostics, migrate, transport, eht, inference, figures, workflow, config | 旧工作流各阶段；`workflow.py` 的 RAPTOR 分支由 `scripts/analysis_templates.py` 取代 |
| `eht_axion/research/` | background_validation, joint_days, state_ensemble, time_sampling | 已删除研究脚本的支撑模块 |
| `eht_axion/producer/config.py` | 整个文件 | 生产配置校验并入 `eht_axion/experiment.py`；硬编码副本（GRID、alpha 表、spin、build 默认值）删除 |
| `tests/` | 40 个测试文件 | 只测试已删除代码的测试；混合测试只删除相关部分 |

保留模块中删除的死代码（静态可达性分析：从全部保留脚本出发，`reports/` 中无此分析的
副本，方法见本节末）：`limits.calculate_from_products / calculate_limits` 及其 legacy 辅助函数
（只被旧 inference 使用；`research/observed.py` 只用到 `LimitOptions` 与 `_validate_request`）、
`observations.load_1d/2d_observations / select_independent_1d_indices`、
`model.load_fit_parameters`（历史 pickle 读入）、`statistics` 中三个未用函数、
`polarization.circular_window / integrated_evpa / radial_integration_matrix`、
`polarimetry` 中四个未用函数、`angles.circular_mean / unwrap_around`、
`eht_official` 中三个别名、`eht_optimizer.solve_fixed_p`、
`products.build_observations / load_observations`、`layout` 的 legacy 读入与 EHT pickle 路径、
`common` 中的 `stage.json`/“imported product” 兼容逻辑。
另外删除了库函数签名里**隐藏的科学默认值**（`LimitOptions`、`evpa_templates`、`omega`、
`cloud_amplitude_scaling`、`calculate_observed` 的 alpha 网格、常量 `ALPHAS / AXION_NORM /
TIME_CONVERSION`）——它们是 YAML 数值的副本；现在所有值都由调用者从 `experiment.yaml` 传入。

### 1.2 行数（`.py` 与 `.yaml`）

| 目录 | 整理前（main 5b1807b） | 整理后 | 变化 |
|---|---|---|---|
| `eht_axion/` | 42 个文件，10 055 行 | 29 个文件，5 266 行 | −48 % |
| `scripts/` | 60 个文件，9 578 行 | 17 个文件，3 048 行 | −68 % |
| `tests/` | 75 个文件，13 036 行 | 36 个文件，5 977 行 | −54 % |
| `configs/` → `experiment.yaml` | 11 个文件，886 行 | 1 个文件，133 行 | −85 % |

测试：整理前 1715 个（eht-rebuild 环境），整理后 658 个，全部通过；默认环境 647 通过、
11 个因缺少 ehtim/astropy 跳过。`ruff`（E,F,W，行宽 110）无告警（main 上有 49 条）。

---

## 2. `experiment.yaml` 的结构

一个过程只能读 `shared` 和它之前各过程的段。

| 段 | 一句话说明 |
|---|---|
| `shared` | 不止一个过程用到的值：M87* 质量与距离、快照标签步长、alpha 网格、主量子数 n，以及公共天空网格（64 像素、120 µas）和 19 µas 匹配模糊（EHT 的 display19 产品与模拟图像用的是**同一个**物理量：天空平面上对 I/Q/U/V 施加的圆高斯核 FWHM，已核对 ehtim `blur_circ` 与 `pipeline.raptor.preprocess_frame`） |
| `raptor` | 过程 1：相机宽度 46.2261 rg、build 选项、等离子体参数、快照区间、响应/背景分辨率与轴子归一化、径向函数表及 spin、alpha→表的对应、BHAC 网格 |
| `eht_imaging` | 过程 2：原 `eht-official.yaml` 去掉移入 `shared` 的 npix/fov/display_blur 以及四个机器路径；保留 uvfits 文件名 + sha256（数据身份）、权重、轮数、重启方案、D 项统计、抽样数与种子 |
| `analysis` | 过程 3：旋转角、孔径、特征选择规则（强度/偏振阈值 0.1/0.2、像素步长 5、10° 扇区间隔与 4° 偏移、最小半径 2 rg、mask 模式、resultant ≥ 0.7、centering native）、云振幅 b(α) 参数、背景时滞与方法、推断网格/先验/时间 |

设计决定（均为“每个值只写一次”）：

* `cloud.reference_norm` 与 `producer.response.reference_norm` 是同一个量（模板振幅以 RAPTOR
  轴子归一化为单位）→ 只写在 `raptor.response.reference_norm`。
* `background.source_resolution` 与 `producer.background.resolution` 相同 → 只写在 `raptor.background`。
* `background.alpha: 0.3` 只是无轴子背景帧的目录名 → 变成 `pipeline/layout.py` 的常量
  `BACKGROUND_DIRECTORY_ALPHA`（不是科学参数）。
* RAPTOR 的 16 个相位（k·π/8）由 RAPTOR 内核决定 → `producer/plan.py` 的 `PHASE_COUNT`，不在 YAML 中。
* 网格/先验检查（phase 128 & 240 点；phase 128 & 260 点 & coupling_max 1e4）不再是 YAML 副本，
  而是 `analysis_check_grid.py` 的命令行覆盖 `--phase-count/--coupling-points/--coupling-max`，
  覆盖值记录在该检查的 manifest 里。`analysis_limits.py --centering centroid` 只用于质心对比。
* 机器路径不进 YAML：RAPTOR 路径在 `.env.raptor`，EHT 成像的 uvfits 目录、输入收据、馈源旋转表、
  官方脚本在 `.env`（`EHT_IMAGING_*`）或命令行；它们仍按原键名记录进成像计划的 `parameters`。
  uvfits 文件名与 sha256 属于数据身份，留在 `eht_imaging.inputs`。
* 原先 `eht_sample`（批量 D 项蒙特卡洛）与 `eht_summarize`（系综汇总）中硬编码的抽样数、种子、
  fsys、D 项统计、站点、网格、显示模糊、孔径半径/宽度、resultant 阈值，现在都从 `experiment.yaml` 读取。
* 名称 `shared` 保留（它同时包含天空网格与模糊核，不只是物理常数）。

---

## 3. 按顺序运行三个过程

```bash
pixi install --locked && pixi run check
# 过程 1：RAPTOR（环境 raptor-build；路径来自 .env.raptor）
pixi run -e raptor-build raptor-build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json
pixi run -e raptor-build raptor-generate
pixi run -e raptor-build raptor-merge --plan plan.json --job-root RUN_A --job-root RUN_B --output-root MERGED
# 过程 2：EHT 成像（环境 eht-rebuild；uvfits 位置来自 .env）
pixi run -e eht-rebuild eht-reconstruct --stage I --day 95 --output-root I95 --review-file REVIEW
pixi run -e eht-rebuild eht-reconstruct --stage fiducial --day 95 --i-root I95 --output-root PD95 --review-file REVIEW
pixi run -e eht-rebuild eht-sample --bundle-root BUNDLE --output-root MC
pixi run -e eht-research eht-summarize --dataset-root MC_COPY --output-root ENSEMBLE
# 过程 3：分析（默认环境）
pixi run analysis-templates --raptor-input-root MERGED/data --output-root TEMPLATES
pixi run analysis-limits --templates TEMPLATES --ensemble ENSEMBLE --output-root LIMITS
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root GRID --phase-count 128 --coupling-points 240
pixi run analysis-check-grid --source LIMITS --templates TEMPLATES --output-root PRIOR --phase-count 128 --coupling-points 260 --coupling-max 10000
pixi run analysis-check-covariance --source LIMITS --templates TEMPLATES --output-root COV
pixi run analysis-summarize --primary LIMITS_CENTROID --native LIMITS --gridcheck GRID --priorcheck PRIOR --covcheck COV --output-root SUMMARY
```

脚本一览（过程顺序）：`raptor_build / raptor_generate / raptor_merge`，
`eht_reconstruct / eht_refine / eht_refine_d / eht_sample / eht_summarize`，
`analysis_templates / analysis_limits / analysis_check_grid / analysis_check_covariance / analysis_summarize`，
以及 `verify_products`。`analysis_templates.py` 取代了原来的 `start_eht_120_templates.py` + `_stage.py` +
`workflow.run_stage` 的 RAPTOR 分支，产生同名阶段与相同参数（`raptor:preprocess`、`raptor:fit`、
`raptor:summarize`、`raptor:background-preprocess`、`raptor:background-unmasked`）。

溯源（provenance）保持不变：每个阶段目录有 `manifest.json`，记录参数、输入/输出 SHA-256、代码哈希、
运行环境；`pixi run verify DIR --recursive` 重新检查整条血统。只删除了纯“副本一致性”检查。

---

## 4. 不兼容旧产品；作业 id 与一次性迁移

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

## 5. 等价性证据

所有证据 JSON 同时保存在 `reports/cleanup-20261106/` 与
`/home/canoziia/GRMHD_DATA/analysis/cleanup-verify-opus-20261106/evidence/`。

RESULTS_PLACEHOLDER

---

## 6. 未验证或遗留的问题

OPEN_PLACEHOLDER
