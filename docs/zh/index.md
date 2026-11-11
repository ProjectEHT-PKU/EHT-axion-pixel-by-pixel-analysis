# M87* EHT 轴子分析流水线文档

[English](../en/index.md) | 简体中文

本仓库实现了一套完整的科学计算流水线，用于在事件视界望远镜（Event Horizon Telescope, EHT）2017 年 M87\* 偏振观测中搜寻由轴子引起的电矢量位置角（electric-vector position angle, EVPA）振荡。流水线将偏振广义相对论光线追踪（ray tracing）模拟（`RAPTOR`）、含仪器偏振泄漏（D 项，D-term）蒙特卡罗（Monte Carlo）的 EHT Paper VII 偏振成像，以及针对无量纲耦合常数 $c$ 的相位边缘化（phase-marginalized）贝叶斯推断（Bayesian inference）结合起来。代码以直接运行的命令行脚本形式组织，由 [Pixi](https://pixi.prefix.dev/) 管理运行环境，不打包为可安装的 Python 程序包。

## 快速开始

### 环境要求

- [Pixi](https://pixi.prefix.dev/)（根据 `pixi.lock` 锁定并配置 Python `3.12.*` 及全部科学计算与编译依赖）和 `git`。
- **平台支持**：`default`、`eht-summary` 和 `eht-imaging` 环境支持 `linux-64`、`linux-aarch64`、`osx-64`、`osx-arm64` 和 `win-64`；编译与运行 C 语言光线追踪核心程序 `RAPTOR`（`raptor-build` 环境）仅支持 `linux-64`。
- **洁净的 Git 工作树**：写出任何阶段输出时，仓库受版本控制的文件不得含有未提交的修改（`git status --untracked-files=no`），因为仓库当前的 `HEAD` 提交号是唯一的代码溯源（provenance）身份标识。

### 安装环境与配置路径

克隆仓库并安装所需的锁定 Pixi 环境：

```bash
git clone <repository-url> eht-script-pipeline
cd eht-script-pipeline

pixi install --locked                  # default environment: analysis, lint, core tests
pixi install --locked -e eht-imaging   # EHT imaging + full test suite (adds ehtim 1.3.2 and astropy)
pixi install --locked -e raptor-build  # Linux only: RAPTOR C toolchain and job orchestration
```

将 `.env.example` 复制为 `.env`，并设置本机的三个路径根目录：

```bash
cp .env.example .env
# Edit .env and set:
# EHT_RAPTOR_SOURCE_ROOT=/path/to/raptor
# EHT_INPUT_ROOT=/path/to/m87-axion-inputs
# EHT_DATA_ROOT=/path/to/m87-axion
```

在 `.env` 中设置 `EHT_DATA_ROOT`（或在 shell 中导出，或通过 `--data-root` 传入）后，所有脚本都会按[默认数据布局](data-layout.md#默认目录布局)自动解析输入和输出目录，并在 `EHT_INPUT_ROOT` 的只读[输入目录布局](data-layout.md#输入目录布局)下定位全部外部输入。所有科学参数集中存放在 [`configs/experiment.yaml`](../../configs/experiment.yaml) 中；机器本地路径存放在唯一的 `.env` 中（见[配置参考](configuration.md)）。

### 运行测试与代码检查

所有单元测试和集成测试均在临时目录中使用合成数据运行，既不需要观测数据文件，也不需要编译好的 `RAPTOR` 可执行文件：

```bash
pixi run --locked check                                       # default-env tests + Ruff linter
OPENBLAS_NUM_THREADS=2 pixi run --locked -e eht-imaging test  # full test suite (incl. ehtim & astropy)
pixi run --locked lint                                        # Ruff check only
```

### 输入数据来源

`RAPTOR` 模拟与 EHT 成像是两条相互独立的上游流水线，其产出共同作为下游分析的输入：

1. **`RAPTOR` Stokes 图像数据立方**（`$EHT_DATA_ROOT/raptor/images/`）：可以直接放入已有的模拟归档数据，也可以由广义相对论磁流体力学（general-relativistic magnetohydrodynamics, GRMHD）BHAC 快照（`data001600.dat`–`data003000.dat`）、两张无表头 $(\alpha, r, R)$ 径向轮廓 CSV 表，以及处于提交 `20056006820b6536c96de790dd5fdaad5ada4eae` 的洁净 `github.com/projecteht-pku/raptor` 工作副本依次运行 `raptor-build` 和 `raptor-generate` 生成（由 `.env` 中的 `EHT_INPUT_ROOT` 与 `EHT_RAPTOR_SOURCE_ROOT` 定位；见 [RAPTOR 模拟](raptor-simulation.md)）。
2. **汇总后的 EHT 偏振系综**（`$EHT_DATA_ROOT/eht/ensemble/`）：可以直接放入已有的成像归档数据，也可以基于 2017 年 M87\* 四个观测日（第 `95`、`96`、`100`、`101` 天）的官方低频段 HOPS `zbl-dtcal+selfcal` `uvfits` 文件以及官方 EHT 成像脚本，依次运行 `eht-reconstruct`、`eht-refine`、`eht-sample` 和 `eht-summarize` 重建生成（位于 `$EHT_INPUT_ROOT/eht/`；见 [EHT 成像](eht-imaging.md)）。

### 在已有数据产品上运行分析链

当 `$EHT_DATA_ROOT/raptor/images` 与 `$EHT_DATA_ROOT/eht/ensemble` 就绪后，可在 `default` 环境中运行完整的统计推断：

```bash
pixi run analysis-expected-signal
pixi run analysis-limits
```

每条命令都会在 `$EHT_DATA_ROOT/analysis/` 下发布一个新的阶段目录（`expected-signal/` 和 `limits/`，其中包含 `limits/limits.csv`），并拒绝覆盖任何已存在的输出目录。

若还需在 `$EHT_DATA_ROOT/analysis/paper/` 中生成论文插图、LaTeX 数值宏（`numbers.tex`、`numbers.json`）及表格（`table-limits.tex`、`table-chi2.tex`）——该步骤额外读取 `$EHT_DATA_ROOT/raptor/plan.json`、位于 `$EHT_RAPTOR_SOURCE_ROOT` 的 `RAPTOR` 工作副本中的 `src/axion.c`，以及 `$EHT_INPUT_ROOT/radial/` 中的径向轮廓 CSV 表——请运行：

```bash
pixi run analysis-paper
```

## 文档目录

| 文档页面 | 内容概要 |
|---|---|
| [概览与基本概念](overview.md) | 流水线整体架构、Mermaid 数据流图、物理尺度、偏振符号约定及核心设计准则 |
| [配置参考](configuration.md) | `configs/experiment.yaml` 全部键的逐项说明、校验规则与派生设置，以及 `.env` 中的三个路径根目录与执行设置 |
| [RAPTOR 模拟](raptor-simulation.md) | 固定版本 `RAPTOR` 编译、BHAC 快照与径向表校验、任务标识符、作业调度、CPU 亲和性与图像数据立方 |
| [EHT 成像](eht-imaging.md) | EHT Paper VII 偏振成像复现、凸优化 PDHG 精化、D 项蒙特卡罗、系综汇总，以及与 Paper VII 的已知差异 |
| [分析](analysis.md) | 图像预处理、扇形孔径算子、预期信号、特征选择、协方差模型、贝叶斯似然与后验上限及论文输出 |
| [数据布局与溯源机制](data-layout.md) | 默认 `$EHT_DATA_ROOT` 目录树、只读 `$EHT_INPUT_ROOT` 输入目录布局、各阶段文件与数组规范、原子化 `stage_directory` 清单与非清单回执 |
| [开发指南](development.md) | Pixi 环境与任务、测试套件与代码检查、模块与脚本索引、架构不变量，以及新增参数或分析阶段的扩展指南 |
