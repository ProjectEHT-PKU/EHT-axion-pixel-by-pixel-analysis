# 来源记录

## 用途

本页说明每个科学产品如何记录自己的来历：带 `manifest.json` 的阶段目录，记录参数、代码的 git 提交，
以及读写的每个数据文件的 SHA-256；这些目录如何以原子方式创建；以及如何校验一个阶段及其
整条来源链（lineage）。本页还说明两个改用自有回执的地方（EHT 打包和 D 项蒙特卡罗批次）。

## 输入

- 任何已完成的阶段目录，用于校验。
- 写出时：由调用脚本提供的阶段名、可 JSON 序列化的参数记录（`config`）和输入文件列表。

## 输出

- 每个阶段目录中的 `manifest.json`。
- 仅在创建过程中存在：目标旁边的 `<stage>.lock` 和临时目录 `<stage>.building-*`，二者在调用结束
  时都会被删除。

## 运行方法

```bash
pixi run verify PATH/TO/STAGE              # check this stage's outputs
pixi run verify PATH/TO/STAGE --recursive  # also re-hash every input and upstream stage
```

`scripts/verify_products.py` 在一个校验会话中完成检查（每个文件只计算一次哈希，无论有多少个阶段
声明了它），检查成功后打印清单，且从不重新计算数据。

## 细节

### 阶段目录

`stage_directory(destination, stage, config, input_paths, *, verify_lineage=True)`
（`eht_axion/provenance.py`）是一个上下文管理器，每个写出科学产品的脚本都使用它：

- 它拒绝经由目录符号链接到达的目标，也拒绝包含其自身某个输入的目标。
- 它以独占方式创建 `<destination>.lock`；第二个写入者直接失败，而不是等待。
- 如果目标已经存在，就校验它（见[输出根目录](#输出根目录)），并把它记录的阶段名、`config` 和
  `inputs` 与新的比较：相同则复用该阶段（跳过主体，上下文返回 `None`），不同则报错。代码提交不参与
  比较，因此可续跑的 RAPTOR 运行在不改变其任务的提交之后仍能继续。已有的阶段从不被替换。
- 只有当仓库中受版本控制的文件没有改动时才写出新阶段（`code_identity`；未跟踪的文件不影响），
  否则拒绝。
- 否则主体写入同一父目录下的一个临时目录。主体结束时，对每个文件计算哈希，最后写入
  `manifest.json`，再把临时目录重命名为目标。失败时删除临时目录。

### 清单

| 字段 | 内容 |
|---|---|
| `stage` | 阶段名，例如 `analysis:fit` |
| `schema` | 清单模式版本（1） |
| `config` | 该阶段的参数记录 |
| `inputs` | 每个输入文件的绝对路径 `path`、`size`、`sha256` |
| `implementation` | `git_commit`：仓库的完整提交号 |
| `runtime` | Python 版本和实现、平台、机器，以及 numpy、scipy、h5py、matplotlib 和 PyYAML 的版本字符串（不计算哈希） |
| `complete` | `true` |
| `outputs` | 每个输出文件的相对路径 `path`、`size`、`sha256`；`.npy` 还记录形状和数据类型，`.npz` 记录其中每个数组的形状和数据类型 |
| `elapsed_seconds` | 主体的墙钟时间 |

代码只由 git 提交标识。由于只能在受版本控制的文件没有改动时写出阶段，提交就固定了每个脚本、每个
`eht_axion` 模块、`configs/experiment.yaml`、`pixi.toml` 和 `pixi.lock`；这些文件都不计算哈希，也不列为输入。
第三方程序包由 `pixi.lock` 固定；依赖特定程序包的阶段记录其版本（例如 EHT 成像计划中的
`ehtim_version`）。git 不保存的文件都计算哈希：从仓库之外读取的数据输入（uvfits、GRMHD 快照、径向
表、官方 EHT 脚本、上游阶段的文件、编译得到的 RAPTOR 可执行文件）以及每个阶段输出。

### 校验

`verify_stage(directory, recursive=False)` 检查清单是否完整，以及每个声明的输出是否以记录的大小和
SHA-256 存在于目录之内。`recursive=True` 时，它还重新计算每个输入的哈希，并对名为 `manifest.json`、
属于某个已完成阶段的输入递归校验该阶段；每个上游阶段只访问一次。

脚本读取另一个阶段的文件时总是使用 `require_declared_outputs(directory, paths)`：该文件必须是那个
阶段声明的输出，而不能是随意多出来的文件。

`verification_session()` 在每次调用中对每个输入文件只计算一次哈希，并在结束时检查没有任何已计算
哈希的文件发生变化（设备、inode、大小、修改时间和状态改变时间）。

默认情况下，`stage_directory` 在开始之前校验每个输入阶段的完整来源链。小型结果步骤
`scripts/analysis/check_grid.py`、`scripts/analysis/check_covariance.py`、`scripts/analysis/summarize.py` 和 `scripts/analysis/paper.py` 传入
`verify_lineage=False`：它们校验所读取阶段的输出，并把这些清单作为已计算哈希的输入记录下来，但不再
重新计算整条原始来源链的哈希。

### 输出根目录

所有脚本适用同一条规则：产生结果的脚本拒绝已经存在的输出根目录。例外是流水线中可续跑的部分：
`scripts/raptor/build.py`、`scripts/raptor/generate.py` 和 `scripts/raptor/merge.py` 复用其运行根目录中已完成的阶段（阶段名、
配置和输入数据都相同的阶段被复用；同一位置上不同的阶段报错），`scripts/eht/sample.py` 在批次身份不变时于同一输出根目录中续跑
批次。`scripts/analysis/templates.py` 在一个新的根目录下写出多个阶段。

`stage_directory` 遇到已存在的目标时，用与新阶段相同的来源链设置校验它：默认递归校验，
`verify_lineage=False` 时只校验输出。

### 阶段名

阶段名是记录在清单中、由使用方检查的固定标识符。其形式为 `<part>:<step>`，部分为 `raptor`、`eht`
或 `analysis`。清单之外的 JSON 记录带有形式为 `<part>-<record>-v<N>` 的模式字符串
（`raptor-generated-v2`、`raptor-merge-v2`、`eht-bundle-v2`、`eht-sample-batch-v2`、
`eht-sample-task-v2`、`eht-sample-finish-v2`、`eht-ensemble-v2`、`analysis-limits-v2`、
`analysis-summary-v2`、`analysis-paper-v1`、`analysis-paper-numbers-v1`）。

| 阶段 | 写出者 |
|---|---|
| `raptor:build` | `scripts/raptor/build.py` |
| `raptor:job` | `scripts/raptor/generate.py`（每个任务一个） |
| `raptor:delivery` | `scripts/raptor/generate.py`、`scripts/raptor/merge.py`（每个产品一个） |
| `raptor:merge` | `scripts/raptor/merge.py` |
| `eht:i`、`eht:fiducial` | `scripts/eht/reconstruct.py` |
| `eht:refine` | `scripts/eht/refine.py` |
| `eht:ensemble` | `scripts/eht/summarize.py` |
| `analysis:preprocess`、`analysis:fit`、`analysis:templates`、`analysis:background` | `scripts/analysis/templates.py` |
| `analysis:limits` | `scripts/analysis/limits.py` |
| `analysis:check-grid` | `scripts/analysis/check_grid.py` |
| `analysis:check-covariance` | `scripts/analysis/check_covariance.py` |
| `analysis:summary` | `scripts/analysis/summarize.py` |
| `analysis:paper` | `scripts/analysis/paper.py` |

### 阶段目录之外的回执

- `scripts/eht/export_bundle.py` 写出 `bundle.json` 而不是清单：列出每个复制文件的 SHA-256、大小和原始
  路径、科学约定以及 git 提交（`code`）。复制之前校验源阶段（完整来源链），复制得到的文件设为只读。
- `scripts/eht/sample.py` 写出 `batch.json`，其 `identity` 是批次计划的 SHA-256（打包快照的哈希、科学取值、
  git 提交以及 ehtim、numpy、scipy、astropy 和 threadpoolctl 的版本，以及线程数）。每次尝试写出一份
  带有其输出哈希的 `taskreceipt.json`；`finish.json` 列出每次抽样获证书尝试的回执。
  `scripts/eht/summarize.py` 在构建系综阶段之前重新检查这些回执。
- `scripts/raptor/generate.py` 和 `scripts/raptor/build.py` 在运行根目录的 `requests/<uuid>.json` 中追加一条请求
  记录，并把实验文件的副本保存为 `config.yaml`。

## 相关页面

- [文档索引](index.md)
- [数据布局](data-layout.md)
- [RAPTOR 模拟](raptor-simulation.md)
- [EHT 成像](eht-imaging.md)
- [分析](analysis.md)
