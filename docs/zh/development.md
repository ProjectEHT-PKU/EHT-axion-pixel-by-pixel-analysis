# 开发

## 用途

本页说明仓库的软件环境：Pixi 环境和任务，如何运行测试和代码检查工具（linter），如何通过 `.env`
文件给出本机设置，以及新代码必须遵守的规则。本仓库不是可安装的程序包；脚本直接运行或通过
`pixi run` 运行。

## 输入

- `pixi.toml` 和 `pixi.lock`：环境和任务。
- `.env` 和 `.env.raptor`：本机路径和执行设置，由 `.env.example` 和 `.env.raptor.example` 复制而来
  （两个文件都被 git 忽略）。

## 输出

无；测试只写入临时目录。

## 运行方法

```bash
pixi install --locked
pixi run check                       # tests and lint in the default environment
pixi run -e eht-imaging test         # all tests, including those that need ehtim and astropy
pixi run lint
pixi run verify PATH/TO/STAGE --recursive
```

## 细节

### 环境

| 环境 | 增加的内容 | 用于 |
|---|---|---|
| `default` | Python 3.12、numpy、scipy、matplotlib、h5py、PyYAML、pytest、ruff、python-dotenv | 分析脚本、测试、代码检查、`verify` |
| `raptor-build` | C 编译器、make、GSL、HDF5、pkg-config、rsync（仅 linux-64） | `raptor-build`、`raptor-generate`、`raptor-merge` |
| `eht-summary` | astropy | `eht-summarize` |
| `eht-imaging` | astropy 和 ehtim 1.3.2 | `eht-reconstruct`、`eht-refine`、`eht-export-bundle`、`eht-sample` |

所有环境都设置 `MPLBACKEND=Agg`、`OPENBLAS_NUM_THREADS=2`、`OMP_NUM_THREADS=2` 和
`PYTHONUNBUFFERED=1`。

### 任务

| 任务 | 命令 |
|---|---|
| `raptor-build`、`raptor-generate`、`raptor-merge` | `scripts/raptor/*.py` |
| `eht-reconstruct`、`eht-refine`、`eht-export-bundle`、`eht-sample`、`eht-summarize` | `scripts/eht/*.py` |
| `analysis-templates`、`analysis-limits`、`analysis-check-grid`、`analysis-check-covariance`、`analysis-summarize`、`analysis-paper` | `scripts/analysis/*.py` |
| `verify` | `scripts/verify_products.py` |
| `test` | `pytest -q -ra`，带 `testpaths=tests` 和 `pythonpath=.` |
| `lint` | `ruff check --isolated --target-version py312 --line-length 110 --select E,F,W eht_axion tests scripts` |
| `check` | `test` 和 `lint` |

### 测试

测试只使用合成输入；从不读取观测数据，也不运行 RAPTOR。在默认环境中，需要 ehtim 或 astropy 的测试
被跳过；`eht-imaging` 环境运行全部测试。有些测试不是某个函数的单元测试，而是对仓库规则的检查：

- `test_experiment.py` 检查 `configs/experiment.yaml` 每个节的校验，以及 `eht_axion` 中没有任何公共函数带
  数值默认值（一个简短的允许列表涵盖执行设置和数组布局参数）。
- `test_docs_links.py` 检查 `docs/en/` 和 `docs/zh/` 两棵文档树：每一页都能从 `README.md` 经相对
  链接到达，每个相对链接和锚点都能解析，每一页恰有一个 H1 标题，每一页都链接回本语言的 `index.md`，
  两棵树的文件名集合相同、对应页面的标题层级结构相同。

### .env 文件

`eht_axion/environment.py` 以数据方式读取 dotenv 文件：不做 shell 展开，不做变量插值，不在调用者的
工作目录中查找。相对路径相对于仓库根目录。未知的键和重复的键都是错误。优先级依次为：命令行选项、
已导出的环境变量、`.env.raptor`（仅对 `scripts/raptor/build.py` 和 `scripts/raptor/generate.py`）、`.env`。
`--env-file PATH` 替换默认文件，`--no-env-file` 禁用它们。

| 键 | 读取者 | 含义 |
|---|---|---|
| `M87_AXION_DATA_ROOT` | 所有脚本 | 默认布局的数据根目录（[数据布局](data-layout.md#默认布局)）；可被 `--data-root` 覆盖 |
| `EHT_IMAGING_INPUT_ROOT` | `scripts/eht/reconstruct.py` | HOPS uvfits 文件所在目录 |
| `EHT_IMAGING_INPUT_RECEIPT` | `scripts/eht/reconstruct.py` | 输入回执 JSON |
| `EHT_IMAGING_FEED_ROTATION` | `scripts/eht/reconstruct.py` | 馈源旋转表 JSON |
| `EHT_IMAGING_OFFICIAL_SCRIPT` | `scripts/eht/reconstruct.py` | 官方 EHT 成像脚本 |
| `EHT_RAPTOR_SOURCE_ROOT`、`EHT_GRMHD_ROOT`、`EHT_RADIAL_LOW`、`EHT_RADIAL_HIGH`、`EHT_RAPTOR_RUN_ROOT`、`EHT_RAPTOR_THREADS`、`EHT_RAPTOR_WORKERS`、`EHT_RAPTOR_CPU_AFFINITY`、`EHT_RAPTOR_CPU_LIST`、`EHT_RAPTOR_TIMEOUT_SECONDS` | RAPTOR 脚本 | 见 [RAPTOR 模拟](raptor-simulation.md) |

分析脚本、`scripts/eht/refine.py`、`scripts/eht/export_bundle.py`、`scripts/eht/sample.py` 和 `scripts/eht/summarize.py` 只从 `.env`
读取 `M87_AXION_DATA_ROOT`；它们的其他路径是命令行选项或数据根目录之下的默认值。

### 代码修改规则

- 科学取值应写在 `configs/experiment.yaml` 中，只写一次，并在 `eht_axion/experiment.py` 中校验；库函数以显式
  参数接收它。
- 机器路径和执行设置应放在命令行或 `.env` 中。
- 写出科学产品的脚本使用 `stage_directory`，并把它从仓库之外读取的每个数据文件声明为输入；仓库中的
  文件由提交标识，从不作为输入（[来源记录](provenance.md)）。运行前先提交：受版本控制的文件有未提交
  改动时，结果会被拒绝写出。
- 产品只按当前代码写出的格式读取；改变了的产品需要新的输出根目录。除了设计上可续跑的地方，脚本都
  拒绝已有的输出根目录（[来源记录](provenance.md#输出根目录)）。

## 相关页面

- [文档索引](index.md)
- [包结构](package-map.md)
- [来源记录](provenance.md)
- [实验配置](experiment-config.md)
- [RAPTOR 模拟](raptor-simulation.md)
