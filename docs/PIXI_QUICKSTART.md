# 用 Pixi 按顺序运行科研流水线

不需要 Nix/NixOS，也不需要预先安装 Python、Conda、Jupyter 或 LaTeX。
生产流程不执行 Notebook。

## 1. 安装 Pixi、进入代码目录

按 [Pixi 官方安装说明](https://pixi.prefix.dev/latest/installation/) 安装适合操作系统的版本。
打开一个新的终端，确认：

```text
pixi --version
```

进入仓库根目录（包含 `pixi.toml` 的目录），例如本机：

```bash
cd /home/canoziia/MyProject/PHYS/eht-script-pipeline
pixi install --locked
pixi run help
pixi run check
```

Pixi 会下载 Python 和依赖，并将本项目安装为 editable package。
所有命令前加 `pixi run`，不需要激活另一个虚拟环境，也不要使用仓库旧的 ARM `.pixi` 环境。
首次安装需要网络；之后不会每次重装。`--locked` 保证不擅自修改锁文件。

支持生成锁文件的平台：Linux x86-64/ARM64、macOS Intel/Apple Silicon、Windows x86-64。
真实安装和执行验证以报告为准；其他平台锁定成功不代表已在这些机器上跑过。
当前锁文件的系统基线为Linux glibc≥2.28、macOS≥13、Windows≥10；很旧的Linux
服务器（如glibc2.17）不能假定可直接安装，应使用更新的计算节点或容器。

## 2. 准备输入与新的输出目录

输入结构：

```text
analysis-inputs/draft/
  raptor-stokes/res-100/alpha-*.npy
  raptor-stokes/res-400/alpha-0.3/stokes-t*-t*.npy
  eht-image-ensemble/230_{95,96,100,101}.pkl
  raptor-transport-diagnostics/...
```

必须是真实 NumPy/pickle 文件，不是 Git LFS 的几行指针文本。
首次从远程克隆时需额外安装 Git LFS 并下载数据：

```bash
git lfs install
git lfs pull --include='shared_data/analysis-inputs/draft/**'
```

本机已有输入：

```text
/home/canoziia/MyProject/PHYS/EHT_ZC/shared_data/analysis-inputs/draft
```

输出选一个新的目录，例如：

```text
/home/canoziia/MyProject/PHYS/eht-pipeline-runs/my-pixi-run
```

原始输入只读，不需要复制进代码目录，也不能将原始数据目录作为输出目录。

## 3. 按科学顺序运行

Linux/macOS Bash 示例（可以放在 tmux 中执行）：

```bash
export INPUT=/path/to/analysis-inputs/draft
export OUTPUT=/path/to/my-pixi-run
export CONFIG=configs/research.toml

# ① RAPTOR图像处理、相位响应拟合、时间分布
pixi run raptor --config "$CONFIG" --input-root "$INPUT" --output-root "$OUTPUT"

# ② EHT候选图匹配、topset、定位、观测统计
pixi run eht --config "$CONFIG" --input-root "$INPUT" --output-root "$OUTPUT" --trust-pickle

# ③ 无轴子时间演化背景误差
pixi run background --config "$CONFIG" --input-root "$INPUT" --output-root "$OUTPUT"

# ④ 选点、日差、协方差
pixi run observations --config "$CONFIG" --output-root "$OUTPUT"

# ⑤ 1D/2D × 三种误差模式的六条限制
pixi run inference --config "$CONFIG" --output-root "$OUTPUT"

# ⑥ 结果图、CSV表
pixi run figures --config "$CONFIG" --output-root "$OUTPUT"
```

这些变量只是避免重复输入，程序不依赖它们。Windows PowerShell 可以直接写路径，例如：

```powershell
pixi run raptor --config configs/research.toml --input-root "D:/EHT-data/analysis-inputs/draft" --output-root "D:/EHT-runs/my-run"
pixi run eht --config configs/research.toml --input-root "D:/EHT-data/analysis-inputs/draft" --output-root "D:/EHT-runs/my-run" --trust-pickle
```

再按同样方式依次执行 background、observations、inference、figures。
不要使用 Bash 的 `export` 或 shell脚本作为 Windows 必需步骤。

`--trust-pickle` 意味着你信任 EHT pickle 的来源；pickle可能执行代码，不要对未知输入开启。

## 4. 一次执行整条主链

```bash
pixi run pipeline --config configs/research.toml \
  --input-root /path/to/analysis-inputs/draft \
  --output-root /path/to/my-pixi-run --trust-pickle
```

默认顺序就是上述①至⑥。重复执行时，配置/输入/代码一致且输出校验通过的阶段会显示
`CACHE HIT`。只要修改了配置或代码，可能收到 `Stale cache`，不是程序故障。
初次使用建议换新的输出目录；确认影响范围后再对所需阶段使用 `--rebuild`。

不同 Python/NumPy/SciPy 环境可能有微小数值差异。不要把此前 Nix 环境的输出宣称为
Pixi 重算结果；Pixi 验证记录将单独保留。生产stage目前记录输入/配置/源码，环境精确
依赖由 `pixi.lock` 固定；切换环境时建议新输出目录。

## 5. 看中间量与结果

```text
01_raptor/           处理后的IQUV
02_response/         每个时间片的A、φ、baseline
03_templates/       A、φ时间分布
04_eht/              EHT选择、IQUV、EVPA、mask
05_background/       背景散布/直方图
06_observations/     日差与协方差
07_inference/        六曲线数值、likelihood/CDF、运行配置
08_figures/          limits.png/pdf/svg、limits.csv
```

检查完整性与生成诊断图：

```bash
pixi run verify /path/to/my-pixi-run/04_eht
pixi run diagnostics --run-root /path/to/my-pixi-run --output-dir /path/to/my-pixi-run/09_diagnostics
pixi run convergence --run-root /path/to/my-pixi-run --output /path/to/convergence.json --modes std cor cor_grmhd
```

图4诊断旁路和HDF5转换：

```bash
pixi run transport --help
pixi run convert-hdf5 --help
```

默认 `configs/research.toml` 是零日差/null配置；实际观测日差用`configs/observed.toml`，
并使用独立输出，不要覆盖null结果。已有NumPy cubes时无需HDF5转换。

## 6. 历史工具环境（可选）

新主链不依赖statsmodels。只有旧工具测试/历史库需要时：

```bash
pixi install --locked -e legacy
pixi run -e legacy test
```

默认 `pixi run test` 会明确跳过legacy-utils模块，而不是隐藏失败；`legacy`环境运行完整测试。
本机首次环境磁盘占用约1.5GiB；default和legacy两套环境合计约3GiB（不含原始科学数据）。
通常合作者只需default；`.pixi/`不进入Git。

## 7. 维护规则

- 提交 `pixi.toml` 和 `pixi.lock`，不提交 `.pixi/`。
- 改依赖后由维护者执行 `pixi install` 更新锁，合作者优先 `pixi install --locked`。
- Python环境由Pixi管理；Git LFS仍需系统单独安装。
- 可查看全部任务：`pixi task list`。
- 停止/中断时以stage.json为完成标记，勿手动将临时目录视为有效缓存。
- 具体科学定义与待确认假设见 [SCRIPT_PIPELINE.md](SCRIPT_PIPELINE.md)。
