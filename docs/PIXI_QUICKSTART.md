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

Pixi只下载Python和依赖，**不会安装本项目为包**。任务直接执行`scripts/`下的脚本，
脚本调用本地可复用模块；不注册`eht-pipeline`总命令，也不要求pip install或editable安装。
通过`pixi run`执行任务，不需要额外激活交互环境，也不要使用仓库旧的ARM `.pixi`环境。
首次安装需要网络；之后不会每次重装。`--locked` 保证不擅自修改锁文件。

本文面向Linux x86-64/ARM64与macOS Intel/Apple Silicon。
真实安装和执行验证以报告为准；其他平台锁定成功不代表已在这些机器上跑过。
当前锁文件的系统基线为Linux glibc≥2.28、macOS≥13；很旧的Linux
服务器（如glibc2.17）不能假定可直接安装，应使用更新的计算节点或容器。

## 2. 直接使用已提供的数据，指定结果保存位置

**本机输入已经准备好，不需要你下载、搬运、改名或重新整理数据。**
在上一步的代码目录中，直接使用：

```text
data/inputs/
```

它已经包含全部八组RAPTOR图像、res400分片、四天EHT候选图和传输诊断数据的访问入口。
其中大文件链接到原来的共享数据，因此不要删除原数据。这里不是要求用户建立一套
“推荐结构”，而是在说明我们已经提供好的输入位置。

`--output-root`只是告诉程序**把本次新生成的中间量和结果存到哪里**，不是要求你准备
另一份输入数据。例如下文统一使用：

```text
../eht-pipeline-runs/my-pixi-run
```

从当前代码目录解析后，它是：

```text
/home/canoziia/MyProject/PHYS/eht-pipeline-runs/my-pixi-run
```

这个目录不需要手动创建，程序会自动创建。第一次从头运行选一个尚未使用的名字，
是为了不覆盖已有的`research-001`等结果；后续同配置运行可以继续用同一目录复用缓存。

**你只需安装环境、直接执行下面的命令。**另一台机器首次取得数据的步骤见文末附录，
本机使用时跳过。

## 3. 按科学顺序运行：推荐使用 .env

### 3.1 一次性填写本机设置

在代码根目录中复制我们提供的模板（若已有`.env`，不要覆盖它）：

```bash
cp .env.example .env
```

打开`.env`，确认：

```dotenv
EHT_CONFIG=configs/research.toml
EHT_INPUT_ROOT=data/inputs
EHT_OUTPUT_ROOT=../eht-pipeline-runs/my-pixi-run
EHT_TRUST_PICKLE=false
EHT_REBUILD=false
EHT_REUSE_IMPORTED=false
```

本机数据已经提供，输入路径通常不需修改。输出目录选择本次运行的名字。
确认使用的是我们提供的可信EHT数据后，**手动将`EHT_TRUST_PICKLE`改为`true`**。
模板默认false，不因加载环境变量就默认授权执行pickle。
`EHT_REBUILD`和`EHT_REUSE_IMPORTED`通常保持false；需要时用命令行临时开启。

`.env`是本机配置，不进入Git；官方只提交`.env.example`。物理方法和数值参数仍放在
`configs/research.toml`，`.env`不替代科学配置文件。

### 3.2 加载到当前 shell

在项目根目录中加载`.env`。Bash、Zsh使用以下命令（也可在tmux中执行）：

```bash
set -a
source .env
set +a
```

`set -a`让变量自动导出给Python/Pixi子进程；仅`source .env`而不export，Python看不到普通shell变量。
每个新终端需要加载一次，修改`.env`后也要重新加载。关闭自动导出不会取消已有变量。
**只source可信的本机文件**：shell会执行文件内容，Python脚本不会自动查找或执行`.env`。

### 3.3 顺序执行，不再重复填写路径

```bash
pixi run raptor        # ① RAPTOR图像、相位响应、时间分布
pixi run eht           # ② EHT topset、定位、观测统计
pixi run background    # ③ 无轴子背景散布
pixi run observations  # ④ 选点、日差和协方差
pixi run inference     # ⑤ 六条限制曲线
pixi run figures       # ⑥ 图片与表格
```

脚本按以下优先级选择设置：**显式命令行参数 > 当前shell中的`EHT_*`环境变量**。
环境变量没有值且必要参数未指定时会报出缺少的名称，而不是猜测输入或输出。
参数覆盖仅影响这次命令，不会修改`.env`或shell：

```bash
pixi run inference --config configs/observed.toml --output-root ../eht-pipeline-runs/observed-run
pixi run figures --rebuild
pixi run eht --no-trust-pickle
pixi run inference --no-rebuild --no-reuse-imported
```

布尔环境变量支持true/false、1/0、yes/no、on/off，不区分大小写；非法值明确报错。
`--no-*`可以关闭环境中已设为true的选项。对同一次运行的后续步骤使用同一配置/输出，
例如observed计算后画图也应指定相同的observed配置和输出，或统一修改`.env`并重新加载。

## 4. 一次执行整条主链

```bash
# 已按第3节加载并确认.env后：
pixi run pipeline
```

`pixi run pipeline`只是`python scripts/run_all.py`的任务别名，不是安装的统一CLI。
其他步骤同样由Pixi直接运行对应脚本，无需项目安装或交互环境激活。
脚本参数可临时覆盖.env导出的配置/路径/安全开关；计算与缓存机制仍在可复用模块中。

默认顺序就是上述①至⑥。重复执行时，配置/输入/代码一致且输出校验通过的阶段会显示
`CACHE HIT`。只要修改了配置或代码，可能收到 `Stale cache`，不是程序故障。
初次使用建议换新的输出目录；确认影响范围后再对所需阶段使用 `--rebuild`。

不同 Python/NumPy/SciPy 环境可能有微小数值差异。不要把此前 Nix 环境的输出宣称为
Pixi 重算结果；Pixi 验证记录将单独保留。生产stage记录输入/配置/源码，调用环境写入`invocations/`；精确依赖由
`pixi.lock`固定。切换环境时建议新输出目录。

## 5. 看中间量与结果

```text
config.toml                             本次配置快照
manifest.json                           稳定的run身份说明
invocations/*.json                      每次请求的配置与环境，追加记录
raptor/processed/alpha-0.30/stokes.npy     处理后IQUV，可mmap
raptor/response/alpha-0.30/response-*.npz 每时间片A/φ/baseline具名数组
raptor/templates/{1d,2d}/templates.npz   时间分布
raptor/background/stokes.npy            无轴子背景图像
background/evpa-difference.npz          背景散布/直方图
background/summary.json                 单位与真实时间差
eht/2017-04-05/{stokes.npy,evpa-1d.npy,summary.npz}
eht/{selection.npz,mask.npy}            EHT选择与区域
observations/{1d,2d}/observations.npz    日差与协方差
inference/limits-{1d,2d}.npz            每维度三模式，cuts轴为(mode,alpha)
inference/posterior-{mode}-{dim}.json   likelihood/CDF
figures/limits.png/pdf/svg,limits.csv   图与表
```

检查完整性与生成诊断图：

```bash
pixi run verify ../eht-pipeline-runs/my-pixi-run/eht
pixi run diagnostics --run-root ../eht-pipeline-runs/my-pixi-run --output-dir ../eht-pipeline-runs/my-pixi-run/diagnostics
pixi run convergence --run-root ../eht-pipeline-runs/my-pixi-run --output ../eht-pipeline-runs/my-pixi-run/convergence.json --modes std cor cor_grmhd
```

图4诊断旁路和HDF5转换：

```bash
pixi run transport --help
pixi run convert-hdf5 --help
```

默认 `configs/research.toml` 是零日差/null配置；实际观测日差用`configs/observed.toml`，
并使用独立输出，不要覆盖null结果。已有NumPy cubes时无需HDF5转换。

## 6. 一个环境即可

```bash
pixi run test
pixi run lint
```

依赖和检查参数统一在`pixi.toml`中，没有额外的pyproject/pytest.ini/ruff.toml。
请通过以上Pixi任务运行检查，裸执行pytest或ruff不会自动取得这些参数。
所有当前测试都在默认环境运行，不需要额外的`legacy`环境或statsmodels。
旧Notebook、旧工具库和对应测试已从此worktree清除，原始`EHT_ZC`中仍完整保留。
数值回归用的兼容代码和合成测试仍在，清理没有改变科学公式。
本地`.pixi/`是实际运行所需的环境，不是科学数据，也不进入Git。

## 7. 维护规则

- 提交 `pixi.toml` 和 `pixi.lock`，不提交 `.pixi/`。
- 改依赖后由维护者执行 `pixi install` 更新锁，合作者优先 `pixi install --locked`。
- Python环境由Pixi管理；仅在数据交付使用Git LFS时才需另行安装Git LFS。
- 可查看全部任务：`pixi task list`。
- 停止/中断时以各阶段`manifest.json`为完成标记；run根目录manifest只是运行描述，不是阶段完成标记。
- 旧`stage.json`可校验；新文件不再用它命名。
- 迁移旧结果：`pixi run migrate run --source /path/to/old-run --destination /path/to/new-run`。
  大NPY复用源文件，响应三列精确转为named NPZ；没有重新计算科学结果。
- 使用迁移产物时明确加`--reuse-imported`。它校验配置/来源/字节后复用，不冒称本版本重算。
  要验证新代码实际计算行为，请使用新的输出目录完整运行。
- 具体科学定义与待确认假设见 [SCRIPT_PIPELINE.md](SCRIPT_PIPELINE.md)。

## 附录：合作者在另一台机器首次取得数据（本机不需要）

数据由我们提供，合作者不需要自行生成或整理科学数据。代码和数据分开交付；
`.pixi/`环境和本机`data/inputs/`的链接视图不随代码Git提交。

合作者取得我们提供的**规范输入数据目录**后，在`.env`中填写该目录的位置，例如：

```dotenv
EHT_INPUT_ROOT=/home/researcher/eht-data/inputs
```

提供的数据目录应直接包含：

```text
inputs/
  raptor/res-100/alpha-0.30/stokes.npy
  raptor/res-400/alpha-0.30/stokes-t1600-t1690.npy
  eht/2017-04-05/image-ensemble.pkl
  transport/snapshot-2110-res-400/...
```

其他alpha和日期同样已由我们整理好；合作者只需填写实际路径，再按第3节加载`.env`。
数据放在代码目录的`data/inputs/`时，可以保留模板的`EHT_INPUT_ROOT=data/inputs`。

**交付注意：**本机的规范输入目前是链接到原始归档的视图，不是独立可搬运的数据包。
向另一台机器交付时，我们需要提供真实数据文件，或在目标机器建立有效链接，不能只复制
指向本机绝对路径的软链接。新的独立规范数据包尚未在本次文档修改中发布；实际下载地址
和交付方式应由我们随数据另行提供，不能把代码checkout当作已经包含全部数据。

如果数据通过Git LFS交付，下载指令随相应数据仓库提供，不要求合作者理解旧归档命名。
历史路径兼容与迁移工具仅供维护者使用，不是合作者从头运行的必需步骤。
具体命名结构见[数据命名说明](DATA_LAYOUT_PLAN.md)。
