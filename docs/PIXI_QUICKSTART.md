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
pixi run check  # 安装后可先检查一次
pixi shell
python scripts/run_all.py --help
```

Pixi只下载Python和依赖，**不会安装本项目为包**。任务直接执行`scripts/`下的脚本，
脚本调用本地可复用模块；不注册`eht-pipeline`总命令，也不要求pip install或editable安装。
推荐用`pixi shell`一次进入项目环境，后面直接执行`python scripts/…`。
不需要再套一层其他虚拟环境，也不要使用仓库旧的ARM `.pixi`环境。
退出项目shell用`exit`；重新打开终端时再次进入`pixi shell`。
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

### 3.2 在 Pixi shell 中加载 .env

先完成第1节的`pixi shell`，然后在项目根目录加载`.env`。
如果使用tmux，建议先进入tmux，再执行`pixi shell`和下面的加载命令。
Bash、Zsh使用：

```bash
set -a
source .env
set +a
```

`set -a`让变量自动导出给Python子进程；仅`source .env`而不export，Python看不到普通shell变量。
`pixi shell`只激活Python环境，不会自动加载`.env`，这两步都需要执行。
每个新终端需要加载一次，修改`.env`后也要重新加载。关闭自动导出不会取消已有变量。
**只source可信的本机文件**：shell会执行文件内容，Python脚本不会自动查找或执行`.env`。

### 3.3 顺序执行，不再重复填写路径

```bash
python scripts/process_raptor.py       # ① RAPTOR图像、相位响应、时间分布
python scripts/process_eht.py          # ② EHT topset、定位、观测统计
python scripts/estimate_background.py  # ③ 无轴子背景散布
python scripts/build_observations.py   # ④ 选点、日差和协方差
python scripts/calculate_limits.py     # ⑤ 六条限制曲线
python scripts/plot_limits.py          # ⑥ 图片与表格
```

脚本按以下优先级选择设置：**显式命令行参数 > 当前shell中的`EHT_*`环境变量**。
环境变量没有值且必要参数未指定时会报出缺少的名称，而不是猜测输入或输出。
参数覆盖仅影响这次命令，不会修改`.env`或shell：

```bash
python scripts/calculate_limits.py --config configs/observed.toml --output-root ../eht-pipeline-runs/observed-run
python scripts/plot_limits.py --rebuild
python scripts/process_eht.py --no-trust-pickle
python scripts/calculate_limits.py --no-rebuild --no-reuse-imported
```

布尔环境变量支持true/false、1/0、yes/no、on/off，不区分大小写；非法值明确报错。
`--no-*`可以关闭环境中已设为true的选项。对同一次运行的后续步骤使用同一配置/输出，
例如observed计算后画图也应指定相同的observed配置和输出，或统一修改`.env`并重新加载。

## 4. 一次执行整条主链

```bash
# 已进入pixi shell，并按第3节加载和确认.env后：
python scripts/run_all.py
```

在`pixi shell`里，直接执行脚本即可，不必反复加`pixi run`。
`raptor`、`inference`等只是Pixi任务名，不是shell命令；应使用上面的脚本名。
如果不进入交互shell，仍可用`pixi run pipeline`或`pixi run raptor`等已有任务。
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
python scripts/verify_products.py ../eht-pipeline-runs/my-pixi-run/eht
python scripts/plot_diagnostics.py --run-root ../eht-pipeline-runs/my-pixi-run --output-dir ../eht-pipeline-runs/my-pixi-run/diagnostics
python scripts/check_convergence.py --run-root ../eht-pipeline-runs/my-pixi-run --output ../eht-pipeline-runs/my-pixi-run/convergence.json --modes std cor cor_grmhd
```

图4诊断旁路和HDF5转换：

```bash
python scripts/transport_diagnostics.py --help
python scripts/raptor_tools.py convert --help
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
- Python环境由Pixi管理；Git LFS仍需系统单独安装。
- 可查看全部任务：`pixi task list`。
- 停止/中断时以各阶段`manifest.json`为完成标记；run根目录manifest只是运行描述，不是阶段完成标记。
- 旧`stage.json`可校验；新文件不再用它命名。
- 迁移旧结果：`pixi run migrate run --source /path/to/old-run --destination /path/to/new-run`。
  大NPY复用源文件，响应三列精确转为named NPZ；没有重新计算科学结果。
- 使用迁移产物时明确加`--reuse-imported`。它校验配置/来源/字节后复用，不冒称本版本重算。
  要验证新代码实际计算行为，请使用新的输出目录完整运行。
- 具体科学定义与待确认假设见 [SCRIPT_PIPELINE.md](SCRIPT_PIPELINE.md)。

## 附录：合作者在另一台机器首次取得数据（本机不需要）

数据由我们提供，但`.pixi/`环境和`data/inputs/`本地链接视图不随Git提交。
合作者拿到新的代码checkout后，仍需取得我们发布的数据文件，不是自行生成或整理科学数据。

如果共享数据已经以真实文件提供，直接把`--input-root`指向收到的目录即可，旧发布结构也支持。
如果只有Git LFS指针，则先在包含已发布共享数据的checkout中执行：

```bash
# 在提供数据的原项目checkout中执行，不是在这个已清理的代码worktree中：
git lfs install
git lfs pull --include='shared_data/analysis-inputs/draft/**'
```

下载完成后，可直接使用那个`shared_data/analysis-inputs/draft`作为输入，或在代码目录中
一次性建立规范访问入口：

```bash
pixi run migrate inputs --source /path/to/downloaded/shared_data/analysis-inputs/draft --destination data/inputs
```

迁移只是建链接并核对数据，不是重算，也不复制21GB数组。
也可以不迁移，直接读取收到的旧共享目录。
具体命名结构见[数据命名说明](DATA_LAYOUT_PLAN.md)。
