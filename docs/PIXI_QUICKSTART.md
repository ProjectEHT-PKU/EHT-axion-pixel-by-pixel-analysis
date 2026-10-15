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
另一份输入数据。直接使用已有的结果目录，不额外建运行名称子目录：

```text
../eht-pipeline-runs
```

从当前代码目录解析后，它是：

```text
/home/canoziia/MyProject/PHYS/eht-pipeline-runs
```

运行后，`raptor/`、`eht/`、`inference/`、`figures/`等直接位于此目录下。
同配置重复运行会复用校验通过的缓存；配置改变时会明确报过期缓存，不会静默覆盖。
只有要同时保存另一套独立分析时，才按需指定其他输出位置。

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
EHT_CONFIG=configs/research.yaml
EHT_INPUT_ROOT=data/inputs
EHT_OUTPUT_ROOT=../eht-pipeline-runs
EHT_TRUST_PICKLE=false
EHT_REBUILD=false
EHT_REUSE_IMPORTED=false
```

本机数据已经提供，输入路径通常不需修改。输出直接保存到`eht-pipeline-runs`，无需另取运行名。
确认使用的是我们提供的可信EHT数据后，**手动将`EHT_TRUST_PICKLE`改为`true`**。
模板默认false，不因加载环境变量就默认授权执行pickle。
`EHT_REBUILD`和`EHT_REUSE_IMPORTED`通常保持false；需要时用命令行临时开启。

`.env`是本机配置，不进入Git；`.env.example`只包含分析/共享设置。
可选生成端另用`.env.raptor.example`和本机`.env.raptor`，见[RAPTOR生成端指南](raptor-generation.md)。
物理方法和数值参数仍放在
`configs/research.yaml`，`.env`不替代科学配置文件。科学配置使用YAML；只有Pixi自身的
依赖和任务配置仍是`pixi.toml`，那是Pixi规定的格式，不需要日常修改。
代码仅支持观测EVPA日差固定为0的假设，配置无需`observation`选项；
观测标准差、协方差仍然参与计算，不是把观测误差设为0。
不要在同一mapping中重复写同一个参数名；读取器会拒绝重复键，避免后值静默覆盖前值。

### 3.2 脚本自动读取本机配置

六个主分析阶段和`run_all.py`每次启动都会自动读取**仓库ROOT下的`.env`**。
首次复制并审阅模板后即可直接`pixi run`；新终端或修改文件后都不需要手动导出或加载。
脚本不向cwd或其父目录搜索`.env`。生成端则依次读取`ROOT/.env`、`ROOT/.env.raptor`，
后者叠加前者；分析入口不读取`.env.raptor`。

解析使用`python-dotenv`且`interpolate=False`，只读取变量值，**不执行shell、不展开变量**。
路径请填写实际值，不要依赖`$HOME`、`${VAR}`或命令替换。文件提供的相对路径统一以仓库ROOT为基准，
即使自定义env文件放在别处也一样；CLI和外部已导出的环境变量中的相对路径仍以当前cwd为基准。

可用`--env-file PATH`替代默认文件选择（生成端也不会再叠加默认两文件），
`PATH`本身相对当前cwd。`--no-env-file`禁用文件加载，不会清除外部环境变量；
要独立于本机env文件运行，请显式给出所需CLI输入，例如：

```bash
pixi run pipeline --env-file /path/to/local.env
pixi run pipeline --no-env-file --config configs/research.yaml \
  --input-root data/inputs --output-root ../eht-pipeline-runs --no-trust-pickle
```

上例保持pickle授权关闭；若输入含pickle，仍须确认可信后显式使用`--trust-pickle`。

### 3.3 顺序执行，不再重复填写路径

```bash
pixi run raptor        # ① RAPTOR图像、相位响应、时间分布
pixi run eht           # ② EHT topset、定位、观测统计
pixi run background    # ③ 无轴子背景散布
pixi run observations  # ④ 选点、日差和协方差
pixi run inference     # ⑤ 六条限制曲线
pixi run figures       # ⑥ 图片与表格
```

脚本按以下优先级选择设置：**显式命令行参数 > shell已导出的`EHT_*`环境变量 > env文件值**。
旧终端中已导出的值会遮住文件修改；遇到配置不符时先检查这些外部变量。
各来源都未提供必要参数时会报出缺少的名称，而不是猜测输入或输出。
参数覆盖仅影响这次命令，不会修改`.env`或shell：

```bash
pixi run inference --output-root /另一个结果位置
pixi run figures --rebuild
pixi run eht --no-trust-pickle
pixi run inference --no-rebuild --no-reuse-imported
```

布尔环境变量支持true/false、1/0、yes/no、on/off，不区分大小写；非法值明确报错。
`--no-*`可以关闭环境中已设为true的选项。对同一次运行的后续步骤使用同一配置/输出，
例如临时换输出位置计算后，画图也应指定相同输出，或统一修改`.env`，下一次调用自动读取。

## 4. 一次执行整条主链

```bash
# 已按第3节填写并确认.env后，脚本自动读取：
pixi run pipeline
```

`pixi run pipeline`只是`python scripts/run_all.py`的任务别名，不是安装的统一CLI。
其他步骤同样由Pixi直接运行对应脚本，无需项目安装或交互环境激活。
脚本参数可临时覆盖自动读取的env文件或外部环境中的配置/路径/安全开关；计算与缓存机制仍在可复用模块中。

默认顺序就是上述①至⑥。重复执行时，配置/输入/代码一致且输出校验通过的阶段会显示
`CACHE HIT`。只要修改了配置或代码，可能收到 `Stale cache`，不是程序故障。
默认直接使用上述结果目录；遇到过期缓存时，确认影响范围后再对所需阶段使用`--rebuild`。

不同 Python/NumPy/SciPy 环境可能有微小数值差异。不要把此前 Nix 环境的输出宣称为
Pixi 重算结果；Pixi 验证记录将单独保留。生产stage记录输入/配置/源码，调用环境写入`invocations/`；精确依赖由
`pixi.lock`固定。切换环境时建议新输出目录。

## 5. 看中间量与结果

```text
config.yaml                             本次配置快照
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
pixi run verify ../eht-pipeline-runs/eht
pixi run diagnostics --run-root ../eht-pipeline-runs --output-dir ../eht-pipeline-runs/diagnostics
pixi run convergence --run-root ../eht-pipeline-runs --output ../eht-pipeline-runs/convergence.json --modes std cor cor_grmhd
```

图4诊断旁路和HDF5转换：

```bash
pixi run transport --help
pixi run convert-hdf5 --help
```

`configs/research.yaml`是唯一科学配置。推断只使用零观测日差，不能切换到实际日差。
已有NumPy cubes时无需HDF5转换。

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

其他alpha和日期同样已由我们整理好；合作者只需填写实际路径，再按第3节直接`pixi run`，脚本自动读取`.env`。
数据放在代码目录的`data/inputs/`时，可以保留模板的`EHT_INPUT_ROOT=data/inputs`。

**交付注意：**本机的规范输入目前是链接到原始归档的视图，不是独立可搬运的数据包。
向另一台机器交付时，我们需要提供真实数据文件，或在目标机器建立有效链接，不能只复制
指向本机绝对路径的软链接。新的独立规范数据包尚未在本次文档修改中发布；实际下载地址
和交付方式应由我们随数据另行提供，不能把代码checkout当作已经包含全部数据。

如果数据通过Git LFS交付，下载指令随相应数据仓库提供，不要求合作者理解旧归档命名。
历史路径兼容与迁移工具仅供维护者使用，不是合作者从头运行的必需步骤。
具体命名结构见[数据命名说明](DATA_LAYOUT_PLAN.md)。
