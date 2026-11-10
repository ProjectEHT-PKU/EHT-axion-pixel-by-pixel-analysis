# RAPTOR 模拟

## 用途

本页说明 RAPTOR 输入如何产生：编译固定版本的 RAPTOR 核心程序，规划光线追踪任务，把任务作为
按内容寻址的作业运行，校验每个作业的输出，把输出转换为 Stokes 数据立方，并为分析部分交付（或
合并）这些数据立方。RAPTOR 模拟不读取任何 EHT 成像产品。

## 输入

- RAPTOR 仓库的一个本地工作副本（`--source-root` 或 `EHT_RAPTOR_SOURCE_ROOT`）：没有未提交的修改，
  恰好处于 `raptor.source.commit`，并且有一个 git 远程仓库指向 `raptor.source.repository`。
- `raptor.snapshots` 所列标签的 BHAC GRMHD 快照 `data<label>.dat`（`--fluid-root` 或
  `EHT_GRMHD_ROOT`）。
- 两张无表头的 CSV 径向表，列为 $\alpha, r, R$（`--radial-low`、`--radial-high` 或
  `EHT_RADIAL_LOW`、`EHT_RADIAL_HIGH`）。
- `experiment.yaml`：`shared` 中的物理键和 `raptor` 节（[实验配置](experiment-config.md)）。

## 输出

所有输出都在一个运行根目录之下（`--run-root` 或 `EHT_RAPTOR_RUN_ROOT`；默认布局使用
`raptor/runs/<name>/`，见[数据布局](data-layout.md#默认布局)）：

- `build/`：编译好的 RAPTOR 及其回执（阶段 `raptor:build`）。
- `jobs/<task id>/`：每个任务一个阶段（阶段 `raptor:job`）；失败的尝试保存在
  `attempts/<task id>/` 中。
- `data/raptor/res-100/alpha-<α>/stokes.npy`：每个 $\alpha$ 一个交付的数据立方；以及
  `data/raptor/res-400/background/t<A>-t<B>/stokes-t<A>-t<B>.npy`：无轴子的时间分片
  （阶段 `raptor:delivery`）。
- `config.yaml`（实验文件的副本）和 `requests/<uuid>.json`（每次调用一条记录）。

`raptor_merge.py` 把同样的 `data/raptor/...` 目录树写入一个新的输出根目录，另外写出
`inventory/`（阶段 `raptor:merge`）。文件细节见[数据布局](data-layout.md)。

## 运行方法

使用 `raptor-build` 环境（仅限 Linux）。路径依次取自命令行、已导出的环境变量、`.env.raptor`、
`.env`（[开发](development.md)）。设置了数据根目录（`M87_AXION_DATA_ROOT` 或 `--data-root`）时，下面每个根目录选项都有取自
[默认布局](data-layout.md#默认布局)的默认值；这里的命令显式给出了这些根目录。

```bash
pixi run -e raptor-build raptor-build                       # build into RUN_ROOT/build
pixi run -e raptor-build raptor-generate --plan --plan-output plan.json   # list tasks only
pixi run -e raptor-build raptor-generate                    # run jobs, then deliver cubes
pixi run -e raptor-build raptor-generate --convert-only     # deliver already complete jobs
pixi run -e raptor-build raptor-merge --plan plan.json \
    --job-root RUN_A --job-root RUN_B --output-root MERGED  # merge runs of several hosts
```

`raptor_generate.py` 的选项：`--max-jobs N`（只运行前 N 个任务；交付要等其余任务完成）、
`--threads`（每个 RAPTOR 进程的 OpenMP 线程数）、`--workers`（并发进程数）、
`--cpu-affinity none|physical|logical`、`--cpu-list`、`--timeout-seconds`。编译和生成都接受
`--config`、`--source-root`、`--run-root`、`--env-file` 和 `--no-env-file`。

## 细节

### 编译

`build_raptor`（`eht_axion/raptor/build.py`）检查源码工作副本没有未提交的修改（包括未跟踪的文件）、
HEAD 恰好是 `raptor.source.commit`，并且至少有一个 git 远程仓库地址在规范化之后（去掉协议、用户、端口
和末尾的 `.git`，`git@host:` 与 `https://host/` 视为相同，转为小写）等于 `raptor.source.repository`；然后
把受版本控制的源码字节复制到 `build/source/`，并用由 `raptor.build` 派生的参数（`--code`、
`--metric`、`--int`、`--grid`、`--rad`、`--sf`、`--bflip`、`--sigma-cut`、`--uniform`、
`--thetae-min`、`--thetae-max`）运行官方的 `setup.sh`。它记录编译器命令行和本次编译的预处理宏
（`macros.txt`），并核对这些宏与 YAML 设置一致。编译记录包含这些设置、RAPTOR 仓库与提交和工具链（源码由提交标识，不再逐文件计算哈希），以及编译得到的可执行
文件的 SHA-256，每个作业都把它记录为输入；径向表和快照不是编译输入。

### 任务计划

`plan_generation`（`eht_axion/raptor/plan.py`）只读取小型输入：

- 对每个标签，它读取 `data<label>.dat` 的 BHAC 尾部信息（小端序块布局、以 $t_g$ 为单位的时间、
  绝热指数、自旋），并检查：自旋等于 `raptor.fluid_spin`（0.9375），块大小能整除
  `raptor.grid.nxlone`，文件长度不短于尾部信息所隐含的长度，以及相邻尾部时间之差等于
  `label_step_tg` × 标签差。
- 每张径向表被载入，检查它是否为覆盖使用它的那些 $\alpha$ 值的矩形递增 $(\alpha, r)$ 网格，并计算
  哈希值。振幅不做重新缩放。

对每个快照，计划包含：

- 一个分辨率为 `response.resolution`（100）的无轴子（“baseline”）任务和一个分辨率为
  `background.resolution`（400）的无轴子任务；
- 对每个 $\alpha$ 和每个相位序号 $k = 0,\dots,15$，一个 100 像素的有轴子（“response”）任务，其中
  `AXION_ALPHA` $=\alpha$，`AXION_OMEGA` $=\Omega(\alpha)=\alpha(1-\alpha^2/(2n^2))$，
  `AXION_NORM` $=$ `reference_norm`，`AXION_PHASE` $=k$（核心程序施加的相位为 $k\pi/8$），并使用
  `raptor.alpha_table` 指定的径向表。

按当前取值共有 $141 \times 130 = 18\,330$ 个任务。每个任务带有 19 个 `model.in` 取值（以 $M_\odot$
为单位的 `MBH`、以 kpc 为单位的 `DISTANCE`、以 g 为单位的 `M_UNIT`、`R_LOW`、`R_HIGH`、以度为单位的
`INCLINATION`、图像尺寸、以 $r_g$ 为单位的相机尺寸、一个频率、`STEPSIZE`、四个轴子取值、
`AXION_TABLE`、`MAX_LEVEL` = 1）及其输出槽位（baseline 为 0，相位 $k$ 为 $k+1$）。

### 任务标识符

任务标识符是一个规范化 JSON 文档（`task_identity`）的 SHA-256，该文档包含：

- `schema_version`（2）；
- 快照文件名及其尾部元数据；
- 除表路径之外的 19 个 `model.in` 取值；
- 对有轴子任务，表的自旋和表内容的 SHA-256；
- `raptor.grid` 和 `raptor.build`。

天空网格、EHT 成像和分析设置从不进入任务标识符，因此改变它们永远不会使光线追踪结果失效。

### 作业

`run_jobs`（`eht_axion/raptor/jobs.py`）用固定数目的工作槽位运行任务。启动之前，它校验编译阶段，检查编译的仓库和提交等于计划中的 `raptor_source`（由 `raptor_generate.py`
从 `raptor.source` 复制而来），
（包括其来源链），重新读取每个快照尾部和径向表并与计划比对（径向表哈希总是比对；计划中含有快照
哈希时也比对快照哈希），把编译目录的网格文件和宏与计划比对，并检查内存余量。每个快照和每张表的
SHA-256 记录在读取它的每个作业的清单中。每个任务在它自己的阶段目录 `jobs/<task id>/` 中运行，
输入为编译清单、可执行文件、`grid_mks.in`、`definitions.h`、快照和表。运行根目录不得位于某个阶段
目录之内（任何上级目录都不得含有 `manifest.json`）。RAPTOR 在一个临时目录中以
`RAPTOR model.in <snapshot> <output slot>` 的形式调用；若它成功退出，则校验其输出（下一节）并发布
该阶段。失败的尝试保存在 `attempts/<task id>/attempt-*/` 下（日志、回执和不完整的工作目录），
从不发布。重新运行时会续跑：只有当一个已发布作业的清单、成功回执、已校验的输出元数据、`model.in`、
任务回执和网格文件全部一致时才复用它。取消运行会终止所有进程组并保留它们的尝试记录。

### 作业输出校验

`validate_job_output`（`eht_axion/raptor/outputs.py`）检查一个作业：

- HDF5 文件 `output/img_data_<slot>.h5` 在单一频率下恰好包含数据集 I、Q、U、V、tau、tauF、alpha、
  beta，每个的形状为 $((N/10)^2, 100)$（10×10 像素的块），只用硬链接、不用外部存储；
- 每个块的相机坐标 `alpha`、`beta` 等于预期的像素中心（块和像素均按 x 优先顺序）；
- 谱文件 `output/spectrum_<slot>_<inclination>.dat` 只有一行（频率以及以 Jy 为单位的 I、Q、U、V），
  其流量等于 HDF5 强度（Jy/sr）乘以像素立体角 $\Delta\Omega = (\Delta x\,\Delta y)\,(r_g/D)^2$，
  所用常数取自编译目录的 `definitions.h`，相对容差为 $10^{-10}$。

结果（`output.json`）还记录 $I<0$ 的像素数和 $\sqrt{Q^2+U^2+V^2} > I$ 的像素数作为诊断量；
它不做任何科学上的接受判断。

### 格式转换

`convert_hdf5`（`eht_axion/raptor/convert.py`）读取一个由作业 HDF5 文件组成的矩形“时间 × 状态”
网格，写出一个形状为 $(T, S, N, N, 4)$ 的 float64 数组（Stokes I、Q、U、V，单位 Jy/sr；第 2 轴为
相机 x，第 3 轴为相机 y）。每个 10×10 的块放在其记录的原点处；这些块必须恰好铺满整幅图像。

### 交付与合并

`deliver_stokes`（`eht_axion/raptor/delivery.py`）把每个计划中的任务恰好分配给一个作业根目录，
校验每个作业阶段及其来源链，检查每个作业是否声明了其任务的编译和源输入，并检查作业的成功回执。
作业的 HDF5 通过作业清单中的 SHA-256 值与已校验的 `output.json` 绑定，不再重新校验。然后它写出：

- 对每个 $\alpha$：`res-100/alpha-<α>/stokes.npy`，形状为 (141, 17, 100, 100, 4)；状态 0 是无轴子
  基线，状态 $k+1$ 的轴子相位为 $k\pi/8$；
- 对无轴子背景：`res-400/background/t<A>-t<B>/stokes-t<A>-t<B>.npy`，每个分片最多 10 个快照，
  形状为 $(T, 1, 400, 400, 4)$。

每个产品目录含有 `generated.json`：光线追踪配置、$\alpha$（背景为 null）、分辨率、快照标签、以
$t_g$ 为单位的尾部时间、状态数、$\Omega$、轴顺序和单位。只有当这条记录等于当前 `experiment.yaml`
会交付的内容（`generated_inputs`）时，分析部分才接受该产品；帧间隔取自记录的尾部时间
（`input_cadence_hours`）。

`raptor_merge.py` 用多个 `--job-root` 参数和一个新的 `--output-root` 调用同一个交付函数，校验冻结
计划的光线追踪配置，并写出一个列出各产品的清单（inventory）阶段。原有的作业目录不会被修改。

### CPU 亲和性

`allocate_cpu_sets`（`eht_axion/raptor/affinity.py`）给每个工作进程分配一组互不相交的、共
`threads` 个 CPU。`physical` 每个物理核只用一个逻辑 CPU（依据
`/sys/devices/system/cpu/*/topology`），`logical` 使用任何允许的 CPU，`none` 不做绑定。绑定使用
`taskset`；`--cpu-list` 限制允许使用的 CPU。启用绑定时，CPU 集合和 OpenMP 绑定方式属于作业身份的
一部分。`check_memory` 要求每个工作进程有足够的内存余量（至少 2.5 GiB，在 `raptor_generate.py`
中还要求至少为计划给出的每进程下限的两倍），并考虑 cgroup-v2 限制。

### .env.raptor

`raptor_build.py` 和 `raptor_generate.py` 读取 `.env.raptor` 和 `.env`（由 `.env.raptor.example`
复制而来）。键：`EHT_RAPTOR_SOURCE_ROOT`、`EHT_GRMHD_ROOT`、`EHT_RADIAL_LOW`、`EHT_RADIAL_HIGH`、
`EHT_RAPTOR_RUN_ROOT`、`EHT_RAPTOR_THREADS`（默认 2）、`EHT_RAPTOR_WORKERS`（默认 1）、
`EHT_RAPTOR_CPU_AFFINITY`（默认 `none`）、`EHT_RAPTOR_CPU_LIST`、`EHT_RAPTOR_TIMEOUT_SECONDS`。
源码根目录和运行根目录没有默认值。运行根目录不得与源码目录、流体数据目录或径向表目录重叠；一个
运行根目录始终使用它启动时的实验文件（`config.yaml`），换用不同的文件需要新建一次运行。

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [实验配置](experiment-config.md)
- [来源记录](provenance.md)
- [数据布局](data-layout.md)
- [开发](development.md)
