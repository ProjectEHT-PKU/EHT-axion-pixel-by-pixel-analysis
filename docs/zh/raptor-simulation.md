# RAPTOR 模拟

[← 返回文档索引](index.md)

`RAPTOR` 模拟流水线编译固定版本的偏振广义相对论光线追踪（ray tracing）核心程序 `RAPTOR`，校验广义相对论磁流体力学（general-relativistic magnetohydrodynamics, GRMHD）BHAC 快照和轴子径向轮廓表，以有界并发和可选 CPU 绑定方式运行 $18\,330$ 个内容寻址（content-addressed）光线追踪作业，逐一校验 HDF5 与谱文件，并为下游分析交付 float64 精度的斯托克斯（Stokes）数据立方。本部分不依赖任何 EHT 成像产出。

## 前置条件与输入数据

运行 `RAPTOR` 脚本需要 Linux（`linux-64`）系统及 `raptor-build` Pixi 环境（含 `c-compiler`、`make`、`gsl`、`hdf5`、`pkg-config`、`rsync`）。共需四项外部输入（通过 `.env` 中的三个根目录 `EHT_RAPTOR_SOURCE_ROOT`、`EHT_INPUT_ROOT` 与 `EHT_DATA_ROOT` 定位，或用命令行选项显式指定；见[配置参考](configuration.md#本机环境变量env)）：

- **`RAPTOR` 源码工作副本**（`--source-root` / `EHT_RAPTOR_SOURCE_ROOT`）：一个洁净的 Git 工作副本（不含已修改或未跟踪文件），其 `HEAD` 恰好处于 `raptor.source.commit`（`20056006820b6536c96de790dd5fdaad5ada4eae`），且至少有一个规范化后等于 `raptor.source.repository`（`github.com/projecteht-pku/raptor`）的 Git 远程仓库地址。
- **BHAC GRMHD 快照**（`--fluid-root`，默认 `$EHT_INPUT_ROOT/grmhd`）：二进制快照文件 `data<label:06d>.dat`（从 `data001600.dat` 到 `data003000.dat`，步长为 `10`，共 141 个快照）。解析后的文件路径必须是不含空白字符且不超过 255 个 UTF-8 字节的单一路径标记。
- **径向轮廓表**（`--radial-low`、`--radial-high`，默认 `$EHT_INPUT_ROOT/radial/result_a=0.9375_l1m1.txt` 与 `$EHT_INPUT_ROOT/radial/result_a=0.99_l1m1.txt`）：分别对应 `low`（`spin = 0.9375`）与 `high`（`spin = 0.99`）的两张无表头三列 CSV 文件 $(\alpha, r, R)$，路径同样不超过 255 个 UTF-8 字节且不含空白字符。
- **实验参数文件**（`--config`，默认 `configs/experiment.yaml`）：读取 `shared` 物理标量及 `raptor` 节。

运行根目录（`--run-root`，默认 `$EHT_DATA_ROOT/raptor/runs/main`）不得与本仓库根目录、`RAPTOR` 源码目录、快照目录或径向表目录重叠，也不得位于任何已有阶段目录之内。

## 运行方法

将 `.env.example` 复制为 `.env` 并设置三个根目录 `EHT_RAPTOR_SOURCE_ROOT`、`EHT_INPUT_ROOT` 与 `EHT_DATA_ROOT`。完整流程如下：

```bash
# 1. Compile the pinned RAPTOR binary into $EHT_DATA_ROOT/raptor/runs/main/build
pixi run -e raptor-build raptor-build

# 2. Validate inputs and freeze the task plan without launching C processes
pixi run -e raptor-build raptor-generate --plan --plan-output $EHT_DATA_ROOT/raptor/plan.json

# 3. Execute ray-tracing jobs (resumable) and deliver Stokes cubes in the run root
pixi run -e raptor-build raptor-generate

# Deliver cubes from already completed jobs in one run root without launching RAPTOR:
pixi run -e raptor-build raptor-generate --convert-only

# 4. Merge completed job roots from one or more hosts/runs into $EHT_DATA_ROOT/raptor/delivered
pixi run -e raptor-build raptor-merge
```

若要显式指定多台主机的作业根目录进行合并（而不扫描 `$EHT_DATA_ROOT/raptor/runs/*`）：

```bash
pixi run -e raptor-build raptor-merge --plan /path/to/plan.json \
    --job-root /path/to/runs/host-a --job-root /path/to/runs/host-b \
    --output-root /path/to/delivered
```

**`scripts/raptor/generate.py` 的执行选项：**
- `--max-jobs N`：仅运行计划中的前 `N` 个任务（仍有未完成任务时跳过数据立方交付）。
- `--workers W`（`EHT_RAPTOR_WORKERS`，默认 `1`）与 `--threads T`（`EHT_RAPTOR_THREADS`，默认 `2`）：并发运行的 `RAPTOR` 进程数及每个进程的 OpenMP 线程数。
- `--cpu-affinity {none,physical,logical}`（`EHT_RAPTOR_CPU_AFFINITY`，回退值 `none`）与 `--cpu-list RANGE`（`EHT_RAPTOR_CPU_LIST`，如 `0-47`）：使用 `taskset` 对工作进程做 CPU 亲和性（CPU affinity）绑定。
- `--timeout-seconds SEC`（`EHT_RAPTOR_TIMEOUT_SECONDS`）：单个进程的墙钟超时限制。

`scripts/raptor/build.py` 只需要 `RAPTOR` 工作副本（`--source-root`）和运行根目录（`--run-root`，或 `--data-root` / `EHT_DATA_ROOT`）；它既不要求也不解析输入根目录，因此纯编译机器无需挂载快照或径向表。同一数据根目录下的第二个运行根目录（如 `raptor/runs/background`）通过 `--run-root` 指定。

## 编译固定版本的 RAPTOR 核心

`build_raptor`（`eht_axion/raptor/build.py`）在 `<run_root>/build/` 发布阶段 `raptor:build`：

1. 校验 `--source-root` 为处于 `raptor.source.commit` 的洁净 Git 仓库，且至少有一个远程地址规范化后等于 `raptor.source.repository`。
2. 提取受版本控制的源码字节（`git ls-tree` / `git cat-file`）到 `build/source/`，清理外部 `make` 和头文件环境变量（`VPATH`、`CPATH`、`MAKEFLAGS`、`RAPTOR` 等），并解析编译器 `h5cc`（`HDF5_CC`）。
3. 按 `raptor.build` 派生的参数（`--code`、`--metric`、`--int`、`--grid`、`--rad`、`--sf`、`--bflip`、`--sigma-cut`、`--uniform`、`--thetae-min`、`--thetae-max`）运行官方 `setup.sh`。
4. 通过 `make -nB all` 记录编译器命令行（`compiler-flags.txt`），通过 `h5cc ... -dM -E` 导出生效的预处理器宏（`macros.txt`），并断言 `metric`、`int_method`、`AMR`、`SMR`、`POL`、`SFC`、`BPOL`、`SIGMA_CUT`、`UNIF`、`THETAE_MIN` 和 `THETAE_MAX` 与 `raptor.build` 一致。
5. 写出 `receipt.json`，记录工具链信息、生效宏取值以及 `RAPTOR` 可执行文件的 SHA-256（每个下游作业都将其声明为输入）。

## 任务规划与内容寻址标识符

`plan_generation`（`eht_axion/raptor/plan.py`）仅检查轻量级元数据，不启动任何光线追踪计算：

- **快照尾部检查**（`read_snapshot_footer`）：读取每个 `data<label:06d>.dat` 末尾 40 字节的小端序 BHAC 尾部及紧邻其前的 `nx`、`neqpar` 记录；核对三维块结构、块尺寸 `nx` 整除 `raptor.grid.nxlone`、黑洞自旋等于 `raptor.fluid_spin`（`0.9375`，绝对容差 $10^{-12}$）、文件字节长度不小于尾部隐含的数据载荷下界、绝热指数 $\gamma$ 恒定，以及相邻快照尾部时间差等于 `shared.label_step_tg * step`（每 10 个标签间隔 $1.0\,t_g$）。
- **径向表检查**（`validate_radial_table`）：载入每张无表头 CSV，要求有限的 $(\alpha, r, R)$ 且 $R \ge 0$，各行按严格递增的 $\alpha$ 分块排列并共用同一严格递增的 $r$ 网格（$\ge 2$ 个点），核对其覆盖 `raptor.alpha_table` 分配给它的全部 $\alpha$ 值，并计算文件 SHA-256（不对振幅做重新归一化）。

对 141 个快照中的每一个，计划生成 130 个任务（总计 $18\,330$ 个任务）：
- **2 个 `baseline` 任务**（无轴子：`AXION_ALPHA = 0`、`AXION_OMEGA = 0`、`AXION_NORM = 0`、`AXION_PHASE = 0`、`AXION_TABLE = disabled`、`output_slot = 0`），分辨率分别为 `response.resolution`（`100`）和 `background.resolution`（`400`）；
- **128 个 `response` 任务**（分辨率 `100` 像素），遍历全部 8 个 $\alpha$ 取值与 16 个相位序号 $k \in \{0, \dots, 15\}$（`AXION_ALPHA` $= \alpha$、`AXION_OMEGA` $= \Omega(\alpha) = \alpha(1 - \alpha^2/(2n^2))$、`AXION_NORM` $= 0.1$、`AXION_PHASE` $= k$ 对应核心相位 $k\pi/8$、`output_slot` $= k + 1$，并使用 `raptor.alpha_table` 指定的径向 CSV 表）。

每个任务携带全部 19 个 `model.in` 参数（以 $M_\odot$ 为单位的 `MBH`、以 $\mathrm{kpc}$ 为单位的 `DISTANCE`、以 $\mathrm{g}$ 为单位的 `M_UNIT`、`R_LOW`、`R_HIGH`、以 $\mathrm{deg}$ 为单位的 `INCLINATION`、`IMG_WIDTH`、`IMG_HEIGHT`、以 $r_g$ 为单位的 `CAM_SIZE_X` 与 `CAM_SIZE_Y`、`FREQS_PER_DEC = 1`、以 $\mathrm{Hz}$ 为单位的 `FREQ_MIN`、`STEPSIZE`、`AXION_ALPHA`、`AXION_OMEGA`、`AXION_NORM`、`AXION_PHASE`、`AXION_TABLE`、`MAX_LEVEL = 1`）。

**任务标识符（`task_identity`）：** 规范化 JSON 文档的 SHA-256，包含快照文件名及尾部元数据、除 `AXION_TABLE` 外的 19 个 `model.in` 取值、径向表的 `{spin, sha256}`（`baseline` 为 `null`）、`raptor.grid` 和 `raptor.build`。天空网格、EHT 成像和分析参数均不进入任务标识符，因此修改下游设置永远不会使已完成的光线追踪作业失效。

## 作业执行、CPU 亲和性与输出校验

在启动作业前，`scripts/raptor/generate.py` 与 `run_jobs`（`eht_axion/raptor/jobs.py`）会：
- 检查剩余磁盘空间（$\ge 4 \times \text{output\_stokes\_lower\_bound\_bytes} + 2\,\mathrm{GiB}$）与可用内存余量（$\ge \text{workers} \times \max(2.5\,\mathrm{GiB},\,2 \times \text{per\_process\_memory\_lower\_bound\_bytes}) + 2\,\mathrm{GiB}$，同时考虑 Linux `MemAvailable` 与 cgroup-v2 的 `memory.max - memory.current`）；
- 校验 `raptor:build` 阶段，核对 `raptor_source`、`grid_mks.in` 及必需的 C 宏（`POL=1`、`num_frequencies=1`、`num_pixels_1d=10`、`NSPIN=3`、`SMR=1`、`AMR=0`、`UNIF=0`、`metric=MKSBHAC`、`DF=TH`、`IMGFILE=1`、`SPECFILE=1`），并重新校验实时快照尾部与径向表哈希；
- 在单个 `verification_session` 中对共享输入闭包预先计算一次哈希。

**调度与 CPU 亲和性（`eht_axion/raptor/affinity.py`）：** 第 `i` 个任务确定性地分配给工作槽位 `i % workers`。`allocate_cpu_sets` 为每个工作进程分配互不相交的 `threads` 个 CPU（`physical` 依据 `/sys/devices/system/cpu/*/topology` 在每个物理核 `(socket, core)` 上只选一个逻辑 CPU；`logical` 使用任意允许的逻辑 CPU；`none` 不做绑定）。启用绑定时，进程通过 `taskset --cpu-list ...` 启动并设置 `OMP_PROC_BIND=close`、`OMP_PLACES=cores`，这些绑定参数会写入 `raptor:job` 的阶段配置，保证断点续跑时缓存身份完全一致。

**执行、续跑与失败现场保留：** 每次实际执行都会创建 `<run_root>/attempts/<task_id>/attempt-*/`（记录 `input.json`、`stdout.log`、`stderr.log`、`result.json`），并在临时阶段目录中调用 `RAPTOR model.in <snapshot> <output_slot>`。若 `RAPTOR` 非零退出、超时或被中断（`SIGINT`/`SIGTERM` 会终止所有活跃进程组），未完成的临时工作目录会被复制到 `attempts/<task_id>/attempt-*/work/`，不会发布为完成阶段。重新运行时，已存在的 `jobs/<task_id>/` 仅在其阶段清单（manifest）、`success.json`、`output.json`、`model.in`、`task.json` 和 `grid_mks.in` 全部校验通过时才会被复用。

**单作业输出校验（`eht_axion/raptor/outputs.py` 中的 `validate_job_output`）：**
- 检查 `output/img_data_<slot>.h5` 仅包含内部硬链接数据集 `{I,Q,U,V,tau,tauF}<freq:e>`、`alpha` 和 `beta`，每个数据集均为有限值且形状为 `((N/10)^2, 100)`（按 x 优先排列的 $10 \times 10$ 像素块）；
- 核对每个块的相机坐标 `alpha`、`beta` 与预期像素中心（单位 $r_g$）一致；
- 核对单行谱文件 `output/spectrum_<slot>_<inclination:.2f>.dat`（`frequency` 及以 $\mathrm{Jy}$ 为单位的 $I, Q, U, V$）与 $\sum \text{HDF5}_{IQUV}\,\Delta\Omega$ 在相对容差 $10^{-10}$ 内一致，其中像素立体角 $\Delta\Omega = (\Delta x\,\Delta y)(r_g / D)^2$ 使用从 `build/definitions.h` 解析出的字面常数 `GGRAV`、`MSUN`、`SPEED_OF_LIGHT` 和 `KPCTOCM`；
- 在 `output.json` 中记录 $I < 0$ 或 $\sqrt{Q^2 + U^2 + V^2} > I$ 的像素数作为诊断信息，不做数据过滤。

## Stokes 数据立方转换、交付与多机合并

`convert_hdf5`（`eht_axion/raptor/convert.py`）将矩形 `time × state`（时间 × 状态）网格上的作业 HDF5 分块拼接为形状 `(T, S, N, N, 4)` 的 float64 数组（`[time, state, x, y, IQUV]`，单位 $\mathrm{Jy\,sr^{-1}}$），并写出 `convert.json`。

`deliver_stokes`（`eht_axion/raptor/delivery.py`）在一个或多个作业根目录中将每个计划任务唯一对应到一个已完成作业阶段，校验其来源链与 `success.json`，并在 `<target>/data/raptor/` 下发布 `raptor:delivery` 阶段：
- `res-100/alpha-<alpha:.2f>/stokes.npy`，形状为 `(141, 17, 100, 100, 4)`（共 8 个阶段：状态 `0` 为无轴子基线，状态 `k+1` 对应轴子相位 $k\pi/8$）；
- `res-400/background/t<A>-t<B>/stokes-t<A>-t<B>.npy`，按最多 10 个连续快照分片（共 15 个阶段：`t1600-t1690` … `t3000-t3000`），形状为 `(T, 1, 400, 400, 4)`。

每个交付阶段同时写出 `generated.json`（模式 `raptor-generated`，记录 `ray_tracing_config`、`alpha`、`resolution`、`snapshot_labels`、BHAC 尾部时间 `times_tg`、`state_count`、`omega`、坐标轴与单位）和 `sources.json`（任务到源作业根目录的映射）。下游分析阶段通过 `generated_inputs` 强制核对 `generated.json` 与当前 `configs/experiment.yaml` 完全一致。

当作业分散在多台主机或多个运行目录中完成时，`scripts/raptor/merge.py` 通过 `deliver_stokes` 将各作业根目录合并交付到互不相交的新目录 `--output-root`（默认 `$EHT_DATA_ROOT/raptor/delivered/`），并写出 `inventory/index.json`（阶段 `raptor:merge`，布局 `raptor-merge`），不修改任何源作业目录。
