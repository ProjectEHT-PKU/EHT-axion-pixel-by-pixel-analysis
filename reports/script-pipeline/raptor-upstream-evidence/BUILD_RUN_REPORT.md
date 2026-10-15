# RAPTOR 构建、运行与源码身份：独立只读审查

## 范围、证据与限制

本报告专注构建入口、运行配置、表选择、debug/实验变体、源码/二进制 provenance；不重复核心算法和 fluid 模型的完整物理审查。

路径简写：

- `L` = `/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel`
- `G` = `/tmp/eht-github-raptor-GckNId`；core 在 `G/raptor`，项目运行目录在 `G/run/raptor_run`
- `U` = `/tmp/eht-raptor-upstream-audit-lVJjW7/upstream`
- `E` = `/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT`

四处主要审查对象为 `L/raptor`、`L/raptor_run`、`G/raptor`、`G/run/raptor_run`，另以完整上游 git 为参照。补充检查严格限定在明确指定的 `L/raptor_run_test`、`L/raptor_run_test_eq_fakevalue`、`E/raptor_run`、`E/raptor_run_gkx`、`E/test_raptor_difference` 的必要入口、头文件、少量 src 和保存对象/二进制；未扩大扫描 vmdk。

没有修改仓库、挂载数据或远端；没有执行任何保存的 RAPTOR binary、下载数据或重跑科学计算。仅在 `/tmp` 写报告/证据，执行独立的 read_model 解析探针、GCC 预处理和 toy make dry-run。当前审查环境未安装 `h5cc`，也未设置 `RAPTOR`，所以**没有完成 HDF5/GSL 科学程序的真实编译、链接或运行验证**。

补充文件：

- `/tmp/raptor-build-run-audit-probe.py`
- `/tmp/raptor-build-run-audit-probe.log`
- `/tmp/raptor-build-run-probe-hmwlewy9/`（解析探针、空 hdf5.h 预处理占位、toy make 环境；不是科学构建）
- `/tmp/raptor-build-run-static-evidence.json`（必要文件 SHA256 与 get_fluid_params 反汇编）
- `/tmp/raptor-main-vs-test-pol-static.diff`

## 一、结论与重跑阻断项

1. **主 fiducial 当前源码没有启用“关闭 plasma effects／仅赤道面／假发射值”实验块。** 这些有意义的注释是实验设计线索，不是活跃代码。不能把它们当成主模型正在关闭吸收和 Faraday 效应的证据。
2. **本地主 core 模板 `SIGMA_CUT=1000` 不等于主运行目录/主 binary 使用 1000。** `L/raptor_run/definitions.h` 是 1；按实际 include 路径预处理仍是 1。静态反汇编也确认保存的主 binary 的磁化截断比较常量为 1.0。
3. **高 α 表不自动选择。** 活跃源码只 `fopen("NormR.csv","r")`。主 `NormR.csv` 覆盖 0.09–0.31；同目录的未跟踪 `NormR_0.99.csv` 覆盖 0.30–0.44、只有 15 个块，而当前 reader 固定读 23×2400。直接添加 α=0.35/0.4/0.44 或仅把高 α 文件改名，不是安全的重跑方案。
4. **项目自带的 `model.in` 与当前 axion reader 不兼容。** 它缺三项 axion 字段；解析会把末尾 `MAX_LEVEL=1` 吞成 axion_norm=1，余下 omega/phase/max_level 留为零。不能直接沿用 README/run.sh 作为本项目的可复现生产命令。
5. **make 不追踪头文件或编译环境变化。** 修改 `definitions.h`、切换 RAPTOR 路径、改变编译参数后，普通 `make all` 可以保持旧对象/旧 binary。必须全新隔离构建并保存 compile manifest。
6. **批处理不检查 RAPTOR/绘图退出状态，断点判断只看 H5 文件存在。** 失败后可能移动旧输出、误写 done；并发进程共享 model.in/output/figures 必然有冲突风险。旧结果目录不可作为新配置的缓存继续复用。

## 二、四处身份与上游关系

| 对象 | HEAD / 状态 | 与其他处关系 |
|---|---|---|
| `L/raptor` | `b6df33c69e650f3617845f4324ba642bdb878582`；dirty 两文件 | `model/bhac/definitions.h` 从 sigma 1 改为 1000；`src/pol_rte_integrator.c` 字节不同，但去注释后词法 token 与 G core 完全相同（7669 token） |
| `G/raptor` | 同一 `b6df33c`；clean | 可作为该项目提交的干净源码基线，不是 upstream main |
| `L` / `L/raptor_run` | 项目 HEAD `7d3ecac555e325e5ab6c160bfa00ff2af15e578e`；主 Run.py 改动，model.in 已删除 | Run.py 实际只扫描 `[0.1]`；其他跟踪 runtime 的 makefile、run.sh、headers/model.c、grid、NormR.csv、plotter 与 G runtime 相同 |
| `G/run/raptor_run` | 项目同一 `7d3ecac`；clean | Run.py 扫 `[0.3,0.25,0.2,0.15,0.1]`；含旧式 model.in；submodule gitlink 精确固定到 `b6df33c` |

公开仓库身份：项目为 `ProjectEHT-PKU/EHT-axion-pixel-by-pixel`，core 为 `ProjectEHT-PKU/raptor`；上游为 `jordydavelaar/raptor`。仅输出了公开 owner/repo，未输出含认证信息的原始 remote 配置。

上游 git：

- upstream main = `b91e937b7b746eb24992b6c1dd8831c5dea0cef1`
- project-main = `b6df33c69e650f3617845f4324ba642bdb878582`
- merge base = `cdc2547b868ed8d6ae0c9a5428bc1ea86fce901a`
- `main...project-main` 分别有 13 / 5 个独有提交；**双方分叉，项目不是当前上游 main 加几个补丁后的直接后代。** 项目提交时间较新不表示吸收了上游修复。

本领域控制文件差异：项目相对 mergebase 的 makefile 增加 `-lgslcblas`；项目 setup.sh 与 mergebase 相同，而 upstream main 增加了 `-n=*|--nfreq=*` 处理。README 在上游新增 spectra 教程；项目 README 基本保留原版，未说明 axion 三字段、NormR 数据依赖和 wrapper。上游 main 没有本项目的 axion 表工作流；不能用上游二进制代替项目版来“原样复现”。是否吸收上游 emission/fluid 修复，应与另外两份算法/模型审查合并后决定。

## 三、实际构建：source/header 如何选择

### 3.1 makefile 的确切规则

`raptor/run/makefile:1–28` 与 runtime makefile 相同：

```make
CC = h5cc -I$(RAPTOR)/src -I$(PWD)
VPATH=$(RAPTOR)/src:
CPATH=$(RAPTOR)/src:
SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c pol_emission.c tetrad.c model.c constants.c camera.c
$(OBJECTS): $(OBJDIR)/%.o: %.c
```

- `.c` prerequisite 先找工作目录，再由 VPATH 找 `$RAPTOR/src`。主 runtime 目前只有 `model.c` 这一个本地 .c；因此它用本地 BHAC model.c，其他 12 个 .c 来自所选 RAPTOR core 的 src。运行目录放入同名 core.c/pol_rte_integrator.c 会覆盖 VPATH 源码，不能只按 core git HEAD 判身份。
- 所审查主 core 的 `src` 没有 definitions.h/model_*.h，只有 functions.h/global_vars.h。因此 `definitions.h`、model_definitions.h/model_functions.h/model_global_vars.h 来自 runtime。`model/bhac/definitions.h` 是 setup 的模板，不是这次 make 的 include 路径。
- `#include "…"` 先查包含它的文件所在目录，再查 `-I`。core 文件通常用 core/src 的 functions.h/global_vars.h，找不到的 runtime 头通过第二个 `-I` 找到；本地 model.c 先找自身所在 runtime 的头。若另一个 RAPTOR/src 有同名配置头，它可能优先于第二个 `-I`，甚至形成不同 translation unit 配置不一致。
- GNU make 的 VPATH 与 C 预处理器的 `-I` 不是同一机制。makefile 中 `CPATH=...` 并不能无条件保证给编译器导出；若外部已有 exported CPATH，它又可能参与标准 include 搜索。重跑应显式记录/清理 CPATH、C_INCLUDE_PATH 等相关变量，而不是打印全部环境。
- `$(PWD)` 依赖环境/命令行，绝不是写死的 runtime；`make -C` 或外部设置可能让它与实际工作目录不同，建议改成基于 `CURDIR` 的显式绝对 include 路径。命令行还可覆盖 RAPTOR、CC、CFLAGS 等。

已用 `gcc -E -H -dM` 对挂载主 core.c 与 runtime model.c 做只读预处理，保留同一 `-I core/src -I runtime` 顺序，结果：

```text
SIGMA_CUT (1)
metric (MKSBHAC)
POL (1)
DEBUG (1)
LIGHT_TRANSPORT (1)
RAD_TRANS (1)
```

仅为预处理补了 `/tmp` 的空 hdf5.h，占位不代表可链接/可运行。

### 3.2 漏依赖与编译并发

- 对象只依赖 `.c`，没有 `.h` prerequisite，也没有 `-MMD -MP` 和 `.d` include。只改 sigma/header，`make all` 不重编译。toy probe 将头时间设为晚于 binary，仍只打印目录检查。
- 更换 RAPTOR 到内容不同但时间更旧的完整 core，toy probe 仍不重编译。make 不跟踪源码内容 hash、CC/CFLAGS/-I/VPATH 或链接库变化。
- create_directories 是 all 的并列 prerequisite，对象没有 order-only 依赖它；干净环境用 `make -j` 可能在 build/ 创建前编译。也没有 `.PHONY` 防止同名文件压掉辅助目标。
- link 命令把 `-lm -lgsl -lgslcblas` 放在 objects 前，对 static archive 或 `--as-needed` 情况并不稳健；旧 binary 能链接不证明新环境能链接。应把库作为 LDLIBS 放在 objects 后。
- 保存 binary 的 RPATH 指向 `/home/yuxin/anaconda3/envs/eht/lib`，动态依赖包括 libgsl.so.27/libgslcblas.so.0/libgomp.so.1 等；挂载路径不是自动重定位后的运行环境。未执行 ldd 或 binary。

挂载 home 的 `.bashrc`、`.profile` 仅做 RAPTOR 相关定向检查，大小写不敏感也没有命中 RAPTOR；不能从这两个文件恢复历史 RAPTOR export，也不能排除历史 shell/job/其他脚本临时设置。当前审查 shell 的 RAPTOR 未设置。

## 四、setup.sh/readme/run.sh：不能当成可靠的生产配置来源

`setup.sh` 的实际默认值是 `cks / rk4 / bhac / smr / pol / nosfc`。主 runtime 则是 `MKSBHAC / RK4 / POL=1 / SMR=1 / SFC=1 / BPOL=PLUS / NSPIN=3`。

- 裸跑 setup 默认改成 CKS、SFC=0，**不是当前 MKS+SFC1 fiducial 的原样重建**。应先确认数据格式，通常重建当前配置需要显式 mks/bhac/pol/sfc；不能只因为 setup 默认如此就改变物理配置。
- 两次 `rsync -u` 按 mtime 保留较新的 runtime 文件，不按 git/hash；可能保留旧实验 model.c/header。`rsync` 的日期顺序不等于 provenance。
- setup 最后会 `make clean && make all`，所以它会清理对象，但清理前到底复制了哪个 source/header 仍不确定。其行为与普通 `make all` 不同。
- README 写 `-g=amr/smr`，脚本判断却是大写 `AMR/SMR`；默认小写 smr 也不触发替换，仅碰巧模板本来就是 SMR1。`-g=amr` 可能仍得到 SMR，且日志照样说使用 amr。此问题项目与当前上游都存在。
- README 写 HAMR，脚本接受小写 hamr；参数实际必须是 `-m=mks`/`--metric=mks` 这种带等号格式。未知/不合法参数多数静默忽略，没有完整 validation。
- 未设置 RAPTOR、model 选错、编译失败等分支在 echo 后执行裸 `exit`，通常退出码为 0；copy/sed 失败也不统一 fail-fast。日志“成功”或 shell 退出 0 不是配置/编译成功的充分证据。
- `setup.sh` 不给项目 runtime 安装 Run.py 或 NormR.csv；纯 core 新目录按旧教程初始化后仍没有轴子表。
- 所有主 run.sh 都仍是旧 Slurm 样例：`./RAPTOR model.in harm3d.txt 0`。它只设 `OMP_STACKSIZE=20m`，没有显式 OMP_NUM_THREADS；`-n 24` 是 task 数请求，不等于可靠的每任务 24 OpenMP 线程配置。推荐单 task、明确 cpus-per-task/线程/绑定，并按内存测试调整 stack。
- 本地主 `20231114.sh` 是未跟踪旧流水线：先 Run.py、再 analy_grmhd.py、再 cd 另一个 M87-image/eht 路径；没有 `set -e`/失败停止。其“gen img_data_0_high_res”注释与当前 `RUN_HIGH_RESOLUTION=False` 不一致，不能把注释当执行结果。
- README 对 MBH 写 gram，但 core.c:102 实际 `MBH *= MSUN`；运行输入为太阳质量。DISTANCE 实际输入 kpc，乘 KPCTOCM。README 没解释 axion 参数。

## 五、Run.py/model.in 的实际优先级、命令和表

### 5.1 当前活跃配置

`G/run/raptor_run/Run.py` 与 `L/raptor_run/Run.py` 唯一跟踪内容差异是扫描 omega list：G 为 `[0.3,0.25,0.2,0.15,0.1]`，L 为 `[0.1]`。字典中的 `AXION_OMEGA=0.4` 被 main 循环覆盖，不是当前生产批处理的有效默认 α。

其余活跃配置：MBH=6.2e9 Msun、DISTANCE=16800 kpc、M_UNIT=3e25 g、R_LOW=R_HIGH=1、inclination=163°、100×100、CAM_SIZE=46.2261 Rg、230 GHz、stepsize=.005、MAX_LEVEL=1；每 snapshot 先 norm=0 baseline，再 norm=.1、phase=0…15 共 16 个相位；高分辨率开关 False。

优先级：

1. 编译期 definitions/model headers 决定 metric/POL/SIGMA/SFC/num_frequencies 等；model.in 不能覆盖它们。
2. 批处理 main 的循环覆盖 PARAMS_JSON 中 omega 等字段，然后 `Run` 每次创建/覆盖当前工作目录 model.in。
3. core/read_model 按固定位置读取 17 个三元组，不校验标签；第三个 CLI 参数读到 TIME_INIT，用作输出 index，不是轴子相位。
4. 轴子相位来自 model.in，活跃公式为 `phase_index/8 * pi`，因此 16 相位均匀覆盖一周期。baseline 输出编号0；相位 k 输出编号 k+1。
5. Run 不编译、不检查 RAPTOR 环境；它只执行工作目录的 `./RAPTOR`。运行时换 RAPTOR export **不会**改变这个已有 binary，但下次 make/setup 会改变源码来源。

命令确切为：

```text
./RAPTOR model.in <absolute-GRMHD-file> 0
./RAPTOR model.in <absolute-GRMHD-file> <phase+1>
python plotter-example.py <same-output-index>
```

Run.py 没有 CLI 配置解析，也不主动 cd 到脚本所在目录；从错误 cwd 调它会读/写错误的 model.in、NormR.csv、output/，或找不到 binary。硬编码 snapshot 路径仍是旧 `/home/yuxin/...`，在审查宿主该路径不存在，挂载数据有另一个前缀。必须明确重跑用原系统、容器内路径还是宿主挂载路径。

文件选择是 `sorted(os.listdir)` + 所有普通文件 + `[100:]`，不是按 simulation time 筛选，也不验证 .dat。此刻定向目录元数据为 241 个 .dat，截掉100后为141个，首个 data001600.dat，末个 data003000.dat。若目录插入其他文件，选择集合会静默改变；应固定 snapshot manifest。

### 5.2 model.in 的两个重要区分

- Run 的 f.write **没有换行**，但不能据此宣称它读错。数字后紧跟下一个标签，现有 fscanf 数值转换会在字母处停止，再由后续 `%s` 读标签。用抽出的真实 read_model 和重建写入表达式做 `/tmp` probe，100像素/norm=.1/omega=.3/phase15/maxlevel1 全部正确。仍建议改成每项一行并严格校验，避免依赖这种脆弱分词。
- 相反，仓库静态 `raptor/run/model.in` 和 `G/run/raptor_run/model.in` 只有上游旧14项，缺 AXION_NORM/OMEGA/PHASE。同一 probe 确认结果是 norm=1、omega=0、phase=0、maxlevel=0；不是安全的“无轴子默认值”。`L/raptor_run/model.in` 已删除不妨碍 wrapper，因其每次生成后删除，但使旧 run.sh 更不能直接使用。
- main/read_model 不验证 argc、不检查 fscanf/sscanf 返回值，`%s` 无宽度；短参数、带空格或超长路径有失败/截断/溢出风险。argv3 用 double 读入，文件命名再 cast int，不是严格整数验证。

### 5.3 高 α 表和边界

活跃 reader 在 `src/pol_rte_integrator.c:654–705`：固定 `omega_count=23`、`r_count=2400`，只打开 CWD 的 `NormR.csv`；不根据 alpha、spin、文件名或环境自动选表。Run.py/makefile/setup.sh/run.sh 也没有高 α 选择逻辑。

| 表 | 内容 | SHA256 |
|---|---|---|
| 主/G runtime NormR.csv | 55200行；23块，α=.09,.10,…,.31，每块2400；r=1…240.9 | `0a352558f8b846f2481ab2acf35c2d6f8b30ce1f83f2e6f5196660f23433c353` |
| `L/raptor_run/NormR_0.99.csv`（未跟踪） | 36000行；15块，α=.30,.31,…,.44，每块2400；r相同 | `0ff8b6c816c997f32ef4d2029f52fd6cabd22990dab85362a54e35cbc62f3815` |
| 旧 E runtime 两处 NormR.csv | 38×490行；第一列约 .074… .444 | `e079a959bdd763a68594b0d4014b8845a98fc8ec5ecef30c3caa90ed33017100` |

风险：

- 主表超过 .31 的 omega 搜索没有上界保护，会越界；低于最小值又会外推而非拒绝。Run 默认字典 .4 单独调用时不安全，虽然当前 main 覆盖到安全扫描集合。
- 即便把高 α 文件改名 NormR.csv，固定23块 reader 仍会读过 EOF；fscanf 未检查。这不是完整适配。
- 文件打不开、列数错误、过短、非有限值、α/r非单调都没充分验证。即使 norm=0，仍进入表初始化/索引逻辑，不能认为 baseline 不需要有效表。
- 首次读表在 OpenMP 像素循环内部，用共享 static flag/table、没有 once/critical/barrier 同步，存在多线程重复初始化、数据竞争或读到不完整表的风险。应在并行区前完成唯一一次读取/验证，再只读共享。此项是运行可靠性，不在本报告展开数值算法。
- `.99` 只是文件名线索，不能独立证明生成时 spin=.99；必须追回表生成参数与归一化定义，尤其文章 fiducial spin=.9375 / alpha=.3 与旧文 .4/.99 设置不同。

## 六、保存 binary 与实验目录：静态可确认的边界

### 6.1 sigma 不再只靠猜测

对保存 ELF 的 get_fluid_params 反汇编，定位 `(rho+1e-20)` 之后的浮点除法、紧随比较及清 n_e 分支，并由 ELF LOAD 段把常量地址映射到文件偏移解码 double：

- 主 binary：`0x412286 divsd`，`0x41228b comisd [0x5d9370]`，常量 **1.0**；大于阈值跳到 `0x412378` 清 n_e 槽。旁边三个比较常量为 RT cutoff=50、theta_e max=100/min=.001，吻合源代码结构。
- test binary：相同结构，sigma 常量 **1.0**。
- eq_fakevalue binary：同一比较引用 `0x5d8c00`，解码 **1000.0**。
- 旧 E/raptor_run 与 E/raptor_run_gkx binary：sigma 常量 **1.0**。

这是对所列 SHA256 保存文件中特定磁化截断机器码的静态确认，不是执行结果，更不是其历史所有输出的配置证明。

| 保存文件 | binary SHA256 | sigma 静态常量 |
|---|---|---:|
| `L/raptor_run/RAPTOR` | `7137ba3b79e7120f45f81a643e6a888139011680abe9456dc799920f00e6412d` | 1 |
| `L/raptor_run_test/RAPTOR` | `69294d7762c45e69dc0d8096fcaf721f38138b4cf88ba66ae5864aed67be34ee` | 1 |
| `L/raptor_run_test_eq_fakevalue/RAPTOR` | `bd738b9eed66a2e26f5ae83d53dd972815ce2200c5725952861f60e3463ca5e3` | 1000 |
| `E/raptor_run/RAPTOR` 与 `E/raptor_run_gkx/RAPTOR` | `2726c3401103dccda4f80543ac61bc50a37ca8275ec9117b07a1ec04bd4ed721` | 1 |

### 6.2 主源码和真正活跃实验不能混同

- 主/G core 的 `turn off plasma effects` 注释块（G:766起）和赤道面/fake j 块（G:776起）没有生效。本地 core 与 G 去注释 token 相同。DEBUG=1 确实打开 BHAC 的 NaN/异常诊断/退出；它**不是**控制 plasma effects 的开关。LIGHT_TRANSPORT、RAD_TRANS 虽定义为1，在所审查当前 src/BHAC model 没有相应使用点，不能靠改名字就宣称关闭了传输。
- `L/raptor_run_test`：Run 只选 `[100:101]` 一个snapshot、α=.3、100²、输出到 raptor_output_test；所有保存对象相对主版仅 pol_rte_integrator.o 不同。其 binary 的 pol_integration_step 确有额外分支：计算径向量×角度余弦后 `cvttsd2si` 做整数截断，再检查非零，并清零连续11个系数槽（发射/吸收/Faraday）。这支持“赤道附近筛选类实验”存在，**不能从缺失的历史 C 源完整恢复准确条件/厚度**。尤其整数截断特征与 `abs(double)` 误用相符；若原意是 .01/.1 的高度阈值，机器码不支持如此精细的浮点阈值，应核对 `fabs` 与实验定义。
- `L/raptor_run_test_eq_fakevalue`：当前 Run 选 `[151:152]`、α=.3、400²；头sigma=1000。保存的 **pol_rte_integrator.o 与主版逐字节相同**；两 binary 的 pol_integration_step、radiative_transfer_polarized、RK4/trapezoid 函数机器码逐字节相同，且全部13个对象仅 model.o 不同。因此当前保存文件不能被称为已启用“赤道面＋fake emissivity”版本；目录名、GRMHD_cut 残留、输出名都不是活跃配置证据。
- 该 eq 输出目录却有18个 snapshot，而当前脚本只选一个，已经足以说明目录积累了不同运行时期结果。不能把整个目录绑定到当前脚本/binary。
- 主/test/eq binary 的 Axion_Radial symbol 大小为 `23*2400*3*8`；旧 E 两 binary 为 `50*500*3*8`。旧 E/raptor core 必要 src 确认仍按38×490读取旧表，与新项目23×2400不相同；source/table不能跨版本随便替换。
- `E/raptor_run_gkx/Run.py` 主 snapshot list 为空，后面真正执行 `for name in [40]`、固定 data001700.dat、omega=.4；不是当前项目的141-snapshot批处理。旧 `E/raptor_run/Run.py` 活跃循环是 omega=.1，其尾部 `for name in []` 不执行。
- `E/test_raptor_difference` HEAD `79d7605`，大量 dirty。必要 core/pol src 已是 DM_mass/DM_cross_section/channel/pikachu 分支、电子分布表读取、额外 GRMHD_cut append；不是轴子主流程的另一个 clean 构建。其模板sigma/theta max均1000。若 RAPTOR误指这里，可能发生ABI/输入字段混用；不应作为生产 core。

### 6.3 与文章/图4的关系

`main.tex:355,367` 明确主数值包含真实 MAD plasma 的发射、吸收、Faraday；同时为分解 photon groups 做了“关掉赤道面外发射”、with/without absorption/Faraday，以及含/不含 gaussian smearing 的比较。`old_images/old.tex_copy:189` 四列设计是 toy → GR geodesic → 厚度 washout → 完整MAD平均（最后一列包含所有plasma effects）。因此这些实验开关有明确科学目的，不能一律删除或认定是生产污染。

但应区分：

1. fiducial：完整实际 GRMHD＋plasma effects；
2. w/o plasma：明确关闭吸收/Faraday，发射保持模型；
3. equatorial-only：限制发射位置；是否同时关吸收/Faraday须按实验定义单独记录；
4. toy/fake emissivity：替代 j 的受控验证；
5. LP/group/smearing 处理：分别说明在哪一层实施。

草稿图号尚不稳定：main.tex:408 的字面“Fig. 4”是 uncertainty 占位文字，开关意图实质在355/367和old189。不能只根据“图4”编号把某实验目录硬绑定到某幅生产图，也不能把 commented 实验当活跃。文章的3%/5%、snapshot数量和部分caption仍有待修改标记；本审查未验证这些数值。

## 七、批处理完整性、错误处理与并发

证据：Run.py:5–16,31–61,87–140；io.c:29–40,79,324；plotter-example.py:输出标题及保存路径。

- `os.system` 返回值从不检查；C exit非零、信号退出、plot失败都不会自动触发 Python except。
- 每个 snapshot 的 try 包住整个baseline+16相位，一次 Python异常会跳过后续相位；error_list只记异常字符串，无snapshot/alpha/phase上下文。
- done文件每个alpha循环后无条件创建；即使没有snapshot、部分失败也写done。空 log_run_error.txt不等于成功。
- 缓存只用目标 `img_data_i.h5` 是否存在，不核验尺寸、dataset、finite、完整写入、source/binary/table/config hash。改sigma/M_UNIT/表/步长后仍可能静默跳过旧文件。
- H5 使用 `H5F_ACC_TRUNC`，输出只含 I/Q/U/V/tau/tauF/alpha/beta 数据集，没有 git/binary/table/model/snapshot provenance attrs；Run还会删除 model.in。
- `move_file` 搬走 output 下所有普通文件、figures/* 全部图，并覆盖目标同名文件。失败后有残留时可能把另一相位/另一snapshot旧结果当本轮产物；未按严格本轮产物白名单归档。
- 同cwd并发写同一 model.in、output/img_data_0…16.h5、spectrum_<index>_<inclination>.dat、uniform_img_<frequency>_<index>.dat、figures/img_<index>.png，以及done/log。存在先检查后运行竞态、H5截断互踩、错误归档和model删除冲突。不能仅靠目标数据目录不同来并发。
- 多线程读表竞争另见5.3；Slurm/OpenMP线程数、栈、绑定、资源日志都应固定。
- 最后 `os.rmdir("output")` 无条件执行；若所有结果都缓存命中而output不存在，或有残留子目录，会在末尾异常。done可能已写，不能据此宣布完整成功。
- 壳命令不引用路径，snapshot路径含空格/特殊字符会分词甚至被shell解释；改用参数列表 subprocess.run(...,check=True,cwd=...)。
- plotter 固定 M=6.2e9、d=16800，没有从本轮config读；标题把输出index当 `t=ind*10*Tunit [days]`，既不是当前phase扫描的物理时刻，也未转换秒到天，是实际provenance/标注错误。
- 未来多频时 io外层每个频率重复写同一H5（writer内部又写所有频率）；不是本次单频主配置问题，但属于上游继承的低效/错误处理风险。

### 仅文件名层面的现有完整性

明确限定在列出的output根下两个目录层级，未读科学图像内容：

- 主 raptor_output 的 α=.1,.15,.2,.25,.3,.35,.4,.44 各有141个snapshot目录；每目录均有 img_data_0…16.h5，未见这些H5零字节。
- test α=.3：1个snapshot，17名齐全。
- eq_fakevalues α=.3：18个snapshot，17名齐全。

这只是**文件名/大小层面完整**。高α文件存在与当前低α表、当前Run list并不一致，说明历史运行步骤/表切换没有被源码和输出充分记录；不能据此证明高α结果错误，也不能证明正确，必须追回历史table/binary/config或重新在可追溯配置下生成。

## 八、重跑前的建议（本次未实施）

### A. Clean-source

- 在新的非挂载工作区建立专用 clean checkout，固定项目 `7d3ecac` 与core `b6df33c`，保存submodule状态；或另开显式审查后的upstream融合分支，不能静默换main。
- 保留有意义实验注释/历史补丁；把 fiducial、w/o plasma、equatorial-only、toy 分为独立命名配置和构建输出。不要用目录名猜开关。
- 明确sigma=1还是1000；主源码基线对应1，不能把dirty模板1000自动带入。将待确认物理修改做成可审阅patch，不与构建修复混在一起。
- 排除build/旧RAPTOR/ignored实验目录，明确RAPTOR绝对路径和cwd；对13个实际编译source及所有实际include头记录路径+hash。

### B. Compile manifest

- 全新build目录；修好目录order-only依赖、`.PHONY`、`-MMD -MP`、链接库顺序；每次变更配置形成独立binary，避免依赖mtime猜更新。
- 保存make -n完整命令、compiler/h5cc -show版本与wrapper展开、CFLAGS/LDFLAGS/LDLIBS、VPATH/include路径、相关环境白名单、预处理macro dump和include tree/依赖文件、git SHA/dirty patch、source/header/table SHA256、binary SHA256、ELF库依赖。
- 表先由串行初始化读取并严格校验，尺寸来自数据/元数据而非写死；高alpha选择应是显式run config/table manifest，按支持区间拒绝外推/越界。

### C. Run config 与可靠调度

- 使用显式JSON/YAML/逐行model配置，C端验证argc、字段名/顺序、scanf返回和取值范围；保存实际生成model.in，记录完整argv和cwd。
- 每个 `(snapshot, alpha, norm, phase, variant, resolution)` 独立临时工作目录/输出目录，只共享只读input与binary；参数列表调用、check=True、每任务stdout/stderr和退出码。
- H5先写临时文件，检查预期dataset/shape/finite/基本数值质量后原子归档；本轮白名单移动；done只在manifest所有任务validated后写。
- 缓存key含source/binary/table/config/snapshot身份hash，禁止只靠img_data_i存在来跳过。保留历史结果只读，新的重跑输出用全新run ID。
- 固定snapshot清单与排序/物理时间，明确OpenMP线程数/栈/绑定；纠正图标题把phase编号当时间的问题。
- 完成构建/解析/小型受控回归验证后，再由用户批准科学重跑；本次不执行。

### D. 需要用户确认的物理选项

1. fiducial是完整plasma effects，还是文章的w/o plasma/equatorial/toy对照；equatorial厚度与是否关absorption/Faraday要分别确认。
2. sigma cut=1或1000、theta_e cut、M_UNIT、R_HIGH/R_LOW与0.6Jy归一化要求。
3. GRMHD spin=.9375与高alpha .99表的生成spin是否一致；alpha与真实mode omega如何对应、径向/角向归一化，以及是否保留25Rg taper（当前活跃）/1022时间原点。
4. 扫描alpha集合、snapshot物理时间范围、16phase/norm=.1是否足够、100或400像素、stepsize/MAX_LEVEL、LP/group和smearing选项。
5. 是否以项目b6df33c复现，或接受上游更新后作为新模型版本；任何物理差异必须进入run ID和论文方法说明。

**最终判断：现有四处足以识别“干净项目源码＋真实主运行配置”的大部分链条，但不足以把历史高α/实验结果逐一绑定到当时binary/table/config。建议先解决构建依赖、旧model.in兼容、显式表选择及批处理provenance，再重跑；不要按当前目录名、done文件或dirty模板推定物理开关。**
