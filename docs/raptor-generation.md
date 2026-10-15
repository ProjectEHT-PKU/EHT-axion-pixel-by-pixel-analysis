# RAPTOR生成端：官方构建、Python调度、现有分析

## 分工与边界

核心源码仍在独立 `../raptor` 仓库。pipeline新增的是生成端编排，
**复用核心的官方setup.sh和makefile**，不维护第二份头文件模板。
Python入口调用官方脚本、生成model.in、调度/验证输出和恢复任务。
不会调用旧run.sh、Run.py、plotter或Notebook。

当前源码锁定在 `raptor-source.lock.json`。要求本地checkout干净且commit一致；
修复分支尚未推送，不自动clone/pull，也不假设GitHub已经包含这些提交。

默认 `pixi run pipeline`仍只分析已有数据，绝不自动build或光追。
本次移植仅验收了工程和合成数据，**没有运行真实快照**。

## 环境与路径

从仓库根目录执行。所有依赖/任务在pixi.toml，生成环境首版仅支持Linux x86_64：

```bash
pixi install --locked -e raptor-build
```

`-e raptor-build`选择包含C编译器、GSL、HDF5、make和rsync的生成环境。
默认分析环境没有新增编译工具依赖。

不要覆盖已有.env。可把以下项加入自己的.env，审阅后显式导出：

```bash
set -a
source .env
set +a
```

| 变量 | 当前机器可用示例 |
|---|---|
| EHT_PRODUCTION_CONFIG | configs/production.yaml |
| EHT_RAPTOR_SOURCE_ROOT | ../raptor |
| EHT_GRMHD_ROOT | /home/canoziia/GRMHD_DATA/a094-high |
| EHT_RADIAL_LOW | ../radial-function/reference/result_a=0.9375_l1m1.txt |
| EHT_RADIAL_HIGH | ../radial-function/reference/result_a=0.99_l1m1.txt |
| EHT_RUN_ID | benchmark-211-m65 |
| EHT_OUTPUT_ROOT | ../eht-pipeline-runs |
| EHT_RAPTOR_THREADS | 2 |

如果现有.env还没有生成端选项，不要覆盖它；本机可先在终端显式设置：

```bash
export EHT_PRODUCTION_CONFIG=configs/production.yaml
export EHT_RAPTOR_SOURCE_ROOT=../raptor
export EHT_GRMHD_ROOT=/home/canoziia/GRMHD_DATA/a094-high
export EHT_RADIAL_LOW=../radial-function/reference/result_a=0.9375_l1m1.txt
export EHT_RADIAL_HIGH=../radial-function/reference/result_a=0.99_l1m1.txt
export EHT_RUN_ID=benchmark-211-m65
export EHT_RAPTOR_THREADS=2
```

也可直接指定 `--run-root`/EHT_PRODUCTION_ROOT。CLI优先于环境；相对路径以当前cwd为准。
Python不加载或执行.env。生成阶段不用EHT的pickle输入，不需要trust-pickle。

## 目录

```text
../eht-pipeline-runs/production/<run-id>/
  config.yaml                    不覆盖的科学配置快照
  requests/                      每次调用的追加记录
  build/                         官方setup产物、binary、源码快照和receipt
  jobs/<task-id>/                 成功任务：model/grid、HDF/spectrum、log、manifest
  attempts/<task-id>/<attempt>/   失败/超时/取消记录与部分输出
  data/raptor/res-100/alpha-0.xx/stokes.npy
  data/raptor/res-400/alpha-0.30/tSTART-tEND/stokes-tSTART-tEND.npy
  analysis/                      另行调用现有分析产生
```

`alpha-0.30`背景目录是旧分析接口的兼容定位，不代表无轴子背景属于某个云质量。
背景metadata明确alpha=None/norm=0；每快照每分辨率baseline只计算一次。
run根只是容器，不能当成可整目录rebuild的stage。

## 先规划（不执行C）

导出上述变量后：

```bash
pixi run -e raptor-build generate-raptor --plan
```

它读取快照footer和两个小CSV，检查真实自旋/时间、alpha覆盖、参数与路径，
打印任务数和资源下界；**不完整读取97GB快照、不创建run、不build、不运行模拟**。
实际执行时再完整hash所消费输入，每个调用内共享同一快照的任务只hash一次。

生产profile全范围：141快照，18,048响应任务，282baseline（100和400各141），
共18,330任务。高分辨率仅background state0，不生成16相位。
所有任务串行调度，每个RAPTOR使用指定OpenMP线程；未做pilot前不承诺耗时/RAM。
计划中的RAM/空间是下界，不是实测；执行前另留磁盘余量。

## 构建（不运行快照）

```bash
pixi run -e raptor-build build-raptor
```

Python检查锁定核心，保存Git tracked源快照，在私有cwd调用其setup.sh。
官方setup负责复制BHAC模板、应用明确选项和make；pipeline核对最终宏。
构建一次用于全部相位，已完成构建必须与源码/配置/环境完全匹配才可复用。
不依赖旧.o或源目录mtime。

核心setup有独立工程小适配：严格失败退出、grid选项大小写、显式sigma/温度窗/UNIF、
可靠模板复制；make分离CC和include flags，并适配当前HDF wrapper。
这些没有更新原版C求解器或辐射系数。

## 执行与恢复（会进行真实科学计算，请明确启动）

```bash
# 获得pilot授权后，最多执行任务清单的前N项：
pixi run -e raptor-build generate-raptor --max-jobs 1 --timeout-seconds 1800

# 明确决定全量生成后：
pixi run -e raptor-build generate-raptor
```

前N项按规划顺序选，不表示完成某个alpha的16相位；先列plan核对。
默认顺序每快照先100/400baseline，再各alpha16相位。
恢复使用同一命令、同一run和配置；已成功任务校验后跳过，失败任务重新尝试。
旧任务存在但配置/binary/source变了则拒绝，不自动覆盖或重签。
任务有独占锁；不要在未检查进程前手删锁。

每任务独立cwd、保存19字段model.in；第三argv是0..16输出slot，
**不是GRMHD物理时间**。真实时间来自footer，另行记录。
失败/非零返回/超时/取消不会写完成manifest，保留attempt日志和部分HDF；
进程组超时/取消会终止子进程。运行时不搬动、不删除源快照和表。

## 转换和分析

完整任务执行后自动交付；也可只转换已完成任务、不启动binary：

```bash
pixi run -e raptor-build generate-raptor --convert-only
```

转换入口同样检查build/plan/job依赖，不能给错误模型的旧任务签发新metadata。
验证IQUV/tau/tauF/alpha/beta、frequency、finite、逐元素坐标顺序、spectrum flux。
复用现有convert_hdf5逐帧写mmap，不擅自转置/翻转/U翻号，不按HDF key排序猜通道。
输出强度保持Jy/sr，不乘像素面积；用面积积分对照spectrum Jy只是序列化验收，
不等于完成0.6Jy校准。

交付res100为(T,17,100,100,4)，state0 baseline、state1..16 phase0..15；
res400为(T,1,400,400,4)时间分片，保留单state轴。缺任务不零填充。
原HDF/spectrum永久保留为可验证祖先，本流程不自动删除大数据。

随后独立分析：

```bash
pixi run pipeline --config configs/production.yaml \
  --input-root data/inputs \
  --raptor-input-root ../eht-pipeline-runs/production/benchmark-211-m65/data \
  --output-root ../eht-pipeline-runs/production/benchmark-211-m65/analysis \
  --trust-pickle
```

EHT输入仍是已有data/inputs；新RAPTOR root显式覆盖只影响模拟和background。
指定的新root缺数据不会回退旧输入。生产配置不接受无生成身份的旧NPY；
生成数据也不能用analysis-only profile绕过一致性检查。
labels传给fit，不使用默认标签；后续inference/plot lineage保持到binary/fluid/table。

## 科学配置与缓存注意

production.yaml默认质量为6.5e9 Msun、距离16.8Mpc，与2025 EHT 2018 Paper II主模型库一致。
相机FOV仍明确固定46.2261 Rg（不是固定角FOV），因此质量变化会改变角宽。
低五云表spin.9375、高三.99，背景仍.9375；其它后处理几何沿用旧值，
不表示beam/EHT角尺度已科学校准。

alpha到omega只调用 `eht_axion.model.omega(alpha, 2)`；C的alpha与omega输入独立。
后续更换频率来源须统一生成/分析并用新run，不能混用旧响应。
reference norm固定.1、16phase，MUNIT3e25/R1/1/cut1为明确benchmark，
不宣称已完成平均0.6Jy校准。Dolan扒图频率不参与本次生成。

新增pixi环境和shared代码改变lock/source身份，旧cache可能stale。
**不重签、删除或自动重算旧结果**。默认research.yaml和已有输出不改。

## 测试

```bash
pixi run test
pixi run lint
```

测试使用合成footer/HDF、fake executable、官方setup假工具链，涵盖错误输入、
失败/超时/恢复、锁、源码/宏/输入变化、坐标/单位/完整性、baseline共享、
生成身份防降级、convert-only绑定和新数据接入真实后处理链。
完整实际C构建单独验收，未用真实fluid执行生成。
