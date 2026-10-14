# GitHub RAPTOR仓库与实际生产工程对照

本轮只读检查GitHub、临时克隆及此前实际生产文件核查证据。未修改远端、未推送、未执行旧Run.py/RAPTOR，也未改生产数据。

## 1. 组织里的实际结构

GitHub公开API与git ls-remote确认：

| 仓库 | main提交 | 职责 |
|---|---|---|
| ProjectEHT-PKU/raptor | b6df33c69e650f3617845f4324ba642bdb878582 | 核心C、model模板、编译/绘图工具；已经含轴子修改，不是纯上游原版 |
| ProjectEHT-PKU/EHT-axion-pixel-by-pixel | 7d3ecac555e325e5ab6c160bfa00ff2af15e578e | 运行工程；raptor_run/参数、脚本、径向表，raptor/为submodule |

没有独立名为raptor-run的公开仓库；所说运行仓库是第二个。第二仓库.gitmodules指向第一仓库的SSH URL，gitlink固定b6df33c。两个远端当前都只公布main，没有高alpha分支/tag。

临时克隆：`/tmp/eht-github-raptor-GckNId/{raptor,run}`，没有在PHYS新建正式RAPTOR工作仓库。

## 2. 比较能力边界

此前挂载实际运行目录P=`vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel`的父/子仓库ref分别同为7d3ecac/b6df33c。

但当时实际子模块有未提交修改：model/bhac/definitions.h以及src/pol_rte_integrator.c，后者git diff --stat为1030insert/1016delete。不能因HEAD相同就说编译源码和remote完全一样，也不能仅凭行数说都是格式调整。

**当前vmdk为空，无法重新执行mounted-vs-remote全文件diff。**以下“实际生产”一栏基于先前实际读取、已保存审计/数值身份，不假装本轮重新验证所有文件。需恢复挂载才能完成源码增量和binary build链核验。

## 3. 配置和数据差异

| 项目 | GitHub运行仓库 | 此前读取的实际生产工程 | 裁决 |
|---|---|---|---|
| 核心ref | b6df33c | b6df33c，但dirty | 同基线，非保证同源码 |
| 运行ref | 7d3ecac | 7d3ecac，实际Run.py有后续修改 | 同基线，非完整生产快照 |
| M、D | Run.py:65–66，6.2e9/16800kpc | 相同 | 采用值一致 |
| camera、倾角 | 46.2261Rg、163° | 相同 | 一致 |
| 分辨率 | 100；high=400，默认highFalse | 相同默认，实际有400产物 | remote默认不重现所有产物 |
| 辐射参数 | M_UNIT3e25g，R_LOW=R_HIGH1，230GHz，step.005 | 当时读取同值 | 参数基准一致 |
| 参考强度/相位 | 0及.1、16相位p0..15 | 同 | 一致 |
| 扫描列表 | [.3,.25,.2,.15,.1] | 实际最后保存Run.py仅[.1]，产物有八alpha | remote只覆盖低五点；实际文件是末次配置非完整历史 |
| 径向表 | 仅NormR.csv：.09–.31，23×2400 | 另有NormR_0.99.csv：.30–.44 | 高三点来源未纳入remote |
| table identity | NormR SHA0a352558... | 本地EHT_ZC/NormR及新版Radial结果同SHA | **字节一致** |
| 高alpha流程 | 无.35/.4/.44、无高表选择 | 有对应cube；分析明确高三点用.99表 | 当前remote不能完整复现 |
| 输出转NPY | 未跟踪h5tonpy.py | 实际生产目录有此脚本 | 缺交付步骤 |
| transport variants | 核心文件有注释掉的实验段 | 实际有独立test/common/jonly/eqfake目录 | 当前remote不记录各variant完整配置 |

两张径向表对应及已生成八cube/四pickle身份见MOUNTED_CONSTANTS_COMPARISON、RADIAL_ARCHIVES_AUDIT，不在本轮重新哈希大数据。

## 4. remote内部需要整理的地方

### 运行脚本与示例配置冲突

remote `raptor_run/model.in`还是SgrA*式示例：M4.14e6、D8.127kpc、200²、camera40、inclination50、step.01，且没有AXION三个字段。

实际`Run.py`运行时重写model.in为M87参数；因此通过Run.py运行不是采用示例值。但`run.sh`直接调用该旧model.in和harm3d.txt，当前C解析器按固定顺序要求AXION字段且不检查读取返回值；示例MAX_LEVEL的1会错位读成axion_norm=1，之后EOF留下未正确设置的其它参数，不是“默认无轴子”。这条直接入口不应当作可用的M87生产入口。

Run.py生成配置没有换行，但%lf能在下一字段字母前结束，因此不能断言必然解析失败；它是脆弱格式。

### 重复model文件与构建选择

本轮逐字节比较remote：raptor_run/{definitions.h,model.c,model_definitions.h,model_functions.h,model_global_vars.h}与raptor/model/bhac/对应文件全部一致。

run makefile通过外部`RAPTOR`环境变量找src，而不是必然使用固定submodule。即使gitlink已固定，环境指向别处仍会编译其它源码。需要明确唯一构建源并记录hash。

当前对象规则只依赖.c而不追踪headers；改definitions.h后make all不保证重编译。因此历史dirty header也不证明保存binary实际使用了该配置。

### 高alpha查表不是自动支持

核心写死NormR.csv、omega_count23、r_count2400，没有0.99表选择，也没有完整alpha上界检查。不能只把Run.py列表扩到.44就认为能重现高alpha。

norm=0时也没有跳过查表，因此直接使用参数字典默认.4也有超出低alpha表的风险，主循环覆盖low5才避开。径向轴j有端点限制，外部半径是端区间外推，不应与alpha轴i无上界混称。

并行像素调用下静态flagg和首次表加载没有同步，存在初始化数据竞争；本轮没有据此断言历史结果损坏。

### 输出存在不是完整性证明

Run.py以img_data_i.h5是否存在跳过，os.system不检查返回码；move_file会替换同名目标，done可在错误后照写。整理时应将配置、表身份、源码、构建、命令返回码与完整输出验证一起记录，不能把done文件当科学完成标记。

### 轴子定义是共同遗留行为，不是本轮版本差异

remote同样用查表参数axion_omega直接乘time，含N/(2π)、1022、exp[-(r/25)^15]；α/ω分离尚未修。相机t_init0，TIME_INIT是输出编号。faraday/absorption实验开关目前是注释块，不能据注释推断实际variant。

### 另一个需要正式核验的核心索引问题

remote pol_rte_integrator.c:955–962从path_counter=steps−1开始，却读取X_u_next=lightpath[path_counter+1]。
gr_integrator.c:382–427只填0…steps−1，最后再将steps加1。因此首轮读取的是未初始化尾后记录，若填满max_steps还会越过数组范围。父已复读两处源码确认路径；未运行实际光追，未判断应改为哪一端。该问题应与后续alpha/omega修复分别测试，而非整理时直接凭猜测改索引。

## 5. 建议的整理边界

1. 保留两个仓库的职责：核心算法与运行配置/数据调度，不把分析pipeline或radial生成器重新塞进核心。
2. 核心以remote b6df33c为基线，恢复挂载后将真正使用的dirty差异作为明确提交；先不把unknown binary当可复现构建。
3. 运行工程固定构建到submodule，去掉错误的示例入口歧义，明确低/高alpha表与节点、自旋、分辨率、快照范围。
4. 补h5tonpy等实际生产交付步骤；独立variant用明确配置而非临时注释改源码。
5. 原径向表和科学产物保持不变；后续频率修正/方法对照与纯组织整理分开提交，不把二者混成一个差异。

本轮未开始正式移动/清理；主要阻塞是实际生产挂载目前不可读。GitHub仓库本身均可访问。
