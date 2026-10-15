# RAPTOR四处代码相对原上游审查：重跑前决策表

## 0. 结论

用户要求由subagent审查本地core/run、云端core/run，分清轴子实现、debug和故意关闭效应的实验，并结合文章判断重跑参数。本轮三个GPT-6.1 Sol分别完成偏振/轴子、流体/模型、构建/运行审查。模型agent因断网中断后恢复完成。主助手核查报告关键源码、Git差分和一手论文，独立复核保存ELF的sigma比较常量，并重新编译运行纯合成solver探针。

**没有修改生产源码、挂载文件、GitHub远端或科学数据，没有执行保存RAPTOR binary，没有重跑真实图像。**完整科学程序尚未在当前环境编译；合成探针不是科学计算验收。

主要判断：

1. 主代码不是偷偷关闭plasma的jonly：真实thermal发射、吸收、Faraday rotation/conversion及GR输运保留。注释中的关效应/假发射有文章诊断意图，不应按垃圾删除或自动激活。
2. 当前主run的有效sigma cut可确认是1，不是core模板后来改的1000。1适合有条件fiducial；10优先systematic；1000为有标签的压力测试/另立物理模型，不是天然更准确。
3. 自定义轴子链存在活跃端点、fluid gate、表初始化/边界和物理约定问题，必须先验证修复，不能直接复用当前wrapper启动所有cube。
4. 电子温度闭合、旧thermal Faraday fit是继承/版本选择，不是本地debug。改变它们应与算法修复分开、重新校准flux。
5. 保存binary/对象已帮助识别实验构建，但不能把今天这些文件身份追溯成每个历史cube的完整来源。

## 1. 对照基线和四处路径

```text
L = /home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel
LC = L/raptor                       本地挂载核心
LR = L/raptor_run                   本地挂载主运行工程
CC = ProjectEHT-PKU/raptor           GitHub核心
CR = ProjectEHT-PKU/EHT-axion-pixel-by-pixel/raptor_run
```

| Git身份 | commit |
|---|---|
| 原版上游main：jordydavelaar/raptor | b91e937b7b746eb24992b6c1dd8831c5dea0cef1 |
| 我们fork main | b6df33c69e650f3617845f4324ba642bdb878582 |
| 共同祖先 | cdc2547b868ed8d6ae0c9a5428bc1ea86fce901a，2022-11-15 |
| 我们run main | 7d3ecac555e325e5ab6c160bfa00ff2af15e578e |

共同祖先→fork为我们自己的改动（5个提交）；共同祖先→上游main为之后13个上游提交。不能把二者混为“我们相对今日上游故意删掉的物理”。

### 四处实际差异

| 项目 | LC | LR | CC | CR |
|---|---|---|---|---|
| SIGMA_CUT源码 | **1000，dirty模板** | **1** | **1** | **1** |
| model.c、model_definitions、grid等 | 四处各自对应文件字节一致 | 同 | 同 | 同 |
| pol_rte active token | 与CC一致；多排版/注释 | 编译时取core/src，无独立源码 | fork active基线 | submodule固定fork |
| Run.py扫描 | 不适用 | 当前最后配置[.1] | 不适用 | [.3,.25,.2,.15,.1] |
| model.in | core旧示例 | 运行后已删除 | core旧示例 | 旧SgrA示例，无AXION字段 |
| 低alpha NormR | 不提供run表 | .09–.31/23×2400 | 不提供run表 | 同LR字节 |
| 高alpha NormR_0.99 | 不提供run表 | .30–.44/15×2400，额外文件 | 不提供run表 | 未跟踪 |
| high-res/导出 | 不适用 | 有400产物及h5tonpy额外脚本 | 不适用 | high默认False，缺导出步骤 |

四处inventory确认core59个tracked文件仅两个本地字节不同；run14个源/参数tracked文件仅Run.py变化和model.in消失，另有23个小型未跟踪run文件。两份CSV均与radial-function的reference字节匹配。

## 2. 相对原版上游增加/改变了什么

下表用fork相对共同祖先，不把上游后续修正归为本项目修改。

| 类别 | 我们的修改 | active？ | 对重跑的判断 |
|---|---|---|---|
| 轴子输入 | 新增norm、名为omega的质量查表参数、phase；经core/functions传入偏振积分 | 是 | 保留能力；alpha和omega应分开、参数严格校验 |
| 局部轴子旋转 | RK4/trapezoid完成plasma后，再用端点角场差旋转Q/U | 是 | 科学核心，非debug；须修端点/符号合同并验证分裂误差 |
| 径向表 | 静态23×2400、读工作目录NormR、两维线性插值 | 是 | 保留查表；加载/并发/高alpha维数/域检查需要修 |
| 场表达式 | N/(2π)sinθ R cos[-φ+(1022+t)alpha+phase]及25Rg指数窗 | 是 | 定义需显式；timealpha、BL/KS、人工window分别决策 |
| 调用gate | 上游valid-fluid AND半径域，改为OR | 是 | 有让真空轴子传播的合理意图，但当前强行计算plasma不安全；应分离 |
| 光线next端点 | 新增lightpath[p+1]供场差分 | 是 | 首轮未初始化/可能越界且与正向段错位；优先修 |
| block/pixel | 新增位置编号传递 | 是，但主物理不用 | 用于诊断；不是额外物理 |
| BHAC适配 | SFC0→1、根网格128×64×64→96×48×48、hslope.25→0、域边界变化 | 是 | 与当前流体格式相关，不能盲目恢复上游示例 |
| 模型物理式 | model.c新增调试注释；温度、units、sigma式未被fork改写 | 注释 | inherited模型选择，不归为我们debug |
| metric.c | 排版空白 | 否 | 未发现Kerr metric自定义公式改动 |
| 诊断输出 | camera printf；io路径buffer64→128 | 是但不改物理 | 可保留日志/工程修正 |
| 绘图工具 | M87尺度、tau图方向/色标/log、双频RM读图及除Δλ² | 仅绘图 | 不改变生产IQUV；图像手性/标签仍需一致 |
| USERSPEC/Camera字段 | USERSPEC0→1、移除未用pdf/avg | 是，谱函数体被注释 | 不等于开启非热EMISUSER；不应误判 |

### 关闭效应的注释应怎样理解

| 注释块 | 当前主模型状态 | 有意义用途 | 主full是否启用 |
|---|---|---|---|
| rQ/rU/rV=0 | 未启用 | 去掉Faraday rotation/conversion | 否 |
| aI/aQ/aU/aV=0 | 未启用 | 去掉吸收，配合前项形成jonly | 否 |
| 赤道高度外j/a/r=0 | 未启用 | 厚度/多穿越诊断 | 否，独立eq模式 |
| 常数或高斯环j替换真实源 | 未启用 | fake源/总导数/几何校验 | 否 |
| X0<-10023关闭发射 | 未启用 | 时间gate/特定光线片段实验 | 否；不是通用direct-only |
| 均匀toy轴子场、替代径向窗 | 未启用 | 解析回归 | 否 |
| rV+=轴子项 | 未启用 | 早期另一种实现尝试 | split方案下不启用，防双计 |
| 小场阈值或pixel编码覆盖Stokes | 未启用 | 历史debug | full关闭 |
| GRMHD_cut/phi_cut并行append | 未启用 | 轨迹/场诊断 | 可具名开关，有序写而非无锁混写 |

文章main.tex:355,367和旧稿的toy→GR→厚度→full顺序支持这些诊断。作者宏中的3%/5%/30%及待改数字不能代替本轮科学验证。

## 3. SIGMA_CUT：哪一个适合主重跑

### 当前实际状态已经确认

主保存binary SHA `7137ba3b...`中磁化比值比较的ELF double为**1.0**；eq_fakevalue SHA `bd738b9e...`为**1000.0**。主助手独立从divsd后comisd地址解码复核，没有执行binary。

主源码比较sigma≈b²/rho，相当于CGS的B²/(4πρc²)。cut条件设置n_e=0，正常有限状态下关掉该cell全部plasma j/a/r，不是Gaussian σ，也不是混合模型sigma_min。现caller OR仍继续轴子；这有物理意图但需安全真空分支。

| 选项 | 建议定位 | 理由 |
|---|---|---|
| **1** | **当前thermal/full MAD的有条件fiducial** | 延续主run；EHT V/VIII及项目引用Mizuno2021排除sigma>1 floors影响区 |
| **10** | **首要systematic** | EHT VIII有直接全模型补充库先例，但其评分稳定不证明我们像素级轴子量稳定 |
| **1000** | **有标签的压力测试/特殊模型** | 放行大量高sigma、可能floor控制的区域；不等于更真实，也不必是“错误值” |

阈值不是物理发光开关。若目标包含真实漏斗电子对或非热辐射，单纯把cut提高不等于实现了相应电子物理。

改变cut应分两类实验：固定MUNIT诊断新增区域的贡献；各模型单独校准同一时间窗平均flux后比较科学结果。不能把输出乘个系数替代RT，也不应逐snapshot强制0.6Jy而声称保留原flux变率。

## 4. 主重跑参数建议：事实、基准与待决定分开

| 项目 | 当前主run | 重跑建议 | 证据/缺口 |
|---|---|---|---|
| M/D | 6.2e9/16.8Mpc | 保持，统一下游单位 | 已证实本批生产尺度 |
| 流体自旋、NSPIN | .9375、3 | 当前数据保持3 | footer index3=.9375，index5=20，不能盲改5 |
| SFC/grid/hslope | 1/96×48×48/0 | 先保留数据适配 | 应补生产参数证明；不能默认上游示例正确 |
| POL/DF | 1/TH | 保留 | 正式thermal偏振链 |
| j、a、Faraday | 全开；B-aligned frame jU/aU/rU=0 | full全部保留 | local U=0不是关闭天空U |
| SIGMA_CUT | 1 | fiducial1、systematic10 | 1000另立压力测试 |
| R_LOW/R_HIGH | 1/1 | 历史benchmark先保留，不称最佳M87模型 | 等于constant R=1；可再加冷电子分支 |
| 电子温度闭合 | (u/rho)(mp/me)/[3(1+R)] | 明确保留历史式或另立EHT式对照 | EHT式为2/[3(2+R)]，R1时当前温度为其.75倍；继承而非fork新增 |
| M_UNIT | 3e25g | 只作校准初值 | 代码没求解0.6Jy；文章说法仍需平均flux实证 |
| thermal rhoV fit | 旧Dexter型 | 历史基准显式保留；新fit另做对照 | 不能静默同步上游造成模型变化 |
| 频率 | 单色230GHz | 初始保持；精确观测频段差异另列 | FITS约226.19GHz；不宣称宽带 |
| cameraFOV/倾角 | 46.2261Rg/163° | 保持真实输入几何 | output角尺度和blur属后处理另定 |
| rcam/RT_OUTER | 1e4Rg/50Rg | 历史基准保留，后测外边界 | RT50不等于相机距离或全云积分范围 |
| theta_e窗 | .001–100 | 初始保留、记录删点/贡献 | 不是clamp，有完整plasma cut后果 |
| horizon margin | .01相对r+ | 保留并检查收敛 | 实际r_stop=1.01r+ |
| STEPSIZE | .005 | 必须小步长对照再决定 | 含轴子源项split不因RK4名字就四阶 |
| 径向自旋节点 | 低五.9375，高三.99 | 显式记录两表和节点 | 非同spin完全模型；已由用户及分析调用确认 |
| alpha/omega | 同一变量 | 必须分开 | 质量查表与时间频率明确不同 |
| 场坐标 | 直接KS t/phi | 与文章BL约定统一或明确另立toy | 径向依赖phase不能由常数1022消除 |
| norm/phase | .1、16等分 | 保留小信号测试、明确N=bc与约定 | b只计入一次；norm0安全旁路 |
| t_ref1022 | 固定 | 配置化并统一full/screen | 整体phase可被nuisance吸收，非自动bug |
| 云窗 | exp[-(r/25)^15]开启 | **科学待决定**：legacy保留，测试无窗/扩大窗与真空补齐 | 不能因像debug就删除，也不能未经验证当物理fiducial |

上述建议不是一次性改动清单。尤其温度公式、Faraday fit、云窗不是“顺手修bug”的内容，必须分开模型对照。

## 5. 重跑前必须处理/验证的活跃问题

| 优先 | 问题 | 确认证据 | 正确处理边界 |
|---|---|---|---|
| P0 | X_u_next=p+1首轮读lightpath[steps] | 存储仅0..steps−1；达到maxsteps还可越界 | forward端点、delta与旋转符号成套，不机械只改索引 |
| P0 | valid-fluid AND→OR半径gate | 失败/陈旧modvar仍可进pitch/tetrad/coeff | valid fluid算plasma；真空/sigma-cut/pitchaligned仍独立算轴子 |
| P0 | 静态NormR首次加载race/未检查/越界 | OpenMP首次读写flag无同步，synthetic .4报index23 | 并行区前一次load、完整验证、动态尺寸/显式表选择、norm0旁路 |
| P0 | alpha用于cos时间项、BL/KS及full/screen不统一 | producer质量表、原场表达式及文章BL分离式 | 先定义唯一模型，验证坐标变换和场包络近似 |
| P0 | 旧model.in字段不匹配 | parser probe读成norm1、omega0等 | 标准逐行输入、字段/返回值校验，禁旧入口 |
| P0 | build/header/table/mode来源可漂移 | env VPATH、header无依赖、mtime对象复用 | 新clean build+manifest；每任务隔离，记录hash/macro/library |
| P0 | 批处理只看文件存在、忽略返回码 | Run.py os.system/done/搬移策略 | 不复用旧输出、新runID、退出码/完整文件验证及事务发布 |
| P1 | 轴子plasma Lie分裂整体一阶 | synthetic constant-source步数加倍误差减半，父复跑确认 | 选正确端点后的受控split/联合解，做收敛，旧法显式legacy |
| P1 | 继承stiffness/clip等保护解释与实际不同 | 零系数仍STIFF，complex sqrt/单位问题 | 常系数解析测试、统计solver/clip触发；不要归为项目故意关物理 |
| P1 | flux与温度/R配方未验收 | 无0.6Jy求解，闭合不等EHT公式 | 每模型同时间窗一个MUNIT校准，保留flux变率 |
| P1 | 独立test来源混杂 | 对象机器码与目录命名不一致 | 诊断模式显式版本化，不从目录名恢复历史配置 |

没有证明：旧cube受影响多少、哪个修复使上限变强/变弱、25Rg窗唯一作者原意、所有现有径向表错误、3%/5%声明已被本轮否定。当前.1径向边界修复重算差1e-14，八节点归一化通过，不能将未使用slice问题重列为当前主数据缺陷。

## 6. 诊断变体应该保留，但定义和来源必须清楚

| 模式 | 真实发射 | 吸收/Faraday | 薄层限制 | 轴子/GR | 用途 |
|---|---|---|---|---|---|
| full/common | 同一GRMHD thermal | 全开 | 无 | 全ray/全GR | 主基准 |
| jonly | 同一GRMHD thermal | 全关 | 无 | 全ray/全GR | 隔离传播plasma，不等于screen |
| eq_jonly | 同一GRMHD thermal | 全关 | 明确赤道发射mask | 全ray/全GR | 厚度/穿越组诊断；轴子不随发射mask消失 |
| fake | 显式常数/高斯环toy | 明确选择，通常全关 | 独立选择 | 明确选择 | 总导数/几何回归，不是full |
| screen | baseline图后处理 | 不重新联合传播 | 单发射映射 | 映射已含GR | 近似对照 |

### 保存实验对象实际提供的证据

- raptor_run_test的pol.o有active z=r cosθ转整数再清零11个系数。对应int abs语义，约|z|<1Rg，不是注释中.01/.1Rg薄层。
- eq_fakevalue的pol.o与主run字节相同，13个对象仅model.o不同；其binary阈值1000。不能凭名字认证它当前使用fake发射或全域jonly。
- 二者无独立pol.c源，外部RAPTOR路径不固定；今天objects不证明每份历史common/jonly/eq_jonly/fake cube的来源。
- 另test_raptor_difference是DM/pikachu分支，不是可直接替代的轴子clean core。

因此文章相关百分比应在同snapshot、same cutoff、same normalization/table/frequency/basis等统一条件下重新计算。不要把不同cut的图差别全归为厚度或plasma传播。

## 7. 上游后续修复：哪些会影响当前thermal主链

| 上游变化 | 当前影响 | 建议 |
|---|---|---|
| 添加EMISUSER偏振调用分支 | 当前0无影响；旧fork只single | 未来启用必须真的实现调用，非只改宏 |
| thermal jQ/rhoQ中间翻号 | 后来撤回，最终与fork一致 | 不摘中间commit“修正”符号 |
| KAPPA/POWER符号、get_w等 | 当前TH不活跃 | 非热研究另同步整套约定 |
| thermal rhoV新Shcherbakov fit、DEXTER宏 | **当前TH会受影响** | historical oldfit与newfit显式对照，禁静默拷贝 |
| model spin越界报错 | 健壮性改善 | 加入，但不能代替数据格式正确性 |
| setup --nfreq/README/license等 | 工程/交付改善 | 保留上游身份许可，独立于科学模型更新 |

## 8. 交付及后续顺序

完整子报告：

- [偏振/轴子核心](raptor-upstream-evidence/AXION_POL_REPORT.md)
- [流体/模型/参数及一手论文](raptor-upstream-evidence/FLUID_MODEL_REPORT.md)
- [构建/运行/debug provenance](raptor-upstream-evidence/BUILD_RUN_REPORT.md)

审查证据包括fork-vs-ancestor/upstream-later diffs、四处文件SHA、编译/解析/表边界/solver合成探针及binary静态记录。大型GRMHD、保存binary、文献PDF不作为代码提交复制。

建议下一步：**先建立独立、保留GitHub历史的core/run工作区，归档四处差异和真正有效的主模型参数；再做工程与轴子数值修复；物理模型变化分别对照；小规模验收后才全量重跑。**这份报告不构成自动执行或全量重跑授权。
