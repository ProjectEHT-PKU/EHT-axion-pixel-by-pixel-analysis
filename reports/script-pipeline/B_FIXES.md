# B类工程问题修复记录

## 范围

基于`SUBAGENT_REVIEW_VERIFIED.md`中B类问题修复。科学YAML值、频率态、耦合归一化、FOV/beam、mask默认值、null-only中心及likelihood模型不变。未做真实科学重算，未修改原EHT_ZC、输入链接或用户的eht-pipeline-runs。

## 已处理

| 问题组 | 实现及回归 |
|---|---|
| 输出安全 | common统一拒绝输出root/祖先/stage/文件链接，包括悬空链接；workflow在resolve前检查；RAPTOR/EHT/products/transport和辅助写入口同样保护。完整source_closure检查间接来源，拒绝将raw所在目录作为替换目标。 |
| 来源断边 | response manifest进入templates inputs，template manifest进入observations，EHT manifest进入background。验证递归且传播imported权限；inference及辅助读取检查必要规范边，包括inference最后一层。 |
| 声明输出覆盖 | require_declared_outputs保证实际消费的曲线、模板、观测、背景、diagnostic summaries和convergence run.json均由被验证stage声明；目录旁置文件不能冒充已验证产品。 |
| 严格YAML | config.py安全解析并完整schema检查，不做自动类型转换、补值或科学选择；错误bool、计数、非有限数值、unknown/missing字段等fail-fast。 |
| 时间分片 | 数值排序、步长10、包含端点帧数、连续性/重叠、非时间shape/dtype及重复源检查；res400双来源拒绝歧义。 |
| 收敛工具 | 递归验证+当前run绑定，先同base重算并匹配保存曲线再加密；默认覆盖三种mode；旧模型不冒充网格误差。 |
| 数值边界 | shared finite_evpa统一有效性；float64加权并拒绝溢出；全NaN padding与部分未知Stokes区分；机器精度退化圆均值无效，不引入任意物理浓度cut。 |
| legacy prior | legacy离散CDF仅接受其实际log-uniform测度，拒绝uniform虚标；corrected算法保持。 |
| 迁移来源 | 独立origin_inputs持续校验原manifest/产品，与新layout inputs分开；原科学identity和迁移identity都保存；旧摘要收据不能冒充来源闭环。 |
| 迁移原子性 | 整run preflight→同filesystem私有root→拓扑重绑定最终路径及hash→单次rename。失败清理本次构建及锁，原始source不改，可重试。 |
| 标签metadata | 保存已有snapshot标签及验证状态；缺失才用显式unverified默认标签，拒绝非法长度、重复/逆序。 |
| 辅助工具 | 完整run严格验证；独立canonical stage按其声明来源验证，不强制整个run；历史裸六文件标记科学lineage unavailable，不假称完整来源验证。 |
| 运行环境身份 | stage signature绑定Python/平台、科学依赖版本及pixi.lock；已有结果不自动重签。 |
| 历史null兼容 | 迁移新manifest只对options.observed严格False规范化，True拒绝，原receipt/source不改。 |

## 验证方法

- 所有新增输入均为pytest临时目录中的微型synthetic cubes、EHT候选图或旧layout fixtures。
- 完整真实脚本入口生成微型六曲线，固定旧主助手探针的12个数值作为回归oracle，要求绝对误差≤1e-12 dex。
- stale raw、imported祖先、旧缺边manifest、foreign inference、未声明旁置产品、输出链接、间接raw包含及运行环境变化均有拒绝测试。
- 迁移覆盖multi-alpha、完整旧run、原输入不变、临时路径重绑定、late failure、磁盘/rename等注入异常、重试、标签来源、live origin变化、真实inference opt-in reuse。
- 配置/shard专项318项、数值专项160项、迁移专项42项及独立复审；最终全套结果以本报告下节为准。
- 独立只读复审的边界发现均已处理，包括间接raw、inference最后一层、未声明产品、独立bundle与历史False。

## 最终验证结果

源码冻结后，后台tmux执行`pixi run test`：**742 passed，78.18s**。
随后`pixi run lint`及`git diff --check`通过。合法微型六曲线的12个cut
与修复前固定值在**1e-12 dex**容差内一致。测试/静态检查日志随本报告保存。
此次测试没有读取真实科学cube/pickle，也没有运行真实pipeline。

## 缓存与兼容注意

本次改变了源码签名、运行环境身份及必要来源边。旧缓存可能被明确拒绝；这不是自动判定旧数值错误，更不是授权删除或重建用户结果。请保留旧产物；重新计算或重新迁移要使用明确的目标，不能只修改manifest使它通过。

migration采用live origin契约，必须保留原归档。它验证旧来源**声明过的**图，不保证旧manifest从未记录的历史输入都已恢复。

输出检查防止正常误操作，不宣称抵御攻击者并发替换目录的完整TOCTOU沙箱。

## 未处理的物理事项

M/D及图像角标定、Stokes基底、n4/211、0.1与b/2π、超辐射解释、超参筛选概率含义、观测/模型相关性、幅相联合分布、Stokes曝光及1d背景匹配等，仍按复核报告等待生产记录或独立方法研究。本次不以通过工程测试宣称这些已解决。
