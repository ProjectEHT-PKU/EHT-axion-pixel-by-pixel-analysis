# 工程修复后的输入、输出与缓存契约

这次修复不更改科学参数、null中心、观测量定义或likelihood公式。它使不安全、失效或歧义输入明确失败，并补齐验证范围。

## 本机env文件与参数来源

- 主分析阶段及`run_all.py`自动读取仓库`ROOT/.env`；生成入口依次读取`ROOT/.env`和`ROOT/.env.raptor`，后者覆盖前者的文件值。分析不读取生成端文件，也不从cwd或祖先目录暗中查找文件。
- 优先级为**CLI > shell已导出的环境变量 > 文件值**；自动读取不会覆盖外部已导出的值，也不要求用户先手动导出或shell加载。修改文件后下一次调用重新读取；外部旧值可能遮住修改。
- 使用`python-dotenv`、`interpolate=False`解析变量值，不执行shell、不做变量展开或命令替换。不要将env文件当作脚本；路径应填写实际值。
- 文件提供的相对路径锚定仓库ROOT（包括自定义env文件中的路径）；CLI和外部环境中的相对路径仍相对cwd。
- `--env-file PATH`替代默认文件选择，生成端也不再读取默认两文件；选择器PATH相对cwd。`--no-env-file`禁用文件loader但不清除外部环境，独立运行应显式传入所需CLI配置和输入/输出路径。
- `.env.example`仅含分析/共享设置；可选生成端首次复制`.env.raptor.example`为本机`.env.raptor`并审阅即可，不覆盖已有文件。本机env文件不进入Git，物理方法和数值仍由科学YAML定义。
- 布尔guard保持不变：接受大小写不敏感的true/false、1/0、yes/no、on/off，拒绝非法值。pickle信任、重建和导入复用不会因自动加载而默认开启；显式`--no-trust-pickle`、`--no-rebuild`、`--no-reuse-imported`仍可关闭对应选项。
- 自动加载仅改变参数来源，不绕过输出保护、来源校验或RAPTOR源码锁：核心checkout仍须干净且commit匹配。它不启动生成，也不构成真实科学计算授权。

## 输出路径

- 写入不允许经过输出根、祖先、stage目录或文件软链接，包括悬空链接；只读科学输入仍可使用软链接。
- 目录替换还检查完整来源图：即使目标内的raw只被外部manifest间接引用，也拒绝删除/替换。
- 即使不使用`--rebuild`，也不能通过链接创建外部目标。
- 如果使用旧兼容目录别名作为输出，请改用其真实目录路径；不要删除链接目标。
- 低层RAPTOR、EHT、观测量、transport及绘图/诊断入口同样拒绝危险输出链接。
- 独立绘图、诊断和收敛输出应使用新的文件/目录，拒绝覆盖已有产品。主stage重建仍由事务目录替换处理。
- 这些检查防止常见误操作；不是在攻击者可并发替换目录的环境下提供完整OS沙箱。

## 科学 YAML

`configs/research.yaml`的数值未改变。schema现在拒绝未知或缺失字段、非有限数值、错误类型和非法范围。

```yaml
historical_mask: false       # 正确：布尔值
# historical_mask: "false"   # 错误：字符串，不再被静默当作真
```

计数必须为整数，不能用bool或浮点数冒充。alpha必须正、唯一、严格递增且符合两位小数的规范身份。配置检查不自动补值、排序、改正或增加物理选择。

## 时间分片

自动选择res400输入时：

- `stokes-tSTART-tEND.npy`按整数START排序；END为包含端点的标签。
- 当前标签约定步长10，帧数须与区间一致，非时间shape/dtype一致。
- 拒绝时间缺口、重叠、重复文件来源，以及完整文件/shards或canonical/legacy双来源歧义。
- 不自动插值、填补缺帧或选择某一个冲突来源。
- 检查通过header/只读mmap完成，不为预检复制科学cube。

这些检查确认命名/结构一致性，不独立证明标签与真实物理时间的映射。

## 来源与缓存

必要的依赖manifest现纳入来源图：response→templates、templates→observations、EHT→background。

- 推断、绘图和收敛不仅递归校验已有声明，还核对必需的规范依赖边，拒绝缺边的旧来源图。
- imported身份沿manifest传播；使用导入产品仍须明确`--reuse-imported`。
- 缓存签名新增Python/平台、NumPy/SciPy等版本及`pixi.lock`身份。软件环境改变，不再静默视为同一重算缓存。
- 旧缓存可能因此被拒绝。这不证明旧数值错误，也不会自动删除、补签或重算它。
- 建议保留旧输出，使用新的输出根重新计算需要的stage；若确实希望重建现有根，先检查目录内容与当前任务，再显式使用`--rebuild`。
- 仅更新旧manifest来使验证通过是不允许的做法。

## 迁移

迁移仍不做科学重算，也不修改原数据。当前契约采用**持续校验原来源**：

- `inputs`表示新规范布局依赖；独立的`origin_inputs`表示实际读取的原manifest和产品。
- 原计算implementation/signature与迁移implementation分别保存。
- 旧null metadata中的`options.observed: false`仅在新manifest配置中规范化为当前无选择器形式；true始终拒绝，原来源记录和文件不改。
- 旧静态收据（有manifest摘要却没有可校验origin inputs）不能被当作已验证来源；应保留旧产物并迁移到新目的地。
- 完整preflight后，在目的地同filesystem的临时root中构建，成功才一次rename发布。失败清理本次临时构建，不删除未知用户目录。
- 发布前按拓扑顺序重绑定最终路径和hash，不能留下指向临时root的依赖。
- 由于大数组仍链接原文件，原归档必须保留；可复制小产品并不意味着可以删除原科学来源。
- 有原snapshot labels时保留并说明验证状态；缺失时默认1600+10*i，明确标记unverified，不冒充已知时间。

## 数值边界与兼容

- Q/U非有限或线偏振为零时EVPA无定义，返回NaN，而不是利用atan2偶然得到有限角。
- 积分前转float64并检查加权溢出；仅四通道全NaN空间padding当作零flux。部分未知I/Q/U明确拒绝，不自动补零。
- EHT resultant在机器精度上退化时设统计量无效；没有新增人为的物理浓度cut。
- legacy等log网格CDF仅支持其实际log-uniform测度；显式uniform请求报错。corrected主链仍支持两种prior。

## 独立诊断工具

- canonical曲线要求有效stage来源，实际读取文件必须在manifest声明的outputs中；额外放入目录的文件不能冒充已验证产品。
- stage身份为`inference:six-curves`时要求完整run来源闭环；独立的规范曲线bundle可按其自身stage验证，不强制提供RAPTOR/EHT整run。
- 历史六文件直接绘图仅验证数值schema，输出明确标记科学lineage unavailable。
- 比较允许有意的算法差异，但必须记录各自身份，不能混合一组曲线内的prior/variant。
- 收敛工具先复算base并匹配保存曲线，再改变网格；默认覆盖std/cor/cor_grmhd。`--reuse-imported`只授权来源，不跳过baseline一致性。
- 此次没有加入物理收敛接受阈值，也没有用数值测试替代生产配置、基底及统计coverage验证。
