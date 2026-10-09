# 两份 GPT-6.1 Sol 审查的逐项复核

## 范围、方法及证据边界

审查基线：`refactor/script-pipeline`，`b2cd143`。本报告是主助手自行读源码、公式与一手论文后，对两份subagent报告的复核，而不是转录其结论。

本轮未修改任何执行代码、科学配置、真实输入、原始EHT_ZC或用户运行输出；只新增本报告及小型证据文件。没有启动真实科学重算，没有重新执行全部pytest；subagent重新执行的297 tests通过不能代替本次独立复核。

主助手在两个独立tmux任务中执行：

- 自己新写的边界／数值探针；
- 从现有微型输入生成器创建全新synthetic run，再测试stale raw、imported response、root symlink与convergence；
- 阅读迁移探针源码后，在新的临时目录重新执行，不消费subagent已经修改过的结果；
- 自己编写24维dense Cholesky oracle，与加速算法逐点比较；
- 重新推导频率、超辐射条件、phase积分、Stokes平均和截断正态例子。

证据：同目录 `parent-verification.log`、`parent-extra-check.log`、两个 `parent-*.py`。日志中的临时绝对路径仅为审查证据，不是生产路径。

分类：**确认缺陷**表示存在实现/接口/验证问题，不表示真实数据已触发；**条件性风险**表示近似在特定条件下会失效；**待生产证据**不能据现资料决定替代参数；**非bug**表示有意定义、历史兼容或资源选择。

## A. 确定实现问题：逐项结论

| ID | 问题与当前位置（相对项目根） | 主助手核验与裁决 |
|---|---|---|
| B1a | 模板inputs只有response NPZ，没有response manifest。`pipeline/workflow.py:165–176` | **确认**。`verify_stage`只沿manifest递归，来源图断边。独立最小图中模板通过recursive检查，直接验证响应却拒绝stale raw。 |
| B1b | inference额外检查processed/response时不递归。`pipeline/inference.py:98–104,133–136` | **确认**。新synthetic run中改α=.1 raw后，`run_inference(rebuild=True)`仍成功。只改raw不改变已保存模板的数值，但违反完整来源验证承诺。 |
| B1c | imported response祖先没有传播。`common.py:249–264` | **确认**。新run中将response标imported，无reuse授权仍可推断。新模板没有得到imported_ancestors。 |
| B2a | workflow先resolve输出root，使root symlink保护失效。`workflow.py:108,121–132` | **确认**。自己新建run别名后figures rebuild改变了真实目标manifest时间戳。可能是用户有意用别名，但与当前“拒绝通过共享目录链接重建”的保护语义冲突。 |
| B2b | `stage_directory`仅rebuild时检查目录symlink。`common.py:170–188` | **确认**。非rebuild通过悬空stage链接在外部目标创建产品。祖先链接同理。 |
| B2c | 低层输出writer仅exists，不拒绝悬空文件链接。`raptor.py:46–56` | **确认**。`out/stokes.npy`指向不存在的外部文件时，预处理创建外部文件。专家入口同样需要保护，不能只靠父stage临时目录。 |
| B3 | 科学YAML无完整严格类型schema。`workflow.py:34–54` | **确认**。两个`"false"`都加载成str且truthy。默认YAML正确；缺陷在自定义配置可能静默改方法。下层已有部分范围检查，但不是完整fail-fast schema。未知字段、计数bool、alpha精度等应统一验证。 |
| B4 | EHT圆均值R≈0仍输出有限中心与covariance。`eht.py:144–158` | **确认退化边界缺陷**。独立重现R约6.1e-17，一列整体加1e-12rad，cross-covariance从+.02945243变−.02945243。只证明退化输入问题，不能外推真实选点已受影响。宽分布阈值应另定，不能随手加任意物理cut。 |
| B5a | res400 glob字典排序不是数值时间排序。`layout.py:34–40` | **确认**。t1000排在t900前；当前固定四位文件可能不触发。 |
| B5b | shard没有校验重叠、缺口、文件帧数与标签步长。`layout.py`、`raptor.py:591–597` | **确认验证缺口**。background把拼接位置当等间隔，缺口可被当连续lag；没有真实缺帧证据。低层显式caller order不是bug，自动调度应承担排序验证。 |
| B5c | canonical单文件被legacy shards覆盖。`layout.py:34–40` | **确认优先级行为；条件性设计缺陷**。独立复现，但没有书面规则时不能说所有情况下选shards一定错。存在两个不同来源时，拒绝歧义或明确统一优先级更可靠。 |
| B6 | convergence只验证inference本地文件，读取当前未验证templates/observations。`convergence.py:25,35–36` | **确认**。自己将amplitude_mean乘2，工具仍成功并报告.29094/.27303dex“网格变化”。应核对来源并先匹配base重算。 |
| B7a | 迁移保存origin摘要，但verifier不处理该字段。`migrate.py:97–114`、`common.py:139–149` | **行为确认，需缩小原报告措辞**。重新运行迁移后改旧response manifest，新递归检查仍成功。复制的数值没有因此改变；静态origin收据本身不是错误。缺陷是验证只能保证新布局来源，不能声称已持续验证原完整科学来源。可用不可变原lineage快照，不必强制永远依赖可变旧目录。 |
| B7b | 迁移无条件生成1600+10*i标签，复制旧fit.json未明确标记新增标签来源。`migrate.py:139` | **确认metadata缺陷**。适用于标准完整输入的默认约定，不适用于任意subset。主fit_stokes已标default/unverified，migration没有同等说明。标签当前不进入拟合均值，但会误导时间身份解释。 |
| B8 | run迁移无整run事务／恢复协议。`migrate.py:82–295` | **确认恢复性缺陷**。缺08_figures后留下10个完整stage，无migration.json；同命令重试拒绝existing destination。逐stage事务仍有效，原输入未损坏。不应描述为结果已被静默篡改。 |
| B9a | RAPTOR/transport EVPA不检查QU finite，Inf可给有限角。`raptor.py:339–342`、`transport.py:218–221` | **确认输入边界缺陷**。自己的Q=Inf,U=1得到0rad。不能由最终角finite判断原Stokes有效。并非所有预处理路径都能让Inf存活；直接API已足以证明问题。 |
| B9b | integrate_stokes将每通道NaN独立变0，混淆padding与未知量。`polarization.py:190–200` | **确认输入语义缺陷**。I=1,Q=NaN,U=1被当45°。全通道padding当零通量合理，不应一律禁止NaN。EHT主入口会拒绝非有限raw pixels，因此此例不证明正常EHT主链已有坏像素。 |
| B9c | integer I*Q/I*U在转float前可溢出。`polarization.py:199–200` | **确认**。uint8 I200,Q2,U1给.47338664rad，float64给.23182380rad。主调度float64不触发，通用数值API仍有错。 |
| B10 | legacy posterior接受uniform但忽略prior。`limits.py:88–91,285–289` | **确认兼容接口缺陷**。自己验证两个prior都返回log10cut=1.8915783。历史等log网格CDF作为回归算法并非bug；应限制或诚实标注支持的测度。默认corrected主链不受影响。 |
| B11 | compare/diagnostics直接读产品，验证弱于canonical stage。`scripts/compare_original_limits.py:43–61`、`pipeline/diagnostics.py:14起` | **确认验证边界，不是全部绘图计算错误**。canonical load_curves缺finite/alpha单调校验；compare旧六文件分支仅检查observation，不统一核对prior/variant；diagnostics无manifest验证。比较本来可有意比较不同算法，不能要求variant必须相同，应显式标注差异。 |
| B12 | source signature不绑定运行库/lock版本。`common.py:35–105`、`workflow.py:289–291` | **确认复现性缺口**。NumPy/SciPy版本只进invocation元数据，不进stage缓存身份；库升级后旧cache仍可能命中。不等于已证明当前缓存数值错误。 |

## B. 物理、统计与常数问题：全部去重清单

### B.1 单位、图像、坐标与模拟来源

| ID | 问题 | 核验与裁决 |
|---|---|---|
| S01 | 当前M6.5e9/D16.8是否等于本批RAPTOR生产参数 | **待生产证据**。EHT-I采用值有出处，不能替代生产记录；另一EHT模拟库用6.2e9/16.9。质量、距离不确定度当前未传播，这是固定基准选择。 |
| S02 | FOV46.2262、32.6868、sigma2.16871是否统一物理标定 | **差异确认；实际正确值待查**。32.6868Rg当前换算124.8307µas，beam19.5033µas。比值高度符合128/20/2.355；历史raw相机给46.2261和另一M/D。不因匹配数值就确认当前pickle实际FOV。 |
| S03 | 108°、17°/163°、±.9375是否与天空/WCS/tetrad一致 | **待生产证据**。108与jet288反方向有关，但PA east-of-north不自动等于数组主动旋转角。自旋符号可随反射约定改变，不应单看负号判错。 |
| S04 | EHT transpose反射是否漏掉P共轭及V翻号 | **条件性风险，未证实bug**。空间map确为transpose，代码P'=-P；若是同一物理平面基底的反射应为−P*。但若只是存储坐标重排，sky Stokes可不这样变。缺WCS/Stokes基底前不可直接修。 |
| S05 | 标签1600+10i、每标签.1tg、141帧是否为真实时间 | **默认来源和形状确认；本批时间待查**。当前主fit metadata诚实标unverified。141点首末只有140间隔；M6.5下约51.877d，不是141tg。 |
| S06 | 旧N−1、.518838与当前N/FOV半像素网格不同 | **历史约定差异确认，非当前主链bug**。旧blur=False/早期Notebook混用中心跨度和边界宽度；不能只调M/D解释。 |
| S07 | 输入频率、flux单位、单色230GHz、fast-light、.6Jy与分辨率是否同批 | **待生产证据**。文件名/草稿不是运行配置；本轮未证明主cube和背景、transport同源。 |

### B.2 轴子态、耦合与超辐射

| ID | 问题 | 核验与裁决 |
|---|---|---|
| S08 | n4与211主量子数n2冲突 | **冲突确认**。model.py omega默认4；已发表氢样式主量子数为nr+l+1，211应2；旧result用(2+1+1)²。只证明代码/文字不一致，生产态可能不同；不要只改推断频率而不查空间模板。 |
| S09 | α=.35/.4/.44与自旋.9375,m1超辐射条件 | **条件不满足确认；解释条件化**。自己算ΩH=.34774119，当前ω4分别.34866016/.398/.437338，大于mΩH。不是禁止计算这些模板，可能是另一spin的云、预存衰减云或外推；不能统称该spin下增长211云。必要条件也不证明长寿命/占据量。 |
| S10 | AXION_NORM=.1与c、b、2π关系 | **历史意图可确认，实际注入待查；不得断言漏2π**。主助手在原EHT_ZC HEAD:result.ipynb cell36找到注释，声称ga=axion_norm/(2π)，解释除.1后乘bc。若producer确如此，现放缩可一致；subagent把“无生产证据”近似说成“无归一化定义”不完整。旧cell26也有check it，注释不是实际执行公式证据。 |
| S11 | 旧稿c含b而前作c不含b | **文字定义冲突确认**。必须决定纵轴到底是c=2πgf还是bc，不可仅靠文件名c95。S10历史注释支持不含b的c解释。 |
| S12 | b=.03(α/.05)^(7/4)适用范围 | **出处确认、外推风险确认**。自己读Baryakhtar Eq80，针对特定最大场幅、fBC附近、两能级/小α标度；不是任意时间、任意f、整网格的无条件精确云幅度。现值需作为基准假设，不证明“b>1所以必错”。 |
| S13 | 参考注入响应线性外推到c1000 | **条件性风险**。prior上界大不证明信号仍弱；需检查后验实际质量区间、多注入强度、高谐波及Stokes退偏振。不能仅凭1000就断言真实cut已非线性。 |
| S14 | producer16相位顺序、零点、时间方向、snapshot绝对相位 | **待生产证据**。共同绝对相位可吸收，但空间依赖或时间推进方向/已经嵌入的snapshot相位未必可吸收。DFT数学正号不能证明producer正号。 |

### B.3 观测样本、mask与协方差

| ID | 问题 | 核验与裁决 |
|---|---|---|
| S15 | 4800→1572不是EHT-IV原Top Set规则 | **确认**。自己核对EHT-IV是37500→1572，当前按EVPA中心性排名，不用visibility χ²。1572的灵感有出处，概率解释未成立。 |
| S16 | 先筛中心样本再估covariance的选择偏差 | **机制确认，真实偏差未量化**。自己复算单变量正态中心保留32.75%后σ比.24116；不可把此例当本项目高维covariance缩小倍数。需noise/calibration injection coverage与topset敏感性。 |
| S17 | 超参网格等权散布不自动是测量噪声后验 | **统计解释缺口确认**。它可作为方法散布代理，但不是天然iid posterior；不应改ddof或除sqrt1572假装解决。 |
| S18 | 观测跨日协方差设0 | **假设确认，真实性待估计**。一般差分C=C1+C2−C12−C21。共享校准、Dterms、共同筛样本可能相关；身份匹配既不证明相关也不证明独立。 |
| S19 | 模型误差空间对角、两个daypair独立 | **假设确认，真实性待估计**。与观测跨日独立是两件事；beam/湍流/慢变会产生相关，不能预设忽略一定保守或一定乐观。 |
| S20 | 历史mask r≥2、10%I、5列/矩形排除 | **来源明确、科学选择理由待查**。EHT-VII的10%显示cut不等于本统计mask；r是屏幕半径，不是BL事件视界；不能将人为mask当观测常数。 |
| S21 | 固定两阶段offset、1d不继承2d mask、offset4 | **实现/历史产物确认，非自动bug**。需要配准诊断和敏感性。10°非重叠扇区不保证beam/成像误差独立。 |
| S22 | ridge0与高维covariance条件数 | **稳定性风险**。有Cholesky失败保护和保存condition numbers，不是完全没检测；缺实际敏感性/近奇异稳定性证明。 |
| S23 | 圆周单branch Gaussian适用于宽、多峰/低P分布吗 | **条件性风险**，与B4数值退化分开。当前同branch计算正确，但不能说明宽分布高斯；应诊断resultant、近零偏振及branch稳定性。 |
| S24 | EHT-VII Eq11印刷QI/UI与当前atan2(UI,QI)不同 | **原文差异确认，不据此判当前错**。标准P=Q+iU给atan2(U,Q)/2；单个印刷式不能取代真实基底/象限核对，也可能是原文笔误。 |

### B.4 模型分布、时间平均与背景

| ID | 问题 | 核验与裁决 |
|---|---|---|
| S25 | A/phi delta propagation未含Cov(A,phi) | **假设/近似确认**。源码仅两项，完整一阶传播有交叉项。需联合分布，不一定只补一个标量即可解决空间/跨时相关。 |
| S26 | 宽相位下用meanA与meanphi构造nominal signal | **近似偏差机制确认**。自己normal例E[sinphi]≈.25972，而sin(meanphi)=.38942；解析衰减e^(−s²/2)与MC一致。数据实际散布待诊断。 |
| S27 | 141时间片散布当模型误差，样本自相关/不同时刻有效样本不一致 | **分析定义及风险确认**。用scatter不是mean standard error；不能直接除sqrt141。A/phi计数不同、近零幅度相位不稳定需诊断。 |
| S28 | 平均标量EVPA不等于平均Stokes/visibilities后取角 | **不等价确认，条件性风险**。自己例平均Stokes角.50848rad、平均EVPA .49000rad。小旋转稳定P时一阶可成立；不证明真实posterior中已有此量级偏差。 |
| S29 | 1d average(I P)与Ibar Pbar不等价 | **代数确认、实际重要性待查**。当前日图observable与动态模拟时间平均需同定义。 |
| S30 | 6h矩形窗口、24h/120h名义间隔不等于实际scan weights | **近似来源确认**。前作使用6h；实际EHT离散扫描，April10覆盖更少。需scan级时间和权重，不应说sinc公式错。 |
| S31 | background lag3=26.6797h，不是24h | **差异确认，非当前隐藏bug**。代码和metadata已经说明；需要结构函数/lag敏感性，不能无条件线性缩放或除sqrt2。 |
| S32 | 背景没有同样6h曝光平均 | **观测量匹配风险确认**。与S31分开，需Stokes时间平均。 |
| S33 | 同一个2d pooled σ用于1d/2d | **近似确认**。1d应可在aperture层另估；pooling异方差、非Gaussian、空间相关也未由标量表达。可作phenomenological选择，非算术错。 |

### B.5 后验与收敛

| ID | 问题 | 核验与裁决 |
|---|---|---|
| S34 | log-uniform有限prior边界影响null上限 | **确认方法属性，非bug**。1/137、1000、95%有前作出处；零响应仍给prior-driven有限cut。需边界敏感性，不把边缘mass当证明。1/137与精确αEM小差不是主要问题。 |
| S35 | 64 phase不普遍收敛、一次doubling不够 | **自己复现条件性风险**。合法α=.23329159，64/128/连续相位logcut=−.047991/−.023442/−.014219。不是默认八点，不能据此宣布真实曲线错。 |
| S36 | convergence入口默认不包括cor_grmhd | **确认覆盖缺口**。模块默认三mode，但scripts/check_convergence.py:13默认std/cor；没有coupling/phase联合加密、多个phase原点、积分网格/beam/mask敏感性或接受阈值。B6修好前已有报告更不能作为完整收敛凭据。 |

### B.6 transport旁路、历史预处理与显示

| ID | 问题 | 核验与裁决 |
|---|---|---|
| T01 | common/jonly/eq_jonly/fake真实开关未验证 | **待生产证据**。代码已标inferred；不能定量宣称plasma/LP百分比。 |
| T02 | 用r_coord圈选像素不等于沿geodesic关闭发射项 | **物理区别确认**。mask先后顺序回答不同问题，不可互换解释。 |
| T03 | 先求场幅再插值与先插值/模糊Stokes后求幅不交换 | **自己合成验证确认**。相邻同幅反相，插值幅度图为.08，插值Stokes响应为0；所以unblurred ratio≠1也不一定是传输物理。现测试常场不能覆盖。 |
| T04 | 单位spin2 phasor插值不是带偏振强度的Stokes插值 | **定义差异确认**。用于角标签合理，不能称物理重成像完全等价。 |
| T05 | 默认legacy RMS包括高谐波、16态中心off-by-one | **确认且有意兼容**。自己纯cos例比真幅1.00443461，另相位1.0。不能当精确第一谐波，DFT已另存；不应破坏legacy回归来掩盖。 |
| T06 | photon公式选正spin下较大逆行赤道轨道，不是一般倾角critical curve | **公式确认，实际分支待几何约定**。a=.9375为3.9441，小顺行1.4345；旧负spin/reflection不能忽略。 |
| T07 | BL4比值颜色贴photon-reference曲线 | **确认proxy，已标注，非偷偷bug**。不能说是该曲线的实测幅值。 |
| T08 | crossing first-outward只扫描有限网格，可能漏切触/窄双根 | **数值方法边界确认**。并非所有物理r_coord都触发；分支是否对应所需赤道交点待查。 |
| T09 | radius geometry NaN→0、旧field NaN→0/失败(200,0,0)可能已烘焙进数组 | **历史代码行为确认，当前产品污染未证实**。需validity mask和生产日志；不能把零场和失败混同。 |
| T10 | 1022、.29915625、observer r500及取最后theta<π/2的位置 | **历史表达式确认，物理身份待查**。.29915625=omega4(.3)；1022原因未知。r500为几何单位非距离Mpc。最后True不自动等于第一次赤道交点。 |
| T11 | analytical screen用common[0]，common_noaxion另作背景，两者可能不同批 | **代码选择确认**。已保存差值但未证明可互换，需要producer解释；本轮未读取大数组比较。 |
| T12 | source图固定±16.3434、比值色标.9–1.1、强度−4到0 | **显示选择确认**。自定义FOV时图范围仍固定；可裁图/压色，不代表科学效应≤10%。 |

## C. 报告中不是问题的数值，以及应纠正的措辞

1. IAU nominal GM、c、parsec、hbar换算正确。当前主链时间系数派生为.1124448628，旧.112411只在低层默认/legacy中；两者差约.0301%，不是当前主链仍硬编码错单位。
2. 64²是匹配网格；16相位/17态是生产形状；64 nuisance phases是积分网格，不能互相混称。/8是DFT 2/N，不是额外物理常数。
3. 10°窗口和1°移动有EHT-VII支持，强度加权I(Q+iU)不是必然重复加权错误。
4. null=0是用户要求的条件性研究中心，保留观测covariance正确，不应被“修回实测日差”。不等于EHT测得每像素日差为零，也不等于median expected limit。
5. ddof0可定义经验分布/MLE，不必机械改N−1；time scatter不应机械除sqrtT；背景日差σ不应机械除sqrt2。
6. 数值容差、hash块8MiB、batch64/128、LRU8、BLAS2线程、DPI/颜色等多为工程设置，不是论文物理常数；应测数值敏感性而非逐个寻找天文出处。
7. 背景旧1.975261°是历史histogram拟合产物，当前主链重新估计；旧阈值.2、.053、Lilliefors门槛不是当前MLE都在使用的gate。
8. 精确eigen变换没有发现公式bug。自己24维oracle、covariance谱跨度2e4，最大log-likelihood差3.64e-12。它确认此例数值一致，不证明任意770维病态矩阵都稳定。
9. 固定covariance省略共同logdet正确；变化covariance保留logdet正确；共同初相位下两组likelihood先乘后边缘化正确。
10. 旋转RAPTOR spin2的内部代数正确，但与真实天空基底对齐还需要S03/S04。
11. EHT-VII引用位置需纠正：120µas/64²为Appendix C.2，窗口为§5.2；早期文档写D.1/§5.1不准确。不能把Eq11印刷分式反向就直接认定当前atan2错。
12. 旧草稿999 snapshots/1234days/888pairs是明确占位，不是可用的数据事实；不应传播到结果说明。

## D. 性能与可维护性建议（不是正确性bug）

| 建议 | 主助手判断 |
|---|---|
| std/cor用2×2 sin/cos Gram与I0e连续边缘化phase | **数学成立**，自己独立用它作phase oracle。固定covariance、null、共同均匀相位才可直接用；不能原样套phase-dependent cor_grmhd。是精度/性能增强，不是当前离散法必错。 |
| null利用phi+π对称，只算一半phase | **成立**，模板反号、方差不变。限无重复端点的匹配偶数grid；不能无条件套所有legacy覆盖开关。 |
| batch累计同branch covariance | **可行方向**，应保留同center、同有效样本，避免改变missing-data规则。 |
| invocation内复用hash验证 | **可行但有并发条件**。不能mtime替代checksum；输入可变时要避免TOCTOU，不能无限信任_seen。 |
| transport复用blur/crossings/analytic maps | **可行**，mask-before-blur不能错误地和blur-before-mask交换；sigma0的结果可共用。 |
| RAPTOR缓存几何coords，背景ring buffer | **可行**，目前多遍计算存在重复；只优化不改观测量。 |
| 记录lock/runtime identity | **有必要，见B12**，当前invocation有部分版本但stage身份不绑定。 |
| 隔离legacy兼容函数和直接诊断入口 | **维护建议**，未被主链调用不等于无用。_legacy_quadratic和marginalize_log_likelihood没有调用；oracle函数有保留价值；transport/HDF5是有效专家入口。 |

## E. 全部测试缺口（不是声称真实数据已出错）

1. 完整raw→processed→response→template来源，imported中间节点与重算祖先传播。
2. 输出root/祖先/叶目录链接、悬空文件链接及只读cache与写入分开保护。
3. migrate_run/migrate_inputs全链测试；late failure恢复、磁盘满、rename失败、源变化。
4. YAML严格类型、未知字段、bool/int边界、alpha身份精度和数据范围。
5. shard数值排序、缺帧/重叠/长度和双layout歧义。
6. resultant0/近退化、空ensemble、部分NaN/Inf、整数与近零P。
7. 独立WCS/sky spin2反射oracle，不用同一转换逆生成输入自证。
8. 多coupling线性域、DFT高谐波残差、producer相位/时间方向、A/phi联合散布与时间相关。
9. Stokes层曝光、真实scan weights、1d背景observable匹配。
10. topset noise-injection coverage、跨日systematics、mask/topset/ridge敏感性。
11. 中高维病态covariance、连续phase oracle、不同phase原点、prior边界扩展、全部六曲线收敛。
12. transport空间变幅/变相位/变偏振强度与真实开关控制，而非仅常field。
13. compare/diagnostics非法产品、方法身份差异、stale来源拒绝测试。

## F. 下一步优先级（本轮未实施）

先修B2输出安全、B1来源断边、B3严格schema、B5时间排序/连续性、B6收敛身份。B7/B8需要先定义迁移来源“静态归档”与“持续依赖”的契约，再完善验证/恢复，不能随意删除用户半成品。

原文件到手后先核实S10/S11耦合纵轴、S08频率态和S01–S04真实图像基底，再评估筛样本概率解释、幅相联合传播、曝光/背景/相关性。

结论不是“两个subagent报告全对”，也不是“所有问题都能由改常数解决”：大多数确定实现缺陷自己可复现；物理项主要是已确认的约定冲突和有条件近似，仍没有证据断言本批真实六曲线已偏移多少。
