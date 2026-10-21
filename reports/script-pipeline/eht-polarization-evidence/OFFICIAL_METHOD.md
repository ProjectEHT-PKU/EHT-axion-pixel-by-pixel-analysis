# M87* 2017 EHT 官方极化图与 Top Set：独立原文核验

## 1. 范围、证据与执行边界
- 本次先联网查官方论文与 EHT 官方代码组织；只读现项目源码作静态对照，不运行成像、拟合、科学统计，不加载原始 pickle，不修改项目源码或原 data。
- 主证据：[VII, arXiv:2105.01169](https://arxiv.org/abs/2105.01169)，正式刊 DOI [10.3847/2041-8213/abe71d](https://doi.org/10.3847/2041-8213/abe71d)，ApJL 910:L12，48 页。
- definition 对照：[IV, arXiv:1906.11241](https://arxiv.org/abs/1906.11241)，正式刊 DOI [10.3847/2041-8213/ab0e85](https://doi.org/10.3847/2041-8213/ab0e85)，ApJL 875:L4，52 页。
- 已先 curl HTML/source，均设置 25–60 s timeout；VII HTML 无可用正文，arXiv `/src` 和 `/e-print` 实际返回 PDF，未把它们假称 TeX。
- PDF 用 `pdftoppm -scale-to 1800 -png` 渲染后直接读页图，未用 pdf2text、OCR 或自造 legacy extractor。
- 原件与全部渲染页：`/tmp/eht-official-method/{PaperVII.pdf,PaperIV.pdf,pageimages/VII-01.png…VII-48.png,pageimages/IV-01.png…IV-52.png}`。
- 本报告页码均为正式刊页脚，也与本地渲染页号一致；未深读 VIII，因为成像与筛选问题在 VII/IV 已有一手证据。

## 2. 必须先分清的目标图：mean 不是“第五天”
|目标|准确图号/位置|caption 与方法含义|
|---|---|---|
|各 pipeline × 四天主极化图|VII Fig.6，p.11；§5.1 p.10|五行依次 eht-imaging、polsolve、LPCAL、DMC、THEMIS；四列 April 5/6/10/11；low-band fiducial images。|
|每日官方保守平均图|VII Fig.7，p.12；§5.1 p.10|对每一天的五种已等效模糊 fiducial 图做简单平均；I/Q/U 分别平均。四列仍是四天，不是四天合成一个 mean。|
|方位分布与日日差|VII Fig.8，p.13；§5.2 pp.10–12|April 5 对 April 11，五种方法、每方法 1000 个 D-term/posterior 样本的方位 P、EVPA 分布，不是 Top Set EVPA 中心筛选图。|
|积分 EVPA 全四天|VII Fig.10，p.16；§6 p.15|low/high band 各方法净 EVPA，4/10 约束差，不能把误差条当像素噪声。|
|Stokes I 官方 consensus|IV Fig.15，p.21；Fig.14 同页|每日三种 fiducial I 图在等效分辨率下平均；是 VII 选模糊核与定位所用参照。|
- Fig.7 上排为 tick plot；下排为橙色总强度上的 EVPA streamline 可视化。caption 明确 streamline 的长度和 opacity 按 polarized intensity 的平方缩放，灵感来自 LIC；不是直接磁场线，也不是额外观测量。
- 本项目若增加“跨四天 mean”面板，应标明自定义 temporal mean，不能称为 VII Fig.7 的官方第五面板。

## 3. 图像量、mask、tick 和单位：可执行的显示合同
- VII §2 Eq.(1)，p.3：复线偏振 `𝒫=Q+iU=I|m|exp(2iχ)`；标量 `P=|𝒫|=sqrt(Q²+U²)`，`|m|=P/I`，`χ=0.5 arg(Q+iU)`。
- χ 从北向东测量，模 180°；实现必须 `0.5*atan2(U,Q)`，不能普通 atan 单参数、更不能把 P 当有符号 Stokes。
- Fig.6：底图是灰度 I；tick 方向为 EVPA；长度为 P；颜色为 `|m|`，色条 0–0.30；P 等值线实/虚/点分别为 20/10/5 μJy μas⁻²。
- 绝对 tick 标度：10 μas 的 tick 长度约对应 30 μJy μas⁻²；caption 明确 **不对每张图的最大 P 重新归一化**。每日 peak-normalized tick 会掩盖真正幅度变化。
- 显示 FOV 120 μas；标尺 50 μas。Fig.7 上排沿用 Fig.6 的 I、fraction、EVPA 显示规则，I 色标为亮温，约 0–5×10⁹ K。
- caption 原句：`Cuts were made to omit all regions in the images where Stokes I < 10% of the peak brightness and P < 20% of the peak polarized brightness.`
- **mask 细节 unknown**：这里文字使用 and；仅靠 caption 不能严格确定作图代码究竟使用“低 I 且低 P 才去掉”，还是常见的“任一低于阈值即去掉”。复现须保留该原句，找到原 plotting code 前不要静默把 AND 改 OR；可并列展示两种实现并明确差异。
- 以上是显示 cut，不是 visibility 拟合筛选条件，也不是计算全图积分 polarization fraction 的默认掩膜。
- 若 FITS 是 Jy/pixel，先除以像素面积 μas²，再乘 10⁶ 得 μJy μas⁻²。若为 Jy/beam，应先除以高斯 beam 面积 `π bmaj bmin/(4 ln2)`，不能按 pixel 面积误换。
- I 的亮温可由 Rayleigh–Jeans `Tb=c² Iν/(2kν²)` 换算，注意 Iν 必须是 Jy/sr；显示与对比要记录实际 band ν，而非凭图片颜色调整 flux。

## 4. beam / resolution 与正确平均次序
- VII §5.1 p.10：对各方法的 I 图与 IV April 11 consensus I 图做 normalized cross-correlation，选择使相关最大之圆高斯卷积核。
- **VII 最终核 FWHM**：eht-imaging、DMC、THEMIS 各 19 μas；LPCAL 20 μas；polsolve 23 μas。同一方法对 I/Q/U 使用同一核（VII §5.3 p.14）。
- 这是一组“匹配有效分辨率”的核，不是“给所有图加相同 20 μas beam”，更不是像素大小等于分辨率。约 20 μas 是结果的有效尺度描述。
- VII Fig.21 p.30 是 preliminary：eht-imaging 17.1 μas、CLEAN 20 μas；不能冒充 Fig.6 的最终 19/20/23 μas。
- IV Fig.14 p.21 的 eht-imaging 17.1、SMILI 18.6、DIFMAP 20 μas 也不等于 VII 的最终极化规则。
- 正确 Fig.7：先统一网格、坐标/位置并按方法卷积各 Stokes，再 `Ibar=ΣI/5,Qbar=ΣQ/5,Ubar=ΣU/5`，最后派生 `Pbar=sqrt(Qbar²+Ubar²),mbar=Pbar/Ibar,χbar=0.5 atan2(Ubar,Qbar)`。
- 不可先平均 χ，也不可平均 `sqrt(Q²+U²)` 代替平均后的 P。即便是 `arg mean exp(2iχ)/2`，也只是等权方向平均，不等于 amplitude-weighted 的平均 Stokes 方向。
- 若是 single-pipeline hyperparameter ensemble，仍可同样做 Stokes 平均，但必须标成该 ensemble mean；它不是“五个独立方法各一票”的 Fig.7。
- 已被模糊的 FITS 不可再盲加上述核；先查 native/restored provenance。高斯近似下新增核满足 `badd²=btarget²-bexisting²`，但这不能代替论文的 cross-correlation 匹配。

## 5. centering / alignment 必须同作用于 I/Q/U
- VII §4.3 p.8 Eq.(8)–(9)：synthetic 对比先平移图像，最大化 Stokes I 相关；I 是减均值 normalized correlation，复偏振相关取实部且不减均值。
- VII §5.2 p.10：方位分析中心由 I 图与自校准使用的 representative M87 I 图的相关峰估计；不能由 EVPA 中心性决定图像空间位置。
- VII §5.3 p.14：β₂ 的中心通过 blurred I 与 IV April 11 consensus I 的相关最大像素位移确定；论文称位移通常小，不定位时 β₂ 结果也无显著改变。
- IV §6.3.1 p.16：因绝对 visibility phase 不知，逐一离散平移计算 ρNX，取最大值；CLEAN 用 20 μas 恢复后参与该比较。
- 工程 recipe：先以明确 sky/WCS convention 重网格；只通过 I 求 `(dx,dy)`，同一偏移施于所有 Stokes，保存 shift、参照、边界损失。禁止分别对 Q/U 配准或为抹平日日差旋转 EVPA。
- 图像翻转/旋转不是纯数组操作：Q/U 是 spin-2，必须验证所用天球基底；用已知北向/东向 EVPA 的小型坐标单元测试，不用科学结果来“调对方向”。

## 6. Paper IV 的 Top Set：1572 的真正含义
- IV §6.3 pp.14–17：Top Set 是 **参数组合集合**，需同时通过真实数据 fit 与 synthetic recovery 两阶段要求；不是按真实图像 EVPA 接近 ensemble 中心排序。
- 数据门槛：scan-averaged low-band，closure phases 与 log closure amplitudes；只从起始时间同步的平均 visibilities 构造 closure。任一天任一要求的 χ²≥2，则该参数组合被排除（p.15）。
- DIFMAP 特例：门槛评估前给每条 visibility 热噪声另加 10% systematic，quadrature 合成；不应把它套到所有 RML 或 polarimetric pipeline。
- IV §2.1 pp.3–4：这些 χ² 是 mean squared standardized residual；没有扣 effective image degrees of freedom，**不严格是 reduced χ²**，不能用常规 χ² 分布自由度作置信度解释。
- synthetic 四个几何训练源：ring/crescent/disk/double；每模型取跨日 ρNX 中位数，再用各模型的 Gaussian-smearing correlation 曲线换算 equivalent FWHM α。
- 取四源 α 的平均得到 αchar；Top Set 需 `αchar < 1/|u|max`，即优于 nominal EHT resolution，约 25 μas。Top Set 内 αchar 最佳的一组参数定义每方法 fiducial。
- IV Table 3 p.17：DIFMAP 1008→30；eht-imaging 37500→1572；SMILI 10800→529。**1572 是 2019 Stokes I 特定参数 survey 的结果，不是极化图需固定保留的样本数。**
- 同表百分比指 Top Set 内含某参数值的比例，不是“该值在全部候选中通过的概率”；fiducial 由联合表现最优选择，不是逐参数取众数。
- IV Fig.10 p.16、§6.3.2 p.17：留一模型检验与 GRMHD 测试用于检查选择偏好；也承认训练源有限、compact extent/细结构/flux 假设的局限。

## 7. Paper VII 的极化 RML fiducial：不能机械复制 IV Top Set
- VII §4.3 p.8：每方法在 synthetic 上选最佳参数，再固定该参数重建真实四天；其对应图称 fiducial images。五方法为 eht-imaging、polsolve、LPCAL、DMC、THEMIS。
- 六套 polarized synthetic 使用 **April 11 low-band (u,v)**；两个 GRMHD、一个 filled disk、三个不同 EVPA pattern 的 ring（Fig.4 p.9）。不能偷换成仅一个环模板或真实图自身。
- synthetic visibilities 加 thermal noise、phase/gain offsets、D-terms；SMT/LMT/PV 的 mock leakage 参考初步观测量级；ALMA/APEX/JCMT/SMA 加 1% residual D-term normal scatter，再用 SMILI fiducial I 自校准，模拟真实处理（p.8）。
- VII App.C.2 pp.22–23：scan average；I 重建用 IV eht-imaging fiducial，FOV 120 μas、64² pixels、总 compact flux 0.6 Jy；重新自校准后极化阶段去掉 zero baselines，I 在每个 polarized imaging step 固定。
- 目标函数是数据项 αP/αm 与 HW、polarimetric TV 正则的加权组合（Eq.C4/C5）；同时使用 polarized visibility 与 visibility-domain polarization ratio。不能把图像 EVPA dispersion score 称作该 likelihood。
- HW 偏好 `m<0.75`，但 **0.75 不是 hard cut**；另用变量变换使 `Q²+U²<I²`。visibility ratio 幅度可大于 1，不可用 image-domain m≤1 将其裁掉（§2 p.3）。
- imaging 与完整 Jones D-term 解交替，不线性化 leakage；假设 V=0，常需 50–100 轮稳定（p.23）。重启可用前一轮图经 20 μas 模糊作初始值；这是优化步骤，不是最终显示 beam。
- VII App.G.1/Table 6 p.31：固定 αm=1；αP、βHW、βTV 各扫描 `[0.01,0.1,1,10,100]`；niter `[1,10,50]`；fsys `[0,0.002,0.005,0.01]`。
- **最终极化 fiducial**：αP=1、αm=1、βHW=100、βTV=1、niter=50、fsys=0.005；fsys×|I-visibility| 与 thermal error quadrature 合成。
- fidelity score 用 **标量 P=sqrt(Q²+U²)** 重建与 truth 的 ρNX（G1）；D-term score 用 SMT/LMT/PV 平均复平面 ℓ2 distance dD，`sdterm=1-erf[dD/(sqrt(2) dtol)]`，dtol=5%（G2）。
- 对六 synthetic 把 `sfid,a × sdterm,a` 相乘（G3），取总分最高参数；dtol=5% 是 soft score 尺度，不是每个 D-term 都必须小于 5% 的拒绝阈值。
- 在这些已核读 VII 极化 RML 方法段中，**没有定义“按 IV 的 χ²<2 加 1572 个极化 Top Set”**；不可把 IV 的 I 筛选和 VII 的 polarized synthetic score 合并后声称官方。
- 正则权重依赖实现中的 normalization、data term、pixel/FOV 和软件版本；数字相同并不足以实现同一个 objective。Eq.C5 的符号应与代码 entropy/TV 定义配套核对，不凭重排印刷公式猜符号。

## 8. acceptable polarization rates 与误差样本到底是什么
- VII Eq.(12)–(13) p.13：`|m|net=|Σ(Q+iU)|/ΣI`；`<|m|>=Σsqrt(Q²+U²)/ΣI`。前者是 coherent net，后者分辨图的 intensity-weighted fraction，不能混称平均偏振率。
- beam 增大会降低 `<|m|>`；理想全图 flux-conserving convolution 不改 `|m|net`。裁剪、mask、边界会破坏该不变性。
- VII Table 2 p.15 给模型评分范围：`|m|net=1.0–3.7%`，`<|m|>=5.7–10.7%`，`|β₂|=0.04–0.07`，`arg β₂=-163°…-127°`。
- 这些范围来自 April 5/11 × 五方法，各均值±1σ 的最小/最大包络；**是后续 GRMHD 评分约束，不是 RML 像素 cut、Top Set 训练门槛或每一天必须满足的先验。**
- 约 15% 是约 20 μas 图像西南部 resolved peak fraction（§5.1），非 net fraction；Fig.6 色条到 30% 也不表示 ring 普遍 30%。
- VII §5.2、App.H pp.33–34：每方法/每天 1000 图；非 MCMC 方法随机抽 D-terms 后固定该校准再重建 Q/U，不同时再解 leakage；DMC/THEMIS 抽 posterior。
- D-term sampling 的跨站与 L/R 独立假设是 worst-case 风格测试，不是完整 posterior；eht-imaging 保留同一 D-term 实虚部 covariance，另两 imaging 方法没有。不能把同一参数 identity 跨日配对解释成官方 posterior temporal covariance。
- Fig.8 方位窗宽 10°、每 1° 滑动；PA 从北逆时针；χ 需对 Q/U 做强度加权 coherent 合成以避 ±90° wrap，不是算术平均角度。

## 9. 官方对跨日变化的措辞与可复现诊断
- VII §5.1 p.10：4/5–6 对 4/10–11，南侧 polarized emission 后期较弱；变化较小，局部细节不易与覆盖差异分离。
- §5.2 p.12：P 峰方位从 4/5 `199°±11°` 到 4/11 `244°±10°`（跨方法平均）；西南高亮区 EVPA 结果一致性较好，北部弱光区方法差异大；中心的小偏移也可能影响方位 χ 对比。
- §6 p.15：变化与 VLBI 数据中的 intrinsic variability signature 一致，不太可能仅由每日不同 (u,v) 覆盖导致；但不等于每个像素差异都有物理解释。
- App.B.4 Fig.14 pp.20–21：closure traces 给 station-based calibration-insensitive 的 evolving source structure 证据；可能主要追踪 I。
- App.B.5 Fig.15：polarization-sensitive conjugate closure trace product 有跨日演化迹象，但显著 quadrangle 太少，结论 **suggestive at best**。不能写成校准无关地高显著检出完整 polarization movie。
- IV Fig.23 p.26：对 4/6 与 4/11 仅重叠 (u,v) 重建，I 变化仍在；这支持 intrinsic I evolution，不自动证明所有 Q/U 差图。
- 执行匹配：统一坐标与 beam 后，分别输出 ΔI/ΔQ/ΔU、由各日 Stokes 推导的 ΔP、`Δχ=0.5 arg[(Q2+iU2)(Q1-iU1)]`；仅在明确共同有效显示区解释 Δχ，并同时标出低 P 区。
- 日日差须和各自方法/校准样本散布并列；single pipeline 只能给条件性差异与敏感性，不能覆盖五方法系统误差。4/10 稀疏覆盖导致长尾和弱约束，应单独标警，不强行套 4/5/11 区间。

## 10. 当前项目 TopSet 的只读静态诊断（不采信现处理报告为 oracle）
- 直接读取：`eht-script-pipeline/eht_axion/pipeline/eht.py` 与 `configs/production.yaml`；未依赖既有审查报告的结论，未运行模块。
- `FIELDS=(zbl,prior_fwhm,sys_noise,MEM,TV,T2V,l1)` 是 I 型 hyperparameter identity；源码先核四日 identity 集合一致，然后每像素全 ensemble χ 的 branch/mean/std，跨四日累加 `Σ[(χsample-χcenter)/σχ]²`，`argsort` 后固定取前 `topset_size=1572`（约 L430–483）。
- 这里 score 衡量 **EVPA ensemble 中心性**，没有观测 visibility residual、closure χ² 或 synthetic truth recovery。故它不是 IV Top Set，也不是 VII fiducial selection；即使最终数目相同也没有定义等价性。
- `RANK_OFFSETS` 与 `FINAL_OFFSETS` 为两套硬编码偏移；report 明写 `fixed; never refitted`，不是本文所查官方 I cross-correlation centering。
- `build_mask`（约 L254–293）是四天 I>10% peak、finite spread 的交集，加 r≥2 Rg、前 5 列及历史矩形排除；没有 Fig.6 的 P 相对阈值。它是研究分析 mask，不能直接宣称官方显示 mask。
- `summarize_evpa` 与最终 `summary.npz` 的 EVPA mean 是角分布摘要；Fig.7 要先平均 I/Q/U，再算角。源码存了 per-image relocated Stokes，可在另行授权的绘图工作中构造正确的 single-pipeline Stokes mean，无需把现 EVPA summary 当 oracle。
- 选择靠中心样本，再用所选样本估 dispersion/covariance，存在条件选择压窄散布的机制；这里只能指出方法风险，未量化本项目实际偏差，不能宣称官方 1σ coverage。
- 必须补 provenance：上游 Q/U 如何生成、D-term correction stage、完整参数表与版本、是否 polarization synthetic trained、beam 是否已加、单位/WCS、每图 calibration identity。此次不打开大 ensemble，不假设它一定是 4800 候选或来自官方公开包。

## 11. 官方下载、license 与已经核实/尚未核实的边界
- 官方入口 `https://eventhorizontelescope.org/for-astronomers/data` 本次返回 HTTP 403（含 www/带查询重试也失败）；不能由入口受阻推断“官方没有公开”。
- **确认可下载 polarized UVFITS**：[EHT 官方 2023-D01-01](https://github.com/eventhorizontelescope/2023-D01-01)，下载 README/tree/INVENTORY/LICENSE，tree commit `b71399714fdfe13579c741075d757b4ce65fd3d1`。
- README 明确 2017 M87 四天×high/low×三处理版本=24 数据集：HOPS `ALMArot` 无 D-term correction；HOPS `zbl-dtcal+selfcal` 已校 ALMA/APEX/SMA/JCMT leakage、并以 SMILI I 自校准；CASA 已校全部 D-terms，未 selfcal，R/L gain 假设 V=0。
- 所有版本已有 ALMA absolute EVPA 与 parallactic rotation。与 VII RML 起点最接近的是 HOPS zbl+selfcal；不要对 CASA 重复做相同 full D-term correction。
- 小文件例（181440 bytes，仅核 inventory/tree，未下载科学数据）：`hops_data/April10/SR2_M87_2017_100_lo_hops_zbl-dtcal+selfcal.uvfits`；可经 raw.githubusercontent.com 的该仓库/固定 commit/路径下载。
- LICENSE.txt 首段：data 为 **ODC-PDDL 1.0**，Python/shell 为 **GPLv3**。README 的 IX DOI 仍是占位符，引用 VII/IV 应使用上面的已核 DOI，不能照抄占位 DOI。
- **2019 UVFITS 不能恢复极化**：[2019-D01-01 README](https://github.com/eventhorizontelescope/2019-D01-01)：RL/LR 及误差置零，RR/LL 均为 I；这不是原始 polarimetric data。
- **确认公开 I imaging scripts**：[2019-D01-02](https://github.com/eventhorizontelescope/2019-D01-02)，DIFMAP/eht-imaging/SMILI；README 提醒公开 I-only 数据与论文 dual-pol、依赖版本会使数值略异。不能称为 VII 完整极化五 pipeline release。
- VII p.6 脚注直接链接 [achael/eht-imaging](https://github.com/achael/eht-imaging)。已查 tree/示例与 LICENSE：含 polarimetric imaging、pol_cal、example_pol.py、tutorial；GPLv3。**通用库/教程不是 Fig.6 的完整参数 survey 与作图代码**，也不是五种方法所有内部脚本。
- **官方 Fig.6/7 image FITS、五方法 fiducial/1000-sample ensemble、原 plotting code：本次未取得可验证下载清单**；不能把 UVFITS 当 image FITS，也不能把本文 PNG 当 FITS 真值。
- `ehtpublic2017pol`：Harvard Dataverse 搜索无此记录，alias API 返回 not found；测试的 achael 同名 GitHub 返回 404；搜索引擎受 challenge/异常结果影响。因此本次 **未核实该名称的实际服务/URL/license**，不把猜测路径写成可用官方资源。
- CyVerse 猜测路径只返回通用 Angular shell，无清单；不足以证明 dataset 存在。下一步应从恢复可访问的官方数据入口或作者给定精确 URL 获取 tiny README/manifest/FITS header，再验证 provenance；不要继续下载猜测目录的大包。
- 论文 VII p.1 明示 **CC BY 4.0**：复用原图应保留作者、标题、journal citation 与 DOI；不要把论文图 license、代码 GPL 与数据 PDDL 混为一个许可。

## 12. 下一轮最小 matching recipe / 验收清单（本次未执行科学跑）
1. 确认目标：先匹配 VII Fig.6 的 **eht-imaging 行**，再谈五方法 Fig.7；每日四幅与自定义 temporal mean 分开命名。
2. 先拿官方 image manifest/header（若找到），记录每幅 day/band/method/Stokes/unit/WCS/beam/processing；只做 metadata 验证，不读大数组。若只能拿 UVFITS，则明确这是重新成像任务，非仅重绘 FITS。
3. 输入 provenance 无误后，使用论文最终核规则匹配有效分辨率；通过 I 相关平移并保存偏移；对所有 Stokes 同操作，禁止以 EVPA 好看为准调位置/角度。
4. 先输出未掩膜 I/Q/U/P 与 flux、unit/beam 自检，再按原 caption 做显示 cuts；mask AND/OR 未决处显式标记双版本，不把研究历史 mask 当官方。
5. 单 pipeline ensemble 的 Stokes mean 与 EVPA direction mean 各自命名；绝对 P tick 统一比例，fraction 色条固定，不按每日 peak 缩放。
6. 筛选表分开报告：真实 closure χ² gate、synthetic recovery score、D-term score、是否只是 EVPA-centrality selection；没有相应证据就不叫 official Top Set。样本数不作验收依据。
7. 日日差同时给 mask/beam/alignment/calibration sensitivity；优先 April 5/11，April 10 加弱约束标记；不把肉眼差异或单 ensemble std 写成 intrinsic variability 显著性。
8. 完成上述图与定义核验后，才能另行授权完整 synthetic/visibility scientific rerun；本次没有做该跑，也没有修改项目处理选择。
