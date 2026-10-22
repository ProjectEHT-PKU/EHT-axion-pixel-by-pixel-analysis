# Paper VII → installed ehtim 1.3.2：focused objective audit

## 主判 / 给父任务的阻塞项
- **不能认证 exact objective mapping。没有证据证明 VII 的 β 对应 unnormalized，也没有证据证明对应当前 `norm_reg=True`。** 表 6 数值相同不够；切换 norm 而保留 β 是换目标，不是等价实现。
- 确证额外项：旧 `pol_cal.leakage_cal` 始终加零中心 D 二次 penalty；`.1` 不是论文认证的 fiducial，不是硬上限。
- 确证实现差别：VII p.23 描述 Jones 正向污染模型；driver 用 Jones 逆校正数据、随候选 D 传播对角误差，不能仅凭“full Jones”认定 likelihood/weights 相同。
- VII p.22 的 HW 描述提及 `.75`；当前 HW 没有 `.75` 参数/阈值，只在 m→1 处斜率发散。软偏好同类不等于公式已认证。
- **working-current-model 冻结建议**：保留 `norm_reg=True`、现行 fsys helper、旧 `leakage_cal(inverse=False)`（父 runner 显式指定）为对照；承认 paper exact norm 未知，不盲改。旧重建/科学结果全部冻结；本审计不重新成像、不优化、不修改 owner code。仅新写本报告与 `/tmp/eht-objective-audit/`。

## 证据与范围
- 论文仅直接读 rendered images：`/tmp/eht-official-method/pageimages/VII-{22,23,31}.png`；未用 PDF text/OCR。
- p.22 C4 是数据项与正则加权和；C5 按 C4 和上下文理解为 **两数据项相加**，不将排版异常解释为相减/相等。
- p.31：αP=1、固定 αm=1、βHW=100、βTV=1、niter=50、fsys=.005；加噪幅度 `fsys*|Ĩ|` 与 thermal quadrature 合成，但未指定每个 Stokes/cross-hand 分量约定。
- 当前库根 `L=/home/canoziia/MyProject/PHYS/eht-script-pipeline/.pixi/envs/eht-rebuild/lib/python3.12/site-packages/ehtim`。
- `D=/home/canoziia/MyProject/PHYS/eht-script-pipeline/scripts/rebuild_eht_pol_warmstart.py`；helper 为同目录 `rebuild_eht_fiducial.py:138–163`；配置 `configs/eht-fiducial.yaml`。
- 有界联网：仅固定官方 commit `c8c7a89ea96da0c3f0217e95c41ce8d3e7ca8765` 的三个 pol 源文件和 `pyproject.toml`，每次≤20s、≤200KB；保存在 `/tmp/eht-objective-audit/`。
- 该 commit 的 `pol_imager_utils.py`、旧 `pol_cal.py` 与安装态逐字相同；这是当前源码证据，**不是 2021 原执行脚本证据**。未作 1.2.7 大 diff。
- SHA256：installed `imager.py` `ce505c2bcacfacccda9fffe2c986e5c7bcef740fdb58c8b998e4b7dca5b62b87`；pol utils `4a44b878d968ff12059ab63c61cebaaa2027e283360bf0934e408dcec3807948`；旧 pol_cal `52b39c1fffa88d568c6bb88d7ca7a6317f31197f9d5f010f3785c146bf4456f2`。

## 公式：当前 P imaging 实际最小化什么
令像素值为积分 Jy/pixel，`Pj=Qj+iUj=Ij*mj*exp(iφj)`，φ=2EVPA；`F=.6 Jy`，Δ=1.875 μas，b=`obs.res()`=1/max|uv|，不是恢复核 19 μas 或 restart 核 20 μas。
令 A 为带 image pulse 的 direct Fourier matrix，dP 为观测 P，dm=dP/dI。
```text
CP = [Σk |dP,k − (AP)k|² / σP,k²] / (2 NP)
Cm = [Σk |dm,k − (AP)k/(AI)k|² / σm,k²] / (2 Nm)
h(m) = ((1+m)/2) ln((1+m)/2) + ((1−m)/2) ln((1−m)/2)
Hraw = −Σj Ij h(mj)                              [positive entropy]
Traw = Σj sqrt(|Pj+ŷ−Pj|² + |Pj+x̂−Pj|²)         [positive TV]
Rhw = −Hraw/F
Rptv = Traw / (F Δ/b)
Jcurrent = (CP−1) + (Cm−1) + 100 Rhw + Rptv
         = CP + Cm − 2 − 100 Hraw/F + Traw*b/(F Δ).
```
- `pol_imager_utils.py:578–635,1413–1481`：分母明确 **2×返回的复数样本数 N**，不减拟合参数、不是 residual Ndof；P/m 均 `unpack(...,conj=True)`。
- 原观测 n 行时返回 N=2n：负 uv 的 P 是另一 cross-hand 的共轭，不可简单称为重复同一 P。Stokes 表示的 sigma 对正负 uv 相同。
- `imager.py:1572–1613`：默认 `chisq_transform=False`，数据项各减 1；总 offset 为 −2，与最优点/梯度无关。不能把负 J 当负 χ²。
- `pol_imager_utils.py:449–488,976–1032`：wrapper 返回 **负 entropy** 与 **正 TV**；TV 是复数梯度的联合 Euclidean norm，不是分别求 Q/U TV 后相加。
- TV 在右/下方向取 forward differences，图外补零，计入末行/末列边界；无额外 1/2、无 epsilon。Δ、b 只在 normalization 中出现。
- `imager.py:196–201`：b 自动取输入观测 nominal resolution；换 uv 子集可能改变正则相对权重。
- 相对于同一 raw H/T，βHW,raw=100/F=166.6667，βTV,raw=b/(FΔ)=22.2319；二者不是同一共同缩放，不能整体 rescale J 消除。
- **这些是当前 norm-on/off 的确定换算，不是已证实的 VII→current β 换算。** VII 指定页未给 reduced-χ² 分母、flux/beam/grid normalization 或离散边界，仍需原脚本/公式定义。

## prior、mcv、`.75`
- `polregularizer` 与 gradient 明注 `priorarr is currently not used`；HW 依赖 **当前 m 与固定 I**，并非只依赖 I；不依赖 prior 的 Q/U。
- driver 每轮 `Imager(current, init, init, .6, ...)` 不产生相对上一轮 P 的 KL/距离先验，不“重新锚定 HW entropy”。
- prior 的 I 决定 `I>clipfloor(0)` mask、几何/pulse；fixed I 全正、flux 显式 .6、`norm_init=False`，这些逐轮不变；init 的 Q/U 则确实决定起点。
- `pol='P'` 的 solve slots 是 `[0,1,1,0]`：I 固定、V固定为0、优化 m 与 φ。`mcv` 在 V=0 时 `m=.5+atan(κ/.5)/π`，确保 0<m<1；φ=2EVPA 只是重参数化。
- `mcv` 无 prior 或 Jacobian log-density penalty；不能把 variable transform 当额外 Bayesian prior。
- 当前 `d(−Hraw/F)/dmj=(Ij/F)*atanh(mj)`：从0连续增罚；没有 `.75` 特征点。VII 的软偏好描述与这个族兼容，但不足以证明作者用了同一缩放/改写；不建议擅加 .75 clip。
- driver 的 20 μas QU blur/projection 是 restart 起点操作，不是目标中的新正则；指定轮次 2,12,22,32,42 未由论文认证。

## fsys=.005、factor 2 与实际权重
记 `s=.005*|dI|`。helper 每次从 fresh base 取数据，向 `σQ²,σU²` 各加 `s²/2`，不向 `σI²,σV²` 加噪。
```text
σP² = σQ²+σU²                    → ΔσP²=s²
Stokes→circ: σRL²=σLR²=σQ²+σU²   → ΔσRL²=ΔσLR²=s²
σm² = [σQ²+σU²+σI²|dm|²]/|dI|² → Δσm²=fsys²   (before Jones)
```
- 依据 `obsdata.py:256–326,878–894`、`obs_helpers.py:686–700`；`psigma` 是代码所用复数残差单分量误差约定，不应另乘√2再保留同一 χ² 公式。
- **这里没有实际 RL 方差×2 的 bug。** 若另一个实现向 Q/U 各加 s²，则其 P/RL 方差增加 **2s²**；等价到本 helper 要 fsys×√2。若把 s 解释为复残差总 RMS，分量又需另一个约定；论文指定页未选定。
- circ→Stokes→circ 把不等 RL/LR thermal variances 对称化为 `(σRL²+σLR²)/2`，丢弃相关性；并非无损恢复原 cross-hand weights。
- P imager 的 polarimetric branch 不传 `systematic_noise`；helper 外部注入确有作用，不能再靠构造器参数重复注入。
- P 和 D 收到的是同一 fresh-base inflation **规则**，不是全过程相同权重：round1 P 用未 D 校数据，以后 P 用上一轮校正结果；D 每次从 base 重解绝对 D。

## D objective：旧 API 的确切形式
`L/calibrating/pol_cal.py:43–46,192–239`：默认 `pol_fit=['RL','LR']`, `const_fpol=False`, `inverse=False`。
```text
Vc(D) = J1(D)⁻¹ Vobs J2(D)⁻ᴴ
σc,ab²(D) = Σcd |(J1⁻¹)ac|² |(J2⁻¹)bd|² σcd²
Ch(D) = Σk |Vc,h,k(D)−Vimage,h,k|²/σc,h,k²(D) / (2 n)   h=RL,LR
JD(D) = [CRL(D)+CLR(D)]/2 + Σstations (|DR/.1|²+|DL/.1|²).
```
- 三待解站共6 complex/12 real D；每个 |D|=.1 贡献1 penalty，|D|=.01 仍贡献.01，全部6个=.1 时总 penalty=6；无截断、无 tolerance “死区”。
- prior 与 **mean reduced χ²** 相加，不随 n 缩小；若改为 raw residual sum，必须同时乘 **4n** 到 prior 才等价。这是 data/prior 相对权重，不是无害 offset。
- 每轮 D 起点均零，不是 previous D warm start；校正非线性完整，逆 Jones 含 `1/(1−DL*DR)`，不是 const_fpol 线性近似分支。
- `obs_simulate.py:1109–1189`：若任一 SEFD≤0，用输入 sigma；否则重算 SEFD thermal sigma，会覆盖外部 inflation。本次旧输入所有 sefdr=0，确认走输入 sigma 分支。
- 传播只保存各相关产品 marginal variance，不保留 full covariance；候选 D 改变权重，无 logdet 项。只用 RL/LR 的 inverse residual 一般不同于正向 RL/LR 固定权重拟合。
- VII p.23 写 χ² between gain-calibrated data and Jones-corrupted model、V=0；**没有写 D prior，不等于证实历史代码没有 prior**。当前 prior 是“超出论文形式的已确认项 / 历史是否存在未知”。
- 当前 commit 的 `pol_cal_new.py:199–253` 用包含二阶 D 的正向污染模型与固定 sigma，**更贴近 VII p.23 文字，适合作独立 controlled case，而非直接替换 baseline**；仍有 D prior，`rescale_leakage_tol=True` 除以12个实参数。本 driver 未调用它，不能换 API 就称历史官方等价。

## 无优化 objective 分解：只读冻结 snapshot
读取 `/home/canoziia/GRMHD_DATA/analysis/eht-official-PD-warmstart-apr11-v132-v1/round-01-raw-IQU.fits` 与原 I-gaincal UVFITS；未调用 minimize/make_image_P/leakage_cal。
- n=150，NP=Nm=300，F=.6，Δ=1.875 μas，b=25.0108782844 μas；HW divisor=.6，TV divisor=.0449804276046。
- P：CP=8.20900031111936，Cm=8.25174180148320，Rhw=−.675137046415785，Rptv=1.38879295900950。
- **J=8.2090003111+8.2517418015−2−67.5137046416+1.3887929590=−51.66416956996645**，复现旧 steps objective −51.664169569966525。
- 同 snapshot 与 saved D：CRL=5.50528574659，CLR=4.48130160958，mean=4.99329367808，prior=.639984881608；**JD=5.633278559692**，复现旧 D objective 5.633278559692。
- prior/data≈12.82%，不能忽略其对解的作用（数值比例本身不测量解偏移）；n=150 的等价 full residual sum 写法须将 prior 同乘600，当前相当于 raw-sum penalty≈383.991，不能保留.640后称同一目标。
- 此 D 后 RL sigma/输入 sigma 范围 [.9984058924,1.0001274590]，证明权重确随 D 改变；P helper 的 post-D pvis=4.99357934322 不等于 D mean，不能把两日志量强行互换。
- 父 probe `/tmp/eht-objective-audit/component-values.json` final 图：Hnorm=.6747202937574、Hraw=.4048321762545；Tnorm=.672062131404、Traw=.0302296420474。
- 故 **HW=1.6666667×raw、TV=22.2318918×raw，不是44.46×**；weighted reg cost norm-on −66.7999672443、off −40.4529879834。该文件无 χ²，不能据此给 final total J。
- 本次可重读脚本、数值、日志在 `/tmp/eht-objective-audit/{decompose.py,decomposition.json,decomposition.log}`；不声称检验新收敛/科学结论。

## 映射判定表与后续（由 parent 执行）
| 项目 | 判定 | 置信与行动 |
|---|---|---|
| P/m 两数据项，α=1/1；HW/PTV β=100/1 | 名称/数值一对一 | 高；不代表归一化一对一 |
| 固定 I、V=0、m/EVPA、m<1 | 一对一（V=0） | 高；当前 φ=2EVPA 不改目标 |
| χ²=Σ/(2N)，−1 offsets | 当前确证；VII normalization未知 | 高；禁止以 Ndof 替换，offset不影响极值 |
| norm_reg flux/beam/pixel | 当前明确缩放；VII对应未知 | 高；无证据主张应关 norm，也无证据认证应开 |
| HW `.75`、TV离散边界 | 公式细节未认证 | 高（当前），低（历史）；需原公式/脚本 |
| fsys=.005 | 幅度值一致、分量约定未知 | 不默认把×2方差选项认作修复 |
| Full Jones，D prior与weights | full nonlinear一致；objective非认证 | 当前高；论文未明述 prior、正向/逆向不同 |
| prior=init | 非逐轮 polarization prior | 高；固定 I 下只改起点，不改 HW/PTV函数 |
- 最小 parent controlled cases：冻结同一 I/obs/image，打印 CP、Cm、Hraw、Traw、F、Δ、b、N、D data/prior；不同 norm 仅做敏感性，不标“论文修正”。
- D 比较需分别控制 prior 有/无、正向/逆向、固定/传播权重；勿一次换 API+prior+noise 再归因。保留旧运行 manifest 与结果，不覆盖。
- exact certification blocker：需作者/官方 VII executable script/version、HW/PTV normalization、fsys 表示约定、D prior/likelihood 定义；当前 tag 不是历史替代证据。
