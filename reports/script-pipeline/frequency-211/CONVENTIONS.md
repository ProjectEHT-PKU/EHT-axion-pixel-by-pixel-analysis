# Kerr massive-scalar 准束缚态实频率与 n 约定：只读核查报告

## 0. 结论先行

在自然单位 \(G=c=\hbar=1\)、\(\alpha=M\mu\) 下，massive scalar Kerr 准束缚态的领先实频率是

\[
\boxed{\frac{\omega_R}{\mu}=1-\frac{\alpha^2}{2N^2}+O(\alpha^4),\quad
N=n_r+l+1,\quad n_r=0,1,2,\ldots.}
\]

因此，无量纲角频率为

\[
\boxed{M\omega_R=\alpha\left(1-\frac{\alpha^2}{2N^2}+O(\alpha^4)\right).}
\]

- 本地文章明确采用 **principal n**，并明确全文研究 **211**；该态是 \(N=2,l=m=1,n_r=0\)，领先分母应为 \(2^2=(0+1+1)^2\)。
- 原 Notebook 的 `(2+1+1)**2` 等于 \(4^2\)。它**可以**来自合法公式中的 **radial overtone \(n_r=2\)、\(l=1\)、加常数 1**；此时若 \(m=1\)，对应 principal 标记 **411**，不是文章的 211。
- 不能把 \(N=n_r+l+1\) 中的 \(n_r\) 直接替换为文章的 principal \(n=2\)：那会重复加上 \(l+1\)。
- \(m\) **不进入领先氢样能级**。`2+1+1` 本身不足以证明作者在做 `n+l+m`，也不足以证明其实际意图是 \(n_r=2\)。两种误读都可能，但现有证据只查明了表达式、文稿约定及两者不一致，未查明作者心理过程。
- 径向仓库频率表低 \(\alpha\) 数值更接近 \(N=2\) 而非 \(N=4\)，但表有自旋插值且缺少 overtone 元数据；**数值接近不是完整本征身份、精度或边界条件的证明**。
- \(\alpha\simeq0.44\) 不能把 leading 小 \(\alpha\) 公式当精确 Kerr 频率。文献有明确的高阶结构及数值比较，详见第 6 节。

## 1. 范围与取证方法

只读项目；未修改任何项目文件或配置，未执行 Notebook 单元、径向生成器、科学求解器、编译或测试。仅下载一手论文至本目录，使用 `pdftotext`/页面图片读论文，以及做所要求的小表算术与现有 CSV 两点线性插值。没有重算本征值、径向函数或科学结果。

原 Notebook **通过 `git show HEAD:result.ipynb`** 读取，未使用工作树版本代替：

- `/home/canoziia/MyProject/PHYS/EHT_ZC` HEAD：`74ebae19eae6bb3f443f62b96abb4654d248b67a`。
- `/home/canoziia/MyProject/PHYS/radial-function` HEAD：`80ca3bd62b45cda49155d3e26c8864a01cdf59b9`；README、源码和表读取的是当前工作树。
- 文章读取 `/home/canoziia/MyProject/PHYS/EHTaxion_imaging_pipeline_constraint/main.tex` 当前工作树。

先联网读取官方 arXiv 摘要页/原始 PDF，再跟随其参考文献读取 Detweiler 的 APS 原始论文，以及 Baumann 2019 引用的 Baumann–Chia–Porto 论文。所有 curl 均使用 `--connect-timeout 10 --max-time 45`，无批量文献抓取。

官方重新下载的 Dolan PDF 与已有缓存 `/tmp/eht-constant-sources/0705.2880.pdf` SHA256 相同：`5589ca3bcc21a2675d18759aba9f1bac596d5557957823de77287ed4f09175ad`。故正文检索可可靠复用该缓存的 `.txt`。

## 2. 一手论文证据：三套 n 不能混用

以下页码优先给论文印刷页；Baumann PDF 另列文件页码以避免封面/目录偏移。只摘录必要短句。

### 2.1 Dolan 2007：n̄ 才是 principal，n 是径向整数

Sam R. Dolan, *Instability of the massive Klein-Gordon field on the Kerr spacetime*, Phys. Rev. D **76**, 084001 (2007)。

- 官方记录：<https://arxiv.org/abs/0705.2880v2>
- 原始 PDF：<https://arxiv.org/pdf/0705.2880v2>
- DOI：<https://doi.org/10.1103/PhysRevD.76.084001>
- **p.10，Eq. (32)**：氢样谱以 \(\bar n\) 为分母，紧接原句是 “**n̄ = n + l + 1 is the principal quantum number**”。
- 同段说明 scalar 也有该谱，并说最低阶不依赖黑洞旋转或场的自旋。
- **Appendix A，p.23，Eqs. (A8)–(A9)**：给出非相对论 Schrödinger 方程，\(E_{NR}=\omega-\mu\)，引力 Coulomb 势为 \(-M\mu/r\)，并回指能级 Eq. (32)。

对应关系：\(\bar n=N\)，其 Eq. (32) 中的 \(n=n_r=N-l-1\)。不要因为同一个字母 n 出现在不同文章中，就把它等同本地文稿的 principal n。

也不要将 Dolan Eq. (33) 的级数求和指标 n，或 Appendix B 的 **QNM** overtone n，误作正在讨论的准束缚态标签；这些是不同上下文。源码的递推求和 n 同样不是选定本征态的径向 overtone。

### 2.2 Detweiler 1980：n 非负，ν 在氢样极限是 principal

Steven Detweiler, *Klein-Gordon equation and rotating black holes*, Phys. Rev. D **22**, 2323–2326 (1980)。此文早于 arXiv，直接读取出版商原文。

- DOI：<https://doi.org/10.1103/PhysRevD.22.2323>
- APS PDF：<https://journals.aps.org/prd/pdf/10.1103/PhysRevD.22.2323>
- **p.2324，Eqs. (9)–(10)**：\(k^2=\mu^2-\omega^2\)、\(\nu=M\mu^2/k\)。
- **p.2324，Eq. (14) 后**：“**ν = l + 1 + n where n is non-negative**”；随后明确 \(\nu\) 是 principal quantum number。
- **Eq. (15)**：\(\nu-l-1=n+\delta\nu\)，说明实际 Kerr 准束缚态的 \(\nu\) 可有小复修正；不是对所有 \(\alpha\) 都严格等于整数。
- **p.2325，Eq. (26)**：
  \[
  \mu^2-\sigma^2=\mu^2\left(\frac{\mu M}{l+1+n}\right)^2,
  \quad\omega=\sigma+i\gamma.
  \]
  展开实部即得到领先 \(\sigma/\mu\simeq1-\alpha^2/[2(l+1+n)^2]\)。这不是声称 Eq. (26) 的平方根包含可靠的全部高阶 Kerr 修正。
- **p.2325，§III，Eq. (29) 前**：最快增长模写作 “**l = 1, m = 1, and ν = 2**”，并称为氢原子 **2p** 的对应态。因此其径向整数 n 是 **0**，不是 2。

APS 文本抽取存在 ν/v、σ 等 OCR 混淆；已用原 PDF 的 p.2324、p.2325 页面图逐项核对。保存在 `Detweiler-page2324.png`、`Detweiler-page2325.png`。

### 2.3 Baumann et al. 2019：n 直接是 principal

Daniel Baumann, Horng Sheng Chia, John Stout, Lotte ter Haar, *The Spectra of Gravitational Atoms*, JCAP **12** (2019) 006。

- 官方记录：<https://arxiv.org/abs/1908.10370v1>
- PDF：<https://arxiv.org/pdf/1908.10370v1>
- DOI：<https://doi.org/10.1088/1475-7516/2019/12/006>
- **p.3（PDF p.5），Notation and conventions**：\(|nlm\rangle\) 的 n 明确定义为 principal，并说 “**Following the convention in atomic physics, we have n ≥ l + 1**”。
- **p.8（PDF p.10），Eq. (2.18)**：径向函数含 \(L^{(2l+1)}_{n-l-1}\)。其阶数 \(n-l-1\) 就是弱耦合氢样极限的径向整数 \(n_r\)，而非 principal n 本身。
- **p.11（PDF p.13），Eq. (2.28)**：scalar 实能级第一项束缚修正为 \(-\mu\alpha^2/(2n^2)\)，另有 \(\alpha^4\)、\(\tilde a m\alpha^5\) 修正。这里的 n 已经是 principal，不能再加 \(l+1\)。
- **p.8，Fig.3**：211 在 n=2 层，311 在 n=3 层，直观展示该标记体系。

此文引用的前作也已跟读：Baumann, Chia, Porto, *Probing Ultralight Bosons with Binary Black Holes*, Phys. Rev. D **99**, 044001 (2019)，<https://arxiv.org/abs/1804.03208v2>，<https://arxiv.org/pdf/1804.03208v2>，DOI <https://doi.org/10.1103/PhysRevD.99.044001>。

- **p.6（PDF p.8），Eqs. (2.7)–(2.8)**：再次明确定义 principal n，给出 \(\omega\simeq\mu[1-\alpha^2/(2n^2)]\)。
- **p.8（PDF p.10），Eq. (2.14) 前的段落**：\(l=m=1\) 有一整塔 \(|n_g11\rangle\)、\(n_g\ge3\) overtones，而最快增长为 211。这直接说明 **l1m1 不唯一确定 overtone**。
- **p.13（PDF p.15），§3.2**：将 “fundamental mode (n = 2)” 与 “overtones (n ≥ 3)” 分开；这里仍是 principal 标签，不是把 overtone 从 2 开始计数。

## 3. n=0 ground 与 n=2 principal：完整对照

统一记 \(N\) 为 principal，\(n_r\) 为从零起的 radial overtone：

| 物理态/系列 | \(n_r\) | \(l\) | \(m\) | \(N=n_r+l+1\) | 本文/Baumann 标记 |
|---|---:|---:|---:|---:|---|
| 全谱氢样基态 | 0 | 0 | 0 | 1 | 100 |
| 固定 l=m=1 的径向基模；主要超辐射态 | 0 | 1 | 1 | 2 | 211 |
| l=m=1 第一径向激发 | 1 | 1 | 1 | 3 | 311 |
| l=m=1 第二径向激发 | 2 | 1 | 1 | 4 | 411 |

- **“n=0 ground”** 只有在 n 指 radial overtone 时才成立；固定 l 后，该系列的最低径向态为 \(n_r=0\)。全谱最低氢样态是 \(n_r=l=m=0\)、\(N=1\)。211 是 l=m=1 系列的径向基模，并非全谱最低能级；100 的 m=0 不满足正频率 scalar 超辐射条件。
- **“n=2”** 在本地文章和 Baumann 中指 principal \(N=2\)。当 l=1 时必有 \(n_r=0\)。principal n 不存在 n=0 的合法氢样态。
- Detweiler 的 \(n=2\) 或 Dolan Eq. (32) 的 \(n=2\)，若 l=1，则是 \(N=4\)，不能与 principal 211 中的 2 等同。
- 一些其他作者可能采用不同起点的 overtone 计数；本报告的 \(n_r\) 明确从 0 开始。跨论文必须先读定义。

## 4. 本地文稿与原 Notebook：查明事实和未查明意图

### 4.1 文稿明确约定，不是推测

`/home/canoziia/MyProject/PHYS/EHTaxion_imaging_pipeline_constraint/main.tex`：

- **L151–152**：n、l、m 分别为 principal、orbital、magnetic，类比氢原子。
- **L162**：再次称 “principal quantum number n”。
- **L179**：“**(n,l,m) = (2,1,1), throughout this work**”。
- 后续 **L362–363**、**L389–391** 继续使用 \(R_{211}\)、\(\omega_{211}^r\)，并用于 EVPA 时间相位与时间平均。
- **L181–185**：给出自旋依赖的超辐射范围，benchmark 自旋为 0.9375，fiducial \(\alpha=0.3\)。

所以不能把 L159 的 “Following Dolan” 解读成文稿也自动沿用 Dolan 的 radial n；文稿自己的显式定义已采用 principal n。

### 4.2 Notebook HEAD 中的事实

`git -C /home/canoziia/MyProject/PHYS/EHT_ZC show HEAD:result.ipynb`，JSON 的 **cell index 25（从零计）**：

```python
def _omega(alpha):
    return alpha * (1 - alpha**2 / 2 / (2 + 1 + 1) ** 2)
```

这只是读取 Git blob/解析 JSON，没有执行单元。该 cell 没有说明三个数字各指什么。其他调用单元把 `_omega(alpha)` 放入 EVPA 的时间差、时间平均正弦因子及相位中；其前面的 `alpha` 表明自然物理解读是 **Mω**，不是直接的 \(\omega_R/\mu\)。此处不额外认证 `0.112411` 的单位换算。

三个层次必须区分：

1. **已查明**：函数实现等价于 \(N=4\) 的 leading \(M\omega_R\)。
2. **可成立的公式解释**：三个数为 \(n_r=2,l=1,+1\)，而非 n、l、m 三个量；若 m=1，则其 leading 能级对应 411。
3. **未查明**：作者究竟有意选择径向第二激发，还是混用了 principal 与 radial n，或把 n+l+m 相加。代码字面数字和现有注释不足以判定。

如果坚持文稿全文 211 的定义，这个 `N=4` leading 表达式就不一致。仅由频率的 leading 数值也无法唯一认定 411，因为同一 N 的多个 l,m 在该阶简并；“411”还依赖 l=m=1 的额外前提。

## 5. 径向仓库：单位、第一列、数值小表及证据边界

路径：`/home/canoziia/MyProject/PHYS/radial-function/`。

### 5.1 get_omega 实际读取什么

- `source/a=0.9375_l1m1_real.txt` 第一行 `304` 是声明的行数，不是一个量子数。随后的每行第一列是质量参数 \(\alpha=M\mu\)，第二列按实现解释为 \(q=\omega_R/\mu\)。例如第一条数据为 `0.0749,0.9991870006389632`。
- `Axion_Radial.c:L149–154` 读行数和两列；**L174** 将第二列的质量方向线性插值结果乘 `mu`，得到 `omega_real`。**L100** 对虚部做相同的乘法，**L176** 合成复频率。
- `Function_R` **L234** 设 `M=1`，故代码变量 `mu` 数值上就是 \(\alpha\)，返回的频率数值相当于 \(M\omega\)。不能直接把表第二列当 \(M\omega\)，也不能把 `get_omega` 返回值再次误当 \(\omega/\mu\)。
- `Axion_Radial.h:L14–17` 的 `FLAG_L_M=11`、`AZIMUTHAL_M=1`、`LAMBDA_L_M=2` 只涉及角向选择及近似，**没有 principal/radial overtone 参数**。
- `Axion_Radial.c:L223` 的 n 是 1000 项递推级数指标，不能把循环 n=0 当作“已选择径向基态”的证据。这里读取既有频率，未见独立重新搜索完整本征谱的过程。

### 5.2 原 README 对精度与来源的限制

`README.txt:L28–32` 明说：

- a=0.95 的实部是 a=0.99 与 a=0.9 的插值；
- a=0.9375 的实部及虚部是 a=0.95 与 a=0.9 的插值；原句为 “**real and image part ... is an interpolation**”；
- l=2,m=2 还沿用 l=1,m=1 的实部数据。源码 **L133–145** 与这一说明一致。

`README.md` 也明确 `source/` 不是新求解的精确本征值。README 提及 Dolan，但没有提供足够的原始数据点、提取流程或 overtone 元数据以完整追踪每张表的本征分支。

### 5.3 指定小表

表列出 \(q_N=1-\alpha^2/(2N^2)\)。\(q_{tab}\) **只是在既有 a=0.9375 表上按 get_omega 的两点线性规则取值**，不运行 get_omega 或科学程序。它是来源不确定的插值参考值，不是真值。误差先定义为

\[
\epsilon_N^{freq}=\frac{q_N-q_{tab}}{q_{tab}}
=\frac{(M\omega)_N-(M\omega)_{tab}}{(M\omega)_{tab}}.
\]

下面百分数均为这一带符号的 **频率相对差**，不是束缚能相对误差。

| α | q，N=2（211 leading） | q，N=4（原函数/411 leading） | q，a=.9375 表 | Mω，表 | ε₂ freq | ε₄ freq |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.998750000000 | 0.999687500000 | 0.998547049499 | 0.099854704950 | +0.020325% | +0.114211% |
| 0.2 | 0.995000000000 | 0.998750000000 | 0.994854808317 | 0.198970961663 | +0.014594% | +0.391534% |
| 0.3 | 0.988750000000 | 0.997187500000 | 0.987730244214 | 0.296319073264 | +0.103242% | +0.957474% |

复核用的原表插值括号（文件行号含第一行的计数）：

- α=.1：**L7–8**，x=.0999、.1049；
- α=.2：**L27–28**，x=.1999、.2049；
- α=.3：**L47–48**，x=.2999、.3049。

需要 leading Mω 时，将 q 乘 α：N=2 依次为 **0.099875、0.199000、0.296625**；N=4 依次为 **0.09996875、0.19975000、0.29915625**。

另一种比较是原函数相对“211 leading”的差：\((q_4-q_2)/q_2\)，在上述三点依次为 **+0.093867%、+0.376884%、+0.853350%**。这与相对插值表的误差不是同一个分母。

若比较 **束缚能**，应取 \(B/\mu=1-q\)，并明确

\[
\epsilon_N^{bind}=\frac{(1-q_N)-(1-q_{tab})}{1-q_{tab}}.
\]

| α | ε₂ bind（相对表） | ε₄ bind（相对表） |
|---:|---:|---:|
| .1 | −13.968% | −78.492% |
| .2 | −2.822% | −75.705% |
| .3 | −8.311% | −77.078% |

仅比较两条 leading 公式时，N=4 的束缚能恰为 N=2 的 **1/4**，相对 N=2 是 **−75%**；绝不能把它称作频率误差 75%。

**解读限度**：这些低 α 点明显更支持 N=2 的 leading 分支，而不是 N=4；但表的自旋插值、潜在提取误差与有限 α 高阶效应混在一起，尤其 .1 点的束缚能相对差并不小。不能单凭这些数字给整张表盖上“精确 211 本征解”印章。完整鉴定还需原始频率数据 provenance、弱耦合分支连续性、径向阶数/节点对应关系以及两端边界条件和收敛/本征残差的验证；此次未进行这些科学检验。

## 6. α≈.44：leading 适用范围与 spin corrections

### 6.1 一手来源真正支持的精度声明

- Detweiler **p.2324，Eq. (8) 前**要求 \(\omega M\ll1\)、\(\mu M\ll1\)。
- Dolan **p.10，Eq. (32) 附近**给出非相对论/低耦合氢样近似；**pp.16–18，Figs.4–5** 给出随自旋变化的数值实频率，而不是把 leading 谱当全耦合精确公式。
- Baumann et al. 2019 **p.12（PDF p.14），Fig.6** 比较 scalar 数值谱和 **含高阶项的 Eq. (2.28)**，图注明确称其为 “**an excellent approximation for α ≲ 0.2**”。该图自旋为 **ã=.5**。

这里的 .2 不是普适硬截止，也不是单独 leading 项的保证；原图比较的是含 fine/hyperfine corrections 的式子，且它的 ratio 面板针对 **ν**。不能把图中 ν 的百分差未经换算就当作 ω 的百分差，更不能外推成 ã=.99、α=.44 的精度保证。

Baumann–Chia–Porto 2018/2019 **p.8** 提到 α<.5 时某解析结果与数值有合理一致，但那段指的是 **增长率 Eq. (2.12)**，不是实频率 leading Eq. (2.8)。超辐射允许的 α 上限同样不是氢样 leading 实频率的精度上限。

### 6.2 已可靠核对的高阶结构

Baumann et al. 2019 **Eqs. (2.28)–(2.29)，pp.11–12**：

\[
\frac{\omega_R}{\mu}
=1-\frac{\alpha^2}{2N^2}
-\frac{\alpha^4}{8N^4}
+\frac{f_{Nl}}{N^3}\alpha^4
+\frac{h_l}{N^3}\tilde a\,m\alpha^5+\cdots,
\]
\[
f_{Nl}=-\frac6{2l+1}+\frac2N,
\qquad h_l=\frac{16}{2l(2l+1)(2l+2)},\qquad\tilde a=a/M.
\]

以上是直接转写经核对的一手结果，不自行推导高阶谱。它还与前作 **Eq. (2.10)，p.6** 的写法相符。h 的该表达式用于 l>0；l=0 时 m=0，不应直接代入一个形式上发散的 h。原文 §3.1.3/脚注亦说明这一点。

由此可可靠得出：

- leading 束缚修正是 \(O(\alpha^2)\) in **ω/μ**，在同 N 内不含 l、m、自旋；l 只通过 N 的允许值约束出现。
- 在 **ω/μ** 中，\(O(\alpha^4)\) 给出相对论/精细结构并解除 l 简并；第一项显式旋转相关的 m 劈裂是 \(O(\tilde a m\alpha^5)\)。这不表示更高阶没有其它自旋结构。
- 在 **Mω** 中整体多乘一个 α，故上述阶次分别变为 leading 质量项 α、束缚项 α³、精细项 α⁵、显式旋转/m 项 α⁶。说“spin correction 是 α⁵”时必须说明指 ω/μ 的展开。
- m 可以影响高阶实频率、吸收/增长率及超辐射条件 \(\omega_R<m\Omega_H\)，不能因此把 m 塞进 leading 分母。

仅作上述可靠公式的代值展示（**不是强耦合精度认证**）：211 有

\[
\omega_R/\mu=1-\alpha^2/8-(17/128)\alpha^4
+(\tilde a m/12)\alpha^5+\cdots.
\]

α=.44 时，leading 211 给 **q₂=.9758、Mω₂=.429352**；原 N=4 式给 **q₄=.99395、Mω₄=.437338**。同点上述 α⁴ 项为 **−.00497794**，取 ã=.99、m=1 的 α⁵ 自旋项为 **+.00136056**。这些高阶项已不能在需要精密相位时无条件忽略；更不能把截断和精确解的误差认定为“恰好等于下一项”。未求可靠数值本征值，故不报告 α=.44 的真实相对误差。

### 6.3 自旋与 α=.44 的额外限制

本地 README 明示 .9375 使用的虚部表上限为 .316，get_omega 检查实、虚表的共同支持范围；**不能用这张 .9375 表对 .44 无提示外推**。a=.99 的 l1m1 支持上限说明为 .446，也不能把它当作近似精度保证。

\(M\Omega_H=\tilde a/[2(1+\sqrt{1-\tilde a^2})]\)：ã=.9375 时约 **.347741**，ã=.99 时约 **.433804**。α=.44 的 leading 211 Mω≈.429 因而靠近高自旋 .99 的超辐射阈值，却不属于 .9375 的 211 增长范围。此阈值比较只是提示敏感性，不替代精确谱判定；“束缚态存在”与“该态仍在超辐射增长”也不是同一命题。

## 7. ω、μ、Mω 与建议版本（仅建议，不实现）

在本报告自然单位中 μ 为场质量，也可作为自由场静止角频率；\(\omega_R\) 是带束缚修正的模式角频率，一般小于 μ，二者不能在精确相位中无条件互换。若恢复 SI 且轴子质量为 \(m_a\)：

\[
\mu_{ang}=m_ac^2/\hbar,\qquad t_g=GM_{BH}/c^3,
\qquad\alpha=t_g\mu_{ang},
\]
\[
\omega_R=\mu_{ang}\left(1-\frac{\alpha^2}{2N^2}+\cdots\right),
\qquad \widehat\omega=t_g\omega_R=\alpha\left(1-\frac{\alpha^2}{2N^2}+\cdots\right).
\]

相位是 \(\omega_R t_{SI}=\widehat\omega\,(t_{SI}/t_g)\)。若使用每秒周期数，则另有 \(f=\omega_R/(2\pi)\)，不要混用角频率和周期频率。原函数的返回形式对应 \(\widehat\omega=M\omega_R\)。

建议未来审定后采用以下**显式数学接口之一**，不混用两套 n：

1. **principal 版**：\(\widehat\omega_{leading}(\alpha;N_{principal})=\alpha[1-\alpha^2/(2N_{principal}^2)]\)，并要求 \(N\ge l+1\)。
2. **radial 版**：\(N=n_r+l+1\)，再代入上式，明确 \(n_r\ge0\) 从零起。
3. 本地文章 211 应明确 **Nprincipal=2**，或 **nr=0,l=1,m=1**；leading 表达式是 **α(1−α²/8)**。m 不作为 leading 能级分母输入。
4. 函数名/文档必须标明输出 `Mω_R` 还是 `ω_R/μ`，以及 `leading_small_alpha` 近似级别。α≈.44 若需要精密相位，应先取得可追踪的、同一 (N,l,m)、同一自旋的 Kerr 数值谱或经过误差验证的近似。
5. 频率的时间振荡、径向表使用的频率、文稿态标签应一致审计；只换分母或只换径向表，均不等于整个 pipeline 已修正。

**最终判定**：现有证据足以指出“原函数 N=4 leading 与本地文稿 principal 211 不一致”，也足以给出正确的 211 leading 版本；不足以断定原作者一定把 n+l+m 相加，或足以认证径向表为精确 211 本征解。此次只提出建议，没有实现任何变更。
