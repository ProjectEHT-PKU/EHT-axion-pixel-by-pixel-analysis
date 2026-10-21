# EHT 本地 producer / 候选集 / selection 来源审计

## 结论与边界

本报告只审计本地代码、metadata、sidecar 和小样本 FITS；不重复官方 paper 研究，不宣称重建官方 Top Set。未跑 reconstruction、EVPA/covariance/science 重算或大数据 hash，未修改原文件。机器总清单：`/tmp/eht-local-provenance/inventory.json`。父 driver 应使用标准 trusted `pickle.load` 或 astropy/ehtim；`/tmp` 中一次性 metadata-prefix/FITS 审计 reader **不是项目 runtime API**。

**已查明：**当前 raw 四天候选各4800且 identity 唯一，完整网格保留；五个底层 data-fit 数值已在 pkl 与 txt sidecar 中保存。当前1572排名不是这些 visibility χ² 的排名。可以立即用现有 raw 图像及 metadata 重做本地 selection 敏感性，无需重新成像。

**重要缺口：**确切 ehtim 库版本/自定义修改、χ² normalization 与数据绑定细节未找回；最终 FITS 是 restoration 后图像，而所存 χ² 是之前的 optimizer-stage 数值。缺少 `amp`、rho、偏振 closure 等新指标时，应从 existing FITS + UVFITS 做有限 data-fit 重评估，不可把旧指标改名冒充。

## 1. 本地来源与备份

代码根（下文 C）：
`/home/canoziia/GRMHD_DATA/source-backup/EHT/M87-image/eht`

数据根（下文 D）：
`/home/canoziia/GRMHD_DATA/EHT_DATA/EHT/M87-image/eht`

已读 `EHT_DATA/records/COMPLETE.json`：备份 source 为原 `PHYS/vmdk/home/yuxin/Work/Work/EHT`，scopes 为 M87-image、eht、eht-imaging、eht-imaging (copy)，72739文件、70623281765 bytes，`complete=true, sha256_verified=true`。本次采用既有验证记录，没有再次 hash。代码备份 `source-backup/records/summary.json` 记录132660个 selected files，包含 Git/source/scripts/notebooks/configuration/CSV 等。原 `PHYS/vmdk` 当前为空，未依赖它读取原机 install。

`D/poloutput/230_{95,96,100,101}`：每目录4800 FITS，每个149760 bytes，identity 全覆盖，无 extra/missing。`D/pol_sorted/poloutput/230_*` 也各4800，但**不是相同 payload 的排序副本**；16对同名小样本文件全部不相同，Q/U 全部不同，部分 I 也不同。另有13张历史试跑 FITS，不应并入四天 grid。

## 2. 实际 candidate grid 与 data-term 权重分离

代码 `C/eht-running_Intenisty.py:22–37` 和 `C/eht-running.py`，并由全部4800 metadata identities 与 FITS 文件名交叉确认：

| metadata identity key | 候选值 |
|---|---|
| `zbl` | 0.4, 0.5, 0.6, 0.7, 0.8 Jy |
| `prior_fwhm` | 40, 50 µas |
| `sys_noise` | 0, 0.01, 0.02, 0.05 |
| `MEM` (`simple`) | 10, 100 |
| `TV` (`tv`) | 0, 1, 10 |
| `T2V` (`tv2`; CLI `TSV`) | 0, 1, 10, 100 |
| `l1` | 0, 1, 10, 100, 1000 |

乘积 `5×2×4×2×3×4×5=4800`。频率标签230；四天95/96/100/101。当前 Intensity runner 文件只启用 day95，pol runner 只启用101；完整四天实际产物存在，不能把保存时的 loop snapshot 当作四天执行日志。

`C/eht-imaging_pipeline_pol_Intensity.py:86–104, 302–340, 351–372` 的固定 optimization weights：

- I round1：`amp:0.2, cphase:1, logcamp:1`。
- I round2：`vis:2, cphase:10, logcamp:10`。
- I round3/4：`vis:10, cphase:20, logcamp:20`；final I rounds 使用 `systematic_noise=0.01`，前两轮使用 station SEFD error budget（LM额外容差）。
- I regularization 为上表 grid，另有固定 `flux:1e4`，`norm_reg=True`。
- P 两轮：先 `pvis:1,m:1`，再 `pvis:5,m:5`；`hw:1,ptv:1` 固定，`transform='mcv'`。

**这些是本地 producer 的权重，grid 的 `l1=1000` 不是官方 Dterm=1000。**没有为本地4800候选找到扫描 Dterm=1000 的证据。官方 acceptance 与 weight 定义由父侧 paperagent 独立核对。

## 3. I / polarization / self-cal 与 χ² 保存阶段

### 可直接由 pipeline source 确认

输入为 `D/data/casa_data/April{05,06,10,11}/SR2_M87_2017_{095,096,100,101}_lo_casa_EVPA_rotation+dcal+ampscale+netcal+10s-avg+RLgain-amp.uvfits`。实际 lo/hi 共8份均存在，runner 不传 `-i2`，所以执行配置指向 lo，而非自动 lo+hi。

数据路径：scan coherent averaging → AA–AP median 估 total flux → 短基线 `<0.1e9` 的 Stokes vis/sigma 按 compact zbl 缩放 → station SNR ordering → reverse taper 5 µas → candidate `sys_noise`。保存 `obs_sc_init` 后进行 LM 短基线 Gaussian model selfcal（0.6 Jy、60 µas、`method='both'`），随后多轮 phase / LM amplitude / all-listed-site both selfcal。最后两轮 imaging 后，各有一次从 `obs_sc_init` 出发的 gain solve/apply。

I uses `amp` initially, then complex `vis` plus `cphase,logcamp`；polarization uses `pvis,m`。源注释把 `m` 注成 closure phases，但这里只能安全识别其为库的 polarimetric-ratio term，**不能称作偏振 closure χ²**；精确 m/rho 关系需要缺失库实现确认。source 没有独立 `rho` term 或输出。

### 保存时机（最容易混淆）

1. 最后 I round imaging/converge 结束后，脚本执行最后一次 selfcal，得到局部变量 `obs_sc`。
2. 随后 `chi2_I_term_dict = imgr.chisq_dict`（Intensity版第351行）。它是 **最后 I imager 的诊断**，不是显式以 restoration 后 FITS 重新调用 `Obsdata.chisq`。最后 selfcal 的新 `obs_sc` 没有在可见代码中显式重新绑定 `imgr`；是否存在库内部绑定，无法由缺失的库证明。
3. 两轮 `make_image_P` 后，`chi2_term_dict = imgr.chisq_dict`，将 `m,pvis` 及缓存的 I `cphase,logcamp,vis` 写 sidecar。
4. **之后**才 `im_out = im_out.blur_circ(5 µas)` 并 `save_fits`。

因此 stored χ² 不是 final-restored-FITS data-fit 的自动保证；I/P 本身也对应不同阶段。仅凭脚本无法核定 `chisq_dict` 是 reduced χ² 的哪种 normalization、closure covariance/minimal-set具体实现，以及 polarization 使用的 observation 精确状态。也不能把 optimization data weight 与 sidecar 中的单项值相乘后当既有排序。

未在 examined target tree 找到按 candidate 保存的 gains/caltables 或 selfcal UVFITS。不能用最终图像逆推出历史 gains。历史阶段精准复现可能需要原库/阶段 checkpoint，不只是 final FITS。

## 4. χ² 在哪里：实际 pkl / logs / FITS

### raw/current input（可供 driver 直接使用）

`D/result/230_{day}.pkl` 四份各630542192 bytes。此前项目已记录 top-level：`meta_data,I_data,Q_data,U_data,theta_data`；I/Q/U/theta 原轴为 `(64,64,N)`。本次只读 metadata 前缀，未复制大数组。

每份 `meta_data` **4800行，4800唯一 identity，全部数值有限**；keys 恰为：

```text
zbl, prior_fwhm, sys_noise, MEM, TV, T2V, l1,
m, pvis, cphase, logcamp, vis, chi2
```

`C/poloutput/meta/230_{day}_log_chi2.txt` 四份每份4800行，无 parse error/重复/缺格；去掉 pkl 派生 `chi2` 后，日志与 pkl metadata multiset 全相等。

`/home/canoziia/GRMHD_DATA/EHT_DATA/analysis-inputs/230_{day}.pkl` 实际解析到 `EHT_DATA/EHT/eht-imaging/result/230_{day}.pkl`，本次四份 metadata 与上述 raw **逐行完全相同**。没有对全部大数组做相等性或 hash 验证。

### 原 producer sorting

`C/eht-running_Intenisty.py`：

```python
chi2 = (vis + cphase + logcamp) / 3.0
```

四天实际 raw pkl 的4800行均逐项验证精确符合，且按 `chi2` 升序。**这里三项等权，不是 optimization 最后阶段的10:20:20；不是按数据点数加权。**源代码先 stable-sort `m`，再 `pvis`，最后 `chi2`；所以排序等价于 `(chi2,pvis,m,原日志行顺序)` 的 tie-breaking，不是三次加权组合。代码保留所有候选，没有在 producer 输出阶段截为1572或设 acceptance threshold。

`C/eht-running.py` 的另一旧分支定义 `chi2=(m+pvis)/2`。实际 `D/pol_sorted/result/230_{day}.pkl` 符合该定义，**没有** `vis,cphase,logcamp`：

| day | metadata行 | 唯一identity | 额外重复行 |
|---|---:|---:|---:|
| 95 | 4800 | 4800 | 0 |
| 96 | 4800 | 4800 | 0 |
| 100 | **6742** | 4800 | **1942** |
| 101 | 4800 | 4800 | 0 |

day100 sidecar 同样6742行，和该 pkl multiset 相符。不要无声截断/按row对齐/任意去重：不同重跑的日志指标未必属于当前同名 FITS。该885156976-byte pkl 本次未整份反序列化，只读约1.47MB metadata 前缀。`pol_sorted` 不是“把raw原图按另一指标排序”而已；同名图像样本的Q/U不同。

### FITS 与 sidecars

`fits-headers.json` 记录32个 FITS 的全HDU header（raw/pol_sorted各16，跨四天，抽取χ²首/中/末及固定grid样本）。每个有 I,Q,U,V 四个64×64 HDU；**第3号HDU是V，不是EVPA theta**。producer runner 把第3HDU塞入 `theta_data`，所以 driver 应从Q/U计算EVPA，不应信其变量名。无扩大坐标/flip/spin2审计。

这32份 sample 的全部 HDU header 未出现 χ²、rho、pvis、Dterm、selfcal、HISTORY/version 等匹配项；常规 WCS/FREQ/MJD/BUNIT/STOKES 信息存在。这里是**样本结论，不声称穷举19200 headers**。可见 writer 仅 `im_out.save_fits`，无将χ²显式加header逻辑。

相关 metadata txt 位于代码备份，不在 D/poloutput。未在 examined target tree 发现承担 fit table 的 CSV；已有 txt 和 pkl 足够，无需要求不存在的 CSV。默认保存 display PDF；`--imgsum` 默认False，代码中可选 summary 若运行会另做 selfcal/data comparison，不能拿默认display PDF当存储过独立fit的证据。

## 5. 底层指标统计与 threshold 诊断

以下 raw 指标各天均4800个finite。范围是 stored producer diagnostic，不等于最终FITS的新fit。

| day | vis范围 | cphase范围 | logcamp范围 | m范围 | pvis范围 | mean-I chi2范围 |
|---|---|---|---|---|---|---|
|95|0.43912–8.44748|0.10272–0.90595|0.23896–1.31871|0.25115–4.37101|0.27381–4.58875|0.26179–3.16506|
|96|1.03353–10.39120|1.02215–2.22283|0.48317–1.72369|0.31309–13.99834|0.36190–14.50373|0.84928–4.34114|
|100|0.97280–9.80341|0.35275–1.43221|0.82314–2.13551|0.31300–3.94199|0.35582–4.31245|0.71657–4.10307|
|101|0.76286–11.35637|0.19473–1.72363|0.33581–1.96323|0.25997–4.23829|0.29422–4.48180|0.43637–4.58780|

**以下全为本地 sensitivity conditions，不是官方 acceptance。**比较符为 `<=`；joint 指同一7-key identity 在四天都通过。顺序95/96/100/101。

| 条件 | ≤1.5 每日 | joint | ≤2 每日 | joint |
|---|---|---:|---|---:|
| vis |2875/833/1186/2006|608|3445/2890/3009/2992|2695|
| cphase |4800/1263/4800/4743|1263|4800/4664/4800/4800|4664|
| logcamp |4800/4753/3532/4750|3531|4800/4800/4787/4800|4787|
| I closures 各自 |4800/1263/3532/4724|1263|4800/4664/4787/4800|4654|
| I三项各自 |2875/428/1158/2006|**290**|3445/2890/3009/2992|**2695**|
| I三项均值 chi2 |4273/2778/3272/3638|**2736**|4559/4120/4228/4326|**4093**|
| m,pvis各自 |3600/3600/3600/3600|3600|3600/3600/3600/3600|3600|
| 五项各自 |2230/428/850/1547|290|2608/2124/2148/2228|1989|

≤1 时 day96 `vis,cphase` 均无通过项，所以I三项joint=0；但mean-I chi2 joint仍217。这再次说明均值条件不能替代各项条件。完整阈值1/1.5/2/3/5、分sys_noise数量、各selection交集、分位数在 `selection-diagnostics.json`。

## 6. 项目1572、producer旧排序与索引契约

当前 `eht-script-pipeline/eht_axion/pipeline/eht.py:25–26, 426–485`：先按7-key identity建立四天共同字典序轴，再用各天已配准 EVPA residual：

`score(identity) = sum_days sum_common_mask_pixels [circular_residual(theta, ensemble_mean) / ensemble_std]^2`

稳定升序取前1572。它是 **image-space EVPA ensemble 中心性**，既不是 visibility likelihood，也不是官方 Top Set。std来自全部4800 ensemble；grid等权散布不自动是测量噪声posterior。旧 notebook 为histogram拟合mean/std，当前pipeline采用其ensemble summary实现；不应混淆版本。

已检查当前档：
`/home/canoziia/GRMHD_DATA/analysis/m65-211-input-derived-v1/eht/selection.npz`

| key | shape / 含义 |
|---|---|
| `fields` | (7,) identity fields |
| `identity` | (4800,) JSON encoded common identity |
| `selected_identity` | (1572,) |
| `sorted_indices` | (4800,) score排序，common轴 |
| `selected_indices` | (1572,) common轴 |
| `raw_sorted_indices` | (4,4800)，common轴→各天pkl原axis索引 |
| `selected_raw_indices` | (4,1572) |
| `score` | (4800,) EVPA score |
| `rank_mask` | (64,64) |
| `days` | (4,) =95,96,100,101 |

**本次已用实际metadata验证** `identity,days,raw_sorted_indices` 全相符。四天 raw row0 指向不同 identity，不能拿同一原行号跨日配对。

历史档 `EHT_ZC/audit/topset_ranking_comparison.npz` keys：`legacy_score,corrected_score,legacy_order,corrected_order` 各(4800,)，`selected_pixels` (933,)；order 是按7-key排序后的共同轴，前1572为对应历史选择。不是pkl raw轴。`EHT_ZC/scripts/validate_topset_ranking.py` 为其producer代码；本次只读已有archive，没有重跑。

当前EVPA1572和以下集合的交集（均映射到共同identity轴后检查）：

| 对照 | 交集 /1572 |
|---|---:|
| historical legacy EVPA1572 |1252|
| corrected histogram EVPA1572 |1400|
| 同raw候选四天 mean-I chi2和前1572 |**680**|
| 同raw候选四天 (m+pvis)/2 和前1572 |643|
| day95 mean-I chi2前1572 |627|
| day96 mean-I chi2前1572 |729|
| day100 mean-I chi2前1572 |550|
| day101 mean-I chi2前1572 |797|

后两种 fourday 排名只是metadata级对照，不是发现了原 producer 自动执行的四日 selection。`/home/canoziia/GRMHD_DATA/source-backup/EHT/eht-imaging/Gaussian_fit.ipynb`存在早期按四天chi2和排序的notebook探索，而当前主pipeline最终按EVPA中心性选择，需以运行artifact为准。680交集显示选择规则不同确实足以显著换样本；**不能从索引交集直接断言daydiff漂移已被证明由TopSet造成**。

## 7. 版本/install/Git 缺口（停止深挖后的明确状态）

- pipeline docstring 引用 eht-imaging v1.1.0；脚本 `__version__='1.0'` 是pipeline版本，**不是运行库版本pin**。
- 代码备份 `EHT/eht-imaging/README` 明言：安装 ehtim 后要用本地 `image.py,imager.py,io/load.py` 替换。但 selected/excluded records 中未找到这些库文件，也未在所查EHT目录找到ehtim package/其Git。
- targeted Git backups 为2023-D01-01、EHT-axion-pixel-by-pixel、raptor、test_raptor_difference；remote/HEAD tree未定位这些ehtim实现。未遍扫整个home，也未跑全Git历史取证。
- 当前系统Python无 numpy/astropy/ehtim；项目default Python有numpy，无astropy/ehtim。所查另一checkout/default/eht-research package目录未定位astropy/ehtim。没有安装新库或下载源码。
- 故“具体版本、commit、chi2归一化、ratio定义、P阶段精确obs绑定”均列未查明，不能凭现代ehtim的常见行为回填成历史事实。

## 8. 给父 driver 的可运行最小输入与下一步

### 不需要重建/新datafit即可做

1. 以 `D/poloutput/230_*/*.fits` 为 **raw/base sensitivity case**；四天现有4800足够。或使用父已有trusted-pickle driver读取相应raw I/Q/U。不要换成 `pol_sorted` 后还沿用raw指标。
2. metadata筛选、mean-I排序、I闭合条件、m/pvis条件和joint/perday条件可直接使用已保存指标；明确标 `stored_optimizer_stage`。
3. 保留 original/legacy/current EVPA1572 baseline，对比不同条件及N，不把新条件称官方TopSet。identity交集后分别读各天对应raw索引；可自然诊断 selection 对daydiff的影响，而不是先预设结论。
4. `/tmp/eht-local-provenance/raw-{day}.indexed.json` 每行包含13原keys + `raw_index,common_identity_index`；`identity-index.json` 提供字段、共同identity与四天映射。父runtime不需引入审计reader。

### 仅缺指标时做少量重评估

- **需要 I visibility-amplitude `amp`：**metadata无该key；不能把 `vis` 改名。用少量 existing raw FITS + 对应 lo UVFITS，明确 averaging、zbl、noise、taper、selfcal convention，计算amp及需要的closure指标，先验证约定再扩展。
- **需要 rho 或真正 polarization closure 指标：**当前metadata没有，FITSheader样本也没有；不能以 `m` 或 `pvis` 替代。恢复可信ehtim版本/API后，在少量图上新增该datafit计算。
- **需要 final-restored-FITS 的一致χ²：**重新forward-evaluate final FITS（必要时限定gain solve），不是重做image optimization；记录这是新评估，并明确与历史stage诊断不可直接等同。
- **需要历史原阶段χ²精确重现：**final FITS不一定足够，缺原未restored图/obs/gains和自定义库；不能承诺仅靠少量recompute可精确还原。这与“从现有图像建立一致的新fit指标”是两个任务。
- `pol_sorted` day100重复日志不宜未经消歧进入science。先用原raw分支即可，不必为此重跑全候选。

## 9. 机器文件说明

目录 `/tmp/eht-local-provenance/`：

- `inventory.json`：总清单，包含全部grid/summary/counts/缺口。
- `metadata-audit.json`：12份pkl metadata summary、全部日志核对和8组FITS目录完整性；只有小txt日志计算SHA。
- `raw-*.rows.json`, `raw-*.indexed.json`：原metadata及driver映射；`analysis-*`/`polsorted-*`为分支核查。
- `selection-diagnostics.json`：底层分位数、各阈值perday/joint、noise分组、selection交集。
- `identity-index.json`, `*.common-indices.json`：共同轴映射与候选集合，可被标准JSON读入。
- `fits-headers.json`：32份小FITS所有HDU headers、16对小文件payload对照及8份lo/hi UVFITS主header。

一次性extract只在`meta_data`之后、第一大数组之前停止，绝大多数进程RSS约21–25MB；逐文件处理，不mmap pickle、不复制数组、不整读day100 885MB pkl。首个raw95记录的进程ru_maxrss约542MiB为异常高水位（不据此宣称21MB保证），仍未超过约600MiB级单文件限制；其实际读取metadata前缀仅1396108 bytes。所有临时reader只服务本报告，**不得移入driver/runtime**。
