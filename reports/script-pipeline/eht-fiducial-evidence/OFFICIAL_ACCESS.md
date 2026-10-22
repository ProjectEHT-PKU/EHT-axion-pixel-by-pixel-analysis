# April 11 polarized fiducial：官方访问与可执行 API 交接

## 立即下一步（给父任务）

**无需等大包：官方 April 11 low-band、部分 D-term 校正且已 Stokes-I 自校准的 HOPS UVFITS 已取到。**

- 本地：`/tmp/eht-fiducial-rebuild/access/Apr11-lo.uvfits`
- 官方固定 URL：<https://raw.githubusercontent.com/eventhorizontelescope/2023-D01-01/b71399714fdfe13579c741075d757b4ce65fd3d1/hops_data/April11/SR2_M87_2017_101_lo_hops_zbl-dtcal+selfcal.uvfits>
- 实测 HTTP 200，524160 bytes；Git blob SHA1 `97521a035216a685e15c82fd93ee78d246296726` 与官方 tree 一致。
- SHA256 `5cf936dce926b5b13a0f39065096da1f836ab52ca95fdfffef581057029e8d5b`。
- 只读 FITS header：`DATE-OBS=2017-04-11`，`OBJECT=M87`，`MJD=57854`，`CRVAL4=227070703125 Hz`；random-groups，5877 groups，STOKES -1,-2,-3,-4（RR/LL/RL/LR），**不是 image FITS**。
- 数据 README：<https://raw.githubusercontent.com/eventhorizontelescope/2023-D01-01/b71399714fdfe13579c741075d757b4ce65fd3d1/README.md>。该输入已校 ALMA/APEX/JCMT/SMA 泄漏和相位/振幅；剩余 SMT/LMT/PV 为主要待解站，站码 `AZ,LM,PV`。不要对 CASA 全 D-term 校正产品再次套同一求解。

**结论边界：未找到已验证的 Paper VII Fig.6 April11 eht-imaging image FITS，也未找到可逐行认证为 VII Appendix C/G 完整 survey / D-term alternating pipeline 的公开脚本。当前可立即做的是“采用官方数据和论文参数的重新实现”，不能称已取得官方原图或 exact pipeline。**

## 新版 ehtim：上游已经核实，安装态由父任务确认

- 官方仓库：<https://github.com/achael/eht-imaging>；GitHub repo API 确认默认分支为 **main**，不是 master。
- 本次 main tree：`7acf089bfe57f1463219a8e9400c7a18f106d9f1`，与之前下载 tree 一致。
- PyPI API <https://pypi.org/pypi/ehtim/json> 当前 `info.version=1.3.2`。
- GitHub release API <https://api.github.com/repos/achael/eht-imaging/releases?per_page=5> 首项 `v1.3.2`，published `2026-06-16T20:54:56Z`。
- 固定 release 元数据 <https://raw.githubusercontent.com/achael/eht-imaging/v1.3.2/pyproject.toml> 确认 `version="1.3.2"`；main 同为 1.3.2，故只看版本字符串不能区别 main 与 release。
- main 与 v1.3.2 的 `ehtim/imager.py` 本次均取到；这里不安装、不重复父任务 pilot，也不声称当前父环境已经升级成功。
- 父任务可在其实际运行解释器执行以下检查（本任务未执行）：

```python
import sys, importlib.metadata as md, inspect, ehtim
from ehtim.calibrating.pol_cal import leakage_cal
from ehtim.calibrating.pol_cal_new import leakage_cal_new
print(sys.executable, md.version('ehtim'), ehtim.__file__)
print(md.distribution('ehtim').read_text('direct_url.json'))
print(inspect.signature(ehtim.imager.Imager.make_image_P))
print(inspect.signature(leakage_cal))
print(inspect.signature(leakage_cal_new))
```

## 官方代码：链接与 App.C/G → 实际 API

固定代码根：`https://raw.githubusercontent.com/achael/eht-imaging/7acf089bfe57f1463219a8e9400c7a18f106d9f1/`

1. **官方极化 imaging 示例**：<https://github.com/achael/eht-imaging/blob/7acf089bfe57f1463219a8e9400c7a18f106d9f1/examples/example_pol.py>
   - `eh.imager.Imager(...)`、`make_image_P()`、`dat_term_next`、`reg_term_next`、`transform_next` 均有实际调用。
   - 示例没有 leakage solve/alternating loop；其注释写 ratio 但某些步骤实设 `pvis`，不能照注释误认为已使用 `m`。
2. **官方 polarization notebook**：<https://github.com/achael/eht-imaging/blob/7acf089bfe57f1463219a8e9400c7a18f106d9f1/tutorials/ehtim_tutorial_polarization.ipynb>
   - 有 `{'m':1}`、`{'hw':1,'ptv':1}`、`make_image_P()` 和 `make_image_IP()`；也没有完整 VII D-term loop。
3. **旧/历史泄漏 API**：<https://github.com/achael/eht-imaging/blob/7acf089bfe57f1463219a8e9400c7a18f106d9f1/ehtim/calibrating/pol_cal.py#L43>
   ```python
   from ehtim.calibrating.pol_cal import leakage_cal
   calibrated = leakage_cal(
       base_obs, im=current_pol_image, sites=['AZ','LM','PV'],
       leakage_tol=.1, pol_fit=['RL','LR'], dtype='vis',
       const_fpol=False, inverse=False, ttype='direct', show_solution=False)
   ```
   这是**官方 API 调用示意，不是找到的 VII 原脚本**。`leakage_tol=.1` 是库默认，不是本次认证的论文 fiducial。`sites=[]` 会解全部站，不能无意重解已校四站。该路径在 `const_fpol=False,inverse=False` 时调用 `apply_jones_inverse`，不是 const_fpol 的线性近似分支。完整 Jones inverse 含 `1/(1-dL*dR)`，见 `observing/obs_simulate.py:894`。
4. **新泄漏 API（2024 实现）**：<https://github.com/achael/eht-imaging/blob/7acf089bfe57f1463219a8e9400c7a18f106d9f1/ehtim/calibrating/pol_cal_new.py#L52>
   ```python
   from ehtim.calibrating.pol_cal_new import leakage_cal_new
   calibrated = leakage_cal_new(
       base_obs, current_pol_image, sites=['AZ','LM','PV'],
       leakage_tol=.1, rescale_leakage_tol=False,
       pol_fit=['RL','LR'], dtype='vis', ttype='direct',
       use_grad=True, show_solution=False, apply_solution=True)
   ```
   `base_obs.ampcal` 与 `.phasecal` 必须真；有 `apply_solution`，旧 API 没有此参数。源码说明 `rescale_leakage_tol=True` 会修正旧版站数缩放问题，所以不能把新实现默认为论文数值完全相同。新函数不接收 `const_fpol` / `inverse`。

### 论文 fiducial 参数的明确映射

已有论文原文核验见项目只读报告 `reports/eht-polarization-research/OFFICIAL_METHOD.md` §7。

| Paper VII | 当前公开 API | 注意 |
|---|---|---|
| αP=1, αm=1 | `data_term={'pvis':1.,'m':1.}` 或 `imgr.dat_term_next=...` | 不是只用 pvis |
| βHW=100, βTV=1 | `reg_term={'hw':100.,'ptv':1.}` | 是 polarimetric TV `ptv`，非 I 的 `tv` |
| I 固定的极化步 | `imgr.make_image_P(show_updates=False)` | 不用 `make_image_IP` |
| 物理 polarization transform | 构造器 `transform=['mcv']` | 当前构造器会 flatten 成数组；后续赋值优先 `imgr.transform_next=['mcv']` |
| 64²、120 μas、0.6 Jy | `eh.image.make_square(obs,64,120*eh.RADPERUAS)`；固定 I 应来自相应 I 重建 | 方格并不提供官方 I prior |
| scan average | `obs.add_scans(); obs=obs.avg_coherent(0,scan_avg=True)` | 实际 scan 边界需记录；自动识别不等于原发布 scan table |
| 去 zero baselines | `obs.flag_uvdist(uv_min=...)` | exact cutoff 未获官方脚本，不凭空给阈值 |
| niter=50 | 外层 image ↔ leakage alternating loop | **不是** `Imager(maxit=50)`，也不能直接等同 `make_image_P(niter=50)` |
| 20 μas restart blur | `im.blur_circ(20*eh.RADPERUAS)` 等 | 再保持固定 I 的合同；`blur_circ` 第二参数涉及极化相对模糊，不能盲套 |
| Fig.6 最终显示 | eht-imaging 19 μas kernel | 别拿 restart 的 20 μas 替代最终核 |

### 必须马上处理的 0.5% noise API 陷阱

- `Imager.__init__` 确实接受 `systematic_noise=0.005`，但 `imager.py:1123–1149` 的 polarimetric 分支调用 `polchisqdata` **没有传它**；只有普通 I data branch 传递此参数。
- `pol_imager_utils.py:1412+` 的 `pvis/m` 数据准备直接读 `psigma/msigma`，无 0.5% 注入。故单写 `systematic_noise=.005` 不等于实现 VII 极化 fsys。
- 公开现成 API：`obs_with_noise = obs.add_fractional_noise(0.005, debias=False)`（`obsdata.py:1937`），源码用 Stokes I visibility amplitude 给四种 sigma quadrature 加噪。应在明确的处理阶段**只加一次**，重建新的 Imager 或确保清空/重建 data tuple cache。
- 仍需核对原论文到底对哪一 polarization representation 的误差分量增加该量及 ratio 的传播：该 API 是可用落点，不是已经证明与论文每个 factor 相等。特别不要同时外部 inflate 和 Imager 再 inflate；D-term 校准 likelihood 是否用相同 inflated errors 也需父任务明确。
- 交替时从同一个部分校正 `base_obs` 重新求绝对剩余 D-terms；不要每轮对上一轮已应用过 D-term 的数据盲目累计纠正。旧 API 返回 `calibrated.tarr['dr'/'dl']` 和已经校正的数据，注意保留基准。

## image FITS / exact pipeline 的 not-found 边界

- 官方 `2023-D01-01` 完整 tree（固定 b713997…）只有 24 UVFITS 及说明/license，无 image FITS，无极化重建脚本。
- EHT GitHub org API 返回的公开仓库中没有专门 Paper VII polarization pipeline；`2021-D02-01` 是多波段 campaign，`2021-D03-01` 是 Centaurus A，不是目标图。
- `2017-april` root 是旧的 er* calibration 目录；递归 tree 超本任务 2MB 响应限而中止，未声称穷尽这个仓库。
- 官方 data portal 本次仍 HTTP 403：<https://eventhorizontelescope.org/for-astronomers/data>。不绕过身份，也不由 403 推论没有公开数据。
- Zenodo API `M87 AND polarization` 命中 104 个，核看前 5 个无官方目标；`"M87" AND "fiducial"` 4 条也没有目标。前者不是穷尽检索；未下载返回列表中 >100MB 的 simulation 包，且这些并非已验证目标。
- Caltech API `M87 AND polarization` 返回 0；宽松 query 返回无关记录。不是整个 Caltech 档案的不存在证明。
- CyVerse 旧猜测路径本次重定向循环，无权威 listing；不把 Angular shell/重定向当数据存在证据，不继续猜目录。
- 因此现阶段**可验证 small FITS URL = 上述官方 UVFITS；可验证官方 image FITS URL = 未找到**。真正逐像素验证 fiducial 仍需作者/官方给出图像 manifest 或原脚本、版本、prior/calibration provenance。

## 本任务执行边界与证据

- 仅写 `/tmp/eht-fiducial-rebuild/access/*` 与本报告；未改项目源、data 或 outputs；未安装、不跑成像、不使用 Goal。
- 每个 curl 均 `--max-time 20`，单请求 `--max-filesize` 1.5–2MB；总获取为小文件级（远小于 100MB），没有下载大包。递归 tree 超限的截断响应已保留但不用于完整性推论。
- URL 请求表：`access/urls{,2,3,4}.json`；repo/PyPI/release metadata、代码和 header 在同目录。旧已有 README/tree/示例位于 `/tmp/eht-official-method/`。
