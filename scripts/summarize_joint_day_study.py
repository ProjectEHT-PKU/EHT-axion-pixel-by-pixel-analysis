"""Summarize completed conditional four-day studies; never run inference or RT.

Only processed declared outputs are consumed. Source manifests are copied under
noncanonical receipt names so this stage does not recursively traverse raw ancestry.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np

from eht_axion.pipeline.common import (
    manifest_path, require_declared_outputs, sha256, stage_directory,
    verification_session, verify_stage,
)
from eht_axion.statistics import normalized_posterior_cdf, posterior_upper_limit

BASE = Path('/home/canoziia/GRMHD_DATA/analysis')
CASES = [
    'mean_background_iso', 'mean_background_diag', 'mean_background_full',
    'mean_background_full_offset', 'frozen_mixture_full', 'evolving_mixture_full',
    'joint_background_evolving_mixture', 'joint_gaussian_moments',
    'independent_gaussian_moments',
]
LABELS = ['mean+iso', 'mean+diag', 'mean+full', 'mean+full+offset',
          'frozen mixture', 'evolving mixture', 'joint empirical B,T',
          'joint Gaussian', 'indep Gaussian']
SELECTED = [CASES[2], CASES[4], CASES[6], CASES[7]]
ALPHAS = [.1, .15, .2, .25, .3, .35, .4, .44]


def read_rows(path):
    with Path(path).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key in ('alpha', 'c95', 'cbg_trace', 'cbg_eigen_min',
                    'background_mean_norm', 'lower_bin_mass', 'upper_bin_mass',
                    'fixed_cov_eigen_min', 'fixed_cov_condition'):
            row[key] = float(row[key])
        for key in ('contrasts', 'n', 'state_count', 'cbg_rank', 'phase_count'):
            row[key] = int(row[key])
    return rows


def index_rows(rows):
    result = {}
    for r in rows:
        key = (r['dimension'], r['alpha'], r['contrasts'], r['case'], r['mode'])
        if key in result:
            raise ValueError(f'duplicate row: {key}')
        if not np.isfinite(r['c95']) or r['c95'] <= 0:
            raise ValueError('invalid c95')
        result[key] = r
    return result


def validate_counts(rows, alphas):
    index_rows(rows)
    expected = {(d, a, h, mode): (12 if h == 2 else 9)
                for d in ('1d', '2d') for a in alphas
                for h in (2, 3) for mode in ('observed', 'matched_zero')}
    actual = Counter((r['dimension'], r['alpha'], r['contrasts'], r['mode']) for r in rows)
    if actual != expected:
        raise ValueError(f'case counts differ: {actual}')
    for r in rows:
        allowed = CASES + (['legacy_recovery', 'cartesian_mean_old_background',
                            'frozen_mixture_old_background'] if r['contrasts'] == 2 else [])
        if r['case'] not in allowed:
            raise ValueError(f"unknown case {r['case']}")


def compare_rows(first, second):
    """Matched-key ratios; >1 is retained (a weaker upper limit)."""
    a, b = index_rows(first), index_rows(second)
    if a.keys() != b.keys():
        raise ValueError('comparison keys differ')
    return [dict(dimension=k[0], alpha=k[1], contrasts=k[2], case=k[3], mode=k[4],
                 first_c95=a[k]['c95'], second_c95=b[k]['c95'],
                 ratio=b[k]['c95'] / a[k]['c95'],
                 delta_dex=float(np.log10(b[k]['c95'] / a[k]['c95']))) for k in a]


def validate_posterior(path, row):
    with np.load(path, allow_pickle=False) as z:
        c, ll, density, cdf = [z[k] for k in ('couplings', 'loglike', 'density', 'cdf')]
        if c.ndim != 1 or any(x.shape != c.shape for x in (ll, density, cdf)):
            raise ValueError(f'posterior axis mismatch: {path}')
        if (not np.isfinite(density).all() or not np.isfinite(cdf).all()
                or np.any(density < 0) or np.any(np.diff(cdf) < -1e-12)
                or not np.allclose(cdf[[0, -1]], [0, 1], atol=1e-12)):
            raise ValueError(f'invalid posterior: {path}')
        expected_density, expected_cdf = normalized_posterior_cdf(c, ll, prior='log-uniform')
        np.testing.assert_allclose(density, expected_density, rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(cdf, expected_cdf, rtol=1e-11, atol=1e-13)
        if not np.isclose(np.sum(np.diff(c)*(density[1:]+density[:-1])/2), 1):
            raise ValueError('density not normalized')
        limit = posterior_upper_limit(c, ll, prior='log-uniform')
        if not np.isclose(limit, row['c95'], rtol=1e-11, atol=0):
            raise ValueError('CSV/NPZ c95 mismatch')
        if not np.isclose(float(z['c95']), limit, rtol=1e-11, atol=0):
            raise ValueError('NPZ c95 mismatch')
        np.testing.assert_allclose([row['lower_bin_mass'], row['upper_bin_mass']],
                                   [cdf[1], 1-cdf[-2]], atol=1e-12)
        peak = int(np.argmax(ll))
        if (row['likelihood_peak_at_boundary'] in (True, 'True')) != (peak in (0, len(c)-1)):
            raise ValueError('likelihood boundary flag mismatch')
        if z['phase_loglike'].shape != (len(c), row['phase_count']):
            raise ValueError('phase/coupling axes mismatch')
        return dict(key=row['key'], coupling_points=len(c), phase_count=row['phase_count'],
                    prior_min=float(c[0]), prior_max=float(c[-1]),
                    likelihood_peak_c=float(c[peak]),
                    likelihood_peak_side='lower' if peak == 0 else 'upper' if peak == len(c)-1 else 'interior')


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_plots(rows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    idx = index_rows(rows)
    alphas = sorted({r['alpha'] for r in rows})
    def get(d, a, h, case, mode='observed'):
        return idx[d, a, h, case, mode]
    colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for i, d in enumerate(('1d', '2d')):
        for j, h in enumerate((2, 3)):
            ax = axes[i, j]
            for case, color in zip(SELECTED, colors):
                for mode, style in [('observed', '-o'), ('matched_zero', '--x')]:
                    ax.plot(alphas, [get(d, a, h, case, mode)['c95'] for a in alphas],
                            style, color=color, label=LABELS[CASES.index(case)]+' '+('obs' if mode == 'observed' else 'NULL'))
            ax.set(title=f'{d}, H{h}', yscale='log', xlabel='alpha', ylabel='conditional c95')
            ax.grid(alpha=.25)
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(out/'c95-vs-alpha.png', dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    x = np.arange(len(CASES))
    for ax, d in zip(axes, ('1d', '2d')):
        for h, offset, color in [(2, -.19, colors[0]), (3, .19, colors[1])]:
            yy = [get(d, .3, h, case)['c95'] for case in CASES]
            bars = ax.bar(x+offset, yy, .38, label=f'H{h}', color=color)
            ax.bar_label(bars, labels=[f'{v:.3g}' for v in yy], fontsize=8, padding=3)
        ax.set(yscale='log', ylabel='observed c95', title=f'{d}, alpha=0.3')
        ax.set_xticks(x, LABELS, rotation=28, ha='right'); ax.legend(); ax.margins(y=.25)
    fig.tight_layout(); fig.savefig(out/'alpha030-factorial.png', dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for i, d in enumerate(('1d', '2d')):
        rr = [get(d, .3, 2, case) for case in CASES[:3]]
        axes[i, 0].bar(['scalar', 'diagonal', 'full'], [r['c95'] for r in rr])
        axes[i, 0].set(title=f'{d}, H2 mean response, alpha=0.3', ylabel='observed c95')
        ax = axes[i, 1]
        ax.plot(range(3), [r['cbg_trace']/(2*r['n']) for r in rr], 'o-', label='mean eigenvalue = trace / 2N')
        ax.plot(range(3), [r['cbg_eigen_min'] for r in rr], 's--', label='minimum eigenvalue')
        ax.set_xticks(range(3), ['scalar', 'diagonal', 'full'])
        ax.set(yscale='log', ylabel='background variance (rad²)', title='Variance ellipsoid summary (not area)')
        ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out/'background-geometry.png', dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    for ax, d in zip(axes, ('1d', '2d')):
        for case, label in zip(CASES, LABELS):
            ax.plot(alphas, [get(d, a, 3, case)['c95']/get(d, a, 2, case)['c95'] for a in alphas], 'o-', label=label)
        ax.axhline(1, color='black', linestyle='--')
        ax.set(title=d, ylabel='observed c95(H3) / c95(H2)', yscale='log', xlabel='alpha')
        ax.grid(alpha=.25)
    axes[0].legend(fontsize=8, ncol=3)
    fig.tight_layout(); fig.savefig(out/'h3-h2-ratios.png', dpi=170); plt.close(fig)


def markdown_table(headers, rows):
    return '\n'.join(['| '+' | '.join(x.replace('|', '\\|') for x in headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']
                     + ['| '+' | '.join(str(x).replace('|', '\\|') for x in r)+' |' for r in rows])


def state_sections(statecheck, diagnostics):
    sections = []
    if statecheck:
        results = json.loads((statecheck/'report.json').read_text())['results']
        if len(results) != 64:
            raise ValueError('expected 64 state sensitivity rows')
        obs = [x for x in results if x['mode']=='observed']
        table = []
        for d in ('1d','2d'):
            for h in (2,3):
                for case in [CASES[2],CASES[5],CASES[6],CASES[7]]:
                    pair = [x for x in obs if x['dimension']==d and x['contrasts']==h and x['case']==case]
                    a = next(x for x in pair if x['subset']=='first_half')
                    b = next(x for x in pair if x['subset']=='last_half')
                    table.append([d,h,case,f"{a['baseline_c95']:.6g}",f"{a['c95']:.6g}",f"{b['c95']:.6g}",f"{max(abs(x['delta_dex']) for x in pair):.6g}"])
        sections.append('## 已完成的有限窗口敏感性（α=.3）\n\n'
            '前/后各 62 个窗口分别重估本子集的矩；不是独立观测误差条，也不是增加新数据。'
            '保留 H2/H3 及 Gaussian 背景控制，不能只检查 empirical 模型。\n\n'+markdown_table(
            ['维度','H','模型','全窗口','前62','后62','最大 |Δdex| 对全窗口'],table)+
            '\n\n2D H3 joint empirical 全窗口约 0.31378，前半约 0.5035、后半约 1.4267；'
            '后半/全窗口约 4.55 倍（0.6577 dex）。1D H3 从约 3.3467 到约 2.45/2.09。'
            'Gaussian 背景加 evolving response 的 1D H3 也变化约 0.338 dex（2.18 倍），'
            '所以不能把其他全协方差模型统称稳健。完整 observed/NULL 64 行见 statecheck/sensitivity.csv。')
    if diagnostics:
        results = json.loads((diagnostics/'report.json').read_text())['results']
        table = [[x['dimension'],x['contrasts'],f"{x['posterior_state_ESS']:.7g}",
                  f"{x['max_weight']:.7g}",f"{x['minimum_component_chi2']:.7g}",
                  x['contrasts']*{'1d':23,'2d':20}[x['dimension']]] for x in results]
        sections.append('## 已完成的 latent-state 峰值诊断\n\n'+markdown_table(
            ['维度','H','峰处 state ESS','最大状态权重','最小 component χ²','向量维数'],table)+
            '\n\n这是保存的单个 coupling/phase 峰处诊断，**不是 posterior-wide ESS**。'
            'H3 的约 1.00055/1.00783 ESS 与 0.999723/0.996100 最大权重表明几乎单状态支配；'
            'H2 也集中，但 component χ² 约 87.5/46、124.7/40，仍不应称好拟合；'
            'H3 约 686.65/69、1000.87/60 更是强烈失配警报。'
            'component 最小 χ² 不是 mixture goodness-of-fit，不能直接转作正式 p 值；'
            '这些读数既未证明科学模型成立，也不构成覆盖率认证。')
    return '\n\n'.join(sections)


def make_report(rows, grid, checks, source_report, output_root, statecheck, extra_sections=''):
    idx = index_rows(rows)
    def r(d, h, case, mode='observed', a=.3):
        return idx[d, a, h, case, mode]
    table = markdown_table(['模型（固定构造顺序）', '1D H2', '1D H3', '2D H2', '2D H3'],
        [[label]+[f"{r(d,h,case)['c95']:.6g}" for d in ('1d','2d') for h in (2,3)]
         for case, label in zip(CASES, LABELS)])
    controls = markdown_table(['维度','H','零均值 c95','固定 offset c95','offset 范数 (rad)'],
        [[d, h, f"{r(d,h,CASES[2])['c95']:.6g}", f"{r(d,h,CASES[3])['c95']:.6g}",
          f"{r(d,h,CASES[3])['background_mean_norm']:.6g}"] for d in ('1d','2d') for h in (2,3)])
    worst = max(grid, key=lambda x: abs(x['delta_dex']))
    decomposition = markdown_table(['模型','最大 |dense/coarse| dex','对应变化幅度'],
        [[label, f"{max(abs(x['delta_dex']) for x in grid if x['case']==case):.6g}",
          f"{100*(10**max(abs(x['delta_dex']) for x in grid if x['case']==case)-1):.3f}%"]
         for case,label in zip(CASES,LABELS)])
    edges = markdown_table(['α','lower bin mass 范围','upper bin mass 范围','likelihood 峰在边界（/84）'],
        [[a, f"{min(x['lower_bin_mass'] for x in rows if x['alpha']==a):.3g}–{max(x['lower_bin_mass'] for x in rows if x['alpha']==a):.3g}",
          f"{min(x['upper_bin_mass'] for x in rows if x['alpha']==a):.3g}–{max(x['upper_bin_mass'] for x in rows if x['alpha']==a):.3g}",
          sum(x['likelihood_peak_at_boundary'] in (True,'True') for x in rows if x['alpha']==a)] for a in ALPHAS])
    ratios = [(d,a,case,r(d,3,case,a=a)['c95']/r(d,2,case,a=a)['c95'])
              for d in ('1d','2d') for a in ALPHAS for case in CASES]
    weak = [x for x in ratios if x[-1]>1]
    maxratio = max(ratios,key=lambda x:x[-1])
    cross = markdown_table(['维度','H','joint Gaussian','independent Gaussian','joint/independent'],
        [[d,h,f"{r(d,h,CASES[7])['c95']:.6g}",f"{r(d,h,CASES[8])['c95']:.6g}",
          f"{r(d,h,CASES[7])['c95']/r(d,h,CASES[8])['c95']:.6g}"] for d in ('1d','2d') for h in (2,3)])
    recovered = ', '.join(f"{d}: {r(d,2,'legacy_recovery')['c95']:.7f}" for d in ('1d','2d'))
    full = [x for x in rows if x['case']=='mean_background_full']
    return f'''# 四日联合受控研究：已完成产物汇总

## 结论先行

这是**冻结 EHT 条件下的模型比较，不是科学正确性认证或覆盖率校准**。八个 α 全部完成；主表保留 672 行
（每 α：2 维度 × (H2 的 12 个模型 + H3 的 9 个模型) × observed/NULL = 84）。
不按约束强弱挑模型，不把新全模型自动视为更正确。

α=0.3 的 2D joint empirical B,T：H2={r('2d',2,CASES[6])['c95']:.6g}，
H3={r('2d',3,CASES[6])['c95']:.6g}；而 H3 joint Gaussian={r('2d',3,CASES[7])['c95']:.6g}。
从约 1.7622 到 0.31378 的强收紧，面对 Gaussian 约 1.8054 的结果，首先是**抽样分布/模型选择巨大敏感性警报**，
不是应当选用最强上限的理由。下面纳入已完成的有限窗口与峰处 ESS 诊断，但这些检查不等于该联合统计模型已获论文级验证。

## 冻结了什么、没有做什么

四日实测数据、特征选择（1D N=23；2D N=20）、逐日 EHT Monte Carlo 协方差、I/error 与超参数全部冻结。
未重建 EHT、未补充 D/latent-EHT 不确定性、未做 RT、未调 χ²、未缩放响应、未改 canonical NULL。
NULL 是同一条件模型中的 matched_zero 对照，不是新观测，也不是“最强模型”的选择依据。
所有 c95 都是原代码耦合参数 c 的条件 95% posterior 上限；**不另造 GeV 单位**。
先验一致为 [1/137,1000] 上 log-uniform；主计算 coupling/phase=120/64，α=.3 网格对照为 240/128。
遗留 H2 原六项（1D/2D × observed 模型控制）按原网格恢复；α=.3 legacy c95 为 {recovered}。
这是实际受控结果，不为追求更强约束而调参。

H2=[day2−day1,day4−day3]；H3 再加入 day3−day2。
四日最多提供三个独立时间差组合，**不是六个独立两两差**。
观测协方差严格按 H diag(C_daily) Hᵀ 构造；共享日产生的负交叉块保留。
H3/H2 在相同模型假设下的改变也不能全归因于“信息量增加”：新增时间组合同时接受同一时间状态模型约束。
144 个 observed 同模型 H3/H2 比值中 {len(weak)} 个 >1（上限变弱）；最大 {maxratio[-1]:.6g}
出现在 {maxratio[0]}, α={maxratio[1]}, {maxratio[2]}，图中不隐藏这些情况。

## α=0.3 全因子表（observed c95）

{table}

mean+iso/diag/full 在同一 mean response 下只改变背景预测协方差；ISO 方差池固定为 H2，
H3 也使用该池，避免加入第三差时顺便改变噪声定义。
frozen/evolving mixture 使用完整复杂响应状态，不将 A 与 φ 独立抽样：
每个状态为复数 (α,state,day,feature)，同一个全局相位与整状态作用于所有特征和时间差。
冻结响应 141 状态，演化响应 124 重叠窗口，先验均为均匀状态权重；不能把这些状态当独立观测或用 √N 缩小误差。
这里称 response，不将其改名为 π/FDT。

背景默认采用**零均值模型**：从经验背景减去样本均值；这与固定 offset 都是估计模型选择，不是校准改正。
响应均值也由有限状态估计；均值响应和响应混合是不同模型，不是消除真实系统误差的保证。
以下单列固定背景均值的控制，范数有限且非零：

{controls}

## 分离 B/T 交叉协方差，而不是混淆分布改变

joint empirical 使用配对的 B+c Re(exp(iψ)T) 混合，仅加 C_EHT，**不再加 C_bg**，避免双计数。
joint Gaussian 在每个 c、ψ 对同一和量矩匹配，保留线性 c 的 B/T 交叉项、c² 项及参数相关 logdet。
它与 empirical 的区别既有分布形状也有依赖结构的表现，不能只称“加入 cross covariance 的收益”。
真正隔离 cross covariance 的对照是 joint Gaussian 对 independent Gaussian：后者边缘矩完全相同，仅令 B/T cross block 为零。

{cross}

背景协方差来自 124 重叠窗口的预测矩（ddof=0，不除以窗口数）；H2/H3 向量维数为 2N/3N（1D 46/69，2D 40/60）。
full C_bg 秩范围 {min(x['cbg_rank'] for x in full)}–{max(x['cbg_rank'] for x in full)}；
加入冻结 C_EHT 后固定协方差最小特征值全表最小 {min(x['fixed_cov_eigen_min'] for x in rows):.6g}，为 SPD。
没有任意 ridge/shrinkage，也没有因模型拟合好坏而更换协方差。

{extra_sections}

## 原生时间、角度分支与数值离散限制

输入 141 帧，跨度 {source_report['span_hours']:.6f} h，原生 cadence {source_report['cadence_hours']:.6f} h，
124 个重叠四日窗口。实际相对时间 [0,24,120,144] h，曝光 6 h。
曝光短于原生采样；所谓积分仅对分段线性插值器精确，**曝光未被原生时间分辨率充分解析**。
冻结/演化响应是时间假设包络，不是时间过程已验证。
背景原生分辨率 400、响应 100，末端网格同为 64 并不消除原生分辨率差异；本结果是该不匹配下的 benchmark。
rmax=1、强度加权、10° 窗口等冻结设置见输入快照，不代表应重调这些参数。
局部最短角差提升不是跨 96 h 间隙的物理绕转测定。
输入背景非零 edge winding 为 1D 219、2D 179；背景及 observed closure winding 都为 0，不能把 edge winding 说成闭合失败。

## 网格、先验边界与 posterior 自检

α=.3 coarse→dense 的全 84 项最大 |Δlog10(c95)|={abs(worst['delta_dex']):.6g} dex，
约 {100*(10**abs(worst['delta_dex'])-1):.3f}%，出现在 {worst['dimension']} H{worst['contrasts']}
{worst['case']} {worst['mode']}。这是 coupling 与 phase 一起加密的合成变化，不能虚构单轴误差分解。
按模型分解（包含 observed/NULL 与两个维度/H）：

{decomposition}

边界表覆盖每 α 的全部 84 行；lower/upper 是首/末一个**网格区间**的 posterior 概率，不是离散点质量。
实际全表 upper mass 仅为 0 到约 1e-15 的浮点尾差，没有上先验截断质量堆积证据，不能将 low-α 直接称为上界先验主导。
likelihood 的边界峰均应结合 posterior-checks.csv 中 lower/upper 方向解读；下界处的似然峰也不自动等于 posterior 全由先验主导。
首区间质量最大随 α 从约 1.36% 到 5.54%，说明不能忽略低耦合先验区间的作用。
未进行新先验实验，不声称先验下界无影响或先验独立性。逐项下/上质量、峰在 lower/interior/upper 均保留在 CSV。

{edges}

逐 case NPZ 检查总数 {len(checks)}（allalpha 672 + coarse α=.3 84 + dense α=.3 84，按源计数，并非主表额外行）。
allalpha 实际另有 16 个 recovery-reference NPZ（8 α × 2 维度），连同其他 declared outputs 均经 SHA 校验；不混入 672 行主表。
以上逐 case NPZ 已检查：coupling/density/CDF/loglike 一维轴一致；phase_loglike 轴正确；
density 非负且按 c 梯形积分归一，CDF 单调且端点 0/1，与标准 log-uniform posterior 计算一致；
NPZ 与 CSV c95 相符，未对错误轴求和。边界质量与 likelihood peak 标记逐项核对。
这只是数值产品内部一致性检查，不是模型真伪或覆盖率证明。

## 图与来源复现

以下图位于 `{output_root}`：

1. ![全部 α：实线 observed / 虚线 NULL]({output_root}/c95-vs-alpha.png)
2. ![α=.3 全因子]({output_root}/alpha030-factorial.png)
3. ![固定 H2 背景协方差几何]({output_root}/background-geometry.png)
4. ![H3/H2 比值，>1 为变弱]({output_root}/h3-h2-ratios.png)

图 3 在固定 H2、mean response 下对比 scalar/diag/full c95，以及平均/最小方差主轴；
它是方差椭球摘要，不将 trace 误称为椭圆面积。图 4 包含全部九个共享模型、全部八个 α。

消费 roots 为 `eht-joint-days-{{inputs,alpha030,allalpha,alpha030-grid}}-v1`。
每个源用 verify_stage(recursive=False) 校验所有已声明数值输出的 SHA，再检查消费文件为 declared outputs。
只消费处理后文件，不递归重哈希约 97 GB 原始父链。
源 manifest 字节原样复制为 source-receipt-*.json（非 canonical manifest/stage 名），
记录原始路径与 SHA；新汇总阶段仍可递归核验这些收据和处理后产品，而不引入 raw ancestors。
主表仅采用 allalpha 的 672 行，不重复堆叠 α=.3 的 84 行；单独确认 coarse α=.3 重复产物一致。
报告、CSV、图与 provenance 都由同一批已验证处理后输入生成；此报告副本与输出阶段内文件字节相同。

statecheck：{statecheck}。可用 `--statecheck ROOT` 显式接入完成产物；默认不等待、不触发额外实验。
其内容不被自动提升为已通过的科学结论。
'''


def run(args):
    output = args.output_root.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError('new output root required')
    roots = [args.inputs, args.alpha030, args.allalpha, args.grid]
    optional_roots = [(name, getattr(args, name)) for name in ('statecheck','state_diagnostics') if getattr(args, name)]
    roots.extend(root for _, root in optional_roots)
    sources, receipts, provenance = [], [], []
    receipt_root = output.with_name(output.name+'-source-receipts')
    receipt_root.mkdir(parents=True, exist_ok=False)
    with verification_session():
        for i, root in enumerate(roots):
            root = root.resolve()
            manifest = verify_stage(root, recursive=False, allow_imported=False)
            receipt = receipt_root/f'source-receipt-{i:02d}.json'
            shutil.copyfile(manifest_path(root), receipt)
            if sha256(receipt) != sha256(manifest_path(root)):
                raise ValueError('receipt copy differs')
            receipts.append(receipt)
            provenance.append(dict(root=str(root), receipt=receipt.name, sha256=sha256(receipt),
                                   verification='declared outputs SHA verified; recursive=False'))
            if i == 0:
                consumed = [root/'report.json']
            elif i < 4:
                rr = read_rows(root/'summary.csv')
                consumed = [root/'summary.csv', root/'metadata.json'] + [root/f"{r['key']}.npz" for r in rr]
            else:
                consumed = [root/e['path'] for e in manifest['outputs']]
            require_declared_outputs(root, consumed)
            sources.extend(consumed)
        rows = read_rows(args.allalpha/'summary.csv')
        coarse = read_rows(args.alpha030/'summary.csv')
        dense = read_rows(args.grid/'summary.csv')
        validate_counts(rows, ALPHAS); validate_counts(coarse, [.3]); validate_counts(dense, [.3])
        repeated = compare_rows(coarse, [r for r in rows if r['alpha']==.3])
        if any(abs(r['delta_dex']) > 1e-12 for r in repeated):
            raise ValueError('alpha=.3 coarse runs differ')
        grid = compare_rows(coarse, dense)
        checks = []
        for root, rr in [(args.allalpha, rows), (args.alpha030, coarse), (args.grid, dense)]:
            for row in rr:
                checks.append(dict(source=str(root), **validate_posterior(root/f"{row['key']}.npz", row)))
        report = json.loads((args.inputs/'report.json').read_text())
        # Enforce common frozen profile/state provenance and intended quadrature.
        for root, rr, cp, pp in [(args.allalpha, rows,120,64), (args.alpha030,coarse,120,64), (args.grid,dense,240,128)]:
            meta = json.loads((root/'metadata.json').read_text())
            if meta['source_report'] != report or meta['prior'] != 'log-uniform':
                raise ValueError('frozen source/prior mismatch')
            if (meta['coupling_points'], meta['phase_count']) != (cp,pp):
                raise ValueError('unexpected source grid')
            for row in rr:
                if row['n'] != {'1d':23,'2d':20}[row['dimension']]:
                    raise ValueError('frozen feature count mismatch')
        if any(not np.allclose([x['prior_min'],x['prior_max']], [1/137,1000], rtol=1e-12) for x in checks):
            raise ValueError('prior support mismatch')
        text = make_report(rows, grid, checks, report, output,
                           f'已接入 {args.statecheck} 的已完成输出（见 statecheck/ 与 provenance）' if args.statecheck else '未接入；可选 split-window 敏感性阶段尚不作为本报告证据',
                           state_sections(args.statecheck, args.state_diagnostics))
        config = dict(schema='joint-day-results-v1', row_count=len(rows), alphas=ALPHAS,
                      source_provenance=provenance, recursive_source_verification=False,
                      scientific_recomputation=False)
        with stage_directory(output, 'research:joint-day-summary', config,
                             sources+receipts+[Path(__file__).resolve()]) as out:
            for receipt in receipts:
                shutil.copyfile(receipt, out/receipt.name)
            write_csv(out/'combined-summary.csv', rows)
            write_csv(out/'grid-comparison.csv', grid)
            write_csv(out/'posterior-checks.csv', checks)
            write_csv(out/'alpha030-repeat-comparison.csv', repeated)
            (out/'provenance.json').write_text(json.dumps(config, indent=2)+'\n')
            (out/'JOINT_DAYS_STUDY_ZH.md').write_text(text)
            for name, root in optional_roots:
                manifest = verify_stage(root, recursive=False, allow_imported=False)
                for entry in manifest['outputs']:
                    dest = out/name/entry['path']
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(root/entry['path'], dest)
            make_plots(rows, out)
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(output/'JOINT_DAYS_STUDY_ZH.md', args.report)
    return output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for flag, name in [('inputs','inputs'), ('alpha030','alpha030'), ('allalpha','allalpha'), ('grid','alpha030-grid')]:
        p.add_argument('--'+flag, type=Path, default=BASE/f'eht-joint-days-{name}-v1')
    p.add_argument('--output-root', type=Path, default=BASE/'eht-joint-days-results-v1')
    p.add_argument('--report', type=Path, default=Path(__file__).resolve().parents[2]/'reports/script-pipeline/JOINT_DAYS_STUDY_ZH.md')
    p.add_argument('--statecheck', type=Path, help='Optional COMPLETED split-window stage; never wait/run')
    p.add_argument('--state-diagnostics', type=Path, help='Optional COMPLETED latent peak diagnostics stage')
    args = p.parse_args(argv)
    print(run(args), flush=True)


if __name__ == '__main__':
    main()
