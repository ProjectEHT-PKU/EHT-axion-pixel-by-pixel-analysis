# 背景模型重新评估：完成与决定

用户补充原始统一sigma的有限141帧动机后，本次不是只为旧方案补充警告，而是重新作了预测验证并重算较低复杂度候选。

## 明确纠正

此前全样本full covariance能计算、加入EHT误差后能求逆，以及背景含相关性，都不足以支持将原始经验full作为更可靠的主背景模型。此次留出预测表明该升级不受支持，撤回其优先地位。旧结果、源码和清单保留，不重签或覆盖。

同理，124个离散背景／响应联合状态并不充分代表联合分布；此前峰值ESS近1和分半大变化已构成不采用其最强上限作主结论的理由。本次不再重跑这个脆弱方案。

保留原方案中正确的部分：四天三个独立时间差、共享日EHT误差交叉块、共同轴子相位、完整快照幅相结构及正确Gaussian归一化。

## 实际重做

冻结EHT均值、每日MC协方差、23/20特征选择、物理参数和先验。只用无轴子GRMHD背景评分，不以EHT的chi2或c95选模型。

- 5个连续测试块、3块敏感性；从训练集剔除任何与测试块共享原始帧的窗口，包括曝光插值端点。
- 训练窗口数5折为82/65/65/65/83，3折为65/49/66。没有把剩余窗口当独立观测。
- 所有候选训练期中心化估计散布；共同固定预测均值0；不减测试均值或全样本均值。
- 另加24小时隔离敏感性，不把无共享帧误称为时间独立。
- 固定候选：全局scalar；时间小矩阵⊗空间单位矩阵；0.5池化的逐特征方差；0.75时间目标+0.25经验空间协方差。权重不根据分数或限制再优化。
- 未正则full与未池化diagonal只作诊断；full奇异时直接标记不可评分，不加epsilon抢救。
- 全8alpha重算上述四个背景、平均／完整冻结响应、H2/H3、observed／matched-zero，共512个条件后验。另64个alpha=.3加密网格对照。

## 预测证据支持的结论

1. 原始经验full在留出时段NLL/d可达数十至约116000，而scalar约−2；部分折秩亏。支持撤回原始full的优先地位，而不是宣称真实背景无相关性。
2. H2两个24h差分中，时间小矩阵相对scalar无平均优势；部分池化异方差在1D反而每折更差，2D却改善。原统一sigma的稳定性动机有实际意义，不能全盘否定。
3. H3加入96h差分后，时间小矩阵相对一个scalar，在5折／3折均改善NLL/d约0.35；额外24h隔离后保留优势。当前模拟对应24h约1.9deg、96h约5.9–6.0deg，不能把空间池化的理由用来主张所有lag共享同一sigma。
4. 固定收缩的空间模型在5折、3折以及额外24h隔离中各fold均改善scalar的背景预测；可以保留为正则化候选，但不是lambda=.25唯一最优、5次独立统计显著或真实天空覆盖率认证。

因此新的工作建议是：H2保留统一sigma基线；H3以少参数时间结构为主要简单候选，并列给出固定空间收缩的敏感性结果。不再以未正则经验full作为默认升级。现有单条轨迹不足以断言唯一正确模型。

## alpha=.30：使用完整冻结响应后的实际条件上限

|背景|1D H2|1D H3|2D H2|2D H3|
|---|---:|---:|---:|---:|
|scalar|1.819|1.786|1.420|1.576|
|池化时间矩阵|1.746|1.262|1.399|1.589|
|0.5部分池化方差|1.844|1.743|1.517|1.488|
|0.75时间目标+0.25空间经验C|1.714|1.012|1.374|1.532|

这是比较模型，不按最小c95挑选。H3 scalar这次按三差分池化，区别于上一轮H2固定方差池，已记录：不能将该差异全部归因于增加数据。32项H2 scalar平均响应同旧研究逐项恢复，最大差0dex。

## 数值、独立复核与完成状态

- 纯GRMHD192个fold/model拟合记录；额外24h隔离192个记录。
- 新512个后验与64个dense后验内部一致性检查；alpha=.3加密变化最大.002946dex（约0.68%），不外推成所有alpha网格已独立收敛。
- 独立代码审查未发现所查路径的致命数据泄漏；未将其升级成完整谱系／科学正确性认证。
- 全回归1659passed227.23s；新增最终focused27passed；本次新增源Ruff E/F/W全部通过。
- 四个新完成stage全部声明产物SHA验证通过；运行所绑定Python源码与当前字节一致，证据在background-reassessment-evidence/FINAL_VERIFIED.json。
- 首次v1因CSV dimension字段重复在推断前失败，保留日志；v2完成，未抹去失败。
- EHT、RT和旧canonical NULL未重建、未覆盖；没有push或新增远程任务。

## 结果与复现

主计算：`/home/canoziia/GRMHD_DATA/analysis/eht-background-reassessment-v2`

汇总：`/home/canoziia/GRMHD_DATA/analysis/eht-background-reassessment-results-v1`

图：`cv-nll.png`、`raw-full-failure.png`、`alpha030-c95.png`、`allalpha-c95.png`。

完整说明：同目录`BACKGROUND_REASSESSMENT_ZH.md`以及本仓库对应报告；CSV与后验保留。

```bash
pixi run python scripts/reassess_background_study.py \
  --old-summary /home/canoziia/GRMHD_DATA/analysis/eht-joint-days-allalpha-v1/summary.csv \
  --output-root /new/unused/output
```

输出必须是未占用新目录。额外检查分别使用`check_background_embargo.py`、`check_background_reassessment_grid.py`；汇总用`summarize_background_reassessment.py`。没有引入包安装、运行notebook或独立依赖清单。
