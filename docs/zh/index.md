# 文档索引

## 用途

本页是 M87\* EHT 轴子分析中文文档的入口页。本仓库使用两个相互独立的输入：带轴子场的
GRMHD 吸积流经 RAPTOR 光线追踪得到的图像，以及 M87\* 的 EHT 2017 偏振图像。分析部分把二者
结合起来，给出无量纲耦合常数 $c$ 的贝叶斯上限。下列每一页都描述本仓库中代码的现状，并且每一页
都链接回本页。英文文档见 [English documentation](../en/index.md)。

## 页面

| 页面 | 内容 |
|---|---|
| [概览](overview.md) | 三个部分、各部分产品之间的关系，以及数据流图 |
| [实验配置](experiment-config.md) | `experiment.yaml` 的每个键：含义、单位、取值以及读取它的代码 |
| [RAPTOR 模拟](raptor-simulation.md) | 任务计划与任务标识符、编译、作业、作业输出、格式转换、交付与合并、CPU 亲和性、`.env.raptor` |
| [EHT 成像](eht-imaging.md) | 按实现描述的 Paper VII 方法：I 阶段、共享增益、基准 P/D 交替、精化、打包导出、D 项蒙特卡罗、系综汇总、已知偏差 |
| [分析](analysis.md) | 图像预处理、响应拟合、模板、孔径观测量、特征选择、GRMHD 背景、协方差模型、似然、后验、检验与汇总 |
| [来源记录](provenance.md) | 阶段清单（manifest）、哈希值、原子化阶段目录、校验 |
| [数据布局](data-layout.md) | 每个脚本写出的目录和文件，以及数组形状和单位 |
| [包结构](package-map.md) | `eht_axion` 的每个模块和每个脚本，各附一句用途说明 |
| [开发](development.md) | Pixi 环境与任务、测试、代码检查（lint）、`.env` 文件 |
| [术语表](glossary.md) | 这些页面中使用的符号和术语 |

## 从哪里开始

- 理解实验：先读[概览](overview.md)，再读[实验配置](experiment-config.md)。
- 审查某个科学步骤：读 [EHT 成像](eht-imaging.md)或[分析](analysis.md)，并结合[来源记录](provenance.md)。
- 运行脚本：读 [RAPTOR 模拟](raptor-simulation.md)、[EHT 成像](eht-imaging.md)和[分析](analysis.md)
  中的“运行方法”一节；环境说明见[开发](development.md)。
- 查找函数或文件：读[包结构](package-map.md)和[数据布局](data-layout.md)。

## 本文档的约定

- 文件名、键名、函数名和阶段名写成代码格式，例如 `analysis.features`。
- YAML 键带上所在节的路径，例如 `eht_imaging.refine.gap_atol`。
- 所有数组中的电矢量位置角（EVPA）均以弧度为单位；EVPA 以 $\pi$ 为模定义。EHT 产品中的天空
  角度从北向东量度。
- 天空上的长度以微角秒（µas）或引力半径 $r_g = GM/c^2$ 为单位；时间以小时或 $t_g = GM/c^3$
  为单位。这些定义中的光速 $c$ 不是耦合常数 $c$；见[术语表](glossary.md)。

## 相关页面

- [仓库 README](../../README.md)
- [English documentation](../en/index.md)
- [概览](overview.md)
- [术语表](glossary.md)
