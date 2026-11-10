# 包结构

## 用途

本页列出 `eht_axion` 包的每个模块和 `scripts/` 中的每个脚本，各附一句说明及其主要公共函数，使读者
能够从某个科学步骤找到实现它的代码。包的结构与脚本对应：`eht_axion/raptor` 服务于 `raptor_*`
脚本，`eht_axion/eht` 服务于 `eht_*` 脚本，`eht_axion/analysis` 服务于 `analysis_*` 脚本；被多个
部分使用的模块位于顶层。

## 细节

### 顶层模块

| 模块 | 用途 | 主要函数和类 |
|---|---|---|
| `eht_axion/experiment.py` | 载入并校验 `configs/experiment.yaml`；派生设置 | `read_experiment`、`validate`、`validate_ray_tracing`、`physical_scale`、`ray_tracing_config`、`analysis_geometry`、`eht_imaging_parameters`、`ANGULAR_BINS`、`EHT_DAYS` |
| `eht_axion/provenance.py` | 阶段清单、输入和输出的 SHA-256、原子化阶段目录、校验 | `stage_directory`、`verify_stage`、`require_declared_outputs`、`verification_session`、`sha256`、`manifest_path`、`source_closure`、`contained_path` |
| `eht_axion/physics.py` | 物理常数、M87\* 的单位尺度、轴子 211 态频率 | `PhysicalScale`、`omega` |
| `eht_axion/evpa.py` | 线偏振的 EVPA 和以 $\pi$ 为周期的差值，由 `eht/` 和 `analysis/` 共用 | `finite_evpa`、`circular_difference` |

### RAPTOR 模拟：eht_axion/raptor

| 模块 | 用途 | 主要函数 |
|---|---|---|
| `plan.py` | 只读的快照尾部和径向表检查；任务列表和任务标识符 | `plan_generation`、`task_identity`、`read_snapshot_footer`、`validate_radial_table` |
| `build.py` | 用官方 `setup.sh` 编译固定版本的 RAPTOR 核心程序 | `build_raptor` |
| `jobs.py` | 并发、逐任务原子化、可续跑的作业运行器 | `run_jobs`、`serialize_parameters`、`parameter_values` |
| `affinity.py` | 每个工作进程互不相交的 CPU 集合和内存余量检查 | `allocate_cpu_sets`、`check_memory`、`launch_prefix` |
| `outputs.py` | 校验一个作业的 HDF5 和谱文件 | `validate_job_output`、`read_core_constants`、`hdf_mapping` |
| `convert.py` | 把作业 HDF5 块转换为 float64 Stokes 数据立方；小型数组文件辅助函数 | `convert_hdf5` |
| `delivery.py` | 数据立方的交付与合并；按实验检查交付的数据立方 | `deliver_stokes`、`generated_inputs`、`input_cadence_hours`、`validate_processed_source` |
| `layout.py` | 交付数据立方的目录布局（`alpha-<α>/` 和 `background/`） | `alpha_label`、`product_directory`、`input_raptor_files` |

### EHT 成像：eht_axion/eht

| 模块 | 用途 | 主要函数 |
|---|---|---|
| `calibration.py` | 基于 ehtim 1.3.2 的 Paper VII 校准基本操作：原生数据、系统噪声、固定 D 项改正、零先验前向 D 拟合、四天 D 分布及抽样 | `load_native`、`add_fsys`、`make_imaging_obs`、`apply_fixed_d`、`fit_forward_d`、`d_distribution`、`sample_dterms` |
| `convex.py` | 固定 I、固定 D 的精确凸 P 求解器（PDHG），附对偶间隙证书 | `exact_cartesian_problem`、`solve_problem`、`CartesianProblem` |
| `ensemble.py` | 原生 FITS 读取器、Paper VII 图像平均量、按天的系综汇总 | `load_stokes_fits`、`global_diagnostics`、`summarize` |

### 分析：eht_axion/analysis

| 模块 | 用途 | 主要函数 |
|---|---|---|
| `templates.py` | 模拟图像预处理、响应拟合、模板汇总及读取、模板根目录布局 | `preprocess_frame`、`preprocess_cube`、`fit_phase_dft`、`fit_processed`、`summarize_fits`、`load_templates`、`RunLayout` |
| `aperture.py` | 强度加权的孔径 EVPA（1D 观测量）和扇区选择统计量，使用同一个求积算子 | `aperture_operator`、`integrate_stokes`、`sector_brightness` |
| `model.py` | 轴子 EVPA 信号模型 | `cloud_amplitude_scaling`、`interpolate_fit_parameters`、`evpa_templates` |
| `observations.py` | 从选定的 EHT 特征得到逐日矩和日对差 | `ObservationData`、`moments`、`build_observations`、`branch_risk`、`covariance_diagnostics` |
| `likelihood.py` | 推断选项、高斯似然、相位边缘化、上限 | `LimitOptions`、`MODES`、`calculate_observed` |
| `posterior.py` | 网格上的后验密度、累积分布函数（CDF）和上限 | `normalized_posterior_cdf`、`posterior_upper_limit` |

### 脚本

| 脚本 | 用途 |
|---|---|
| `raptor_build.py` | 把固定版本的 RAPTOR 核心程序编译到 `RUN_ROOT/build` |
| `raptor_generate.py` | 规划任务、运行作业并交付一个运行根目录的 Stokes 数据立方 |
| `raptor_merge.py` | 把多台主机作业根目录中的数据立方交付到一个新的根目录 |
| `eht_reconstruct.py` | EHT 成像阶段 `I` 和 `fiducial`，每个输出根目录一天 |
| `eht_refine.py` | 对基准阶段做精确凸 P 精化和 D 重拟合 |
| `eht_export_bundle.py` | 把四个 I 阶段和四个精化阶段打包为蒙特卡罗输入 |
| `eht_sample.py` | D 项蒙特卡罗：每次 D 抽样求一幅精确 P 图像 |
| `eht_summarize.py` | 蒙特卡罗图像的逐日系综汇总 |
| `analysis_templates.py` | 由交付的 RAPTOR 数据做预处理、响应拟合并生成模板 |
| `analysis_limits.py` | 特征选择、观测量、GRMHD 背景、似然和上限 |
| `analysis_check_grid.py` | 在另一个相位或耦合常数网格上重新计算保存的观测量 |
| `analysis_check_covariance.py` | 由抽样的两半得到的协方差敏感性 |
| `analysis_summarize.py` | 汇总表、比较和图 |
| `analysis_paper.py` | 论文的图、引用数字（`numbers.tex`、`numbers.json`）和表格 |
| `_paper_figures.py` | 论文图的样式和绘图函数 |
| `verify_products.py` | 校验一个已完成的阶段（`pixi run verify`） |
| `_bootstrap.py` | 把仓库根目录加入 `sys.path`，以便直接执行脚本 |
| `_environment.py` | 以数据方式读取 `.env` 和 `.env.raptor`，只允许一组固定的键 |
| `_raptor_options.py` | `raptor_build.py` 和 `raptor_generate.py` 共用的命令行选项 |
| `_data_root.py` | `--data-root` / `M87_AXION_DATA_ROOT` 以及数据布局中的默认根目录 |

### 测试

`tests/` 含有合成数据的单元测试和集成测试，按所测试的模块或脚本命名（例如 `test_raptor_plan.py`、
`test_eht_convex.py`、`test_likelihood.py`、`test_analysis_limits.py`）。`test_docs_links.py` 检查
本文档的链接。见[开发](development.md)。

## 相关页面

- [文档索引](index.md)
- [概览](overview.md)
- [数据布局](data-layout.md)
- [开发](development.md)
