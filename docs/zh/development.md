# 开发指南

[← 返回文档索引](index.md)

本仓库由公共模块目录（`eht_axion/`）、直接执行的命令行入口脚本（`scripts/`）以及合成数据测试套件（`tests/`）组成，通过 [Pixi](https://pixi.prefix.dev/)（`pixi.toml` 与 `pixi.lock`）实现可复现的软件环境管理。本仓库不打包为 wheel，也不安装进 `site-packages`。

## Pixi 环境与任务

全部四个 Pixi 环境在激活时都会统一设置环境变量 `MPLBACKEND=Agg`、`OPENBLAS_NUM_THREADS=2`、`OMP_NUM_THREADS=2` 和 `PYTHONUNBUFFERED=1`。

| 环境名 | 支持平台 | 在 `default` 基础上新增的依赖 | 用途 |
|---|---|---|---|
| `default` | `linux-64`、`linux-aarch64`、`osx-64`、`osx-arm64`、`win-64` | Python `3.12.*`、`numpy`、`scipy`、`matplotlib`、`h5py`、`pyyaml`、`pytest`、`ruff`、`python-dotenv` | 分析脚本、代码检查（`lint`）及不依赖 `ehtim`/`astropy` 的测试 |
| `raptor-build` | `linux-64` | `c-compiler`、`make`、`gsl`、`hdf5`、`pkg-config`、`rsync` | 编译 `RAPTOR`、运行光线追踪作业、写出图像数据立方 |
| `eht-summary` | 全部 5 个平台 | `astropy` | 将仪器偏振泄漏（D 项，D-term）蒙特卡罗（Monte Carlo）FITS 抽样汇总为 `eht:ensemble` |
| `eht-imaging` | 全部 5 个平台 | `astropy`、`ehtim == 1.3.2`（PyPI） | EHT 重建、凸优化精化、蒙特卡罗抽样以及完整测试套件 |

| 任务名 | 环境 | 执行命令 / 目标脚本 |
|---|---|---|
| `raptor-build`、`raptor-generate` | `-e raptor-build` | `scripts/raptor/{build,generate}.py` |
| `eht-reconstruct`、`eht-refine`、`eht-sample` | `-e eht-imaging` | `scripts/eht/{reconstruct,refine,sample}.py` |
| `eht-summarize` | `-e eht-summary`（或 `-e eht-imaging`） | `scripts/eht/summarize.py` |
| `analysis-expected-signal`、`analysis-limits`、`analysis-paper` | `default` | `scripts/analysis/{expected_signal,limits,paper}.py` |
| `test` | `default` 或 `-e eht-imaging` | `python -m pytest -q -ra -o "testpaths=tests" -o "pythonpath=."` |
| `lint` | `default` | `ruff check --isolated --target-version py312 --line-length 110 --select E,F,W eht_axion tests scripts` |
| `check` | `default` | 依次运行 `test` 与 `lint` |

## 测试与代码检查

在审查或提交代码改动前，请在 `eht-imaging` 环境中运行完整测试套件（含 `ehtim` 与 `astropy` 相关测试）并运行隔离配置的 Ruff 代码检查工具（linter）：

```bash
OPENBLAS_NUM_THREADS=2 pixi run --locked -e eht-imaging test
pixi run --locked lint
```

在 `default` 环境中（`pixi run --locked check`），导入 `ehtim` 或 `astropy` 的测试会通过 `pytest.importorskip` 自动跳过。所有测试仅在 `pytest` 的 `tmp_path` 临时目录中使用合成数据运行，不读取观测 `uvfits` 文件、快照或 `RAPTOR` 可执行文件。

除常规单元测试外，以下测试模块用于强制检查全仓库级别的架构约束：
- `tests/test_experiment.py`：检验 `configs/experiment.yaml` 的结构与取值范围校验，并扫描 `eht_axion/` 的抽象语法树（AST，`test_library_signatures_have_no_numeric_default_values`），禁止公共函数带有数值型默认参数（仅保留五项并发与数组轴白名单）。
- `tests/test_data_root.py`：检验 `eht_axion.environment.LAYOUT` 的每个条目均出现在 `docs/en/data-layout.md` 与 `docs/zh/data-layout.md` 的 `<data root>/` 目录树中、`INPUT_LAYOUT` 的每个条目均出现在 `<input root>/` 目录树中，并验证全部脚本的命令行与 dotenv 路径优先级。
- `tests/test_docs_links.py`：检验 `docs/en/` 与 `docs/zh/` 具有完全相同的文件集合与标题层级序列，每页均可由 `README.md` 经相对链接到达并链接回对应语言的 `index.md`，且所有相对链接与 `#anchors` 锚点均可解析。
- `tests/test_stage_safety.py` 与 `tests/test_provenance.py`：检验符号链接拒绝、输入/输出路径重叠防护、脏工作树拒绝以及异常时的原子回滚。

## 代码结构与模块索引

### 顶层公共模块（`eht_axion/`）

| 模块 | 职责说明 | 主要公共接口 |
|---|---|---|
| `eht_axion/experiment.py` | `configs/experiment.yaml` 的严格加载器、模式/范围校验器及派生配置视图 | `read_experiment`、`validate`、`validate_ray_tracing`、`physical_scale`、`ray_tracing_config`、`analysis_geometry`、`eht_imaging_parameters`、`ANGULAR_BINS`、`EHT_DAYS` |
| `eht_axion/environment.py` | 唯一 `.env` 的安全解析器、三个根目录（`EHT_RAPTOR_SOURCE_ROOT`、`EHT_INPUT_ROOT`、`EHT_DATA_ROOT`）的解析，以及输入与产物目录布局 | `resolve_environment`、`add_environment_options`、`add_data_root_option`、`add_input_root_option`、`data_root`、`input_root`、`raptor_source_root`、`layout_path`、`input_path`、`fill_defaults`、`fill_input_defaults`、`LAYOUT`、`INPUT_LAYOUT` |
| `eht_axion/provenance.py` | 原子化阶段目录、`manifest.json`、Git 提交身份及输入文件 SHA-256 | `stage_directory`、`verify_stage`、`require_declared_outputs`、`sha256`、`manifest_path`、`source_closure`、`contained_path`、`code_identity` |
| `eht_axion/physics.py` | SI/IAU 物理常数、M87\* 引力时间与角度尺度、211 束缚态频率 | `PhysicalScale`、`omega` |
| `eht_axion/evpa.py` | 线偏振电矢量位置角（EVPA）与 $\pi$ 周期圆统计角度差 | `finite_evpa`、`circular_difference` |

### 子包模块

| 模块 | 职责说明 | 主要公共接口 |
|---|---|---|
| `eht_axion/raptor/plan.py` | BHAC 快照尾部与径向 CSV 校验、任务列表及按参数命名的任务 ID | `plan_generation`、`task_id`、`read_snapshot_footer`、`validate_radial_table` |
| `eht_axion/raptor/build.py` | 洁净源码工作副本核验与固定版本 `RAPTOR` 的官方 `setup.sh` 编译 | `build_raptor` |
| `eht_axion/raptor/jobs.py` | 带确定性工作槽位与断点续跑的有界并发 `RAPTOR` 作业运行器 | `run_jobs`、`serialize_parameters`、`parameter_values` |
| `eht_axion/raptor/affinity.py` | 互不相交的 `taskset` CPU 分配及 Linux/`cgroup-v2` 内存余量检查 | `allocate_cpu_sets`、`check_memory`、`launch_prefix` |
| `eht_axion/raptor/outputs.py` | 单作业 HDF5 像素块/坐标及积分流量谱校验 | `validate_job_output`、`read_core_constants`、`hdf_mapping` |
| `eht_axion/raptor/convert.py` | 将分块作业 HDF5 拼接为 float64 `(T, S, H, W, 4)` `stokes.npy` | `convert_hdf5` |
| `eht_axion/raptor/images.py` | 由作业 HDF5 写出 Stokes 图像数据立方及下游实验配置身份核对 | `write_images`、`generated_inputs`、`input_cadence_hours` |
| `eht_axion/raptor/layout.py` | 图像数据立方的目录命名规则（`alpha-<α>/`、`background/`） | `alpha_label`、`product_directory`、`image_cube` |
| `eht_axion/raptor/options.py` | `RAPTOR` 脚本共用的命令行/dotenv 解析器 | `options` |
| `eht_axion/eht/calibration.py` | 原生圆偏振 HOPS 载入、`fsys`、逆琼斯 D 改正、零先验前向 D 拟合与四天 D 抽样 | `load_native`、`add_fsys`、`make_imaging_obs`、`apply_fixed_d`、`fit_forward_d`、`d_distribution`、`sample_dterms` |
| `eht_axion/eht/convex.py` | 带对偶间隙（duality gap）证书的固定 $I$/固定 $D$ 精确凸优化原始–对偶混合梯度（PDHG）$P$ 求解器 | `exact_cartesian_problem`、`solve_problem`、`CartesianProblem` |
| `eht_axion/eht/ensemble.py` | 严格四扩展单元（4-HDU）Stokes FITS 读取器、Paper VII 全局诊断量与逐像素 EVPA 图 | `load_stokes_fits`、`global_diagnostics`、`global_values`、`map_evpa` |
| `eht_axion/analysis/results.py` | 上限类阶段的结果文件：数组存入 `results.npz`，其余存入 `report.json` | `write_results`、`read_results` |
| `eht_axion/analysis/expected_signal.py` | 高斯模糊与天空旋转重采样、逐快照 16 相位离散傅里叶拟合、跨快照极大似然统计及读取 | `preprocess_frame`、`fit_phase_dft`、`snapshot_response`、`alpha_responses`、`expected_signal`、`processed_background`、`load_expected_signal` |
| `eht_axion/analysis/aperture.py` | 用于 1D 扇区 EVPA 和扇区亮度的压缩稀疏行（CSR）双线性求积算子 | `aperture_operator`、`integrate_stokes`、`sector_brightness` |
| `eht_axion/analysis/model.py` | 轴子云振幅标度 $b(\alpha)$、$\alpha$ 节点圆插值及双历元日对 EVPA 信号模型 | `cloud_amplitude_scaling`、`interpolate_fit_parameters`、`evpa_templates` |
| `eht_axion/analysis/observations.py` | 逐日圆统计矩、总体 D 散布协方差、日对差值以及分支割线/特征值诊断 | `ObservationData`、`moments`、`build_observations`、`branch_risk`、`covariance_diagnostics` |
| `eht_axion/analysis/likelihood.py` | 三种协方差模式（`std`、`cor`、`cor_fixed_background`）、高斯 $\chi^2$ 与相位边缘化 | `LimitOptions`、`MODES`、`calculate_observed` |
| `eht_axion/analysis/posterior.py` | 梯形法则后验归一化与分段线性累积分布函数（CDF）解析求逆上限 | `normalized_posterior_cdf`、`posterior_upper_limit` |
| `eht_axion/analysis/paper_figures.py` | 论文全部 8 幅确定性矢量 PDF 图的样式与绘制函数 | `save`、`radial_function`、`two_photon`、`signal_maps`、`template_maps`、`background`、`observational_uncertainty`、`data_vs_model`、`limits` |

### 入口脚本（`scripts/`）

每个入口脚本开头都带有两行 `sys.path` 引导代码，定位最近的含 `eht_axion/` 的上级目录作为仓库根目录，因此可在任意当前工作目录下直接运行 `python scripts/<part>/<name>.py`。测试文件命名与模块和脚本一一对应（例如 `scripts/<part>/<name>.py` 对应 `tests/test_<part>_<name>.py`）。

## 扩展流水线

### 架构不变量

1. **所有科学参数在 `configs/experiment.yaml` 中只写一次**，在 `eht_axion/experiment.py` 中严格校验，并以无默认值的显式参数传给 `eht_axion/` 函数。
2. **机器路径与并发开关绝不写入 YAML**：统一通过命令行选项、`add_data_root_option` / `fill_defaults` 或在 `eht_axion/environment.py`（`PATH_KEYS`、`VALUE_KEYS`）中注册的环境变量传入。
3. **只读取已校验且在清单中声明的阶段输出**：对上游阶段调用 `verify_stage` 与 `require_declared_outputs`；写出结果时置于 `stage_directory(...)` 内，将上游 `manifest_path(dir)` 及实际读取的数据文件列入 `input_paths`（仓库源码文件不列入）。
4. **保持每行不超过 110 字符**，且无任何 Ruff `E`、`F`、`W` 规则告警（目标版本 `py312`）。

### 添加新的科学参数

1. 在 `configs/experiment.yaml` 对应的节（`shared`、`raptor`、`eht_imaging` 或 `analysis`）下添加新键及其取值。
2. 在 `eht_axion/experiment.py` 的 `SCHEMA` 中登记该键名，并在 `_validate_shared`、`_validate_raptor`、`_validate_eht_imaging` 或 `_validate_analysis` 中补齐类型与区间检查。
3. 若该参数决定 `RAPTOR` 光线追踪结果，确保 `ray_tracing_config`（若影响单任务还需纳入任务参数）纳入该键；若影响 EHT 成像阶段，确保 `eht_imaging_parameters`（`_EHT_PARAMETER_ORDER`）或求解器计划字典纳入该键。
4. 从入口脚本把该取值显式传入 `eht_axion/` 函数签名（不设数值默认值），并在 `tests/test_experiment.py` 中补充合法与非法取值的测试用例。
5. 在 `docs/en/configuration.md` 与 `docs/zh/configuration.md` 中同步记录该参数。

### 添加新的分析阶段

1. 在 `eht_axion/environment.py` 的 `LAYOUT` 中新增布局键（如 `"analysis.my_check": "analysis/my-check"`），并在 `docs/en/data-layout.md` 与 `docs/zh/data-layout.md` 的 `<data root>/` 目录树 `analysis/` 下同步添加对应行（由 `tests/test_data_root.py` 强制检查）。
2. 在 `eht_axion/analysis/` 中实现不含文件 I/O、不带科学默认值的纯数值计算逻辑。
3. 新建 `scripts/analysis/my_check.py`，使用标准两行 `sys.path` 引导、`add_data_root_option` 和 `fill_defaults`；用 `verify_stage` + `require_declared_outputs` 校验上游阶段，拒绝已存在的 `--output-root`，并在 `stage_directory(destination, "analysis:my-check", identity, inputs)` 中写出文件。
4. 在 `pixi.toml` 中注册 `analysis-my-check` 任务，在 `tests/test_analysis_my_check.py` 与 `tests/test_data_root.py` 中添加合成数据测试，并在 `docs/en/` 与 `docs/zh/` 中补充文档。
