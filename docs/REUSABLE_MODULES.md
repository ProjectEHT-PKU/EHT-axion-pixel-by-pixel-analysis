# 可复用模块与研究脚本的边界

本项目不是可安装的软件包，不注册统一CLI；`eht_axion/`只是本地源码的组织名称。
这并不妨碍单独调用已经验证的函数。也没有要求每位使用者都执行整条研究流程。

## 1. 哪些内容适合单独复用

| 内容 | 源文件 | 示例函数 |
|---|---|---|
| 周期角度统计 | `eht_axion/angles.py` | `circular_difference`、`circular_mean`、`unwrap_around` |
| HDF5显式Stokes转换 | `eht_axion/pipeline/raptor.py` | `convert_hdf5` |
| 图像模糊/旋转/重采样 | 同上 | `preprocess_frame`、`preprocess_stokes` |
| 相位余弦拟合 | 同上 | `fit_phase_dft` |
| 孔径积分 | `eht_axion/polarization.py` | `aperture_operator`、`integrate_stokes` |
| 概率与后验数值 | `eht_axion/statistics.py` | `gaussian_log_likelihood`、`posterior_upper_limit` |
| 缓存/校验 | `eht_axion/pipeline/common.py` | `stage_directory`、`verify_stage` |
| 结果图 | `eht_axion/pipeline/figures.py` | `plot_results` |

这些函数不包含命令行参数解析。它们的输入shape、单位、坐标和假设仍须遵守文档。
例如默认角度周期是π（EVPA），普通相位要显式指定2π；HDF5必须明确数据集名；
绘图读取我们的具名限制schema，并不是适配任意数据的通用绘图库。

`workflow.py`、`inference.py`以及配置中的topset/配准/噪声模型，则是我们研究专用的逻辑，
不宣称已经形成普适标准流程。

## 2. 在代码目录中单独调用

```bash
pixi shell
python
```

```python
import numpy as np
from eht_axion.angles import circular_difference

# EVPA周期为π；89°减去-89°应取最近的周期差-2°。
difference = circular_difference(np.deg2rad(89), np.deg2rad(-89))
print(np.rad2deg(difference))
```

没有pip install；Python从当前源码目录读取这些模块。
在其他目录，不会自动全局找到本项目，这是预期行为，不是环境安装失败。

## 3. 写自己的脚本

如果脚本放在本项目`scripts/`中，可以这样写：

```python
import _bootstrap  # 将相邻源码目录加入本脚本的搜索路径，不安装项目
import numpy as np
from eht_axion.angles import circular_mean

angles = np.deg2rad([89, -89])
print(np.rad2deg(circular_mean(angles)))
```

然后：

```bash
# 已进入pixi shell后：
python scripts/my_analysis.py
```

脚本通过自身路径寻找源码，所以从其他cwd用绝对脚本路径运行也可以。
若复用到另一个项目，可以复制经过确认的模块及对应测试，或显式管理源码路径；
不需要为此把整个M87分析发布为可安装库。

`pixi shell`激活环境后可直接运行Python脚本；Pixi任务名本身不会变成shell命令。
如果不进入交互环境，也可以继续使用`pixi run`任务。

## 4. 各研究脚本

| Pixi任务 | 直接执行的文件 |
|---|---|
| `pipeline` | `scripts/run_all.py` |
| `raptor` | `scripts/process_raptor.py` |
| `eht` | `scripts/process_eht.py` |
| `background` | `scripts/estimate_background.py` |
| `observations` | `scripts/build_observations.py` |
| `inference` | `scripts/calculate_limits.py` |
| `figures` | `scripts/plot_limits.py` |
| `diagnostics` | `scripts/plot_diagnostics.py` |
| `convergence` | `scripts/check_convergence.py` |
| `verify` | `scripts/verify_products.py` |
| `migrate` | `scripts/migrate_data.py` |
| `transport` | `scripts/transport_diagnostics.py` |
| `convert-hdf5` | `scripts/raptor_tools.py convert` |

`raptor_tools.py`、`process_eht_ensemble.py`和`plot_saved_limits.py`是单独处理已有数据的
低级工具，不负责全run的缓存调度；主研究步骤使用上表带cache的脚本。
参数解析留在脚本，复用模块不依赖它。没有`eht-pipeline`console命令、`python -m
 eht_axion.pipeline`统一入口、editable安装或项目distribution元数据。

依赖、测试搜索路径及代码检查规则全部保存在`pixi.toml`中。
不再保留`pyproject.toml`、`pytest.ini`或`ruff.toml`。请使用`pixi run test`和
`pixi run lint`，确保应用这里统一声明的参数。
