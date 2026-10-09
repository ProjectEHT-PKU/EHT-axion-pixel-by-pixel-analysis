# 图像视场与 Gaussian 核：来源专项核查

**后续直接证据更新：**[挂载原运行目录的常数重算](../reports/script-pipeline/MOUNTED_CONSTANTS_COMPARISON.md)。
当前八cube/四pickle已与原目录SHA256匹配；FITS确认128µas视场，RAPTOR生产M6.2e9/D16.8、camera46.2261Rg。
本报告早期“实际FOV未确认”的结论已更新；20µas匹配核的科学依据仍未闭环，数值参数未改。

## 结论先行

本轮只核查来源，不修改配置、算法或科学结果。

当前三个数值 **不是互不相关的独立常数**：

```text
initial_fov_rg = 46.2262
final_fov_rg   = 32.6868
gaussian_sigma_rg = 2.16871
```

- **直接确认：**历史 RAPTOR HDF5 的图像坐标给出输入FOV=46.2261，与代码46.2262一致到约2.2ppm；该FOV表示正方形边长，不是半径。
- **强数值关联：**initial≈sqrt(2)×final；sigma≈final×(20/128)/2.355，相对误差不到1e-7。
- **文献确认：**EHT-IV的eht-imaging采用64²、128μas；EHT-VIII模拟比较采用20μas FWHM。EHT-VII偏振成像则采用64²、120μas，对eht-imaging图像匹配核为19μas。
- **尚未确认：**当前4800图pickle的实际角FOV、主RAPTOR cubes生成时的M/D/角标定。它们不能仅凭论文或数字吻合决定。
- **新发现：**保留的历史RAPTOR HDF5+uniform输出反推theta_g=3.643272μas/Rg，与旧RAPTOR常数下M6.2e9/D16.8完全吻合；并非128/32.6868隐含的3.915954μas/Rg。

因此：**可以高置信度恢复三个数值的候选设计关系，但还不能宣布整批实际输入已使用统一物理标定。**

## 1. 先区分三个物理量

设角引力半径为：

```text
theta_g = (GM/c²)/D          [μas / Rg]
```

则：

```text
FOV_angle = FOV_rg × theta_g
pixel_rg = FOV_rg / N
FWHM_angle = sqrt(8 ln2) × sigma_rg × theta_g
sigma_input_pixels = sigma_rg × N_input / initial_fov_rg
```

`scipy.ndimage.gaussian_filter`的第二参数是**标准差σ，以像素为单位**。
旧代码`ImageS.Filter`使用`2.16871 × N_input/46.2262`，故2.16871是**σ[Rg]**，不是20μas本身，也不是FWHM。

FOV控制覆盖范围，Gaussian控制分辨率；二者物理作用不同，但都需要同一theta_g换算为角尺度。
不能把“图像128μas宽”理解成“Gaussian128μas宽”。

## 2. 旧源码的直接证据

原始项目仍在`../EHT_ZC`，未修改。

### 2.1 三个数值共同进入同一处理过程

`EHT_ZC/GRMHD/Gaussian Blur.ipynb` 第2个cell（零基1，ID a6439ddc，In[25]）：

```python
self.res = (self.shape - 1) / 46.2262
xgrid = np.linspace(-23.1131, 23.1131, self.shape, endpoint=False)
gausian_kernal = 2.16871  # * 1.558722, 20uas = 2.16871
rotation_angle = 0.6 * np.pi
# 最终坐标网格 ±16.3434
```

后来的`EHT_ZC/utils/libanaly.py:360–386,411–417,440–470`将它们明确写成参数：

```python
init_fov=46.2262
final_fov=32.6868
final_shape=64
self.gaussian_kernel = 2.16871
self.res = self.shape / self._init_fov
gaussian_filter(Q0, self.gaussian_kernel * self.res, mode="nearest")
```

这能证明同一软件处理链和作者“20uas”的意图，不能证明原作者具体使用的M/D。
注释中的`*1.558722`没有实际执行，未找到来源，不能擅自纳入换算。

### 2.2 旧网格有N和N−1约定差异

早期Notebook还用`63×0.518838=32.686794`作为输出首末中心跨度，而之后标记坐标的
FOV=32.6868、endpoint=False应给中心跨度`32.6868×63/64`。早期输入也用(N−1)/FOV，
后来的正常blur分支使用N/FOV。

因此早期实现混有“中心到中心跨度”和“边界到边界宽度”；这不是换一个物理质量就能解释的差别。
当前新实现沿用后者的半像素中心约定；此次未改变它。

## 3. 两个FOV为什么差sqrt(2)

直接计算：

```text
46.2262 / 32.6868 = 1.414216136
sqrt(2)          = 1.414213562
相对差           = 1.82e-6
```

一个旋转后的正方形若要完整采样内侧目标正方形，所需外框边长为：

```text
L_input >= L_output × (|cos θ| + |sin θ|)
```

任意旋转角度的最坏值是sqrt(2)，因此这组比例**高度符合“放大输入视场，旋转后裁出最终图像”设计**。
当前θ=108°时最小几何比例约1.260074；sqrt(2)是更保守的余量，而不是108°专属精确值。

还需区分：

- 这解释的是几何覆盖关系，不自动确定最终角FOV是120还是128μas。
- 额外边框不保证Gaussian边界效应为零；旧代码用nearest边界延拓，必须另行评估边缘。
- 未找到原作者写明“因为旋转所以乘sqrt(2)”的推导文本，不能把强数值证据说成原始配置文件已经找到。

## 4. Gaussian与最终FOV的更精确关联

取128μas视场、20μas FWHM、旧常用近似FWHM=2.355σ：

```text
theta_g = 128 / 32.6868 = 3.915953841 μas/Rg
sigma_rg = 20 / (2.355 × theta_g) = 2.168710191
代码 sigma_rg                         = 2.168710000
相对差                               ≈ 8.8e-8
```

等价地：

```text
2.16871 × 2.355       = 5.10731205 Rg
32.6868 × 20 / 128    = 5.10731250 Rg
```

这是比“约10像素”更紧的关系。使用精确sqrt(8ln2)=2.354820045：

```text
FWHM / final_pixel = 9.999234979 pixels
若finalFOV=128μas，FWHM=19.998469958μas
```

所以**候选原设计链**很可能是：

```text
final_fov_rg = 128μas / theta_g
initial_fov_rg ≈ sqrt(2) × final_fov_rg
sigma_rg = 20μas / (2.355 × theta_g)
```

但在已查源码/草稿中，没有发现128/theta_g或2.355的原计算表达式。它是高度精确的反推，不是当前pickle FOV的直接证据。

## 5. EHT论文能明确给出什么

### EHT-IV（2019）

[First M87 EHT Results IV](https://arxiv.org/abs/1906.11241)，[PDF](https://arxiv.org/pdf/1906.11241)，§6.2.2，PDF p13：

> using a 64 × 64 pixel grid with a 128 μas FOV

因此64²/128μas有直接文献依据，像素尺度为2μas。
同段还说明eht-imaging重建最终有**5μas FWHM**卷积（与inverse taper配套）；
不能把本文中CLEAN恢复beam20μas误认为每一张eht-imaging图原生都经过同样的20μas卷积。

### EHT-VII（2021，偏振）

[First M87 EHT Results VII](https://arxiv.org/abs/2105.01169)，[PDF](https://arxiv.org/pdf/2105.01169)：

- Appendix D.1，PDF p22：
  > fix the image field of view at 120 μas and solve for a grid of 64 × 64 pixels
- §5.1，PDF p10：通过与2019共识强度图匹配，选择Gaussian FWHM：
  **eht-imaging/DMC/THEMIS=19μas，LPCAL=20μas，polsolve=23μas**。

120μas/64=1.875μas每像素。这里的19μas是附加的匹配卷积核，不应不加区分地叫作
每个原始重建图唯一精确PSF。当前pickle是否已经执行这种匹配blur也未知，需生成代码。

### EHT-VIII（2021，模拟比较）

[First M87 EHT Results VIII](https://arxiv.org/abs/2105.01173)，[PDF](https://arxiv.org/pdf/2105.01173)：

- §4.4，PDF p9：模拟偏振指标比较使用**20μas FWHM圆Gaussian**。
- PDF p8：库标定采用**M=6.2e9太阳质量、D=16.9Mpc**；模拟库原始图为320²/160μas。

这支持“用20μas圆Gaussian模拟EHT有限分辨率”的方法，但不能证明我们的RAPTOR数据采用该库全部参数。

## 6. 历史RAPTOR实际文件中的坐标证据

本地保留了：

```text
EHT_ZC/GRMHD/data_in/data001600.dat/img_data_0.h5
EHT_ZC/GRMHD/data_in/data001600.dat/uniform_img_2.30e+11_0.dat
```

HDF5中的alpha、beta各有100个唯一像素中心：

```text
中心范围：[-22.8819195, +22.8819195]
间隔：0.462261
边界FOV：100×0.462261 = 46.2261
```

HDF5没有质量、距离或单位属性；从公开RAPTOR源码的camera/输出约定，可解释为Rg坐标。
这份HDF5的I按坐标重排后，与历史`GRMHD/data_in/grmhdData1600.npy` state0 **逐元素完全相同**。
它与当前主cube α=.3、t1600、state0的I并不相同（差/峰值≈0.192），因此不能把它当成主cube同批的原文件。

### 角尺度从uniform输出反推

项目旧稿引用的RAPTOR源码版本b91e937...：

- [io.c](https://github.com/jordydavelaar/raptor/blob/b91e937b7b746eb24992b6c1dd8831c5dea0cef1/src/io.c#L317-L354)：
  `arcsec_factor=206265*R_GRAV/SOURCE_DIST`，uniform输出前两列为x/y×此因子。
- [camera.c](https://github.com/jordydavelaar/raptor/blob/b91e937b7b746eb24992b6c1dd8831c5dea0cef1/src/camera.c#L53-L77)：HDF5的alpha/beta是相机几何坐标。
- [definitions.h](https://github.com/jordydavelaar/raptor/blob/b91e937b7b746eb24992b6c1dd8831c5dea0cef1/model/bhac/definitions.h#L147-L157)：
  G=6.674e-8、MSUN=1.989e33、KPCTOCM=3.086e21、c=2.99792458e10。
- [core.c](https://github.com/jordydavelaar/raptor/blob/b91e937b7b746eb24992b6c1dd8831c5dea0cef1/src/core.c#L97-L98)：质量/距离转换。

将本地uniform坐标与HDF5坐标对应：

```text
theta_g = 3.6432720895784816 μas/Rg
原始角FOV ≈ 168.41426μas
```

用上述旧常数、M=6.2e9、D=16.8Mpc计算，theta_g与它**数值完全一致**。
这只确定M/D及相应常数约定，不能唯一确定M和D各自的值；也不能证明实际执行程序就是该Git版本。

## 7. 不同标定下的后果

保持当前三个代码参数不动：

| 角尺度依据 | theta_g (μas/Rg) | final角FOV (μas) | Gaussian FWHM (μas) |
|---|---:|---:|---:|
| 候选128μas设计 | 3.91595384 | 128.0000 | 19.9985 |
| 假定120μas图像 | 3.67120673 | 120.0000 | 18.7486 |
| 当前配置M6.5e9/D16.8，IAU单位 | 3.81899325 | 124.8307 | 19.5033 |
| 历史raw文件反推 | 3.64327209 | 119.0869 | 18.6059 |

因此绝对尺度不同，但代码的**beam/最终FOV比值基本固定为20/128**。
把质量改成6.5或距离改成16.8而不重新推导FOV/σ，不能声称最终图就自动变成120或128μas。

## 8. 实际EHT候选图仍缺什么

本轮读取四天pickle并逐个检查全部4800条metadata字段并集：

```text
MEM,T2V,TV,chi2,cphase,l1,logcamp,m,prior_fwhm,pvis,sys_noise,vis,zbl
```

四天均没有fov、psize、pixel_size、mass、distance字段；顶层只有meta_data和I/Q/U/theta数组。
其中prior_fwhm是成像先验宽度，**不是最终视场或恢复beam**，不能据其数值恢复像素角尺度。

现有EHT处理Notebook是pickle消费者，直接硬编码±16.3434或Rg/pixel=0.51073125；
未找到生成这些pickle的原重建脚本、FITS header或FOV设置。

必须取得至少一项：

1. 实际生成4800图的脚本中fov/npix/psize设置；
2. 同批原始FITS的CDELT/NAXIS等图像头；
3. 本次主RAPTOR运行的CAM_SIZE、M/D配置及HDF5坐标。

然后才适合决定用120、128或其他角FOV，以及19/20μas还是其他匹配核。

## 9. 建议（本轮未实施）

先固定“真实输入图像的角坐标”和“希望匹配的角分辨率”，再派生Rg参数：

```text
final_fov_rg = actual_eht_fov_uas / theta_g
sigma_rg = target_fwhm_uas / (sqrt(8ln2) × theta_g)
initial_fov_rg = actual_raptor_camera_fov_rg  # 优先直接读取坐标，不凭裁切关系猜值
```

FOV是输入数据几何，不能为满足公式随意改；裁切关系应做覆盖检查而非覆盖真实元数据。
原始camera若只有46.2261Rg，直接把initialFOV改成另一个数只是在错误地重标图像。
当前名义配置没有随本轮研究改动，结果也没有重算。
