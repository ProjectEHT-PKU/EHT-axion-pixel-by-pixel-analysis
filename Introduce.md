## GRMHD
141个吸积盘，时间间隔$1R_g$（$141*1R_g$），omega=0.4， axion_norm=0.1
  - 每个吸积盘跑了1个无轴子的图和16个不同相位有轴子的图
    - 141个无轴子的图的EVPA -> 天文delta EVPA uncertainty（所有满足条件的像素delta EVPA一起Gaussian Fit得到一个数）（条件：$2<r<9$，$P_{noaxion}$>0.1e9）

    - 同一时间16个不同轴子相位的图 -> cosine 拟合出A，φ（部分点会拟合失败）。选取某个稳定时的吸积盘（n=51），以该吸积盘的A和φ为准去计算轴子一天引起的EVPA变化的大小和uncertainty（pixel by pixel 但是无cov，对角阵）：
	  - 用该吸积盘的A和φ得到一天引起的EVPA变化
	  - 对141张图的A和φ高斯拟合，得到轴子效应的A与φ的uncertainty（部分点会拟合失败） -> 误差传递得到最后轴子引起的delta EVPA的uncertainty （还是用n=51的A，φ）
	  - 这里计算的变化和变化的std，都是基于omega=0.4，axion_norm=0.1。为了估计出不同omega和axion_norm的影响，依据photon ring处的轴子场强做放缩

	  最终得到不同omega下的每个pixel的理论预言的delta EVPA和其uncertainty

Note：
  - 黑洞、吸积盘角动量三维朝向确定
  - 拟合出来的A在部分pixel比axion norm得到的理论最大值要大：高斯模糊导致的




## EHT-Imaging
topset大小为1572，按照chi2排序最好的1572张，每张图可以都得到EVPA。对这1572张图再进行二次筛选，方法：

对于每个（$2<r<9$）的pixel，EVPA值都随不同图像有一个分布；对于每个图像，其中（$2<r<9$）的像素都离平均值有个相对偏差（$\Delta/\sigma$）；把一张图所有这些像素的相对偏差求平方和得到这张图的偏离指数，选取偏离指数最小的前一半。

这样就得到了大小为786的topset。此时把所有图像平均会得到一个黑洞图像，将其重定位到图像中心，并记住这个移动的大小。在这个基础上，可以高斯拟合求出新的每个pixel的EVPA平均值和std（部分点会拟合失败），以及不同像素间的关联。

Note：
  - 重新定位黑洞中心曾经有两种考虑，1. 对图像的平均重定位 2. 对每张图像先重定位再求平均。最后选择了1，因为EHT IV给出的每张图像并没有重定位
  - EVPA的定义的理解：在一张图的不同pixel中，EVPA的0度方向始终是一致沿x轴。

## Result
最终考虑的因素：uncertainty: GRMHD、Axion、Observation; Axion prediction; 
不考虑的因素：认为 
  - Observation=0
  - GRMHD prediction=0

Note:
  - 所有在过程中拟合失败的pixel都被丢弃了
  - 1d的选点有自由度。360度角由于积分时包括了周围10度所以有36个相隔10度的独立数据点，起始角度可以是0-9任意。

	

## 要画的图

Fig. 1 Radial Function of axion field for different ωs and rs.

Fig. 2 2D-plot of axion-signal: Amplitude, relative phase—-Fix omega

Fig. 3 1D-plot
of axion-signal: Amplitude, relative phase—-comparing ωs, comparing analytic results and
raptor’s result.

Fig. 4 2D-plot of uncertainty: GRMHD, axion-signal: Amplitude and phase.

Fig.5 Topset- distribution, Fiducial Image.

Fig.6 Observation-uncertainty-sigma:2D-plot. covariance?.

Fig.7 Constraints.