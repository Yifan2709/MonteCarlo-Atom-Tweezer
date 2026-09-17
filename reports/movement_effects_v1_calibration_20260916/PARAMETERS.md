# 四组实际参数和物理约束

| 参数 | 修正前 | 200 us | 400 us | 600 us | ML400 us | 单位 | 来源 | 操作 / 范围 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T0 harmonic proposal | 15 | 19 | 19 | 19 | 19 | uK | 物理工作域内未测标定量 | 联合调整；4.3–28；微波阶段温度锚点与T<=U/5工作域，非实验置信区间 |
| SLM/AOD nominal waist | 1.17 | 1.11 | 1.11 | 1.11 | 1.11 | um | SLM论文测量1.17(6)；AOD匹配入瞳假设 | 共同束腰联合调整；1.11–1.23；采用报告括号不确定范围，未拟合独立AOD腰 |
| AOD RIN PSD | 1e-08 | 5.75e-09 | 5.75e-09 | 5.75e-09 | 5.75e-09 | 1/Hz | 弱白噪声有效模型，尚未装置测量 | 联合调整；0–1e-8；0–100kHz积分rms<=3.2%弱扰动工作域 |
| SLM RIN PSD | 1e-08 | 0 | 0 | 0 | 0 | 1/Hz | 论文132s静态寿命约束；短时安静近似 | 结构修正，非拟合；固定0；另1e-12敏感性 |
| Tref | 15 | 不使用 | 不使用 | 不使用 | 不使用 | uK | 旧线性化单模参考能量 | 在物理扩散模型中删除；只在历史/系数修正对照保留15 |
| Heating reference frequency | 36007 | 不使用固定频率 | 不使用固定频率 | 不使用固定频率 | 不使用固定频率 | Hz | 原配置显式采用标称AOD径向解析频率 | 局部力自动确定响应；旧固定功率仍用36007；新模型不再整段套同一个谐振频率 |
| Constant total heating | 485.631567 | 不使用固定功率 | 不使用固定功率 | 不使用固定功率 | 不使用固定功率 | kB*uK/s | 旧单模参量+指向+固定反冲功率 | 由局部力扩散替代；历史值保留，修正系数对照的精确值见运行JSON |
| Recoil scattering | 4.23 | 随局部强度计算 | 随局部强度计算 | 随局部强度计算 | 随局部强度计算 | 1/s | 旧SLM+AOD两能级估计；非直接测量 | 改为局部反冲扩散；新模型采用Cs D2两能级近似，另关闭对照 |
| SLM depth | 140 | 140 | 140 | 140 | 140 | uK | 论文Fig6 Methods直接报告 | 固定；未报告可信校准误差，未擅自拟合 |
| AOD target depth | 280 | 280 | 280 | 280 | 280 | uK | 论文Fig6 Methods直接报告 | 固定；ML矢量回放终点按论文280；原回放280.5167另存 |
| Transport distance | 2.4 | 2.4 | 2.4 | 2.4 | 2.4 | um | 论文直接报告 | 固定；ML原始4nm偏置修正来自PDF控制点，非存活拟合 |
| Wait | 100 | 100 | 100 | 100 | 100 | us | 论文直接报告 | 固定；100 us；具体SLM关断边沿未报告 |
| SLM wavelength | 1061 | 1061 | 1061 | 1061 | 1061 | nm | 论文直接报告 | 固定；Gaussian Rayleigh range |
| AOD wavelength | 1055 | 1055 | 1055 | 1055 | 1055 | nm | 论文直接报告 | 固定；Gaussian Rayleigh range |
| Cs mass | 132.90545196 | 132.90545196 | 132.90545196 | 132.90545196 | 132.90545196 | u | 原配置Cs-133质量 | 固定；不移植Rb/Yb |
| SLM depth relative sigma | 0.114 | 0.114 | 0.114 | 0.114 | 0.114 | fraction | 论文全阵列实测；子阵列未单独给出 | 固定；未逐组调散布 |
| AOD depth relative sigma | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | fraction | 无装置测量假设 | 不拟合；零散布对照；无可用测量范围 |
| AOD waist relative sigma | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | fraction | 无装置测量假设 | 不拟合；零散布对照 |
| Alignment sigma | 50 | 50.0 | 50.0 | 50.0 | 50.0 | nm | 无装置测量假设 | 不拟合；固定及关闭敏感性 |
| Segment alignment jitter sigma | 10 | 10.0 | 10.0 | 10.0 | 10.0 | nm | 无带宽依据的原模型假设 | 不拟合；固定及关闭敏感性；不等同连续PSD |
| AOD pointing PSD | 1e-24 | 1e-24 | 1e-24 | 1e-24 | 1e-24 | m²/Hz | 原安静光路敏感性假设 | 不拟合；固定与关闭对照；只沿运输x轴 |
| SLM pointing PSD | 1e-24 | 0 | 0 | 0 | 0 | m²/Hz | 短时安静近似 | 结构修正，非拟合；0；另1e-26敏感性 |
| Lensing effective vs | 0.5391689206547424 | 0 | 0 | 0 | 0 | m/s | 旧severity情景；论文SI与器件参数支持额外0.412425情景 | 模型情景比较，不连续拟合；0在代码中表示关闭，物理上vs趋于无穷；非零须测光束填充率 |
| Parametric coefficient | pi²/4 | pi² | pi² | pi² | pi² | dimensionless | Savard Eq12及本次独立推导 | 确认修正；PSD仍为单边每Hz，定义未更改 |
| Pointing coefficient | 1/8 | 1/4 | 1/4 | 1/4 | 1/4 | dimensionless | Savard Eq15谱转换及本次独立推导 | 确认修正；同一PSD定义 |
| Heating implementation | P=Gamma kB Tref | Gaussian force diffusion | Gaussian force diffusion | Gaussian force diffusion | Gaussian force diffusion | model | Ito力噪声推导，谐振与实际核基准通过 | 物理结构修正；没有额外运动加热；SLM/AOD噪声分别规定但四组共用 |
| Preparation | site-bound harmonic | site-bound harmonic | site-bound harmonic | site-bound harmonic | site-bound harmonic | model | 已修复的原初始化，保持 | 不调整；每个实际位点E<0；没有逐组截断/热尾/预筛选 |
| Observation | E>=0 absorbing | E>=0 absorbing | E>=0 absorbing | E>=0 absorbing | E>=0 absorbing | model | 实验末次SLM成像的近似 | 未拟合阈值；等待末AOD，放下末SLM；另允许再捕获对照 |
| Manual timing | 48/52 | 48/52 | 48/52 | 48/52 | 48/52 | percent | 论文直接报告 | 固定；二次升深/单段三次移动 |
| ML interpolation | raster CubicSpline | PDF-marker CubicSpline | PDF-marker CubicSpline | PDF-marker CubicSpline | PDF-marker CubicSpline | model | 独立PDF几何，未使用存活率 | 波形来源修正；未夹紧导数、未单调平滑、位置/深度同钟 |

操作时长分别为200/400/600/400 us，每30往返总时长为15/27/39/27 ms；这两行是实验允许的差异。主要影响与退化关系详见CSV的effect和degeneracy列。所有数值均为实际传给积分器的共同配置，Tref字段仅在历史对照使用。
