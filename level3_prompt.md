# Codex 执行 Prompt：Level 3 三维 AOD 长距离运输与柱面透镜效应

请在当前工作区中，在已经完成并验证的 Level 0 静态势阱、Level 1 一维顺序转移和 Level 2 同步转移波形优化基础上，实际完成“Level 3：三维 AOD 长距离运输与柱面透镜效应”。

不要只给方案或代码片段。必须检查并复用已有项目，创建必要模块、配置、测试、数据和图表，实际运行预检查、参数扫描和独立蒙特卡洛验证，并诚实报告所有模型假设、未通过项和无法由论文唯一确定的硬件参数。

如果工作区存在 `sources/`，其中所有内容均为只读，不得编辑、移动、重命名或删除；如果不存在，直接继续。任何输出路径解析后都不得位于 `sources/`。

论文依据：arXiv:2403.12021v4，尤其是 Fig. 5、Extended Data Fig. 10、Methods “Atom transport” 和 Supplementary Information 第 III 节。论文报告 280 μK AOD 光镊中 610 μm、1.6 ms 的对角运输，并指出直线运输受 AOD 柱面透镜效应限制；论文还比较了 constant-jerk 与 adiabatic-sine 轨迹。官方版本：https://arxiv.org/abs/2403.12021v4

## 一、Level 3 的定位

Level 2 已经解决短距离 2.4 μm 的 AOD-SLM drop-off 波形问题。Level 3 暂不把 transfer 与 long move 组合，而是隔离研究“AOD 单阱已经装有原子后，如何长距离移动”。

核心问题：

> 对一个初始处于 280 μK AOD 光镊中的有限温度铯原子，在 270-610 μm 的运输中，运输方向、轨迹形状、运输时间和 AOD 柱面透镜强度如何影响经典留阱率、末态激发和轴向运动？

Level 3 必须完成：

1. 把现有一维模型扩展为两个横向坐标加一个光轴坐标的三维经典模型；
2. 在零扫频极限恢复标准三维高斯光镊及 Level 0 的径向结果；
3. 实现论文 SI 启发的 crossed-AOD 柱面透镜势；
4. 比较直线运输与对角运输；
5. 比较 piecewise constant-jerk、adiabatic-sine 和 minimum-jerk/control 轨迹；
6. 研究 270、510 和 610 μm 距离及多个运输时长；
7. 使用热初态蒙特卡洛统计经典留阱率和加热；
8. 识别直线运动中瞬时势阱由单阱变为双阱的临界区；
9. 检查显式时变三维势中的功-能关系与时间步长收敛；
10. 使用无量纲透镜强度扫描，避免伪造论文未提供的唯一硬件校准；
11. 将模型结果与论文实验锚点作定性比较，但不得声称复现量子相干保真度。

本 Level 3 的“survival/retention”仅表示经典末态仍被 AOD 势束缚。它不包括真空损失、光子散射、量子态退相干、XY4、随机基准测试或成像误差，因此不得与论文的 99.953(2)% 瞬时量子通道保真度直接等同。

## 二、坐标与静态三维高斯势

使用以下坐标命名，避免与论文图中的实验室轴混淆：

- \(x_1,x_2\)：焦平面内两条正交的 AOD 横向轴；
- \(z\)：光轴方向；
- \(\boldsymbol c(t)=[c_1(t),c_2(t)]\)：移动光镊在焦平面的中心；
- \(\xi=x_1-c_1(t)\)，\(\eta=x_2-c_2(t)\)。

默认忽略重力。如果以后启用重力，必须明确实际竖直轴和势能 \(mgq\)，不得在轴定义不清时加入。

波长和 Rayleigh 长度：

\[
\lambda=1055\ \mathrm{nm},
\qquad
z_R=\frac{\pi w_0^2}{\lambda}.
\]

零扫频、无透镜效应时，三维高斯势应化为

\[
U_0(\xi,\eta,z)
=-\frac{D}{1+(z/z_R)^2}
\exp\!\left[
-\frac{2(\xi^2+\eta^2)}{w_0^2[1+(z/z_R)^2]}
\right].
\]

中心的小振幅频率为

\[
\omega_r=\sqrt{\frac{4D}{mw_0^2}},
\qquad
\omega_z=\sqrt{\frac{2D}{mz_R^2}}.
\]

默认使用：

- \(^{133}\mathrm{Cs}\) 精确同位素质量；
- AOD 深度 \(D=k_B\times280\ \mu\mathrm K\)；
- \(w_0=1.17\ \mu\mathrm m\)，标注为与论文 SLM 后孔径匹配所作的 AOD 假设；
- \(\lambda=1055\ \mathrm{nm}\)；
- 初态温度 5 μK；
- 转移全程命令深度默认恒定为 280 μK；
- 运输结束后保持 100 μs，检查末态留阱与静态能量误差。

## 三、论文 SI 启发的 crossed-AOD 透镜势

快速扫动某条 AOD 轴时，该轴对应的柱面焦点沿光轴偏移。定义两个高斯因子：

\[
G_1(\xi,z;z_{s,1})=
\frac{1}{\sqrt{1+[(z-z_{s,1})/z_R]^2}}
\exp\!\left[
-\frac{2\xi^2}{w_0^2\{1+[(z-z_{s,1})/z_R]^2\}}
\right],
\]

\[
G_2(\eta,z;z_{s,2})=
\frac{1}{\sqrt{1+[(z-z_{s,2})/z_R]^2}}
\exp\!\left[
-\frac{2\eta^2}{w_0^2\{1+[(z-z_{s,2})/z_R]^2\}}
\right].
\]

三维移动势为

\[
U_{\rm AOD}(x_1,x_2,z,t)
=-D_{\rm cmd}(t)G_1G_2.
\]

当 \(z_{s,1}=z_{s,2}=0\) 时，必须严格恢复上一节的标准三维高斯势。

对两条 AOD 轴分别采用

\[
\frac{z_{s,j}}{z_R}=2\,\sigma_j\frac{\dot c_j}{v_s},
\qquad j\in\{1,2\},
\]

其中 \(v_s\) 是单轴 characteristic shift velocity，\(\sigma_j=\pm1\) 表示具体 AOD 安装和衍射方向。默认选择与论文“对角运动把柱面透镜变为球面透镜”相符的同号方向，并将符号写入配置。

论文 SI 给出近似关系

\[
v_s\approx C\frac{w_0}{T_a},
\]

其中 \(T_a\) 是声波穿过输入光束的 access time，简单 Airy 估计给出 \(C\approx5.3\)，论文装置表现更接近 7-8。由于论文没有提供足够数据在本项目中唯一确定 \(T_a\) 和有效 \(C\)，不得硬编码一个值并称为“论文参数”。

必须支持三种校准模式：

1. `lensing_off`：\(z_{s,1}=z_{s,2}=0\)，作为理想模型；
2. `absolute_critical_velocity`：用户明确给出 \(v_s\) 及来源；
3. `normalized_reference_scan`：以一个固定参考运输的单轴峰值速度为基准，扫描 \(v_{\max,\mathrm{axis}}/v_s\)。这是默认模式。

默认参考运输：610 μm 的对角 adiabatic-sine 运输，时长 1.6 ms。默认扫描：

```text
reference_axis_vmax_over_vs = [0.0, 0.50, 0.75, 1.00, 1.25]
```

对每个非零比值先由参考运输计算一个固定的 \(v_s\)，然后将这个相同 \(v_s\) 用于所有距离、方向、时长和轨迹比较。不得为每个候选重新缩放 \(v_s\)，否则比较不再代表同一硬件。

### 单阱/双阱诊断

对直线单轴运动，在横向中心处扫描瞬时轴向势 \(U(z)\)：

- 当 \(|z_{s,1}|<2z_R\) 时通常为一个展宽的局部最低点；
- 当 \(|z_{s,1}|>2z_R\) 时可能分裂为两个轴向最低点；
- 临界条件等价于单轴速度达到 \(v_s\)。

不得只依赖上述解析判断。必须用数值极值搜索检查每个时间点的局部最低点数量、位置、深度和势垒，并用解析临界关系作回归测试。

对同号、等速的理想对角运输，两条柱面焦点共同移动，模型应保持单一球面对称焦点，而不是产生直线运输式的双阱。必须通过数值极值测试验证这一点。

## 四、运输几何

将配置中的 `distance_um` 定义为焦平面中的欧氏路径长度 \(L\)。

### 直线运输

\[
c_1(t)=Lq(t),
\qquad
c_2(t)=0.
\]

只有第一条 AOD 轴扫频，因此只产生 \(z_{s,1}\)。

### 对角运输

\[
c_1(t)=\frac{L}{\sqrt2}q(t),
\qquad
c_2(t)=\frac{L}{\sqrt2}q(t).
\]

两条 AOD 的单轴速度均为路径速度的 \(1/\sqrt2\)。报告速度时必须区分路径速度和单轴 AOD 速度。

路径必须满足：起点、终点准确；端点速度为零；除明确的 piecewise jerk 跳点外，位置、速度和加速度连续；运动方向单调；无过冲。

## 五、必须实现的轨迹族

所有轨迹以 \(s=t/T\in[0,1]\) 和归一化路径 \(q(s)\in[0,1]\) 表示，并返回 \(q,\dot q,\ddot q\)，必要时返回 jerk。

### A. Adiabatic-sine

使用与论文 Methods 等价的归一化形式：

\[
q_{\rm sine}(s)
=s-\frac{\sin(2\pi s)}{2\pi}.
\]

验证其端点位置和速度，并检查在简谐、无透镜、小激发极限下，加热趋势与

\[
\Delta N\propto\frac{L^2}{\omega^5T^6}
\]

一致。只检查缩放趋势，不要求未知比例常数。

### B. Piecewise constant-jerk S-curve

使用四个等时段、jerk 符号依次为 \(+J,-J,-J,+J\) 的对称轨迹，从零速度和零加速度开始并结束。对总距离 \(L\)、总时长 \(T\)，采用

\[
J=\frac{32L}{T^3}.
\]

分段解析积分得到位置、速度和加速度，不得把 quintic minimum-jerk 误称为 constant jerk。检查其简谐加热趋势与

\[
\Delta N\propto\frac{L^2}{\omega^3T^4}
\]

一致。

### C. Minimum-jerk/control

使用

\[
q_{\rm minjerk}(s)=10s^3-15s^4+6s^5.
\]

它作为常见平滑控制轨迹，但必须正确命名，不得声称这是论文的 constant-jerk 轨迹。

三种轨迹在相同 \(L,T\) 下必须使用同一初态样本比较。

## 六、深度模式与补偿边界

默认 `commanded_depth_mode: constant`，即 \(D_{\rm cmd}=280\ \mu\mathrm K\)。这只表示命令光功率恒定；透镜模型中的瞬时有效最低点深度仍会随速度变化。

必须实现并区分：

- `ideal_constant_depth`：关闭透镜效应，标准移动高斯阱；
- `lensing_uncompensated`：命令深度恒定，允许有效阱深改变；
- `exploratory_lensing_compensation`：可选探索模式，根据数值瞬时最低点深度提高命令深度，但设定最大功率倍率。

探索补偿不是论文已经实现的 lensing compensation，不得如此表述。论文说明其沿位置进行了 RF 功率校准，同时仍把更完整的 AOD lensing compensation 留作改进。

补偿模式必须输出所需最大功率倍率、饱和时间比例和实际有效阱深；若达到功率上限，不得继续声称“保持恒定阱深”。

## 七、初态三维热采样

默认温度 5 μK。为保持与前级一致，使用三维简谐近似提议分布：

\[
x_1,x_2\sim\mathcal N(0,\sigma_r^2),
\qquad
z\sim\mathcal N(0,\sigma_z^2),
\]

\[
v_1,v_2,v_z\sim\mathcal N(0,k_BT/m),
\]

其中

\[
\sigma_r^2=\frac{k_BT}{m\omega_r^2},
\qquad
\sigma_z^2=\frac{k_BT}{m\omega_z^2}.
\]

用完整静态三维高斯势检查初态总能量

\[
E_i=\frac12m(v_1^2+v_2^2+v_z^2)+U_0(x_1,x_2,z)<0,
\]

只接受束缚样本。保存提议数、接受率、六个相空间分量的均值/协方差和随机种子。

必须明确：这是简谐提议分布加完整势束缚拒绝，不是完整有限深高斯阱中严格的三维正则系综。不同运输条件使用 common random numbers。

至少分离：

- `scan_pool`：粗参数扫描；
- `validation_pool`：预先声明的关键条件独立验证；
- `convergence_pool`：步长和模型回归；

三套种子不得相同。看到 validation 结果后不得改变预先声明的主条件；若修改，必须生成新的 validation 种子并记录原因。

## 八、三维推进与显式时变功-能核算

扩展现有时变 velocity-Verlet 为三维向量形式：

\[
\boldsymbol r_{n+1}=\boldsymbol r_n+\boldsymbol v_n\Delta t
+\frac12\boldsymbol a_n\Delta t^2,
\]

\[
\boldsymbol v_{n+1}=\boldsymbol v_n
+\frac12(\boldsymbol a_n+\boldsymbol a_{n+1})\Delta t.
\]

对显式时变势，它是二阶推进方法；不得仅因为静态 velocity-Verlet 是辛的，就声称三维显式时变实现具有普通自治系统的能量守恒。

总机械能：

\[
E(t)=\frac12m|\boldsymbol v|^2+U_{\rm AOD}(\boldsymbol r,t).
\]

应满足

\[
\frac{\mathrm dE}{\mathrm dt}
=\left.\frac{\partial U_{\rm AOD}}{\partial t}\right|_{\boldsymbol r}.
\]

时间依赖至少来自 \(c_1,c_2,z_{s,1},z_{s,2}\) 和可能的 \(D_{\rm cmd}\)。实现解析链式导数或经过严格收敛验证的固定坐标数值偏导，积分得到外部功，并检查

\[
R_W(t)=E(t)-E(0)-W_{\rm ext}(t).
\]

至少将外部功分解为：横向中心移动项、焦点轴向偏移项、命令深度项。若数值偏导无法稳定分解，可以将总功作为主要验收量，但必须报告分解误差而不能伪造精度。

## 九、末态留阱与加热定义

运输结束时 \(\dot c_1=\dot c_2=0\)，故 \(z_{s,1}=z_{s,2}=0\)。随后保持 100 μs，势阱保持静止。

在运输终点的静态 AOD 势中定义

\[
E_f=\frac12m|\boldsymbol v_f|^2
+U_0[x_1-c_1(T),x_2-c_2(T),z].
\]

严格判据：

\[
E_f<0\Rightarrow\text{classically retained},
\qquad
E_f\ge0\Rightarrow\text{not retained}.
\]

还要在 100 μs 保持结束时再次判定。主留阱率使用预先指定的“运输结束判据”；保持结束结果作为一致性诊断。两者不一致时必须保存样本并解释数值或定义原因。

每个 shot 至少记录：

- 初始和末态总能量；
- 初始/末态相对阱底激发 \(\varepsilon=E+D\)；
- \(\Delta\varepsilon\)；
- 径向和轴向最大位移；
- 最大速度；
- 是否进入瞬时双阱时段；
- 双阱期间跟随的轴向分支及是否发生分支切换；
- 总外部功和功-能残差；
- 运输结束与保持结束的留阱标签；
- 是否接近 \(E_f=0\) 分类阈值。

“进入双阱”和“分支切换”必须给出明确、可测试的数值定义。无法可靠分支分类的 shots 标为 `branch_ambiguous`，不得强行归类。

## 十、参数扫描与独立验证

采用分阶段计算预算，避免在巨大笛卡尔积上盲目运行高精度大样本。

### 阶段 A：确定性势阱与轨迹诊断

不使用 MC，扫描：

- 轨迹：adiabatic-sine、constant-jerk、minimum-jerk；
- 几何：straight、diagonal；
- 距离：270、510、610 μm；
- 时长：0.4、0.6、0.8、1.0、1.2、1.6、2.0 ms；
- 参考 lensing severity：0、0.50、0.75、1.00、1.25。

保存峰值路径速度、峰值单轴速度、最大 \(|z_s|/z_R\)、瞬时最低点深度、最低点数、双阱时间比例、最大加速度和 jerk。

### 阶段 B：粗蒙特卡洛相图

默认每个格点 128 shots、步长 0.10 μs。至少覆盖：

- straight 510 μm；
- diagonal 610 μm；
- adiabatic-sine 与 constant-jerk；
- 时长 0.6、0.8、1.0、1.2、1.6、2.0 ms；
- lensing severity 0、0.75、1.00、1.25。

使用完全相同的 `scan_pool` 比较所有格点。输出留阱率、Wilson 95% 区间、末态激发中位数和功-能残差。粗扫描用于找趋势，不得用 128 shots 的小差异声称高显著性。

### 阶段 C：预先声明的高精度验证

在运行 validation 前冻结以下主条件：

1. `paper_anchor_diagonal`：610 μm、1.6 ms、280 μK、diagonal、adiabatic-sine；
2. 同一条件的 constant-jerk；
3. `straight_control`：510 μm、1.6 ms、280 μK、straight、adiabatic-sine；
4. 同一条件的 constant-jerk；
5. 每个条件分别运行 ideal lensing-off 和一个预先选择的非零 severity；默认选择参考比值 1.00，若改动必须在打开 validation 前写入配置。

每个主条件至少 2,000 shots、步长 0.05 μs，并使用同一个独立 `validation_pool`。报告：

- 留阱率与 Wilson 95% 区间；
- 同初态配对 bootstrap 的留阱率差；
- 配对成败四格计数；
- 末态激发、轴向最大偏移和双阱经历；
- near-threshold 比例；
- 功-能残差；
- 0.025 μs 子集复算结果。

论文 610 μm、1.6 ms 只作为距离/时间锚点。不得把本模型得到的留阱率拟合或描述成论文的相干运输保真度。

### 阶段 D：敏感性

冻结主模型后，至少检查：

- 温度：3、5、10 μK；
- 深度：160、280、920 μK；
- \(w_0\)：1.17 μm 及其 ±0.06 μm；
- lensing severity：0.50、0.75、1.00、1.25；
- AOD lensing signs：同号与反号，反号只作模型敏感性，不冒充论文配置。

默认每个敏感性条件 500 shots。采用 one-factor-at-a-time 为主；若做联合最坏情形，必须单独标记。

## 十一、建议配置

在不破坏已有配置系统的前提下新增 Level 3 配置，例如：

```yaml
level3:
  name: 3d_aod_transport_lensing
  model_dimension: 3

atom:
  isotope: Cs-133
  mass_u: 132.90545196

aod_trap:
  wavelength_nm: 1055.0
  waist_um: 1.17
  depth_uK: 280.0
  waist_source: assumed_matched_to_slm_back_aperture

initial_ensemble:
  temperature_uK: 5.0
  sampler: harmonic_3d_rejection
  scan_shots: 128
  validation_shots: 2000
  sensitivity_shots: 500
  scan_seed: 31001
  validation_seed: 31002
  convergence_seed: 31003

integration:
  scan_dt_us: 0.10
  validation_dt_us: 0.05
  convergence_dt_us: 0.025
  post_transport_hold_us: 100.0

transport:
  geometries: [straight, diagonal]
  profiles: [adiabatic_sine, constant_jerk, minimum_jerk]
  distances_um: [270.0, 510.0, 610.0]
  durations_us: [400.0, 600.0, 800.0, 1000.0, 1200.0, 1600.0, 2000.0]

lensing:
  calibration_mode: normalized_reference_scan
  axis_signs: [1, 1]
  reference_geometry: diagonal
  reference_profile: adiabatic_sine
  reference_distance_um: 610.0
  reference_duration_us: 1600.0
  reference_axis_vmax_over_vs: [0.0, 0.50, 0.75, 1.00, 1.25]
  absolute_critical_velocity_m_per_s: null
  access_time_us: null
  access_time_speed_factor: null

depth_control:
  commanded_depth_mode: constant
  enable_exploratory_lensing_compensation: false
  max_power_multiplier: 2.0

classification:
  retained_energy_threshold_J: 0.0
  near_threshold_fraction_of_depth: 0.001
  escape_extent_waist: 20.0
  escape_extent_rayleigh: 20.0

output:
  directory: outputs/level3_3d_transport_lensing/demo
  save_full_validation_shots: true
  save_representative_trajectories: true
  representative_shots_per_class: 5
```

`escape_extent_*` 仅用于检测积分域是否过小和数值溢出，不得默认把越界立即当成物理损失，除非另有经过验证的吸收边界模型。主留阱判据仍是静态末能量 \(E_f<0\)。

## 十二、代码结构与复用

先检查当前仓库并沿用已有打包方式。必须复用现有常数、单位换算、随机种子管理、统计、序列化、图片检查和通用积分器逻辑；不得复制一套独立实现绕开 Level 0-2 回归。

建议职责，可与现有文件合理合并：

```text
gaussian_3d.py             # 标准三维高斯势、力和频率
aod_lensing.py             # G1/G2、焦点偏移、极值与分裂诊断
transport_profiles.py      # sine、constant-jerk、minimum-jerk
integrators_3d.py          # 向量化三维时变推进
thermal_sampling_3d.py     # 三维初态采样和拒绝
transport_work_energy.py   # 外部功和残差
transport_statistics.py    # Wilson、配对 bootstrap、四格计数
level3_simulation.py       # 扫描与 validation 编排
level3_visualization.py
level3_cli.py
```

所有公开函数包含简短中文 docstring。三维力应使用解析导数或高效稳定实现，并通过有限差分/复步长梯度测试。不得在每个时间步对每个 shot 调用通用标量优化器；瞬时极值诊断与轨迹推进应分离。

提供清晰 CLI，例如：

```bash
python -m level0_static_trap.level3_cli \
  --config configs/level3_3d_transport.yaml \
  --stage all \
  --output-dir outputs/level3_3d_transport_lensing/demo
```

若现有包名不同，使用实际包名。`--stage` 至少支持：

```text
preflight
potential_scan
coarse_mc
validate
sensitivity
all
```

不得使用临时 `PYTHONPATH` 掩盖打包问题。必须实际验证 editable/标准安装、现有测试和模块入口。

## 十三、预检查：失败则停止大规模 MC

按顺序运行并写入 `preflight.json`：

1. 运行全部 Level 0-2 测试；
2. 验证三维标准高斯势在 \(z=0,x_2=0\) 切片恢复 Level 0 一维势和力；
3. 验证静态三维数值曲率恢复 \(\omega_r,\omega_z\)；
4. 验证 \(z_{s,1}=z_{s,2}=0\) 时 lensing 势等于标准三维高斯势；
5. 验证三维解析力与势梯度一致；
6. 验证三类路径的端点、速度、加速度、距离和单调性；
7. 验证 constant-jerk 的四段解析积分和 \(J=32L/T^3\)；
8. 验证直线运动的数值单阱/双阱临界点接近 \(|z_s|=2z_R\)；
9. 验证同号等速对角模型保持单一轴向最低点；
10. 静态三维积分器能量误差有界并随步长二阶收敛；
11. 理想无透镜平移的功-能残差随步长收敛；
12. lensing 模型的功-能残差随步长收敛；
13. 理想中心、零速度初态完成 610 μm/1.6 ms 对角运输；
14. 对至少 128 个固定热初态做 0.10/0.05/0.025 μs 配对收敛；
15. 100 μs 静态保持段的能量误差有界；
16. 检查所有配置值有限、单位正确、主 validation 条件已在运行前冻结。

任何关键检查失败时停止大规模扫描，输出失败原因。不得跳过失败预检查后继续生成“完整结果”。

## 十四、输出文件

默认至少生成：

```text
outputs/level3_3d_transport_lensing/demo/
├── config_used.yaml
├── preflight.json
├── potential_minima_scan.csv
├── trajectory_profiles.csv
├── coarse_scan.csv
├── validation_shots.csv
├── sensitivity.csv
├── representative_trajectories.csv
├── metrics.json
├── image_validation.json
├── summary.md
├── static_3d_trap_slices.png
├── lensing_potential_bifurcation.png
├── trajectory_profile_comparison.png
├── straight_vs_diagonal.png
├── survival_phase_diagram.png
├── heating_comparison.png
├── axial_dynamics.png
├── final_energy_distributions.png
├── work_energy_balance.png
├── timestep_convergence.png
└── sensitivity_panels.png
```

### `potential_minima_scan.csv`

至少包含：geometry、axis velocities、\(v_j/v_s\)、\(z_{s,j}/z_R\)、时间、最低点数量、各最低点位置/深度、势垒、是否 split、有效阱深和补偿倍率。

### `trajectory_profiles.csv`

至少包含：profile、geometry、distance、duration、time、path fraction、centers、velocities、accelerations、jerks、焦点偏移和命令深度。

### `coarse_scan.csv`

每个参数格点一行，包含完整条件、shot 数、留阱率、Wilson 区间、near-threshold 数、激发统计、轴向偏移、split 时间比例、功-能残差和运行状态。

### `validation_shots.csv`

每个 shot 和每个冻结主条件一行，至少包含：共同初态 ID、六维初态、条件、末态六维状态、初末能量、retained 标签、near-threshold、激发变化、最大横向/轴向偏移、split/branch 指标、外部功、功-能残差及保持段结果。

### `metrics.json`

至少包含：

- 代码/配置/随机种子信息；
- 全部预检查；
- 物理常数、\(z_R\)、解析频率和模型假设；
- lensing 校准模式、由参考比值生成的每个固定 \(v_s\)；
- 所有冻结 validation 条件；
- 各主条件留阱率与 Wilson 区间；
- straight/diagonal、sine/constant-jerk 的配对差异和区间；
- 加热、轴向运动、split/branch 统计；
- 功-能残差和步长收敛；
- 温度、深度、束腰和 lensing severity 敏感性；
- 与论文可比较和不可比较的量；
- 图片程序化检查和视觉检查状态。

JSON 中不得出现 NaN 或 Infinity。

## 十五、诊断图

1. `static_3d_trap_slices.png`：无透镜静态势的横向和轴向切片，标出数值/解析曲率。
2. `lensing_potential_bifurcation.png`：直线模型在多个 \(z_s/z_R\) 下的轴向势，显示单阱到双阱；并给出对角模型对照。
3. `trajectory_profile_comparison.png`：三类 \(q,\dot q,\ddot q\) 和 jerk。
4. `straight_vs_diagonal.png`：相同物理路径长度/时长下，两种几何的单轴速度、焦点偏移、最低点位置和有效阱深。
5. `survival_phase_diagram.png`：留阱率随时长与 lensing severity 的相图，straight/diagonal 分面。
6. `heating_comparison.png`：仅 retained shots 的 \(\Delta\varepsilon\) 和末态激发分布；未留阱样本另行表示。
7. `axial_dynamics.png`：代表性 retained、lost、split、branch-ambiguous 轨迹的 \(z(t)\) 与瞬时最低点/分支。
8. `final_energy_distributions.png`：末态静态能量分布及 \(E_f=0\) 阈值。
9. `work_energy_balance.png`：\(\Delta E,W_{\rm ext},R_W\)；理想与 lensing、直线与对角均至少一例。
10. `timestep_convergence.png`：连续量误差、功-能残差和标签不一致率随步长。
11. `sensitivity_panels.png`：温度、深度、束腰和 lensing severity 敏感性。

使用非交互后端，至少 200 dpi，不调用 `plt.show()`，保存后关闭 figure。

## 十六、图片和数据验收

对每张 PNG 做程序化检查并写入 `image_validation.json`：存在、可重新打开、尺寸非零、不是全白/全黑/近乎单色、数据和坐标有限、图例数量合理、figure 已关闭。

- 若环境提供图像查看工具，逐张打开并记录 `visual_review_performed: true`；
- 若没有，记录 `false`，明确说明只完成程序化检查；
- 未实际看图时不得声称布局清晰或坐标轴合理。

对 CSV 检查 schema、行数、参数格点完整性、共同初态 ID 对齐和有限值。对 JSON/YAML/Markdown 重新读取验证。不得只检查文件存在。

## 十七、测试要求

保留全部 Level 0-2 测试，并至少新增：

1. 标准三维高斯势的中心值、对称性、远场极限；
2. 三维力与势梯度一致；
3. 静态径向/轴向曲率和频率正确；
4. lensing 零偏移极限严格恢复标准势；
5. 直线/对角中心和速度分解正确，欧氏路径长度均为 \(L\)；
6. adiabatic-sine 解析导数正确；
7. piecewise constant-jerk 分段连续、端点正确、总距离正确；
8. minimum-jerk 命名和公式正确；
9. 直线临界速度附近的最低点数量变化；
10. 对角同号模型无直线式双阱分裂；
11. lensing sign 反转按配置生效；
12. 三维热采样协方差、束缚拒绝和固定种子复现；
13. 三维静态积分能量误差有界；
14. 无透镜移动势的功-能关系；
15. lensing 势的功-能关系；
16. 步长减半的二阶收敛；
17. \(E_f<0\) 留阱、\(E_f=0\) 不留阱；
18. near-threshold 样本不被删除；
19. 主条件之间 common random numbers 对齐；
20. scan/validation/convergence 种子隔离；
21. validation 条件在运行前冻结；
22. Wilson 和配对 bootstrap 在构造数据上正确；
23. CLI 各 stage smoke test；
24. 严格 JSON、CSV schema、图片程序化检查；
25. 输出位于 `sources/` 时被拒绝。

测试不得要求模型必须得到论文的实验保真度，也不得要求 diagonal 必须在所有无透镜条件下优于 straight。应测试方程、几何、统计和数据隔离是否正确。

## 十八、执行顺序

1. 检查现有 Level 0-2 代码、配置、测试和输出；
2. 运行全部已有测试并记录基线；
3. 扩展三维势、力、轨迹、采样和推进；
4. 实现 lensing 模型、固定 \(v_s\) 场景生成和最低点诊断；
5. 实现功-能核算和统计；
6. 编写 Level 3 测试；
7. 运行 preflight，失败则修复并重跑；
8. 运行确定性 potential/trajectory scan；
9. 运行粗 MC 相图；
10. 在打开 validation 前冻结主条件和非零 severity；
11. 运行独立 2,000-shot validation；
12. 冻结主模型后运行敏感性；
13. 生成全部数据、报告和图片；
14. 程序化检查所有产物；
15. 若有图像查看工具，逐图视觉检查并修复展示问题；
16. 重新运行完整测试和冻结主条件，确保结果对应最新代码；
17. 检查没有修改 `sources/`；
18. 给出可核查的最终回复。

## 十九、最终回复必须包含

- 创建/扩展的主要模块和实际成功命令；
- 测试通过、失败和跳过数量；
- preflight 状态；
- \(z_R,\omega_r,\omega_z\) 的数值；
- lensing 校准模式及所有使用的固定 \(v_s\)，明确哪些是无量纲场景而非论文测量；
- 粗扫描覆盖范围和失败格点；
- 冻结主条件的留阱率、Wilson 区间和配对差异；
- straight/diagonal 与 sine/constant-jerk 的趋势；
- 单阱/双阱、轴向运动和加热统计；
- 功-能残差和步长收敛；
- 温度、深度、束腰、lensing severity 敏感性；
- 哪些结果只能定性对照论文，为什么不能等同量子保真度；
- 图片程序化检查及是否真正视觉检查；
- 输出目录和关键文件绝对路径；
- 建议 Level 4 才加入的 pick-up—长距离运输—drop-off 组合序列、重复往返和量子内态退相干。

不得只汇报“最佳留阱率”，不得隐藏不利的 straight、constant-jerk 或高 lensing severity 结果，也不得用实验保真度为当前经典模型背书。
