# Codex 执行 Prompt：Level 2 同步波形优化与鲁棒蒙特卡洛验证

请在当前工作区中，在已经完成并验证的 Level 0 静态势阱框架和 Level 1 顺序式 AOD→SLM 转移蒙特卡洛模拟基础上，实际完成“Level 2：同步移动-降深波形优化与鲁棒验证”。

不要只给方案或代码片段。必须检查现有项目、复用 Level 0/Level 1 的公共实现、创建必要代码与配置、运行测试和模拟、保存中间数据与最终结果，并诚实报告所有未通过项。

如果工作区存在 `sources/`，将其中内容视为只读，不得编辑、移动、重命名或删除；如果不存在，直接继续，不要创建空目录。任何输出路径解析后都不得落入 `sources/`。

## 一、Level 2 的定位

Level 1 已经实现如下顺序式 drop-off 基线：

1. 一个有限温度的铯原子初始被 280 μK AOD 光镊束缚；
2. 静止的 140 μK SLM 光镊位于目标位置；
3. 总时长 600 μs；
4. 前 52% 时间将 AOD 中心从距离 SLM 2.4 μm 处移动到 SLM 中心；
5. 后 48% 时间将 AOD 深度按二次曲线从 280 μK 降到 0；
6. 用时变 velocity-Verlet 推进；
7. 用末态 SLM 单阱能量 (E_f<0) 判定是否被 SLM 俘获；
8. 对 5 μK 热初态进行蒙特卡洛统计。

Level 2 不立刻进入二维、量子动力学或长距离运输，而是回答以下更聚焦的问题：

> 在相同的一维经典模型和初态分布下，允许 AOD 的位置与深度同时变化，能否把单次 drop-off 缩短到约 400 μs，同时保持或提高 SLM 俘获率，并降低末态激发？

Level 2 必须完成：

- 原样复现 Level 1 顺序波形，作为不可篡改的基线；
- 实现位置与深度可以重叠变化的受约束波形族；
- 使用固定训练样本优化波形；
- 使用从未参与优化的独立蒙特卡洛样本进行最终验证；
- 比较 600 μs Level 1 基线、压缩到 400 μs 的顺序基线和 400 μs 优化同步波形；
- 检查显式时变势中的功-能关系，而不是错误要求总机械能守恒；
- 进行时间步长、初态温度和末端对准误差的鲁棒性检查；
- 输出可复现数据、统计区间、诊断图和中文报告。

优化后未优于 Level 1 不等于实现失败。只要优化流程、独立验证和数值验收正确，应如实报告“未观察到显著提升”，不得重用验证集调参、删除不利样本、修改俘获判据或无依据放宽容差。

## 二、物理模型与坐标约定

继续使用 Level 0/Level 1 已验证的一维焦平面径向高斯势。令 SLM 中心为实验室坐标原点：

\[
x_{\rm SLM}=0.
\]

SLM 势为

\[
U_{\rm SLM}(x)
=-D_{\rm SLM}\exp\!\left[-\frac{2x^2}{w_{\rm SLM}^2}\right],
\qquad D_{\rm SLM}=k_B\times140\ \mu{\rm K}.
\]

AOD 势为

\[
U_{\rm AOD}(x,t)
=-D_{\rm AOD}(t)
\exp\!\left[-\frac{2[x-c_{\rm AOD}(t)]^2}{w_{\rm AOD}^2}\right].
\]

总势与总力为

\[
U_{\rm tot}(x,t)=U_{\rm SLM}(x)+U_{\rm AOD}(x,t),
\qquad
F_{\rm tot}(x,t)=-\frac{\partial U_{\rm tot}}{\partial x}.
\]

默认参数：

- 原子：(^{133}\mathrm{Cs})，质量使用现有项目中的精确同位素质量；
- (D_{\rm AOD}(0)=k_B\times280\ \mu\mathrm K)；
- (D_{\rm SLM}=k_B\times140\ \mu\mathrm K)；
- (w_{\rm SLM}=1.17\ \mu\mathrm m)；
- (w_{\rm AOD}=1.17\ \mu\mathrm m)，继续标记为“与 SLM 后孔径匹配所作的假设”，不是论文单独测量值；
- (c_{\rm AOD}(0)=2.4\ \mu\mathrm m)；
- (c_{\rm AOD}(T)=0)；
- (D_{\rm AOD}(T)=0)；
- 初态温度默认 5 μK；
- 转移结束后增加 100 μs 的纯 SLM 保持段，用于检查末态判据和静态能量守恒。

本 Level 2 仍然不模拟：二维/轴向运动、AOD 柱面透镜效应、光子散射、技术噪声、量子内态、相干性、加热与损失的随机过程、多个原子的相互作用、反向 pick-up、连续多次转移或 610 μm 长距离运输。

## 三、显式时变势中的功-能检查

转移期间机械能不守恒，因为外部控制改变了 AOD 的中心和深度。不得把 (E(t)) 的变化当成积分器“能量漂移”。

定义

\[
E_{\rm tot}(t)=\frac12mv^2+U_{\rm SLM}(x)+U_{\rm AOD}(x,t).
\]

沿轨迹应满足

\[
\frac{\mathrm dE_{\rm tot}}{\mathrm dt}
=\left.\frac{\partial U_{\rm AOD}}{\partial t}\right|_x.
\]

若将 AOD 深度记为 (D(t)>0)、中心记为 (c(t))、(q=x-c(t))，则

\[
\frac{\partial U_{\rm AOD}}{\partial t}
=-\dot D(t)e^{-2q^2/w^2}
+F_{\rm AOD}(x,t)\dot c(t).
\]

数值计算外部功

\[
W_{\rm ext}(t)=\int_0^t
\frac{\partial U_{\rm AOD}}{\partial t'}\,\mathrm dt',
\]

并检查

\[
R_W(t)=E_{\rm tot}(t)-E_{\rm tot}(0)-W_{\rm ext}(t).
\]

必须分别记录深度变化功、中心移动功、总外部功和功-能残差。转移后的 100 μs 保持段中 AOD 已关闭，此时 SLM 单阱能量应恢复为有界振荡误差而无 secular drift。

## 四、初态采样与公平比较

为保持与 Level 1 可比，默认严格复用 Level 1 的初态采样器：

1. 在 AOD 中心附近按简谐近似热分布采样 (x_0,v_0)；
2. 默认温度 5 μK；
3. 按完整 AOD 高斯势的 (E_{\rm AOD}(x_0,v_0)<0) 做拒绝；
4. 保存接受率、实际样本均值、标准差和随机种子；
5. 明确说明这是“简谐提议分布加束缚拒绝”，不是完整高斯阱中的严格正则系综。

比较不同波形时必须使用 common random numbers：同一组候选波形使用完全相同的初态数组。不得为每个候选重新抽样，否则蒙特卡洛噪声会掩盖波形差异。

至少设置三套互不重叠的随机样本：

- `optimization_pool`：只用于搜索候选波形；
- `selection_pool`：只用于从少量入围波形中选择最终波形；
- `validation_pool`：只做一次最终报告，禁止根据其结果继续调参。

随机种子必须写入配置和输出。若在看到 validation 结果后修改波形，原 validation 集即作废，必须换用新的、此前未生成过的种子并在报告中说明。

## 五、必须实现的三类波形

所有波形都必须返回 AOD 中心 (c(t))、深度 (D(t)) 及其一阶导数；为功-能检查和光滑度诊断，还应能稳定计算二阶导数，必要时计算三阶导数。

### A. Level 1 原始基线：`sequential_600us`

严格复现现有 Level 1：

- 总转移时间 600 μs；
- 前 52%：AOD 中心从 2.4 μm 按 Level 1 原有 smoothstep 移至 0，深度保持 280 μK；
- 后 48%：中心保持 0，深度按 Level 1 原有 ((1-s)^2) 形式从 280 μK 降至 0；
- 不得为了让 Level 2 更容易胜出而改变 Level 1 的实现、步长、初态或判据。

### B. 压缩顺序基线：`sequential_400us`

保持完全相同的 52%/48% 形状，只把总转移时间缩放为 400 μs。它用于区分“同步波形带来的改进”和“单纯缩短时间带来的退化”。

### C. 同步受约束波形：`overlapped_optimized_400us`

默认总转移时间 400 μs，允许移动和降深区间重叠。核心实现应采用一个低维、可解释、受约束的参数化族，以便在合理计算预算内优化。

定义五次 smootherstep：

\[
S(z)=10z^3-15z^4+6z^5,\qquad z\in[0,1].
\]

定义窗口函数：区间开始前为 0，结束后为 1，区间内为 (S[(s-a)/(b-a)])，其中 (s=t/T\in[0,1])。

至少优化以下参数：

```text
move_start_fraction
move_end_fraction
ramp_start_fraction
ramp_end_fraction
ramp_power
```

构造例如：

\[
c(t)=d\,[1-W_x(s)],
\qquad
D(t)=D_0\,[1-W_D(s)]^{p_D}.
\]

参数约束至少包括：

- (0\le a_x<b_x\le1)；
- (0\le a_D<b_D\le1)；
- 每个有效区间宽度不小于配置下限；
- (0.5\le p_D\le4)；
- (c(t)) 单调不增且始终位于 ([0,d])；
- (D(t)) 单调不增且始终位于 ([0,D_0])；
- 端点严格满足 (c(0)=d,c(T)=0,D(0)=D_0,D(T)=0)；
- 不允许样条过冲、负深度或超出起止位置；
- 端点速度和加速度应为零或数值上足够接近零；若深度幂导致端点导数奇异，候选必须判为无效。

### 可扩展的论文启发样条

代码结构应允许未来把位置和深度各自扩展为 14 个控制点。Level 2 默认可以使用 8 个含端点控制点的 PCHIP/单调样条做第二阶段精修，但不得要求在高维 14+14 参数空间中进行无计算预算限制的盲目全局搜索。

若实现样条精修：

- 由低维最优波形初始化；
- 默认 8 个控制点，配置可改为 14；
- 使用单调参数化或 PCHIP，严禁普通三次样条造成过冲；
- 把样条结果与低维波形分别报告，不得只保留较好结果而隐藏失败候选。

## 六、优化策略与防止蒙特卡洛过拟合

不得把 2,000 个验证 shots 放入每次优化迭代。采用分阶段预算：

### 阶段 1：候选探索

- 使用固定 `optimization_pool`，默认 128 shots；
- 使用 Latin hypercube、受约束随机搜索、SciPy differential evolution 的有限预算版本，或等价可复现的无导数方法；
- 默认最多评估 128-256 个有效候选；
- 优化积分步长可用 0.10 μs，但必须在后续阶段用正式步长复算；
- 所有候选及其参数、约束状态、得分和运行状态写入 `candidates.csv`，不得只保存最优候选。

### 阶段 2：精修与选择

- 取阶段 1 的前 5-10 个不同候选；
- 使用独立 `selection_pool`，默认 512 shots；
- 使用正式步长 0.05 μs；
- 可以进行小范围局部精修，但选择集不得变成无限次反复调参的隐藏训练集；
- 选择唯一最终波形，并在打开 validation 结果前冻结为 `best_waveform.yaml`。

### 阶段 3：独立验证

- 使用 `validation_pool`，默认至少 2,000 shots；
- 对三类主波形使用同一 validation 初态数组；
- 不再改变波形参数；
- 报告每种波形的 Wilson 95% 置信区间；
- 对波形间俘获结果使用配对 bootstrap 或等价配对方法，报告俘获率差的 95% 置信区间；
- 报告两波形对同一初态的四种配对计数：都成功、仅基线成功、仅优化成功、都失败；
- 不得仅凭点估计更高就声称有显著改进。

### 优化目标

优化以俘获为首要目标，以末态激发和波形光滑度为次要目标。可以采用加权软目标或字典序，但必须保存每个组成项，避免一个不可解释的总分。

建议优先级：

1. 最大化固定训练池中的俘获数或平滑俘获代理；
2. 最大化俘获裕量的低分位数；
3. 降低已俘获原子的末态归一化激发；
4. 惩罚过大的速度、加速度、jerk 和非法端点；
5. 不用验证集分数参与优化。

若使用平滑俘获代理，最终所有统计仍必须使用严格判据 (E_f<0)，不得用代理概率替代真实俘获率。

## 七、末态判据与物理量定义

在转移结束、AOD 深度严格为 0 时，定义 SLM 单阱末能量

\[
E_f=\frac12mv_f^2+U_{\rm SLM}(x_f).
\]

严格判据：

\[
E_f<0\Rightarrow\text{captured},
\qquad
E_f\ge0\Rightarrow\text{not captured}.
\]

所有 shots 都保留，不得删除接近阈值的样本。另行报告数值边界样本比例，例如

\[
|E_f|/D_{\rm SLM}<10^{-3},
\]

用于解释步长敏感性，但这些样本仍按严格符号计入俘获率。

对每个 shot 至少计算：

- 初始 AOD 束缚能量；
- 初始相对 AOD 阱底激发
  \(arepsilon_i=E_{\rm AOD,i}+D_{\rm AOD}(0)\)；
- 末态 SLM 单阱能量 (E_f)；
- 末态相对 SLM 阱底激发
  \(arepsilon_f=E_f+D_{\rm SLM}\)；
- 归一化俘获裕量
  \(M=-E_f/D_{\rm SLM}\)；
- 归一化末态激发
  \(e_f=\varepsilon_f/D_{\rm SLM}=1-M\)；
- 激发能变化 (Delta\varepsilon=\varepsilon_f-\varepsilon_i)；
- 深度变化功、移动功、总功和功-能残差；
- 100 μs 保持段中的 SLM 能量振荡幅度。

激发能统计必须明确是对“所有 shots”还是“仅已俘获 shots”。对于未俘获原子，不得把 (e_f>1) 与俘获态激发混在一起求一个难以解释的均值。

## 八、预检查：不过关不得启动优化

按顺序执行并保存 `preflight.json`：

1. 运行现有 Level 0 和 Level 1 测试；
2. 检查时变积分器在静态极限下与 Level 0 积分器一致；
3. 用固定初态和种子复现 Level 1 `sequential_600us`；如存在历史 Level 1 输出，比较主要指标；
4. 检查三类波形的端点、范围、单调性、连续性和导数；
5. 对理想初态 (x_0=c_{\rm AOD}(0),v_0=0) 运行三种波形；
6. 对理想轨迹检查功-能平衡；
7. 对理想轨迹将步长从 0.05 μs 减半到 0.025 μs，检查末态位置、速度、能量和功-能残差收敛；
8. 从固定初态池抽取至少 128 shots，对 0.05/0.025 μs 做配对步长检查，报告末能量差、俘获标签不一致率和阈值附近样本；
9. 在 100 μs 纯 SLM 保持段检查 SLM 能量误差幅度；
10. 任一关键检查失败时停止优化，输出失败原因，不得继续生成看似完整的“最优结果”。

建议默认数值验收尺度，可放入配置但不得无依据放宽：

- 理想轨迹最终位置和速度随步长表现出二阶收敛；
- 理想轨迹末态功-能残差除以 `max(D_AOD, D_SLM)` 小于 (5\times10^{-4})；
- MC 子集功-能残差的中位数小于 (10^{-3})，并报告最大值和 95% 分位数；
- 步长减半后的俘获标签不一致率不高于 2%，且不一致样本应主要集中在 (E_f\approx0) 附近；
- 保持段的 `peak_to_peak_energy_error_over_slm_depth` 小于配置阈值。

## 九、鲁棒性与敏感性检查

冻结最终优化波形后，至少进行以下检查，不得重新优化：

### 温度敏感性

在 3、5、7、10 μK 下分别用独立样本比较三类主波形。每个温度至少 500 shots，并报告 Wilson 区间。

### 末端对准误差

把 AOD 最终中心相对 SLM 设置为

```text
-0.10 μm, -0.05 μm, 0, +0.05 μm, +0.10 μm
```

其他参数保持不变，每个偏差至少 500 shots。这里的偏差只用于敏感性评估，不参与优化。

### 时间步长敏感性

对 nominal validation 初态中的固定子集比较 0.10、0.05 和 0.025 μs。必须区分连续量收敛和阈值分类跳变。

可选但推荐：对 SLM/AOD 深度各自做 ±10% 标定敏感性。若执行，必须清楚区分“一次只改变一个参数”和“联合最坏情形”。

## 十、建议配置

在不破坏现有配置系统的前提下新增 Level 2 配置，例如：

```yaml
level2:
  name: joint_transfer_optimization

atom:
  isotope: Cs-133
  mass_u: 132.90545196

traps:
  slm:
    depth_uK: 140.0
    waist_um: 1.17
    center_um: 0.0
    waist_source: paper_1061nm_slm_measurement
  aod:
    initial_depth_uK: 280.0
    waist_um: 1.17
    initial_center_um: 2.4
    waist_source: assumed_matched_to_slm_back_aperture

initial_ensemble:
  temperature_uK: 5.0
  sampler: level1_harmonic_rejection
  optimization_shots: 128
  selection_shots: 512
  validation_shots: 2000
  optimization_seed: 21001
  selection_seed: 21002
  validation_seed: 21003

integration:
  optimization_dt_us: 0.10
  validation_dt_us: 0.05
  convergence_dt_us: 0.025
  post_transfer_hold_us: 100.0

waveforms:
  sequential_baseline_duration_us: 600.0
  compressed_baseline_duration_us: 400.0
  optimized_duration_us: 400.0
  level1_move_fraction: 0.52
  level1_ramp_fraction: 0.48
  min_segment_fraction: 0.10
  ramp_power_bounds: [0.5, 4.0]

optimization:
  method: latin_hypercube_then_local_refine
  max_candidates: 192
  finalists: 8
  enable_monotone_spline_refinement: true
  spline_control_points: 8
  capture_weight: 1.0
  excitation_weight: 0.05
  smoothness_weight: 0.001

classification:
  capture_energy_threshold_J: 0.0
  near_threshold_fraction_of_slm_depth: 0.001

robustness:
  temperatures_uK: [3.0, 5.0, 7.0, 10.0]
  final_alignment_offsets_um: [-0.10, -0.05, 0.0, 0.05, 0.10]
  shots_per_condition: 500

output:
  directory: outputs/level2_joint_transfer
  save_candidate_table: true
  save_representative_trajectories: true
  representative_shots_per_outcome: 5
```

若现有 Level 1 已有不同但等价的配置结构，应合理整合而不是建立第二套互不兼容的配置系统。最终 `config_used.yaml` 必须包含 CLI 覆盖后的实际值和所有随机种子。

## 十一、代码结构与复用要求

先检查当前仓库。若 Level 0/Level 1 已位于同一个可安装包内，直接扩展现有包；若已分层为多个包，则保持现有风格。不得复制一套新的高斯势、单位常数、采样器或积分器来绕开复用。

至少需要清晰承担以下职责；文件名可与现有结构整合：

```text
waveforms.py       # 顺序、同步和样条波形及导数
objectives.py      # 候选评分、约束惩罚和组成项
optimizer.py       # 分阶段搜索、检查点和候选表
work_energy.py     # 外部功与功-能残差
statistics.py      # Wilson 区间、配对 bootstrap、配对计数
robustness.py      # 温度、对准和步长敏感性
level2_simulation.py
level2_visualization.py
level2_cli.py
```

所有公开函数使用简短中文 docstring。优化器必须支持中途保存检查点和从检查点继续；候选运行异常应记录为失败候选，不能让整轮搜索静默丢失。

必须提供清晰 CLI，例如：

```bash
python -m level0_static_trap.level2_cli \
  --config configs/level2_joint_transfer.yaml \
  --stage all \
  --output-dir outputs/level2_joint_transfer/demo
```

如果现有包名不同，使用现有包名，但 README 和最终回复必须给出实际成功执行的命令。`--stage` 至少支持：

```text
preflight
optimize
validate
robustness
all
```

不得依赖临时 `PYTHONPATH` 掩盖打包问题。运行现有项目要求的 editable install 或标准安装命令，并验证模块入口真实可用。

## 十二、输出文件

默认输出目录至少包含：

```text
outputs/level2_joint_transfer/demo/
├── config_used.yaml
├── preflight.json
├── candidates.csv
├── optimization_history.csv
├── optimization_checkpoint.json
├── best_waveform.yaml
├── waveforms.csv
├── validation_shots.csv
├── robustness.csv
├── metrics.json
├── image_validation.json
├── summary.md
├── waveform_comparison.png
├── potential_evolution.png
├── representative_trajectories.png
├── final_phase_space.png
├── capture_rate_comparison.png
├── excitation_and_margin.png
├── optimization_history.png
├── work_energy_balance.png
├── timestep_convergence.png
└── robustness_heatmaps.png
```

### `candidates.csv`

每个候选一行，至少包含：

- 唯一候选 ID；
- 全部波形参数；
- 是否满足约束；
- 失败原因；
- 使用的样本池、步长和随机种子；
- 严格俘获率；
- 平滑代理值（若使用）；
- 裕量分位数；
- 已俘获样本的末态激发统计；
- 最大速度、加速度和 jerk；
- 目标函数各组成项与总分；
- 是否入围和是否最终选中。

### `waveforms.csv`

对三个主波形使用足够密的公共时间轴，至少保存：

```text
waveform
time_s
time_us
aod_center_m
aod_center_um
aod_depth_J
aod_depth_uK
center_velocity_m_per_s
center_acceleration_m_per_s2
depth_rate_J_per_s
```

### `validation_shots.csv`

每个 validation shot 和每个主波形一行，至少保存：

- shot ID 和初态；
- 波形名称；
- 末态位置、速度、SLM 能量；
- captured 布尔值；
- near-threshold 布尔值；
- capture margin；
- initial/final excitation 和变化；
- move work、depth work、total work；
- 功-能残差；
- 保持段能量误差幅度。

### `metrics.json`

至少包含：

- 版本、配置、种子和代码运行信息；
- 全部预检查结果；
- 三个主波形的端点和光滑度指标；
- 优化预算、有效/失败候选数、训练与选择结果；
- 冻结的最终参数；
- validation 中每个波形的俘获率和 Wilson 区间；
- 波形间配对俘获率差及置信区间；
- 配对成功/失败计数；
- 俘获裕量和激发分布统计；
- 功-能残差统计；
- 步长收敛和边界标签不一致率；
- 温度和对准鲁棒性；
- 是否观察到统计上明确的改进；
- 图片程序化检查及是否真正做过视觉检查。

JSON 不得包含 NaN 或 Infinity。

## 十三、必须生成的诊断图

1. `waveform_comparison.png`：三种主波形的 AOD 中心和深度随时间变化；明确标出移动/降深重叠区。
2. `potential_evolution.png`：优化波形下若干关键时刻的总势剖面，或时间-位置势能热图；标出 SLM 与 AOD 中心。
3. `representative_trajectories.png`：成功、失败和阈值附近的代表轨迹，叠加 AOD 中心；不得只挑成功轨迹。
4. `final_phase_space.png`：末态 (x,v) 散点，按 captured 着色，并画出 (E_{\rm SLM}=0) 分界线。
5. `capture_rate_comparison.png`：validation 俘获率及 Wilson 95% 区间；训练结果与独立验证结果视觉上区分。
6. `excitation_and_margin.png`：已俘获原子的末态激发、俘获裕量和 (Delta\varepsilon) 分布。
7. `optimization_history.png`：候选评估顺序、得分组成、最佳值演化和无效候选。
8. `work_energy_balance.png`：理想轨迹及至少一个普通 MC 轨迹的 (Delta E)、移动功、深度功、总功和残差。
9. `timestep_convergence.png`：连续量误差与俘获标签不一致率随步长变化。
10. `robustness_heatmaps.png`：温度与末端对准偏差下的俘获率或分别绘制两组敏感性图。

绘图使用非交互后端，至少 200 dpi，保存后关闭 figure，不调用阻塞式 `plt.show()`。

## 十四、图片与文件验收

对每张 PNG 必做程序化检查并写入 `image_validation.json`：文件存在、可重新打开、尺寸非零、像素不是全白/全黑/近乎单色、数据和坐标范围有限、图例系列数量合理、figure 已关闭。

这些检查不能证明标签不重叠或曲线美观：

- 若环境提供图像查看工具，逐张打开最终图并记录 `visual_review_performed: true`；
- 若没有，记录 `visual_review_performed: false`，并明确写“完成程序化检查，未完成视觉审阅”；
- 未真正查看图片时不得声称“图片清晰可读”或“坐标轴合理”。

CSV 必须检查列名、行数、有限值、shot ID 对齐和三波形使用相同初态。JSON/YAML/Markdown 必须重新读取验证。不得只检查文件是否存在。

## 十五、测试要求

在保留所有 Level 0/Level 1 测试的基础上，至少新增：

1. 三类波形的端点严格正确；
2. 同步波形位置和深度单调且不越界；
3. 非法参数、负深度、过冲和奇异端点导数被拒绝；
4. 波形解析导数与有限差分一致；
5. `sequential_600us` 与 Level 1 固定轨迹一致；
6. 时变积分器静态极限与 Level 0 一致；
7. 显式时变势满足功-能平衡并随步长收敛；
8. AOD 关闭后的纯 SLM 保持段能量误差有界；
9. common-random-number 比较中三波形初态逐行完全相同；
10. optimization、selection、validation 种子互不相同；
11. 最终 validation 数据未被优化器访问；可通过接口分层或测试 spy 验证；
12. Wilson 区间在已知小样本案例上正确；
13. 配对 bootstrap 或配对统计在构造数据上正确；
14. (E_f<0) 判俘获，(E_f=0) 不俘获；
15. near-threshold 样本仍被保留并按严格符号计数；
16. 步长减半时末态连续量收敛，标签变化集中在阈值附近；
17. 单波形运行失败被写入候选表，不会被静默忽略；
18. CLI 各 stage 的 smoke test；
19. CSV/JSON schema、严格 JSON 和图片程序化检查；
20. 输出目录位于 `sources/` 时被拒绝。

测试不得要求优化波形一定优于基线，因为物理结果不是单元测试。应测试优化流程可复现、约束正确、冻结与独立验证没有数据泄漏。

## 十六、执行顺序

必须按以下顺序实际执行：

1. 检查当前仓库、现有 Level 0/Level 1 代码、配置、测试和输出；
2. 运行现有测试并记录基线；
3. 扩展配置和代码，复用已有常数、势、采样器和积分器；
4. 实现三类波形、导数和约束；
5. 实现功-能核算和统计函数；
6. 编写 Level 2 测试；
7. 运行 `preflight`，失败则修复并重跑；
8. 运行候选探索并保存全部候选和检查点；
9. 用 selection pool 精修、选择并冻结最终波形；
10. 用独立 validation pool 一次性验证三个主波形；
11. 冻结参数后运行温度、对准和步长敏感性；
12. 生成全部数据、报告和图；
13. 程序化检查所有输出；
14. 若有图像查看工具，逐图视觉检查并修复展示问题；
15. 重新运行完整测试和 nominal validation，确保最终输出对应最新代码；
16. 检查未修改 `sources/`；
17. 给出简洁、可核查的最终回复。

## 十七、最终回复必须包含

- 实际扩展或创建的主要模块；
- 实际成功执行的安装、测试和 CLI 命令；
- 测试通过/失败/跳过数量；
- preflight 是否全部通过；
- 评估候选数、无效候选数和冻结的最佳参数；
- 600 μs 顺序基线、400 μs 顺序基线、400 μs 优化同步波形的 validation 俘获率和 Wilson 95% 区间；
- 优化波形相对两个基线的配对俘获率差及置信区间；
- 末态激发与俘获裕量结果；
- 功-能残差和步长收敛结果；
- 温度与末端对准敏感性；
- 是否存在统计上明确的提升；若没有，直接说明；
- 程序化图片检查结果及是否真正执行视觉检查；
- 输出目录和关键文件的绝对路径；
- Level 2 仍未包含、应留给后续 Level 的物理过程。

不得只汇报“最优俘获率”，也不得用 optimization/selection 数据替代独立 validation 结果。
