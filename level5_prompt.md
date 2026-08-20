# Codex 执行 Prompt：Level 5 有限脉冲、自旋噪声与多站点交错随机基准测试

请在当前工作区中，在已完成并验证的 Level 0-4 基础上，实际完成“Level 5：沿端到端原子运动轨迹的有限脉冲单比特动力学、噪声谱和多站点参考/交错随机基准测试”。

不要只给方案或代码片段。必须检查并复用 Level 0-4 的代码和冻结产物，建立可验证的量子-经典混合模拟，实际运行预检查、解析基准、单通道表征、随机基准测试、多站点统计、消融与敏感性分析，并生成结构化数据、图表和中文报告。

如果工作区存在 `sources/`，其中内容全部只读，不得编辑、移动、重命名或删除；如果不存在，直接继续。任何输出路径解析后都不得位于 `sources/`。

论文依据：arXiv:2403.12021v4 的 Fig. 4-6、Extended Data Fig. 7-10、Methods “Characterizing the atomic qubits”“Single-qubit gate randomized benchmarking”“Randomized benchmarking of coherent transport”“Combined atom transfer and move”。官方版本：https://arxiv.org/abs/2403.12021v4

## 一、Level 5 的定位

Level 4 已经计算了：

- 完整 pick-up—long transport—return—drop-off 三维经典轨迹；
- 逐阶段 survival、heating 和功-能账本；
- 由差分光移积分得到的半经典相位；
- 理想瞬时 echo/XY4 toggling 下的 ensemble contrast。

Level 5 将内部态模型升级为有限时长、受噪声驱动的单量子比特动力学，并执行可复现的随机基准测试。

Level 5 必须完成：

1. 沿 Level 4 的三维经典轨迹传播 \(^{133}\mathrm{Cs}\) 超精细钟态量子比特；
2. 实现旋转框架下的有限时长微波脉冲；
3. 实现 bare rotations 和 SCROFULOUS 复合脉冲；
4. 实现轨迹相关差分光移、磁场噪声和强度噪声；
5. 支持从 PSD 或时域模型生成相关噪声；
6. 实现理想、有限脉冲 spin echo、XY4 和 XY16；
7. 构造并验证完整 24 元素单比特 Clifford 群；
8. 执行 reference RB 和 transport-interleaved RB；
9. 区分原子损失、条件自旋返回概率和 survival-weighted 返回概率；
10. 通过 Pauli transfer matrix/process characterization 独立表征单次运输通道；
11. 在 47/195 个独立原子站点上加入共同噪声、局域噪声和参数不均匀性；
12. 比较站点平均与站点分辨结果；
13. 检查拟合可辨识性、非指数衰减和数据泄漏；
14. 明确区分模拟量与论文实验 IRB fidelity。

本 Level 5 仍不模拟 Rydberg 相互作用、两比特门、纠缠、多体波函数、测量反馈、量子纠错、原子间碰撞或光学串扰。多站点模拟默认是条件独立的单原子动力学加相关经典噪声，不能称为多体量子模拟。

## 二、冻结上游数据和两种耦合模式

运行前发现并冻结：

1. Level 4 的 `upstream_manifest.json`；
2. `frozen_validation_protocols.yaml`；
3. Protocol A/B 的三维经典 checkpoint 和代表/validation 轨迹；
4. Level 4 survival、first-failure、phase 和功-能结果；
5. Level 3 lensing calibration；
6. Level 2 transfer 波形。

写入新的 `upstream_manifest.json`，包含绝对路径、SHA-256、schema/version、生成时间和关键参数。

支持两种模式：

### A. Frozen-trajectory spin replay

复用 Level 4 已保存的经典轨迹，在其上重复不同量子脉冲和噪声 realization。默认使用此模式进行大量 RB，因为它避免重复昂贵的三维动力学。

必须：

- 对位置、势能和深度使用经过验证的插值；
- 在脉冲边界、DD 翻转点和 protocol segment 边界显式插入时间点；
- 不改变 Level 4 的 survival 标签；
- 记录轨迹插值误差；
- 允许同一经典轨迹配对比较多种脉冲方案。

### B. Joint classical-spin propagation

对预检查和小规模验证，经典位置/速度与量子态在同一时间网格推进。默认不含量子态对机械运动的反作用。

必须用该模式验证 frozen replay 的 phase 和量子结果；若两者不一致超过阈值，停止主 RB。

默认采用 weak-state-dependence approximation：两个钟态共享同一经典势，差分极化率只进入相对相位。若以后启用 state-dependent force，必须作为单独敏感性模型，不能混入 nominal。

## 三、单比特 Hamiltonian

使用

\[
|0\rangle\equiv|6S_{1/2},F=3,m_F=0\rangle,
\qquad
|1\rangle\equiv|6S_{1/2},F=4,m_F=0\rangle.
\]

在微波载频旋转框架和 rotating-wave approximation 下，对第 \(i\) 个原子：

\[
H_i(t)=\frac{\hbar}{2}
\left[
\Delta_i(t)\sigma_z
+\Omega_i(t)
\{\cos[\phi_i(t)]\sigma_x+\sin[\phi_i(t)]\sigma_y\}
\right].
\]

总失谐：

\[
\Delta_i(t)=
\Delta_{\rm drive}
+\Delta_{{\rm LS},i}(t)
+\Delta_{B,i}(t)
+\Delta_{{\rm site},i}.
\]

轨迹相关差分光移：

\[
\Delta_{{\rm LS},i}(t)=
\frac{\eta_{\rm SLM}U_{\rm SLM}[\boldsymbol r_i(t),t]
+\eta_{\rm AOD}U_{\rm AOD}[\boldsymbol r_i(t),t]}{\hbar}.
\]

nominal 使用 \(\eta=1.3\times10^{-4}\)，并扫描论文理论量级 \(1.5\times10^{-4}\)。两个波长的 \(\eta\) 必须独立可配置。

默认微波 Rabi 频率使用论文锚点

\[
\Omega_0=2\pi\times24.611\ {\rm kHz}.
\]

该值来自全阵列平均实验结果。站点局部 Rabi 频率、幅度误差和空间梯度必须另行建模，不得把平均值当作每个站点精确值。

### 模型边界

- 驱动相位和幅度为经典控制；
- 默认采用 RWA；
- 不含 Bloch-Siegert shift，除非在敏感性中显式加入；
- 不含真实多能级泄漏，除非提供经验证的扩展模型；
- 不含自发 Raman scattering，除非以独立 Lindblad rate 明确配置；
- clock frequency 的大载频不直接进入步进，只在旋转框架失谐中体现。

## 四、量子态推进

至少支持：

### A. State-vector stochastic propagation

对每个经典噪声 realization 传播归一化两分量态矢。时间步内 Hamiltonian 视为常数时，使用解析 SU(2) 指数或稳定矩阵指数：

\[
|\psi(t+\Delta t)\rangle
=\exp[-iH(t+\Delta t/2)\Delta t/\hbar]|\psi(t)\rangle.
\]

### B. Density-matrix/Lindblad propagation

用于可选的 \(T_1\)、Markovian pure dephasing 或 leakage surrogate：

\[
\dot\rho=-\frac{i}{\hbar}[H,\rho]
+\sum_k\left(L_k\rho L_k^\dagger
-\frac12\{L_k^\dagger L_k,\rho\}\right).
\]

默认 nominal 可以关闭 Lindblad 项。启用时每个 rate 必须有配置和 provenance。

每一步检查：

- state-vector norm；
- density-matrix trace、Hermiticity、最小 eigenvalue；
- 传播器 unitarity；
- NaN/Infinity；
- 时间步和脉冲边界。

量子步长可以小于经典轨迹步长，但必须通过插值和收敛测试。不得用脉冲面积正确作为唯一数值验证。

## 五、脉冲库

### Bare rotation

定义

\[
R_\phi(\theta)
=\exp\left[
-\frac{i\theta}{2}
(\cos\phi\,\sigma_x+\sin\phi\,\sigma_y)
\right],
\]

并用有限时长矩形、平滑边沿或配置的 envelope 实现。至少支持：

- \(X_{\pm\pi/2},X_\pi\)；
- \(Y_{\pm\pi/2},Y_\pi\)；
- virtual \(Z\) rotation；
- 任意 equatorial rotation。

### SCROFULOUS

按论文 Methods 和其引用的标准定义实现对脉冲长度误差鲁棒的三段对称复合脉冲。必须：

- 明确 `sinc` 和 `arcsinc` 的定义、归一化方式及所选根分支；
- 数值求解不存在或分支不唯一时给出清晰错误；
- 验证目标 unitary；
- 扫描静态幅度误差并与 bare pulse 比较；
- 正确累加总脉冲面积和实际持续时间；
- 特殊处理 identity 和接近零角度；
- 不得仅硬编码论文平均面积 2.02π 而不构造实际 pulse list。

### Pulse envelopes

至少实现：

- square；
- cosine-edge square；
- 可选 Gaussian/DRAG-like envelope，但没有多能级模型时不得声称 DRAG 可抑制 leakage。

脉冲标记必须包含 start/end、phase、nominal angle、amplitude、site scaling、composite sequence ID 和 Clifford ID。

## 六、DD 序列与有限脉冲

至少实现：

- none；
- spin echo；
- XY4；
- XY8；
- XY16；
- custom pulse schedule。

Level 4 的 toggling-function 结果应在以下极限恢复：

- 脉冲时长趋于零；
- pulse error 为零；
- 只存在 \(\sigma_z\) 纯退相位；
- RWA 有效。

有限脉冲下必须实际传播 Hamiltonian，不能在脉冲期间仅翻转 \(y(t)\)。报告：

- free-evolution phase；
- pulse-during-noise contribution；
- pulse area/phase error；
- DD block 的总时间；
- 在每个 transport move 内 pulse placement。

实现 transport-aware XY4：outbound 和 inbound 各一个对称 block。可选实现 transformed-Clifford-frame DD，但只有在 Clifford 共轭、pulse axes 和最终 inverse 均经过群论测试后才能启用。

不得因实现了 XY4 就声称复制论文的 transformed Clifford frame 技术；两者必须作为不同模式。

## 七、噪声模型与 PSD

总 detuning/amplitude noise 至少支持：

### A. Quasistatic

- shot-to-shot detuning；
- shot-to-shot Rabi amplitude scaling；
- site depth/Rabi offset。

### B. Ornstein-Uhlenbeck

给定 RMS 和 correlation time 的 detuning 或 fractional intensity noise。

### C. PSD-synthesized Gaussian noise

支持：

- white noise；
- \(1/f^\alpha\)；
- Lorentzian；
- 60 Hz 窄带线及 harmonics；
- 用户给定 tabulated one-sided PSD。

使用 FFT 或经验证的滤波方法生成时域噪声。必须验证：

- Hermitian symmetry 和实数时序；
- one-sided/two-sided PSD 约定；
- Hz 与 rad/s；
- sample rate、Nyquist 和 record length；
- 生成样本的积分方差与目标 PSD 一致；
- Welch PSD 的统计区间；
- 不同随机 seed 复现/独立。

### D. Common/local spatial correlations

第 \(i\) 个站点噪声至少支持

\[
\delta_i(t)=
\sqrt{\rho}\,\delta_{\rm common}(t)
+\sqrt{1-\rho}\,\delta_{{\rm local},i}(t),
\qquad 0\le\rho\le1.
\]

也可支持基于站点距离的相关核。必须验证经验相关矩阵。

### 参数来源

nominal 噪声参数若没有实验数据，默认设为零。敏感性中的 RMS、谱指数、相关时间和 60 Hz 线强度必须标记为 `assumed`。

可以建立独立的 paper-coherence anchor 配置，对照论文：

- ensemble \(T_2^*=14.0(1)\) ms；
- site-averaged \(T_2^*=25.5\) ms；
- 约 4.3 μK 的温度估计；
- reduced trap depth 55 μK；
- XY16 下 \(T_2=12.6(1)\) s。

该 anchor 与 Level 4 的 transport 协议参数不同。不得用同一组自由噪声参数反复拟合所有指标后再称为独立预测。若做校准，必须分离 calibration 和 validation observables，并报告不可辨识性。

## 八、单次通道的过程表征

在进入 RB 前，分别表征：

- static idle same duration；
- transfer-only；
- transport-only；
- full roundtrip；
- Level 4 Protocol A/B；
- none/echo/XY4/XY16；
- lensing on/off；
- noise on/off。

对输入态

\[
\{|0\rangle,|1\rangle,|+\!x\rangle,|-\!x\rangle,
|+\!y\rangle,|-\!y\rangle\}
\]

传播并重建 survival-conditioned Pauli transfer matrix（PTM）。至少输出：

- affine PTM；
- trace preservation/unitality diagnostics；
- average state fidelity；
- average gate/channel fidelity；
- unitarity；
- coherent rotation component；
- dephasing/contraction axes；
- survival probability；
- leakage/loss 单独通道。

不能把非 trace-preserving loss map 强行归一化后只报一个 fidelity。必须并列报告：

1. conditional spin channel；
2. survival；
3. survival-weighted performance；
4. loss-detected convention 下的结果。

如果使用

\[
F_{\rm avg}=\frac12+\frac{\operatorname{Tr}T}{6}
\]

必须确认这是 trace-preserving、unital qubit channel 的 \(3\times3\) Bloch block \(T\)；不满足条件时使用一般公式或报告不适用。

用非参数 bootstrap 跨 classical/noise realizations 给出区间。

## 九、Clifford 群和序列

构造完整 24 元素单比特 Clifford 群。每个元素必须包含：

- 规范 unitary，忽略 global phase 后唯一；
- 到 primitive rotations 的分解；
- bare 和 SCROFULOUS pulse realization；
- inverse Clifford；
- multiplication table 或可验证 composition。

测试：

- 恰有 24 个唯一元素；
- 群封闭；
- 每个 inverse 正确；
- 对 Pauli 算符的共轭仍为带符号 Pauli；
- ideal pulse realization 与规范 unitary 一致；
- 随机序列与计算出的 inverse 在无噪声时返回初态；
- Clifford sampler 均匀。

若使用论文 SCROFULOUS 分解，报告平均 pulse area，并与论文的约 2.02π 锚点比较；不一致时检查 decomposition，而不是强行归一到 2.02π。

保存所有随机 Clifford 序列和 seed，确保拟合可重放。

## 十、Reference RB

先在 static idle 条件下实现标准 reference RB：

1. 对每个 sequence length \(n\) 均匀抽取 \(n\) 个 Clifford；
2. 计算精确 inverse；
3. 使用有限脉冲和指定噪声传播；
4. 测量回到初态的概率；
5. 对 \(|0\rangle\) 和可选多个输入态重复；
6. 保存每条序列而非只保存平均。

默认 lengths：

```text
[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000]
```

开发阶段每点 8 条序列；正式阶段每点至少 60 条，默认 72 条，以呼应论文实验设计。多噪声 realization 数由计算预算配置，必须报告。

拟合模型：

\[
P(n)=A p^n+B.
\]

报告 SPAM 参数、\(p\)、置信区间、残差和拟合窗口。若噪声具有强时域相关或 gate dependence，单指数可能不适用；必须比较 stretched exponential、分段/非参数趋势或直接 PTM 结果，不得无条件给出 RB fidelity。

## 十一、Transport-interleaved RB

实现论文启发的固定总门数方案：

- 固定总 Clifford 数 \(N=80\)；
- 在前 \(M\) 个 Clifford 之间各插入一次 transport/roundtrip channel；
- 之后施加剩余 \(N-M\) 个 Clifford，使不同 \(M\) 的总 gate 数相同；
- 最后施加总 inverse；
- reference 与 transported 使用相同 Clifford strings、初态和尽可能配对的噪声 realization；
- 默认
  \[
  M\in\{0,1,2,4,8,12,20,30,40,60,80\}.
  \]

明确配置 interleaved operation 的含义：

- `one_way_transport`；
- `aod_roundtrip`；
- `full_pickup_transport_dropoff_roundtrip`；
- `static_idle_matched_duration`。

不得把这些不同 channel 混成一个“move”。

对每个 \(M\) 默认至少 72 个随机 Clifford strings。保存：

- atomic survival；
- conditional spin return probability；
- survival-weighted return；
- loss-detected return convention；
- site/noise/sequence resolved values；
- pulse count、total duration 和 DD count。

### 非指数和损失处理

必须分别拟合 survival 与 conditional spin signal。允许比较论文启发的 clipped-Boltzmann × depolarization 模型，但不得强制。

至少报告：

- 早期 \(M\) 的模型无关配对斜率；
- exponential fit；
- 非指数诊断；
- instantaneous ratios 及其 bootstrap 区间；
- 尾部 risk-set size；
- fit parameter covariance/identifiability。

当 survival 接近零或相邻比值不稳定时，不输出误导性的 instantaneous fidelity。

### IRB fidelity

只有在 reference/interleaved RB 的 gate-independence、Markovianity、拟合质量和 loss convention 满足前提时，才输出标准 IRB estimate，并同时给出 process/PTM estimate。

如果前提不成立，明确写 `IRB estimate not interpretable under standard assumptions`，仍保留原始曲线和通道表征。

## 十二、多站点并行模拟

默认模拟独立原子站点集合：

- 开发：8 sites；
- nominal validation：47 sites；
- 扩展：195 sites；
- 6100 sites 只允许使用经过验证的降阶/分布抽样，不得直接盲目全量传播。

每个 site 可具有：

- SLM/AOD depth offset；
- temperature offset；
- Rabi amplitude scale；
- static detuning；
- \(\eta\) variation；
- alignment/lensing variation；
- local noise；
- common noise。

论文报告全阵列 trap-depth 标准差约 11.4%，可作为单独 paper-inspired scenario；不得自动当作所有 Level 4 协议的精确分布。

至少输出：

- site-resolved survival；
- site-resolved reference/IRB decay；
- site-resolved channel metrics；
- array average；
- between-site variance；
- within-site noise variance；
- common-mode vs local contribution；
- worst decile 和 best decile；
- 与站点参数的相关性。

所有站点使用相同 Clifford strings，以隔离空间不均匀性；noise realizations 按 common/local 结构生成。

不得把 47 或 195 个独立站点的平均称为包含原子间相关或并行 gate crosstalk 的模拟。

## 十三、统计设计和数据隔离

至少分离：

- `unit_test_seed_space`；
- `noise_calibration_pool`；
- `channel_validation_pool`；
- `reference_rb_pool`；
- `interleaved_rb_pool`；
- `site_validation_pool`；
- `sensitivity_pool`。

如果执行噪声校准，校准 observables 和最终 validation observables 不得相同。

使用 hierarchical bootstrap：

1. 重采样 Clifford sequence；
2. 在 sequence 内重采样 site；
3. 在 site 内重采样 classical/noise realization。

报告 bootstrap 层级。不得把所有样本当独立 Bernoulli 导致置信区间虚假变窄。

对确实由逐 shot 的 0/1 survival 或 loss 计数形成的基础比例，同时保存分子、分母并报告 Wilson 95% 区间，作为有限样本描述性统计。Wilson 区间不能替代上述 hierarchical bootstrap，也不能用于把有站点/序列/噪声相关的数据伪装成独立 Bernoulli；conditional return 的分母必须是实际 surviving risk set。若模拟直接输出连续概率而非二项观测，不得套用 Wilson 区间，应使用匹配采样层级的 bootstrap。

对多重条件比较报告 effect size 和区间；如进行显著性检验，说明多重比较处理。

## 十四、计算预算和缓存

Level 5 计算量很大，必须分阶段：

### Stage A：analytic/preflight

无噪声、小规模，验证 Hamiltonian、脉冲、Clifford、DD、PTM 和 Level 4 极限。

### Stage B：single-channel characterization

至少 2,000 条 frozen classical trajectories 或经批准的配对子样本；对多个输入态和 DD 模式批量传播。

### Stage C：RB development

8 sites、8 sequences/length、少量 noise realizations，验证全管线和拟合。

### Stage D：formal reference/IRB

- 47 sites；
- 72 sequences/condition；
- 预先冻结的 sequence lengths 和 \(M\)；
- noise realizations 数根据预算预注册；
- 使用向量化 SU(2)/Bloch 推进、channel caching 或批处理。

缓存只有在物理上允许时使用：

- Markovian、独立重复 channel 可以缓存 superoperator；
- 长相关时间噪声跨 moves 时不能把每次 move 当独立缓存；
- heating 改变后续经典轨迹时不能复用固定单次 channel；
- 使用近似缓存必须与连续时间小样本比较。

### Stage E：195-site extension

使用冻结模型，只改变 site 数。若采用降阶 channel surrogate，必须在 47-site/full propagation 上验证。

所有阶段保存 checkpoint，可恢复运行。不得因资源不足自动减少样本而不在输出中记录实际数量。

## 十五、建议配置

在不破坏现有配置系统的前提下新增 Level 5 配置，例如：

```yaml
level5:
  name: finite_pulse_irb_array

upstream:
  level4_output: auto_discover
  require_sha256_manifest: true
  trajectory_mode: frozen_replay
  joint_validation_shots: 128

qubit:
  basis: cs133_clock_states
  drive_rabi_frequency_kHz: 24.611
  drive_detuning_Hz: 0.0
  eta_slm: 1.3e-4
  eta_aod: 1.3e-4
  rotating_wave_approximation: true
  include_state_dependent_force: false

quantum_integration:
  method: analytic_su2_midpoint
  max_dt_us: 0.05
  convergence_dt_us: 0.025
  norm_tolerance: 1.0e-10
  density_matrix_min_eigenvalue_tolerance: -1.0e-10

pulses:
  default_family: scrofulous
  envelope: square
  cosine_edge_us: 0.0
  include_amplitude_error: true
  include_phase_error: false

dd:
  modes: [none, spin_echo, xy4, xy16]
  transport_mode: xy4_per_long_move
  transformed_clifford_frame: false
  pulses_are_finite: true

noise:
  nominal_enabled: false
  detuning_model: none
  amplitude_model: none
  psd_model: none
  common_fraction: 0.0
  noise_seed: 51999
  provenance: assumed_unless_calibrated

paper_coherence_anchor:
  enabled: true
  trap_depth_uK: 55.0
  temperature_uK: 4.3
  ensemble_t2star_ms: 14.0
  site_t2star_ms: 25.5
  xy16_t2_s: 12.6
  use_for_calibration: false

channel_characterization:
  classical_trajectories: 2000
  bootstrap_samples: 2000
  input_states: [z_plus, z_minus, x_plus, x_minus, y_plus, y_minus]

reference_rb:
  lengths: [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000]
  sequences_per_length: 72
  sequence_seed: 51001

interleaved_rb:
  total_cliffords: 80
  interleaved_counts: [0, 1, 2, 4, 8, 12, 20, 30, 40, 60, 80]
  sequences_per_count: 72
  operation: full_pickup_transport_dropoff_roundtrip
  sequence_seed: 51002

array:
  development_sites: 8
  validation_sites: 47
  extension_sites: 195
  site_seed: 51003
  trap_depth_relative_std: 0.0
  rabi_relative_std: 0.0
  static_detuning_std_Hz: 0.0
  spatial_noise_correlation: 0.0

statistics:
  bootstrap_samples: 2000
  hierarchical_bootstrap: true
  confidence_level: 0.95

output:
  directory: outputs/level5_finite_pulse_irb/demo
  save_sequence_resolved_data: true
  save_site_resolved_data: true
  save_noise_realizations_for_replay: true
```

nominal noise 默认关闭是为了避免将假设参数伪装成实验输入。正式 noise scenario 必须通过单独配置启用并注明来源。

## 十六、代码结构与复用

必须复用 Level 4 的轨迹、survival、协议、差分光移和输出系统。建议职责：

```text
qubit_hamiltonian.py        # rotating-frame Hamiltonian
quantum_propagators.py      # SU(2)、Bloch、Lindblad
pulse_library.py            # bare、SCROFULOUS、envelopes
clifford_group.py           # 24 Clifford、分解、inverse
dd_finite_pulses.py         # echo、XY4/8/16、frame transform
noise_psd.py                # PSD 和时域噪声生成
spatial_noise.py            # common/local/site correlation
trajectory_replay.py        # Level 4 轨迹插值与边界
channel_tomography.py       # PTM 和 channel metrics
reference_rb.py
interleaved_rb.py
rb_fitting.py
array_ensemble.py
level5_statistics.py
level5_visualization.py
level5_cli.py
```

所有公开函数包含简短中文 docstring。不得为不同 pulse family 或 protocol 复制传播器。

沿用 Level 0 已建立的可安装包结构。若代码采用 `src/` layout，仓库根目录必须已有或补齐有效的 `pyproject.toml`（或等价标准构建配置），且包发现规则必须覆盖新增 Level 5 模块。不得创建第二套互相冲突的打包配置。运行 Level 5 前必须从项目根目录实际执行 editable/标准安装，随后在不设置临时 `PYTHONPATH` 的新进程中验证包导入和模块入口；若出现 `ModuleNotFoundError`，preflight 必须失败而不是跳过。

CLI 示例：

```bash
python -m level0_static_trap.level5_cli \
  --config configs/level5_finite_pulse_irb.yaml \
  --stage all \
  --output-dir outputs/level5_finite_pulse_irb/demo
```

实际包名不同则沿用现有名称。`--stage` 至少支持：

```text
preflight
channel
reference_rb
interleaved_rb
array_extension
sensitivity
all
```

不得使用临时 `PYTHONPATH` 掩盖安装问题。

## 十七、preflight：失败则停止正式 RB

按顺序执行并写入 `preflight.json`：

1. 运行 Level 0-4 全部测试；
2. 从项目根目录验证构建配置、实际安装、干净进程导入和 `python -m <package>.level5_cli --help`，不得依赖临时 `PYTHONPATH`；
3. 校验 Level 4 轨迹、hash、time grids、survival 和 schema；
4. frozen replay 与 joint propagation 的差分光移和 phase 一致；
5. 常 Hamiltonian 的解析 Rabi oscillation；
6. 常 detuning 下的解析 Ramsey；
7. quasistatic detuning 下 echo 极限；
8. state-vector norm 和 density-matrix trace/positivity；
9. SU(2) propagator unitarity；
10. bare pulse 的目标 unitary；
11. SCROFULOUS 的目标 unitary 和幅度误差鲁棒性；
12. pulse envelope 面积和持续时间；
13. none/echo/XY4/XY16 的 pulse order；
14. Level 4 instantaneous toggling 极限；
15. PSD normalization、Welch spectrum 和 correlation；
16. common/local 经验相关矩阵；
17. 完整 24 Clifford 群、inverse 和均匀 sampler；
18. ideal RB 所有随机序列返回概率为 1；
19. 已知 depolarizing/dephasing channel 的 PTM 和 \(F_{\rm avg}\)；
20. 构造 reference RB 数据的参数回收；
21. 构造 interleaved channel 的参数回收；
22. loss 和 conditional channel 口径分离；
23. 量子步长 0.10/0.05/0.025 μs 收敛；
24. pulse boundary/segment boundary 无重复或漏积分；
25. seed 空间隔离和 sequence replay；
26. validation 条件冻结且未被 calibration 访问。

任何关键项失败时停止正式 RB，保存失败原因。不得跳过失败后输出 fidelity。

## 十八、输出文件

默认至少生成：

```text
outputs/level5_finite_pulse_irb/demo/
├── config_used.yaml
├── upstream_manifest.json
├── frozen_rb_design.yaml
├── preflight.json
├── pulse_library.csv
├── clifford_table.csv
├── noise_psd.csv
├── noise_validation.json
├── channel_ptm.csv
├── channel_metrics.json
├── rb_sequences.jsonl
├── reference_rb_sequence_results.csv
├── interleaved_rb_sequence_results.csv
├── rb_fit_metrics.json
├── site_resolved_metrics.csv
├── array_summary.csv
├── sensitivity.csv
├── metrics.json
├── image_validation.json
├── summary.md
├── pulse_response.png
├── scrofulous_robustness.png
├── noise_psd_validation.png
├── ramsey_echo_dd.png
├── channel_ptm.png
├── channel_bloch_sphere.png
├── reference_rb_decay.png
├── interleaved_rb_decay.png
├── survival_vs_spin_return.png
├── rb_residuals.png
├── site_resolved_array.png
├── coherence_vs_round.png
├── timestep_convergence.png
└── sensitivity_panels.png
```

### Sequence-resolved 数据

每行至少包含 sequence ID、seed、Clifford list、inverse、site、classical trajectory ID、noise realization、pulse family、DD、survival、conditional return、joint return、duration、pulse count 和 fit group。

不得只保存 averaged RB curve。

### `channel_metrics.json`

至少包含每个 channel/DD/noise 条件的 PTM、survival、conditional fidelity、survival-weighted metric、unitarity、coherent rotation、bootstrap interval 和适用性诊断。

### `rb_fit_metrics.json`

至少包含：

- fit model；
- 数据口径；
- 使用的 lengths/\(M\)；
- 参数和协方差；
- hierarchical bootstrap 区间；
- residual diagnostics；
- AIC/BIC；
- early-window sensitivity；
- loss convention；
- standard RB/IRB assumptions 是否满足；
- 与 PTM estimate 的比较；
- 不可解释时的明确状态。

### `metrics.json`

汇总：

- 上游版本；
- preflight；
- pulse/Clifford/noise 配置；
- channel metrics；
- reference/IRB 结果；
- survival 与 spin 分解；
- site-resolved/array statistics；
- paper-coherence anchor；
- 收敛、消融和敏感性；
- 哪些值来自论文、校准或假设；
- 图片检查状态。

严格 JSON，不含 NaN/Infinity。

## 十九、诊断图

1. `pulse_response.png`：bare pulse 的 Rabi、detuning 和 envelope 响应。
2. `scrofulous_robustness.png`：bare 与 SCROFULOUS 对幅度/失谐误差的 unitary infidelity。
3. `noise_psd_validation.png`：目标 PSD、生成 PSD 和置信带。
4. `ramsey_echo_dd.png`：Ramsey、echo、XY4、XY16 的 contrast。
5. `channel_ptm.png`：idle/transfer/transport/full roundtrip PTM。
6. `channel_bloch_sphere.png`：至少显示三条 Bloch 轴经过条件通道后的变化。
7. `reference_rb_decay.png`：sequence-resolved 散点、均值和多模型拟合。
8. `interleaved_rb_decay.png`：survival、conditional return、joint return 分面。
9. `survival_vs_spin_return.png`：明确显示损失和自旋误差的不同贡献。
10. `rb_residuals.png`：拟合残差、非指数性和 sequence variance。
11. `site_resolved_array.png`：site metrics 分布/空间图和参数相关性。
12. `coherence_vs_round.png`：与 Level 4 相位模型对照。
13. `timestep_convergence.png`：unitary、phase、return probability 和 fidelity 收敛。
14. `sensitivity_panels.png`：Rabi、detuning、\(\eta\)、PSD、相关性、pulse errors。

使用非交互后端，至少 200 dpi，不调用 `plt.show()`，保存后关闭 figure。

## 二十、图片和数据验收

对每张 PNG 程序化检查并写入 `image_validation.json`：存在、可重新打开、尺寸非零、不是全白/全黑/近乎单色、数据和坐标有限、图例数量合理、figure 已关闭。

- 有图像查看工具时逐张打开并记录 `visual_review_performed: true`；
- 没有时记录 `false`，明确只完成程序化检查；
- 未看图时不得声称布局或标签已经视觉确认。

检查 CSV/JSONL 的 schema、键唯一性、sequence 完整性、site/trajectory/noise 对齐、有限值和实际样本数。JSON/YAML/Markdown 必须重新读取验证。

## 二十一、测试要求

保留 Level 0-4 全部测试，并至少新增：

1. Pauli matrices、commutators 和 rotation conventions；
2. 常 Hamiltonian SU(2) analytic propagation；
3. Rabi/Ramsey/echo analytic limits；
4. state norm、density trace、Hermiticity、positivity；
5. bare rotation unitary；
6. pulse envelope area；
7. SCROFULOUS target unitary、branch handling、amplitude-error improvement；
8. DD pulse order 和 finite-pulse propagation；
9. instantaneous-pulse Level 4 limit；
10. trajectory interpolation和边界；
11. differential light shift与 Level 4 phase；
12. OU mean/variance/autocorrelation；
13. PSD normalization、Nyquist 和 reproducibility；
14. common/local spatial correlations；
15. 24 Clifford uniqueness、closure、inverse、Pauli conjugation；
16. Clifford sampler uniformity；
17. ideal sequence + inverse return；
18. known PTM reconstruction；
19. \(F_{\rm avg}\) applicability checks；
20. loss/conditional/joint metrics separation；
21. synthetic reference RB fit recovery；
22. synthetic IRB fit recovery；
23. non-exponential/ill-conditioned fit safe failure；
24. sequence/site/noise hierarchical bootstrap；
25. common random sequence pairing；
26. calibration/validation seed isolation；
27. frozen RB design无数据泄漏；
28. quantum timestep convergence；
29. cached channel与 continuous-noise 小样本一致；
30. 47→195 site extension不改变冻结模型；
31. CLI 各 stage smoke test；
32. CSV/JSON/JSONL schema 和图片程序化检查；
33. 输出位于 `sources/` 时拒绝。

测试不得要求 SCROFULOUS 在所有 detuning error 下优于 bare pulse，也不得要求 RB 必须单指数或 IRB 必须得到论文数值。

## 二十二、执行顺序

1. 检查 Level 0-4 代码、测试、配置和冻结产物；
2. 运行全部已有测试；
3. 建立上游 manifest 和 frozen RB design；
4. 实现 Hamiltonian、传播器、脉冲和 DD；
5. 实现 PSD/noise 和空间相关；
6. 实现 Clifford 群；
7. 实现 PTM/channel characterization；
8. 实现 reference RB、interleaved RB 和稳健拟合；
9. 实现多站点 ensemble；
10. 编写 Level 5 测试；
11. 运行 preflight，失败则修复并重跑；
12. 运行 single-channel characterization；
13. 运行小规模 RB development；
14. 在打开正式结果前冻结 sequences、sites、noise 和拟合窗口；
15. 运行 47-site formal reference/IRB；
16. 运行 195-site extension；
17. 冻结主结果后运行消融和敏感性；
18. 生成全部数据、报告和图片；
19. 程序化检查产物；
20. 若有图像工具，逐图视觉检查并修复；
21. 重跑完整测试和冻结主设计；
22. 检查没有修改 `sources/`；
23. 给出可核查的最终回复。

## 二十三、最终回复必须包含

- 创建/扩展模块和实际成功命令；
- 测试通过、失败、跳过数量；
- 上游路径、hash 和 replay/joint 模式；
- preflight 状态；
- Hamiltonian、pulse family、Rabi frequency、DD 和噪声配置；
- SCROFULOUS 与 bare pulse 的验证结果；
- PSD 验证和所有假设参数；
- idle/transfer/transport/full-roundtrip 的 PTM、survival 和条件通道指标；
- reference RB 拟合、残差和适用性；
- interleaved RB 的 survival、conditional return、joint return 和拟合；
- 标准 IRB estimate 是否可解释；若不可解释，明确说明；
- PTM 与 RB estimate 的一致/不一致；
- 47/195-site 分布、共同/局域噪声贡献；
- Level 4 半经典 contrast 与 Level 5 finite-pulse 结果的差异；
- 收敛、消融和敏感性；
- 哪些数值来自论文、校准或假设；
- 为什么仍不能视为完整实验复现；
- 图片程序化检查及是否真正视觉检查；
- 输出目录和关键文件绝对路径；
- 建议 Level 6 才加入的 Rydberg 两比特门、多体串扰、测量/SPAM、泄漏态和逻辑层电路。

不得只报告一个“fidelity”。必须同时报告 survival、conditional spin channel、joint return、loss convention、拟合假设和区间。
