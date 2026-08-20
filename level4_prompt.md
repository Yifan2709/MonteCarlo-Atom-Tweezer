# Codex 执行 Prompt：Level 4 端到端转移-运输往返与半经典相干性

请在当前工作区中，在已完成并验证的 Level 0 静态势阱、Level 1 一维顺序转移、Level 2 同步转移优化和 Level 3 三维 AOD 长距离运输基础上，实际完成“Level 4：端到端 pick-up—split—long transport—return—merge—drop-off 往返协议与半经典相干性”。

不要只提供方案或代码片段。必须检查并复用 Level 0-3 的实现和冻结产物，创建必要配置、协议编排、测试、模拟、统计、数据和图表，实际运行预检查、单轮验证、重复往返和敏感性分析。

如果工作区存在 `sources/`，其内容均为只读，不得编辑、移动、重命名或删除；如果不存在，直接继续。所有输出路径解析后均不得位于 `sources/`。

论文依据：arXiv:2403.12021v4 的 Fig. 5、Fig. 6、Extended Data Fig. 10、Methods “Atom transport”“Atom transfer between SLM and AOD tweezers”“Combined atom transfer and move”及相关 Supplementary Information。官方版本：https://arxiv.org/abs/2403.12021v4

## 一、Level 4 的目标和边界

前三级分别隔离验证了：

- Level 0：单个静态光镊；
- Level 1：一维顺序 AOD→SLM drop-off；
- Level 2：短距离同步 transfer 波形；
- Level 3：三维 AOD-only 长距离运输与柱面透镜效应。

Level 4 首次将这些组件组合成一个完整往返协议：

1. 原子初始位于静态 SLM；
2. AOD 与 SLM 对齐并 pick up 原子；
3. AOD 与 SLM 分离 2.4 μm；
4. AOD 携带原子作约 375 μm 的对角长距离运输；
5. 在远端保持一个操作时间；
6. AOD 沿反向轨迹返回；
7. AOD 与原始 SLM 合并；
8. 原子 drop off 回 SLM；
9. 重复多个 round trips，统计协议成功、加热和相干对比度。

Level 4 必须完成：

- 组合 Level 2 的冻结 transfer 波形和 Level 3 的冻结三维运输模型；
- 实现论文时序启发的完整 round-trip baseline；
- 对每个阶段设置明确 checkpoint，定位首次失败发生在哪个阶段；
- 区分最终留阱、全协议成功、临时失阱后再俘获；
- 统计多轮往返中的存活和加热累积；
- 实现由轨迹相关差分光移引起的半经典纯退相位；
- 比较无解耦、spin echo 和理想瞬时 XY4 的相干对比度；
- 进行无透镜、无噪声、无差分光移等消融；
- 检查分段协议的连续性、功-能平衡和步长收敛；
- 将经典成功概率、相干对比度和论文的 IRB fidelity 严格区分。

本 Level 4 不模拟完整两能级或多能级 Schrödinger 方程，不模拟 Rabi 驱动脉冲的有限时长和脉冲误差，不实现 Clifford twirling、随机基准测试拟合、Raman 散射、真空碰撞、成像误差或多原子相互作用。因此结果不得称为论文的“transport fidelity”或“IRB fidelity”。

## 二、冻结输入与来源追踪

Level 4 不得重新优化 Level 2 或重新挑选 Level 3 的最优条件。运行前必须发现、加载并冻结：

1. Level 2 的最终 transfer 波形，例如 `best_waveform.yaml`；
2. Level 2 的配置、代码版本、随机种子和 validation 摘要；
3. Level 3 的三维势、AOD lensing 配置和已选定的固定 \(v_s\) 场景；
4. Level 3 的 long-transport profile、步长和 validation 摘要。

对每个输入文件记录绝对路径、SHA-256、生成时间和关键参数，写入 `upstream_manifest.json`。

如果冻结产物不存在：

- 先搜索现有工作区，不要假定固定相对路径；
- 可以从明确的 Level 2/3 配置重建，但必须记录 `reconstructed_from_config: true`；
- 不得暗中重新优化；
- 若无法确定所需波形或 lensing 场景，停止主 validation，输出缺失项和可运行的预检查结果。

Level 4 结果必须能回答：每个 transfer/transport 参数究竟来自论文锚点、Level 2/3 冻结结果，还是本 Level 的显式假设。

## 三、三维总势

继续使用 Level 3 的坐标：

- \(x_1,x_2\)：焦平面两条正交 AOD 轴；
- \(z\)：光轴；
- \(\boldsymbol r=(x_1,x_2,z)\)。

总势为

\[
U_{\rm tot}(\boldsymbol r,t)
=U_{\rm SLM}(\boldsymbol r,t)
+U_{\rm AOD}(\boldsymbol r,t).
\]

### SLM 势

使用标准三维高斯势：

\[
U_{\rm SLM}
=-\frac{D_{\rm SLM}(t)}{1+(z/z_{R,\rm SLM})^2}
\exp\!\left[
-\frac{2[(x_1-X_1)^2+(x_2-X_2)^2]}
{w_{\rm SLM}^2[1+(z/z_{R,\rm SLM})^2]}
\right].
\]

默认 SLM 固定在 \((X_1,X_2)=(0,0)\)，波长 1061 nm，束腰 1.17 μm。SLM 深度随协议变化。

### AOD 势

严格复用 Level 3 的 crossed-AOD lensing 势：

\[
U_{\rm AOD}=-D_{\rm AOD}(t)G_1G_2,
\]

其中 AOD 中心、两轴焦点偏移 \(z_{s,1},z_{s,2}\)、characteristic velocity \(v_s\) 和 lensing signs 均沿用冻结 Level 3 模型。

不得把 Level 2 的一维势与 Level 3 的三维势简单相加两次。短距离 transfer 也必须在同一个三维势实现中执行；在 \(x_2=z=0\)、lensing off 的切片上回归 Level 2。

默认忽略重力。若启用，必须声明实际竖直轴并对所有阶段一致使用。

## 四、完整 round-trip 协议

使用分段协议对象，每段至少返回：

- 起止时间；
- SLM 深度及导数；
- AOD 深度及导数；
- AOD 三维中心、速度、加速度；
- AOD lensing 焦点偏移；
- dynamical-decoupling toggling state；
- checkpoint 类型。

时间轴必须严格递增，每个分段连接处只能保留一个端点样本。除物理上明确的 jerk 跳点外，势、中心和速度必须连续。

### Protocol A：论文时序启发的 anchor

默认完整一轮：

1. `initial_slm_hold`：SLM 180 μK，AOD 关闭，保持 100 μs；
2. `pickup`：850 μs，AOD 从 0 升到 280 μK，同时 SLM 从 180 μK 降到 60 μK，中心重合；
3. `split`：400 μs，AOD 携原子从 SLM 分离 2.4 μm；
4. `outbound_long_move`：1.45 ms，沿对角方向运输 375 μm；
5. `remote_operation_hold`：54 μs，仅模拟在 AOD 中保持，不模拟真实单比特门；
6. `inbound_long_move`：1.45 ms，严格反向返回；
7. `merge`：400 μs，与原始 SLM 合并；
8. `dropoff`：850 μs，AOD 从 280 μK 降到 0，SLM 从 60 μK 升回 180 μK；
9. `final_slm_hold`：SLM 180 μK，保持 100 μs。

这些时长来自论文 Extended Data Fig. 10e 的组合操作锚点。论文没有给出所有深度曲线的逐点数据；默认 ramp 形状必须标为模型假设，不得称为实验精确波形。

默认：

- pickup/dropoff 深度曲线采用有零端点斜率的 monotone smootherstep；
- split/merge 使用 Level 3 已验证的 piecewise constant-jerk 或配置指定轨迹；
- long move 使用 adiabatic-sine；
- split 方向与 long-move 对角方向一致；
- outbound 终点为 AOD 相对起始 SLM 的 \(2.4+375\) μm 路径位置；
- inbound、merge、dropoff 是严格时间反向的控制序列，但原子微观状态不会被人为反转。

### Protocol B：Level 2 + Level 3 composed

- pickup 使用冻结 Level 2 最佳 dropoff 波形的时间反向版本；
- dropoff 使用冻结 Level 2 最佳波形；
- long move 使用冻结 Level 3 指定 profile、lensing 场景和时长；
- split/merge 与 Level 2 波形中已有短移动不得重复计算。

组合前必须检查 Level 2 波形语义：

- 如果 Level 2 dropoff 已包含 2.4 μm merge，则 Protocol B 不再额外添加 merge；
- 时间反向 pickup 必须正确交换起止深度和中心；
- 不得仅反转数组顺序而遗漏导数符号；
- 如果 Level 2 SLM 深度为恒定 140 μK，而 Protocol A 使用 180→60 μK，必须保留为两个不同协议，不能混合后称为同一个 baseline。

### Protocol C：消融控制

至少支持：

- `ideal_no_lensing`；
- `no_differential_light_shift`；
- `no_intensity_noise`；
- `no_dynamical_decoupling`；
- `transport_only`；
- `transfer_only`；
- `static_idle_same_duration`。

每个消融只改变目标机制，其余初态、时间轴和参数保持一致。

## 五、阶段 checkpoint 和成功定义

至少在以下时刻记录 checkpoint：

1. pickup 结束；
2. split 结束；
3. outbound long move 结束；
4. remote hold 结束；
5. inbound long move 结束；
6. merge 结束；
7. dropoff 结束；
8. final hold 结束。

### 局部束缚判据

在单阱主导的 checkpoint，分别计算原子相对 AOD 或 SLM 静态瞬时势阱的能量：

\[
E_{\rm trap}
=\frac12m|\boldsymbol v-\boldsymbol v_{\rm trap}|^2
+U_{\rm local}(\boldsymbol r)-U_{\rm continuum}.
\]

在 trap 静止的 checkpoint，\(\boldsymbol v_{\rm trap}=0\)。只有 \(E_{\rm trap}<0\) 才判为局部束缚。

在 AOD/SLM 重叠且两个势阱都显著存在时，不能用“分别对两个阱算负能量”简单决定归属。必须通过以下一种经过测试的方法分类：

- 找出总势的局部最低点和分隔势垒，将原子分配给相空间上可达的 basin；
- 或在短暂冻结控制的诊断副本中推进若干振动周期确定归属。

分类必须输出 `bound_to_slm`、`bound_to_aod`、`shared_or_ambiguous`、`unbound`。无法稳定归类时标为 ambiguous，不得强行判成功。

### 三类协议结果

- `final_retained`：final hold 后 \(E_{\rm SLM}<0\)；
- `roundtrip_success`：所有必需 checkpoint 均满足预期阱归属，且 final_retained；
- `absorbing_success`：一旦某 checkpoint 失败，该 shot 在后续轮次中永久记为失败，即使动力学上偶然再俘获；
- `recaptured`：中间 checkpoint 失败但最终又回到 SLM。

多轮 survival 曲线主指标使用 `absorbing_success`，同时单独报告非吸收的最终留阱和 recapture。不得静默删除失败粒子或在每轮重新热化/重新采样。

### 首次失败归因

为每个 shot 记录 first-failure stage：

`pickup`、`split`、`outbound_transport`、`remote_hold`、`inbound_transport`、`merge`、`dropoff`、`final_hold` 或 `none`。

归因是基于模型 checkpoint 的操作性定义，不等于实验中的真实损失机制。

## 六、初态和重复往返

### 初态

原子初始位于 180 μK SLM。使用 Level 3 的三维 harmonic-proposal + full-potential bound rejection 方法，按默认 5 μK 采样六维初态。

保存接受率、均值、协方差和种子。明确说明该采样不是有限深高斯阱的严格正则系综。

### Common random numbers

所有协议和消融使用相同初态 ID。至少隔离：

- `development_pool`：预检查和调试；
- `single_round_validation_pool`；
- `repeated_round_pool`；
- `sensitivity_pool`。

### 多轮推进

一个 shot 在完成一轮后，必须把其真实 final state 作为下一轮 initial state；不得重新抽样、重置到阱中心或人为冷却。

默认：

- 单轮 validation：2,000 shots；
- 重复往返：500 shots，连续推进至 10 rounds；
- 长尾探索：128 shots，较粗步长推进至 40 rounds；
- 在第 1、2、5、10、20、40 轮保存 checkpoint；如果最大轮数较小，只保存可用项。

一次运行到最大轮数并在中间保存前缀结果，不得为每个 \(N\) 从头重复积分。

对 absorbing-failed shots 可以停止详细积分以节省计算，但必须保留其失败阶段和时间；若还要估计 recapture，则需用单独非吸收子样本继续积分，不能混淆两种数据。

## 七、经典加热和功-能账本

对每个分段计算：

- 初末机械能；
- SLM 深度变化功；
- AOD 深度变化功；
- 横向中心移动功；
- AOD lensing 焦点偏移功；
- 总外部功；
- 功-能残差；
- 相对当前目标阱底的激发。

显式时变阶段满足

\[
\Delta E=W_{\rm ext}+R_W.
\]

静态 hold 段应只有有界数值能量误差。

对完整 round trip 输出 segment-resolved energy ledger。总残差不能仅因正负分段抵消而通过；每段和全周期都必须检查。

加热至少用以下指标：

- checkpoint 激发 \(\varepsilon=E-E_{\min}\)；
- 每轮结束 SLM 激发；
- retained/success shots 的均值、中位数、90% 和 95% 分位数；
- 径向与轴向能量代理；
- 随轮数的 heating slope，但只有在线性拟合残差合理时才报告线性率。

未留阱粒子的正能量不得和留阱态激发混在一个均值中。

## 八、半经典差分光移与相干对比度

### 模型定位

只追踪每条经典轨迹上两个超精细钟态之间的相对相位，不求解完整自旋动力学。定义差分光移角频率

\[
\delta\omega_i(t)
=\frac{\eta_{\rm SLM}U_{\rm SLM}[\boldsymbol r_i(t),t]
+\eta_{\rm AOD}U_{\rm AOD}[\boldsymbol r_i(t),t]}{\hbar}.
\]

默认 nominal：

\[
\eta_{\rm SLM}=\eta_{\rm AOD}=1.3\times10^{-4},
\]

并以 \(1.5\times10^{-4}\) 作为理论/敏感性对照。符号和两个波长是否完全相同必须设为可配置；不得用一个无依据的精确共同值掩盖波长差异。

### Dynamical-decoupling toggling function

用 \(y(t)=\pm1\) 表示理想瞬时 \(\pi\) 脉冲对纯退相位项的符号翻转：

\[
\phi_i(T)=\int_0^T y(t)\,\delta\omega_i(t)\,\mathrm dt.
\]

至少实现：

- `none`：\(y(t)=1\)；
- `spin_echo`：协议中点一个理想 \(\pi\) 脉冲；
- `xy4_per_long_move`：每个 outbound/inbound long move 内放置一个对称 XY4 序列；
- `custom_pulse_times`：由配置给出。

在仅有理想纯退相位和瞬时脉冲的模型中，XY4 的 X/Y 轴顺序对 toggling function 等价；不得因此声称已经模拟了 XY4 对脉冲误差或高阶噪声的完整抑制。

脉冲时刻不得恰好落在重复时间样本或分段边界的歧义位置。积分必须精确处理 toggling discontinuity，必要时把时间网格在脉冲处拆分。

### 对比度

对一组 shots 定义条件相干对比度

\[
C_{\rm cond}
=\left|\frac{1}{N_{\rm succ}}\sum_{i\in\rm success}e^{i\phi_i}\right|.
\]

同时报告：

- `contrast_all_final_retained`；
- `contrast_roundtrip_success`；
- `phase_mean_direction`；
- `circular_phase_std`；
- `usable_coherent_fraction=P(roundtrip_success)\times C_cond`。

最后一个量只是工程代理，不是量子过程 fidelity。

### 可选技术噪声

默认 nominal 运行不加入随机技术噪声。敏感性模式可加入：

- 每个 shot 的 quasistatic fractional intensity offset；
- Ornstein-Uhlenbeck fractional intensity noise；
- SLM 与 AOD 独立或相关噪声。

所有 RMS、相关时间和相关系数必须来自配置并标记 `assumed`，不得声称由论文唯一确定。噪声使用独立 seed，且与初态 seed 分离。

不包含磁场噪声、Raman scattering、pulse-area error、有限脉冲时长或 Clifford gate error。remote 54 μs 只作为相位积累 hold，不得输出 gate fidelity。

## 九、重复轮次统计和拟合

每轮至少输出：

- risk set 大小；
- absorbing success 数和概率；
- non-absorbing final retention；
- first-failure stage 分布；
- retained shots 的激发分布；
- 条件相干对比度；
- usable coherent fraction；
- near-threshold 和 ambiguous 比例。

使用 Wilson 95% 区间估计二项成功概率。协议比较使用同初态配对 bootstrap，并报告四格配对计数。

可以拟合论文启发的经验生存曲线，例如 clipped-Boltzmann 或 exponential × clipped-Boltzmann，但必须：

- 至少有足够的非平凡轮次；
- 同时比较简单指数和非参数结果；
- 报告参数区间、残差和 AIC/BIC；
- 不把拟合参数直接解释成实验温度或每次量子通道 fidelity；
- 如果数据几乎全为 0 或 1、参数不可辨识或协方差奇异，跳过拟合并说明。

不得把 \(S_{n+1}/S_n\) 在小样本尾部的剧烈波动称为可靠瞬时 fidelity。

## 十、运行矩阵

### Stage A：preflight

只使用理想初态和小样本，完成所有回归和连续性检查。

### Stage B：single-round validation

运行前冻结：

1. Protocol A，Level 3 nominal nonzero lensing；
2. Protocol A，lensing off；
3. Protocol B，冻结 Level 2+3 参数；
4. static idle same total duration；
5. 每个协议分别运行 no DD 与 XY4-per-long-move。

每个条件至少 2,000 个共同初态、正式步长 0.05 μs。使用独立 validation pool，不得根据结果重新优化。

### Stage C：repeated rounds

- Protocol A 和 B；
- nominal lensing；
- none 与 XY4；
- 500 shots、10 rounds、0.05 或经过收敛批准的 0.10 μs；
- 128-shot 长尾至 40 rounds。

如果使用 0.10 μs 做多轮主结果，必须在固定子集上证明其 survival、heating、phase 和 work residual 与 0.05/0.025 μs 一致到预设容差。

### Stage D：消融

至少比较：

- transfer-only；
- transport-only；
- no lensing；
- no differential light shift；
- no intensity noise；
- static idle；
- 删除 remote hold；
- long move 改用 constant-jerk。

消融用于归因，不得从消融中挑选一个更好协议后重新命名为预注册主协议。

### Stage E：敏感性

冻结主协议后检查：

- 初态温度：3、5、10 μK；
- SLM 初始/最终深度：180 μK ±10%；
- SLM 低深度：60 μK ±10%；
- AOD 深度：280 μK ±10%；
- SLM/AOD 相对对准：横向每轴 0、±0.05、±0.10 μm；
- lensing severity：Level 3 固定 nominal 的相邻场景；
- \(\eta\)：1.2、1.3、1.5、1.6 × \(10^{-4}\)；
- intensity noise RMS 和 correlation time 的明确假设网格；
- DD pulse timing offset。

默认每个条件 500 shots，以 one-factor-at-a-time 为主；联合最坏情形必须另行标记。

## 十一、建议配置

在不破坏现有配置系统的前提下新增 Level 4 配置，例如：

```yaml
level4:
  name: end_to_end_roundtrip_coherence

upstream:
  level2_best_waveform: auto_discover
  level3_metrics: auto_discover
  require_sha256_manifest: true
  allow_reconstruction_without_reoptimization: true

atom:
  isotope: Cs-133
  mass_u: 132.90545196

slm:
  wavelength_nm: 1061.0
  waist_um: 1.17
  initial_depth_uK: 180.0
  transport_depth_uK: 60.0
  center_um: [0.0, 0.0]

aod:
  wavelength_nm: 1055.0
  waist_um: 1.17
  transport_depth_uK: 280.0
  split_distance_um: 2.4
  inherit_lensing_from_level3: true

protocol_anchor:
  initial_hold_us: 100.0
  pickup_us: 850.0
  split_us: 400.0
  long_distance_um: 375.0
  outbound_us: 1450.0
  remote_operation_hold_us: 54.0
  inbound_us: 1450.0
  merge_us: 400.0
  dropoff_us: 850.0
  final_hold_us: 100.0
  long_geometry: diagonal
  long_profile: adiabatic_sine
  split_profile: constant_jerk

initial_ensemble:
  temperature_uK: 5.0
  sampler: harmonic_3d_rejection
  single_round_validation_shots: 2000
  repeated_round_shots: 500
  long_tail_shots: 128
  single_round_seed: 41001
  repeated_round_seed: 41002
  long_tail_seed: 41003
  sensitivity_seed: 41004

integration:
  validation_dt_us: 0.05
  repeated_dt_us: 0.10
  convergence_dt_us: 0.025
  max_rounds: 10
  long_tail_max_rounds: 40

coherence:
  enabled: true
  eta_slm: 1.3e-4
  eta_aod: 1.3e-4
  dd_modes: [none, spin_echo, xy4_per_long_move]
  pulses_are_instantaneous: true
  include_magnetic_noise: false
  include_pulse_errors: false

intensity_noise:
  enabled_nominal: false
  model: ornstein_uhlenbeck
  slm_fractional_rms: 0.0
  aod_fractional_rms: 0.0
  correlation_time_us: null
  cross_correlation: 0.0
  noise_seed: 41999
  provenance: assumed_sensitivity_only

classification:
  energy_threshold_J: 0.0
  near_threshold_fraction_of_depth: 0.001
  ambiguous_basin_tolerance: 0.001
  absorbing_protocol_success: true

output:
  directory: outputs/level4_end_to_end_roundtrip/demo
  save_all_checkpoint_states: true
  save_representative_trajectories: true
  representative_shots_per_failure_stage: 3
```

如已有 Level 0-3 使用不同等价结构，应整合到现有配置，不要新建互不兼容的系统。`config_used.yaml` 必须保存解析后的上游路径、hash、CLI 覆盖、所有随机种子和冻结条件。

## 十二、代码结构与复用

必须复用现有三维势、积分器、Level 2 波形、Level 3 lensing、热采样、统计、输出和图片检查。建议职责，可合理合并：

```text
protocol_segments.py          # 分段时序和连续拼接
upstream_artifacts.py         # 发现、hash 和冻结 Level 2/3 产物
trap_assignment.py            # 重叠势 basin 与 checkpoint 归属
roundtrip_simulation.py       # 单轮和多轮向量化推进
energy_ledger.py              # 分段功-能账本
dephasing.py                  # differential shift、toggling 和 phase
dd_sequences.py               # none、echo、XY4、custom pulse times
noise_models.py               # 可选 quasistatic/OU 强度噪声
roundtrip_statistics.py       # Wilson、paired bootstrap、生存曲线
level4_visualization.py
level4_cli.py
```

所有公开函数包含简短中文 docstring。协议编排与底层动力学必须解耦；不得为 Protocol A/B 复制两套积分器。

CLI 示例：

```bash
python -m level0_static_trap.level4_cli \
  --config configs/level4_roundtrip.yaml \
  --stage all \
  --output-dir outputs/level4_end_to_end_roundtrip/demo
```

如果实际包名不同，使用现有包名。`--stage` 至少支持：

```text
preflight
single_round
repeated_rounds
ablations
sensitivity
all
```

不得使用临时 `PYTHONPATH` 掩盖打包问题。

## 十三、preflight：失败则停止大样本

按顺序执行并写入 `preflight.json`：

1. 运行全部 Level 0-3 测试；
2. 发现并校验上游冻结产物、hash 和 schema；
3. 检查每个 segment 的起止值、导数、单位和时间连续性；
4. 检查完整 round-trip 的 AOD/SLM 最终状态严格回到初始控制状态；
5. 在三维切片和 lensing-off 极限复现 Level 2 transfer；
6. 单独运行 long move 并复现 Level 3 冻结条件；
7. 检查 Protocol B 没有重复 split/merge；
8. 对理想中心零速度初态完成整个 Protocol A/B；
9. 在无噪声、lensing-off 条件下进行控制时间反向测试：将终态速度反转并反演控制，回到初态误差随步长收敛；
10. 对每个 segment 和全周期检查功-能残差；
11. 在静态 hold 段检查能量误差有界；
12. 检查 overlap basin 分类在构造势上正确，并能返回 ambiguous；
13. 检查 checkpoint 和 first-failure 状态机；
14. 检查单轮末态正确传给下一轮，不重采样；
15. 用常数差分光移验证相位积分解析值；
16. 对对称 quasistatic shift 验证 ideal echo/XY4 toggling cancellation；
17. 检查 DD pulse discontinuity 的时间网格处理；
18. 检查 ensemble contrast 的全同相、均匀相位和已知两点相位案例；
19. 对至少 128 shots 做 0.10/0.05/0.025 μs 的状态、成功标签、功-能和 phase 收敛；
20. 冻结主 validation 条件并确认 validation 数据未被调参代码访问。

任何关键项失败时停止主 validation，输出原因。不得带着失败 preflight 继续生成“完整 Level 4”结论。

## 十四、输出文件

默认至少生成：

```text
outputs/level4_end_to_end_roundtrip/demo/
├── config_used.yaml
├── upstream_manifest.json
├── frozen_validation_protocols.yaml
├── preflight.json
├── protocol_timeline.csv
├── single_round_shots.csv
├── checkpoint_states.csv
├── repeated_round_metrics.csv
├── failure_events.csv
├── energy_ledger.csv
├── phase_records.csv
├── ablations.csv
├── sensitivity.csv
├── representative_trajectories.csv
├── metrics.json
├── image_validation.json
├── summary.md
├── protocol_timeline.png
├── trap_path_and_depths.png
├── representative_roundtrips.png
├── checkpoint_survival.png
├── failure_stage_breakdown.png
├── survival_vs_round.png
├── heating_vs_round.png
├── phase_and_contrast.png
├── dd_toggling_and_phase.png
├── energy_ledger.png
├── ablation_comparison.png
├── timestep_convergence.png
└── sensitivity_panels.png
```

### `protocol_timeline.csv`

保存 time、segment、SLM/AOD 深度、AOD 中心/速度/加速度、lensing shifts、DD toggling state 和 pulse markers。

### `single_round_shots.csv`

每个 shot/协议/DD 条件一行，至少包含初态、各 checkpoint 归属、first-failure、final_retained、roundtrip_success、recaptured、末态、激发、功-能残差、phase 和 contrast grouping keys。

### `checkpoint_states.csv`

长格式保存每个 shot 在每个 checkpoint 的六维状态、局部能量、basin、near-threshold、ambiguous 和 cumulative work。

### `repeated_round_metrics.csv`

每轮保存 risk set、absorbing/non-absorbing success、Wilson 区间、failure-stage counts、激发分位数、contrast 及 bootstrap 区间。

### `phase_records.csv`

每个 shot/round/DD 条件保存各 segment phase、总 phase、survival conditioning 和噪声 realization ID。不得只保存 ensemble 平均而丢失相位分布。

### `metrics.json`

至少包含：

- 上游 manifest 和 hash；
- 全部 preflight；
- Protocol A/B 每段参数与来源；
- 冻结 validation 条件；
- 单轮 final retention、roundtrip success、recapture 和 Wilson 区间；
- first-failure 分布；
- 协议/DD 的配对差异与 bootstrap 区间；
- 多轮 survival、heating 和相干对比度；
- phase 分布、usable coherent fraction；
- 功-能账本和步长收敛；
- 消融和敏感性；
- 经验生存拟合的适用性和诊断；
- 可与论文比较/不可比较的量；
- 图片程序化与视觉检查状态。

JSON 不得包含 NaN 或 Infinity。

## 十五、诊断图

1. `protocol_timeline.png`：完整分段时间线与 checkpoint。
2. `trap_path_and_depths.png`：SLM/AOD 深度、AOD 路径、速度和 lensing shifts。
3. `representative_roundtrips.png`：成功、各主要失败阶段、recaptured 和 ambiguous 代表轨迹；不得只画成功轨迹。
4. `checkpoint_survival.png`：每个 checkpoint 的条件/累计成功概率。
5. `failure_stage_breakdown.png`：first-failure stage 堆叠图。
6. `survival_vs_round.png`：absorbing 与 non-absorbing survival，Wilson 区间和可选拟合。
7. `heating_vs_round.png`：仅成功/留阱 shots 的激发分布随轮次。
8. `phase_and_contrast.png`：相位分布、none/echo/XY4 contrast 随轮次。
9. `dd_toggling_and_phase.png`：代表轨迹的 \(y(t)\)、差分光移和累计相位。
10. `energy_ledger.png`：各分段外部功、能量变化和残差。
11. `ablation_comparison.png`：transfer、transport、lensing、DD、idle 消融。
12. `timestep_convergence.png`：状态、功-能、success、phase 和 contrast 收敛。
13. `sensitivity_panels.png`：温度、深度、对准、\(\eta\)、噪声和 pulse timing。

使用非交互后端，至少 200 dpi，不调用 `plt.show()`，保存后关闭 figure。

## 十六、图片和数据验收

对每张 PNG 做程序化检查并写入 `image_validation.json`：存在、可重新打开、尺寸非零、不是全白/全黑/近乎单色、绘图数据有限、坐标范围非零、图例数量合理、figure 已关闭。

- 若环境提供图像查看工具，逐张打开并记录 `visual_review_performed: true`；
- 若没有，记录 `false`，明确说明只做程序化检查；
- 未真正看图不得声称布局清晰或标签无重叠。

对 CSV 检查 schema、行数、shot/round/checkpoint 键唯一性、共同初态对齐、协议格点完整性和有限值。对 JSON/YAML/Markdown 重新读取验证。

## 十七、测试要求

保留 Level 0-3 全部测试，并至少新增：

1. 协议分段端点与导数连续；
2. 完整控制 round trip 回到初始 SLM/AOD 状态；
3. Protocol B 不重复 split/merge；
4. Level 2 transfer 三维回归；
5. Level 3 long-move 回归；
6. 控制时间反向和速度反转测试；
7. 总势、总力和各组件力一致；
8. overlap basin 分类和 ambiguous；
9. 各 checkpoint 局部束缚判据；
10. first-failure 状态机；
11. recaptured 与 absorbing failure 区分；
12. 多轮状态连续传递且不重采样；
13. 每段及总周期功-能关系；
14. 静态 hold 能量误差有界；
15. constant differential shift 的解析 phase；
16. toggling function 在 none/echo/XY4 下正确；
17. pulse 时刻网格拆分无重复或漏积分；
18. ensemble contrast 构造案例；
19. survival conditioning 不混入失败 shots；
20. quasistatic/OU noise 的均值、方差、相关时间和 seed 复现；
21. initial-state seed 与 noise seed 隔离；
22. common random numbers 对齐；
23. validation 冻结和无数据泄漏；
24. Wilson、paired bootstrap、circular statistics；
25. 生存拟合在不可辨识数据上安全失败；
26. 0.10/0.05/0.025 μs 状态、work、phase 二阶收敛；
27. near-threshold 和 ambiguous 不被删除；
28. CLI 各 stage smoke test；
29. CSV/JSON schema 和图片程序化检查；
30. 输出位于 `sources/` 时被拒绝。

测试不得要求 Protocol B 必须优于 Protocol A，也不得要求 XY4 在所有噪声模型下必然改善。应测试实现和统计是否正确，而不是把预期物理结论写死。

## 十八、执行顺序

1. 检查现有 Level 0-3 代码、测试、配置和冻结输出；
2. 运行全部已有测试；
3. 建立 upstream manifest 并冻结输入；
4. 实现协议分段、连续拼接和 checkpoint；
5. 实现总势、basin 分类和多轮状态机；
6. 实现分段功-能账本；
7. 实现 differential-shift phase、DD toggling 和可选噪声；
8. 编写 Level 4 测试；
9. 运行 preflight，失败则修复并重跑；
10. 在打开 validation 前写入 `frozen_validation_protocols.yaml`；
11. 运行独立 2,000-shot 单轮 validation；
12. 运行 500-shot/10-round 重复协议；
13. 运行 128-shot/40-round 长尾探索；
14. 冻结主结果后运行消融和敏感性；
15. 生成数据、报告和图；
16. 程序化检查全部产物；
17. 若有图像工具，逐图视觉检查并修复；
18. 重新运行完整测试和冻结主协议，确保结果对应最新代码；
19. 检查未修改 `sources/`；
20. 给出可核查的最终回复。

## 十九、最终回复必须包含

- 实际创建/扩展的模块和成功命令；
- 测试通过、失败和跳过数量；
- upstream 文件路径、hash 和是否重建；
- preflight 状态；
- Protocol A/B 的完整时序及假设；
- 单轮 final retention、roundtrip success、recapture、Wilson 区间；
- first-failure stage 分布；
- Protocol A/B 及消融的配对差异；
- 多轮 survival 和 heating；
- none/echo/XY4 的条件相干对比度、phase 分布和 usable coherent fraction；
- 功-能残差和步长收敛；
- 温度、深度、对准、lensing、\(\eta\)、噪声和 pulse timing 敏感性；
- 哪些输入来自论文、Level 2/3 或假设；
- 为什么结果不能等同论文 IRB fidelity；
- 图片程序化检查和是否真正视觉检查；
- 输出目录和关键文件绝对路径；
- 建议 Level 5 才加入的有限脉冲、自旋哈密顿量、噪声谱、Clifford/IRB 和多原子并行操作。

不得只汇报最终成功率；不得隐藏中间失败、recapture、ambiguous、相位后选择效应或不利的消融结果。
