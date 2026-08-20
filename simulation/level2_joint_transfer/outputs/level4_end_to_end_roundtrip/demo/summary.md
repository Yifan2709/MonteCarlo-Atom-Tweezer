# Level 4 端到端往返协议与半经典相干性

组合 Level 2 冻结 transfer 波形与 Level 3 冻结三维 AOD 运输，完成 pick-up—split—长运输—返回—merge—drop-off 往返协议的经典蒙特卡洛验证与半经典差分光移相位分析。

- Protocol A 总时长 5654 μs；Protocol B 总时长 3954 μs。
- nominal v_s = 0.5392 m/s（Level 3 severity 1.0 冻结场景）。


## 单轮 validation（validation_pool，2000 shots，dt=0.05 μs）

| 条件 | roundtrip success | Wilson 95% CI | final retained | recaptured |
|---|---:|---|---:|---:|
| protocol_a_nominal_lensing | 2000/2000 (1.0000) | [0.9981, 1.0000] | 1.0000 | 0 |
| protocol_a_lensing_off | 2000/2000 (1.0000) | [0.9981, 1.0000] | 1.0000 | 0 |
| protocol_b_frozen | 2000/2000 (1.0000) | [0.9981, 1.0000] | 1.0000 | 0 |
| static_idle_same_duration | 2000/2000 (1.0000) | [0.9981, 1.0000] | 1.0000 | 0 |

### 配对比较（同初态）
- lensing_off_vs_nominal: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0）。
- protocol_b_vs_a: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0）。
- static_idle_vs_a: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0）。

### 相干对比度（条件于 roundtrip success）
- protocol_a_nominal_lensing|none: C=0.8797，mean direction=-0.364 rad，circular std=0.5062。
- protocol_a_nominal_lensing|spin_echo: C=1.0000，mean direction=-0.001 rad，circular std=0.0091。
- protocol_a_nominal_lensing|xy4_per_long_move: C=0.9695，mean direction=+0.450 rad，circular std=0.2488。
- protocol_a_lensing_off|none: C=0.8798，mean direction=-0.366 rad，circular std=0.5062。
- protocol_a_lensing_off|spin_echo: C=1.0000，mean direction=-0.000 rad，circular std=0.0024。
- protocol_a_lensing_off|xy4_per_long_move: C=0.9695，mean direction=+0.449 rad，circular std=0.2490。
- protocol_b_frozen|none: C=0.9308，mean direction=+1.126 rad，circular std=0.3788。
- protocol_b_frozen|spin_echo: C=1.0000，mean direction=-0.001 rad，circular std=0.0068。
- protocol_b_frozen|xy4_per_long_move: C=0.9954，mean direction=+1.913 rad，circular std=0.0958。
- static_idle_same_duration|none: C=0.9165，mean direction=+2.228 rad，circular std=0.4177。
- static_idle_same_duration|spin_echo: C=1.0000，mean direction=-0.000 rad，circular std=0.0015。
- static_idle_same_duration|xy4_per_long_move: C=0.9165，mean direction=+2.228 rad，circular std=0.4177。

## 功-能账本（Protocol A 代表批）

- initial_slm_hold: ΔE 均值 -154.897 μK，W_ext 均值 -154.896 μK，max|R| 2.65e-03 μK。
- pickup: ΔE 均值 +58.252 μK，W_ext 均值 +58.252 μK，max|R| 2.53e-03 μK。
- split: ΔE 均值 -0.017 μK，W_ext 均值 -0.017 μK，max|R| 2.22e-03 μK。
- outbound_long_move: ΔE 均值 -0.000 μK，W_ext 均值 -0.000 μK，max|R| 1.96e-03 μK。
- remote_operation_hold: ΔE 均值 +0.016 μK，W_ext 均值 +0.015 μK，max|R| 2.35e-03 μK。
- inbound_long_move: ΔE 均值 -58.250 μK，W_ext 均值 -58.250 μK，max|R| 3.17e-03 μK。
- merge: ΔE 均值 +154.895 μK，W_ext 均值 +154.895 μK，max|R| 3.98e-03 μK。
- dropoff: ΔE 均值 +0.000 μK，W_ext 均值 +0.000 μK，max|R| 1.67e-03 μK。
- final_slm_hold: ΔE 均值 +0.000 μK，W_ext 均值 +0.000 μK，max|R| 0.00e+00 μK。

## 诊断图清单

- `protocol_timeline.png`：Protocol A/B 分段时间线、checkpoint 与 DD toggling/脉冲
- `trap_path_and_depths.png`：Protocol A/B 深度、AOD 中心/速度与 lensing 焦点偏移
- `initial_ensemble.png`：初态六维采样分布与接受率/温度/种子元数据
- `upstream_provenance.png`：上游冻结输入（路径、SHA-256、时间）一览
- `preflight_checks.png`：20 项 preflight 状态与关键数值诊断
- `representative_roundtrips.png`：成功/各失败阶段/recaptured 代表轨迹（位置+能量+功）
- `checkpoint_survival.png`：四个条件的逐 checkpoint 条件/累计成功概率
- `basin_classification.png`：各 checkpoint basin 归属、重叠区 E_SLM/E_AOD 散点
- `checkpoint_energy_evolution.png`：checkpoint 总能量分位数带随协议阶段演变
- `failure_stage_breakdown.png`：单轮 first-failure 阶段堆叠（按条件）
- `segment_phase_breakdown.png`：各 DD 模式的分段相位贡献分解
- `paired_comparisons.png`：配对四格计数与 bootstrap 差异区间
- `survival_vs_round.png`：10 轮 absorbing/非吸收生存与拟合
- `long_tail_survival.png`：128-shot 长尾（至 40 轮）absorbing 生存曲线
- `failure_events_by_round.png`：重复往返每轮新增失败事件的阶段堆叠 （本次未生成：全程无 first-failure 事件，无可绘数据）
- `heating_vs_round.png`：留阱 shots 的 SLM 激发分位数随轮次
- `phase_and_contrast.png`：相位分布与各 DD 模式对比度随轮次
- `usable_coherent_fraction.png`：usable_coherent_fraction = P(success)×C_cond 随轮次
- `dd_toggling_and_phase.png`：代表轨迹的 y(t)、δω 分解与三种 DD 累计相位
- `energy_ledger.png`：分段功分解、ΔE–W_ext 对照与残差
- `noise_model_diagnostic.png`：assumed 强度噪声示例轨迹与统计
- `ablation_comparison.png`：消融的成功率、激发与 C(none)/C(xy4)
- `timestep_convergence.png`：状态/功/相位/成功率随步长收敛
- `sensitivity_panels.png`：温度、深度、对准、lensing、η、噪声与脉冲时刻敏感性

## 与论文的关系与边界

- 本级输出的是经典协议成功概率（checkpoint basin 归属 + 末态E_SLM<0）与半经典轨迹相位对比度，**不得**称为论文的 transport fidelity 或 IRB fidelity：未包含 Clifford twirling、随机基准拟合、Raman 散射、真空碰撞、脉冲误差与成像误差。
- usable_coherent_fraction = P(roundtrip_success)×C_cond 只是工程代理，不是量子过程保真度。
- remote 54 μs 仅作为相位积累 hold，不输出任何 gate fidelity。
- 深度 ramp 形状（smootherstep）与 η=1.3e-4 为模型假设；噪声参数标记 assumed。

## Level 5 建议

- 有限时长/含脉冲误差的 Rabi 驱动与完整自旋哈密顿量；
- 真实噪声谱（磁场、强度 OU 谱）与 pulse-area error；
- Clifford/IRB 随机基准协议拟合；
- 多原子并行操作的相互作用与成像误差。

图片完成程序化检查，并经逐图视觉审阅（23 张，zcode-visual-review）：全部通过。
