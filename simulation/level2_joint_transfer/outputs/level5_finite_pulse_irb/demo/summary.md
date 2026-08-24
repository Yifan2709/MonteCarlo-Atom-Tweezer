# Level 5 有限脉冲单比特动力学、噪声谱与多站点 RB/IRB

在 Level 4 冻结的经典轨迹与协议之上，把内部态升级为旋转框架下有限时长微波脉冲驱动的 ^133Cs 钟态量子比特，完成单通道 PTM 表征、reference RB、transport-interleaved RB 与多站点统计。

- Rabi 锚点 Ω₀ = 2π×24.611 kHz（论文 Fig. 4a）；η_SLM=η_AOD=0.00013；nominal 噪声关。
- 脉冲族 scrofulous（Clifford 门）/ bare（DD π 脉冲）；DD transport 模式 xy4_per_long_move。
- Clifford 群 24 元素；平均 SCROFULOUS 面积 2.280π（论文平均锚点 2.02π，分解不同不作归一）。
- preflight：26 项全部通过（1220 s）。

## 单通道表征（survival-conditioned PTM）

| 条件 | DD | survival | F_avg(cond) | 95% CI | F_avg(joint) |
|---|---|---:|---:|---|---:|
| static_idle | none | 1.0000 | 0.476391 | [0.47275, 0.48013] | 0.476391 |
| static_idle | spin_echo | 1.0000 | 0.333523 | [0.33352, 0.33352] | 0.333523 |
| static_idle | xy4 | 1.0000 | 0.999998 | [1.00000, 1.00000] | 0.999998 |
| static_idle | xy16 | 1.0000 | 0.999996 | [1.00000, 1.00000] | 0.999996 |
| transfer_only | xy4 | 1.0000 | 0.993590 | [0.99358, 0.99360] | 0.993590 |
| transport_only | none | 1.0000 | 0.806914 | [0.79780, 0.81607] | 0.806914 |
| transport_only | xy4_per_long_move | 1.0000 | 0.562426 | [0.56011, 0.56474] | 0.562426 |
| transport_only | xy16 | 1.0000 | 0.999989 | [0.99999, 0.99999] | 0.999989 |
| protocol_a | none | 1.0000 | 0.943190 | [0.94051, 0.94561] | 0.943190 |
| protocol_a | spin_echo | 1.0000 | 0.333380 | [0.33338, 0.33338] | 0.333380 |
| protocol_a | xy4_per_long_move | 1.0000 | 0.957305 | [0.95536, 0.95915] | 0.957305 |
| protocol_a | xy16 | 1.0000 | 0.995996 | [0.99599, 0.99600] | 0.995996 |
| protocol_b | none | 1.0000 | 0.795773 | [0.78595, 0.80536] | 0.795773 |
| protocol_b | xy4_per_long_move | 1.0000 | 0.553874 | [0.55138, 0.55630] | 0.553874 |
| protocol_a_no_lensing | xy4_per_long_move | 1.0000 | 0.956368 | [0.95217, 0.96037] | 0.956368 |
| protocol_a | xy4_per_long_move | 1.0000 | 0.939439 | [0.93636, 0.94228] | 0.939439 |
| protocol_a | xy4_per_long_move | 1.0000 | 0.956798 | [0.95473, 0.95871] | 0.956798 |

## Reference RB（47 sites × 72 序列）

- 单指数 p = 0.9967100824486766；首选模型 stretched；早期斜率 -0.0035880310423875997。

## Interleaved RB（N=80，47 sites × 72 strings）

- IRB 判定：IRB estimate interpretable under standard assumptions
- 标准 IRB 每移动 F_avg = 1.001650（p_ref=0.996710, p_comb=1.000000）
- survival（M=80）= 1.0000；conditional（M=80）= 0.8785；joint（M=80）= 0.8785。
- 相邻比值稳定性（conditional）：True；尾部 risk set = None。

## sites_47

- array conditional = 0.87465；between-site std = 0.00106；within-site std = 0.11385；worst decile = 0.87332。
- 站点参数相关性：{}。

## sites_195

- array conditional = 0.87832；between-site std = 0.00065；within-site std = 0.10732；worst decile = 0.87739。
- 站点参数相关性：{}。

## 敏感性（one-factor-at-a-time，assumed）

- baseline = nominal: F_avg(cond)=0.957410
- rabi_scale = 0.98: F_avg(cond)=0.959493
- rabi_scale = 1.02: F_avg(cond)=0.955322
- delta_drive_hz = -50.0: F_avg(cond)=0.962949
- delta_drive_hz = 50.0: F_avg(cond)=0.747883
- eta = 0.00012: F_avg(cond)=0.726995
- eta = 0.00015: F_avg(cond)=0.717968
- eta = 0.00016: F_avg(cond)=0.445098
- quasistatic_detuning_rms_hz = 10.0: F_avg(cond)=0.957386
- quasistatic_detuning_rms_hz = 50.0: F_avg(cond)=0.957436
- ou_rms_hz_tau200us = 5.0: F_avg(cond)=0.957591
- ou_rms_hz_tau200us = 20.0: F_avg(cond)=0.954825
- pulse_amp_error = 0.01: F_avg(cond)=0.956366
- pulse_amp_error = 0.03: F_avg(cond)=0.954280
- envelope = cosine_edge_1us: F_avg(cond)=0.957412
- dd_mode = none: F_avg(cond)=0.940858
- dd_mode = spin_echo: F_avg(cond)=0.333382
- dd_mode = xy4: F_avg(cond)=0.501324
- dd_mode = xy16: F_avg(cond)=0.995994
- dd_mode = xy4_per_long_move: F_avg(cond)=0.957410
- paper_trap_depth_std_114pct = site_depth_spread: F_avg(cond)=0.958185

## 与 Level 4 半经典模型对照

- Level 4 contrast_none: 10 轮对比度 0.0568。
- Level 4 contrast_xy4: 10 轮对比度 0.0496。
- Level 4 contrast_spin_echo: 10 轮对比度 0.0351。
- Level 5 none: 10 轮有限脉冲对比度 0.2568。
- Level 5 xy4_per_long_move: 10 轮有限脉冲对比度 0.7582。

## 与论文的关系与边界

- 本级所有数值为模拟量：survival 是经典轨迹标签的经验均值，conditional spin channel 是冻结轨迹上的幺正系综平均；两者都**不是**论文实验 IRB fidelity。
- 论文的 transformed Clifford frame XY4 未启用（配置transformed_clifford_frame=false）；本实现的 xy4_per_long_move 是不同的、更简单的模式。
- nominal 噪声关闭；噪声敏感性参数全部标记 assumed。paper-coherence anchor 仅作对照，不用于校准。
- 47/195 站点是条件独立单原子 + 可配置共同/局域经典噪声；不包含原子间相互作用或并行门串扰。
- 经典轨迹池冻结重放：站点阱深不均匀只通过势能线性缩放进入量子通道，经典 survival 未按站点重新生成。
- 未包含：Rydberg 相互作用、两比特门、测量/SPAM 误差、多能级泄漏与自发 Raman 散射。

图片完成程序化检查；未进行人工视觉审阅。
