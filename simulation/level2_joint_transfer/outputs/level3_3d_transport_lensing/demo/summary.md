# Level 3 三维 AOD 长距离运输与柱面透镜效应

隔离研究 AOD 单阱长距离搬运：三维经典模型 + 论文 SI 启发的 crossed-AOD 柱面透镜势；留阱率仅指经典末态仍被束缚，不代表量子保真度。

- z_R = 4.076 μm，ω_r/2π = 36.01 kHz，ω_z/2π = 7.31 kHz。
- preflight：16 项检查全部通过（96 s）。

## 加热标度律（§5，简谐、无透镜、小激发极限）

- adiabatic_sine: T 上包络斜率 -6.64（拟合 35 点）。
- constant_jerk: T 上包络斜率 -6.09（拟合 35 点）。
- minimum_jerk: T 上包络斜率 -6.60（拟合 35 点）。
- const_accel_reference: T 上包络斜率 -3.74（拟合 35 点）。
- adiabatic_sine: L 斜率 +2.00。
- constant_jerk: L 斜率 +2.00。
- minimum_jerk: L 斜率 +2.00。
- adiabatic_sine: Δε 对 ω 斜率 -4.03（ΔN = Δε/ħω 再减 1）。
- const_accel_reference: Δε 对 ω 斜率 -1.85（ΔN = Δε/ħω 再减 1）。

- **与 §5B/论文声称的 T⁻⁴ 不一致**：四段 constant-jerk 实测 -6.09，与其 m=2（加速度连续）定义自洽；T⁻⁴/ω³ 属于 m=1（分段恒加速、三角速度）族，已由 const_accel_reference 参照轨迹实测复现。差异如实报告，不改数据。

## 冻结主条件独立验证（validation_pool，一次性）

| 条件 | 留阱 | 留阱率 | Wilson 95% CI |
|---|---:|---:|---|
| paper_anchor_diagonal|sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_diagonal|sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_constant_jerk|sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_constant_jerk|sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control|sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control|sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control_constant_jerk|sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control_constant_jerk|sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |

### 配对比较
- sine_vs_constant_jerk_diagonal_sev1.0: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0，不显著）。
- sine_vs_constant_jerk_straight_sev1.0: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0，不显著）。
- straight_vs_diagonal_sine_sev1.0: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0，不显著）。
- lensing_off_vs_sev1.0_paper_anchor: 差 +0.0000，95% CI [+0.0000, +0.0000]（含 0，不显著）。

### 初态采样
- 物理接受率 1.0000（提议 4000 次，D/kT=56 时几乎无拒绝；简谐提议+束缚拒绝，非完整正则系综）。

## 敏感性（one-factor-at-a-time，每点 500 shots）

- temperature = 3.0: 1.0000（Wilson [0.9924, 1.0000]）
- temperature = 5.0: 1.0000（Wilson [0.9924, 1.0000]）
- temperature = 10.0: 1.0000（Wilson [0.9924, 1.0000]）
- depth = 160.0: 1.0000（Wilson [0.9924, 1.0000]）
- depth = 280.0: 1.0000（Wilson [0.9924, 1.0000]）
- depth = 920.0: 1.0000（Wilson [0.9924, 1.0000]）
- waist = 1.11: 1.0000（Wilson [0.9924, 1.0000]）
- waist = 1.17: 1.0000（Wilson [0.9924, 1.0000]）
- waist = 1.23: 1.0000（Wilson [0.9924, 1.0000]）
- severity = 0.5: 1.0000（Wilson [0.9924, 1.0000]）
- severity = 0.75: 1.0000（Wilson [0.9924, 1.0000]）
- severity = 1.0: 1.0000（Wilson [0.9924, 1.0000]）
- severity = 1.25: 1.0000（Wilson [0.9924, 1.0000]）
- axis_signs_opposite = 1.0: 0.9440（Wilson [0.9203, 0.9610]）

## 探索性 lensing 深度补偿（§6，非论文已实现功能）

- 条件：straight 510 μm / 1.6 ms / sine / sev1.0，500 shots 同初态对照。
- 峰值所需功率倍率 2.36，实际上限 2，饱和时间比例 0.259。
- 未补偿最低有效阱深 118.4 μK → 补偿后 236.8 μK（基准 280 μK）。
- 恒定阱深声明：not_maintained_power_cap_reached。
- lensing_uncompensated: 留阱率 1.0000，retained 中位 Δε 15.833 μK。
- exploratory_lensing_compensation: 留阱率 1.0000，retained 中位 Δε 17.647 μK。

## 边界

经典留阱率不可等同于论文的 99.953(2)% 量子通道保真度；lensing severity 为无量纲扫描场景而非论文测量值。未包含：真空损失、光子散射、退相干、pick-up/drop-off 组合、重复往返。

## Level 4 建议

- pick-up—长距离运输—drop-off 组合序列（含深度 ramp 与 SLM 关断）；
- 重复往返运输与损失累计/自洽存活模型；
- 量子内态退相干（XY4、强度噪声 dephasing）与成像误差；
- 论文式 RF 功率-位置校准的 lensing 补偿（本级的补偿仅为探索模式）。

图片完成程序化检查，并经逐图视觉审阅（12 张，ZCode (GLM vision) 逐图视觉审阅 2026-08-19）：全部通过。
