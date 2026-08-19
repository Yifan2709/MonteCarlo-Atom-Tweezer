# Level 3 三维 AOD 长距离运输与柱面透镜效应

隔离研究 AOD 单阱长距离搬运：三维经典模型 + 论文 SI 启发的 crossed-AOD 柱面透镜势；留阱率仅指经典末态仍被束缚，不代表量子保真度。

- z_R = 4.076 μm，ω_r/2π = 36.01 kHz，ω_z/2π = 7.31 kHz。
- preflight：15 项检查全部通过（89 s）。

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

## 边界

经典留阱率不可等同于论文的 99.953(2)% 量子通道保真度；lensing severity 为无量纲扫描场景而非论文测量值。未包含：真空损失、光子散射、退相干、pick-up/drop-off 组合、重复往返。

图片完成程序化检查；未进行人工视觉审阅。
