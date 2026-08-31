# Level 2B: AOD 移动速度边界与变速剖面验证（tag: level2, speed_modify）

本目录是 Level 2 的**独立扩展工程**：`src/level0_static_trap`、
`src/level1_transfer_1d`、`src/level2_joint_transfer` 为原 Level 2 工程的
逐字节副本（校验和一致，作为不可篡改基线），新增 `src/level2b_speed_limit`
实现 [`level2b_speed_prompt.md`](../../level2b_speed_prompt.md) 定义的任务：

> 把 AOD 移动速度本身作为显式扫描变量——实验 1 纯移动跟随上限、
> 实验 2 端到端速度上限（v95/v50 + 双方法 CI）、实验 3 恒速巡航（梯形）
> vs 变速（smootherstep）双口径配对比较，全部带独立复核与诚实边界声明。

原 Level 2 工程目录（`simulation/level2_joint_transfer/`）未被修改。

## 运行

```bash
cd simulation/level2b_speed_modify
.venv/bin/python -m pytest                      # 全部测试（含 level0/1/2 基线）
.venv/bin/python -m level2b_speed_limit.level2b_cli \
    --config configs/level2b_speed_limit.yaml --stage all \
    --output-dir outputs/level2b_speed_limit/demo
```

`--stage` 支持 `preflight` / `scan` / `compare` / `noise` / `maxspeed` / `all`；
`scan` 与 `compare` 依赖 `preflight.json`；`noise`（实验 4）与 `maxspeed`
（实验 5）依赖各自配置段，可独立运行。

实验 5（`--stage maxspeed`）最高速度测算：初始噪声模型（×1 默认幅度）下
按硬件条件（AOD 深度倍数 / 束腰）测端到端 v95，并在存活率陡降区对
"AOD 起动到停止时长相同"的四种移动方案（smootherstep / 梯形 0.30 /
快加速梯形 0.15 / 三角 0.50）做同位移配对比较；产物为
`max_speed_hardware.csv`、`max_speed_profiles.csv`、`max_speed_summary.json`。
调试旗标 `--skip-level2-preflight`、`--skip-existing-tests` 仅用于冒烟测试，
正式运行不得使用。

## 关键设计

- **梯形恒速巡航剖面**（`SpeedWaveform(move_profile="trapezoid")`）：
  smootherstep 加/减速 + 恒速巡航，v_cruise = d/[(1−f)·T_move]，
  峰值速度严格等于巡航速度（校准 <1e-6）；
- **自适应步长**：dt ≤ min(50 ns, T_move/2000, 0.01·w/v_peak)，
  转移段与保持段两段积分均为同一种 velocity-Verlet；
- **失效三分类**：AOD 深度 ≥50%D0 承载段内 |x−c|>2w → 跟丢；
  深度 <50%D0 降深段内逃出 → 降深段激发；全程未逃出但 E_f≥0 → 交接失败；
- **防作弊**：scan_seed(22001) 与 boundary_seed(22002) 互不相同；
  曲线未穿越 95% 水平时 `boundary_covered=false`，只允许声明
  "v95 > 最大扫描速度"；边界点用 boundary 种子独立复核后才算数；
- **复用不复制**：势/常数/积分器/采样器/Wilson 与配对 bootstrap/功-能核算
  全部 import 自 level0/1/2 包；`smootherstep` 模式与
  `OverlappedWaveform` 逐点一致（有测试守护）。

## 实验 4：平均强度噪声加热 × 移动速度（`noise_speed_validation`）

复用 `level2_joint_transfer.noise_heating` 的三通道确定性平均加热模型
（`build_mean_heating`/`scaled_heating`，可 pickle 支持多进程），通过
`--stage noise` 或 `all` 执行，逐行增量写入 `noise_speed_results.csv`，
v95 双方法估计写 `noise_v95.json`，四联诊断图 `noise_speed_validation.png`。
scan 阶段保持无噪声（冻结基线），噪声只做敏感性评估、不回调任何参数。

默认幅度（全部 assumed_sensitivity_only，1061 nm、ν₀ 取 SLM 解析阱频
25.46 kHz）：光子反冲 2·E_r·Γ_sc = 0.355 μK/s（Γ_sc=2.77 s⁻¹）、
参量 Γ_par = (π²/4)ν₀²·S_RIN = 16 s⁻¹（−80 dBc/Hz）、指向
m·ω₀⁴·S_x/8 = 1.31 μK/s（1 pm²/Hz）。

Demo 运行结果（300 shots/点）：

| 子实验 | 结论 |
|---|---|
| 4a 边界移动 | v95(bootstrap)：关 40.8 [40.2, 41.7] ≈ ×1 40.8 [40.2, 41.7]（完全不变）→ ×100 20.8 [6.1, 21.4] mm/s（腰斩）。×1 档与关档逐点俘获率完全一致。 |
| 4b 开/关配对 | 8 个 (速度×温度) 组合配对差全部 = 0（CI [0,0]），裕量退化 ≤4.5×10⁻⁴，单次注入 ~0.04-0.09 μK —— 默认幅度在单次转移内不可分辨。 |
| 4c 幅度敏感性 | v_peak=55 mm/s：×1/×10 俘获率 100%（裕量 0.84）；×100 → 0.50（裕量 0.42）；×1000 → 0.50（裕量 −0.29）。退化在 ×10 与 ×100 之间起步。 |
| 4d 通道分解 | ×100 下：仅参量 0.480 ≈ 三通道 0.477；仅反冲 1.000（注入 2.3×10⁻³ μK）；仅指向 1.000（8.3×10⁻³ μK）。参量通道注入 11.6 μK，是唯一显著通道。 |

- 功-能恒等式扩展 ΔE = W_ext + W_noise + R_W 在全部 61 个噪声点成立：
  残差中位数 8.9×10⁻⁷ ~ 6.6×10⁻⁶（阈值 1×10⁻³）；
- 注入能量列沿用 level2 的 `*_uK` 命名约定（数值单位实为开尔文，×10⁶ 换算 μK）；
- 方法差异如实报告：×100 档 bootstrap 插值 20.8 与 logistic 拟合 15.6 mm/s
  相差明显（俘获率悬崖过陡时 logistic 尾部偏缓），以非参 bootstrap 为主口径。

## 输出

`outputs/level2b_speed_limit/demo/`：`preflight.json`、`speed_scan_results.csv`
（逐点增量写入，断点续跑）、`scan_flags.npz`（逐 shot 俘获旗标，供 bootstrap）、
`critical_speeds.csv`、`scaling_fit.json`、`comparison_results.csv`、
`noise_speed_results.csv`（实验 4 全部子实验）、`noise_v95.json`、
`metrics.json`、`image_validation.json`、`summary.md` 与 9 张诊断 PNG。

## Demo 运行结果（2026-08-23，300 shots/点）

- 预检查 8/8 通过（含 Level 2 十项复跑，456 s）；全部测试 100 项通过
  （98 快 + 2 slow CLI）；基线三包与原工程校验和一致（`8d9ae671cfe5b53f`）。
- 派生尺度：ν_AOD = 36.0 kHz、ν_SLM = 25.5 kHz、绝热速度 v_ad = 264.7 mm/s。
- **实验 1（纯移动）**：v_peak ≤ 63 mm/s（η≤0.24）零跟丢；100 mm/s 跟丢率 75%，
  ≥158 mm/s 全跟丢。激发标度 ε_med ∝ v^1.33（CI [1.08, 1.82]）。
- **实验 2（端到端）**：9 条曲线全部覆盖失效边界，v95(v_avg) = 41–56 mm/s，
  随温度升高而下降（20 μK 时 frozen 41.2 < sequential 53.2 ≈ trapezoid 54.5 mm/s，
  排序与各结构峰值/平均速度比一致：smootherstep 1.875 > sequential 1.5 >
  trapezoid 1.43）。早停在 <90% 处截断部分曲线的 v50（above_scanned，如实记录）。
- **实验 3（双口径）**：同 v_avg 时低峰值者胜（f=0.15 梯形 > f=0.30 梯形 >
  smootherstep，配对差显著）；同 v_peak 时低平均速度（更长时长）者胜
  （smootherstep +0.0875 显著优于 f=0.30 梯形，f=0.15 梯形 −0.18 显著更差）。
  独立复核 30/30 组合 CI 重叠。
- 图片程序化检查 8/8 通过；未做人工视觉审阅（`visual_review_performed: false`）。

## 边界

一维经典模型；不含 AOD 声学器件带宽/RF 扫描速率限制（测的是动力学上限
而非器件规格）、二维/柱面透镜、光子散射（随机部分）、量子内态与多原子。
