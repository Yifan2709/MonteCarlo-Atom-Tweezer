# Level 2: 同步移动-降深波形优化与鲁棒蒙特卡洛验证

本项目从已验证的 Level 1 工程复制而来：`src/level0_static_trap` 与
`src/level1_transfer_1d` 保持原样（基线实现），新增 `src/level2_joint_transfer`
完成 Level 2 任务——在相同一维经典模型与初态分布下，把 600 μs 顺序 drop-off
压缩到 400 μs 并允许 AOD 位置与深度重叠变化，同时保持或提高 SLM 俘获率、
降低末态激发，并用从未参与优化的独立样本做最终验证。

## 运行

在本目录使用项目自带的解释器：

```bash
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest

.venv/bin/python -m level2_joint_transfer.cli \
  --config configs/level2_joint_transfer.yaml \
  --stage all \
  --output-dir outputs/level2_joint_transfer/demo
```

`--stage` 支持 `preflight` / `optimize` / `validate` / `robustness` / `all`；
`validate` 与 `robustness` 依赖 `optimize` 冻结的 `best_waveform.yaml`。

## 防过拟合设计

- `optimization_pool`（128 shots）只用于 LHS 候选探索（dt=0.10 μs）；
- `selection_pool`（512 shots）只用于入围精修与样条精修（dt=0.05 μs）；
- `validation_pool`（2000 shots）在波形冻结后一次性评估三个主波形，
  使用与基线完全相同的初态数组（common random numbers）；
- 波形比较使用配对 bootstrap 与四种配对计数，不凭点估计声称显著；
- 三池种子互不相同且写入配置与 `metrics.json`；`optimizer.py` 源码层
  不引用 validation 池（有测试守护）。

## 物理与数值要点

- 三类波形：`sequential_600us`（严格复现 Level 1）、`sequential_400us`
  （同形状压缩）、`overlapped_optimized_400us`（smootherstep 窗口族 +
  可选 8 控制点单调 PCHIP 样条精修）；
- 显式时变势按 dE/dt = ∂U_AOD/∂t 做功-能核算，分别记录深度功、移动功、
  总功与残差 R_W(t)；转移后 100 μs 纯 SLM 保持段检查静态能量误差；
- 俘获判据严格 E_f < 0（E_f = 0 不俘获）；|E_f|/D_SLM < 1e-3 的样本
  单独标记 near-threshold 但仍按严格符号计数；
- 预检查十项（含运行 Level 0/1 测试、Level 1 复现、步长收敛、功-能平衡、
  保持段有界性）全部通过才允许启动优化。

## 平均强度确定性噪声加热（`noise_mean_heating`）

三类光镊噪声以**系综平均加热功率**确定性注入（无随机性、默认常数，
随时间累积使系统越发不稳定），三个噪声均以函数形式外放为输入参数
（`noise_heating.py`，签名 `(x_m, v_m_s, t_s, e_osc_J) -> W`，可注入任意
时变函数；`build_mean_heating` 构造的常数模型可 pickle、支持多进程）：

| 通道 | 默认平均模型 | 单次 400 μs+100 μs 转移的注入 |
|---|---|---|
| 光子反冲 | P = 2·E_r·Γ_sc ≈ 0.36 μK/s（1061 nm、280 μK，Cs D2 两能级） | ~2×10⁻⁴ μK（可忽略） |
| 强度噪声（参量） | dE/dt = Γ_par·ε_osc，Γ_par = (π²/4)ν₀²S_RIN(2ν₀) ≈ 16 s⁻¹（−80 dBc/Hz） | ε_osc 增长 ~1.3% |
| 指向噪声 | P = m·ω₀⁴S_x(ν₀)/8 ≈ 1.3 μK/s（S_x = 1 pm²/Hz @25.5 kHz） | ~7×10⁻⁴ μK（可忽略） |

- 加热踢沿当前速度方向缩放 v ← sign(v)·√(v²+2ΔE/m)，逐通道累计注入能量；
  功-能恒等式扩展为 ΔE = W_ext + W_noise + R_W，开噪后 R_W 仍处积分截断水平；
- 配置 `noise_mean_heating.enabled: true` 后 validate/robustness 阶段生效
  （优化阶段保持无噪声，不改变已冻结流程）；`validation_shots.csv` 新增
  4 列注入能量（recoil/parametric/pointing/total）；
- 幅度参数全部为 assumed_sensitivity_only；演示见 `demo_mean_heating.py`
  （输出 `outputs/level2_joint_transfer/demo_mean_heating/`）。

**已知既有问题**：`validation_shots.csv` 等处的 `*_uK` 字段数值实际单位为
开尔文（换算 `/1.380649e-29/1e6` = J→K），沿自 Level 2 既有实现；新增噪声
列保持同约定以便跨列运算，跨列比较时无需再换算。

## 输出

`outputs/level2_joint_transfer/demo/` 下包含 `preflight.json`、`candidates.csv`
（全部候选含无效与失败）、`best_waveform.yaml`、`waveforms.csv`、
`validation_shots.csv`、`robustness.csv`、`metrics.json`、`image_validation.json`、
`summary.md` 与 10 张诊断 PNG。

## 边界

仍是一维经典模型：不含二维/轴向运动、AOD 柱面透镜效应、光子散射、技术噪声、
量子内态与相干性、加热/损失随机过程、多原子相互作用、反向 pick-up 与
610 μm 长距离运输。经典俘获率不是实验存活率或转移保真度，不可直接对标
文献 99.81%。

## Level 4：端到端往返协议与半经典相干性

`level4_roundtrip_coherence` 包在同一个三维势实现（SLM 三维高斯 + Level 3
crossed-AOD 透镜势）上编排完整 pick-up—split—375 μm 对角长运输—54 μs 远端
保持—返回—merge—drop-off 往返协议：

```bash
level4-roundtrip --config configs/level4_roundtrip.yaml --stage all \
  --output-dir outputs/level4_end_to_end_roundtrip/demo
```

- **Protocol A**（论文 Extended Data Fig. 10e 时序锚点）：100+850+400+1450
  +54+1450+400+850+100 μs，SLM 180→60→180 μK、AOD 0→280→0 μK，
  smootherstep 深度 ramp（模型假设）；**Protocol B**：冻结 Level 2
  `best_waveform.yaml` 的时间反向 pickup + 正向 dropoff（自带 2.4 μm 短
  移动，不重复 split/merge）+ Level 3 冻结 v_s 长运输，SLM 恒 140 μK；
- 上游产物（Level 2 波形/指标、Level 3 指标/配置）以 SHA-256 冻结进
  `upstream_manifest.json`，参数来源逐项标注（论文锚点/冻结/假设）；
- checkpoint basin 归属（`bound_to_slm/bound_to_aod/shared_or_ambiguous/
  unbound`，连线最低点+势垒+相空间可达性）、first-failure 状态机、
  absorbing 与非吸收生存、recaptured 区分；
- 分段功-能账本（SLM 深度功/AOD 深度功/移动功/lensing 偏移功四通道，
  ΔE=W_ext+R_W 逐段校验）；
- 半经典差分光移相位 φ=∫y(t)·(ηU)/ħ dt 以 A_SLM/A_AOD 分离累积（η 线性
  可重算），理想瞬时 none/spin-echo/XY4 toggling（脉冲自动避开网格与分段
  边界），条件相干对比度与 usable coherent fraction；
- preflight 20 项（含 Level 0-3 全部测试、Level 2/3 三维回归、控制时间
  反演到机器精度、0.10/0.05/0.025 μs 收敛）通过后才开放 2000-shot
  validation；500-shot×10-round 重复、128-shot×40-round 长尾、8 项消融、
  温度/深度/对准/lensing/η/噪声/脉冲时序敏感性。

经典成功概率与 C_cond 不可等同论文 IRB fidelity；未包含有限脉冲、自旋
哈密顿量、Clifford/IRB、Raman 散射与多原子相互作用（Level 5 范围）。


## Level 4：端到端往返协议与半经典相干性

`level4_roundtrip_coherence` 包在同一个三维势实现（SLM 三维高斯 + Level 3
crossed-AOD 透镜势）上编排完整 pick-up—split—375 μm 对角长运输—54 μs 远端
保持—返回—merge—drop-off 往返协议：

```bash
level4-roundtrip --config configs/level4_roundtrip.yaml --stage all \
  --output-dir outputs/level4_end_to_end_roundtrip/demo
```

- **Protocol A**（论文 Extended Data Fig. 10e 时序锚点）：100+850+400+1450
  +54+1450+400+850+100 μs，SLM 180→60→180 μK、AOD 0→280→0 μK，
  smootherstep 深度 ramp（模型假设）；**Protocol B**：冻结 Level 2
  `best_waveform.yaml` 的时间反向 pickup + 正向 dropoff（自带 2.4 μm 短
  移动，不重复 split/merge）+ Level 3 冻结 v_s 长运输，SLM 恒 140 μK；
- 上游产物（Level 2 波形/指标、Level 3 指标/配置）以 SHA-256 冻结进
  `upstream_manifest.json`，参数来源逐项标注（论文锚点/冻结/假设）；
- checkpoint basin 归属（`bound_to_slm/bound_to_aod/shared_or_ambiguous/
  unbound`，连线最低点+势垒+相空间可达性）、first-failure 状态机、
  absorbing 与非吸收生存、recaptured 区分；
- 分段功-能账本（SLM 深度功/AOD 深度功/移动功/lensing 偏移功四通道，
  ΔE=W_ext+R_W 逐段校验）；
- 半经典差分光移相位 φ=∫y(t)·(ηU)/ħ dt 以 A_SLM/A_AOD 分离累积（η 线性
  可重算），理想瞬时 none/spin-echo/XY4 toggling（脉冲自动避开网格与分段
  边界），条件相干对比度与 usable coherent fraction；
- preflight 20 项（含 Level 0-3 全部测试、Level 2/3 三维回归、控制时间
  反演到机器精度、0.10/0.05/0.025 μs 收敛）通过后才开放 2000-shot
  validation；500-shot×10-round 重复、128-shot×40-round 长尾、8 项消融、
  温度/深度/对准/lensing/η/噪声/脉冲时序敏感性。

## Level 5：有限脉冲单比特动力学、噪声谱与多站点 RB/IRB

`level5_finite_pulse_irb` 包把内部态升级为旋转框架下有限时长微波脉冲驱动
的 ^133Cs 钟态量子比特（`H=ħ/2[Δσ_z+Ω(cosφσ_x+sinφσ_y)]`，RWA），在
Level 4 冻结轨迹上执行单通道 PTM 表征、reference RB 与
transport-interleaved RB（论文 Fig. 5c Methods 设计）：

```bash
level5-finite-pulse-irb --config configs/level5_finite_pulse_irb.yaml \
  --stage all --output-dir outputs/level5_finite_pulse_irb/demo
```

- Rabi 锚点 Ω₀=2π×24.611 kHz（论文 Fig. 4a）；η_SLM=η_AOD=1.3e-4；
  nominal 噪声关闭（避免把假设参数伪装成实验输入）；
- 24 元素单比特 Clifford 群（规范型、唯一性/封闭/逆/Pauli 共轭自检），
  Z-X-Z-X-Z 分解，bare 与 SCROFULOUS 两种脉冲实现；SCROFULOUS 按论文
  Methods 公式逐段构造（sinc=sin(x)/x，[π/2,π] 根分支），平均面积
  2.28π/Clifford（论文平均锚点 2.02π）；
- 有限脉冲 DD（none/spin_echo/XY4/XY8/XY16/transport-aware XY4），
  transport XY4 在 outbound/inbound 长移动段各放一个对称块；transformed
  Clifford frame 模式默认关闭；
- 6 条冻结轨迹池（protocol_a/static_idle 各 2000 条，protocol_b/
  transfer/transport/no-lensing 各 1000 条）以 NPZ 缓存，survival 标签
  来自 Level 4 引擎（不修改）；frozen replay 与 joint（联合经典-自旋）
  两种模式，joint 一致性在 preflight 中验证；
- 17 个通道条件的 survival-conditioned PTM（六输入态、bootstrapped CI）；
  reference RB（n=1..1000，72 序列/点 × 47 站点，SCROFULOUS 门）与
  IRB（N=80，M∈{0..80}，72 strings × 47 站点，transport XY4），
  survival/conditional/joint 三口径分开拟合；多模型（指数 vs stretched）
  比较、非指数诊断、IRB 适用性判定；
- 47/195 站点 ensemble（条件独立 + 可配置共同/局域噪声）；paper-coherence
  anchor（T2*=14.0/25.5 ms、T2=12.6 s）仅作对照不校准；
- preflight 26 项（含 Level 0-4 全部测试重跑、干净进程导入、解析
  Rabi/Ramsey/echo、PSD/空间相关、Clifford、理想 RB、PTM 参数回收、
  dt 收敛）全部通过才开放正式 RB；测试 159 项全部通过（Level 0-4 的
  124 + Level 5 的 35）。
