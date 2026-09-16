# movement_effects_v0 验证计划与 A 阶段结果

日期：2026-09-15。验证对象：提交 `0b866199d5b6dd2eaa9393c25133f199ea17af30`（分支 `movement_effects_v0`，已核对一致）。

环境差异记录（任务书 §一.6）：任务书工作目录 `D:/work/MonteCarlo-Atom-Tweezer-movement_effects_v0`；本核查在本机镜像仓库执行，HEAD 与目标 SHA 完全一致。`ACCEPTANCE_PROTOCOL.md` 仍不存在于仓库（与合同登记一致）。

本文件分两部分：①A 阶段已执行的小检查与发现；②B/C 阶段验证计划（PLANNED，未执行）。

## ① A 阶段结果（已执行，独立脚本 `check_a.py`，输出 `stage_a_results.json`）

### 通过项

- **A1 Bell 基元（M1）**：独立构造 |Φ+⟩ 概率表与 Pauli 矩阵（不调用被测函数），得 ⟨ZI⟩=⟨IZ⟩=0、⟨ZZ⟩=⟨XX⟩=1、⟨YY⟩=−1，全部精确通过。
- **A4 入口（修复后）**：`paper_integration.cli list` 正常；`simulation_workflow` 注册 `paper-rb2022`/`paper-yb2026` 正常（见下述缺陷修复）。
- **轨迹端点（R1/Y1 共用）**：两类轨迹端点速度≈0（3×10⁻⁶·D/T 数值噪声级）、峰值速度比 1.728/1.875 与代码声明一致、端点位置精确到 D。

### 发现的缺陷（逐条对应任务书"阅读线索"，全部证实）

| # | 位置 | 发现 | 严重度 | 处置 |
|---|---|---|---|---|
| D1 | `quantum_mc.diag_expectation` | 乘积式实现对纠缠态错误：Bell 态上 "ZZ" 返回 0（真值 1），偏差 1.0。生产路径未使用（grep 确认），测试仅用单比特串 "ZI"（乘积式恰好正确）| 潜在误用面 | 建议删除或改为真实联合期望 |
| D2 | `trajectories.sixth_zero_jerk` | 实为五阶剖面 x=6s³−7s⁴+2s⁵，数值端点 jerk(0)=3.4×10⁶≈36·D/T³≠0；名称与模块 docstring 的"六阶/零 jerk"与实现不符（内联注释 a=6/γ≈1.63 如实）| 文档/命名缺陷；Y1 结构判据基于峰值速度差（1.728<1.875）仍成立 | 更名（如 `low_peak_v`）或修正 docstring |
| D3 | `yb_code.teleportation` L375 | `accepted = ~detectedB & ~lostB…any()` 混入**真实丢失**（隐藏标签）→ "postselected" 口径违反"解码只读检测记录"规则 | 口径违规 | 该口径应改标"理想诊断"；观测后选应仅用 detectedB |
| D4 | `yb_vls` | g_optimal=0.04 由论文上界（6×10⁻⁴/单程）反求 → "optimal 在界内"为循环论证；"1.5×bound" 容差非论文口径 | 循环锚定 | Y2 的 optimal 判定应降级为"构造自洽" |
| D5 | `yb_gate.scan_scaling` | 对预置解析公式采样再拟合 → 指数 4/2 是自洽检查，不是门动力学预测（模块 docstring 已声明有效通道，但报告标签需降级）；C_AR=300 为 assumed | 自洽冒充预测风险 | 报告中 Y4 标签改为"给定误差模型的抽样" |
| D6 | `rb_entangled` | L80-82 构造 BatchedState 后未使用，结果全部来自解析读出；解析读出与态演化的等价性未独立证明（死代码+等价性缺口）| 证据缺口 | B 阶段需独立态演化对照 |
| D7 | `simulation_workflow/cli.py` | **0b86619 上统一入口实际损坏**（注册插入双逗号语法错误，import 即失败；原 42 项测试不覆盖该导入）| 高（但正式结果不经此路径）| **已在工作区修复并验证**（7 项 workflow 测试过）；修复未提交——按任务书 §一.5，不得以修复版宣称原提交通过，原提交该入口 = FAIL |
| D8 | `runs.stage_preflight` | B0 读取历史归档 JSON（含 initial_unbound=0），不执行本版本模型 → "本版回归"实为 NOT_RUN | 口径澄清 | B0 应记"历史核对 PASS、本版回归 NOT_RUN" |
| D9 | `integration.py` | `abs(new−old)>1e-9` 传播判据为平凡灵敏度；且当前是标量统计量耦合（加热率/损失率），非同批原子连续演化 | 与"实际接通"声明的差距 | I1 应标注"标量耦合级接通"；连续演化未实现 |

### 对原报告结论的影响

- 13 项 SELF_CHECKED 结论中：M1 基元、R1/R2/Y1 的**结构关系**、Y3 锚点量级、Y5 序关系仍站得住；但 **Y2"optimal 界内"、Y4"指数复现"、Y6"后选口径"、I1"实际接通（连续）"、B0"本版回归"五处标签需要降级或改口径**（见 D3-D5、D8-D9）。
- 原 verify 工具"13 项 0 问题"仅证明文件存在+哈希一致，不构成内容验收（V1 深度不足，D7 即为其漏网实例）。

## ② B/C 阶段计划（PLANNED，未执行）

按任务书 §四通用标准（三档步长、95% 半宽 ≤0.003 且 ≤容差 1/3、冻结拟合区间等）：

- **B（数值收敛/统计稳定）**：R1/Y1 代表点+损失转折点 h、h/2、h/4；R3/Y5 事件级重算置信区间；Y3 逐轮 hazard 与通道关闭消融；Y6 random/adaptive 同预算配对比较。
- **C（论文对照）**：受限于源数据未获取（Zenodo/期刊 PDF 网络受限），R1/R2/Y1 绝对曲线对照 BLOCKED；R3/Y5/Y6 可对标论文标量锚点（含 SPAM 口径核对）。
- 前置失败阻断：D7 修复需先提交独立修复 commit 后，B/C 方可在修复版上执行并单独标注版本。

矩阵见 `validation_matrix.csv`（状态：A_PASS/A_FINDING/PLANNED/BLOCKED/PLANNED_DEFECT_FIXED）。
