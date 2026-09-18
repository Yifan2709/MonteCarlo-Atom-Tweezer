# 对比报告：`codex/long-distance-aod-transport` 新提交 vs `experiment_toolkit`

日期：2026-09-17。对比对象：

- **A 方（新拉取）**：远端新分支 `codex/long-distance-aod-transport`，6 个提交（`63dc14e`→`60aa6db`），基于 `movement_effects_v1`/`1fed632`，2775 个文件、约 88.6 万行新增（含大型 runs 归档）。已检出到本地同名分支。
- **B 方（本轮实施）**：`experiment_toolkit/`（包 `tweezer_experiment` v1.0.0），目前在 `movement_effects_v1` 工作区的**未跟踪目录**，未提交；A 方分支不含它。

一句话结论：**两者是同一物理血统下的两种正交交付——A 是回答一个物理问题的研究工程（长距离运输的条件性存活速度边界 + 论文物理论证修正），B 是按任务书 F1–F5 产品化的独立离线规划工具。不存在替代关系；物理口径一致、能力互补，各有对方没有的东西。**

## 1. 目标与形态

| 维度 | A：codex 长距离研究 | B：experiment_toolkit |
|---|---|---|
| 回答的问题 | 270/510/610 µm 长距离 AOD 运输在 99%/99.9% 机械存活目标下的**条件性速度/时间边界**；两篇论文运输物理的正确形态 | "输入装置与操作条件 → 可追溯预测、候选方案比较、波形导出回读"的**通用离线工程链**（任务书 F1–F5） |
| 形态 | 研究仓库内研究报告模式：`reports/long_distance_aod_transport_20260917/`（PROTOCOL/AUDIT/PARAMETERS/PHYSICS/BOUNDARIES/REPRODUCE + 冻结队列 + 20 余个 runs 归档）+ 新生产模块 `src/long_distance_transport/model.py`（numba） | 独立 src 布局包：wheel 可装（A01 干净环境验证）、CLI `tweezer-experiment` + Python API、7 个自包含示例、48 项测试、A01–A20 验收 |
| 依赖 | 深嵌研究仓库（import `continuous_transfer.engine.field`、`noisy_survival.beam_noise`、`paper_integration`） | 仅 numpy+PyYAML；零研究仓库 import；论文资源仅经 `adapters.paper` 显式可选 |
| 交付产物 | 边界表（BOUNDARIES.md）、8 份冻结交接控制（handoff_controls）、机制/敏感性/收敛证据、论文页渲染与 Zenodo 源数据核对 | 五项能力软件链 + 验收报告/机读清单 + 证据目录 |

## 2. 物理覆盖对比

| 物理项 | A | B |
|---|---|---|
| 三维有限深高斯场 + 双轴柱面焦移（透镜） | ✓（复用 `engine.field`） | ✓（同源 `physics/field.py`，波长/束腰全参数化） |
| 随机力扩散（RIN 经力矢量、指向经 dF/dx、反冲∝局部深度） | ✓（复用 `beam_noise`，固定 PSD 情景） | ✓（同公式逐项迁移；A16 对归档 Joint 模型四曲线 3.5σ 内复现） |
| 位点差异 + 实际位点束缚初态 | ✓（独立实现，显式波长） | ✓（`sample_site_bound`，A16 验证） |
| 轨迹族 | **论文精确族**：Rb cubic 3s²−2s³、Yb zero-jerk γ=1.5625、Yb min-jerk γ=1.875（含精确系数与 jerk/加速度端点分析） | 6 族：paper_manual（=Rb cubic 同式）、minimum_jerk、constant_jerk（四段 S）、adiabatic_sine、linear；**缺 Yb γ 族** |
| 声响应 τ（焦移滞后） | ✓ τ∈{0,2,5} µs 进边界表 | **✗ 无该参数**（透镜按指令速度即时耦合） |
| 操作协议 | A 纯运输 / B 静→动→静（200 µs 共址线性功率交叉渐变 + 100 µs 静止）/ C 重复 1/10/50 单程 | 任意段/事件编排：static/transfer(双向顺序)/csv/嵌套 repeat/任意观察点；B 型交叉渐变**可经 CSV 表达**但无内置段类型 |
| 观测/判据 | 仅完整单程末静止时 E<0；Clopper–Pearson 双侧 95%；条件损失只对存活队列定义 | 段末单阱 E<0 检查点；吸收/再俘获两语义；Wilson 95%（**无 CP 区间**） |
| 短途拾取/放回（2.4 µm 论文协议） | 明确声明不属于本研究 | ✓（含 ML 回放、时间反演放回，A16 验证） |
| 受约束搜索 / 波形导出 / 回读复算 / 缓存续跑 | ✗ | ✓（F4/F5/A11–A15） |

## 3. 验证体系对比

| 维度 | A | B |
|---|---|---|
| 判据预注册 | PROTOCOL.md 扫描前写入目标与区间口径 | A01–A20 判据先固化再执行；验收模板逐项回填 |
| 统计口径 | Clopper–Pearson 精确区间；跨目标"未确定"带显式标注；合并种子后由总计数重算分位数 | Wilson 区间；逐原子计数可复算（A06）；并列候选"统计不可区分"显式报告 |
| 独立物理参照 | DOP853 独立积分、光学势差分、作者导数数据、逐点收敛记录（benchmarks.json） | 有限差分力核对、解析阱频对比振荡频谱、Verlet 二阶收敛、扩散账本实现=期望、拾取/放回位级时间镜像、A16 分布级回归 |
| 复现机制 | REPRODUCE.md + provenance.json + 源码快照哈希 + 冻结队列 | run 记录（哈希/种子/版本）+ RunStore 缓存失效 + 搜索断点续跑 + wheel 哈希 |
| 诚实边界声明 | 无绝对装置上限、无量子保真度、公开输入不足以定量复现论文 | 分层结论：实验校准/真实设备 NOT_RUN，硬件适配 PARTIAL（合成标定） |

## 4. 各自发现并修复的真实缺陷（无重叠）

- **A 修复**（均在 `paper_integration`/旧链）：Yb 五阶轨迹族整体错误（旧 6,−7,2 系数不代表 zero-jerk，γ=1.875 由错误族反求）；恒加速度段位置不连续；拒绝采样布尔索引长度错误；dt 用传参而非段网格；末能量漏时变 SLM/透镜项；运动判据误用实验室系动能；Yb ED2 交接被双重推进；`hash()` 进程随机种子；Rb Eq.(2) 与 110 µm 误标；汇总层"上轮全灭时条件损失=100%"与"复制代表记录"两个统计错误。
- **B 修复**（验收过程暴露，均已留档）：导出物理表速度列单位换算错误（×1e-6 应为 ×1，由控制保真度指标抓出）；SearchJournal 头部行解析崩溃；RunStore 种子口径；ML 回放放回段误用解析波形（应为 ML 自身时间反演，A16 首跑 ΔS=0.083→修复后 ≤0.018）。

**交叉一致性佐证**：A 的 Rb cubic 修正（3s²−2s³、端点加速度跳跃、与旧四段 S 曲线的区分）与 B 独立实现的 `paper_manual` 族及其"加速度跳跃披露不改"的语义**完全一致**——两条线对论文轨迹的读法独立收敛。B 的 A16 回归通过也意味着 B 与 A 所复用的同一扩散内核在分布意义上等价。

## 5. 差距与整合建议

B 相对 A 的缺口（可作为 toolkit 后续项）：声响应 τ 参数；Yb γ 族精确轨迹（B 的族框架可直接挂载，系数由 A 的 `exact_trajectories.polynomial` 提供）；Clopper–Pearson 区间选项；B 型共址交叉渐变内置段类型；边界表式的"目标存活→最短时间"扫描产出（B 的 F4 目标类型已支持 minimize_duration + capture 下界，语义同源，补一个边界扫描驱动即可）。

A 相对 B 的缺口（研究线不需要、产品线需要）：独立可安装性、通用输入契约、导出/回读、搜索、缓存与续跑、干净环境验收。

整合路径建议：`experiment_toolkit` 提交入库（独立目录或独立分支）；从 A 分支 cherry-pick `exact_trajectories.py` 与 τ 耦合实现进 toolkit 物理层（保持参数化），并在 A16 之外补一组对 A 边界表抽样格点的回归对照。

## 6. Git 状态说明

- 本地当前在 `codex/long-distance-aod-transport`（= origin，`60aa6db`）；`movement_effects_v1` 与 origin 一致（`1fed632`）。
- `experiment_toolkit/` 为未跟踪目录，切分支不受影响、未混入 A 的任何内容；A 分支也不包含它。
- 两线合并顺序不影响各自结论；合并时建议保留双方的证据目录与报告原样。
