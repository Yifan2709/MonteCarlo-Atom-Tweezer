# Codex 执行 Prompt：Level 2B AOD 移动速度边界与变速剖面验证

请在当前工作区中，在已经完成并验证的 Level 0 静态势阱、Level 1 顺序转移和 Level 2 同步波形优化（含平均噪声加热扩展）基础上，实际完成"Level 2B：AOD 移动速度边界与变速剖面验证"。

不要只给方案或代码片段。必须检查现有项目、复用 Level 0/1/2 的公共实现、创建必要代码与配置、运行预检查和全部扫描实验、保存中间数据与最终结果，并诚实报告所有未通过项与未覆盖区间。

如果工作区存在 `sources/`，将其中内容视为只读，不得编辑、移动、重命名或删除；任何输出路径解析后都不得落入 `sources/`。

## 一、Level 2B 的定位

Level 2 的独立验证结论是：5 μK、400 μs 工况下问题处于饱和区，三波形俘获率均为 2000/2000 = 100%，配对差全为 0，无法区分波形优劣。当前 AOD 中心的移动速度量级为：优化波形平均 6.1 mm/s、峰值 11.5 mm/s；压缩顺序基线平均 11.5 mm/s、峰值 17.3 mm/s。距离饱和失效还很远，"多快算太快"完全没有被测量过。

Level 2B 不做新的波形优化，而是把**移动速度本身作为显式扫描变量**，回答三个聚焦的问题：

> 1. **跟随上限**：AOD 阱以多快的速度拖动时，原子开始跟丢或被过度激发（纯移动、深度恒定条件下的速度边界）？
> 2. **端到端上限**：完整移动-降深 drop-off 流程中，SLM 俘获率随平均移动速度如何下降？俘获率 95% 与 50% 对应的临界速度 v_95、v_50 是多少，随温度如何变化？
> 3. **恒速 vs 变速**：在失效边界附近，恒速巡航剖面（梯形速度）与变速剖面（smoothstep 族、Level 2 冻结窗口形状）在相同平均速度与相同峰值速度两种口径下谁更优？

Level 2B 必须完成：

- 原样复现 Level 2 三类主波形作为不可篡改的慢极限锚点；
- 实现以巡航速度为直接参数的梯形恒速剖面族和推广的变速剖面族；
- 执行三组扫描实验（纯移动跟随、端到端速度、恒速 vs 变速配对比较）；
- 用从未参与探索的独立样本复核临界边界点；
- 在扫描覆盖的每个速度点保证步长收敛与功-能平衡；
- 输出临界速度及其置信区间、激发-速度标度律、绝热参数坍缩图和中文报告。

两条诚实底线，违反任何一条即为本 Level 失败：

- **饱和是合法结论但必须显式声明**：若直至配置的最大扫描速度俘获率仍为 100%，必须报告"扫描区间未覆盖失效边界，v_95 > v_max_scanned"，给出已达的最大绝热参数，并说明需要把温度或速度上限推到何处才可能有鉴别力；不得宣称"验证了最大移动速度"。
- **临界速度是估计量不是常数**：v_95/v_50 必须带置信区间（bootstrap 或参数拟合），不得只报点估计；边界点复核必须使用新种子独立样本，不得用扫描样本反复重算后挑好看的结果。

## 二、物理模型与坐标约定

继续使用 Level 0/1/2 已验证的一维焦平面径向高斯势，全部参数与 Level 2 完全一致，不得修改：

- SLM 中心为原点，`U_SLM(x) = -D_SLM·exp(-2x²/w²)`，D_SLM = k_B×140 μK，w_SLM = 1.17 μm（论文 1061 nm 测量值）；
- `U_AOD(x,t) = -D(t)·exp(-2[x−c(t)]²/w²)`，D_AOD(0) = k_B×280 μK，w_AOD = 1.17 μm（与 SLM 后孔径匹配的假设值），c(0) = 2.4 μm；
- 原子 ^133Cs，质量沿用现有项目精确同位素质量；
- 初态采样器沿用 Level 1/2（简谐提议 + E_AOD<0 束缚拒绝），主线温度 5 μK；
- 时变 velocity-Verlet 推进；俘获判据严格 E_f < 0；
- 转移结束后 100 μs 纯 SLM 保持段不变。

除 Level 2 已排除的物理外，本 Level 另外显式排除：**AOD 声学器件带宽、RF 驱动功率与扫描速率等硬件限制**——本 Level 测的是该一维经典模型下的"动力学速度上限"，不是器件规格；也不含变速带来的 AOD 柱面透镜焦距漂移（Level 3 范围）。

### 必须计算并报告的派生量

从现有 `level0_static_trap` 的解析阱频函数计算并写入 `metrics.json`：

- ν_AOD、ν_SLM（280 μK 与 140 μK 高斯阱的解析阱频）及对应振荡周期；
- **绝热速度尺度** v_ad = w_AOD·ω_AOD（w 与 AOD 阱角频率的乘积），所有速度结果除以它做无量纲化；
- **绝热参数** η = v_peak / v_ad，报告每个扫描点的 η 与临界 η_95、η_50；
- **每步位移数** δ = dt·v_peak / w_AOD，作为步长充分性的显式指标。

不得把这些量当成"理论预期正确值"去校准模拟结果；它们只用于无量纲化与量级对照。若实测临界速度与 v_ad 相差 orders of magnitude，应如实报告并讨论，不得反向调参凑数。

## 三、速度与剖面的定义

### 速度口径

任何"移动速度"必须区分三个口径，全部单独报告，禁止混用：

- 平均速度 v_avg = d / T_move（d = 2.4 μm，T_move 为中心实际移动时间跨度）；
- 峰值速度 v_peak = max|ċ(t)|；
- 均方根速度 v_rms。

### 剖面族

所有剖面必须返回 c(t)、D(t) 及一、二阶导数（功-能检查与光滑度诊断需要），复用 `WaveformBase` 接口。

**A. `trapezoid_cruise`（梯形恒速巡航）**——本 Level 的核心新剖面：

- 结构：五次 smootherstep 加速段（时长占比 f_acc）→ 恒速巡航段（速度 v_cruise）→ 五次 smootherstep 减速段（占比 f_acc）；
- 参数：v_cruise（或等价的 T_move）、f_acc ∈ [0.1, 0.45]；
- 巡航段即名义最大速度，必须满足校验 |v_peak − v_cruise|/v_cruise < 1×10⁻⁶（加减速段用 smootherstep 则峰值恰在巡航段）；
- 端点速度为零、加速度为零；总位移严格 d。

**B. `smoothstep_family`（变速剖面族）**：

- Level 2 的 smootherstep 窗口族推广：move/ramp 窗口端点与 ramp 幂次为参数；
- 必须包含 Level 2 冻结波形 `overlapped_optimized_400us` 的窗口形状（move [0.0195, 1.0]、ramp [0.607, 1.0]、power 1.112）作为已验证锚点，只允许整体时间缩放，不允许改窗口形状。

**C. `sequential`（顺序基线）**：Level 1/2 的 52/48 顺序结构，仅整体时间缩放。

**D. `move_only_*`（纯移动剖面）**：深度恒定 280 μK，中心按 A 或 B 的剖面移动，用于实验 1；移动结束后中心保持 0。

### 约束

- c(t) 单调不增、始终位于 [0, d]；端点严格 c(0)=d、c(T)=0；
- 不允许过冲、非物理速度跳变；加减速段端点加速度为零；
- 解析导数与有限差分一致性检查是预检查项；
- 等比例缩放不变性：给定剖面形状，所有导数量随 1/T 缩放正确（单元测试验证）。

## 四、三组实验设计

### 实验 1：纯移动跟随上限（`move_only`）

- AOD 深度全程恒定 280 μK，SLM 全程存在但原子初始束缚于 AOD；
- 中心按梯形剖面移动 2.4 μm，到达后保持 100 μs；
- 扫描 v_cruise：至少 12 个近似对数分布的点，从 6 mm/s（当前慢极限）覆盖到至少 500 mm/s；
- 每个 shot 记录：移动段内 max|x−c|/w、保持段末相对 AOD 阱底激发 ε_AOD/D_AOD、末态位置相对中心偏移；
- **跟随丢失判定**（阈值写入配置，一旦冻结不得事后调整）：移动段或保持段内 |x−c| > 2w，或 ε_AOD/D_AOD > 0.5；
- 输出激发-速度标度律：拟合 ε_med ∝ v^α，报告 α、95% CI 与拟合优度；同时报告跟丢率随速度的上升曲线；不预设理论指数。

### 实验 2：端到端速度上限（完整 drop-off）

- 三套窗口结构并行扫描：(a) 顺序 52/48；(b) Level 2 冻结同步窗口形状；(c) 梯形巡航 + 冻结 ramp 窗口（移动由梯形承担，降深窗口形状与 (b) 相同）；
- 扫描总时长 T（等价扫描 v_avg）：至少 12 个点，从 600 μs 一路压缩，直到出现俘获率 < 90% 的点或到达配置下限（默认 T_min = 20 μs，对应 v_avg = 120 mm/s）；
- 三个温度 5、10、20 μK，各自独立样本（增强鉴别力；5 μK 大概率饱和，10/20 μK 用于把边界拉进可测区间）；
- 每点 ≥ 1000 shots，Wilson 95% CI；
- **临界速度估计**：对每条（结构 × 温度）曲线，用 bootstrap（对扫描点重采样后分段线性插值）或 logistic/probit 拟合估计 v_95 与 v_50 及其 95% CI；两种方法都做，不一致时如实报告差异；
- 失效模式分类（每个失败 shot）：移动段跟丢（max|x−c|/w 超阈）/ 交接失败（跟随良好但 E_f ≥ 0）/ 降深段激发（移动段末激发已高但未跟丢）；三类比例随速度演变单独成图。

### 实验 3：恒速 vs 变速（边界附近配对比较）

- 先从实验 2 的 10 或 20 μK 曲线确定边界区，在其中选 3-5 个速度点（覆盖 v_95 上下约 ±50%）；
- 每个点做两组配对比较，**两种公平口径分别报告，缺一不可**：
  - 口径 1（相同 v_avg）：梯形与 smoothstep 剖面平均速度相同（梯形峰值更低）；
  - 口径 2（相同 v_peak）：峰值速度相同（梯形平均速度更低）；
- 比较对象：梯形恒速 vs Level 2 冻结形状缩放 vs 至少一种其他 smoothstep 窗口；
- common random numbers：同一点内所有剖面对同一初态数组求配对结果；
- 配对 bootstrap 差值 CI + 四格配对计数（都成功/仅A/仅B/都失败）+ 激发与裕量分布对比；
- **边界点复核**：实验 3 选点依据来自实验 2 结果，属探索后选择；因此每个边界点必须用新生成、从未用过的种子独立样本（≥ 1000 shots）复核一次，报告复核值与扫描值是否一致（CI 重叠）；此后不得再改剖面或重跑。

## 五、步长与数值验收（高速工况的关键）

- **自适应步长规则**：dt ≤ min(0.05 μs, T_move/2000, 0.01·w_AOD/v_peak)。即峰值速度下每步中心位移不超过束宽的 1%，且每个剖面至少 2000 步；该规则必须实现为独立函数并有单元测试；
- 在最快与最慢扫描点各做步长减半检查：末态连续量二阶收敛、俘获标签不一致率 ≤ 2% 且翻转样本集中在 E_f≈0 附近；高速端的功-能残差（相对阱深）仍须满足 Level 2 同一阈值（1×10⁻³ 中位数），不得因速度快放宽；
- 保持段能量误差验收同 Level 2；
- 报告每个扫描点实际使用的 dt 与每步位移数 δ，δ > 0.01 的点若因配置上限未能满足规则，必须标记为"数值可疑"并在结论中排除。

## 六、预检查：不过关不得启动扫描

在 Level 2 十项预检查全部复跑通过的基础上，新增：

1. 梯形剖面校验：v_peak == v_cruise（相对误差 < 1×10⁻⁶）、巡航段速度恒定、端点速度与加速度为零、总位移严格 d；
2. 全部新剖面的解析导数与有限差分一致；
3. 等比例缩放不变性：形状不变、T 减半时 v_peak 恰好翻倍、a_max 恰好乘 4；
4. 慢极限回归：把梯形剖面参数化到 v_avg ≈ 6 mm/s（与 Level 2 优化波形同平均速度量级）、相同初态（CRN）时，两者的俘获率配对差 CI 含 0、激发分布无显著差异——证明新剖面族在慢极限退化到已验证结果；
5. 最快扫描点（v_peak 最大处）的 dt 与 dt/2 收敛检查；
6. 静态极限：速度 → 0 时退化为 Level 0 静态阱能量守恒。

结果写入 `preflight.json`；任一关键项失败即停止，不得生成"看似完整"的扫描结果。

## 七、建议配置

在现有配置体系中新增（示例）：

```yaml
level2b:
  name: speed_limit_validation

initial_ensemble:
  main_temperature_uK: 5.0
  scan_temperatures_uK: [5.0, 10.0, 20.0]
  scan_shots_per_point: 1000
  boundary_reshot_shots: 1000
  scan_seed: 22001
  boundary_seed: 22002   # 必须与 scan_seed 不同且从未使用

speed_scan:
  pure_move_cruise_speeds_mm_per_s: [6, 10, 16, 25, 40, 63, 100, 158, 250, 397, 500]  # 近似对数
  end_to_end_durations_us: [600, 400, 300, 220, 160, 120, 90, 65, 45, 32, 25, 20]
  min_duration_us: 20.0
  stop_when_capture_below: 0.90
  trapezoid_accel_fraction: 0.30
  capture_target_levels: [0.95, 0.50]
  critical_speed_method: [bootstrap_interpolation, logistic_fit]

follow_loss:
  max_separation_over_waist: 2.0
  max_excitation_over_depth: 0.5

integration:
  base_dt_us: 0.05
  max_step_displacement_over_waist: 0.01
  min_steps_per_move: 2000
  post_transfer_hold_us: 100.0

comparison:
  points_near_boundary: [0.5, 0.75, 1.0, 1.5, 2.0]   # × v_95
  profiles: [trapezoid_cruise, level2_frozen_shape, smoothstep_variant]
  calibers: [same_average_speed, same_peak_speed]
  bootstrap_samples: 10000

output:
  directory: outputs/level2b_speed_limit/demo
```

`config_used.yaml` 必须记录 CLI 覆盖后的实际值与全部种子。若扫描在到达 T_min 前已出现俘获率 < 90%，可以提前停止并如实记录停止点。

## 八、代码结构与复用要求

先检查现有仓库，一切势、常数、采样器、积分器、统计函数（Wilson、配对 bootstrap）、功-能核算、CSV/图片程序化检查工具都必须从 `level2_joint_transfer`（及其上游 level0/level1 包）复用，不得复制第二套。新增模块建议：

```text
speed_waveforms.py    # 梯形巡航、推广窗口族、纯移动剖面及导数
speed_scan.py         # 实验 1/2 编排：速度/时长网格、温度维度、CRN、失效分类
critical_speed.py     # bootstrap 插值与 logistic 拟合两种临界速度估计 + CI
profile_compare.py    # 实验 3：双口径配对比较 + 边界点独立复核
level2b_cli.py        # --stage preflight / scan / compare / all
level2b_visualization.py
```

所有公开函数用简短中文 docstring；CLI 入口真实可运行（安装后 `python -m ...` 或 console script），不得依赖临时 PYTHONPATH。扫描必须支持断点续跑（每个速度点完成即落盘），单点异常记为失败点并保留，不得让整轮扫描静默丢失。

计算预算：运行前估算总 shot×步数并写入配置与报告（例如 12 点 × 3 温度 × 3 结构 × 1000 shots），预算超限时应降低每点 shots 而不是偷偷放粗步长。

## 九、输出文件

```text
outputs/level2b_speed_limit/demo/
├── config_used.yaml
├── preflight.json
├── speed_scan_results.csv      # 实验1+2全量：每行一个 (实验, 剖面, 温度, 速度/时长, shot统计)
├── critical_speeds.csv         # v_95/v_50 + CI，按 结构×温度×方法
├── scaling_fit.json            # ε∝v^α 拟合与拟合优度
├── comparison_results.csv      # 实验3 双口径配对比较与复核
├── boundary_reshot.csv         # 边界点独立复核
├── failure_classification.csv  # 失效模式分类统计
├── metrics.json
├── image_validation.json
├── summary.md
├── capture_vs_speed.png
├── critical_speed_vs_temperature.png
├── excitation_vs_speed_loglog.png
├── adiabatic_collapse.png
├── trapezoid_vs_smoothstep.png
├── boundary_phase_space.png
├── follow_loss_diagnostics.png
└── timestep_convergence_fast.png
```

`metrics.json` 至少包含：派生量（阱频、v_ad）、预检查结果、全部扫描曲线与 CI、两种方法的临界速度及差异、标度指数 α、失效模式分类比例、双口径比较结论、边界点复核一致性、每个扫描点的 dt 与 δ、是否覆盖失效边界的显式判定（`boundary_covered: true/false`）、饱和时已达的最大 η。JSON 不得含 NaN/Infinity。

## 十、必须生成的诊断图

1. `capture_vs_speed.png`：俘获率 ± Wilson CI 随 v_avg（对数轴），温度分面板，三结构分色；标出 v_95/v_50 及 CI 带；
2. `critical_speed_vs_temperature.png`：v_95、v_50 随温度变化及 CI；
3. `excitation_vs_speed_loglog.png`：实验 1 激发中位数随 v 的双对数图与 α 拟合线；
4. `adiabatic_collapse.png`：俘获率/激发对 η = v_peak/v_ad 的坍缩情况（不同温度是否落到同一曲线）；
5. `trapezoid_vs_smoothstep.png`：两种口径下剖面形状叠加 + 配对结果对比；
6. `boundary_phase_space.png`：边界速度点（v_95 附近）的末态 (x,v) 散点按 captured 着色，叠加 E_SLM=0 分界线；
7. `follow_loss_diagnostics.png`：max|x−c|/w 分布与跟丢标记随速度的演变；
8. `timestep_convergence_fast.png`：最快点的步长收敛与标签翻转率。

绘图规则同 Level 2（非交互后端、≥200 dpi、保存后关闭、程序化检查写 `image_validation.json`、未人工审阅不得声称美观）。

## 十一、测试要求

在保留全部既有测试的基础上新增：

1. 梯形剖面：v_peak==v_cruise 校准、巡航段恒速、端点零速/零加速度、位移严格 d；
2. 缩放不变性：T→T/2 时 v_peak×2、a_max×4、形状归一化后逐点一致；
3. 速度三口径换算自洽（T ↔ v_avg ↔ v_peak ↔ v_rms）；
4. 步长规则函数：给定 v_peak、T_move 返回 dt 上限，含边界情形（超快时由位移约束主导）；
5. 临界速度估计：在构造的已知 logistic 数据上回收 v_95/v_50 误差 < 2%；bootstrap CI 覆盖率抽查；
6. 慢极限回归：新剖面在慢极限与 Level 2 冻结波形配对差 CI 含 0；
7. 失效分类逻辑在构造轨迹上的正确性（跟丢/交接/降深三类）；
8. `boundary_covered` 判定逻辑：全 100% 曲线必须判 false；
9. scan 与 boundary 种子互不相同；断点续跑后结果与一次跑完一致；
10. CLI 各 stage 冒烟测试、CSV/JSON schema、图片程序化检查、输出路径落入 `sources/` 时被拒绝。

测试不得要求"变速一定优于恒速"、不得要求"必须找到失效边界"（饱和是合法物理结论）；应测试的是剖面定义正确、估计量无偏、口径区分清楚、无数据泄漏。

## 十二、执行顺序

1. 检查现有仓库与 Level 0/1/2 代码、配置、测试、输出；运行全部既有测试记录基线；
2. 扩展配置，实现剖面族与导数、缩放不变性与步长规则；
3. 编写 Level 2B 测试并跑通；
4. 运行 preflight（Level 2 十项 + 本级六项），失败则修复重跑；
5. 顺序执行实验 1 → 实验 2（含温度维度与失效分类），每点落盘；
6. 从实验 2 选边界点，执行实验 3 双口径配对比较；
7. 用新种子独立复核全部边界点；
8. 生成全部数据、图与中文报告，程序化检查输出；
9. 若有图像查看工具逐图检查，否则如实记录 `visual_review_performed: false`；
10. 重跑完整测试确保最终输出对应最新代码；检查未修改 `sources/`；
11. 给出简洁、可核查的最终回复。

## 十三、最终回复必须包含

- 实际新增/扩展的模块与实际成功执行的命令；
- 测试通过/失败/跳过数量与 preflight 结果；
- 实验 1：跟丢率与激发随速度的曲线要点、标度指数 α 及 CI；
- 实验 2：每条曲线的 v_95/v_50 及 CI（两种方法）、失效模式组成随速度的演变、是否覆盖失效边界（`boundary_covered`）；若 5 μK 全饱和，明确写出"5 μK 未测得边界，10/20 μK 边界为 …"；
- 实验 3：两种口径下恒速 vs 变速的配对差及 CI，明确结论（谁优/无显著差异）；
- 边界点独立复核与扫描值的一致性；
- 绝热参数：v_ad、η_95、η_50 及其与量级预期的对照讨论；
- 数值可信度：最快点步长收敛、功-能残差、每步位移数 δ 范围；
- 程序化图片检查结果及是否真正执行视觉检查；
- 输出目录与关键文件绝对路径；
- 本 Level 仍未包含、应留给后续 Level 的物理过程（器件带宽、二维/柱面透镜、量子动力学等）。

不得只汇报单条曲线的最好点；不得在看到复核结果后回头改剖面或阈值；不得用扫描样本冒充独立复核。
