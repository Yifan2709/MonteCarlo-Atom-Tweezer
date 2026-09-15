# AOD–SLM Monte Carlo 综合阶段报告 v9

更新日期：2026-09-14。基于代码提交 `5e55bb0` 加本轮工作区修改。
前一版为 18 页 Level 2C 汇报，本版加入 Level 3、4、5，PPT 共 32 页。
本报告读取已有 demo 数值并核查源代码，不声称重新完成大规模物理验收。
所有引用文件及 SHA-256 保存在同名 `.evidence.json`。

配套 [32 页 PPT](AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9_reviewed.pptx)。

## 1. 层级关系与计算对象

| 层级 | 输入及方法 | 输出 | 操作单位 |
|---|---|---|---|
| Level 2C | 一维 SLM+AOD；短转移波形、初态散布与平均加热 | 逐单程存活 S(n)、损失阶段、波形比较 | 2.4 μm 单程，60 次为 30 短往返 |
| Level 3 | 已装入 AOD 的三维热初态；移动高斯阱及柱面透镜 | 长运输留阱、激发能、轴向分岔与敏感性 | 单次 270/510/610 μm |
| Level 4 | Level 2 冻结短转移与 Level 3 长运输模型 | 完整往返成功、连续多轮加热、相位与条件对比度 | 每轮含 375 μm 去程和 375 μm 回程 |
| Level 5 | Level 4 冻结轨迹；有限微波脉冲二能级传播 | PTM、reference RB、固定 N 的 M 扫描、站点统计 | 默认 M 为完整组合往返数 |

Level 2C 最新波形没有自动接入 Level 4/5。Level 4 Protocol B 仍使用旧 Level 2 冻结优化波形。
Level 5 每次 move 独立抽取单轮轨迹，不具有 Level 4 连续重复轨迹的机械加热记忆。

## 2. Level 2C：保留原汇报，并修正解读边界

PPT 第 3–17 页保留原来的论文数字化、手工波形、初态、存活曲线、敏感性、ML 回放与验收。
默认位移 2.4 μm，60 次单程不含任何 375/510/610 μm 长距离运输。
此前 T40 扫描同时改变了代码默认参量加热参考能量，不能完全解释为纯初温效应。
当前参量/指向噪声 PSD 系数、手工 constant-jerk 波形定义和跨丢失后的随机配对仍需单独核查。
本轮不改变这些参数，避免把报告整合变成未标定的物理拟合。

## 3. Level 3：三维长运输与柱面透镜

默认 Cs-133、1055 nm、束腰 1.17 μm、AOD 深度 280 μK、初温 5 μK。
Rayleigh 长度约 4.076 μm，径向频率约 36.01 kHz，轴向约 7.31 kHz。
初态为三维谐振提议加束缚拒绝，传播使用三维 velocity-Verlet，结束后保持 100 μs。
两扫描轴速度决定柱面焦点偏移，模型可出现有效深度下降与轴向双势阱。
对角线总路程 L 对应每轴 L/√2，不能将 610 μm 当成每轴各走 610 μm。

配置支持 straight/diagonal、270/510/610 μm、400–2000 μs、sine/四段 constant-jerk/minimum-jerk。
已有粗扫描 CSV 为 96 个格点的子集，应以 CSV 中实际参数为准。
正式验证每条件 2000 shots；数值步长 0.05 μs，另有 0.025 μs 收敛检查。

| 冻结条件 | 留阱数 | 留阱率 | Wilson 95% 区间 |
|---|---:|---:|---|
| paper_anchor_diagonal / sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_diagonal / sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_constant_jerk / sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| paper_anchor_constant_jerk / sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control / sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control / sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control_constant_jerk / sev0.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| straight_control_constant_jerk / sev1.0 | 2000/2000 | 1.0000 | [0.9981, 1.0000] |

这些结果说明所选主条件在当前经典模型中没有观察到失阱，不等于证明实验零损失。
默认透镜强度是相对参考场景，severity=1 对应 v_s≈0.53917 m/s，是推导场景值而非硬件测量。
四段 constant-jerk 的加热时间标度接近 T⁻⁶，尚未证明与论文轨迹定义一致。

探索性补偿：straight 510 μm / 1.6 ms / severity=1，500 shots。
2 倍功率上限下最低有效深度约 118.4→236.8 μK，所需峰值倍率约 2.36。
留阱率两者均为 1，中位激发能增量 15.833→17.647 μK，因此不能声称补偿自动降低加热。
本级没有 SLM 拾取/放回、连续往返或量子内部态。

## 4. Level 4：组合往返与半经典相干性

Protocol A：初始保持 100、拾取 850、split 400、375 μm 去程 1450、远端保持 54、
375 μm 回程 1450、merge 400、放回 850、最终保持 100 μs，共 5654 μs。
SLM 深度 180→60→180 μK，AOD 峰值 280 μK。远端 54 μs 只表示保持窗口，没有实际两比特门。
Protocol B 使用 Level 2 冻结波形反演作为 pickup，正向作为 dropoff，自带 2.4 μm 短移动，
不再次添加 split/merge，总时长 3954 μs；SLM 恒 140 μK。

引擎包含 SLM 三维高斯势、Level 3 AOD 势、势阱归属、首次失败、再俘获区分、分段功-能账本。
重复运行连续继承原子位置和速度。半经典相位积分 φ=∫y(t)ηU/ħ dt，
支持 none、spin-echo、每个长移动段的 XY4，脉冲均为理想瞬时操作。

| 单轮条件 | 成功数 | 成功率 | Wilson 95% 区间 |
|---|---:|---:|---|
| protocol_a_nominal_lensing | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| protocol_a_lensing_off | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| protocol_b_frozen | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| static_idle_same_duration | 2000/2000 | 1.0000 | [0.9981, 1.0000] |

| Protocol A 的 DD 模式 | 单轮条件相干对比度 |
|---|---:|
| none | 0.879740 |
| spin_echo | 0.999958 |
| xy4_per_long_move | 0.969520 |

重复 CSV：Protocol A 第 10 轮仍为 500/500 成功，
C_none=0.0568，C_echo=0.0351，C_XY4=0.0496。
长尾第 40 轮为 128/128 成功。
这些对比度是当前重复序列与统计定义的历史结果，不能与 Level 5 独立抽取 move 的结果直接比较。
联合可用相干份额定义为 P(success)×C_cond，也不是 IRB 门保真度。

## 5. Level 5：有限脉冲、通道与多站点 RB

沿 Level 4 冻结势能轨迹传播 Cs 钟态二能级：
H/ħ = [Δσ_z + Ω(cosφ σ_x + sinφ σ_y)]/2。
Rabi 频率 Ω/(2π)=24.611 kHz，η_SLM=η_AOD=1.3×10⁻⁴。
包含 bare/SCROFULOUS、有限时长 DD、24 元素 Clifford、六输入态 PTM、reference RB、运输插入与 47/195 站点统计。
Reference RB 扫描门数 1–1000；运输插入固定 N=80，M=0–80（包括 M=60），每点默认 72 strings。
默认技术噪声关闭、站点深度/Rabi/失谐散布为零，噪声谱与站点散布属于可选敏感性功能。

以下单通道数值来自旧 channel 输出，未经过本次发现的 IRB 索引路径：

| 条件 | survival | 相对恒等目标的历史 F_avg |
|---|---:|---:|
| protocol_a / none / nominal | 1.0000 | 0.943190 |
| protocol_a / spin_echo / nominal | 1.0000 | 0.333380 |
| protocol_a / xy4_per_long_move / nominal | 1.0000 | 0.957305 |
| protocol_a / xy16 / nominal | 1.0000 | 0.995996 |

目标门必须与脉冲净操作一致。spin-echo 相对恒等门的约 1/3 不能直接解释为失相干。
transformed Clifford frame 尚未进入默认正式传播。冻结回放、忽略态依赖力及多能级泄漏都限制了与实验的比较。

## 6. Level 5 已修复问题与历史结果处理

1. 多站点 IRB 原索引未乘 n_sites 且遗漏 site_index，导致跨序列串写。本轮修为 sequence×n_sites+site。
2. 插入 move 后原先未累计耗时，本轮将 move 的总时长纳入后续门的回放时间。
3. 固定 N、扫描 M 的拟合参数与 reference p(n) 自变量不同，停止套用标准 p_int/p_ref 保真度公式。
4. 估计函数拒绝非有限参数和 F 超出 [0,1] 的结果，保留不可解释原因。

旧输出 F_avg=1.001650 不可解释。旧 IRB、47/195 阵列与相关图片需要全量重跑，
本版不把旧曲线作为修复后的验证证据。原始文件保留，目录新增 RESULT_STATUS.md，旧 summary 顶部添加状态说明。
新增回归使用解析 bit-flip 通道，覆盖 1/3/47 站点和 M=0/1/2/60 的状态、吸收式损失与耗时。

## 7. 距离与次数口径

| 场景 | 单位 | 累计路程 |
|---|---|---:|
| Level 2C n=60 | 60 次 2.4 μm 单程 | 144 μm = 0.144 mm |
| Level 3 L=610 μm | 一次 AOD 长运输 | 610 μm |
| Level 4 一个组合往返 | 2.4+375+375+2.4 μm | 754.8 μm |
| Level 5 默认 M=60 | 60 个完整组合往返通道 | 45,288 μm = 45.288 mm |

Level 5 的累计路程是操作定义的几何计数，不表示同一个经典原子已被连续推进 60 轮。

## 8. 仓库入口、复现与验证

在仓库根目录使用 `python run_simulation.py list`、`check --json` 查看层级与源文件校验。
`run 3 --stage preflight`、`run 4 --stage single_round`、`run 5 --stage preflight` 调用原模块。
`--dry-run` 仅显示命令。默认使用新输出目录，分阶段续跑需指定同一个 --output-dir。
Level 4/5 使用新上游时需显式配置路径；自动发现不保证选择最新产物。

本轮验证记录：`{"status": "passed_scoped_tests", "level3_level4_level5_regressions": {"passed": 89, "deselected": 8, "elapsed_s": 66.78}, "workflow_tests": {"passed": 6, "elapsed_s": 0.02}, "total_passed": 95, "excluded": "slow integration, CLI and full preflight", "production_mc_or_rb_rerun": false}`。
没有在本轮运行完整 2000-shot/47/195 站点正式仿真，也未重做实验参数校准。

## 9. 后续优先工作

先完成修复后的 Level 5 正式重跑及目标门/拟合设计核查，再对齐论文控制波形和绝对 AOD 标定。
将 Level 2C 波形接入三维 Protocol B 需要单独定义 SLM 开关、短转移和长运输的拼接语义，
并做连续轨迹、能量、损失与相位回归；当前报告整合不等于该物理改造已经完成。

## 10. 来源

论文：仓库 `2403.12021v4.pdf` 与 `2403.12021v4.md`；已有图号说明沿用各模块配置中的来源标注。
Level 3/4 结果取自 demo/metrics.json 与对应 CSV；Level 5 单通道取自 demo/metrics.json。
每份数值输入、三层级配置和源代码文件的路径与 SHA-256 见配套 evidence JSON。
