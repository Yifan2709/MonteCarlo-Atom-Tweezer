# AOD–SLM 光镊 Monte Carlo

主工程位于 `simulation/level2_joint_transfer/`，保留 Level 0–5 与 Level 2C。
2026-09-14 更新补齐统一入口、Level 3/4/5 说明、综合报告，以及 Level 5 IRB 的索引与适用性修复。

2026-09-15：论文 Fig. 6d 参考已改为独立提取 PDF 矢量拟合线、散点和误差条。连续模型新增按实际 SLM 格点阱深准备束缚初态，修复先按名义阱采样、再加入位点差异而混入初始不束缚原子的问题；支持独立设置初始温度与参量加热参考。既有 PPT 使用修复前归档，不能当作本次初始化修复后的结果。600 μs 的完整曲线检查脚本为 `simulation/level2_joint_transfer/tools/study_manual600_initialization.py`，参数扫描属于诊断校准，不代表已经复现实验。

云端版本包含源码、配置、测试、参考数据、报告和轻量结果；`tmp/`、临时图表目录、报告构建缓存及连续模型的逐原子 `.npz` 轨迹仅保留本地。需要逐原子诊断时，先根据对应运行计划重新生成轨迹。

| 层级 | 计算内容 | 一次操作与距离 | 当前依赖 |
|---|---|---|---|
| Level 2C | 一维拾取/放回、噪声平均加热、重复转移存活 | 单程 2.4 μm；60 单程为 30 短往返 | Level 2 引擎 |
| Level 3 | 三维长运输、柱面透镜、机械激发与留阱 | 单次 270/510/610 μm | 三维 AOD 势与热初态 |
| Level 4 | 完整组合往返、连续运动状态、半经典相位及理想 DD | 每轮 2.4+375+375+2.4=754.8 μm | Level 2 冻结波形及 Level 3 透镜模型 |
| Level 5 | 有限微波脉冲、PTM、reference RB、固定 N 的运输插入、47/195 站点 | 默认 M 为完整组合往返数，每个 move 独立抽取冻结轨迹 | Level 4 冻结协议与轨迹池 |

**新增 `full` 连续入口，已接入 Level 2C 的四条波形、三维柱面透镜、重复往返与有限 DD/PTM。**每个原子的六维运动状态及内部态连续保留，失阱后不再返回样本池。原有各层级 CLI 保留，原 Level 5 的独立冻结轨迹 RB 不等于这个连续模型。

`matched` 保持 Fig. 6d 的短转移协议，用于原 Level 2C 验收。`extended` 另外加入 375 μm 去程、54 μs 远端保持及 375 μm 回程，单轮名义路程 754.8 μm。两者的 `n=60` 都表示 30 轮中的 60 次短拾取/放回计数。扩展协议总路程约 22.644 mm，不可直接作为 Fig. 6d 拟合结果。

## 安装和统一运行

需要 Python 3.10 或以上版本。在仓库根目录执行：

```bash
python -m pip install -e "./simulation/level2_joint_transfer[dev,continuous]"
python run_simulation.py list
python run_simulation.py check --json
python run_simulation.py run 3 --stage preflight --dry-run
python run_simulation.py run 3 --stage preflight
python run_simulation.py run 4 --stage single_round --output-dir simulation/level2_joint_transfer/outputs/level4_end_to_end_roundtrip/new_run
python run_simulation.py run 5 --stage preflight
python run_simulation.py run full --stage preflight
python run_simulation.py run full --stage all --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/new_run
```

也可使用安装后的 `simulation-workflow` 命令。未指定输出目录时，每次运行使用新的时间戳目录。
按阶段续跑时应显式使用同一个 `--output-dir`。`--dry-run` 只显示命令，不创建目录或启动仿真。
正式全量计算可选 `--stage all`；Level 5 会生成大型轨迹池，运行成本远高于单元测试。
Level 4/5 的上游文件仍由各自 YAML 指定；使用新上游结果时应显式配置路径，不能依赖自动发现选中最新目录。

| 层级 | 可运行阶段 |
|---|---|
| full | preflight / scan / validation / convergence / all |
| 3 | preflight / potential_scan / coarse_mc / validate / sensitivity / compensation / all |
| 4 | preflight / single_round / repeated_rounds / ablations / sensitivity / all |
| 5 | preflight / channel / reference_rb / interleaved_rb / array_extension / sensitivity / all |

## 报告与结果状态

- [统一初始温度 Fig. 6d 重算与参数审计（2026-09-16）](simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260916_unified_T15/REPORT.md)：四条曲线共同 T0=Tref=15 μK、修复后按实际格点初始化、4096 原子、dt=0.0125 μs、冻结 seed 23003；四条曲线均未通过 ±0.02 容差，逐曲线参数审计与复现命令见报告与 `param_audit.csv`。
- [独立 Level 5 PPT：论文参考数据修正版](reports/level5_paper_review/output/Level5_MonteCarlo_论文对照与全流程_2026-09-15_参考数据修正版.pptx)
- [手工 600 μs 初期损失诊断（修复前归档）](reports/level5_paper_review/output/manual600_early_audit/手工600us初期损失核查.md)
- [连续全流程重算报告 v10](reports/stage_review/output/AOD_SLM_continuous_level2c_345_2026-09-14_v10.md)
- [连续全流程 PPT v10](reports/stage_review/output/AOD_SLM_continuous_level2c_345_2026-09-14_v10_reviewed.pptx)
- [连续仿真配置](simulation/level2_joint_transfer/configs/continuous_level2c_345.yaml)
- [本次逐样本结果](simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914)

v10 包含新的 Monte Carlo 重算。下面的 v9 及独立模块 demo 是此前的历史对照：

- [综合报告 v9](reports/stage_review/output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9.md)
- [综合 PPT v9](reports/stage_review/output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9_reviewed.pptx)
- [数值来源及 SHA-256](reports/stage_review/output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9.evidence.json)
- [报告生成说明](reports/stage_review/README.md)
- [各模块详细说明](simulation/level2_joint_transfer/README.md)
- [Level 5 历史结果状态](simulation/level2_joint_transfer/outputs/level5_finite_pulse_irb/demo/RESULT_STATUS.md)

独立 Level 3/4 的 v9 数值来自已有 demo 输出。v10 的 full 结果使用新的连续仿真输出。
Level 5 的历史 IRB/阵列输出受多站点索引错误影响，修复后需要重新生成。旧输出中 `F_avg=1.001650` 无物理意义，新报告仅将其列为待纠正记录。
Level 5 单通道 PTM 与 IRB 是独立路径；单通道数字仍是相对指定目标门的历史模拟结果，不能直接当作论文实验保真度。

## 验证

```bash
cd simulation/level2_joint_transfer
python -m pytest tests/test_level3.py tests/test_level4.py tests/test_level5.py tests/test_level5_irb_regression.py tests/test_simulation_workflow.py tests/test_continuous_transfer.py tests/test_continuous_initialization.py -m "not slow" -k "not cli and not preflight" -q
```

该命令覆盖数值模型与新增索引/耗时/拟合适用性回归。大规模 CLI、完整 preflight 与正式 47/195 站点仿真另行运行，不能用单元测试通过替代物理校准或实验验收。
