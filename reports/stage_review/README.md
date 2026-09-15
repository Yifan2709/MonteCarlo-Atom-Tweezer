# 综合阶段报告

当前版本为 2026-09-14 v10，43 页。第 33–43 页补充 Level 2C 接入 3/4/5 后的真实连续重算与独立验收。

- [最新 PPT v10](output/AOD_SLM_continuous_level2c_345_2026-09-14_v10_reviewed.pptx)
- [最新文字报告](output/AOD_SLM_continuous_level2c_345_2026-09-14_v10.md)
- [新数值及 SHA-256](output/AOD_SLM_continuous_level2c_345_2026-09-14_v10.evidence.json)

v9 为 32 页历史综合报告，保留备查：

- [PPT](output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9_reviewed.pptx)
- [文字报告](output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9.md)
- [数值来源与 SHA-256](output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9.evidence.json)

## 页码

| 页码 | 内容 |
|---|---|
| 1–18 | Level 2C 原汇报；补充温度与加热参考能量耦合的复核说明 |
| 19 | 各层级的对象、输入、输出与衔接关系 |
| 20–22 | Level 3 三维模型、长运输结果、透镜标定与补偿边界 |
| 23–25 | Level 4 完整协议、单轮成功与相干性、连续重复往返 |
| 26–28 | Level 5 有限脉冲、历史 PTM、RB/IRB 设计及结果状态 |
| 29 | 距离与 n=60 / M=60 的计数口径 |
| 30–32 | 仓库入口、代码修复与验证、物理证据缺口 |
| 33–34 | 连续模型结论、协议与距离口径 |
| 35–38 | 四条波形的原始参考、基线、三维与精细复核曲线 |
| 39–40 | 共同参数扫描、额外长运输结果 |
| 41–43 | 有限脉冲、步长与统计复核、未拟合原因 |

## v10 连续计算与报告生成

在根目录安装项目及可选加速依赖后执行。以下输出目录是本次运行目录，重新计算可更换为其他新目录。

```powershell
python -m pip install -e './simulation/level2_joint_transfer[dev,continuous]'
python run_simulation.py run full --stage preflight
python run_simulation.py run full --stage all --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914
continuous-transfer --stage scan --arms matched --temperatures 10 15 20 25 30 35 --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914
continuous-transfer --stage validation --arms matched --temperatures 15 25 --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914
continuous-transfer --stage validation --shots 128 --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914
python reports/stage_review/continuous_precision.py
python reports/stage_review/analyze_continuous.py
node reports/stage_review/build_continuous_report.mjs
```

`continuous_precision.py` 预先锁定扫描选出的 T=5/15/25、4096 个样本、新 seed=23003 和 dt=0.025/0.0125 μs。`analyze_continuous.py` 核对同一初态池、读取实际计算格点，生成报告、曲线图及结果摘要。不存在的结果会报错，不以历史指标替代。

原 128-shot 粗细步长配对需额外运行 `validation --shots 128`。三维采样器的批大小依赖样本数，不能把 512-shot 文件的前 128 行误认为同一个 128-shot 初态池。

PPT 构建使用下面列出的 `RUNTIME_NODE_MODULES`、`RUNTIME_PYTHON` 和 `PRESENTATIONS_SKILL_DIR` 环境变量。v10 导入 v9，恢复其全部 17 个图表工作簿，并新增 4 个原生存活曲线图表。报告输出和私有渲染路径均与 v9 分开。

## v9 历史报告生成流程

`collect_evidence.py` 使用 Python 标准库读取已有指标及 CSV，不启动仿真。
`build_report.mjs` 使用 `@oai/artifact-tool` 导入 v8，编辑少量说明并添加新页。
原有 14 个图表使用 `restore_chart_workbooks.py` 从 v8 原样恢复工作簿和引用，
新增 3 个原生图表使用实际数据快照。没有虚构工作簿公式或把图表截图当可编辑数据。

在仓库根目录执行：

```powershell
python reports/stage_review/collect_evidence.py
$env:RUNTIME_NODE_MODULES = '<bundled runtime>/dependencies/node/node_modules'
$env:RUNTIME_PYTHON = '<bundled runtime>/dependencies/python/python.exe'
$env:PRESENTATIONS_SKILL_DIR = '<installed presentations skill directory>'
$env:REPORT_PPTX = '<absolute path inside this repository to a new output filename>.pptx'
node reports/stage_review/build_report.mjs
```

使用环境提供的 Node/Python 与 `@oai/artifact-tool`。生成器不覆盖已存在的 PPT，
再次运行应通过 REPORT_PPTX 指定新文件。数字来源按 evidence JSON 核查。
私有构建目录包含逐页渲染与最终检查记录，不属于交付文件。

v9 更新共 95 项指定测试通过，没有重新运行正式 Monte Carlo。v10 的连续重算新增后共 104 项相关检查通过，另有 8 项 slow/CLI/preflight 检查未列入此测试命令。完整 47/195 站点 RB 属于原独立协议，本次未重跑。
Level 5 的旧 IRB/阵列输出已明确标记为需要重跑，PTM 数值也保留历史和目标门的限定。
