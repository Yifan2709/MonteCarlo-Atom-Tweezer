# 两篇论文接入 Monte Carlo：实现与验证报告

**状态总览（先讲匹配程度与失败位置）**：

- 结构/锚点级通过：R1 R2 Y1 Y2 Y3 Y4（预注册结构判据全过或差距<容差）。
- 定性/序关系通过、水平差距如实报告：R3 Y5 Y6（见下表差距列）。
- 基础设施：B0 M1 I1 V1（B0 含本版初态束缚实检 D8；Cs 600μs 历史偏差仍未解决且未被改写）。
- **独立验证 A 阶段修正已并入本版**（validation_plan.md，D1-D9）：Y6 后选口径改只读检测；Y2/Y4 标签降级为自洽；I1 声明标量耦合；sixth_zero_jerk→low_peak_v；统一入口语法修复（原 0b86619 该入口损坏）。
- **不能**声称“完整复现两篇论文”：绝对曲线验收所需的论文源数据（R 图1d/ED2、Y ED2b）在本环境未获取（网络受限），Y1/R1/R2 仅做结构判据；R3 raw 差 0.17、Y6 水平差约 +0.11，归因见 §差距。

| ID | 实现 | 数值 | 论文对照 | 关键数值 |
|---|---|---|---|---|
| B0 | done | ok | n/a(基线回归) | 初态束缚=0 确认；600μs 历史偏差 7.28pp 原样保留（未解决，未改写） |
| M1 | done | ok | n/a(基础设施) | 账本/键控随机流自检通过；24 项单元测试覆盖 |
| R1 | done | ran | structural_pass | C1=True slope=-4.338334681716283 anchor=1.0 |
| R2 | done | ran | structural_pass | platform=0.9526(paper raw 0.948) DD on/off=0.9580/0.9150 |
| R3 | done | ran | qualitative_pass | raw=0.545(paper 0.71, gap 0.165) cond=0.769(paper 0.991) accept=0.673(paper 0.66) |
| Y1 | done | ran | structural_pass | C1=True C2=True C3=True |
| Y2 | done | ran | construction_consistent(循环锚定已声明) | optimal/trip=8.88e-04(bound 6e-4, g_optimal 由界反求→自洽非复现) |
| Y3 | done | ran | anchor_pass | RT1=0.0127(paper 0.011) frac RT1/RT5=0.56/0.07(paper 0.5/0.1) |
| Y4 | done | ran | self_consistent_model_sampling(构成锚点通过) | slope AR/TO=4.00/2.00(预置公式自洽)；erasure AR=0.38(paper 0.38,MC 抽样锚点) |
| Y5 | done | ran | ordering_pass | prep raw/flag/ps=0.947/0.958/0.948(paper 0.981/0.990/0.995); storage t=0 uncond/informed=0.947/0.965(paper 0.874/0.902) |
| Y6 | done | ran | ordering_pass | adaptive/random/ps=0.912/0.894/0.912(paper 0.802/0.771/0.87; ps 已改只读检测口径 D3); 水平差 ~+0.110（见报告归因） |
| I1 | done | ran | n/a(标量耦合接通_pass) | mode=scalar_statistic_cou… yb=True rb=True |
| V1 | done | ok | n/a(验收) | 本文件与 verify_paper_integration.py 输出 |

## 差距与未匹配原因

1. **R3 raw ⟨X̄⟩ 0.54 vs 0.71**：本实现 dark 比特 33%（初态损失 1%+3 层 CZ loss+运输层+检错轮），论文含未公开的线路结构与 SPAM 细节；接受率 0.67 与论文 0.66 一致，条件 ⟨X̄⟩ 0.77 vs 0.991 差距来自读出误差 1%/比特与未检出 X 错。
2. **Y5 制备 0.947/0.960 vs 0.981/0.990**：每门误差取论文 AR CZ 1.6%，但制备线路深度（3 CX）为结构假设；论文线路未公开。
3. **Y6 传送 0.91 vs 0.80**：本实现读出含单擦除恢复；论文含未建模的测量层与更多 SPAM。序关系（自适应>随机、后选最高）与论文一致。
4. **Y1 绝对曲线**：τ_eff 为 assumed（扫描 {0,2,4,8}），论文 τ_AOD 值不可得；结构判据（透镜限制快速运输、含透镜 zero-jerk 更优）成立。
5. **Zenodo/PDF 源数据未获取**：所有参考锚点取自 arXiv 正文/Methods/ED 文字数值，未做曲线数字化。

## 复算入口

```bash
cd simulation/level2_joint_transfer
.venv/bin/python -m paper_integration.cli run --stage preflight --output-dir outputs/paper_integration/recheck
.venv/bin/python -m paper_integration.cli run paper-rb2022 --stage scan --shots 1024 --seed 20260915 --output-dir outputs/paper_integration/recheck
.venv/bin/python -m paper_integration.cli run paper-yb2026 --stage scan --shots 1024 --seed 20260916 --output-dir outputs/paper_integration/recheck_yb
.venv/bin/python -m paper_integration.cli run paper-rb2022 --stage validation --trials 200000 --seed 20260917 --output-dir outputs/paper_integration/recheck
.venv/bin/python -m paper_integration.cli run paper-yb2026 --stage validation --trials 200000 --seed 20260918 --output-dir outputs/paper_integration/recheck_yb
.venv/bin/python -m paper_integration.cli run --stage integration --output-dir outputs/paper_integration/recheck
.venv/bin/python tools/verify_paper_integration.py --contract configs/paper_integration/contract --run-dir outputs/paper_integration/run_formal
.venv/bin/python tools/build_paper_report.py
.venv/bin/python -m pytest tests/test_paper_integration.py -q
```

## 配置体系与参数溯源

cs_existing / rb87_2022_native / yb171_2026_native 三套配置在
`paper_integration/species.py` 登记，逐参数来源见
`configs/paper_integration/contract/acceptance_contract.json` 与
运行目录各 run_manifest.json 的 config 字段。

## 遗留（未完成/受限）

- ACCEPTANCE_PROTOCOL.md 在仓库中不存在，按任务书 §3/§4/§8 执行并登记。
- 工作目录差异：任务书路径为 D:/，本实现基于本机仓库（同 HEAD b417aaa）。
- extension 项（R 其余图态、surface/toric code、混合模拟；Y 全部主图逐项）未实现，属 full_paper 范围。
- 60 μs 600 μs Cs 已知拟合偏差（±0.02 未通过）为既有问题，本任务未触碰。