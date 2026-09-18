# tweezer-experiment 一期验收报告（A20）

- 日期：2026-09-17；执行环境：macOS arm64，Python 3.14.7，numpy 2.5.2（venv）
- 提取来源：`monte-carlo-for-aodslm` 分支 `movement_effects_v1`，HEAD `1fed6329e95c`
- 工具版本：tweezer-experiment 1.0.0；wheel SHA-256(16) `b184f463a21f4c3f`
- 证据目录：`acceptance/evidence/`（本仓库随源码交付，不指向执行机临时路径）

## 分层结论

| 层面 | 状态 | 依据 |
| --- | --- | --- |
| independent_software_and_io_contract | PASS | A01 干净环境 wheel、隔离与空格路径；A02 API/CLI 一致与导入纯净；A03/A04 契约与有效配置；A18 文档命令逐条执行；A19 空实现扫描 |
| numerical_physics_benchmarks | PASS | A05 独立基准（力=负梯度、解析阱频对比、Verlet 二阶、扩散账本实现/期望一致、零噪声位相等）；A06 初态/吸收语义；A09 模型覆盖；A10 步长/种子 |
| prediction_calibration_for_real_experiments | NOT_RUN | 需要独立装置测量（温度、噪声谱、阱参数、RF 标定）；本仓库全部示例参数为 assumed/synthetic，不做实测验收声明 |
| offline_file_adaptation_under_given_hardware | PARTIAL | A13/A14：合成标定下的映射、量化、限值与回读复算全部通过；真实设备协议未接入，标定表来自用户输入前只能保持 PARTIAL |
| real_device_execution | NOT_RUN | 一期范围明确不含设备连接与运行 |

## 未通过项

无（A02–A19 18 项就地检查全部 PASS；A01 通过；A16 见下表）。

## 需求逐项追踪

| ID | 验收摘要 | 状态 | 预定判据 | 关键证据 | 命令 |
| --- | --- | --- | --- | --- | --- |
| A01 | clean-environment wheel, old-repo isolation, spaced cwd | PASS | fresh venv; wheel built from this source; import from site-packages only; runs from a path with spaces; PYTHONPATH clear | {"import_location": "/private/tmp/tk_fresh_mGBT/venv/lib/python3.14/site-packages/tweezer_experiment/__init__.py"} | zsh acceptance/run_a01.sh |
| A02 | API/CLI same core; import purity; help needs no artifacts | PASS | CLI survival == API survival (same seeds/shots); --help/--version exit 0 without research artifacts; import has no side  | {"api_equals_cli": true} | python -m tweezer_experiment.cli --help; python -m tweezer_experiment.cli simulate examples/short_transfer.yaml --shots 128 |
| A03 | units and invalid inputs | PASS | unit/error tests pass with codes and paths | {'pytest_tail': '13 passed in 0.08s', 'log': '/Users/imok/work/Monte Carlo/monte-carlo-for-aodslm/experiment_toolkit/acceptance/evidence/A03_pytest.tx | python -m pytest tests/test_units_and_schemas.py -q |
| A04 | effective config shows defaults/overrides/sources; params reach the core | PASS | config command emits sources and SI values (unit tests verify wavelength/waist/duration reach the physics and compiled g | {'config_exit': 0, 'has_sources': True, 'reach_tests': 'tests/test_units_and_schemas.py::test_parameter_reaches_physics_core'} | python -m tweezer_experiment.cli config examples/static_hold.yaml |
| A05 | independent physics/numerics benchmarks | PASS | physics benchmark tests pass (finite-diff force, analytic trap frequency, Verlet order, diffusion ledger, mirrors) | {'pytest_tail': '8 passed in 3.07s', 'log': '/Users/imok/work/Monte Carlo/monte-carlo-for-aodslm/experiment_toolkit/acceptance/evidence/A05_pytest.txt | python -m pytest tests/test_physics_benchmarks.py -q |
| A06 | initial states, continuity, observation semantics | PASS | initialization/semantics tests pass (bound states, absorbing/recapture, denominators, per-atom recount) | {'pytest_tail': '27 passed, 12 warnings in 1.76s', 'log': '/Users/imok/work/Monte Carlo/monte-carlo-for-aodslm/experiment_toolkit/acceptance/evidence/ | python -m pytest tests/test_simulate_search_export.py tests/test_waveforms_and_protocol.py -q |
| A07 | non-paper short transfer (350/410 us, 3.1 um, 82 us wait, 13 repeats) | PASS | runs from input only (no source edit); 13*2+1=27 checkpoints; denominator=shots; zero initially unbound | {"checkpoints": 27} | python -m tweezer_experiment.cli simulate examples/short_transfer.yaml |
| A08 | CSV import, AOD transport, combined sequence, generic checkpoints | PASS | all three run with exit 0; checkpoint count == observed segments (+final hold); counts differ across protocols (no 2R+1  | {'runs': {'csv': {'exit': 0, 'checkpoints': 2}, 'long': {'exit': 0, 'checkpoints': 5, 'labels': ['after_pickup', 'after_outbound', 'after_inbound', 'a | python -m tweezer_experiment.cli simulate examples/combined.yaml |
| A09 | noise models cover all protocol types; no silent degradation | PASS | every (protocol, model) pair either runs with the requested model reported in the ledger or raises a coded ToolkitError; | {'matrix': 'evidence/a09_model_protocol_matrix.json'} | - |
| A10 | mechanical/noise step and seed checks | PASS | same seed reproduces bit-identically; different seeds differ; dt and noise-dt differences are reported as numbers (not h | {"mech_dt_0.1_vs_0.05_max_dS": 0.02734375, "noise_dt_0.05_vs_0.1_max_dS": 0.0234375} | - |
| A11 | real constrained search with validation pool; no-feasible outcome | PASS | 8 candidates executed; recommended candidate validated on a separate pool; no-feasible case returns recommended=None | {"no_feasible_case": true, "resume_identical": true} | python -m tweezer_experiment.cli optimize examples/search_pickup.yaml |
| A12 | device/preparation fixed during search; bounds enforced; degradation reported | PASS | single device hash across all candidates; out-of-range variable values rejected (E_BAD_VALUE); statistically tied top ca | {'device_hash_constant': True, 'out_of_range_rejected': True, 'undistinguishable_reported': True, 'undist_note': 'top candidates are not statistically | - |
| A13 | export checks: boundaries, overshoot disclosure, clock/quantization/limits | PASS | SLM switch events recorded at segment boundaries; quantization error reported; limit violation produces NO compliant ins | {"switch_events": 2, "violation_blocked": true} | - |
| A14 | synthetic calibration mapping, read-back re-simulation, missing branch | PASS | fine clock: CSV round-trip lossless and RF round-trip reported; coarse clock difference reported as a number; no-calibra | {"fine_csv_roundtrip_dS": 0.0, "coarse_rf_roundtrip_dS": 0.0} | python -m tweezer_experiment.cli export examples/export_demo.yaml <dir> |
| A15 | reproducibility, cache reuse and invalidation, interrupted-search resume | PASS | identical inputs reuse the stored curve; changed seed or device refuses the cache with a reason; search journal resumes  | {"cache_hit": true, "seed_change": "no stored run for this signature", "device_change": "no stored run for this signature"} | - |
| A05b | waveform/protocol/export/search suites | PASS | remaining contract suites pass | {'pytest_tail': '13 passed, 12 warnings in 1.71s', 'log': '/Users/imok/work/Monte Carlo/monte-carlo-for-aodslm/experiment_toolkit/acceptance/evidence/ | python -m pytest tests/test_simulate_search_export.py -q |
| A17 | measured runtime/memory for representative workloads | PASS | elapsed times and peak RSS measured and recorded (no unmeasured performance claims) | {'measurements': 'evidence/a17_resources.json'} | - |
| A18 | every documented command actually executes | PASS | all README commands (selftest/validate/config/simulate/optimize/export/report/compare) exit 0 with example inputs | {'commands': ['tweezer-experiment selftest', 'tweezer-experiment validate examples/static_hold.yaml', 'tweezer-experiment config   examples/static_hol | - |
| A19 | F1-F5 real chains; no stubs/TODO/hidden degradation | PASS | all five public entry points execute in this acceptance run; source scan finds no TODO/FIXME/NotImplementedError markers | {'api_coverage': {'F1': True, 'F2': True, 'F3': True, 'F4': True, 'F5': True}, 'scan_findings': []} | - |
| A16 | optional paper regression vs archived same-model results | PASS | >=90% of common n within 3.5*sigma_mc AND max|dS|<=0.06 per curve | {'shots': 2048, 'dt_us': 0.0125, 'noise_dt_us': 0.05, 'seed': 20260917, 'criteria': '>=90% of common n within 3.5*sigma_mc AND max|dS|<=0.06 per curve | python acceptance/run_a16_paper.py |

## 验收过程中发现并修复的缺陷（保留记录）

1. 物理控制表速度列单位换算错误（m/s↔µm/µs 因子应为 1，误写 1e-6）——由 A14 新增的“控制保真度”指标暴露（最大相对变化 0.9999），修复后 2.1e-16（ULP 级），全套测试与验收重跑通过。
2. SearchJournal 读取头部行缺 index 键（A15 断点续跑暴露），已修复。
3. RunStore.save 默认采用结果自带的 seed/shots 而非 bundle 默认值（A15 暴露），已修复。
4. 论文 ML 回放适配器的放回段误用解析 paper_manual，而归档研究协议是 ML 波形自身的精确时间反演（A16 首跑暴露：ML 曲线 ΔS=0.083）；改为反向 CSV 回放后四条曲线全部落入判据（max|dS|≤0.018）。三条手工曲线在修复前后均通过，证明差异隔离在适配器层而非物理内核。

## 已知边界（如实声明，不因验收通过而隐去）

- 带随机噪声时，任何 ULP 级控制差异都会被混沌放大：回读 ΔS（如 0.008@96 原子）包含该放大，导出链无损性以控制保真度为准（<1e-12 判据通过）。
- 一期噪声为逐原子独立流：MC 区间不表示同批实验的共模波动。
- 论文数据仅通过 `adapters.paper` 作为显式用户资源进入（A16）；删除该资源后其余全部验收不受影响。
- 无真空损失/成像效率模型；`cubic` CSV 插值需要 scipy（缺失时显式报错而非降级）。
