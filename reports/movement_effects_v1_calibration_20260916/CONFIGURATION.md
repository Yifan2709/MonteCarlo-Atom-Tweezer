# 实际运行配置与调用链

本研究入口是此目录的 `run_study.py`，不是仓库旧默认扫描命令。

`final_model.yaml` → `run_study.evaluate` → `continuous_transfer.cli.case_inputs`（显式加载 `configs/continuous_fig6d_unified_T15.yaml`，再加载 `level2c_pickup_survival.yaml`）→ 根据锁定的共同腰重建同一实际位点内的束缚初态 → `waveform_for`/PDF矢量ML → `build_protocol` → `run_noisy_survival` → `propagate`。

- 实际参数优先级为每个case字段覆盖研究入口诊断默认值；最终以 `final_lock.json.parameters` 和每个运行JSON的 `case`、`physics` 为准。`full_base_config`与`full_settings`保留原配置作溯源，不能把其中旧扫描温度数组当成实际温度。
- 四组初态只生成一次，位置、速度、实际SLM位点系数逐元素相同；共用19 uK提议温度和同一拒绝采样规则。不存在逐波形温度、能量截断、热尾或预筛选参数。
- 四组每个原子在每轮结束保持其状态，下一轮继续传播。每个case的所有实际初态及每次检查点的状态均保存在NPZ中；没有重置、补充或冷却。
- 最终噪声由 `case.rin/pointing/slm_rin/slm_pointing`输入，`noise_ledger`保存实际/期望注能。底层旧配置里的Tref虽然仍被读取以计算历史对照功率，但最终 `model=diffusion`不使用这些功率，`Tref_used=false`、`fixed_power_uK_s=null`明确记录此事。
- 共同束腰覆盖SLM及AOD标称值；原2% AOD静态腰散布仍存在。`aod_waist_um`只为敏感性接口预留，最终没有独立调AOD腰。
- 手工持续时间来自曲线名称200/400/600 us；48%二次升深、52%单段三次移动。ML使用 `ml_vector_nodes.json`，400 us、同一时间轴的位置与深度普通三次样条，未经平滑或端点速度夹紧。
- `vs=0`是代码关闭声透镜的约定，物理对应忽略该效应/无限大有效速度极限。不是器件声速测为零。非零情景另存。
- 最终使用三段matched协议：拾取、100 us的AOD等待（SLM关）、反向放下（SLM开）。机械存活计算没有自旋反作用，因此最终 `dd=False`；历史对照保留原DD时间网格但`eta=0`，没有额外力。
- n为单程次数。奇数点判断AOD等待末束缚，偶数点判断放下后SLM束缚；实验散点全部为偶数。最终n=60是30往返。所有图和评分按原论文归一化。
- 运行JSON包含Python/Numpy/Numba版本及生产源码SHA256；确切依赖见 `requirements_exact.txt`。JSON和YAML由各自解析器读取，避免YAML 1.1将JSON中`1e-24`误读成字符串。

正式对照：`formal_baseline`通过显式legacy轨迹与原系数功率复现c4775da行为；`formal_corrected`保留15 uK及固定功率，修正两个系数、手工曲线与ML回放来源；最终采用全高斯随机力。参数表不把这三个模型混为一套同时使用的参数。
