# 复现

入口 `reproduce.ps1`。Windows PowerShell，在仓库根目录执行：

```powershell
& ./reports/long_distance_aod_transport_20260917/reproduce.ps1 -Threads 8
```

它先验证物理，再逐阶段计算和绘图。运行较久，尤其8192原子连续50次及细时间步。队列在各阶段第一次执行前写入`*_queue.json`；`confirmation_selection.json`冻结初始确认条件。后续final_single使用另外两组26091741/42种子检查99%边缘、粗扫描零损失点之外的余量点和慢端。正式物理参数没有依据验证结果重拟合。

断点续算根据配置SHA256的前18位命名，只跳过已经写完JSON的点；同名NPZ先写，JSON最后写。要完全独立从零重算，请复制源码/数据到新的研究目录，保留脚本与论文数据，使用空`runs`目录；不要删除原始结果。恢复时须核对源文件hash，不能把改过物理的源码与旧完成点混合。本次运行中仅扩展队列/分析脚本；正式动力学模型冻结后未更改。

从零重算时也应保留本次冻结的全部`*_queue.json`及`confirmation_selection.json`。最终调度器优先读取这些权威配置，不因后来加入验证种子导致报告边界移动而重新挑选旧阶段条件。若要研究全新条件，使用新阶段名或独立目录，而不覆盖旧队列。最终汇总将相同正式步长下匹配的独立验证种子统一加入已确认条件；其他步长不增计样本，合并前统计另存。这项审查后汇总改进没有改变成功标准。

环境原先有两套本地依赖目录：`tmp/audit_deps`和`tmp/continuous_deps`，不上传、不视作源码。迁移环境可按`requirements_exact.txt`安装依赖，再令PYTHONPATH含`simulation/level2_joint_transfer/src`。Python/NumPy随机实现和Numba版本必须记录；参考`environment.json`。机器线程数不改变每原子显式种子。浮点平台变化仍可能改变混沌边界个体标签。

后续创建了本地 `tmp/long_distance_venv` 并对 `simulation/level2_joint_transfer` 执行可编辑安装，以满足旧CLI子进程清理PYTHONPATH后的导入；默认复现命令使用此环境。迁移时可另建虚拟环境、安装冻结依赖并执行 `python -m pip install --no-deps -e simulation/level2_joint_transfer`，然后把该Python路径传给 `reproduce.ps1 -Python ...`。实际版本（含pypdfium2）也记录在 `provenance.json`。

核心输出：每点JSON保存完整参数、运动学、逐程存活和CP区间、耗时、源文件hash；NPZ保存全部原子的初始状态、实际位点因子、拒绝次数、每程六维状态/能量/存活、每阶段能量及噪声功账本。它们足以用配置、种子与源码重新积分任意原子的详细时间轨迹，未仅保存均值。NPZ文件保持本地完整，不要求把大型数据纳入Git；最终文件清单含SHA256。

作者原始Yb数据：Zenodo 19491381，CC-BY-4.0；zip MD5=4701e05c2bb8abac833a586d8f09e428。保留整个数据包，但仅使用运输相关图。Rb原始图在指定本地PDF第17页，向量散点提取算法与坐标误差见`papers_and_optics.py`、`rb_digitization.json`。图线插值仅用于对照误差分析，绝不用于生产动力学或失阱判定。

`benchmarks.log`、`native_coarse.log`、`stages.log`、`confirmation.log`、`operations_confirmation.log`、`final_single.log`保存实际运行记录。`regression.log`与`regression_retry.log`保留旧CLI嵌套预检查失败；不能只保留最终成功记录。所有统计置信区间都是条件模型内逐点区间，不能消除参数不确定性或变成全扫描同时置信。

最后补充的 `native_steps`、`boundary_steps`、`repeat_steps` 与 `edge_robustness` 均为独立检查，未改生产动力学。`run_study.py --start N --stop M stage` 支持零起点的互不重叠队列分片；部分尾部队列使用单独进程计算，完成的JSON会被原顺序进程跳过。没有选择性删除失败点。全配置JSON是执行记录的依据，不依赖恢复旧版调度脚本。

尾部分片与原顺序进程有少数同条件的重复执行，日志均保留；配置哈希相同的确定性记录共用同一个文件，统计按条件和种子去重，不把重复计算增加为独立样本。物理源码哈希按文件原始字节检查；迁移时若Git改变行尾，使用`source_snapshot.zip`还原精确文件字节，避免把行尾变化误当成可忽略的源码差异。

`verify_artifacts.py` 逐一检查全部原始NPZ与JSON的原子数、束缚标签、跨轮失效状态和冻结物理源码哈希，并原样重放一个正式条件；它证明归档和确定性，不能替代实验验证。`package_evidence.py` 保存源代码ZIP、Git源提交号、原始PDF与作者数据哈希、依赖版本，生成覆盖本目录所有交付文件（除自身清单及缓存）的 `artifact_manifest.csv`。`source_snapshot.zip` 包含生产源码、测试及研究脚本。原始NPZ约2.5GB，不进Git但仍保留在本目录；完整迁移需一并复制，不能只复制Git仓库。

归档核对中的`atom_trials`是逐文件样本数总和，包含配对初态及跨阶段重复计算，**不是相同数量的独立样本**。正式区间只按`independent_confirmation.csv`中的条件/种子去重后计数，`source_files`给出合并所用的全部原始JSON。

条件损失及合并能量分位数以最终CSV为准：前一轮存活数为0时条件损失未定义，不能把冻结原始JSON旧占位值1解释为新增100%损失；合并能量分布从全部对应NPZ重算，不使用代表种子的分位数。该统计层修正未改动任何原子轨迹或存活计数。
