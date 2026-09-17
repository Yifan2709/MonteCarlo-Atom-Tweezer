# 研究状态与证据索引

“实现”“计算”“验证”不等价。下面链接到外部数据、原始数值和具体反例；状态文字本身不作为证明。全部队列的实际完成数由 [queue_completion.csv](queue_completion.csv) 自动读取结果生成。

|研究项|已实现|已计算|已验证范围|仍缺证据 / 不能声称|
|---|---|---|---|---|
|远端与独立分支|保留基线，独立分支，无推送|远端GitHub分支API读取|[远端响应](remote_branch_metadata.json)、Git提交记录|Git内置HTTPS helper缺失，采用官方API核查；未改原分支|
|两篇原文运输协议|Rb单原子55μm、Yb ED2纯AOD；区别完整交接|原文定位与作者数据读取|[原文页](paper_pages/)、[作者数据](zhang_source/dataset/)、[来源说明](paper_data_metadata.json)|Rb数值源数据未获得；Yb无作者MC完整输入|
|精确论文轨迹|三次与准确S2、显式命名|端点/连接、全导数/极值|[独立作者导数残差](benchmarks.json)、源码回归|未知硬件带宽如何平滑实际连接|
|生产调用链审计|逐入口跟踪，单独真实研究入口|旧实现/修正核对|[AUDIT.md](AUDIT.md)、[代码修改测试](new_physics_tests.xml)|不认证旧full固定功率核、旧相干/IRB为本装置正确|
|实际位点束缚初态|先位点、后拒绝；保留所有初态/拒绝次数|所有运行NPZ|原始`initial_pos/vel/site_factors/attempts`，初始化测试|实际三维制备分布、特别热尾未测|
|三维高斯及声透镜力|复用实际field，显式两焦点|生产扫描与独立差分|[benchmarks.json](benchmarks.json) 的力/指向导数残差|理想高斯及成像模型的实装精度|
|运动能量参考系|运动不按实验室K吸收；静止末点判E<0|匀速与阶段/末点能量|[Galilean残差](benchmarks.json)、[阶段诊断](loss_timing_diagnostics.csv)|观测时刻不是不可逆飞走时刻；真实成像损失判据|
|机械频谱|精确加速度傅里叶式|DOP853独立解与Rb近似域|[谐振对照](benchmarks.json)、[ωT域](rb_harmonic_averaging.csv)|线性谐振近似不直接预测有限深三维失阱|
|有限恢复力与势垒|倾斜势稳态根/鞍点|深度/轨迹/加速度对照|[inertial_barrier.csv](inertial_barrier.csv)、[深度图](depth_acceleration.png)|瞬态越过临界不等于立即丢失|
|声透镜机制|τ、轴符号、焦移显式控制|零/单轴/双轴/反号|[lensing_geometry.csv](lensing_geometry.csv)、[存活对照](geometry_survival.png)|τ、真实晶体束径、RF标定未测|
|声学有限传播|完整延迟相位及二次项比较|余项/模态重叠|[acoustic_propagation_check.csv](acoustic_propagation_check.csv)|未将全衍射场接入MC；快速100μs附近瞬时近似不足；不认证实际99.9%|
|随机技术噪声|单边/Hz力扩散，无经验loss/heating|谐振均值/方差、全势功账本|[基准矩](benchmarks.json)、[生产噪声](production_noise_benchmark.json)|实测频谱、SLM噪声和阵列共同噪声相关性缺失|
|完整交接B|明确200μs共址线性交叉渐变|1/10/50及零移动对照|[handoff_controls.csv](handoff_controls.csv)、全程原始状态|不冒充用户实测波形；有限波形族非全局优化|
|跨轮连续性|状态及每原子随机流连续，不重采样|原始状态/阶段能量，双程vs分段|[continuity残差](benchmarks.json)、`runs/*/*.npz`|不包含实际轮间冷却、纠错或检测反馈|
|Rb原生运输对照|原生87Rb；无目标曲线输入|24数字化点、20重合点比较|[逐点偏差](rb_comparison_residuals.csv)、[时间步](native_convergence.csv)|RMSE约0.100；**未定量复现实验**；初温/透镜缺少独立输入|
|Yb原生运输对照|原生171Yb；准确两轨迹/纯运输|作者源数据与τ共同情景比较|[逐曲线误差](paper_comparison_metrics.csv)、[时间步](native_convergence.csv)|部分参数是假设；**未定量复现实验**；不逐曲线拟合|
|Cs第一阶段|固定装置/制备规则；D×T×轨迹×深度|跨稳定、过渡、失效区；粗网格+细化|[全部逐点](all_points.csv)、[距离时间图](distance_time_region.png)|真实扫描范围/动态功率补偿未认证|
|Cs第二阶段|逐机制加入PSD/误差/B|保持无噪声对照|[mechanisms.png](mechanisms.png)、逐点JSON|PSD/误差只是共同候选值|
|Cs第三阶段|重复和预设参数敏感性|270/510/610μm、1/10/50次；边缘扰动|[重复图](repeated_transport.png)、[敏感性](sensitivity.png)、[边缘扰动](edge_robustness.csv)|敏感性范围不是测量置信区间|
|抽样及独立种子|精确二项区间，至少两组新种子|关键点各4096×2或更多|[分种子](seedwise_confirmation.csv)、[合并区间](independent_confirmation.csv)|跨目标标未确定；无多重比较同时保证；不抵消模型误差|
|时间步误差|机械与噪声步单独比较|0.1/0.05/0.025μs机械；0.2/0.1/0.05μs噪声|[机械](mechanical_convergence.csv)、[噪声](noise_step_convergence.csv)|重复边缘部分不精确收敛；个体标签差异不等于整体分布差异|
|边界提取|非单调网格、只认证已算点|每距离/轨迹/重复数的表|[BOUNDARIES.md](BOUNDARIES.md)、[完整支持点](boundaries.csv)、[非单调证据](nonmonotonic_screening.csv)|不是连续区间证明、全局最优或普适速度极限|
|能量与损失时点|每阶段能量、每程状态/分布|逐次条件损失与残余质心/速度|[residual_excitation.csv](residual_excitation.csv)、[能量图](energy_conditional_loss.png)|幸存者变冷可能是筛选；触发阶段不等于逃逸时刻|
|散射/量子态|与机械存活分开；保留Yb来源数据|文献指标及源波形累计量核对|[paper_data_metadata.json](paper_data_metadata.json)、原文相干页|**没有可信Cs量子态速度上限**；不混用Yb寿命或Rb退相干保护|
|复现与失败保留|配置、种子、源码哈希、日志、原始NPZ|全点输出和清单校验|[REPRODUCE.md](REPRODUCE.md)、[全原始记录核对/重放](artifact_verification.json)、[哈希清单](artifact_manifest.csv)、[源码来源](provenance.json)|约2.5GB原始NPZ保留本地并列清单，不放进Git；转移研究必须一并复制|

结论强度：可验证的方法和条件性控制计算已经形成；实验一致性、真实装置参数和量子态成功仍不能认证。这不是将未完成项改名为“通过”。
