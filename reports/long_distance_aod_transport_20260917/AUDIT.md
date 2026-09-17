# 生产链审计及修改

审计基线 `1fed6329…`；原分支、旧结果未覆盖。历史报告只帮助定位，公式/协议由本地原论文页与作者Zenodo源数据复核。

|模块|实际入口与已有能力|已验证范围|确认问题|本研究处理|
|---|---|---|---|---|
|level3_3d_transport_lensing|独立level3 CLI，三维高斯/双柱面势、Verlet、功账本|解析力已用独立数值势梯度核对；零透镜、不同轴符号另测|constant_jerk是四段S曲线，非Rb cubic；v_s为情景等效速度；能量接口是实验室系|保留旧轨迹和历史输出；新核精确命名并只在静止观测点判损失|
|continuous_transfer.engine|统一full→continuous.cli→run_continuous→integrate；连续六维/自旋|跨轮状态保存，实际位点初态逻辑可复用|CLI仍传固定加热功率；该核不同于短途标定调用的随机核；谱未知不能用旧通过认证|新研究入口明确调用long_distance_transport.model，直接复用field与beam_noise；不绕回固定功率入口|
|continuous_transfer.noisy_survival|历史run_study脚本直接调用；局部力RIN/指向扩散|力梯度、谐振平均/方差及新步长核查|接口只允许matched短途；不能把short匹配结果变成长程标定|新控制序列支持源/目标两个静态位置与重复方向，固定同一PSD|
|continuous_transfer.initialization|位点深度先抽，按实际位点拒绝采样|保持原位点ID和束缚初态；本研究独立实现同规则并保存初态|函数硬编码SLM1061nm，不能直接跨物种复用|新通用制备显式波长，按A或B实际初始势抽样，不重新采样存活者|
|paper_integration.trajectories|论文入口实际调用|新公式直接对照作者ED2c数据|旧五阶族错误，6,-7,2不能代表zero-jerk；gamma1.875也由错误族反求；constant_accel位置有跳变|替换为精确S2，显式c2–c5；三角速度恒加速度段修为连续；时间网格必达端点|
|paper_integration.dynamics3d|Rb/Yb论文入口共同积分|修后回归及新独立参照|二次拒绝采样布尔索引长度错误；dt使用传参而非段网格；末能量漏时变SLM/透镜；移动判据用实验室动能|索引限制pending，逐步实际dt；单阱运动用相对速度，移动双阱拒绝判定；末能量按末控制|
|paper_integration.yb_transport|旧Y1 CLI调用|旧结构C1–C3只测试趋势，不是论文复现|默认加入交接，与ED2不同；交接先独立推进，又在总轨迹里推进第二次；初始化把两个阱叠加；hash(traj)种子跨进程不稳定|默认改纯AOD；显式旧交接路径只推进一次；单阱束缚初态；成对固定种子|
|paper_integration.rb_transport|原生Rb单次函数|方法比较与独立谐振计算|把Eq2称常加速度；110um误作单原子距离；解析Nmax误标物理边界；原始ΔN包括热基线|新扫描单原子55um；Eq2文档纠正；不以Nmax/erf制造存活曲线|
|paper_integration.species|参数登记，不是测量数据库|逐项与正文/图/数据对照|Yb腰0.77标原文而无法定位；30kHz从周期调制借用于运输；距离50um原标假设现有源数据支持|PARAMETERS明确降级/更新证据；旧跨逻辑模块数值不全局改动，避免无关系统语义变化|
|level4_roundtrip_coherence|独立CLI，完整时变双阱、相位累积|结构上跨轮连续，能量与相位有历史预检|透镜severity和差分光移/噪声未成为本装置标定；标签检查点与激发时点不同|不把旧自旋结果作为长程量子保真度；新记录每阶段能量|
|level5_finite_pulse_irb|独立CLI或frozen replay / joint验证|有限脉冲数学能力，不证明装置相干模型|默认回放冻结经典池，不会自动响应新机械模型/参数；旧success保持不变|本任务不运行完整IRB/纠错；机械与量子态保持明确分开|

新链：`run_study.py → Config → controls/initial → integrate → field/beam_noise → 每单程静止检查 → 原始NPZ/CP区间`。与历史CLI并列且显式调用；所有正式结果由该链实际计算，不是孤立演示函数。

验证的证据是benchmarks.json中的数值残差、独立DOP853和光学势差分、作者导数数据及逐点收敛记录，不是此表的自填状态。旧回归日志保留失败及原因；旧CLI不得因新研究成功而宣称所有入口均已修复或标定。
