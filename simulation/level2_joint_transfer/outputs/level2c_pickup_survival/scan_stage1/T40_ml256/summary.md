# Level 2C SLM→AOD pick-up 存活率复现（论文 Fig. 6d 顶图）

范围：仅转移存活率（atom survival vs 单程转移次数）；不含 IRB 保真度、XY4、相干性、成像保真度与长距离运输。对标数据见 reference/fig6d_parameters.md。


## ML 轨迹（14+14 控制点三次插值）

- 优化器：LHS + local perturbation refine（替代论文在线 ML 优化器，等价黑盒搜索）；目标：60 次连续转移存活率。
- 冻结候选 ml_0000：训练池存活率 0.6680。
- 独立复核（512 shots）：ML S(60)=0.6426，手工 400 µs S(60)=0.6289；配对 bootstrap 差 +0.0137 CI [-0.0254, +0.0527]。

- 程序化图片检查：2/2 张全部通过。
- 未进行人工视觉审阅；程序化检查不能证明标签无重叠或版面美观。
