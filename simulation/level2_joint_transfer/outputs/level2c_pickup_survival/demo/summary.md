# Level 2C SLM→AOD pick-up 存活率复现（论文 Fig. 6d 顶图）

范围：仅转移存活率（atom survival vs 单程转移次数）；不含 IRB 保真度、XY4、相干性、成像保真度与长距离运输。对标数据见 reference/fig6d_parameters.md。

- 预检查：8 项关键检查全部通过（176 s）。

## 单次 pick-up 基线（保守模型，三池）

| 波形 | 池 | 俘获 | 俘获率 | Wilson 95% CI |
|---|---|---:|---:|---|
| pickup_pickup_200us | optimization | 128/128 | 1.0000 | [0.9709, 1.0000] |
| pickup_pickup_200us | selection | 512/512 | 1.0000 | [0.9926, 1.0000] |
| pickup_pickup_200us | validation | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| pickup_pickup_400us | optimization | 128/128 | 1.0000 | [0.9709, 1.0000] |
| pickup_pickup_400us | selection | 512/512 | 1.0000 | [0.9926, 1.0000] |
| pickup_pickup_400us | validation | 2000/2000 | 1.0000 | [0.9981, 1.0000] |
| pickup_pickup_600us | optimization | 128/128 | 1.0000 | [0.9709, 1.0000] |
| pickup_pickup_600us | selection | 512/512 | 1.0000 | [0.9926, 1.0000] |
| pickup_pickup_600us | validation | 2000/2000 | 1.0000 | [0.9981, 1.0000] |

保守模型下单次 pick-up ≈100% 为预期的无害基线（与 SI “200 µs 单程 >98%”兼容，实验值含本模型未包含的损耗）。

## 保守模型 survival(n≤60)（扫描池）

- manual_200us：p0=0.9503，p=1.00000，b=0.0240
- manual_400us：60 次转移全程无丢失；p 的 Wilson 95% 下界 0.99975（总试验 15360）
- manual_600us：60 次转移全程无丢失；p 的 Wilson 95% 下界 0.99975（总试验 15360）

保守模型 S_n≈1：曲线形状机制正确但数值对不上论文——正常，差距归因进入 survival_scan 阶段。

## ML 轨迹（14+14 控制点三次插值）

- 优化器：LHS + local perturbation refine（替代论文在线 ML 优化器，等价黑盒搜索）；目标：60 次连续转移存活率。
- 冻结候选 ml_0000：训练池存活率 0.9844。
- 独立复核（512 shots）：ML S(60)=1.0000，手工 400 µs S(60)=1.0000；配对 bootstrap 差 +0.0000 CI [+0.0000, +0.0000]。

## 损耗机制归因扫描与论文对比

- 扫描变体：conservative, heating, disorder, full；正式对比变体：full（复核池）。

| 曲线 | 最大偏差 | 容差带内比例 | 是否匹配论文 |
|---|---:|---:|---|
| manual_200us | 0.144 | 0.00 | 否 |
| manual_400us | 0.516 | 0.00 | 否 |
| manual_600us | 0.423 | 0.00 | 否 |
| ml_400us | 0.259 | 0.09 | 否 |

验收标准与诚实差距归因见 metrics.json 的 honesty_note 与 comparison_to_paper.json；未为凑数调整任何噪声/差异参数。

## 2D 横向模型检查（条件触发，prompt 第 12 步）

- 触发原因：1D+加热+位点差异未达论文损失水平，按 prompt 第 12 步升级到 2D 横向模型。
- manual_200us：2D S(60)=0.2578，对论文最大偏差 0.229。
- manual_400us：2D S(60)=0.9668，对论文最大偏差 0.553。
- manual_600us：2D S(60)=0.9980，对论文最大偏差 0.425。
- ml_400us：2D S(60)=0.9590，对论文最大偏差 0.296。

- 程序化图片检查：8/8 张全部通过。
- 未进行人工视觉审阅；程序化检查不能证明标签无重叠或版面美观。
