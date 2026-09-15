# Level 2C 对标参考表：论文 Fig. 6 / Methods "Atom transfer between SLM and AOD tweezers"

> **2026-09-15 更正：本页第 3–4 节记录的是旧像素参考，现已停用。** 旧程序混合追踪拟合线、菱形标记及遮挡边缘，产生伪起伏。“曲线单调性可靠”这一旧声明不成立。旧 CSV 保留用于结果追溯，不再作为默认验收参考。
>
> 新参考见 [fig6d_vector_v2/README.md](fig6d_vector_v2/README.md)。直接读取 PDF 中独立的拟合路径、37 个散点中心和原图误差条；默认配置使用 `published_fit_reference.csv`。±0.02 仅保留为与既有验收可比的容差，不再称为提取误差或论文误差条。

> 来源：arXiv 2403.12021v4（*A tweezer array with 6100 highly coherent atomic qubits*）。
> 本文件是 Level 2C 的固定对标数据，一旦建立不得为迁就模拟结果而修改。
> 生成日期：2026-08-31。

## 1. 实验参数（Methods 原文）

| 参数 | 数值 | 出处 |
|---|---|---|
| SLM 阱深（转移时） | ~k_B × 0.14 mK（从全深绝热降下，无损失） | Methods: "ramped down to ∼kB × 0.14 mK" |
| AOD 阱深 | k_B × 0.28 mK | Methods / Fig. 6 |
| 分裂移动距离 | 2.4 µm | Fig. 6  caption: "picked up and moved away by 2.4 µm" |
| pick-up 后等待段 | 100 µs（期间 SLM 关闭以清除残留原子，再重开） | Fig. 6 caption / Methods |
| 手工轨迹结构 | 前 48% 时长：AOD 深度二次上升；后 52%：constant-jerk 移动 2.4 µm；顺序不重叠 | Methods |
| ML 轨迹结构 | 深度与位置各 14 个控制点、三次插值、ramp 与 move 同时；单程固定 400 µs；优化目标为 60 次连续单程转移的存活率 | Methods |
| 单次转移时长组 | 手工 0.2 / 0.4 / 0.6 ms；ML 0.4 ms | Fig. 6d 顶图 |
| 阵列规模 | 195 个 AOD 阱与 SLM 重合（15×14 AOD 网格，5 倍 SLM 间距）；504 µm × 468 µm | Methods |
| SLM 激光 | 1061 nm（AOD 光路 1055 nm，见 Ext. Data Fig. 2） | Methods / Ext. Data Fig. 2 |
| 对准方式 | 先相机对 PSF，再按转移存活率优化；残余误差即位点间差异来源 | Methods |
| 存活率拟合公式 | S_n = p0 · p^n · (1 − e^(−1/(bn)))，不确定度用参数 bootstrap | Methods |

## 2. SI 文字锚点

- **200 µs 单程转移存活率 > 98%**（SI 重排时间预算一节：
  "a single one-way transfer can be done in 200 µs with a survival probability of more than 98%"）。

## 3. Fig. 6d 顶图数字化曲线

见 `fig6d_survival_digitized.csv`（由 `tools/digitize_fig6d.py` 从 PDF 第 7 页自动生成，
校验叠加图 `fig6d_digitize_check.png`）。覆盖范围与精度声明：

- n = 5..60 整数网格；0.2 ms 曲线只到 n = 31（其后曲线超出面板下界，无数据点）；
- n < 5 处四条曲线互相遮挡无法区分，未数字化；目视均为 0.97~1.00；
- 整体精度 ±0.01~0.02（图像目视级数字化），曲线单调性与走势可靠，单点数值不作真值。

关键读数（供快速引用）：

| n | ML 0.4 ms | 手工 0.6 ms | 手工 0.4 ms | 手工 0.2 ms |
|---|---|---|---|---|
| 5  | 0.990 | 0.971 | 0.968 | 0.925 |
| 20 | 0.940 | 0.889 | 0.824 | 0.524 |
| 40 | 0.815 | 0.726 | 0.582 | — |
| 60 | 0.663 | 0.573 | 0.414 | — |

## 4. 明确的验收标准

Level 2C 的对比对象是**存活率曲线本身**；99.81(3)% 是 IRB 单程保真度，不在范围内。

- 主要判据：对每条轨迹（手工 200/400/600 µs 与 ML 400 µs），模拟 survival(n) 曲线
  在 n = 5..60 范围内落入数字化曲线的 ±0.02 容差带（≈数字化精度）内，
  或拟合出的单次转移存活率 p 与论文曲线拟合的 p 在 bootstrap 置信区间内一致；
- 若达不到：按 prompt 阶段 5 第 14 步**诚实报告缺口并归因**，
  不得调噪声参数强行凑数；
- 次要锚点：保守模型（无噪声）下 200/400/600 µs 单程存活率应 ≈100%，
  与 SI "200 µs 单程 > 98%" 不矛盾（实验含未建模损失）。
