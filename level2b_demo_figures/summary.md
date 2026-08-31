# Level 2B AOD 移动速度边界与变速剖面验证

在 Level 0/1/2 一维经典模型内，把移动速度作为显式扫描变量：实验 1 纯移动跟随上限、实验 2 端到端速度上限（三结构 × 三温度）、实验 3 恒速 vs 变速双口径配对比较。

- 预检查：全部通过（456 s）。
- 扫描：94 点完成、25 点曲线级早停跳过、0 点失败（耗时 642 s）。
- frozen_shape_scaled @ 5 μK：v95 ≈ 54.9 mm/s，边界已覆盖。
- frozen_shape_scaled @ 10 μK：v95 ≈ 47.9 mm/s，边界已覆盖。
- frozen_shape_scaled @ 20 μK：v95 ≈ 41.2 mm/s，边界已覆盖。
- sequential_scaled @ 5 μK：v95 ≈ 56.4 mm/s，边界已覆盖。
- sequential_scaled @ 10 μK：v95 ≈ 54.7 mm/s，边界已覆盖。
- sequential_scaled @ 20 μK：v95 ≈ 53.2 mm/s，边界已覆盖。
- trapezoid_frozen_ramp @ 5 μK：v95 ≈ 55.3 mm/s，边界已覆盖。
- trapezoid_frozen_ramp @ 10 μK：v95 ≈ 55.3 mm/s，边界已覆盖。
- trapezoid_frozen_ramp @ 20 μK：v95 ≈ 54.5 mm/s，边界已覆盖。
- 实验 1 标度律：ε_med ∝ v^1.33（95% CI [1.08, 1.82]，R²=0.826）。
- 实验 3 参考曲线：frozen_shape_scaled @ 20 μK（bootstrap_interpolation，v_ref = 41.2 mm/s）。
- 独立复核：30/30 个 (口径, 速度, 剖面) 组合的 95% CI 与扫描池重叠。
- 图片程序化检查：全部通过；未进行人工视觉审阅（visual_review_performed: false）。

## 边界

本 Level 仍是一维经典模型且不含 AOD 器件带宽/RF 扫描速率限制：测得的是该模型下的动力学速度上限，不是器件规格。
