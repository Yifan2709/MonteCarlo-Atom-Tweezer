# Level 5 论文对照与全流程

最新报告：[Level5_MonteCarlo_论文对照与全流程_2026-09-15_参考数据修正版.pptx](output/Level5_MonteCarlo_论文对照与全流程_2026-09-15_参考数据修正版.pptx)

四曲线图：[PNG](output/Level5_修正参考曲线对照.png)、[SVG](output/Level5_修正参考曲线对照.svg)。数值更正：[论文参考数据更正与重新对照](output/论文参考数据更正与重新对照.md)。

本报告是独立的 24 页 PPT，以论文 Fig. 4、Fig. 5、Fig. 6 和 Extended Data Fig. 10 的实验结果为对照。先讨论过程覆盖、定量匹配和差距，再说明 Monte Carlo 六个阶段的输入、输出与跨轮状态传递，并单独解释三项平均加热过程。

数据来自 2026-09-14 已归档的连续仿真，本次没有重新运行模拟或修改物理模型。2026-09-15 更正旧像素追踪参考，直接从 PDF 矢量对象分别读取拟合路径、37 个散点中心及误差条，并重新评估全部 184 组结果。保留原比较区间和 ±0.02 比较容差，四条冻结参数曲线仍均未通过。散点不连成稠密折线，拟合线不作额外平滑或重拟合。旧参考、旧模拟输出和旧 PPT 均保留用于追溯。

当前累计条件相干指标与论文单次瞬时 IRB 保真度分别呈现。旧 IRB 输出的无效状态没有因本次报告更新而改变。

## 来源

- 论文原文：仓库根目录 `2403.12021v4.pdf`，v4。
- 更正后数值：`output/Level5_vector_reference_v2.evidence.json`。
- 全部 184 组重新对照：`output/Level5_vector_reference_v2.recomparison.json`。
- PDF 矢量参考：`../../simulation/level2_joint_transfer/reference/fig6d_vector_v2/`。
- 历史数值：`../stage_review/output/AOD_SLM_continuous_level2c_345_2026-09-14_v10.evidence.json`。
- 原始连续结果：`../../simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914/`。
- 物理流程：`../../simulation/level2_joint_transfer/src/continuous_transfer/`。
- 历史 RB 结果状态：`../../simulation/level2_joint_transfer/outputs/level5_finite_pulse_irb/demo/RESULT_STATUS.md`。

每页备注包含相关来源。表格、数据曲线与流程图均为原生 PowerPoint 对象，论文原图为局部图像。

## 重新生成

依次使用工作区 Python 运行 `../../simulation/level2_joint_transfer/tools/extract_fig6d_vectors.py`、`recompare_vector_reference.py`、`prepare_paper_figures.py`。可运行 `plot_vector_reference.py` 生成带原图误差条的科学图和源图叠加核验图。

使用带有 `@oai/artifact-tool` 的工作区 Node 运行 `build_level5_review.mjs`。此脚本调用 `add_paper_error_bars.py` 为原生散点系列补齐每个点的非对称误差条，所有数值来自 PDF 中实际误差条的上下端点。

最终文件不允许覆盖。再次生成时设置 `REPORT_PPTX` 为新的绝对输出路径，并用 `REPORT_REV` 指定新的内部构建编号。脚本会校验并渲染全部幻灯片，布局复核文件保存在 `.build` 中。
