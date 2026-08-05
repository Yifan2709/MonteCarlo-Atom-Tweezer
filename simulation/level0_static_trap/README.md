# Level 0 — 一维径向静态高斯光偶极阱自洽检查

本目录实现 **Level 0** 数值框架：在光束焦平面的一个径向方向上建立静态、一维、
经典、单原子的高斯光偶极势，并完成数值自洽检查。它为后续 SLM-AOD 转移与运输
模拟提供经过验证的势函数、积分器、频率与能量分析工具。

> 本任务受 arXiv:2403.12021v4 实验参数启发，但**不是论文复现**。

## 明确的边界（必须读懂）

- `x` 是高斯光束**焦平面内的径向坐标**，**不是**光轴方向；
- **不**模拟论文报告的轴向 5.64 kHz 运动；
- **不**模拟 SLM-AOD 交接、势阱深度随时间变化、长距离移动、AOD 柱面透镜效应、
  二维运动、噪声、量子内态、加热、损失或多原子；
- 在本 Level 0 中，SLM 与 AOD **调用同一个高斯势函数**，二者只通过配置参数
  与参数来源区分；
- 数值曲率与解析曲率的比较是**代码一致性检查**，不是独立实验验证；
- 默认 180 μK SLM 与 280 μK AOD 是**两个独立静态势阱**的 Level 0 比较，
  不代表论文 Fig.6 转移时任意时刻的实际重叠势（论文该处 SLM 深度约 140 μK）；
- `velocity_verlet` **只**对自治势 `F(x)` 保证辛性质；若未来扩展到 `F(x,t)`，
  接口、能量解释与辛结构都需要重新审视，**不得**仅凭“传入函数”就宣称当前
  积分器已正确支持任意时变势。

## 安装

```bash
# 在本目录下（simulation/level0_static_trap/）执行
python -m pip install -e .
python -m pytest
python -m level0_static_trap.cli \
  --config configs/default.yaml \
  --trap both \
  --output-dir outputs/level0_static_trap/demo
```

运行依赖：NumPy、SciPy、Matplotlib、PyYAML、Pillow；开发依赖：pytest。

## 模块

| 模块 | 职责 |
|------|------|
| `constants.py` | k_B、m_u 及单位换算（唯一来源） |
| `config.py` | 读取/验证 YAML，转 SI，输出路径解析 |
| `potentials.py` | 高斯势、力、解析曲率、解析频率 |
| `analysis.py` | 平衡位置、数值曲率、轨迹频率、非线性周期、汇总 |
| `integrators.py` | 静态 `F(x)` 的 velocity-Verlet |
| `simulation.py` | 单势阱流程编排 |
| `visualization.py` | 势/力/轨迹/能量/频率比较图（英文标签，Agg 后端） |
| `io_utils.py` | CSV/JSON/YAML/Markdown 写出，图像程序化检查 |
| `cli.py` | 命令行入口 |

## 误差分类（关键）

- **数值离散化误差**：数值曲率频率与解析小振幅频率的相对差异
  `|f_num_curv − f_analytic| / f_analytic`（目标 `< 1e-4`）。
- **有限振幅非简谐频移**：默认 `0.10w` 振幅下，小振幅频率与非线性积分频率
  的相对差异（约 1–3%，物理软化效应，**不是**积分器误差）。
- **轨迹积分误差**：轨迹估计频率与同能量非线性积分频率的相对差异
  （严格数值容差）。

详细输出说明见 `summary.md` 与 `metrics.json`。
