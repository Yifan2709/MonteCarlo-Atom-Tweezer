"""level0_static_trap: Level 0 数值框架——一维径向静态高斯光偶极阱自洽检查。

本包只建立一个静态、一维、经典、单原子的径向光镊模型，用作后续
SLM-AOD 转移和运输模拟之前的数值自洽检查。
"""

# 在导入包时强制使用非交互后端，避免绘图时弹出窗口或因缺少显示而失败。
import matplotlib

matplotlib.use("Agg")

__version__ = "0.1.0"
