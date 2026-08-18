"""一维焦平面径向高斯光镊势、力和解析小振幅量。"""

from __future__ import annotations

import numpy as np


def gaussian_potential(x, depth: float, waist: float, center: float = 0.0):
    """计算高斯光镊势能；输入采用 SI 单位并支持标量或数组。"""

    coordinate = np.asarray(x)
    return -depth * np.exp(-2.0 * ((coordinate - center) / waist) ** 2)


def gaussian_force(x, depth: float, waist: float, center: float = 0.0):
    """计算指向势阱中心的高斯光镊保守力。"""

    coordinate = np.asarray(x)
    offset = coordinate - center
    return -(4.0 * depth * offset / waist**2) * np.exp(-2.0 * (offset / waist) ** 2)


def gaussian_curvature_analytic(depth: float, waist: float) -> float:
    """返回势阱中心的解析曲率。"""

    return 4.0 * depth / waist**2


def trap_frequency_analytic(depth: float, waist: float, mass: float) -> tuple[float, float]:
    """返回解析小振幅角频率和普通频率。"""

    omega = float(np.sqrt(gaussian_curvature_analytic(depth, waist) / mass))
    return omega, omega / (2.0 * np.pi)

