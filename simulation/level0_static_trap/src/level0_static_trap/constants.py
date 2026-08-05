"""物理常数与单位换算（全包唯一来源）。

其他模块不得重复定义 k_B 或 m_u；一律从这里导入。
所有公开函数/常量附简短中文 docstring。
"""

from __future__ import annotations

import math

import numpy as np

#: 玻尔兹曼常数 (J/K)，CODATA 推荐值量级。
BOLTZMANN_J_PER_K: float = 1.380649e-23  # 精确值（2019 SI 重定义）

#: 原子质量单位 (kg)。
ATOMIC_MASS_UNIT_KG: float = 1.66053906660e-27


def uk_to_j(temperature_uk):
    """把微开尔文温度换算为能量 (J)：U = k_B * T_μK * 1e-6。支持标量与数组。"""
    return BOLTZMANN_J_PER_K * 1.0e-6 * np.asarray(temperature_uk, dtype=float)


def um_to_m(length_um):
    """把微米换算为米。支持标量与数组。"""
    return 1.0e-6 * np.asarray(length_um, dtype=float)


def u_to_kg(mass_u: float) -> float:
    """把原子质量单位 (u) 换算为千克。"""
    return ATOMIC_MASS_UNIT_KG * float(mass_u)


def j_to_uk(energy_j):
    """把能量 (J) 换算回微开尔文温度。支持标量与数组。"""
    return np.asarray(energy_j, dtype=float) / (BOLTZMANN_J_PER_K * 1.0e-6)


def is_finite_number(value) -> bool:
    """判断标量是否为有限实数。"""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(f)
