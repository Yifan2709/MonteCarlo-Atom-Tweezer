"""一维径向高斯光偶极势、力与解析频率。

SLM 与 AOD 共享本模块的同一套实现，仅通过配置参数区分。所有公开函数同时支持
标量与 NumPy 数组输入，并附简短中文 docstring。
"""

from __future__ import annotations

import math

import numpy as np


def gaussian_potential(
    x, depth: float, waist: float, center: float = 0.0
):
    """高斯光偶极势 ``U(x) = -U0 * exp(-2 (x-xc)^2 / w^2)``，单位 J。

    参数
    ----
    x : 标量或数组，位置 (m)
    depth : U0 > 0，势阱深度 (J)
    waist : w > 0，1/e^2 强度半径 (m)
    center : xc，势阱中心 (m)
    """
    if depth <= 0.0:
        raise ValueError(f"depth 必须为正，实际值：{depth:g} J")
    if waist <= 0.0:
        raise ValueError(f"waist 必须为正，实际值：{waist:g} m")
    xa = np.asanyarray(x, dtype=float)
    z2 = ((xa - center) / waist) ** 2
    return -depth * np.exp(-2.0 * z2)


def gaussian_force(
    x, depth: float, waist: float, center: float = 0.0
):
    """高斯势的力 ``F(x) = -dU/dx = -4 U0 (x-xc)/w^2 * exp(-2 (x-xc)^2/w^2)`` (N)。

    正方向定义：力指向势阱中心。
    """
    if depth <= 0.0:
        raise ValueError(f"depth 必须为正，实际值：{depth:g} J")
    if waist <= 0.0:
        raise ValueError(f"waist 必须为正，实际值：{waist:g} m")
    xa = np.asanyarray(x, dtype=float)
    dx = xa - center
    z2 = (dx / waist) ** 2
    return -4.0 * depth * dx / (waist * waist) * np.exp(-2.0 * z2)


def gaussian_curvature_analytic(depth: float, waist: float) -> float:
    """势阱中心解析曲率 ``kappa0 = d^2U/dx^2|_xc = 4 U0 / w^2`` (N/m)。"""
    if depth <= 0.0:
        raise ValueError(f"depth 必须为正，实际值：{depth:g} J")
    if waist <= 0.0:
        raise ValueError(f"waist 必须为正，实际值：{waist:g} m")
    return 4.0 * depth / (waist * waist)


def trap_frequency_analytic(depth: float, waist: float, mass: float) -> float:
    """解析小振幅普通频率 ``f0 = sqrt(kappa0/m) / (2 pi)`` (Hz)。"""
    if mass <= 0.0:
        raise ValueError(f"mass 必须为正，实际值：{mass:g} kg")
    kappa = gaussian_curvature_analytic(depth, waist)
    omega = math.sqrt(kappa / mass)
    return omega / (2.0 * math.pi)
