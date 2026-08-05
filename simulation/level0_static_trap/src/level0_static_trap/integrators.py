"""velocity-Verlet 积分器（针对静态、位置依赖力 ``F(x)``）。

注意（README 同步说明）：本积分器仅对**自治势** ``F(x)`` 保证辛性质。若未来扩展到
时变力 ``F(x, t)``，接口、能量解释与辛结构都需要重新审视——不得仅凭“传入函数”
就宣称当前实现已正确支持任意时变势。
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

ForceFunction = Callable[[float | np.ndarray], float | np.ndarray]


class IntegrationError(ValueError):
    """积分器输入或时间网格不合法。"""


def _check_uniform_times(times: np.ndarray) -> None:
    if times.ndim != 1 or times.size < 2:
        raise IntegrationError(f"times 必须是至少 2 个点的一维数组，实际 shape={times.shape}")
    dt = float(times[1] - times[0])
    if dt <= 0.0:
        raise IntegrationError(f"times 必须严格递增，实际首步 dt={dt:g}")
    diffs = np.diff(times)
    if not np.all(np.isfinite(diffs)):
        raise IntegrationError("times 含非有限值")
    rel = np.max(np.abs(diffs - dt)) / dt
    # 等间隔检查：相对偏差远小于典型双精度
    if rel > 1.0e-10:
        raise IntegrationError(
            f"times 必须等间隔，最大相对偏差 {rel:.3e}（阈值 1e-10）"
        )


def velocity_verlet(
    force_function: ForceFunction,
    mass: float,
    x_initial: float,
    v_initial: float,
    times: np.ndarray,
):
    """对静态 ``F(x)`` 使用 velocity-Verlet 推进轨迹。

    返回 ``(times, x, v, a)``，全部为一维 NumPy 数组（与 ``times`` 同形状）。

    采用 kick-drift-kick 形式，每步恰一次力求值（首步除外）。

    参数
    ----
    force_function : 给定位置返回力（标量或数组 -> 标量）
    mass : 质量 (kg)，必须为正
    x_initial, v_initial : 初始位置 (m) 与速度 (m/s)
    times : 等间隔、严格递增的时间网格 (s)
    """
    if mass <= 0.0:
        raise IntegrationError(f"mass 必须为正，实际值：{mass:g} kg")
    if not math.isfinite(x_initial) or not math.isfinite(v_initial):
        raise IntegrationError("初始位置或速度非有限")
    times = np.asarray(times, dtype=float)
    _check_uniform_times(times)

    n = times.size
    dt = float(times[1] - times[0])
    x = np.empty(n, dtype=float)
    v = np.empty(n, dtype=float)
    a = np.empty(n, dtype=float)

    x[0] = float(x_initial)
    v[0] = float(v_initial)
    a[0] = float(np.asarray(force_function(x[0])).item()) / mass

    for i in range(n - 1):
        # drift（半步位置更新）
        x_new = x[i] + v[i] * dt + 0.5 * a[i] * dt * dt
        a_new = float(np.asarray(force_function(x_new)).item()) / mass
        # kick（速度全步）
        v_new = v[i] + 0.5 * (a[i] + a_new) * dt
        x[i + 1] = x_new
        v[i + 1] = v_new
        a[i + 1] = a_new

    return times, x, v, a
