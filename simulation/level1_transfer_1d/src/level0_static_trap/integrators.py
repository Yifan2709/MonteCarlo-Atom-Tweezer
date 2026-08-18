"""静态位置依赖保守力的 velocity-Verlet 积分器。"""

from __future__ import annotations

from typing import Callable

import numpy as np


def velocity_verlet(
    force_function: Callable[[float], float], mass: float, x_initial: float, v_initial: float, times
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """对自治力 F(x) 执行 velocity-Verlet 并返回 t、x、v、a 数组。"""

    time = np.asarray(times, dtype=float)
    if time.ndim != 1 or time.size < 2 or not np.all(np.isfinite(time)):
        raise ValueError("时间数组必须是一维、有限且至少包含两个点")
    differences = np.diff(time)
    if np.any(differences <= 0.0):
        raise ValueError("时间数组必须严格递增")
    step = float(differences[0])
    if not np.allclose(differences, step, rtol=1.0e-11, atol=max(1.0e-18, abs(step) * 1.0e-13)):
        raise ValueError("时间数组必须等间隔")
    if mass <= 0.0 or not np.isfinite(mass):
        raise ValueError("质量必须是有限正数")

    positions = np.empty_like(time)
    velocities = np.empty_like(time)
    accelerations = np.empty_like(time)
    positions[0] = x_initial
    velocities[0] = v_initial
    accelerations[0] = float(force_function(x_initial)) / mass
    for index in range(time.size - 1):
        positions[index + 1] = positions[index] + velocities[index] * step + 0.5 * accelerations[index] * step**2
        accelerations[index + 1] = float(force_function(positions[index + 1])) / mass
        velocities[index + 1] = velocities[index] + 0.5 * (accelerations[index] + accelerations[index + 1]) * step
    return time, positions, velocities, accelerations

