"""理想瞬时 dynamical-decoupling 脉冲序列与 toggling function。

在仅有理想纯退相位与瞬时 π 脉冲的模型中，X/Y 轴顺序对 toggling function
等价（每个脉冲都使退相位项变号）；不得据此声称已模拟 XY4 对脉冲误差或
高阶噪声的完整抑制。
"""
from __future__ import annotations

import numpy as np


def toggling_function(t, pulse_times):
    """y(t)=(-1)^{N(t)}，N(t) 为 t 之前（严格小于）的脉冲数。"""
    t = np.asarray(t, dtype=float)
    pulses = np.sort(np.asarray(pulse_times, dtype=float))
    n_before = np.searchsorted(pulses, t, side="left")
    return (-1.0) ** n_before.astype(float)


def pulse_indices_in_steps(pulse_times, times):
    """定位每个脉冲落入的时间步区间，返回 (step_index, fraction) 数组。

    脉冲必须严格落在某一步内部（不含端点），否则视为放置错误。
    """
    out = []
    for p in np.sort(np.asarray(pulse_times, dtype=float)):
        idx = np.searchsorted(times, p, side="left") - 1
        if idx < 0 or idx >= len(times) - 1:
            raise ValueError(f"脉冲 {p} 落在时间网格之外")
        t_i, t_ip1 = times[idx], times[idx + 1]
        if not (t_i < p < t_ip1):
            raise ValueError(f"脉冲 {p} 恰好落在时间样本 {t_i} 上")
        out.append((int(idx), float((p - t_i) / (t_ip1 - t_i))))
    return out


def validate_pulse_set(pulse_times, dt, total_duration_s, boundaries_s,
                       margin_factor=0.25):
    """检查脉冲时刻：网格内部、远离分段边界与时间样本。"""
    issues = []
    pulses = np.sort(np.asarray(pulse_times, dtype=float))
    for p in pulses:
        if p <= 0 or p >= total_duration_s:
            issues.append(f"脉冲 {p:.6e} s 超出 (0, T) 内部范围")
            continue
        if abs(p / dt - round(p / dt)) < 1e-9:
            issues.append(f"脉冲 {p:.6e} s 恰好落在时间样本上")
        if len(boundaries_s) and np.min(np.abs(np.asarray(boundaries_s) - p)) \
                < margin_factor * dt:
            issues.append(f"脉冲 {p:.6e} s 过于接近分段边界")
    duplicates = [p for p, c in zip(*np.unique(pulses, return_counts=True))
                  if c > 1]
    if duplicates:
        issues.append(f"存在重复脉冲时刻: {duplicates}")
    return issues


def symmetric_xy4_fraction():
    """对称 XY4 的脉冲分数（段内归一化时间 1/8,3/8,5/8,7/8）。"""
    return (1.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0, 7.0 / 8.0)
