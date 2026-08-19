"""三类运输轨迹族：adiabatic-sine、piecewise constant-jerk、minimum-jerk。"""
from __future__ import annotations

import numpy as np


class AdiabaticSine:
    """q(s) = s - sin(2πs)/(2π)，论文 Methods 的绝热正弦轨迹。"""
    name = "adiabatic_sine"

    def q(self, s):
        return s - np.sin(2 * np.pi * s) / (2 * np.pi)

    def q_dot(self, s):
        return 1.0 - np.cos(2 * np.pi * s)

    def q_ddot(self, s):
        return 2 * np.pi * np.sin(2 * np.pi * s)

    def q_jerk(self, s):
        return (2 * np.pi) ** 2 * np.cos(2 * np.pi * s)

    q_dot_max = 2.0
    q_ddot_max = 2 * np.pi


class MinimumJerk:
    """q(s) = 10s^3-15s^4+6s^5，常用平滑控制轨迹（非论文 constant-jerk）。"""
    name = "minimum_jerk"

    def q(self, s):
        return s ** 3 * (10.0 - 15.0 * s + 6.0 * s ** 2)

    def q_dot(self, s):
        return 30.0 * s ** 2 - 60.0 * s ** 3 + 30.0 * s ** 4

    def q_ddot(self, s):
        return 60.0 * s - 180.0 * s ** 2 + 120.0 * s ** 3

    def q_jerk(self, s):
        return 60.0 - 360.0 * s + 360.0 * s ** 2

    q_dot_max = 1.875
    q_ddot_max = 5.7735


class ConstantJerk:
    """四等分段、jerk 符号 +J,-J,-J,+J 的对称 S 曲线，J = 32L/T^3。

    分段解析积分（段内分数 u，s 导数）：段末常数 q:[0,1/12,1/2,11/12]、
    v:[0,1,2,1]、a:[0,8,0,-8]，端点速度/加速度为零，q(1)=1。
    """
    name = "constant_jerk"

    _Q0 = np.array([0.0, 1 / 12, 1 / 2, 11 / 12])

    def _segment(self, s):
        s = np.clip(np.asarray(s, dtype=float), 0.0, 1.0)
        seg = np.clip((s * 4).astype(int), 0, 3)
        u = np.clip(s * 4.0 - seg, 0.0, 1.0)
        return seg, u

    def q(self, s):
        seg, u = self._segment(s)
        return (self._Q0[seg] + np.select(
            [seg == 0, seg == 1, seg == 2, seg == 3],
            [u ** 3 / 12.0,
             u / 4.0 + u ** 2 / 4.0 - u ** 3 / 12.0,
             u / 2.0 - u ** 3 / 12.0,
             u / 4.0 - u ** 2 / 4.0 + u ** 3 / 12.0]))

    def q_dot(self, s):
        seg, u = self._segment(s)
        return np.select(
            [seg == 0, seg == 1, seg == 2, seg == 3],
            [u ** 2,
             1.0 + 2.0 * u - u ** 2,
             2.0 - u ** 2,
             1.0 - 2.0 * u + u ** 2])

    def q_ddot(self, s):
        seg, u = self._segment(s)
        return np.select(
            [seg == 0, seg == 1, seg == 2, seg == 3],
            [2.0 * u, 2.0 - 2.0 * u, -2.0 * u, -2.0 + 2.0 * u]) * 4.0

    def q_jerk(self, s):
        seg, _ = self._segment(s)
        return np.array([32.0, -32.0, -32.0, 32.0])[seg]

    q_dot_max = 2.0
    q_ddot_max = 8.0


def get_profile(name):
    """按名称返回轨迹实例。"""
    profiles = {"adiabatic_sine": AdiabaticSine, "minimum_jerk": MinimumJerk,
                "constant_jerk": ConstantJerk}
    if name not in profiles:
        raise KeyError(f"未知轨迹族 {name}")
    return profiles[name]()


def geometry_centers(q_value, geometry, distance_m):
    """由归一化路径分数返回 (c1, c2)；straight 单轴，diagonal 各 L/√2。"""
    if geometry == "straight":
        return distance_m * q_value, 0.0
    if geometry == "diagonal":
        return distance_m * q_value / np.sqrt(2.0), distance_m * q_value / np.sqrt(2.0)
    raise KeyError(f"未知几何 {geometry}")


def geometry_velocity_factors(geometry):
    """单轴速度相对路径速度 q̇·L 的系数。"""
    if geometry == "straight":
        return 1.0, 0.0
    return 1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)
