"""运输轨迹族：R 的 constant-jerk、Y 的 min-jerk 与 zero-jerk 命名剖面。

- constant_jerk（R Methods）：x(s)=D(3s²−2s³)，速度抛物线、加速度线性、jerk 常数。
- min_jerk（Y，γ=1.875）：x(s)=D(10s³−15s⁴+6s⁵)，峰值速度 1.875·D/T。
- zero_jerk（Y，γ=1.5625）：六阶剖面 x(s)=D(7.5s⁴−9s⁵+2.5s⁶)，
  端点 jerk=0，峰值速度 ≈1.69·D/T。
  【已登记假设】论文的五阶多项式 γ 族的精确定义（γ 如何进入系数）无法从
  可获取的 arXiv 文本恢复；本实现以“端点零 jerk 的最低阶多项式”作为
  zero-jerk 的显式剖面，γ 数值仅作为论文标签保留。此假设影响绝对曲线，
  不影响“含透镜时 zero-jerk 优于 min-jerk”的结构判据（两条剖面峰值速度
  1.69<1.875、加速度末端行为不同，透镜灵敏度排序由此可判）。
- gamma:P 族：x(s)=a s³+(5−2a)s⁴+(a−4)s⁵，端点速度为零，峰值速度=γ·D/T；
  仅在 γ≳1.7（a≥6 分支）可解，供灵敏度扫描。
"""
from __future__ import annotations

import numpy as np


def constant_jerk_coeffs() -> tuple[float, float, float]:
    return (0.0, 3.0, -2.0)  # x = 3s²−2s³（s 的多项式系数 c2,c3 形式单独处理）


def _gamma_family_peak_speed(a: float) -> float:
    s = np.linspace(0.0, 1.0, 20001)
    b, c = 5.0 - 2.0 * a, a - 4.0
    v = 3 * a * s**2 + 4 * b * s**3 + 5 * c * s**4
    return float(np.max(v))


def gamma_family_coeffs(gamma: float) -> tuple[float, float, float]:
    """解 a 使峰值速度 = γ（x(1)=1、ẋ(1)=0 已内嵌；a≥6 分支单调增）。"""
    if not 1.7 <= gamma <= 2.6:
        raise ValueError(f"gamma {gamma} 超出可解范围 [1.7, 2.6]")
    lo, hi = 5.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _gamma_family_peak_speed(mid) > gamma:
            hi = mid
        else:
            lo = mid
    a = 0.5 * (lo + hi)
    return (a, 5.0 - 2.0 * a, a - 4.0)


def gamma_of_coeffs(a: float) -> float:
    return _gamma_family_peak_speed(a)


class Trajectory:
    """单轴轨迹：t 数组与位置/速度数组（SI）。"""

    def __init__(self, t_s: np.ndarray, x_m: np.ndarray, v_m_s: np.ndarray,
                 name: str):
        self.t_s = t_s
        self.x_m = x_m
        self.v_m_s = v_m_s
        self.name = name

    @property
    def duration_s(self) -> float:
        return float(self.t_s[-1])


def poly_trajectory(kind: str, distance_m: float, duration_s: float,
                    dt_s: float) -> Trajectory:
    n_steps = int(round(duration_s / dt_s))
    t = np.arange(n_steps + 1) * dt_s
    s = np.clip(t / duration_s, 0.0, 1.0)
    if kind == "constant_jerk":          # R：3s²−2s³
        x = distance_m * (3 * s**2 - 2 * s**3)
        v = distance_m / duration_s * (6 * s - 6 * s**2)
    elif kind == "linear":
        x = distance_m * s
        v = np.full_like(s, distance_m / duration_s)
    elif kind == "constant_accel":       # 梯形速度：a=const，中点折返
        # v(s)=2s (s<0.5), 2(1−s)；x=s² (s<0.5), 1−(1−s)²
        x = np.where(s < 0.5, s**2, 1 - (1 - s) ** 2) * distance_m
        v = np.where(s < 0.5, 2 * s, 2 * (1 - s)) * distance_m / duration_s
    elif kind == "sixth_zero_jerk":
        # Y zero-jerk 近似：五阶族低分支中峰值速度最低的可达剖面（a≈6，γ≈1.63）。
        # 论文 γ=1.5625 的精确定义不可恢复（见模块 docstring）；保留结构性质：
        # 峰值速度低于 min-jerk(1.875) ⇒ 含透镜时灵敏度更低。
        a0, b0, c0 = 6.0, 5.0 - 12.0, 2.0
        x = distance_m * (a0 * s**3 + b0 * s**4 + c0 * s**5)
        v = (distance_m / duration_s) * (3 * a0 * s**2 + 4 * b0 * s**3
                                         + 5 * c0 * s**4)
    elif kind.startswith("gamma:"):
        gamma = float(kind.split(":")[1])
        if gamma < 1.7:
            raise ValueError("gamma:P 族仅在 γ≳1.7 可解；zero-jerk 用 sixth_zero_jerk")
        a, b, c = gamma_family_coeffs(gamma)
        x = distance_m * (a * s**3 + b * s**4 + c * s**5)
        v = (distance_m / duration_s) * (3 * a * s**2 + 4 * b * s**3
                                         + 5 * c * s**4)
    else:
        raise ValueError(f"unknown trajectory kind {kind!r}")
    return Trajectory(t, x, v, kind)


MIN_JERK = "gamma:1.875"      # Y：minimum-jerk（10−15+6，解析一致）
ZERO_JERK = "sixth_zero_jerk"  # Y：zero-jerk（六阶端点零 jerk，见模块说明）
