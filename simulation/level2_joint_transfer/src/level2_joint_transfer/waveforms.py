"""Level 2 三类 AOD 波形：顺序基线、压缩基线、同步受约束族与单调样条精修。

所有波形统一返回 AOD 中心 c(t)、深度 D(t) 及一至三阶导数（导数用于功-能核算
与光滑度诊断）。时间输入为秒，支持标量与 NumPy 数组。

Level 2C 追加（pick-up 方向，论文 arXiv:2403.12021v4 Fig. 6b/c 结构）：
- PickupSequentialWaveform：先二次升深后 constant-jerk 移动（手工轨迹）；
  direction='dropoff' 时为其时间反演（merge & drop-off）；
- MLCubicWaveform：深度/位置各 N 个控制点的普通三次插值（ML 轨迹），
  允许非单调，配独立过冲校验。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator


def _out(value, original):
    """标量进标量出，数组进出数组。"""
    return float(value) if np.asarray(original).ndim == 0 else value


def smootherstep(z):
    """五次 smootherstep S(z)=10z^3-15z^4+6z^5，输入被裁剪到 [0,1]。"""
    z = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    return _out(10 * z**3 - 15 * z**4 + 6 * z**5, z)


def _window(z):
    """窗口函数：区间前 0，区间内 S(z)，区间后 1；z 已是归一化局部坐标。"""
    return smootherstep(np.clip(z, 0.0, 1.0))


class WaveformBase:
    """波形公共接口：中心、深度及导数，子类只需实现 _impl 系列。"""

    name = "base"
    duration_s = 0.0

    def center(self, t):
        """返回 AOD 中心 c(t)，单位 m；t 超出 [0,T] 时冻结端点值。"""
        raise NotImplementedError

    def depth(self, t):
        """返回 AOD 深度 D(t)，单位 J；t 超出 [0,T] 时冻结端点值。"""
        raise NotImplementedError

    def center_velocity(self, t):
        raise NotImplementedError

    def center_acceleration(self, t):
        raise NotImplementedError

    def center_jerk(self, t):
        raise NotImplementedError

    def depth_rate(self, t):
        """返回 dD/dt，单位 J/s。"""
        raise NotImplementedError

    def dense_table(self, num_points=2001):
        """在 [0,T] 上给出密集波形表及导数，用于 CSV 与光滑度统计。"""
        t = np.linspace(0.0, self.duration_s, num_points)
        return {
            "time_s": t,
            "aod_center_m": self.center(t),
            "aod_depth_J": self.depth(t),
            "center_velocity_m_per_s": self.center_velocity(t),
            "center_acceleration_m_per_s2": self.center_acceleration(t),
            "center_jerk_m_per_s3": self.center_jerk(t),
            "depth_rate_J_per_s": self.depth_rate(t),
        }

    def smoothness_metrics(self, num_points=2001):
        """返回最大速度、加速度、jerk 与深度变化率的无量纲与原始量。"""

        table = self.dense_table(num_points)
        span = self.initial_center_m
        tau = self.duration_s
        return {
            "max_abs_center_velocity_m_per_s": float(np.max(np.abs(table["center_velocity_m_per_s"]))),
            "max_abs_center_acceleration_m_per_s2": float(np.max(np.abs(table["center_acceleration_m_per_s2"]))),
            "max_abs_center_jerk_m_per_s3": float(np.max(np.abs(table["center_jerk_m_per_s3"]))),
            "max_abs_depth_rate_J_per_s": float(np.max(np.abs(table["depth_rate_J_per_s"]))),
            "dimensionless_accel": float(np.max(np.abs(table["center_acceleration_m_per_s2"])) * tau**2 / span),
            "dimensionless_jerk": float(np.max(np.abs(table["center_jerk_m_per_s3"])) * tau**3 / span),
        }

    def describe(self):
        """返回可 YAML 序列化的波形参数描述。"""
        return {"type": self.name, "duration_us": self.duration_s * 1e6}


class SequentialWaveform(WaveformBase):
    """严格复现 Level 1 的两段式顺序波形：先 smoothstep 移动，后 (1-s)^2 降深。"""

    name = "sequential"

    def __init__(self, duration_s, move_fraction, initial_center_m, initial_depth_j):
        self.duration_s = float(duration_s)
        self.move_fraction = float(move_fraction)
        self.initial_center_m = float(initial_center_m)
        self.initial_depth_j = float(initial_depth_j)
        self.move_duration_s = self.duration_s * self.move_fraction
        self.ramp_duration_s = self.duration_s - self.move_duration_s

    def _move_s(self, t):
        return np.clip(np.asarray(t, dtype=float) / self.move_duration_s, 0.0, 1.0)

    def _ramp_s(self, t):
        t = np.asarray(t, dtype=float)
        return np.clip((t - self.move_duration_s) / self.ramp_duration_s, 0.0, 1.0)

    def center(self, t):
        s = self._move_s(t)
        return _out(self.initial_center_m * (1.0 - (3 * s**2 - 2 * s**3)), t)

    def depth(self, t):
        s = self._ramp_s(t)
        return _out(self.initial_depth_j * (1.0 - s) ** 2, t)

    def center_velocity(self, t):
        s = self._move_s(t)
        inside = (np.asarray(t, dtype=float) > 0.0) & (np.asarray(t, dtype=float) < self.move_duration_s)
        derivative = -self.initial_center_m * (6 * s - 6 * s**2) / self.move_duration_s
        return _out(np.where(inside, derivative, 0.0), t)

    def center_acceleration(self, t):
        s = self._move_s(t)
        inside = (np.asarray(t, dtype=float) > 0.0) & (np.asarray(t, dtype=float) < self.move_duration_s)
        derivative = -self.initial_center_m * (6 - 12 * s) / self.move_duration_s**2
        return _out(np.where(inside, derivative, 0.0), t)

    def center_jerk(self, t):
        scalar = np.asarray(t, dtype=float)
        inside = (scalar > 0.0) & (scalar < self.move_duration_s)
        value = np.where(inside, 12.0 * self.initial_center_m / self.move_duration_s**3, 0.0)
        return _out(value, t)

    def depth_rate(self, t):
        s = self._ramp_s(t)
        t = np.asarray(t, dtype=float)
        inside = (t > self.move_duration_s) & (t < self.duration_s)
        derivative = -2.0 * self.initial_depth_j * (1.0 - s) / self.ramp_duration_s
        return _out(np.where(inside, derivative, 0.0), t)

    def describe(self):
        return {
            "type": self.name,
            "duration_us": self.duration_s * 1e6,
            "move_fraction": self.move_fraction,
            "initial_center_um": self.initial_center_m * 1e6,
            "initial_depth_uK": self.initial_depth_j / 1.380649e-29,
        }


class OverlappedWaveform(WaveformBase):
    """同步受约束波形族：移动与降深窗口可重叠，均为 smootherstep 窗口。"""

    name = "overlapped"

    def __init__(self, duration_s, move_start_fraction, move_end_fraction,
                 ramp_start_fraction, ramp_end_fraction, ramp_power,
                 initial_center_m, initial_depth_j):
        self.duration_s = float(duration_s)
        self.move_start_fraction = float(move_start_fraction)
        self.move_end_fraction = float(move_end_fraction)
        self.ramp_start_fraction = float(ramp_start_fraction)
        self.ramp_end_fraction = float(ramp_end_fraction)
        self.ramp_power = float(ramp_power)
        self.initial_center_m = float(initial_center_m)
        self.initial_depth_j = float(initial_depth_j)

    def _local(self, t, a, b):
        s = np.clip(np.asarray(t, dtype=float) / self.duration_s, 0.0, 1.0)
        return (s - a) / (b - a)

    def _window_terms(self, t, a, b):
        """返回窗口值 W、dW/ds 及区间内标志。"""
        z = np.clip(self._local(t, a, b), 0.0, 1.0)
        window = 10 * z**3 - 15 * z**4 + 6 * z**5
        dwindow_ds = (30 * z**2 - 60 * z**3 + 30 * z**4) / (b - a)
        inside = (np.asarray(t, dtype=float) > a * self.duration_s) & \
                 (np.asarray(t, dtype=float) < b * self.duration_s)
        return window, np.where(inside, dwindow_ds, 0.0), inside

    def center(self, t):
        window, _, _ = self._window_terms(t, self.move_start_fraction, self.move_end_fraction)
        return _out(self.initial_center_m * (1.0 - window), t)

    def depth(self, t):
        window, _, _ = self._window_terms(t, self.ramp_start_fraction, self.ramp_end_fraction)
        # p<1 的非法候选在窗口末端 0**负幂会产生警告；此类候选会被约束拒绝
        with np.errstate(divide="ignore", invalid="ignore"):
            value = self.initial_depth_j * (1.0 - window) ** self.ramp_power
        value = np.where(np.isfinite(value), value, 0.0)
        return _out(value, t)

    def center_velocity(self, t):
        _, dwindow_ds, _ = self._window_terms(t, self.move_start_fraction, self.move_end_fraction)
        return _out(-self.initial_center_m * dwindow_ds / self.duration_s, t)

    def center_acceleration(self, t):
        # d2W/ds2 = S''(z)/(b-a)^2, S''=60z-180z^2+120z^3，仅在区间内非零
        a, b = self.move_start_fraction, self.move_end_fraction
        z = np.clip(self._local(t, a, b), 0.0, 1.0)
        second = (60 * z - 180 * z**2 + 120 * z**3) / (b - a) ** 2
        inside = (np.asarray(t, dtype=float) > a * self.duration_s) & \
                 (np.asarray(t, dtype=float) < b * self.duration_s)
        return _out(-self.initial_center_m * np.where(inside, second, 0.0) / self.duration_s**2, t)

    def center_jerk(self, t):
        a, b = self.move_start_fraction, self.move_end_fraction
        z = np.clip(self._local(t, a, b), 0.0, 1.0)
        third = (60.0 - 360.0 * z + 360.0 * z**2) / (b - a) ** 3
        inside = (np.asarray(t, dtype=float) > a * self.duration_s) & \
                 (np.asarray(t, dtype=float) < b * self.duration_s)
        return _out(-self.initial_center_m * np.where(inside, third, 0.0) / self.duration_s**3, t)

    def depth_rate(self, t):
        window, dwindow_ds, _ = self._window_terms(t, self.ramp_start_fraction, self.ramp_end_fraction)
        complement = 1.0 - window
        # p<1 时端点导数本就奇异（此类候选会被约束检查拒绝）；此处置 0 仅为止住 NaN
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = complement ** (self.ramp_power - 1.0)
        factor = np.where(np.isfinite(factor), factor, 0.0)
        rate = -self.initial_depth_j * self.ramp_power * factor * dwindow_ds / self.duration_s
        return _out(rate, t)

    def params(self):
        return {
            "move_start_fraction": self.move_start_fraction,
            "move_end_fraction": self.move_end_fraction,
            "ramp_start_fraction": self.ramp_start_fraction,
            "ramp_end_fraction": self.ramp_end_fraction,
            "ramp_power": self.ramp_power,
        }

    def describe(self):
        data = self.describe_base()
        data.update(self.params())
        return data

    def describe_base(self):
        return {
            "type": self.name,
            "duration_us": self.duration_s * 1e6,
            "initial_center_um": self.initial_center_m * 1e6,
        }


class SplineWaveform(WaveformBase):
    """单调 PCHIP 样条波形：默认 8 个含端点控制点，端点值固定，无过冲。"""

    name = "spline"

    def __init__(self, duration_s, control_s, center_values_m, depth_values_j,
                 initial_center_m, initial_depth_j):
        self.duration_s = float(duration_s)
        self.control_s = np.asarray(control_s, dtype=float)
        self.center_values_m = np.asarray(center_values_m, dtype=float)
        self.depth_values_j = np.asarray(depth_values_j, dtype=float)
        self.initial_center_m = float(initial_center_m)
        self.initial_depth_j = float(initial_depth_j)
        self._center_interp = PchipInterpolator(self.control_s, self.center_values_m, extrapolate=False)
        self._depth_interp = PchipInterpolator(self.control_s, self.depth_values_j, extrapolate=False)

    def _s(self, t):
        return np.clip(np.asarray(t, dtype=float) / self.duration_s, 0.0, 1.0)

    def center(self, t):
        return _out(self._center_interp(self._s(t)), t)

    def depth(self, t):
        return _out(self._depth_interp(self._s(t)), t)

    def center_velocity(self, t):
        return _out(self._center_interp.derivative(1)(self._s(t)) / self.duration_s, t)

    def center_acceleration(self, t):
        return _out(self._center_interp.derivative(2)(self._s(t)) / self.duration_s**2, t)

    def center_jerk(self, t):
        return _out(self._center_interp.derivative(3)(self._s(t)) / self.duration_s**3, t)

    def depth_rate(self, t):
        return _out(self._depth_interp.derivative(1)(self._s(t)) / self.duration_s, t)

    def control_table(self):
        """返回控制点表，便于冻结与报告。"""
        return {
            "control_s": self.control_s.tolist(),
            "center_values_m": self.center_values_m.tolist(),
            "depth_values_j": self.depth_values_j.tolist(),
        }


def _constant_jerk_terms(u, span, tau):
    """constant-jerk（bang–coast–bang jerk）移动剖面：u∈[0,1] 局部时间。

    jerk 分段 +j0（前 1/4）、−j0（中 1/2）、+j0（后 1/4），j0=32·span/τ³；
    端点速度与加速度严格为零，总位移严格 span。返回 (x, v, a, j)。
    u 越界时冻结端点值（导数为零）。
    """
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    j0 = 32.0 * span / tau**3
    t = u * tau
    q = tau / 4.0
    x = np.empty_like(t); v = np.empty_like(t); a = np.empty_like(t); j = np.empty_like(t)
    seg_a = t < q
    seg_c = t > 3.0 * q
    seg_b = ~(seg_a | seg_c)
    # 段 A：[0, τ/4]，jerk = +j0
    ta = t[seg_a]
    x[seg_a] = j0 * ta**3 / 6.0
    v[seg_a] = j0 * ta**2 / 2.0
    a[seg_a] = j0 * ta
    j[seg_a] = j0
    # 段 B：[τ/4, 3τ/4]，jerk = −j0；段 A 末态 a1=j0q, v1=j0q²/2, x1=j0q³/6
    tb = t[seg_b] - q
    x[seg_b] = j0 * q**3 / 6.0 + (j0 * q**2 / 2.0) * tb + (j0 * q) * tb**2 / 2.0 - j0 * tb**3 / 6.0
    v[seg_b] = j0 * q**2 / 2.0 + j0 * q * tb - j0 * tb**2 / 2.0
    a[seg_b] = j0 * q - j0 * tb
    j[seg_b] = -j0
    # 段 C：[3τ/4, τ]，jerk = +j0；段 B 末态 v2=j0q²/2, x2=11j0q³/6, a2=−j0q
    tc = t[seg_c] - 3.0 * q
    x[seg_c] = 11.0 * j0 * q**3 / 6.0 + (j0 * q**2 / 2.0) * tc - (j0 * q) * tc**2 / 2.0 + j0 * tc**3 / 6.0
    v[seg_c] = j0 * q**2 / 2.0 - j0 * q * tc + j0 * tc**2 / 2.0
    a[seg_c] = -j0 * q + j0 * tc
    j[seg_c] = j0
    # 冻结端点之外
    frozen_end = np.asarray(u, dtype=float) >= 1.0
    if np.any(frozen_end):
        x[frozen_end] = span
        v[frozen_end] = 0.0
        a[frozen_end] = 0.0
        j[frozen_end] = 0.0
    return x, v, a, j


class PickupSequentialWaveform(WaveformBase):
    """论文手工 pick-up 波形：前 ramp_fraction 二次升深（中心停在 SLM 上），
    后 (1−ramp_fraction) 以 constant-jerk 移动 0→d，两段不重叠。

    direction='dropoff' 时为时间反演（merge & drop-off）：先 constant-jerk
    移动 d→0，再按 (1−s)² 二次降深到 0。端点速度/加速度均为零。
    """

    name = "pickup_sequential"

    def __init__(self, duration_s, ramp_fraction, final_center_m, final_depth_j,
                 direction="pickup"):
        if direction not in ("pickup", "dropoff"):
            raise ValueError("direction 必须是 'pickup' 或 'dropoff'")
        self.duration_s = float(duration_s)
        self.ramp_fraction = float(ramp_fraction)
        self.final_center_m = float(final_center_m)
        self.final_depth_j = float(final_depth_j)
        self.direction = direction
        self.ramp_duration_s = self.duration_s * self.ramp_fraction
        self.move_duration_s = self.duration_s - self.ramp_duration_s
        # 与 WaveformBase.smoothness_metrics 的归一化兼容（span/深度尺度）
        self.initial_center_m = self.final_center_m
        self.initial_depth_j = self.final_depth_j

    def _ramp_s(self, t):
        t = np.asarray(t, dtype=float)
        if self.direction == "pickup":
            return np.clip(t / self.ramp_duration_s, 0.0, 1.0)
        return np.clip((t - self.move_duration_s) / self.ramp_duration_s, 0.0, 1.0)

    def _move_u(self, t):
        t = np.asarray(t, dtype=float)
        if self.direction == "pickup":
            return np.clip((t - self.ramp_duration_s) / self.move_duration_s, 0.0, 1.0)
        return np.clip(t / self.move_duration_s, 0.0, 1.0)

    def center(self, t):
        u = self._move_u(t)
        span = self.final_center_m if self.direction == "pickup" else -self.final_center_m
        offset = 0.0 if self.direction == "pickup" else self.final_center_m
        x, _, _, _ = _constant_jerk_terms(u, span, self.move_duration_s)
        return _out(offset + x, t)

    def depth(self, t):
        s = self._ramp_s(t)
        value = self.final_depth_j * s**2 if self.direction == "pickup" \
            else self.final_depth_j * (1.0 - s) ** 2
        return _out(value, t)

    def center_velocity(self, t):
        u = self._move_u(t)
        span = self.final_center_m if self.direction == "pickup" else -self.final_center_m
        _, v, _, _ = _constant_jerk_terms(u, span, self.move_duration_s)
        return _out(v, t)

    def center_acceleration(self, t):
        u = self._move_u(t)
        span = self.final_center_m if self.direction == "pickup" else -self.final_center_m
        _, _, acc, _ = _constant_jerk_terms(u, span, self.move_duration_s)
        return _out(acc, t)

    def center_jerk(self, t):
        u = self._move_u(t)
        span = self.final_center_m if self.direction == "pickup" else -self.final_center_m
        _, _, _, jk = _constant_jerk_terms(u, span, self.move_duration_s)
        return _out(jk, t)

    def depth_rate(self, t):
        s = self._ramp_s(t)
        t = np.asarray(t, dtype=float)
        if self.direction == "pickup":
            inside = (t > 0.0) & (t < self.ramp_duration_s)
            derivative = 2.0 * self.final_depth_j * s / self.ramp_duration_s
        else:
            inside = (t > self.move_duration_s) & (t < self.duration_s)
            derivative = -2.0 * self.final_depth_j * (1.0 - s) / self.ramp_duration_s
        return _out(np.where(inside, derivative, 0.0), t)

    def describe(self):
        return {
            "type": self.name,
            "duration_us": self.duration_s * 1e6,
            "ramp_fraction": self.ramp_fraction,
            "final_center_um": self.final_center_m * 1e6,
            "final_depth_uK": self.final_depth_j / 1.380649e-29,
            "direction": self.direction,
        }


class MLCubicWaveform(WaveformBase):
    """论文 ML 轨迹：深度与位置各 N 个控制点的普通三次插值（not-a-knot）。

    与 SplineWaveform（PCHIP 单调）不同，允许非单调与过冲——这是 ML 优化的
    真实搜索空间；过冲幅度由 validate_pickup_waveform 按配置容差拦截。
    端点值由调用方固定在 (c0, D0)（pickup：0→d 与 0→D0；dropoff 反之）。
    """

    name = "ml_cubic"

    def __init__(self, duration_s, control_s, center_values_m, depth_values_j,
                 initial_center_m, initial_depth_j):
        self.duration_s = float(duration_s)
        self.control_s = np.asarray(control_s, dtype=float)
        self.center_values_m = np.asarray(center_values_m, dtype=float)
        self.depth_values_j = np.asarray(depth_values_j, dtype=float)
        self.initial_center_m = float(initial_center_m)
        self.initial_depth_j = float(initial_depth_j)
        if self.control_s[0] != 0.0 or self.control_s[-1] != 1.0:
            raise ValueError("ML 三次样条控制点必须覆盖 [0,1]")
        self._center_interp = CubicSpline(self.control_s, self.center_values_m)
        self._depth_interp = CubicSpline(self.control_s, self.depth_values_j)

    def _s(self, t):
        return np.clip(np.asarray(t, dtype=float) / self.duration_s, 0.0, 1.0)

    def center(self, t):
        return _out(self._center_interp(self._s(t)), t)

    def depth(self, t):
        return _out(self._depth_interp(self._s(t)), t)

    def center_velocity(self, t):
        return _out(self._center_interp.derivative(1)(self._s(t)) / self.duration_s, t)

    def center_acceleration(self, t):
        return _out(self._center_interp.derivative(2)(self._s(t)) / self.duration_s**2, t)

    def center_jerk(self, t):
        return _out(self._center_interp.derivative(3)(self._s(t)) / self.duration_s**3, t)

    def depth_rate(self, t):
        return _out(self._depth_interp.derivative(1)(self._s(t)) / self.duration_s, t)

    def control_table(self):
        return {
            "control_s": self.control_s.tolist(),
            "center_values_m": self.center_values_m.tolist(),
            "depth_values_j": self.depth_values_j.tolist(),
        }

    def describe(self):
        return {
            "type": self.name,
            "duration_us": self.duration_s * 1e6,
            "control_points": int(self.control_s.size),
            "initial_center_um": float(self.center_values_m[0]) * 1e6,
            "final_center_um": float(self.center_values_m[-1]) * 1e6,
            "initial_depth_uK": float(self.depth_values_j[0]) / 1.380649e-29,
            "final_depth_uK": float(self.depth_values_j[-1]) / 1.380649e-29,
        }


def validate_pickup_waveform(waveform, final_center_m, final_depth_j,
                              direction="pickup", overshoot_fraction=0.0,
                              ml_overshoot_fraction=0.05):
    """pick-up/drop-off（时间反演）波形校验：端点、范围、单调性与过冲。

    - 顺序 pick-up 族：中心与深度必须单调（pickup 不减、dropoff 不增）且零过冲；
    - ML 三次样条：允许非单调，但深度/位置越界幅度不得超过
      ml_overshoot_fraction（论文 ML 轨迹本身允许摆动，仅防病态过冲）；
    - 端点严格：pickup 为 c(0)=0→c(T)=d、D(0)=0→D(T)=D0；dropoff 反之；
    - 端点中心速度/加速度为零（constant-jerk 剖面性质；ML 样条不强制，
      其端点导数参与优化目标约束，由优化器惩罚而非此处拒绝）。
    """
    reasons = []
    table = waveform.dense_table(2001)
    center = table["aod_center_m"]
    depth = table["aod_depth_J"]
    span = abs(final_center_m)
    c_tol = 1e-9 * max(span, 1e-9) + 1e-15
    d_tol = 1e-9 * final_depth_j + 1e-30

    if direction == "pickup":
        endpoints = [(center[0], 0.0), (center[-1], final_center_m),
                     (depth[0], 0.0), (depth[-1], final_depth_j)]
    else:
        endpoints = [(center[0], final_center_m), (center[-1], 0.0),
                     (depth[0], final_depth_j), (depth[-1], 0.0)]
    labels = ["c(0)", "c(T)", "D(0)", "D(T)"]
    for (value, target), label in zip(endpoints, labels):
        scale = span if label.startswith("c") else final_depth_j
        tol = c_tol if label.startswith("c") else d_tol
        if abs(value - target) > tol:
            reasons.append(f"{label}={value:.6e} 不等于端点目标 {target:.6e}")

    if isinstance(waveform, PickupSequentialWaveform):
        dc = np.diff(center)
        dd = np.diff(depth)
        if direction == "pickup":
            if np.any(dc < -1e-9 * span):
                reasons.append("pick-up 中心轨迹非单调不减")
            if np.any(dd < -1e-9 * final_depth_j):
                reasons.append("pick-up 深度轨迹非单调不减")
        else:
            if np.any(dc > 1e-9 * span):
                reasons.append("drop-off 中心轨迹非单调不增")
            if np.any(dd > 1e-9 * final_depth_j):
                reasons.append("drop-off 深度轨迹非单调不增")
        for probe in (0.0, waveform.duration_s):
            if abs(float(waveform.center_velocity(probe))) > 1e-12 * span / waveform.duration_s:
                reasons.append(f"端点中心速度非零 t={probe:.3e}")
            if abs(float(waveform.center_acceleration(probe))) > 1e-9 * span / waveform.duration_s**2:
                reasons.append(f"端点中心加速度非零 t={probe:.3e}")

    if isinstance(waveform, MLCubicWaveform):
        c_over = ml_overshoot_fraction * max(span, 1e-12)
        d_over = ml_overshoot_fraction * final_depth_j
        lo_c, hi_c = (0.0, final_center_m) if direction == "pickup" else (0.0, final_center_m)
        if np.any(center < min(lo_c, hi_c) - c_over) or np.any(center > max(lo_c, hi_c) + c_over):
            reasons.append("ML 样条中心过冲超过容差")
        if np.any(depth < -d_over) or np.any(depth > final_depth_j + d_over):
            reasons.append("ML 样条深度过冲超过容差")
        if float(np.min(depth)) < -d_tol:
            reasons.append("ML 样条出现负深度（物理非法）")
    else:
        if np.any(center < -1e-9 * span - c_tol) or np.any(center > final_center_m + 1e-9 * span):
            reasons.append("中心越出 [0, d] 范围")
        if np.any(depth < -d_tol) or np.any(depth > final_depth_j * (1 + 1e-9)):
            reasons.append("深度越出 [0, D0] 范围")

    return (len(reasons) == 0), reasons


def validate_waveform(waveform, wcfg, final_center_target_m=0.0):
    """检查端点、范围、单调性和（同步族）参数约束；返回 (是否有效, 原因列表)。"""

    reasons = []
    table = waveform.dense_table(2001)
    center = table["aod_center_m"]
    depth = table["aod_depth_J"]
    span = abs(waveform.initial_center_m)
    depth_scale = waveform.initial_depth_j

    endpoint_tolerance = 1e-9 * max(span, 1e-9)
    if abs(center[0] - waveform.initial_center_m) > endpoint_tolerance:
        reasons.append(f"c(0)={center[0]:.6e} 不等于初始中心")
    if abs(center[-1] - final_center_target_m) > 1e-6 * max(span, 1e-9) + 1e-15:
        reasons.append(f"c(T)={center[-1]:.6e} 不等于目标终点")
    if abs(depth[0] - depth_scale) > 1e-9 * depth_scale:
        reasons.append(f"D(0)={depth[0]:.6e} 不等于初始深度")
    if depth[-1] > 1e-12 * depth_scale:
        reasons.append(f"D(T)={depth[-1]:.6e} 不为零")

    if np.any(np.diff(center) > 1e-9 * span):
        reasons.append("中心轨迹非单调不增")
    if np.any(np.diff(depth) > 1e-9 * depth_scale):
        reasons.append("深度轨迹非单调不增")
    if np.any(center < -1e-9 * span) or np.any(center > waveform.initial_center_m + 1e-9 * span):
        reasons.append("中心越出 [0, d] 范围")
    if np.any(depth < -1e-9 * depth_scale) or np.any(depth > depth_scale * (1 + 1e-9)):
        reasons.append("深度越出 [0, D0] 范围")

    if isinstance(waveform, OverlappedWaveform):
        p = waveform.params()
        min_width = wcfg["min_segment_fraction"]
        if not (0.0 <= p["move_start_fraction"] < p["move_end_fraction"] <= 1.0):
            reasons.append("移动窗口分数次序非法")
        if not (0.0 <= p["ramp_start_fraction"] < p["ramp_end_fraction"] <= 1.0):
            reasons.append("降深窗口分数次序非法")
        if p["move_end_fraction"] - p["move_start_fraction"] < min_width - 1e-12:
            reasons.append("移动窗口宽度低于下限")
        if p["ramp_end_fraction"] - p["ramp_start_fraction"] < min_width - 1e-12:
            reasons.append("降深窗口宽度低于下限")
        if not (wcfg["ramp_power_bounds"][0] <= p["ramp_power"] <= wcfg["ramp_power_bounds"][1]):
            reasons.append("ramp_power 超出边界")
        if p["ramp_power"] < 1.0:
            reasons.append("ramp_power<1 导致深度端点导数奇异")
        # 端点速度/加速度：精确端点必须为零；近端点（1e-·T 内）必须远小于峰值
        peak_speed = max(float(np.max(np.abs(table["center_velocity_m_per_s"]))), 1e-30)
        peak_depth_rate = max(float(np.max(np.abs(table["depth_rate_J_per_s"]))), 1e-30)
        eps_t = 1e-3 * waveform.duration_s
        for probe in (0.0, waveform.duration_s):
            if abs(float(waveform.center_velocity(probe))) > 1e-9 * peak_speed:
                reasons.append(f"端点中心速度非零 t={probe:.3e}")
                break
        for probe in (0.0, waveform.duration_s):
            if abs(float(waveform.depth_rate(probe))) > 1e-9 * peak_depth_rate:
                reasons.append(f"端点深度变化率非零 t={probe:.3e}")
                break
        for probe in (eps_t, waveform.duration_s - eps_t):
            if abs(float(waveform.center_velocity(probe))) > 0.05 * peak_speed:
                reasons.append(f"近端点中心速度超过峰值 5% t={probe:.3e}")
                break
        for probe in (eps_t, waveform.duration_s - eps_t):
            if abs(float(waveform.depth_rate(probe))) > 0.05 * peak_depth_rate:
                reasons.append(f"近端点深度变化率超过峰值 5% t={probe:.3e}")
                break

    if isinstance(waveform, SplineWaveform):
        if np.any(np.diff(waveform.center_values_m) > 0.0):
            reasons.append("样条中心控制点非单调不增")
        if np.any(np.diff(waveform.depth_values_j) > 0.0):
            reasons.append("样条深度控制点非单调不增")

    return (len(reasons) == 0), reasons
