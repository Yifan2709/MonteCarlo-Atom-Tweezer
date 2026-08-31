"""Level 2 三类 AOD 波形：顺序基线、压缩基线、同步受约束族与单调样条精修。

所有波形统一返回 AOD 中心 c(t)、深度 D(t) 及一至三阶导数（导数用于功-能核算
与光滑度诊断）。时间输入为秒，支持标量与 NumPy 数组。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator


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
