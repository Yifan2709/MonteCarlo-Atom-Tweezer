"""Level 2B 剖面族：梯形恒速巡航、窗口化变速与纯移动剖面。

全部类复用 level2_joint_transfer.waveforms.WaveformBase 接口与
OverlappedWaveform 的降深窗口机制，只新增"移动侧"的解析剖面：
- trapezoid：窗口内 smootherstep 加速 → 恒速巡航 → smootherstep 减速；
  加速/减速段各贡献 v_c·f·T/2 位移，巡航段贡献 v_c·(1−2f)·T，
  故巡航速度 v_c = d/[(1−f)·T_move]，峰值速度严格等于 v_c；
- smootherstep：与 OverlappedWaveform 完全相同的五次 smootherstep 移动；
- depth_mode="constant"：深度恒定（纯移动实验）。
"""
from __future__ import annotations

import numpy as np

from level2_joint_transfer.waveforms import OverlappedWaveform, WaveformBase


def _out(value, original):
    """标量输入返回 Python float，数组输入原样返回（与 level2 约定一致）。"""
    return float(value) if np.asarray(original).ndim == 0 else value


def _smootherstep(z):
    """五次 smootherstep S(z)=10z^3-15z^4+6z^5，z 裁剪到 [0,1]。"""
    zc = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    return 10 * zc**3 - 15 * zc**4 + 6 * zc**5


def _smootherstep_int(z):
    """I5(z)=∫₀ᶻS(u)du=2.5z⁴-3z⁵+z⁶；I5(1)=0.5。"""
    zc = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    return 2.5 * zc**4 - 3.0 * zc**5 + zc**6


def _s5_prime(z):
    """S'(z)=30z²-60z³+30z⁴。"""
    return 30 * z**2 - 60 * z**3 + 30 * z**4


def _s5_double(z):
    """S''(z)=60z-180z²+120z³。"""
    return 60 * z - 180 * z**2 + 120 * z**3


class SpeedWaveform(OverlappedWaveform):
    """速度研究用统一剖面：移动窗口内 trapezoid 或 smootherstep，降深可选恒定。

    参数与 OverlappedWaveform 兼容，额外要求 move_profile ∈
    {"trapezoid","smootherstep"}、accel_fraction ∈ (0,0.5]、
    depth_mode ∈ {"ramp_window","constant"}。
    """

    def __init__(self, duration_s, move_start_fraction, move_end_fraction,
                 ramp_start_fraction, ramp_end_fraction, ramp_power,
                 initial_center_m, initial_depth_j,
                 move_profile="trapezoid", accel_fraction=0.30,
                 depth_mode="ramp_window"):
        super().__init__(duration_s, move_start_fraction, move_end_fraction,
                         ramp_start_fraction, ramp_end_fraction, ramp_power,
                         initial_center_m, initial_depth_j)
        if move_profile not in ("trapezoid", "smootherstep"):
            raise ValueError(f"未知移动剖面 {move_profile!r}")
        if not (0.0 < accel_fraction <= 0.5):
            raise ValueError("accel_fraction 必须位于 (0, 0.5]")
        if depth_mode not in ("ramp_window", "constant"):
            raise ValueError(f"未知深度模式 {depth_mode!r}")
        self.move_profile = str(move_profile)
        self.accel_fraction = float(accel_fraction)
        self.depth_mode = str(depth_mode)
        self.move_t0 = self.move_start_fraction * self.duration_s
        self.move_t1 = self.move_end_fraction * self.duration_s
        self.move_duration_s = self.move_t1 - self.move_t0
        if self.move_duration_s <= 0.0:
            raise ValueError("移动窗口宽度必须为正")
        self.cruise_speed_m_per_s = self.initial_center_m / (
            (1.0 - self.accel_fraction) * self.move_duration_s)

    # ------------------------------------------------------------------ 移动侧
    def _tau(self, t):
        """相对移动窗口起点的本地时间，裁剪到 [0, T_move]。"""
        return np.clip(np.asarray(t, dtype=float) - self.move_t0, 0.0, self.move_duration_s)

    def _phases(self, t):
        """返回 (inside, accel, decel) 布尔掩码。"""
        t_arr = np.asarray(t, dtype=float)
        tau = self._tau(t)
        f, tm = self.accel_fraction, self.move_duration_s
        inside = (t_arr > self.move_t0) & (t_arr < self.move_t1)
        accel = (tau <= f * tm) & inside
        decel = (tau >= (1.0 - f) * tm) & inside
        return inside, accel, decel

    def center(self, t):
        if self.move_profile == "smootherstep":
            return super().center(t)
        tau = self._tau(t)
        f, tm = self.accel_fraction, self.move_duration_s
        v_c, d = self.cruise_speed_m_per_s, self.initial_center_m
        # 相位按裁剪后的 τ 判定（不排除窗口端点）：三分支公式在边界
        # 连续闭合（τ=0→d，τ=T_move→0），窗口外经裁剪自然冻结端点值。
        accel = tau <= f * tm
        decel = tau >= (1.0 - f) * tm
        # 加速段：c = d − v_c·f·T·I5(τ/(fT))
        c_acc = d - v_c * f * tm * _smootherstep_int(tau / (f * tm))
        # 巡航段：c = d − v_c·(f·T/2 + τ − f·T)
        c_cru = d - v_c * (0.5 * f * tm + (tau - f * tm))
        # 减速段：z' = (T−τ)/(fT)，c = d − v_c·[(1−f)T − f·T·I5(z')]
        z_dec = (tm - tau) / (f * tm)
        c_dec = d - v_c * ((1.0 - f) * tm - f * tm * _smootherstep_int(z_dec))
        value = np.where(accel, c_acc, np.where(decel, c_dec, c_cru))
        return _out(value, t)

    def center_velocity(self, t):
        if self.move_profile == "smootherstep":
            return super().center_velocity(t)
        tau = self._tau(t)
        f, tm = self.accel_fraction, self.move_duration_s
        v_c = self.cruise_speed_m_per_s
        inside, accel, decel = self._phases(t)
        v_acc = -v_c * _smootherstep(tau / (f * tm))
        v_dec = -v_c * _smootherstep((tm - tau) / (f * tm))
        value = np.where(accel, v_acc, np.where(decel, v_dec, np.where(inside, -v_c, 0.0)))
        return _out(value, t)

    def center_acceleration(self, t):
        if self.move_profile == "smootherstep":
            return super().center_acceleration(t)
        tau = self._tau(t)
        f, tm = self.accel_fraction, self.move_duration_s
        v_c = self.cruise_speed_m_per_s
        _inside, accel, decel = self._phases(t)
        a_acc = -v_c * _s5_prime(tau / (f * tm)) / (f * tm)
        a_dec = +v_c * _s5_prime((tm - tau) / (f * tm)) / (f * tm)
        value = np.where(accel, a_acc, np.where(decel, a_dec, 0.0))
        return _out(value, t)

    def center_jerk(self, t):
        if self.move_profile == "smootherstep":
            return super().center_jerk(t)
        tau = self._tau(t)
        f, tm = self.accel_fraction, self.move_duration_s
        v_c = self.cruise_speed_m_per_s
        _inside, accel, decel = self._phases(t)
        j_acc = -v_c * _s5_double(tau / (f * tm)) / (f * tm) ** 2
        j_dec = -v_c * _s5_double((tm - tau) / (f * tm)) / (f * tm) ** 2
        value = np.where(accel, j_acc, np.where(decel, j_dec, 0.0))
        return _out(value, t)

    # ------------------------------------------------------------------ 深度侧
    def depth(self, t):
        if self.depth_mode == "constant":
            return _out(np.full_like(np.asarray(t, dtype=float), self.initial_depth_j), t)
        return super().depth(t)

    def depth_rate(self, t):
        if self.depth_mode == "constant":
            return _out(np.zeros_like(np.asarray(t, dtype=float)), t)
        return super().depth_rate(t)

    def describe(self):
        data = super().describe()
        data.update({
            "move_profile": self.move_profile,
            "accel_fraction": self.accel_fraction,
            "depth_mode": self.depth_mode,
            "move_window_us": [self.move_t0 * 1e6, self.move_t1 * 1e6],
            "cruise_speed_mm_per_s": self.cruise_speed_m_per_s * 1e3,
        })
        return data


def build_speed_waveform(spec: dict, physics: dict):
    """由可序列化 spec 构造 Level 2B 剖面；level2 原有类型委托原实现。"""
    kind = spec["type"]
    if kind in ("sequential", "overlapped", "spline"):
        from level2_joint_transfer.level2_simulation import build_waveform
        return build_waveform(spec, physics)
    if kind == "speed":
        return SpeedWaveform(
            spec["duration_s"], spec["move_start_fraction"], spec["move_end_fraction"],
            spec["ramp_start_fraction"], spec["ramp_end_fraction"], spec["ramp_power"],
            physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
            move_profile=spec.get("move_profile", "trapezoid"),
            accel_fraction=spec.get("accel_fraction", 0.30),
            depth_mode=spec.get("depth_mode", "ramp_window"))
    raise ValueError(f"未知 Level 2B 剖面类型 {kind!r}")


# ---------------------------------------------------------------- 速度口径与步长
def speed_calibers(waveform: WaveformBase, num_points: int = 4001) -> dict:
    """返回剖面的三种速度口径（平均/峰值/均方根）与移动时长。"""
    table = waveform.dense_table(num_points)
    t = table["time_s"]
    c = table["aod_center_m"]
    v = table["center_velocity_m_per_s"]
    peak = float(np.max(np.abs(v)))
    moving = np.abs(v) > 1e-9 * max(1e-30, peak)
    total_displacement = float(c[0] - c[-1])
    if hasattr(waveform, "move_duration_s"):
        move_duration = float(waveform.move_duration_s)
    else:  # sequential 等无窗口属性：用速度非零区间的跨度
        idx = np.flatnonzero(moving)
        move_duration = float(t[idx[-1]] - t[idx[0]]) if idx.size >= 2 else float(t[-1])
    v_avg = total_displacement / move_duration if move_duration > 0 else 0.0
    v_rms = float(np.sqrt(np.mean(v[moving] ** 2))) if np.any(moving) else 0.0
    return {"move_duration_s": move_duration,
            "total_displacement_m": total_displacement,
            "v_avg_m_per_s": v_avg, "v_peak_m_per_s": peak, "v_rms_m_per_s": v_rms}


def max_step_for(v_peak_m_per_s: float, move_duration_s: float, waist_m: float,
                 base_dt_s: float, min_steps_per_move: int = 2000,
                 max_step_displacement_over_waist: float = 0.01) -> float:
    """自适应步长规则：dt ≤ min(base, T_move/N, frac·w/v_peak)。"""
    if v_peak_m_per_s <= 0.0 or move_duration_s <= 0.0 or waist_m <= 0.0:
        raise ValueError("速度、移动时长与束宽必须为正")
    return float(min(base_dt_s,
                     move_duration_s / max(1, int(min_steps_per_move)),
                     max_step_displacement_over_waist * waist_m / v_peak_m_per_s))


def snap_dt(total_s: float, dt_max: float) -> float:
    """把 dt 收缩为 total_s 的整数分之一，保证网格整除。"""
    if total_s <= 0.0 or dt_max <= 0.0:
        raise ValueError("总时长与步长上限必须为正")
    return total_s / max(1, int(np.ceil(total_s / dt_max)))
