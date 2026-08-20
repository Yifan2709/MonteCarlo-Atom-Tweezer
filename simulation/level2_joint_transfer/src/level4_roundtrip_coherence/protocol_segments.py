"""Level 4 分段协议时序：Protocol A/B、消融控制与连续拼接。

约定：路径坐标 p(t) 为 AOD 中心沿 long-move 方向离 SLM 的距离（米），
diagonal 几何下 c1 = c2 = p/√2。深度 ramp 默认 smootherstep（模型假设，
非论文逐点波形）。lensing 焦点偏移由瞬时轴速度按 Level 3 冻结模型计算。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from level3_3d_transport_lensing.gaussian_3d import KB
from level3_3d_transport_lensing.transport_profiles import (get_profile,
                                                            geometry_velocity_factors)
from level2_joint_transfer.waveforms import OverlappedWaveform

U_TO_KG = 1.66053906660e-27

CHECKPOINT_STAGES = {
    "pickup_end": "pickup",
    "split_end": "split",
    "outbound_end": "outbound_transport",
    "remote_hold_end": "remote_hold",
    "inbound_end": "inbound_transport",
    "merge_end": "merge",
    "dropoff_end": "dropoff",
    "final_hold_end": "final_hold",
}


def smootherstep(x):
    """S(x)=10x³-15x⁴+6x⁵，输入裁剪到 [0,1]。"""
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    return 10 * x ** 3 - 15 * x ** 4 + 6 * x ** 5


def smootherstep_prime(x):
    """dS/dx = 30x²(1-x)²（区间内），端点为零。"""
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    return 30.0 * x ** 2 * (1.0 - x) ** 2


@dataclass
class Segment:
    """一个协议分段；controls(tau) 返回路径坐标与深度的向量数组。"""

    name: str
    kind: str
    duration_s: float
    controls: Callable[[np.ndarray], dict]
    checkpoint: str | None = None
    expected_basin: str | None = None


def _constant_path(p):
    """恒定路径坐标控制。"""

    def controls(tau):
        n = np.shape(tau)
        return {"p": np.full(n, float(p)), "pd": np.zeros(n), "pdd": np.zeros(n)}

    return controls


def _profile_move(p0, p1, duration_s, profile_name):
    """p0→p1 的轨迹族移动（q 族归一化分数 s=tau/T）。"""
    prof = get_profile(profile_name)
    dp = p1 - p0

    def controls(tau):
        s = np.clip(tau / duration_s, 0.0, 1.0)
        return {"p": p0 + dp * prof.q(s),
                "pd": dp * prof.q_dot(s) / duration_s,
                "pdd": dp * prof.q_ddot(s) / duration_s ** 2}

    return controls


def _reversed_profile_move(p0, p1, duration_s, profile_name):
    """时间反向移动：p(τ)=p0+Δp(1−q(1−s))，速度同号导出、二阶导变号。"""
    prof = get_profile(profile_name)
    dp = p1 - p0

    def controls(tau):
        s = np.clip(tau / duration_s, 0.0, 1.0)
        u = 1.0 - s
        # p(τ)=p0+Δp(1−q(u))，du/dτ=−1/T：
        # ṗ=+Δp·q̇(u)/T，p̈=−Δp·q̈(u)/T²（时间反向的二阶导符号翻转）
        return {"p": p0 + dp * (1.0 - prof.q(u)),
                "pd": dp * prof.q_dot(u) / duration_s,
                "pdd": -dp * prof.q_ddot(u) / duration_s ** 2}

    return controls


def _depth_ramp(d0_j, d1_j, duration_s):
    """smootherstep 深度 ramp（零端点斜率）。"""

    def ramp(tau):
        x = np.clip(tau / duration_s, 0.0, 1.0)
        return (d0_j + (d1_j - d0_j) * smootherstep(x),
                (d1_j - d0_j) * smootherstep_prime(x) / duration_s)

    return ramp


def _constant_depth(d_j):
    """恒定深度 ramp。"""

    def ramp(tau):
        return (np.full(np.shape(tau), float(d_j)), np.zeros(np.shape(tau)))

    return ramp


def _combine(path_controls, slm_ramp, aod_ramp, override_aod=False):
    """组合路径控制与 SLM/AOD 深度 ramp。

    override_aod=True 时忽略 path_controls 自带的 AOD 深度（如 L2 波形段）。
    """

    def controls(tau):
        local = dict(path_controls(tau))
        d_slm, ddot_slm = slm_ramp(tau)
        local["d_slm"] = d_slm
        local["ddot_slm"] = ddot_slm
        if override_aod or "d_aod" not in local:
            d_aod, ddot_aod = aod_ramp(tau)
            local["d_aod"] = d_aod
            local["ddot_aod"] = ddot_aod
        return local

    return controls


def _l2_waveform_path(spec, d_init_center_m, d_depth_j, duration_s, reverse):
    """Level 2 冻结波形映射到路径坐标控制；reverse=True 为时间反向 pickup。

    正向（dropoff）：p(τ)=c(τ)（d→0），D(τ)=D_L2(τ)（D0→0）。
    反向（pickup）：p(τ)=c(T-τ)（0→d），D(τ)=D_L2(T-τ)（0→D0）；
    一阶导数经链式法则取负、二阶导数不变——不是简单倒序数组。
    """
    wave = OverlappedWaveform(duration_s, spec["move_start_fraction"],
                              spec["move_end_fraction"], spec["ramp_start_fraction"],
                              spec["ramp_end_fraction"], spec["ramp_power"],
                              d_init_center_m, d_depth_j)

    def controls(tau):
        t = duration_s - tau if reverse else tau
        center = np.asarray(wave.center(t), dtype=float)
        depth = np.asarray(wave.depth(t), dtype=float)
        vel = np.asarray(wave.center_velocity(t), dtype=float)
        acc = np.asarray(wave.center_acceleration(t), dtype=float)
        drate = np.asarray(wave.depth_rate(t), dtype=float)
        if reverse:
            vel, drate = -vel, -drate
        return {"p": center, "pd": vel, "pdd": acc,
                "d_aod": depth, "ddot_aod": drate}

    return controls


@dataclass
class Physics4:
    """Level 4 总物理参数（SI 单位）。"""

    mass_kg: float
    slm_waist_m: float
    slm_zr_m: float
    slm_center_m: tuple
    aod_waist_m: float
    aod_zr_m: float
    v_s_m_per_s: float | None
    axis_signs: tuple
    geometry: str

    @property
    def geometry_factors(self):
        """路径速度到两轴速度的系数。"""
        return geometry_velocity_factors(self.geometry)

    def lensing_shifts(self, pd, pdd):
        """由路径速度/加速度返回 (zs1, zs2, zs1d, zs2d)。"""
        f1, f2 = self.geometry_factors
        if self.v_s_m_per_s is None:
            z = np.zeros_like(np.asarray(pd, dtype=float))
            return z, z.copy(), z.copy(), z.copy()
        c1d, c2d = pd * f1, pd * f2
        c1dd, c2dd = pdd * f1, pdd * f2
        zr = self.aod_zr_m
        zs1 = 2.0 * self.axis_signs[0] * zr * c1d / self.v_s_m_per_s
        zs2 = 2.0 * self.axis_signs[1] * zr * c2d / self.v_s_m_per_s
        zs1d = 2.0 * self.axis_signs[0] * zr * c1dd / self.v_s_m_per_s
        zs2d = 2.0 * self.axis_signs[1] * zr * c2dd / self.v_s_m_per_s
        return zs1, zs2, zs1d, zs2d

    def centers_from_path(self, p):
        """路径坐标 → (c1, c2)。"""
        f1, f2 = self.geometry_factors
        return p * f1, p * f2


def physics_from_config(cfg, v_s=None) -> Physics4:
    """由 Level 4 配置构造物理参数。"""
    return Physics4(
        mass_kg=cfg.mass_u * U_TO_KG,
        slm_waist_m=cfg.slm_waist_um * 1e-6,
        slm_zr_m=np.pi * (cfg.slm_waist_um * 1e-6) ** 2
        / (cfg.slm_wavelength_nm * 1e-9),
        slm_center_m=(cfg.slm_center_um[0] * 1e-6, cfg.slm_center_um[1] * 1e-6),
        aod_waist_m=cfg.aod_waist_um * 1e-6,
        aod_zr_m=np.pi * (cfg.aod_waist_um * 1e-6) ** 2
        / (cfg.aod_wavelength_nm * 1e-9),
        v_s_m_per_s=v_s,
        axis_signs=cfg.lensing_axis_signs,
        geometry=cfg.long_geometry,
    )


class ProtocolTimeline:
    """分段协议时间线：连续拼接、网格化控制量与 DD 脉冲时刻。"""

    def __init__(self, name: str, segments: list[Segment], physics: Physics4):
        self.name = name
        self.segments = segments
        self.physics = physics
        self.segment_starts = []
        t = 0.0
        for seg in segments:
            if seg.duration_s <= 0:
                raise ValueError(f"分段 {seg.name} 时长非正")
            self.segment_starts.append(t)
            t += seg.duration_s
        self.total_duration_s = t

    def _segment_index(self, t: float) -> int:
        for i in range(len(self.segments) - 1, -1, -1):
            if t >= self.segment_starts[i] - 1e-12:
                return i
        return 0

    def control_at(self, t: float) -> dict:
        """标量时刻的完整控制量（路径坐标、中心、深度与 lensing 偏移）。"""
        t = float(np.clip(t, 0.0, self.total_duration_s))
        idx = self._segment_index(t)
        tau = t - self.segment_starts[idx]
        local = self.segments[idx].controls(np.array([tau]))
        return self._expand(local, idx, np.array([t]))[0]

    def controls_on_grid(self, times: np.ndarray) -> dict:
        """整条时间网格上的分段控制量数组（分段边界样本归左侧分段）。"""
        n = len(times)
        keys = ("p", "c1", "c2", "c1d", "c2d", "c1dd", "c2dd",
                "zs1", "zs2", "zs1d", "zs2d",
                "d_slm", "ddot_slm", "d_aod", "ddot_aod")
        out = {key: np.zeros(n) for key in keys}
        out["segment_index"] = np.zeros(n, dtype=int)
        covered = np.zeros(n, dtype=bool)
        for i, seg in enumerate(self.segments):
            start = self.segment_starts[i]
            mask = ((times >= start - 1e-12)
                    & (times <= start + seg.duration_s + 1e-12))
            if i > 0:
                mask = mask & ~covered
            if not np.any(mask):
                continue
            local = seg.controls(times[mask] - start)
            rows = self._expand(local, i, times[mask])
            for key in keys:
                out[key][mask] = [row[key] for row in rows]
            out["segment_index"][mask] = i
            covered |= mask
        if not np.all(covered):
            raise ValueError("时间网格存在未被分段覆盖的样本")
        return out

    def _expand(self, local, seg_idx, times_abs) -> list[dict]:
        """把分段本地控制量展开为完整控制字典。"""
        phys = self.physics
        p, pd, pdd = local["p"], local["pd"], local["pdd"]
        c1, c2 = phys.centers_from_path(p)
        zs1, zs2, zs1d, zs2d = phys.lensing_shifts(pd, pdd)
        f1, f2 = phys.geometry_factors
        rows = []
        for i in range(len(p)):
            rows.append({
                "t": times_abs[i], "segment": self.segments[seg_idx].name,
                "p": p[i], "c1": c1[i], "c2": c2[i],
                "c1d": pd[i] * f1, "c2d": pd[i] * f2,
                "c1dd": pdd[i] * f1, "c2dd": pdd[i] * f2,
                "zs1": zs1[i], "zs2": zs2[i], "zs1d": zs1d[i], "zs2d": zs2d[i],
                "d_slm": local["d_slm"][i], "ddot_slm": local["ddot_slm"][i],
                "d_aod": local["d_aod"][i], "ddot_aod": local["ddot_aod"][i],
            })
        return rows

    def build_grid(self, dt: float) -> dict:
        """以步长 dt 构建整条协议的控制网格与 checkpoint 信息。"""
        for seg in self.segments:
            ratio = seg.duration_s / dt
            if abs(ratio - round(ratio)) > 1e-9:
                raise ValueError(
                    f"分段 {seg.name} 时长 {seg.duration_s * 1e6:g} μs 不能被 "
                    f"dt={dt * 1e6:g} μs 整除")
        n_total = int(round(self.total_duration_s / dt))
        times = np.arange(n_total + 1) * dt
        grid = self.controls_on_grid(times)
        grid["times"] = times
        grid["dt"] = dt
        grid["n_steps"] = n_total
        checkpoints = []
        for i, seg in enumerate(self.segments):
            if seg.checkpoint is not None:
                step = int(round((self.segment_starts[i] + seg.duration_s) / dt))
                checkpoints.append({
                    "name": seg.checkpoint, "segment": seg.name,
                    "segment_index": i, "step": step,
                    "time_s": step * dt, "expected_basin": seg.expected_basin,
                    "stage": CHECKPOINT_STAGES[seg.checkpoint],
                })
        grid["checkpoints"] = checkpoints
        grid["segment_names"] = [s.name for s in self.segments]
        grid["segment_kinds"] = [s.kind for s in self.segments]
        grid["segment_steps"] = [int(round(self.segment_starts[i] / dt))
                                 for i in range(len(self.segments))]
        return grid

    def dd_pulse_times(self, mode: str, dt: float, pulse_offset_s: float | None,
                       n_repeats: int = 1, custom_times_s=()) -> np.ndarray:
        """按 DD 模式返回脉冲时刻（秒；n_repeats>1 时按轮重复）。

        脉冲时刻经 _snap_off_grid 处理：保证与最近时间样本距离 ≥ 0.2·dt
        （分段边界都是时间样本，因此也自动远离边界）。
        """
        if pulse_offset_s is not None:
            base = np.asarray(pulse_offset_s, dtype=float)
        else:
            base = None
        if mode == "none":
            pulses = np.empty(0)
        elif mode == "spin_echo":
            t_mid = 0.5 * self.total_duration_s
            if pulse_offset_s is not None:
                t_mid = t_mid + float(pulse_offset_s)
            pulses = np.array([t_mid])
        elif mode == "xy4_per_long_move":
            pulses = []
            for i, seg in enumerate(self.segments):
                if seg.kind in ("outbound_long_move", "inbound_long_move"):
                    t0 = self.segment_starts[i]
                    for k in (1, 3, 5, 7):
                        pulses.append(t0 + seg.duration_s * k / 8.0)
            pulses = np.array(sorted(pulses))
        elif mode == "custom_pulse_times":
            pulses = np.array(sorted(float(t) for t in custom_times_s))
        else:
            raise KeyError(f"未知 DD 模式 {mode}")
        if pulse_offset_s is not None and pulses.size:
            pulses = pulses + base
        if pulses.size:
            pulses = np.array([_snap_off_grid(p, dt) for p in pulses])
        if n_repeats > 1:
            if pulses.size == 0:
                return np.empty(0)
            offsets = (np.arange(n_repeats) * self.total_duration_s
                       + pulses[:, None]).ravel()
            return np.sort(np.array([_snap_off_grid(p, dt) for p in offsets]))
        return pulses

    def validate_pulse_placement(self, pulses: np.ndarray, dt: float) -> list[str]:
        """脉冲不得落在时间样本或分段边界附近；返回问题列表。"""
        issues = []
        boundaries = np.array([0.0] + [self.segment_starts[i] + seg.duration_s
                                       for i, seg in enumerate(self.segments)])
        for t in np.asarray(pulses, dtype=float):
            if t < 0 or t > self.total_duration_s * (1 + 1e-12):
                issues.append(f"脉冲 {t * 1e6:.4f} μs 超出协议时间")
                continue
            if abs(t / dt - round(t / dt)) < 1e-6:
                issues.append(f"脉冲 {t * 1e6:.4f} μs 恰好落在时间样本上")
            if np.min(np.abs(boundaries - t)) < dt / 4:
                issues.append(f"脉冲 {t * 1e6:.4f} μs 过于接近分段边界")
        return issues

    def reversed(self) -> "ProtocolTimeline":
        """控制序列的严格时间反向（preflight 可逆性检查用，非物理协议）。"""
        reversed_segments = []
        for seg in reversed(self.segments):
            inner, duration = seg.controls, seg.duration_s

            def make_controls(inner=inner, duration=duration):
                def controls(tau):
                    local = inner(duration - tau)
                    out = dict(local)
                    out["pd"] = -local["pd"]
                    if "ddot_aod" in local:
                        out["ddot_aod"] = -local["ddot_aod"]
                    if "ddot_slm" in local:
                        out["ddot_slm"] = -local["ddot_slm"]
                    return out
                return controls

            reversed_segments.append(Segment(seg.name, seg.kind, duration,
                                             make_controls()))
        return ProtocolTimeline(f"{self.name}_reversed", reversed_segments,
                                self.physics)


def _uk(depth_uK: float) -> float:
    """μK → J。"""
    return depth_uK * 1e-6 * KB


def _snap_off_grid(t: float, dt: float, margin: float = 0.2) -> float:
    """把脉冲时刻移到与最近时间样本距离 ≥ margin·dt 的位置。

    分段边界均为时间样本，因此该规则同时保证脉冲远离分段边界；
    偏移量 ≤ 0.5·dt，物理上可忽略（理想瞬时脉冲模型）。
    """
    frac = (t / dt) % 1.0
    if min(frac, 1.0 - frac) >= margin:
        return t
    return t + 0.5 * dt


def build_protocol_a(cfg, v_s, long_profile=None, split_profile=None) -> ProtocolTimeline:
    """Protocol A：论文时序启发的完整 round-trip anchor。"""
    d_slm0 = _uk(cfg.slm_initial_depth_uK)
    d_slm_low = _uk(cfg.slm_transport_depth_uK)
    d_aod = _uk(cfg.aod_transport_depth_uK)
    split = cfg.split_distance_um * 1e-6
    far = split + cfg.long_distance_um * 1e-6
    us = 1e-6
    long_prof = long_profile or cfg.long_profile
    split_prof = split_profile or cfg.split_profile
    aod_on = _constant_depth(d_aod)
    slm_low = _constant_depth(d_slm_low)
    segments = [
        Segment("initial_slm_hold", "hold", cfg.initial_hold_us * us,
                _combine(_constant_path(0.0), _constant_depth(d_slm0),
                         _constant_depth(0.0))),
        Segment("pickup", "pickup", cfg.pickup_us * us,
                _combine(_constant_path(0.0),
                         _depth_ramp(d_slm0, d_slm_low, cfg.pickup_us * us),
                         _depth_ramp(0.0, d_aod, cfg.pickup_us * us)),
                checkpoint="pickup_end", expected_basin="shared_or_ambiguous"),
        Segment("split", "split", cfg.split_us * us,
                _combine(_profile_move(0.0, split, cfg.split_us * us, split_prof),
                         slm_low, aod_on),
                checkpoint="split_end", expected_basin="bound_to_aod"),
        Segment("outbound_long_move", "outbound_long_move", cfg.outbound_us * us,
                _combine(_profile_move(split, far, cfg.outbound_us * us, long_prof),
                         slm_low, aod_on),
                checkpoint="outbound_end", expected_basin="bound_to_aod"),
        Segment("remote_operation_hold", "remote_hold",
                cfg.remote_operation_hold_us * us,
                _combine(_constant_path(far), slm_low, aod_on),
                checkpoint="remote_hold_end", expected_basin="bound_to_aod"),
        Segment("inbound_long_move", "inbound_long_move", cfg.inbound_us * us,
                _combine(_reversed_profile_move(far, split, cfg.inbound_us * us,
                                                long_prof), slm_low, aod_on),
                checkpoint="inbound_end", expected_basin="bound_to_aod"),
        Segment("merge", "merge", cfg.merge_us * us,
                _combine(_reversed_profile_move(split, 0.0, cfg.merge_us * us,
                                                split_prof), slm_low, aod_on),
                checkpoint="merge_end", expected_basin="shared_or_ambiguous"),
        Segment("dropoff", "dropoff", cfg.dropoff_us * us,
                _combine(_constant_path(0.0),
                         _depth_ramp(d_slm_low, d_slm0, cfg.dropoff_us * us),
                         _depth_ramp(d_aod, 0.0, cfg.dropoff_us * us)),
                checkpoint="dropoff_end", expected_basin="bound_to_slm"),
        Segment("final_slm_hold", "final_hold", cfg.final_hold_us * us,
                _combine(_constant_path(0.0), _constant_depth(d_slm0),
                         _constant_depth(0.0)),
                checkpoint="final_hold_end", expected_basin="bound_to_slm"),
    ]
    return ProtocolTimeline("protocol_a", segments, physics_from_config(cfg, v_s))


def build_protocol_a_no_remote_hold(cfg, v_s) -> ProtocolTimeline:
    """Protocol A 变体：删除 remote hold（消融）。"""
    timeline = build_protocol_a(cfg, v_s)
    keep = [s for s in timeline.segments if s.kind != "remote_hold"]
    return ProtocolTimeline("protocol_a_no_remote_hold", keep, timeline.physics)


def build_protocol_b(cfg, level2_spec, v_s) -> ProtocolTimeline:
    """Protocol B：冻结 Level 2 波形 pickup/dropoff + Level 3 冻结长运输。

    Level 2 波形语义：dropoff 自带 2.4 μm merge 移动，因此本协议不再添加
    独立 split/merge 段（避免重复计算）；SLM 深度沿用 Level 2 恒定 140 μK。
    """
    d_slm_b = _uk(cfg.slm_protocol_b_depth_uK)
    d_aod = _uk(cfg.aod_transport_depth_uK)
    split = cfg.split_distance_um * 1e-6
    far = split + cfg.long_distance_um * 1e-6
    us = 1e-6
    l2_dur = float(level2_spec["duration_us"]) * us
    pickup_path = _l2_waveform_path(level2_spec, split, d_aod, l2_dur, reverse=True)
    dropoff_path = _l2_waveform_path(level2_spec, split, d_aod, l2_dur,
                                     reverse=False)
    slm_const = _constant_depth(d_slm_b)
    aod_off = _constant_depth(0.0)
    aod_on = _constant_depth(d_aod)
    segments = [
        Segment("initial_slm_hold", "hold", cfg.initial_hold_us * us,
                _combine(_constant_path(0.0), slm_const, aod_off)),
        Segment("pickup", "pickup", l2_dur,
                _combine(pickup_path, slm_const, aod_off, override_aod=False),
                checkpoint="pickup_end", expected_basin="bound_to_aod"),
        Segment("outbound_long_move", "outbound_long_move", cfg.outbound_us * us,
                _combine(_profile_move(split, far, cfg.outbound_us * us,
                                       cfg.long_profile), slm_const, aod_on),
                checkpoint="outbound_end", expected_basin="bound_to_aod"),
        Segment("remote_operation_hold", "remote_hold",
                cfg.remote_operation_hold_us * us,
                _combine(_constant_path(far), slm_const, aod_on),
                checkpoint="remote_hold_end", expected_basin="bound_to_aod"),
        Segment("inbound_long_move", "inbound_long_move", cfg.inbound_us * us,
                _combine(_reversed_profile_move(far, split, cfg.inbound_us * us,
                                                cfg.long_profile), slm_const,
                         aod_on),
                checkpoint="inbound_end", expected_basin="bound_to_aod"),
        Segment("dropoff", "dropoff", l2_dur,
                _combine(dropoff_path, slm_const, aod_off, override_aod=False),
                checkpoint="dropoff_end", expected_basin="bound_to_slm"),
        Segment("final_slm_hold", "final_hold", cfg.final_hold_us * us,
                _combine(_constant_path(0.0), slm_const, aod_off),
                checkpoint="final_hold_end", expected_basin="bound_to_slm"),
    ]
    return ProtocolTimeline("protocol_b", segments, physics_from_config(cfg, v_s))


def build_static_idle(cfg, total_us=None, depth_uK=None) -> ProtocolTimeline:
    """static_idle_same_duration 消融：SLM 恒深、AOD 关闭。"""
    total = total_us if total_us is not None else protocol_a_total_us(cfg)
    depth = _uk(depth_uK if depth_uK is not None else cfg.slm_initial_depth_uK)
    segments = [
        Segment("static_idle", "final_hold", total * 1e-6,
                _combine(_constant_path(0.0), _constant_depth(depth),
                         _constant_depth(0.0)),
                checkpoint="final_hold_end", expected_basin="bound_to_slm"),
    ]
    return ProtocolTimeline("static_idle", segments, physics_from_config(cfg, None))


def build_transfer_only(cfg, v_s) -> ProtocolTimeline:
    """transfer_only 消融：只保留 hold+pickup+dropoff+hold（无任何移动段）。"""
    d_slm0 = _uk(cfg.slm_initial_depth_uK)
    d_slm_low = _uk(cfg.slm_transport_depth_uK)
    d_aod = _uk(cfg.aod_transport_depth_uK)
    us = 1e-6
    segments = [
        Segment("initial_slm_hold", "hold", cfg.initial_hold_us * us,
                _combine(_constant_path(0.0), _constant_depth(d_slm0),
                         _constant_depth(0.0))),
        Segment("pickup", "pickup", cfg.pickup_us * us,
                _combine(_constant_path(0.0),
                         _depth_ramp(d_slm0, d_slm_low, cfg.pickup_us * us),
                         _depth_ramp(0.0, d_aod, cfg.pickup_us * us)),
                checkpoint="pickup_end", expected_basin="shared_or_ambiguous"),
        Segment("dropoff", "dropoff", cfg.dropoff_us * us,
                _combine(_constant_path(0.0),
                         _depth_ramp(d_slm_low, d_slm0, cfg.dropoff_us * us),
                         _depth_ramp(d_aod, 0.0, cfg.dropoff_us * us)),
                checkpoint="dropoff_end", expected_basin="bound_to_slm"),
        Segment("final_slm_hold", "final_hold", cfg.final_hold_us * us,
                _combine(_constant_path(0.0), _constant_depth(d_slm0),
                         _constant_depth(0.0)),
                checkpoint="final_hold_end", expected_basin="bound_to_slm"),
    ]
    return ProtocolTimeline("transfer_only", segments, physics_from_config(cfg, v_s))


def build_transport_only(cfg, v_s) -> ProtocolTimeline:
    """transport_only 消融：AOD-only 长运输（Level 3 式），SLM 全程关闭。

    共同初态数组被解释为起点 AOD 阱内的相对六维态（CRN 映射），用于隔离
    长运输段的损失/加热/相位贡献。
    """
    d_aod = _uk(cfg.aod_transport_depth_uK)
    split = cfg.split_distance_um * 1e-6
    far = split + cfg.long_distance_um * 1e-6
    us = 1e-6
    slm_off = _constant_depth(0.0)
    aod_on = _constant_depth(d_aod)
    segments = [
        Segment("split", "split", cfg.split_us * us,
                _combine(_profile_move(0.0, split, cfg.split_us * us,
                                       cfg.split_profile), slm_off, aod_on),
                checkpoint="split_end", expected_basin="bound_to_aod"),
        Segment("outbound_long_move", "outbound_long_move", cfg.outbound_us * us,
                _combine(_profile_move(split, far, cfg.outbound_us * us,
                                       cfg.long_profile), slm_off, aod_on),
                checkpoint="outbound_end", expected_basin="bound_to_aod"),
        Segment("remote_operation_hold", "remote_hold",
                cfg.remote_operation_hold_us * us,
                _combine(_constant_path(far), slm_off, aod_on),
                checkpoint="remote_hold_end", expected_basin="bound_to_aod"),
        Segment("inbound_long_move", "inbound_long_move", cfg.inbound_us * us,
                _combine(_reversed_profile_move(far, split, cfg.inbound_us * us,
                                                cfg.long_profile), slm_off,
                         aod_on),
                checkpoint="inbound_end", expected_basin="bound_to_aod"),
        Segment("merge", "merge", cfg.merge_us * us,
                _combine(_reversed_profile_move(split, 0.0, cfg.merge_us * us,
                                                cfg.split_profile), slm_off,
                         aod_on),
                checkpoint="merge_end", expected_basin="bound_to_aod"),
        Segment("final_aod_hold", "final_hold", cfg.final_hold_us * us,
                _combine(_constant_path(0.0), slm_off, aod_on),
                checkpoint="final_hold_end", expected_basin="bound_to_aod"),
    ]
    return ProtocolTimeline("transport_only", segments, physics_from_config(cfg, v_s))


def protocol_a_total_us(cfg) -> float:
    """Protocol A 总时长（μs）。"""
    return (cfg.initial_hold_us + cfg.pickup_us + cfg.split_us + cfg.outbound_us
            + cfg.remote_operation_hold_us + cfg.inbound_us + cfg.merge_us
            + cfg.dropoff_us + cfg.final_hold_us)


def protocol_b_total_us(cfg, level2_duration_us: float) -> float:
    """Protocol B 总时长（μs）。"""
    return (cfg.initial_hold_us + 2.0 * level2_duration_us + cfg.outbound_us
            + cfg.remote_operation_hold_us + cfg.inbound_us + cfg.final_hold_us)
