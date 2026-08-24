"""有限脉冲自旋通道的批量传播引擎（RB 与 channel 表征共用）。

一条“移动通道”= 冻结轨迹池中的一条 u_SLM/u_AOD 时序 + DD 有限脉冲。
自由演化段用样本中点矩阵向量积一次性累积 z 相位（与逐步施加在分段
常数假设下严格一致）；脉冲段逐子步解析 SU(2) 推进。静态门（Clifford）
在静态池轨迹回放上施加，虚拟 Z 作为精确零时长帧操作按分解顺序插入。
所有 pulse family（bare/SCROFULOUS）共用同一传播器。大 batch 自动分块
以限制轨迹行物化的峰值内存。
"""
from __future__ import annotations

import numpy as np

from .dd_finite_pulses import dd_pulse_events
from .pulse_library import (PulseSpec, envelope_area_integral,
                            envelope_value, pulse_list_for)
from .quantum_propagators import apply_su2, su2_propagator
from .qubit_hamiltonian import HBAR
from .trajectory_replay import TrajectoryPool

TWO_PI = 2.0 * np.pi
DEFAULT_CHUNK = 4096


def _interp_rows(rows: np.ndarray, pos: np.ndarray):
    """rows [B, n_t]，pos [B, n_sub] → 线性插值 [B, n_sub]。"""
    n_t = rows.shape[1]
    i0 = np.floor(pos).astype(np.int64)
    frac = pos - i0
    v0 = np.take_along_axis(rows, i0, axis=1)
    v1 = np.take_along_axis(rows, np.minimum(i0 + 1, n_t - 1), axis=1)
    return v0 + (v1 - v0) * frac


class MoveChannelEngine:
    """沿冻结轨迹传播一批实例的移动通道（有限脉冲 DD）。"""

    def __init__(self, pool: TrajectoryPool, omega_base: float,
                 dd_mode: str = "xy4_per_long_move",
                 pulse_family: str = "bare", envelope: str = "square",
                 edge_us: float = 0.0,
                 long_move_kinds: tuple = ("outbound_long_move",
                                           "inbound_long_move")):
        self.pool = pool
        self.omega_base = float(omega_base)
        self.dd_mode = dd_mode
        self.edge_us = float(edge_us)
        bounds = pool.segment_boundaries_s
        windows = [(float(bounds[i]), float(bounds[i + 1]))
                   for i, name in enumerate(pool.segment_names)
                   if any(key in name for key in long_move_kinds)]
        self.long_move_windows = windows
        self.events = dd_pulse_events(
            dd_mode, 0.0, pool.total_duration_s, self.omega_base,
            family=pulse_family, envelope=envelope, edge_us=self.edge_us,
            long_move_windows=windows)
        self.n_pulses = len(self.events)

    def propagate(self, psi: np.ndarray, k_idx: np.ndarray,
                  eta_slm, eta_aod, delta_drive=0.0, delta_site=None,
                  noise_qs=None, dt_free_us=0.25, dt_pulse_us=0.1,
                  rabi_scale=None, amp_error=0.0,
                  free_noise_traces=None, chunk=DEFAULT_CHUNK) -> np.ndarray:
        """把 [B,2] 态矢推过整条移动通道（就地修改并返回）。"""
        n = k_idx.size
        eta_slm = np.broadcast_to(np.asarray(eta_slm, dtype=float), (n,))
        eta_aod = np.broadcast_to(np.asarray(eta_aod, dtype=float), (n,))
        base_offset = np.zeros(n)
        base_offset = base_offset + float(delta_drive)
        if delta_site is not None:
            base_offset = base_offset + np.asarray(delta_site, dtype=float)
        if noise_qs is not None:
            base_offset = base_offset + np.asarray(noise_qs, dtype=float)
        # rabi 仅为无量纲 per-instance 比例；峰值 Rabi 由脉冲面积/包络积分给出
        rabi = (np.ones(n) if rabi_scale is None
                else np.asarray(rabi_scale, dtype=float))
        boundaries = [0.0]
        for p in self.events:
            boundaries += [p.start, p.start + p.duration]
        boundaries.append(self.pool.total_duration_s)
        boundaries = np.unique(np.asarray(boundaries))
        for lo in range(0, n, chunk):
            hi = min(lo + chunk, n)
            sl = slice(lo, hi)
            rows_s = self.pool.u_slm[k_idx[sl]].astype(np.float32)
            rows_a = self.pool.u_aod[k_idx[sl]].astype(np.float32)
            noise_rows = (None if free_noise_traces is None
                          else free_noise_traces[sl])
            event_iter = iter(self.events)
            cur = next(event_iter, None)
            for a, b in zip(boundaries[:-1], boundaries[1:]):
                if b - a <= 1e-15:
                    continue
                if cur is not None and abs(a - cur.start) < 1e-12:
                    self._pulse(psi[sl], cur, rows_s, rows_a, a,
                                eta_slm[sl], eta_aod[sl],
                                base_offset[sl], rabi[sl], amp_error,
                                dt_pulse_us)
                    cur = next(event_iter, None)
                else:
                    self._free(psi[sl], a, b, rows_s, rows_a,
                               eta_slm[sl], eta_aod[sl], base_offset[sl],
                               dt_free_us, noise_rows)
        return psi

    def _free(self, psi, a, b, rows_s, rows_a, eta_s, eta_a, offset,
              dt_free_us, noise_rows):
        """自由演化段：中点矩形积分一次性施加 z 旋转。"""
        duration = b - a
        n_sub = max(1, int(np.ceil(duration / (dt_free_us * 1e-6) - 1e-9)))
        dt = duration / n_sub
        sample_dt = self.pool.sample_interval_s
        mids = a + (np.arange(n_sub) + 0.5) * dt
        pos = np.clip(mids / sample_dt, 0.0, rows_s.shape[1] - 1.000001)
        delta = (eta_s[:, None] * _interp_rows(rows_s, pos[None, :])
                 + eta_a[:, None] * _interp_rows(rows_a, pos[None, :])) / HBAR
        delta = delta + offset[:, None]
        if noise_rows is not None:
            s0 = np.clip(((mids - 0.5 * sample_dt) / sample_dt).astype(int),
                         0, noise_rows.shape[1] - 1)
            delta = delta + noise_rows[:, s0]
        phase = delta.sum(axis=1) * dt
        half = 0.5 * phase
        psi[:, 0] *= np.exp(-1j * half)
        psi[:, 1] *= np.exp(1j * half)

    def _pulse(self, psi, pulse: PulseSpec, rows_s, rows_a, a, eta_s, eta_a,
               offset, rabi, amp_error, dt_pulse_us):
        """脉冲段：逐子步解析 SU(2)（envelope 按面积补偿峰值 Rabi）。"""
        duration = pulse.duration
        n_sub = max(1, int(np.ceil(duration / (dt_pulse_us * 1e-6) - 1e-9)))
        dt = duration / n_sub
        env_area = envelope_area_integral(pulse.envelope, duration,
                                          pulse.edge_us)
        peak = pulse.nominal_area * (1.0 + pulse.amp_error + amp_error) \
            / max(env_area, 1e-12)
        sample_dt = self.pool.sample_interval_s
        for k in range(n_sub):
            t_mid = a + (k + 0.5) * dt
            pos = min(max(t_mid / sample_dt, 0.0), rows_s.shape[1] - 1.000001)
            i0 = int(np.floor(pos))
            frac = pos - i0
            us = rows_s[:, i0] + (rows_s[:, min(i0 + 1, rows_s.shape[1] - 1)]
                                  - rows_s[:, i0]) * frac
            ua = rows_a[:, i0] + (rows_a[:, min(i0 + 1, rows_a.shape[1] - 1)]
                                  - rows_a[:, i0]) * frac
            delta = (eta_s * us + eta_a * ua) / HBAR + offset
            omega = rabi * peak * envelope_value(
                pulse.envelope, (k + 0.5) / n_sub, pulse.edge_us,
                (k + 0.5) * dt * 1e6, duration * 1e6)
            u = su2_propagator(delta, omega, pulse.axis_phase, dt)
            apply_su2(psi, *u)

    def phase_diagnostics(self, k_idx, eta_slm, eta_aod):
        """自由段/脉冲段时间与瞬时极限 toggling 相位（诊断用）。"""
        pulse_time = sum(p.duration for p in self.events)
        centers = [p.start + 0.5 * p.duration for p in self.events]
        times = self.pool.times_s
        mids = 0.5 * (times[:-1] + times[1:])
        y = np.ones(mids.size)
        if centers:
            flips = np.searchsorted(mids, np.sort(np.asarray(centers)),
                                    side="left")
            for f in flips:
                y[f:] *= -1.0
        delta = ((eta_slm * self.pool.u_slm[k_idx].astype(np.float64)
                  + eta_aod * self.pool.u_aod[k_idx].astype(np.float64))
                 / HBAR)
        delta_mid = 0.5 * (delta[:, :-1] + delta[:, 1:])
        phi = delta_mid @ (y * self.pool.sample_interval_s)
        return {"free_evolution_time_s": float(self.pool.total_duration_s
                                               - pulse_time),
                "pulse_total_time_s": float(pulse_time),
                "n_dd_pulses": self.n_pulses,
                "instantaneous_limit_phase_rad": phi}


class StaticGateEngine:
    """静态条件下 Clifford 门（含 SCROFULOUS）的批量施加。

    门失谐 = 静态池轨迹回放（每实例一条 trace draw，按序列时间循环
    平铺）+ 站点静态失谐 + quasistatic 噪声 + Δ_drive；虚拟 Z 按分解
    顺序作为精确零时长操作插入。
    """

    def __init__(self, static_pool: TrajectoryPool, omega_base: float,
                 family: str = "scrofulous", envelope: str = "square",
                 group=None):
        self.pool = static_pool
        self.omega_base = float(omega_base)
        self.family = family
        self.envelope = envelope
        self.group = group
        self._ops_cache: dict[int, list] = {}

    def ops_for(self, clifford_index: int) -> list:
        """分解 op 列表：("VZ", angle) 与 ("PULSE", PulseSpec) 按时间序。"""
        if clifford_index not in self._ops_cache:
            ops = []
            for kind, val in self.group.decompose(clifford_index):
                if kind == "VZ":
                    ops.append(("VZ", float(val)))
                else:
                    axis_phase, angle = {
                        "X90": (0.0, 0.5 * np.pi),
                        "Y90": (0.5 * np.pi, 0.5 * np.pi),
                        "X-90": (np.pi, 0.5 * np.pi),
                        "Y-90": (1.5 * np.pi, 0.5 * np.pi)}[val]
                    for p in pulse_list_for(axis_phase, angle, self.family,
                                            self.omega_base, self.envelope,
                                            clifford_id=clifford_index):
                        ops.append(("PULSE", p))
            self._ops_cache[clifford_index] = ops
        return self._ops_cache[clifford_index]

    def clifford_duration_s(self, clifford_index: int) -> float:
        """该 Clifford 的物理脉冲总时长（虚拟 Z 零时长）。"""
        return float(sum(p.duration for k, p in self.ops_for(clifford_index)
                         if k == "PULSE"))

    def apply(self, psi, cliff_indices: np.ndarray, k_trace: np.ndarray,
              t_start, eta_slm, eta_aod, delta_site=None, noise_qs=None,
              rabi_scale=None, amp_error=0.0, n_sub_pulse=4, delta_drive=0.0,
              chunk=DEFAULT_CHUNK):
        """按每实例 Clifford 索引施加门（就地修改并返回）。

        k_trace: [B] 静态池轨迹索引；t_start: [B] 或标量（该实例的序列
        时间起点，trace 按 pool 总时长循环平铺）。
        """
        n = cliff_indices.size
        pool = self.pool
        sample_dt = pool.sample_interval_s
        total = pool.total_duration_s
        eta_slm = np.broadcast_to(np.asarray(eta_slm, float), (n,))
        eta_aod = np.broadcast_to(np.asarray(eta_aod, float), (n,))
        offset = np.full(n, float(delta_drive))
        if delta_site is not None:
            offset = offset + np.asarray(delta_site, float)
        if noise_qs is not None:
            offset = offset + np.asarray(noise_qs, float)
        rabi = (np.ones(n) if rabi_scale is None
                else np.asarray(rabi_scale, float))
        t_clock = np.broadcast_to(np.asarray(t_start, dtype=float),
                                  (n,)).astype(float).copy()
        for lo in range(0, n, chunk):
            hi = min(lo + chunk, n)
            rows_s = pool.u_slm[k_trace[lo:hi]].astype(np.float32)
            rows_a = pool.u_aod[k_trace[lo:hi]].astype(np.float32)
            local = cliff_indices[lo:hi]
            for cliff in np.unique(local):
                sel = np.where(local == cliff)[0]
                gidx = sel + lo
                t_sel = t_clock[gidx].copy()
                for kind, val in self.ops_for(int(cliff)):
                    if kind == "VZ":
                        half = 0.5 * val
                        psi[gidx, 0] *= np.exp(-1j * half)
                        psi[gidx, 1] *= np.exp(1j * half)
                        continue
                    pulse = val
                    duration = pulse.duration
                    n_sub = max(1, int(n_sub_pulse))
                    dt = duration / n_sub
                    env_area = envelope_area_integral(pulse.envelope,
                                                      duration, pulse.edge_us)
                    peak = pulse.nominal_area * (1.0 + pulse.amp_error
                                                 + amp_error) \
                        / max(env_area, 1e-12)
                    for k in range(n_sub):
                        tc = np.mod(t_sel + (k + 0.5) * dt, total)
                        pos = np.clip(tc / sample_dt, 0.0,
                                      rows_s.shape[1] - 1.000001)
                        us = _interp_rows(rows_s[sel], pos[:, None])[:, 0]
                        ua = _interp_rows(rows_a[sel], pos[:, None])[:, 0]
                        delta = ((eta_slm[gidx] * us + eta_aod[gidx] * ua)
                                 / HBAR + offset[gidx])
                        omega = rabi[gidx] * peak * envelope_value(
                            pulse.envelope, (k + 0.5) / n_sub, pulse.edge_us,
                            (k + 0.5) * dt * 1e6, duration * 1e6)
                        u = su2_propagator(delta, omega, pulse.axis_phase, dt)
                        a0 = psi[gidx, 0]
                        b0 = psi[gidx, 1]
                        psi[gidx, 0] = u[0] * a0 + u[1] * b0
                        psi[gidx, 1] = u[2] * a0 + u[3] * b0
                    t_sel = t_sel + duration
                t_clock[gidx] = t_sel
        return psi

    def sequence_duration_s(self, cliff_indices) -> float:
        """一组 Clifford 的物理总时长。"""
        return float(sum(self.clifford_duration_s(int(c))
                         for c in cliff_indices))
