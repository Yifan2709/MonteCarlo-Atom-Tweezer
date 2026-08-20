"""半经典差分光移相位累积、toggling 积分与 ensemble 对比度。

只追踪每条经典轨迹上两个超精细钟态之间的相对相位：
φ_i = ∫ y(t)·[η_SLM·U_SLM(r_i(t),t) + η_AOD·U_AOD(r_i(t),t)]/ħ dt。
积分把 A_SLM=∫y·U_SLM dt 与 A_AOD=∫y·U_AOD dt 分开累积，因此对任意
(η_SLM, η_AOD) 线性可重算（η 敏感性与 no_differential_light_shift 免重跑）。
"""
from __future__ import annotations

import numpy as np

HBAR = 1.054571817e-34


class PhaseAccumulator:
    """按步累积 A_SLM/A_AOD（含 per-segment 分解），精确处理 toggling 跳变。

    每步用与梯形功积分一致的线性插值：步内脉冲把步拆成两段解析积分。
    """

    def __init__(self, modes: dict[str, np.ndarray], times: np.ndarray,
                 n_shots: int, segment_index: np.ndarray):
        self.times = times
        self.n_shots = n_shots
        self.segment_index = segment_index
        self.n_segments = int(segment_index.max()) + 1
        self.state = {}
        for mode, pulses in modes.items():
            steps = _pulses_to_steps(np.asarray(pulses, dtype=float), times)
            self.state[mode] = {
                "pulses_in_steps": steps,
                "y": 1.0,
                "next_pulse": 0,
                "a_slm": np.zeros(n_shots),
                "a_aod": np.zeros(n_shots),
                "a_slm_seg": np.zeros((self.n_segments, n_shots)),
                "a_aod_seg": np.zeros((self.n_segments, n_shots)),
            }

    def snapshot(self):
        """当前各模式的累积量副本（轮次边界保存用）。"""
        return {mode: {"a_slm": st["a_slm"].copy(), "a_aod": st["a_aod"].copy(),
                       "y": st["y"], "next_pulse": st["next_pulse"]}
                for mode, st in self.state.items()}

    def restore(self, snap):
        """恢复快照（多轮重算/回退用）。"""
        for mode, st in self.state.items():
            st["a_slm"] = snap[mode]["a_slm"].copy()
            st["a_aod"] = snap[mode]["a_aod"].copy()
            st["y"] = snap[mode]["y"]
            st["next_pulse"] = snap[mode]["next_pulse"]

    def step(self, i, u_slm_i, u_aod_i, u_slm_ip1, u_aod_ip1):
        """推进一个时间步 [t_i, t_{i+1}] 的相位积分。"""
        dt = self.times[i + 1] - self.times[i]
        seg = self.segment_index[i]
        for mode, st in self.state.items():
            total_slm = 0.5 * (u_slm_i + u_slm_ip1) * dt
            total_aod = 0.5 * (u_aod_i + u_aod_ip1) * dt
            extra = [(s, f) for (s, f) in st["pulses_in_steps"] if s == i]
            contrib_slm = np.zeros_like(u_slm_i)
            contrib_aod = np.zeros_like(u_aod_i)
            if not extra:
                contrib_slm += st["y"] * total_slm
                contrib_aod += st["y"] * total_aod
            else:
                y = st["y"]
                left_slm = np.zeros_like(u_slm_i)
                left_aod = np.zeros_like(u_aod_i)
                for _, frac in extra:
                    part_slm = (u_slm_i * frac
                                + 0.5 * (u_slm_ip1 - u_slm_i) * frac * frac) * dt
                    part_aod = (u_aod_i * frac
                                + 0.5 * (u_aod_ip1 - u_aod_i) * frac * frac) * dt
                    contrib_slm += y * (part_slm - left_slm)
                    contrib_aod += y * (part_aod - left_aod)
                    left_slm, left_aod = part_slm, part_aod
                    y = -y
                contrib_slm += y * (total_slm - left_slm)
                contrib_aod += y * (total_aod - left_aod)
                st["y"] = y
                st["next_pulse"] += len(extra)
            st["a_slm"] += contrib_slm
            st["a_aod"] += contrib_aod
            st["a_slm_seg"][seg] += contrib_slm
            st["a_aod_seg"][seg] += contrib_aod

    def phases(self, eta_slm, eta_aod):
        """按 (η_SLM, η_AOD) 计算各模式的相位 φ（rad）。"""
        return {mode: (eta_slm * st["a_slm"] + eta_aod * st["a_aod"]) / HBAR
                for mode, st in self.state.items()}

    def segment_phases(self, eta_slm, eta_aod):
        """各模式按分段的相位矩阵 (n_segments, n_shots)。"""
        return {mode: (eta_slm * st["a_slm_seg"] + eta_aod * st["a_aod_seg"]) / HBAR
                for mode, st in self.state.items()}


def _pulses_to_steps(pulses, times):
    """脉冲时刻 → (步索引, 步内分数)，要求脉冲严格位于步内部。"""
    out = []
    for p in pulses:
        idx = int(np.searchsorted(times, p, side="left")) - 1
        if idx < 0 or idx >= len(times) - 1:
            raise ValueError(f"脉冲 {p:.6e} s 落在时间网格之外")
        if not (times[idx] < p < times[idx + 1]):
            raise ValueError(f"脉冲 {p:.6e} s 恰好落在时间样本上")
        out.append((idx, (p - times[idx]) / (times[idx + 1] - times[idx])))
    return out


def differential_shift(u_slm_j, u_aod_j, eta_slm, eta_aod):
    """差分光移角频率 δω = (η_SLM·U_SLM + η_AOD·U_AOD)/ħ。"""
    return (eta_slm * u_slm_j + eta_aod * u_aod_j) / HBAR


def conditional_contrast(phases_rad, mask=None):
    """条件相干对比度 C=|mean(exp(iφ))| 与圆统计量。"""
    phi = np.asarray(phases_rad, dtype=float)
    if mask is not None:
        phi = phi[np.asarray(mask, dtype=bool)]
    if phi.size == 0:
        return {"n": 0, "contrast": None, "mean_direction": None,
                "circular_std": None,
                "phase_mean": None, "phase_std": None}
    z = np.exp(1j * phi)
    mean_z = z.mean()
    contrast = float(np.abs(mean_z))
    return {
        "n": int(phi.size),
        "contrast": contrast,
        "mean_direction": float(np.angle(mean_z)),
        "circular_std": float(np.sqrt(max(0.0, -2.0 * np.log(contrast))))
        if contrast > 0 else None,
        "phase_mean": float(np.mean(phi)),
        "phase_std": float(np.std(phi)),
    }


def contrast_bootstrap_ci(phases_rad, mask, n_boot=2000, seed=0):
    """对比度的 bootstrap 95% 区间（success 条件下重采样）。"""
    rng = np.random.default_rng(seed)
    phi = np.asarray(phases_rad, dtype=float)[np.asarray(mask, dtype=bool)]
    if phi.size == 0:
        return {"ci_low": None, "ci_high": None}
    idx = rng.integers(0, phi.size, size=(int(n_boot), phi.size))
    z = np.exp(1j * phi)[idx].mean(axis=1)
    c = np.abs(z)
    return {"ci_low": float(np.percentile(c, 2.5)),
            "ci_high": float(np.percentile(c, 97.5)),
            "n_bootstrap": int(n_boot), "seed": int(seed)}


def analytic_constant_phase(delta_omega, duration_s, pulse_times=()):
    """常数 δω 的解析相位：按 toggling 符号逐段累积。"""
    bounds = np.concatenate(
        ([0.0], np.sort(np.asarray(pulse_times, dtype=float)), [duration_s]))
    phi = 0.0
    y = 1.0
    for a, b in zip(bounds[:-1], bounds[1:]):
        phi += y * delta_omega * (b - a)
        y = -y
    return phi
