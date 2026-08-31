"""显式时变势中的外部功核算与功-能残差检查。"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from level0_static_trap.potentials import gaussian_potential


def work_energy_along_trajectory(time_s, position_m, velocity_m_per_s, waveform, mass, slm_depth_j,
                                  slm_waist_m, slm_center_m, aod_waist_m, noise_work_j=None):
    """沿轨迹计算深度功、移动功、总功与功-能残差 R_W(t)。

    dU_AOD/dt = -Ddot*exp(-2q^2/w^2) + F_AOD*cdot，其中 q = x - c(t)。
    noise_work_j 为可选的确定性噪声注入能量累计数组（J，与轨迹同长）；
    提供时残差额外扣除该通道，恒等式为 ΔE = W_ext + W_noise + R_W。
    """
    t = np.asarray(time_s, dtype=float)
    x = np.asarray(position_m, dtype=float)
    v = np.asarray(velocity_m_per_s, dtype=float)

    center = np.asarray(waveform.center(t), dtype=float)
    depth = np.asarray(waveform.depth(t), dtype=float)
    center_rate = np.asarray(waveform.center_velocity(t), dtype=float)
    depth_rate = np.asarray(waveform.depth_rate(t), dtype=float)

    offset = x - center
    envelope = np.exp(-2.0 * (offset / aod_waist_m) ** 2)
    depth_power = -depth_rate * envelope                       # 深度变化功功率
    aod_force = -(4.0 * depth * offset / aod_waist_m**2) * envelope
    move_power = aod_force * center_rate                       # 中心移动功功率

    work_depth = cumulative_trapezoid(depth_power, t, initial=0.0)
    work_move = cumulative_trapezoid(move_power, t, initial=0.0)
    work_total = work_depth + work_move

    u_slm = gaussian_potential(x, slm_depth_j, slm_waist_m, slm_center_m)
    u_aod = -depth * envelope
    energy_total = 0.5 * mass * v**2 + u_slm + u_aod
    residual = energy_total - energy_total[0] - work_total
    if noise_work_j is not None:
        noise_work = np.asarray(noise_work_j, dtype=float)
        if noise_work.shape != t.shape:
            raise ValueError("noise_work_j 必须与轨迹时间轴同长")
        residual = residual - noise_work
    else:
        noise_work = np.zeros_like(t)

    return {
        "time_s": t,
        "work_depth_J": work_depth,
        "work_move_J": work_move,
        "work_total_J": work_total,
        "work_noise_J": noise_work,
        "energy_total_J": energy_total,
        "residual_J": residual,
        "final_work_depth_J": float(work_depth[-1]),
        "final_work_move_J": float(work_move[-1]),
        "final_work_total_J": float(work_total[-1]),
        "final_work_noise_J": float(noise_work[-1]),
        "max_abs_residual_J": float(np.max(np.abs(residual))),
    }


def hold_segment_energy_error(time_s, position_m, velocity_m_per_s, mass, slm_depth_j,
                              slm_waist_m, slm_center_m, transfer_end_s):
    """AOD 关闭后的纯 SLM 保持段能量误差幅度（应有界，无 secular drift）。"""
    t = np.asarray(time_s, dtype=float)
    x = np.asarray(position_m, dtype=float)
    v = np.asarray(velocity_m_per_s, dtype=float)
    mask = t >= transfer_end_s - 1e-15
    if not np.any(mask) or mask.sum() < 2:
        return {"hold_samples": int(mask.sum()), "peak_to_peak_energy_error_J": None,
                "peak_to_peak_energy_error_over_slm_depth": None}
    energy = 0.5 * mass * v[mask] ** 2 + gaussian_potential(x[mask], slm_depth_j, slm_waist_m, slm_center_m)
    span = float(np.ptp(energy))
    return {
        "hold_samples": int(mask.sum()),
        "peak_to_peak_energy_error_J": span,
        "peak_to_peak_energy_error_over_slm_depth": span / slm_depth_j,
    }


def work_energy_two_trap(time_s, position_m, velocity_m_per_s, mass,
                         aod_depth_j, aod_center_m, aod_depth_rate_j_s, aod_center_rate_m_s,
                         slm_depth_j, slm_waist_m, slm_factor, aod_waist_m):
    """双阱均可时变（含 SLM 开关跳变）的功-能核算（Level 2C 往返序列用）。

    外部功 = AOD 深度功 + AOD 移动功 + SLM 开关跳变功：
    - AOD 两项沿用单阱公式 −Ḋ·exp(−2q²/w²) 与 F_AOD·ċ（数组由调用方
      按解析导数给出，与积分轨迹同一时间轴）；
    - SLM 因子 slm_factor∈[0,1] 的每个跳变（含瞬时开关）贡献
      ΔW = U_SLM(x)·Δfactor = −D_SLM·exp(−2x²/w²)·Δfactor，
      在跳变步索引处按该步位置计入（梯形积分的分布极限）。
    返回功的分量、总功、总能量与残差 R_W(t)=ΔE−W_ext 数组。
    """
    t = np.asarray(time_s, dtype=float)
    x = np.asarray(position_m, dtype=float)
    v = np.asarray(velocity_m_per_s, dtype=float)
    depth = np.asarray(aod_depth_j, dtype=float)
    center = np.asarray(aod_center_m, dtype=float)
    depth_rate = np.asarray(aod_depth_rate_j_s, dtype=float)
    center_rate = np.asarray(aod_center_rate_m_s, dtype=float)
    factor = np.asarray(slm_factor, dtype=float)

    offset = x - center
    envelope = np.exp(-2.0 * (offset / aod_waist_m) ** 2)
    depth_power = -depth_rate * envelope
    aod_force = -(4.0 * depth * offset / aod_waist_m**2) * envelope
    move_power = aod_force * center_rate
    work_depth = cumulative_trapezoid(depth_power, t, initial=0.0)
    work_move = cumulative_trapezoid(move_power, t, initial=0.0)

    u_slm_static = gaussian_potential(x, slm_depth_j, slm_waist_m, 0.0)
    work_slm = np.zeros_like(t)
    jumps = np.diff(factor)
    jump_idx = np.nonzero(jumps)[0]
    for i in jump_idx:
        work_slm[i + 1:] += u_slm_static[i + 1] * jumps[i]

    work_total = work_depth + work_move + work_slm
    energy_total = 0.5 * mass * v**2 + u_slm_static * factor - depth * envelope
    residual = energy_total - energy_total[0] - work_total
    return {
        "time_s": t,
        "work_depth_J": work_depth,
        "work_move_J": work_move,
        "work_slm_switch_J": work_slm,
        "work_total_J": work_total,
        "energy_total_J": energy_total,
        "residual_J": residual,
        "max_abs_residual_J": float(np.max(np.abs(residual))),
        "final_work_total_J": float(work_total[-1]),
        "slm_switch_count": int(len(jump_idx)),
    }
