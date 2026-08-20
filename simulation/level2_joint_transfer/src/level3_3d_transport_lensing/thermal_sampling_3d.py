"""三维热初态采样：简谐提议分布 + 完整三维高斯势束缚拒绝。"""
from __future__ import annotations

import numpy as np

from .gaussian_3d import KB, static_potential


def sample_initial_states_3d(seed, shots, temperature_uK, mass_kg, depth_j, waist_m,
                             z_r_m, omega_r, omega_z, max_attempts=200):
    """简谐近似提议 + E_i<0 拒绝采样，返回六维初态与统计元数据。"""
    rng = np.random.default_rng(seed)
    sigma_r = np.sqrt(KB * temperature_uK * 1e-6 / (mass_kg * omega_r ** 2))
    sigma_z = np.sqrt(KB * temperature_uK * 1e-6 / (mass_kg * omega_z ** 2))
    v_sigma = np.sqrt(KB * temperature_uK * 1e-6 / mass_kg)
    accepted_pos, accepted_vel = [], []
    attempts = 0
    accepted_count = 0
    while accepted_count < shots and attempts < max_attempts * shots:
        need = shots - len(accepted_pos)
        batch = max(need * 2, 64)
        pos = np.stack([rng.normal(0.0, sigma_r, batch),
                        rng.normal(0.0, sigma_r, batch),
                        rng.normal(0.0, sigma_z, batch)], axis=1)
        vel = np.stack([rng.normal(0.0, v_sigma, batch),
                        rng.normal(0.0, v_sigma, batch),
                        rng.normal(0.0, v_sigma, batch)], axis=1)
        energy = (0.5 * mass_kg * np.sum(vel ** 2, axis=1)
                  + static_potential(pos[:, 0], pos[:, 1], pos[:, 2],
                                     depth_j, waist_m, z_r_m))
        keep = energy < 0
        accepted_pos.append(pos[keep])
        accepted_vel.append(vel[keep])
        accepted_count += int(keep.sum())
        attempts += batch
    pos = np.concatenate(accepted_pos)[:shots]
    vel = np.concatenate(accepted_vel)[:shots]
    if len(pos) < shots:
        raise RuntimeError("三维初态拒绝采样未能在预算内凑满样本")
    energy = (0.5 * mass_kg * np.sum(vel ** 2, axis=1)
              + static_potential(pos[:, 0], pos[:, 1], pos[:, 2],
                                 depth_j, waist_m, z_r_m))
    states = {
        "seed": seed, "shots": len(pos), "attempts": attempts,
        # 物理接受率 = 被接受提议数 / 总提议数。批按 2 倍需求超采样，
        # 首批即凑满并截断，被丢弃的多余接受样本不影响该比值。
        "acceptance_rate": accepted_count / attempts,
        "surplus_accepted_discarded": int(accepted_count - len(pos)),
        "sampling_note": (
            "简谐提议 + E_i<0 束缚拒绝；D/kT=56 时几乎无拒绝，物理接受率≈1。"
            "attempts 含 2 倍批超采样；这不是完整有限深高斯阱的三维正则系综。"),
        "pos_m": pos, "vel_m_per_s": vel, "initial_energy_j": energy,
        "sigma_r_m": float(sigma_r), "sigma_z_m": float(sigma_z),
        "mean_pos_m": pos.mean(axis=0).tolist(), "std_pos_m": pos.std(axis=0).tolist(),
        "mean_vel_m_per_s": vel.mean(axis=0).tolist(),
        "std_vel_m_per_s": vel.std(axis=0).tolist(),
        "sampler": "harmonic_3d_proposal_with_bound_rejection",
    }
    return states
