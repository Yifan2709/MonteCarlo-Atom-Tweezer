"""Level 4 初态采样：复用 Level 3 三维简谐提议 + 束缚拒绝方法。"""
from __future__ import annotations

import numpy as np

from level3_3d_transport_lensing.gaussian_3d import KB
from level3_3d_transport_lensing.thermal_sampling_3d import \
    sample_initial_states_3d
from .protocol_segments import physics_from_config


def sample_slm_states(cfg, seed, shots, temperature_uK=None,
                      depth_uK=None):
    """在 SLM 静态阱中按温度采样六维初态（含接受率等元数据）。"""
    physics = physics_from_config(cfg)
    temp = cfg.temperature_uK if temperature_uK is None else temperature_uK
    depth_j = (cfg.slm_initial_depth_uK if depth_uK is None
               else depth_uK) * 1e-6 * KB
    omega_r = float(np.sqrt(4.0 * depth_j
                            / (physics.mass_kg * physics.slm_waist_m ** 2)))
    omega_z = float(np.sqrt(2.0 * depth_j
                            / (physics.mass_kg * physics.slm_zr_m ** 2)))
    states = sample_initial_states_3d(
        seed, shots, temp, physics.mass_kg, depth_j, physics.slm_waist_m,
        physics.slm_zr_m, omega_r, omega_z)
    states["pos_m"] = states["pos_m"] + np.array(
        [physics.slm_center_m[0], physics.slm_center_m[1], 0.0])
    states["temperature_uK"] = float(temp)
    states["trap_depth_uK"] = depth_j / KB / 1e-6
    states["sampler_note"] = (
        "Level 3 三维简谐提议 + 完整势 E<0 束缚拒绝；非有限深高斯阱的严格"
        "正则系综（高斯尾部截断 + 简谐速度提议引入轻微偏差）。")
    return states
