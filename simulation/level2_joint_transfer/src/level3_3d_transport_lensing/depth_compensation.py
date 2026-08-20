"""探索性 lensing 深度补偿（任务书 §6 exploratory 模式，非论文已实现功能）。

思路：命令深度 D_cmd(t) 按数值瞬时最深最低点深度反向放大，
D_cmd(t) = D_base × clip(D_base / d_eff(t), 1, max_power_multiplier)。
补偿只依赖控制路径（轨迹、焦点偏移），与原子状态无关，因此在控制
网格上预计算 D_cmd/Ḋ，积分时线性插值；不逐 shot、逐步调用优化器。

必须诚实报告：所需最大功率倍率、达到上限的饱和时间比例、补偿后的
实际有效阱深。若倍率达到上限，不得声称“保持恒定阱深”。
"""
from __future__ import annotations

import numpy as np

from .aod_lensing import axial_minima
from .gaussian_3d import KB
from .level3_simulation import make_control


def effective_min_depth_j(depth_j, waist_m, z_r_m, zs1, zs2):
    """横向中心处瞬时轴向势最深最低点的深度（正数，J）。"""
    _, depths, _, _ = axial_minima(depth_j, waist_m, z_r_m, zs1, zs2)
    return float(-np.min(depths))


def build_compensation_schedule(physics, profile, geometry, distance_m,
                                duration_s, v_s, signs, max_power_multiplier,
                                n_ctl=321):
    """预计算补偿命令深度表；返回 (schedule dict, 基础 control_at)。

    schedule 含 t_ctl、d_cmd、d_dot、未补偿有效深度 eff_uncomp、补偿后
    实际有效深度 eff_comp 与汇总 report（含饱和比例与功率倍率上限）。
    """
    depth_base = physics["depth_j"]
    base_control = make_control(profile, geometry, distance_m, duration_s, v_s,
                                signs, depth_base, physics["z_r_m"])
    t_ctl = np.linspace(0.0, duration_s, n_ctl)
    eff = np.empty(n_ctl)
    for i, t in enumerate(t_ctl):
        c = base_control(t)
        eff[i] = effective_min_depth_j(depth_base, physics["waist_m"],
                                       physics["z_r_m"], c[4], c[5])
    required = depth_base / eff
    mult = np.clip(required, 1.0, max_power_multiplier)
    d_cmd = depth_base * mult
    d_dot = np.gradient(d_cmd, t_ctl)
    eff_comp = eff * mult
    report = {
        "mode": "exploratory_lensing_compensation",
        "paper_implemented": False,
        "max_required_power_multiplier": float(np.max(required)),
        "max_applied_power_multiplier": float(np.max(mult)),
        "power_cap": float(max_power_multiplier),
        "saturation_time_fraction": float(np.mean(required > max_power_multiplier)),
        "min_effective_depth_uncompensated_uK": float(np.min(eff) / KB * 1e6),
        "min_effective_depth_compensated_uK": float(np.min(eff_comp) / KB * 1e6),
        "depth_base_uK": float(depth_base / KB * 1e6),
        "constant_depth_claim": (
            "maintained" if np.max(required) <= max_power_multiplier
            else "not_maintained_power_cap_reached"),
        "n_control_points": int(n_ctl),
    }
    schedule = {"t_ctl": t_ctl, "d_cmd": d_cmd, "d_dot": d_dot,
                "eff_uncomp": eff, "eff_comp": eff_comp, "mult": mult,
                "report": report}
    return schedule, base_control


def make_compensated_control(base_control, schedule):
    """把基础控制量的 (D_cmd, Ḋ) 替换为补偿表插值，其余分量不变。"""

    def control_at(t):
        c = base_control(t)
        d_cmd = float(np.interp(t, schedule["t_ctl"], schedule["d_cmd"]))
        d_dot = float(np.interp(t, schedule["t_ctl"], schedule["d_dot"]))
        return (*c[:8], d_cmd, d_dot)

    return control_at
