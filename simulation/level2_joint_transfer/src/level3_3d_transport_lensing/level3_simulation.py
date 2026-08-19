"""Level 3 仿真编排：控制量生成、单条件运行、确定性扫描、粗 MC、验证与敏感性。"""
from __future__ import annotations

import numpy as np

from .aod_lensing import (axial_minima, axis_shift_m, vs_from_ratio,
                          reference_axis_vmax)
from .gaussian_3d import KB, physics_defaults, static_potential
from .integrators_3d import run_trajectory, static_hold
from .level3_config import Level3Config
from .thermal_sampling_3d import sample_initial_states_3d
from .transport_profiles import geometry_centers, geometry_velocity_factors, get_profile
from .transport_statistics import (paired_bootstrap_difference, paired_counts,
                                   summarize_retention)


def physics_from_config(cfg: Level3Config, depth_uK=None, waist_um=None) -> dict:
    """由配置生成物理字典（可覆盖深度/束腰用于敏感性）。"""
    depth = cfg.depth_uK if depth_uK is None else depth_uK
    waist = cfg.waist_um if waist_um is None else waist_um
    physics = physics_defaults(depth_j=depth * 1e-6 * KB, waist_m=waist * 1e-6,
                               wavelength_m=cfg.wavelength_nm * 1e-9)
    physics["depth_uK"] = depth
    physics["waist_um"] = waist
    return physics


def resolve_vs(cfg: Level3Config, severity: float):
    """按校准模式返回固定特征速度 v_s（severity=0 → None 表示关闭透镜）。"""
    if severity == 0:
        return None
    if cfg.calibration_mode == "absolute_critical_velocity":
        return cfg.absolute_critical_velocity_m_per_s
    return vs_from_ratio(severity, cfg.reference_profile,
                         cfg.reference_distance_um * 1e-6,
                         cfg.reference_duration_us * 1e-6, cfg.reference_geometry)


def make_control(profile, geometry, distance_m, duration_s, v_s, signs, depth_j, z_r_m):
    """返回 control_at(t) -> (c1,c2,ċ1,ċ2,zs1,zs2,żs1,żs2,D_cmd,Ḋ)。"""
    prof = get_profile(profile)
    f1, f2 = geometry_velocity_factors(geometry)

    def control_at(t):
        s = min(max(t / duration_s, 0.0), 1.0)
        q = float(prof.q(s))
        qd = float(prof.q_dot(s))
        qdd = float(prof.q_ddot(s))
        c1, c2 = geometry_centers(q, geometry, distance_m)
        c1d = distance_m * qd * f1 / duration_s
        c2d = distance_m * qd * f2 / duration_s
        c1dd = distance_m * qdd * f1 / duration_s ** 2
        c2dd = distance_m * qdd * f2 / duration_s ** 2
        if v_s:
            zs1 = float(axis_shift_m(c1d, v_s, signs[0], z_r_m))
            zs2 = float(axis_shift_m(c2d, v_s, signs[1], z_r_m))
            zs1d = float(axis_shift_m(c1dd, v_s, signs[0], z_r_m))
            zs2d = float(axis_shift_m(c2dd, v_s, signs[1], z_r_m))
        else:
            zs1 = zs2 = zs1d = zs2d = 0.0
        return (c1, c2, c1d, c2d, zs1, zs2, zs1d, zs2d, depth_j, 0.0)

    return control_at


def split_fraction(profile, geometry, distance_m, duration_s, v_s, signs, z_r_m,
                   depth_j, waist_m, samples=201):
    """运输期间瞬时轴向势处于双阱状态的时间比例（数值极值统计）。"""
    prof = get_profile(profile)
    f1, f2 = geometry_velocity_factors(geometry)
    counts = 0
    for s in np.linspace(0, 1, samples):
        qd = float(prof.q_dot(s))
        c1d = distance_m * qd * f1 / duration_s
        c2d = distance_m * qd * f2 / duration_s
        zs1 = float(axis_shift_m(c1d, v_s, signs[0], z_r_m)) if v_s else 0.0
        zs2 = float(axis_shift_m(c2d, v_s, signs[1], z_r_m)) if v_s else 0.0
        _, _, _, n = axial_minima(depth_j, waist_m, z_r_m, zs1, zs2)
        counts += (n >= 2)
    return counts / samples


def run_condition(physics, states, profile, geometry, distance_um, duration_us,
                  severity, v_s, signs, dt_us, hold_us, cfg: Level3Config,
                  keep_trajectory_ids=(), trajectory_stride=200):
    """运行一个运输条件的一批初态，返回逐 shot 记录与摘要。"""
    mass = physics["mass_kg"]
    depth_j = physics["depth_j"]
    waist_m = physics["waist_m"]
    z_r = physics["z_r_m"]
    distance_m = distance_um * 1e-6
    duration_s = duration_us * 1e-6
    control_at = make_control(profile, geometry, distance_m, duration_s, v_s, signs,
                              depth_j, z_r)
    n_transport = int(round(duration_s / (dt_us * 1e-6)))
    t_grid = np.arange(n_transport + 1) * dt_us * 1e-6
    pos0 = states["pos_m"].copy()
    vel0 = states["vel_m_per_s"].copy()
    traj_ids = tuple(sorted(keep_trajectory_ids))
    subset_pos = pos0[list(traj_ids)] if traj_ids else pos0[:0]
    subset_vel = vel0[list(traj_ids)] if traj_ids else vel0[:0]
    result = run_trajectory(pos0, vel0, t_grid, control_at, mass, depth_j, waist_m,
                            z_r)
    sub_result = run_trajectory(subset_pos, subset_vel, t_grid, control_at, mass,
                                depth_j, waist_m, z_r,
                                record_stride=trajectory_stride) if traj_ids else None

    # 保持段在终点阱（中心 c_end、无透镜）中静止：先把坐标平移到阱心
    c_end = geometry_centers(1.0, geometry, distance_m)
    hold_pos_rel = result["final_pos"] - np.array([c_end[0], c_end[1], 0.0])
    _, _, hold_pp, hold_e_final = static_hold(
        hold_pos_rel, result["final_vel"], hold_us * 1e-6, dt_us * 1e-6,
        mass, depth_j, waist_m, z_r)
    sf = split_fraction(profile, geometry, distance_m, duration_s, v_s, signs, z_r,
                        depth_j, waist_m)
    depth_uK = physics["depth_uK"]
    near_frac = cfg.near_threshold_fraction_of_depth

    records = []
    eps_i_uK = (states["initial_energy_j"] + depth_j) / KB * 1e6
    for i in range(len(pos0)):
        e_final = result["e_final"][i]
        e_hold = hold_e_final[i]
        retained = e_final < cfg.retained_energy_threshold_J
        retained_hold = e_hold < cfg.retained_energy_threshold_J
        eps_f_uK = (e_final + depth_j) / KB * 1e6
        mean_z = result["mean_z"][i]
        if sf > 0 and abs(mean_z) < 0.1 * z_r:
            branch = "ambiguous"
        else:
            branch = "+" if mean_z >= 0 else "-"
        records.append({
            "shot_id": i, "profile": profile, "geometry": geometry,
            "distance_um": distance_um, "duration_us": duration_us,
            "severity": severity,
            "x0_um": pos0[i, 0] * 1e6, "y0_um": pos0[i, 1] * 1e6,
            "z0_um": pos0[i, 2] * 1e6,
            "vx0_m_s": vel0[i, 0], "vy0_m_s": vel0[i, 1], "vz0_m_s": vel0[i, 2],
            "xf_um": result["final_pos"][i, 0] * 1e6,
            "yf_um": result["final_pos"][i, 1] * 1e6,
            "zf_um": result["final_pos"][i, 2] * 1e6,
            "vxf_m_s": result["final_vel"][i, 0], "vyf_m_s": result["final_vel"][i, 1],
            "vzf_m_s": result["final_vel"][i, 2],
            "e0_uK": result["e0"][i] / KB * 1e6, "e_final_uK": e_final / KB * 1e6,
            "eps_i_uK": eps_i_uK[i], "eps_f_uK": eps_f_uK,
            "delta_eps_uK": eps_f_uK - eps_i_uK[i],
            "retained": bool(retained), "retained_after_hold": bool(retained_hold),
            "near_threshold": bool(abs(e_final) / depth_j < near_frac),
            "max_radial_um": result["max_radial"][i] * 1e6,
            "max_axial_um": result["max_axial"][i] * 1e6,
            "max_speed_m_s": result["max_speed"],
            "split_exposed": bool(sf > 0), "split_fraction": sf,
            "branch": branch, "mean_z_um": mean_z * 1e6,
            "w_total_uK": result["w_total"][i] / KB * 1e6,
            "w_move_uK": result["w_move"][i] / KB * 1e6,
            "w_shift_uK": result["w_shift"][i] / KB * 1e6,
            "w_residual_uK": result["residual"][i] / KB * 1e6,
            "hold_peak_to_peak_uK": hold_pp[i] / KB * 1e6,
        })
    retained_flags = [r["retained"] for r in records]
    summary = summarize_retention(retained_flags, f"{profile}|{geometry}|{severity}")
    summary["split_fraction"] = sf
    summary["median_delta_eps_uK_retained"] = float(np.median(
        [r["delta_eps_uK"] for r in records if r["retained"]])) if any(
        retained_flags) else None
    summary["max_abs_w_residual_over_depth"] = float(np.max(
        np.abs(result["residual"]))) / depth_j
    summary["max_hold_pp_over_depth"] = float(np.max(hold_pp)) / depth_j
    summary["retained_mismatch_with_hold"] = int(sum(
        r["retained"] != r["retained_after_hold"] for r in records))
    out = {"records": records, "summary": summary}
    if sub_result is not None:
        out["trajectories"] = sub_result["trajectory"]
    return out


def sample_states(cfg: Level3Config, physics, seed, shots, temperature_uK=None):
    """按当前物理参数抽样三维热初态。"""
    temp = cfg.temperature_uK if temperature_uK is None else temperature_uK
    return sample_initial_states_3d(
        seed, shots, temp, physics["mass_kg"], physics["depth_j"], physics["waist_m"],
        physics["z_r_m"], physics["omega_r"], physics["omega_z"])


def deterministic_scan(cfg: Level3Config, physics):
    """阶段 A：无 MC 的势阱与轨迹确定性扫描，返回两个行列表。"""
    minima_rows, traj_rows = [], []
    signs = cfg.axis_signs
    z_r = physics["z_r_m"]
    severities = cfg.reference_axis_vmax_over_vs
    for profile in cfg.profiles:
        for geometry in cfg.geometries:
            for distance in cfg.distances_um:
                for duration in cfg.durations_us:
                    for severity in severities:
                        v_s = resolve_vs(cfg, severity)
                        prof = get_profile(profile)
                        f1, f2 = geometry_velocity_factors(geometry)
                        s_grid = np.linspace(0, 1, 41)
                        path_vmax = distance * 1e-6 * prof.q_dot_max / (duration * 1e-6)
                        axis_vmax1 = path_vmax * f1
                        axis_vmax2 = path_vmax * f2
                        zs1_peak = float(axis_shift_m(axis_vmax1, v_s, signs[0], z_r)) \
                            if v_s else 0.0
                        zs2_peak = float(axis_shift_m(axis_vmax2, v_s, signs[1], z_r)) \
                            if v_s else 0.0
                        sf = split_fraction(profile, geometry, distance * 1e-6,
                                             duration * 1e-6, v_s, signs, z_r,
                                             physics["depth_j"], physics["waist_m"])
                        zmin, dmin, barrier, n_min = axial_minima(
                            physics["depth_j"], physics["waist_m"], z_r,
                            zs1_peak, zs2_peak)
                        minima_rows.append({
                            "profile": profile, "geometry": geometry,
                            "distance_um": distance, "duration_us": duration,
                            "severity": severity,
                            "path_vmax_m_s": path_vmax,
                            "axis1_vmax_m_s": axis_vmax1, "axis2_vmax_m_s": axis_vmax2,
                            "v1_over_vs": severity * f1 if severity else 0.0,
                            "v2_over_vs": severity * f2 if severity else 0.0,
                            "zs1_peak_over_zr": zs1_peak / z_r,
                            "zs2_peak_over_zr": zs2_peak / z_r,
                            "split_fraction": sf, "n_minima_at_peak": n_min,
                            "min_depth_uK": float(np.max(dmin)) / KB * 1e6,
                            "barrier_uK": barrier / KB * 1e6,
                        })
                        s_dense = np.linspace(0, 1, 201)
                        qd = prof.q_dot(s_dense)
                        qdd = prof.q_ddot(s_dense)
                        for k, s in enumerate(s_dense):
                            traj_rows.append({
                                "profile": profile, "geometry": geometry,
                                "distance_um": distance, "duration_us": duration,
                                "severity": severity,
                                "time_us": s * duration,
                                "path_fraction": float(prof.q(s)),
                                "q_dot": float(qd[k]), "q_ddot": float(qdd[k]),
                                "axis_speed_m_s": distance * 1e-6 * float(qd[k]) * f1
                                / (duration * 1e-6),
                            })
    return minima_rows, traj_rows


def compare_paired(records_a, records_b, seed, n_boot=2000):
    """同初态配对比较两个条件的留阱结果。"""
    a = np.array([r["retained"] for r in records_a])
    b = np.array([r["retained"] for r in records_b])
    return {
        "paired_counts": paired_counts(a, b),
        "bootstrap": paired_bootstrap_difference(a, b, n_boot, seed),
    }
