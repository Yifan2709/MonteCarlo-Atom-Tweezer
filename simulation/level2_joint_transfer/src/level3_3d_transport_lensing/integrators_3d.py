"""向量化三维时变 velocity-Verlet 推进与显式功-能核算。"""
from __future__ import annotations

import numpy as np

from .aod_lensing import lensing_potential_and_force, partial_t_potential

CONTROL_ZERO = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None, 0.0)


def total_energy(pos, vel, mass_kg, depth_j, waist_m, z_r_m, c1, c2, zs1, zs2):
    """E = ½m|v|² + U（透镜势，zs=0 时即标准三维高斯势），坐标相对阱中心。"""
    u, _, _, _ = lensing_potential_and_force(pos[:, 0] - c1, pos[:, 1] - c2, pos[:, 2],
                                             depth_j, waist_m, z_r_m, zs1, zs2)
    return 0.5 * mass_kg * np.sum(vel ** 2, axis=1) + u


def make_force_fn(depth_j, waist_m, z_r_m, control_at, mass_kg):
    """control_at(t) -> (c1,c2,ċ1,ċ2,zs1,zs2,żs1,żs2,D_cmd,Ḋ)；返回 (accel, power)。

    accel(pos, t) 返回 (N,3) 加速度（力已除以质量）；power(pos, t) 返回功-能分解。
    """

    def accel(pos, t):
        c = control_at(t)
        _, f_xi, f_eta, f_z = lensing_potential_and_force(
            pos[:, 0] - c[0], pos[:, 1] - c[1], pos[:, 2], c[8], waist_m, z_r_m,
            c[4], c[5])
        return np.stack([f_xi, f_eta, f_z], axis=1) / mass_kg

    def power(pos, t):
        c = control_at(t)
        xi = pos[:, 0] - c[0]
        eta = pos[:, 1] - c[1]
        u, f_xi, f_eta, _ = lensing_potential_and_force(
            xi, eta, pos[:, 2], c[8], waist_m, z_r_m, c[4], c[5])
        return partial_t_potential(xi, eta, pos[:, 2], c[8], waist_m, z_r_m,
                                   c[4], c[5], c[2], c[3], c[6], c[7], c[9],
                                   f_xi=f_xi, f_eta=f_eta, g1g2=-u / c[8])

    return accel, power


def run_trajectory(pos0, vel0, t_grid, control_at, mass_kg, depth_j, waist_m, z_r_m,
                   record_stride=0):
    """推进一批初态（velocity-Verlet），返回末态与功-能诊断。

    功率用步中点控制量（二阶精度）；record_stride>0 时按间隔保存轨迹。
    """
    accel, power = make_force_fn(depth_j, waist_m, z_r_m, control_at, mass_kg)
    pos = pos0.astype(float).copy()
    vel = vel0.astype(float).copy()
    acc = accel(pos, t_grid[0])

    def energy(t):
        c = control_at(t)
        return total_energy(pos, vel, mass_kg, depth_j, waist_m, z_r_m,
                            c[0], c[1], c[4], c[5])

    n_steps = len(t_grid) - 1
    e0 = energy(t_grid[0])
    w_total = np.zeros(len(pos))
    w_move = np.zeros(len(pos))
    w_shift = np.zeros(len(pos))
    w_depth = np.zeros(len(pos))
    max_speed = float(np.max(np.linalg.norm(vel, axis=1)))
    max_radial = np.zeros(len(pos))   # 相对瞬时阱中心的横向位移 sqrt(ξ²+η²)
    max_axial = np.zeros(len(pos))    # 相对两柱面焦点中点的轴向位移 |z−(zs1+zs2)/2|
    times_rec, pos_rec, vel_rec = [], [], []
    z_sum = np.zeros(len(pos))        # 用于运输期间 z 的时间平均（分支诊断）
    p_old = power(pos, t_grid[0])
    for i in range(n_steps):
        t = t_grid[i]
        dt = t_grid[i + 1] - t
        pos = pos + vel * dt + 0.5 * acc * dt ** 2
        acc_new = accel(pos, t_grid[i + 1])
        vel = vel + 0.5 * (acc + acc_new) * dt
        acc = acc_new
        p_new = power(pos, t_grid[i + 1])
        # 梯形法则积分外部功（一阶中点/左端点求积会引入 O(dt) 全局残差）
        w_total += 0.5 * (p_old[0] + p_new[0]) * dt
        w_move += 0.5 * (p_old[1] + p_new[1]) * dt
        w_shift += 0.5 * (p_old[2] + p_new[2]) * dt
        w_depth += 0.5 * (p_old[3] + p_new[3]) * dt
        p_old = p_new
        max_speed = max(max_speed, float(np.max(np.linalg.norm(vel, axis=1))))
        c = control_at(t_grid[i + 1])
        radial = np.hypot(pos[:, 0] - c[0], pos[:, 1] - c[1])
        axial = np.abs(pos[:, 2] - 0.5 * (c[4] + c[5]))
        np.maximum(max_radial, radial, out=max_radial)
        np.maximum(max_axial, axial, out=max_axial)
        z_sum += pos[:, 2] * dt
        if record_stride and (i % record_stride == 0):
            times_rec.append(t_grid[i + 1])
            pos_rec.append(pos.copy())
            vel_rec.append(vel.copy())
    e_final = energy(t_grid[-1])
    residual = e_final - e0 - w_total
    result = {
        "final_pos": pos, "final_vel": vel, "e0": e0, "e_final": e_final,
        "w_total": w_total, "w_move": w_move, "w_shift": w_shift, "w_depth": w_depth,
        "residual": residual, "max_residual": float(np.max(np.abs(residual))),
        "max_speed": max_speed, "max_radial": max_radial, "max_axial": max_axial,
        "final_z": pos[:, 2].copy(), "mean_z": z_sum / (t_grid[-1] - t_grid[0]),
    }
    if record_stride:
        result["trajectory"] = {
            "time_s": np.array(times_rec),
            "position_m": np.array(pos_rec),
            "velocity_m_per_s": np.array(vel_rec),
        }
    return result


def static_hold(pos, vel, duration_s, dt, mass_kg, depth_j, waist_m, z_r_m):
    """静态（阱静止、无透镜）保持段，返回 (末态pos, 末态vel, 峰峰值能量误差, 末能量)。"""
    n = int(round(duration_s / dt))
    control = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, depth_j, 0.0)
    accel, _ = make_force_fn(depth_j, waist_m, z_r_m, lambda t: control, mass_kg)
    p, v = pos.copy(), vel.copy()
    acc = accel(p, 0.0)
    energies = [total_energy(p, v, mass_kg, depth_j, waist_m, z_r_m, 0.0, 0.0, 0.0, 0.0)]
    for _ in range(n):
        p = p + v * dt + 0.5 * acc * dt ** 2
        acc_new = accel(p, 0.0)
        v = v + 0.5 * (acc + acc_new) * dt
        acc = acc_new
        energies.append(total_energy(p, v, mass_kg, depth_j, waist_m, z_r_m,
                                     0.0, 0.0, 0.0, 0.0))
    energies = np.array(energies)  # (n+1, N)
    e_final = total_energy(p, v, mass_kg, depth_j, waist_m, z_r_m, 0.0, 0.0, 0.0, 0.0)
    pp = np.max(energies, axis=0) - np.min(energies, axis=0)
    return p, v, pp, e_final
