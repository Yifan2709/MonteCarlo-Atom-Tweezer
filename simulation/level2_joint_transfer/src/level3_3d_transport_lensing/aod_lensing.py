"""论文 SI 启发的 crossed-AOD 柱面透镜势、校准模式与单阱/双阱诊断。"""
from __future__ import annotations

import numpy as np

from .gaussian_3d import physics_defaults


def axis_shift_m(axis_velocity_m_s, v_s_m_s, sigma, z_r_m):
    """柱面焦点轴向偏移 z_s = 2 σ z_R ċ / v_s（v_s=None 表示 lensing_off）。"""
    if v_s_m_s is None or v_s_m_s == 0:
        return np.zeros_like(np.asarray(axis_velocity_m_s, dtype=float))
    return 2.0 * sigma * z_r_m * np.asarray(axis_velocity_m_s, dtype=float) / v_s_m_s


def lensing_potential_and_force(xi, eta, z, depth_j, waist_m, z_r_m, zs1, zs2):
    """U = -D G1(ξ,z;z_s1) G2(η,z;z_s2) 及解析力，全部支持广播。

    返回 (U, F_xi, F_eta, F_z)。
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    z = np.asarray(z, dtype=float)
    dz1 = (z - zs1) / z_r_m
    dz2 = (z - zs2) / z_r_m
    a1 = 1.0 + dz1 ** 2
    a2 = 1.0 + dz2 ** 2
    b1 = 1.0 / np.sqrt(a1)
    b2 = 1.0 / np.sqrt(a2)
    e1 = np.exp(-2.0 * xi ** 2 / (waist_m ** 2 * a1))
    e2 = np.exp(-2.0 * eta ** 2 / (waist_m ** 2 * a2))
    g1 = b1 * e1
    g2 = b2 * e2
    u = -depth_j * g1 * g2
    # ∂G/∂u = G * (-4u/(w^2 A))；∂G/∂z = G * dz/zR * (-1/A + 4u^2/(w^2 A^2))
    dg1_dxi = g1 * (-4.0 * xi / (waist_m ** 2 * a1))
    dg2_deta = g2 * (-4.0 * eta / (waist_m ** 2 * a2))
    dg1_dz = g1 * dz1 / z_r_m * (-1.0 / a1 + 4.0 * xi ** 2 / (waist_m ** 2 * a1 ** 2))
    dg2_dz = g2 * dz2 / z_r_m * (-1.0 / a2 + 4.0 * eta ** 2 / (waist_m ** 2 * a2 ** 2))
    f_xi = depth_j * dg1_dxi * g2
    f_eta = depth_j * g1 * dg2_deta
    f_z = depth_j * (dg1_dz * g2 + g1 * dg2_dz)
    return u, f_xi, f_eta, f_z


def partial_t_potential(xi, eta, z, depth_j, waist_m, z_r_m, zs1, zs2,
                        c1_dot, c2_dot, zs1_dot, zs2_dot, depth_dot,
                        f_xi=None, f_eta=None, g1g2=None):
    """显式时变势的 ∂U/∂t（固定坐标），分解为移动/焦点偏移/深度三项。

    ∂U/∂t = F_ξ·ċ1 + F_η·ċ2 - Σ_j D·(∂G_j/∂z|u)·(other G)·ż_s,j - Ḋ·G1G2
    其中 ∂G_j/∂z_s = -∂G_j/∂z。若已算过力与 G1G2 可传入复用。
    返回 (total, move_term, shift_term, depth_term)。
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    z = np.asarray(z, dtype=float)
    if f_xi is None or f_eta is None or g1g2 is None:
        u, f_xi, f_eta, f_z = lensing_potential_and_force(
            xi, eta, z, depth_j, waist_m, z_r_m, zs1, zs2)
        g1g2 = -u / depth_j
    dz1 = (z - zs1) / z_r_m
    dz2 = (z - zs2) / z_r_m
    a1 = 1.0 + dz1 ** 2
    a2 = 1.0 + dz2 ** 2
    b1 = 1.0 / np.sqrt(a1)
    b2 = 1.0 / np.sqrt(a2)
    e1 = np.exp(-2.0 * xi ** 2 / (waist_m ** 2 * a1))
    e2 = np.exp(-2.0 * eta ** 2 / (waist_m ** 2 * a2))
    g1 = b1 * e1
    g2 = b2 * e2
    dg1_dz = g1 * dz1 / z_r_m * (-1.0 / a1 + 4.0 * xi ** 2 / (waist_m ** 2 * a1 ** 2))
    dg2_dz = g2 * dz2 / z_r_m * (-1.0 / a2 + 4.0 * eta ** 2 / (waist_m ** 2 * a2 ** 2))
    move = f_xi * c1_dot + f_eta * c2_dot
    # ∂U/∂z_s,j = +D·(∂G_j/∂z|_u)·(other G)（因 ∂G_j/∂z_s = -∂G_j/∂z，U 又带负号）
    shift = depth_j * (dg1_dz * g2 * zs1_dot + g1 * dg2_dz * zs2_dot)
    depth = -depth_dot * g1g2
    return move + shift + depth, move, shift, depth


def axial_minima(depth_j, waist_m, z_r_m, zs1, zs2, z_grid=None):
    """横向中心处瞬时轴向势的数值极值诊断。

    返回 (minima_z, minima_depths, barrier, n_minima)。barrier 为相邻最低点间
    的势垒高度（单阱时为 0）。
    """
    if z_grid is None:
        # 网格范围自适应覆盖两个焦点偏移，避免极小值落在域外
        reach = max(8.0, abs(zs1) / z_r_m + 4.0, abs(zs2) / z_r_m + 4.0)
        z_grid = np.linspace(-reach * z_r_m, reach * z_r_m, 4001)
    u, _, _, _ = lensing_potential_and_force(0.0, 0.0, z_grid, depth_j, waist_m, z_r_m,
                                             zs1, zs2)
    dz1 = (z_grid - zs1) / z_r_m
    dz2 = (z_grid - zs2) / z_r_m
    u_center = -depth_j / np.sqrt((1.0 + dz1 ** 2) * (1.0 + dz2 ** 2))
    u = u_center
    inner = u[1:-1]
    is_min = (inner <= u[:-2]) & (inner < u[2:])  # 一侧允许等号：底部平局保护
    idx = np.where(is_min)[0] + 1
    # 合并间隔极近的最低点（数值噪声）
    merged = []
    for i in idx:
        if merged and (z_grid[i] - z_grid[merged[-1]]) < 0.05 * z_r_m \
                and u[i] < u[merged[-1]]:
            merged[-1] = i
        elif merged and (z_grid[i] - z_grid[merged[-1]]) < 0.05 * z_r_m:
            continue
        else:
            merged.append(int(i))
    minima_z = z_grid[merged]
    minima_depths = u[merged]
    order = np.argsort(minima_depths)
    minima_z, minima_depths = minima_z[order], minima_depths[order]
    barrier = 0.0
    if len(merged) >= 2:
        barrier = float(np.max(u[min(merged[0], merged[1]):max(merged[0], merged[1]) + 1])
                        - min(minima_depths[:2]))
    return minima_z, minima_depths, barrier, len(merged)


def reference_axis_vmax(profile, distance_m, duration_s, geometry):
    """参考运输的单轴峰值速度（diagonal 为路径速度的 1/√2）。"""
    from .transport_profiles import get_profile
    q_dot_max = get_profile(profile).q_dot_max
    axis_factor = 1.0 / np.sqrt(2.0) if geometry == "diagonal" else 1.0
    return axis_factor * distance_m * q_dot_max / duration_s


def vs_from_ratio(ratio, profile, distance_m, duration_s, geometry):
    """由参考比值 r = v_max,axis/v_s 反解固定特征速度 v_s。"""
    if ratio <= 0:
        return None
    return reference_axis_vmax(profile, distance_m, duration_s, geometry) / ratio
