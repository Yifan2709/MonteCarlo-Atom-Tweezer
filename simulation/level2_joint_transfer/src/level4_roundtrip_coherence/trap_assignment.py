"""重叠势 basin 分类与 checkpoint 局部束缚判据。

分类输出：bound_to_slm / bound_to_aod / shared_or_ambiguous / unbound。
单阱主导时用 E_trap<0 判据；两阱都显著时沿连线找总势最低点与分隔势垒，
把原子分配到相空间可达的 basin；能量接近势垒或深度接近零时标 ambiguous。
"""
from __future__ import annotations

import numpy as np

from level3_3d_transport_lensing.aod_lensing import lensing_potential_and_force
from level3_3d_transport_lensing.gaussian_3d import static_potential

BASINS = ("bound_to_slm", "bound_to_aod", "shared_or_ambiguous", "unbound")


def slm_potential(pos, d_slm, waist_m, zr_m, center):
    """SLM 三维高斯势（中心静止）。"""
    return static_potential(pos[..., 0] - center[0], pos[..., 1] - center[1],
                            pos[..., 2], d_slm, waist_m, zr_m)


def aod_potential(pos, d_aod, waist_m, zr_m, c1, c2, zs1, zs2):
    """AOD crossed-lens 势（含焦点偏移）。"""
    u, _, _, _ = lensing_potential_and_force(pos[..., 0] - c1, pos[..., 1] - c2,
                                             pos[..., 2], d_aod, waist_m, zr_m,
                                             zs1, zs2)
    return u


def trap_energies(pos, vel, ctrl, phys, v_trap=(0.0, 0.0, 0.0)):
    """分别相对 SLM/AOD 静态瞬时阱的 E_trap（含相对阱速度的动能）。"""
    mass = phys.mass_kg
    rel_v = vel - np.asarray(v_trap)
    ke = 0.5 * mass * np.sum(rel_v ** 2, axis=-1)
    u_slm = slm_potential(pos, ctrl["d_slm"], phys.slm_waist_m, phys.slm_zr_m,
                          phys.slm_center_m)
    u_aod = aod_potential(pos, ctrl["d_aod"], phys.aod_waist_m, phys.aod_zr_m,
                          ctrl["c1"], ctrl["c2"], ctrl["zs1"], ctrl["zs2"])
    return ke + u_slm, ke + u_aod, u_slm, u_aod


def total_potential_line(u_axis, c_slm, c_aod, d_slm, d_aod, phys, n=481):
    """沿 SLM→AOD 连线（z=0 切片）的总势采样。

    采样范围扩展到 [−0.35, 1.35]，使两阱最低点都落在数组内部，
    否则端点处的阱最低点不会被内部极小值搜索发现。
    """
    t = np.linspace(-0.35, 1.35, n)
    line = c_slm[None, :] + (c_aod - c_slm)[None, :] * t[:, None]
    zeros = np.zeros(n)
    u_s = static_potential(line[:, 0] - c_slm[0], line[:, 1] - c_slm[1], zeros,
                           d_slm, phys.slm_waist_m, phys.slm_zr_m)
    u_a = lensing_potential_and_force(line[:, 0] - c_aod[0], line[:, 1] - c_aod[1],
                                      zeros, d_aod, phys.aod_waist_m,
                                      phys.aod_zr_m, 0.0, 0.0)[0]
    return t, u_s + u_a


def line_minima_and_barrier(u_line, depth_scale):
    """连线势的最低点索引、鞍点与势垒高度（相对较深最低点）。"""
    inner = u_line[1:-1]
    is_min = (inner <= u_line[:-2]) & (inner < u_line[2:])
    idx = np.where(is_min)[0] + 1
    # 合并近邻最低点（数值噪声保护，间距 < 2 格）
    merged = []
    for i in idx:
        if merged and i - merged[-1] <= 2:
            if u_line[i] < u_line[merged[-1]]:
                merged[-1] = int(i)
        else:
            merged.append(int(i))
    if len(merged) < 2:
        return merged, None, 0.0
    a, b = merged[0], merged[1]
    saddle = int(a + np.argmax(u_line[a:b + 1]))
    barrier = float(u_line[saddle] - min(u_line[a], u_line[b]))
    return merged, saddle, barrier


def classify_basin(pos, vel, ctrl, phys, ambiguous_tol=1e-3):
    """单个或一批原子在瞬时控制量下的 basin 分类。

    返回 dict：basin（str 或 list）、e_slm、e_aod、near_slm_threshold、
    near_aod_threshold、ambiguous、sep_m、n_line_minima、barrier_j。
    """
    single = np.asarray(pos).ndim == 1
    pos_arr = np.atleast_2d(np.asarray(pos, dtype=float))
    vel_arr = np.atleast_2d(np.asarray(vel, dtype=float))
    e_slm, e_aod, u_slm, u_aod = trap_energies(pos_arr, vel_arr, ctrl, phys)
    d_slm, d_aod = float(ctrl["d_slm"]), float(ctrl["d_aod"])
    c_slm = np.array([phys.slm_center_m[0], phys.slm_center_m[1]])
    c_aod = np.array([float(ctrl["c1"]), float(ctrl["c2"])])
    sep = float(np.linalg.norm(c_aod - c_slm))
    depth_scale = max(d_slm, d_aod, 1e-30)
    near_slm = np.abs(e_slm) < ambiguous_tol * max(d_slm, 1e-30) * 10
    near_aod = np.abs(e_aod) < ambiguous_tol * max(d_aod, 1e-30) * 10

    merged_sep = phys.aod_waist_m  # 距离小于一个束腰视为合并阱
    both_on = (d_slm > 1e-6 * depth_scale) and (d_aod > 1e-6 * depth_scale)
    basins = []
    n_minima = 0
    barrier = 0.0
    u_tot_pos = u_slm + u_aod
    e_tot = 0.5 * phys.mass_kg * np.sum(vel_arr ** 2, axis=1) + u_tot_pos
    # 连线分析只依赖瞬时控制量，批量分类时仅计算一次
    t_line = saddle = None
    if both_on and sep >= merged_sep:
        t_line, u_line = total_potential_line(None, c_slm, c_aod, d_slm, d_aod, phys)
        minima, saddle, barrier = line_minima_and_barrier(u_line, depth_scale)
        n_minima = len(minima)
    for i in range(len(pos_arr)):
        if not both_on:
            if d_aod <= 1e-6 * depth_scale:
                basins.append("bound_to_slm" if e_slm[i] < 0 else "unbound")
            else:
                basins.append("bound_to_aod" if e_aod[i] < 0 else "unbound")
            continue
        if sep < merged_sep:
            # 合并阱：以总势负能量判束缚，归类为 shared
            basins.append("shared_or_ambiguous" if e_tot[i] < 0 else "unbound")
            continue
        # 双阱显著：沿连线找最低点与势垒
        if n_minima < 2:
            # 数值上未分辨出双阱：按更深的那个阱归属
            deep = "aod" if d_aod >= d_slm else "slm"
            e_deep = e_aod[i] if deep == "aod" else e_slm[i]
            basins.append(f"bound_to_{deep}" if e_deep < 0 else "unbound")
            continue
        u_saddle = u_line[saddle]
        tol = ambiguous_tol * depth_scale
        if e_tot[i] >= 0:
            # E_tot ≥ 0 且两势均非正 → 对两阱都不束缚（已逃逸）
            basins.append("unbound")
            continue
        # 相空间可达性：束缚但总能高于势垒 → 两 basin 均可达
        if e_tot[i] > u_saddle - tol:
            basins.append("shared_or_ambiguous")
            continue
        # 位置投影到连线，与鞍点比较
        vec = c_aod - c_slm
        u_proj = float(np.dot(pos_arr[i, :2] - c_slm, vec) / np.dot(vec, vec))
        basin = "aod" if u_proj > t_line[saddle] else "slm"
        e_basin = e_aod[i] if basin == "aod" else e_slm[i]
        basins.append(f"bound_to_{basin}" if e_basin < 0 else "unbound")

    result = {
        "e_slm": e_slm, "e_aod": e_aod, "u_slm": u_slm, "u_aod": u_aod,
        "e_total": e_tot,
        "near_slm_threshold": near_slm, "near_aod_threshold": near_aod,
        "separation_m": sep, "n_line_minima": n_minima,
        "barrier_j": barrier if isinstance(barrier, float) else float(barrier),
        "d_slm": d_slm, "d_aod": d_aod,
    }
    ambiguous = np.array([b == "shared_or_ambiguous" for b in basins])
    result["ambiguous"] = ambiguous
    if single:
        result["basin"] = basins[0]
        for key in ("e_slm", "e_aod", "u_slm", "u_aod", "e_total",
                    "near_slm_threshold", "near_aod_threshold", "ambiguous"):
            result[key] = result[key][0]
    else:
        result["basin"] = basins
    return result


def frozen_dynamics_assignment(pos, vel, ctrl, phys, duration_s, dt, n_periods=3):
    """诊断副本法：短暂冻结控制推进若干振动周期确定归属（测试/校验用）。"""
    from .roundtrip_simulation import integrate_static
    pos_f, vel_f = integrate_static(pos, vel, ctrl, phys, duration_s, dt)
    disp_slm = np.linalg.norm(pos_f[..., :2] - np.array(phys.slm_center_m),
                              axis=-1)
    disp_aod = np.linalg.norm(pos_f[..., :2]
                              - np.array([ctrl["c1"], ctrl["c2"]]), axis=-1)
    return np.where(disp_slm < disp_aod, "bound_to_slm", "bound_to_aod")
