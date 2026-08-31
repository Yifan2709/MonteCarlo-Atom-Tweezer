"""Level 2C 条件触发的 2D 横向模型（prompt 阶段 4 第 12 步）。

当 1D 模型加标定加热与位点差异仍达不到论文损失水平时，把径向一维扩成
x-y 二维：分裂移动沿 x 轴，双阱合并/分离的势垒结构与 1D 不同，且 y 方向
提供额外热能内容（equipartition），束缚判据用二维总能量。柱面透镜按论文
说法在该速度下不构成危害，明确排除在外。

引擎结构与 1D 版一致（同一往返程序、判定节点、早停记账），状态扩为
(x, y, vx, vy)；力为圆对称二维高斯的可分导数。
"""

from __future__ import annotations

import numpy as np

from .survival_engine import DisorderDraws, RoundtripLegs, _leg_steps


def _gaussian_force_2d(x, y, depth, waist, cx, cy=0.0):
    """二维圆对称高斯阱的力分量 (fx, fy)（支持数组逐 shot 参数）。"""
    qx = (x - cx) / waist
    qy = (y - cy) / waist
    envelope = np.exp(-2.0 * (qx**2 + qy**2))
    factor = -(4.0 * depth / waist) * envelope
    return factor * qx, factor * qy


def run_survival_2d(x0_m, y0_m, vx0_m_s, vy0_m_s, mass_kg, legs: RoundtripLegs, dt_s,
                    max_transfers: int, slm_depth_j, slm_waist_m, aod_waist_m,
                    heating=None, disorder: DisorderDraws | None = None,
                    jitter_sigma_m: float = 0.0):
    """2D 版逐 shot 连续转移模拟；记账格式与 1D 版一致。"""
    x0_m = np.asarray(x0_m, dtype=float)
    y0_m = np.asarray(y0_m, dtype=float)
    vx = np.asarray(vx0_m_s, dtype=float).copy()
    vy = np.asarray(vy0_m_s, dtype=float).copy()
    n0 = x0_m.size
    if disorder is None:
        disorder = DisorderDraws.disabled(n0)
    jitter_rng = np.random.default_rng(disorder.jitter_seed) if jitter_sigma_m > 0 else None

    x = x0_m.copy()
    y = y0_m.copy()
    shot_ids = np.arange(n0)
    noise_total = np.zeros(n0)
    transfers_survived = np.zeros(n0, dtype=int)
    lost_leg = np.array([""] * n0, dtype=object)
    lost_transfer = np.zeros(n0, dtype=int)
    judgement_energy_over_depth: list[list[float]] = [[] for _ in range(n0)]
    survival_counts = np.zeros(max_transfers + 1, dtype=int)
    survival_counts[0] = n0

    slm_depths = slm_depth_j * disorder.slm_depth_factor
    aod_waists = aod_waist_m * disorder.aod_waist_factor

    def force(xarr, yarr, center, depth_nominal, slm_on):
        fx, fy = _gaussian_force_2d(xarr, yarr, depth_nominal * disorder.aod_depth_factor,
                                    aod_waists, center)
        if slm_on:
            sx, sy = _gaussian_force_2d(xarr, yarr, slm_depths, slm_waist_m, 0.0)
            fx = fx + sx
            fy = fy + sy
        return fx / mass_kg, fy / mass_kg

    roundtrips = (max_transfers + 1) // 2
    transfer_index = 0
    for rt in range(roundtrips):
        if x.size == 0 or transfer_index >= max_transfers:
            break
        for center_arr, depth_arr, slm_on, judge in _leg_steps(legs):
            if transfer_index >= max_transfers:
                break
            steps = center_arr.size - 1
            # 逐段按当前活跃池尺寸重抽抖动，避免中途收缩后尺寸错位
            jitter = jitter_sigma_m * jitter_rng.standard_normal(x.size) \
                if jitter_rng is not None else 0.0
            offsets = disorder.alignment_offset_m + (jitter if np.ndim(jitter) else 0.0)
            ax, ay = force(x, y, center_arr[0] + offsets, depth_arr[0], slm_on)
            for i in range(steps):
                x = x + vx * dt_s + 0.5 * ax * dt_s**2
                y = y + vy * dt_s + 0.5 * ay * dt_s**2
                ax_new, ay_new = force(x, y, center_arr[i + 1] + offsets, depth_arr[i + 1], slm_on)
                vx = vx + 0.5 * (ax + ax_new) * dt_s
                vy = vy + 0.5 * (ay + ay_new) * dt_s
                ax, ay = ax_new, ay_new
                if heating is not None:
                    powers = heating.channel_powers(x, vx, 0.0, np.zeros_like(x))
                    delta = powers["total"] * dt_s
                    if np.any(delta != 0.0):
                        speed2 = np.maximum(0.0, vx**2 + vy**2 + 2.0 * delta / mass_kg)
                        speed_old = np.sqrt(np.maximum(vx**2 + vy**2, 1e-300))
                        scale = np.sqrt(speed2) / speed_old
                        vx = vx * scale
                        vy = vy * scale
                        noise_total = noise_total + delta
            if judge is None:
                continue
            transfer_index += 1
            if judge == "pickup_wait":
                depth_end = depth_arr[-1] * disorder.aod_depth_factor
                center_end = center_arr[-1] + offsets
                qx = (x - center_end) / aod_waists
                qy = y / aod_waists
                energy = 0.5 * mass_kg * (vx**2 + vy**2) - depth_end * np.exp(-2.0 * (qx**2 + qy**2))
                depth_scale = depth_end
            else:
                energy = 0.5 * mass_kg * (vx**2 + vy**2) \
                    - slm_depths * np.exp(-2.0 * ((x / slm_waist_m) ** 2 + (y / slm_waist_m) ** 2))
                depth_scale = slm_depths
            bound = energy < 0.0
            for pos, sid in enumerate(shot_ids):
                judgement_energy_over_depth[sid].append(float(energy[pos] / depth_scale[pos] + 1.0))
            transfers_survived[shot_ids[bound]] = transfer_index
            survival_counts[transfer_index] = int(bound.sum())
            if not bound.all():
                lost_now = shot_ids[~bound]
                lost_leg[lost_now] = judge
                lost_transfer[lost_now] = transfer_index
                keep = bound
                x, y, vx, vy = x[keep], y[keep], vx[keep], vy[keep]
                shot_ids = shot_ids[keep]
                disorder = DisorderDraws(disorder.slm_depth_factor[keep],
                                         disorder.aod_depth_factor[keep],
                                         disorder.aod_waist_factor[keep],
                                         disorder.alignment_offset_m[keep],
                                         disorder.jitter_seed)
                slm_depths = slm_depth_j * disorder.slm_depth_factor
                aod_waists = aod_waist_m * disorder.aod_waist_factor
            if x.size == 0:
                break

    return {
        "shots": n0,
        "max_transfers": int(max_transfers),
        "transfers_survived": transfers_survived,
        "lost_transfer": lost_transfer,
        "lost_leg": lost_leg.tolist(),
        "survival_counts": survival_counts,
        "survival_fraction": survival_counts / n0,
        "judgement_energy_over_depth": judgement_energy_over_depth,
        "noise_total_energy_uK": noise_total / 1.380649e-29,
        "final_active": int(x.size),
        "model": "2d_transverse",
    }
