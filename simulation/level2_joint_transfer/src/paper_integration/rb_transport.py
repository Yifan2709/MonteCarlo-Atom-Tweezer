"""R1：Rb-87 移动光镊运输（Bluvstein 2022 图1/ED2 的结构复现）。

协议（R Methods）：线路运行时单 AOD 阱深 h×4 MHz 恒定，原子全程在同一阱内；
constant-jerk 轨迹 x(s)=D(3s²−2s³)（与 Cs 仓库的 pick-up/放回不同：无 SLM、
无交接，属 R 论文“线路内同一移动光镊”协议）。

输出：存活率与平均激发量子数 ΔN 随运输时长 T 的扫描 + 解析基准：
- R 等式(2)：常加速度 ΔN = ½(6D/x_zpf/(ω₀T)²)²（constant-jerk 为其 5–10 倍
  加速下界口径）；x_zpf=√(ħ/(2mω₀))。
- erf 留存模型 retention=½(1−erf[(N−Nmax)/σ])，Nmax∈[26,33] 作对照带。
结构判据（绝对曲线不验收，R 图1d/ED2 源数据未获取）：
  C1 存活率在长 T 处 →1，短 T 处下降（存在速度阈值）；
  C2 数值 ΔN 与解析 1/T⁴ 标度一致（对数斜率 −4±0.6）；
  C3 R 图1 锚点：D=110 μm, T=300 μs（b≈0.37 μm/μs < 0.55）存活率高。
"""
from __future__ import annotations

import numpy as np

from .dynamics3d import (Segment, TrapPhysics, propagate,
                         sample_bound_thermal_3d)
from level3_3d_transport_lensing.gaussian_3d import static_potential
from .provenance import REGISTRIES
from .trajectories import poly_trajectory

HBAR = 1.054571817e-34
H_PLANCK = 6.62607015e-34
KB = 1.380649e-23


def rb_physics() -> TrapPhysics:
    rb = REGISTRIES["rb87_2022_native"]
    depth = rb.get("depth_circuit_hz").value * H_PLANCK
    return TrapPhysics(mass_kg=rb.get("mass_u").value * 1.6605390666e-27,
                       depth_j=depth, waist_m=rb.get("waist_m").value,
                       zr_m=rb.get("zr_m").value)


def run_single_transport(distance_m: float, duration_s: float, shots: int,
                         seed: int, dt_s: float, temperature_uK: float = 5.0,
                         traj_kind: str = "constant_jerk") -> dict:
    rb = REGISTRIES["rb87_2022_native"]
    ph = rb_physics()
    pos, vel = sample_bound_thermal_3d(seed, temperature_uK, shots, ph)
    traj = poly_trajectory(traj_kind, distance_m, duration_s, dt_s)
    seg = Segment("move", traj.t_s, traj.x_m, np.zeros_like(traj.x_m),
                  np.full_like(traj.t_s, ph.depth_j), traj.v_m_s)
    res = propagate(pos, vel, [seg], ph, dt_s)
    omega_r = rb.get("omega_r_circuit_rad_s").value
    alive_e = res["final_energy_j"][res["alive"]] - (-ph.depth_j)
    dn = float(np.mean(alive_e)) / (HBAR * omega_r) if res["alive"].any() \
        else float("nan")
    # 初态热激发基线（同批初态，静态阱能量）
    e0 = (0.5 * ph.mass_kg * np.sum(vel ** 2, axis=1)
          + static_potential(pos[:, 0], pos[:, 1], pos[:, 2], ph.depth_j,
                             ph.waist_m, ph.zr_m)) + ph.depth_j
    dn0 = float(np.mean(e0)) / (HBAR * omega_r)
    return {
        "distance_m": distance_m, "duration_s": duration_s,
        "traj": traj_kind, "shots": shots, "seed": seed, "dt_s": dt_s,
        "temperature_uK": temperature_uK,
        "survival": res["survival_final"],
        "alive": res["alive"],
        "final_energy_j": res["final_energy_j"],
        "wilson95_halfwidth": _wilson_half(res["survival_final"], shots),
        "mean_delta_N": dn,          # 含初态热激发
        "mean_delta_N_initial": dn0,
        "mean_delta_N_excess": dn - dn0,
        "loss_stage_counts": {s: int(np.sum(res["loss_stage"] == s))
                              for s in np.unique(res["loss_stage"]) if s},
    }


def analytic_delta_n_const_accel(distance_m: float, duration_s: float,
                                 omega_rad_s: float, mass_kg: float) -> float:
    """Legacy API name: R Eq.(2) is cubic/constant-JERK, frequency-averaged.

    Requires omega*T >> 1 and frequency averaging; not a constant-accel result.
    """
    x_zpf = np.sqrt(HBAR / (2.0 * mass_kg * omega_rad_s))
    ratio = 6.0 * distance_m / (x_zpf * omega_rad_s ** 2 * duration_s ** 2)
    return 0.5 * ratio ** 2


def erf_retention(delta_n: float, n_max: float, sigma: float = 2.0) -> float:
    from math import erf
    return 0.5 * (1.0 - erf((delta_n - n_max) / sigma))


def scan_transports(distance_m: float, durations_s: list[float], shots: int,
                    seed: int, dt_s: float,
                    temperature_uK: float = 5.0) -> list[dict]:
    rows = []
    for i, t in enumerate(durations_s):
        rows.append(run_single_transport(distance_m, t, shots, seed + 101 * i,
                                         dt_s, temperature_uK))
    return rows


def evaluate_structure_criteria(rows: list[dict]) -> dict:
    out: dict = {}
    ts = np.array([r["duration_s"] for r in rows])
    sv = np.array([r["survival"] for r in rows])
    dn = np.array([r["mean_delta_N"] for r in rows])
    out["C1_survival_high_at_long_T"] = bool(sv[-1] > 0.98
                                             and sv[0] < sv[-1])
    # C2：ΔN ∝ T^-4 拟合（仅用干净标度窗口：过剩激发 0.5<N<50，
    # 避开短 T 饱和/幸存者偏置与长 T 热噪声底）
    ex = np.array([r.get("mean_delta_N_excess", float("nan")) for r in rows])
    mask = np.isfinite(ex) & (ex > 0.3) & (ex < 100.0)
    dn_fit = ex
    if mask.sum() >= 3:
        slope = float(np.polyfit(np.log(ts[mask]), np.log(dn_fit[mask]), 1)[0])
        out["C2_deltaN_scaling_slope"] = slope
        out["C2_scaling_ok"] = bool(-4.6 < slope < -3.4)
    rb = REGISTRIES["rb87_2022_native"]
    d_paper = rb.get("bell_move_distance_m").value
    t_paper = rb.get("bell_move_time_s").value
    anchor = [r for r in rows
              if abs(r["distance_m"] - d_paper) < 1e-12
              and abs(r["duration_s"] - t_paper) < 1e-12]
    if anchor:
        out["C3_anchor_55um_per_atom_300us"] = {
            "survival": anchor[0]["survival"],
            "atom_speed_m_per_s": d_paper / t_paper,
            "pair_separation_speed_m_per_s": 2 * d_paper / t_paper,
            "paper_flat_speed": rb.get("bell_flat_speed_m_per_s").value,
            "ok": bool(anchor[0]["survival"] > 0.95
                       and 2 * d_paper / t_paper
                       <= rb.get("bell_flat_speed_m_per_s").value),
        }
    return out


def _wilson_half(p: float, n: int) -> float:
    z = 1.96
    denom = 1 + z * z / n
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(half)
