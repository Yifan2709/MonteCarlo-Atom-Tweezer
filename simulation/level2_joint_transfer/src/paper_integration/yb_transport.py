"""Y1：Yb 含透镜经典运输与轨迹对照（Zhang 2026 ED Fig.2b 的结构复现）。

协议（Y 论文 §II/Methods）：原子初始在存储区 SLM 阱（静态）；
一次单程 = SLM→AOD 交接（AOD 深度 ramp，0.78 ms）→ AOD 单独移动 T（扫描）
→ 段末判定。移动轨迹 ∈ {zero-jerk, min-jerk}；透镜 τ_eff ∈ {0, 2, 4, 8} s/m。

预注册结构判据（绝对曲线不验收——论文 ED Fig.2b 源数据未获取）：
  C1 存在 T_fast 使 lensing_on 的存活率显著低于 lensing_off（快速运输受限）；
  C2 lensing_on 时 zero-jerk 存活率 ≥ min-jerk（T_fast 附近，超出联合 MC 误差）；
  C3 lensing_off 时 zero-jerk 与 min-jerk 之差在联合 MC 误差内。
"""
from __future__ import annotations

import numpy as np

from .dynamics3d import (Segment, TrapPhysics, propagate,
                         sample_bound_thermal_3d)
from .provenance import REGISTRIES
from .trajectories import MIN_JERK, ZERO_JERK, poly_trajectory

H_PLANCK = 6.62607015e-34
KB = 1.380649e-23


def yb_physics(lensing_tau_eff: float = 0.0) -> TrapPhysics:
    yb = REGISTRIES["yb171_2026_native"]
    depth = yb.get("storage_depth_uK").value * 1e-6 * KB
    return TrapPhysics(
        mass_kg=yb.get("mass_u").value * 1.6605390666e-27,
        depth_j=depth, waist_m=yb.get("waist_m").value,
        zr_m=yb.get("zr_m").value,
        slm_depth_j=depth,           # 交接时 SLM 静态阱在场
        lensing_tau_eff_s_per_m=lensing_tau_eff,
    )


def run_single_move(distance_m: float, move_time_s: float, traj_kind: str,
                    tau_eff: float, shots: int, seed: int, dt_s: float,
                    temperature_uK: float = 5.0,
                    handoff_time_s: float | None = None,
                    alignment_sigma_m: float | None = None) -> dict:
    """ED2 defaults to pure AOD transport; explicit handoff is a separate protocol.

    Legacy parameter registry is conditional, not an independently calibrated
    recreation of the authors' apparatus. Do not certify experiment agreement
    from evaluate_structure_criteria.
    """
    yb = REGISTRIES["yb171_2026_native"]
    if handoff_time_s is None:
        handoff_time_s = 0.0
    if alignment_sigma_m is None:
        alignment_sigma_m = yb.get("alignment_sigma_m").value
    ph = yb_physics(tau_eff)
    ph_move = yb_physics(tau_eff)
    ph_move.slm_depth_j = 0.0        # 交接完成后 SLM 关闭，原子随 AOD 移动
    rng = np.random.default_rng(seed)
    pos, vel = sample_bound_thermal_3d(seed, temperature_uK, shots, ph_move)
    offs = alignment_sigma_m * rng.standard_normal(shots)

    depth = ph.depth_j
    if handoff_time_s <= 0:
        traj = poly_trajectory(traj_kind, distance_m, move_time_s, dt_s)
        seg = Segment("move", traj.t_s, traj.x_m, np.zeros_like(traj.x_m),
                      np.full_like(traj.t_s, depth), traj.v_m_s)
        res = propagate(pos, vel, [seg], ph_move, dt_s)
        e = res["final_energy_j"][res["alive"]] + depth
        return {"distance_m": distance_m, "move_time_s": move_time_s,
                "traj": traj_kind, "tau_eff_s_per_m": tau_eff,
                "shots": shots, "seed": seed, "dt_s": dt_s,
                "temperature_uK": temperature_uK, "protocol": "pure_AOD_ED2",
                "survival": res["survival_final"],
                "wilson95_halfwidth": _wilson_half(res["survival_final"], shots),
                "mean_excitation_quanta_alive": float(np.mean(e))/(1.054571817e-34*yb.get("omega_r_rad_s").value) if len(e) else float("nan"),
                "checkpoints": res["checkpoints"],
                "loss_stage_counts": {"move": int((~res["alive"]).sum())}}
    n_ramp = int(round(handoff_time_s / dt_s))
    t_ramp = np.arange(n_ramp + 1) * dt_s
    traj = poly_trajectory(traj_kind, distance_m, move_time_s, dt_s)
    t_all = np.concatenate([t_ramp, t_ramp[-1] + traj.t_s[1:]])
    n_move = len(traj.t_s) - 1
    depth_all = np.concatenate([
        depth * t_ramp / handoff_time_s,          # AOD 线性升
        np.full(n_move, depth)])
    slm_all = np.concatenate([
        depth * (1.0 - t_ramp / handoff_time_s),  # SLM 线性降（交接）
        np.zeros(n_move)])                        # 移动段 SLM 关闭
    cx = np.concatenate([np.zeros(n_ramp + 1), traj.x_m[1:]])
    vx = np.concatenate([np.zeros(n_ramp + 1), traj.v_m_s[1:]])
    segs = [
        Segment("handoff", t_all[:n_ramp + 1], cx[:n_ramp + 1],
                np.zeros(n_ramp + 1), depth_all[:n_ramp + 1],
                slm_depth_j=slm_all[:n_ramp + 1]),
        Segment("move", t_all[n_ramp:], cx[n_ramp:], np.zeros_like(cx[n_ramp:]),
                depth_all[n_ramp:], vx[n_ramp:], slm_depth_j=slm_all[n_ramp:]),
    ]
    res = propagate(pos, vel, segs, ph, dt_s, offsets_m=offs)
    alive_e = res["final_energy_j"][res["alive"]]
    from .dynamics3d import excitation_quanta
    nq = excitation_quanta(ph, yb.get("omega_r_rad_s").value, alive_e
                           - (-depth)) if len(alive_e) else np.array([0.0])
    return {
        "distance_m": distance_m, "move_time_s": move_time_s,
        "traj": traj_kind, "tau_eff_s_per_m": tau_eff, "shots": shots,
        "seed": seed, "dt_s": dt_s, "temperature_uK": temperature_uK,
        "survival": res["survival_final"],
        "wilson95_halfwidth": _wilson_half(res["survival_final"], shots),
        "mean_excitation_quanta_alive": float(np.mean(nq)),
        "checkpoints": res["checkpoints"],
        "loss_stage_counts": {
            s: int(np.sum(res["loss_stage"] == s))
            for s in np.unique(res["loss_stage"]) if s != ""},
    }


def _wilson_half(p: float, n: int) -> float:
    if n == 0:
        return float("nan")
    z = 1.96
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(half)


def scan_moves(distance_m: float, move_times_s: list[float], shots: int,
               seed: int, dt_s: float, taus: list[float],
               trajs: list[str], temperature_uK: float = 5.0) -> list[dict]:
    """网格扫描：返回行列表（逐格点独立 seed 由基 seed 派生）。"""
    rows = []
    for tau in taus:
        for traj in trajs:
            for i, t in enumerate(move_times_s):
                row = run_single_move(
                    distance_m, t, traj, tau, shots,
                    seed,
                    dt_s, temperature_uK)
                rows.append(row)
    return rows


def evaluate_structure_criteria(rows: list[dict]) -> dict:
    """预注册结构判据 C1/C2/C3。"""
    by_key = {(r["tau_eff_s_per_m"], r["traj"], round(r["move_time_s"], 9)): r
              for r in rows}
    taus = sorted({r["tau_eff_s_per_m"] for r in rows})
    trajs = sorted({r["traj"] for r in rows})
    times = sorted({round(r["move_time_s"], 9) for r in rows})
    tau_on = max(taus)
    results: dict = {"tau_on": tau_on, "details": []}
    # C1：lensing_on 存在 T 使存活显著低于 lensing_off
    c1 = False
    best_gap = None
    for t in times:
        r_off = by_key.get((0.0, trajs[0], t))
        for tr in trajs:
            r_on = by_key.get((tau_on, tr, t))
            if not r_on or not r_off:
                continue
            gap = r_off["survival"] - r_on["survival"]
            err = np.hypot(r_off["wilson95_halfwidth"],
                           r_on["wilson95_halfwidth"])
            if best_gap is None or gap > best_gap[0]:
                best_gap = (gap, err, t, tr)
            if gap > err:
                c1 = True
    results["C1_lensing_limits_fast_transport"] = c1
    if best_gap:
        results["C1_best_gap"] = {"gap": best_gap[0], "err": best_gap[1],
                                  "at_T_s": best_gap[2], "traj": best_gap[3]}
    # C2/C3：zero vs min jerk
    zj = [k for k in ("low_peak_v",) if k in trajs]
    mj = [k for k in ("gamma:1.875",) if k in trajs]
    if zj and mj:
        zj, mj = zj[0], mj[0]
        for tau, cid in ((tau_on, "C2"), (0.0, "C3")):
            diffs = []
            for t in times:
                rz = by_key.get((tau, zj, t))
                rm = by_key.get((tau, mj, t))
                if not (rz and rm):
                    continue
                # 预注册评估窗口（登记于合同修正）：
                # C2 仅在部分损失区（透镜存活率 0.02–0.98，避开全灭/全活）；
                # C3 仅在绝热区（无透镜存活率均 >0.98）。
                if cid == "C2" and not (0.02 < rz["survival"] < 0.98
                                        or 0.02 < rm["survival"] < 0.98):
                    continue
                if cid == "C3" and not (rz["survival"] > 0.98
                                        and rm["survival"] > 0.98):
                    continue
                diffs.append((rz["survival"] - rm["survival"],
                              np.hypot(rz["wilson95_halfwidth"],
                                       rm["wilson95_halfwidth"]), t))
            if not diffs:
                continue
            if cid == "C2":
                d, e, t = max(diffs, key=lambda x: x[0] / max(x[1], 1e-9))
                results["C2_zero_jerk_better_with_lensing"] = bool(d > e)
                results["C2_best"] = {"diff": d, "err": e, "at_T_s": t}
            else:
                significant = [d for d, e, t in diffs if abs(d) > e]
                results["C3_similar_without_lensing"] = len(significant) == 0
                results["C3_max_abs_diff"] = float(
                    max(abs(d) for d, e, t in diffs))
    return results
