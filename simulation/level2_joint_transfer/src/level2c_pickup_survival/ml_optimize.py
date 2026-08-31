"""Level 2C ML 风格轨迹优化：14+14 控制点普通三次插值，目标为固定次数连续转移存活率。

与论文的差异（诚实声明）：论文用在线 ML 优化器在实验装置上直接优化存活率；
此处用等价的黑盒优化（LHS 全局采样 + 局部精修）在模拟器上优化同一目标
（固定 transfers_for_objective 次连续单程转移的存活率），优化器形式不同、
目标与搜索空间（深度/位置各 14 控制点、三次插值、ramp-and-move 同时）相同。

防过拟合纪律与 Level 2 一致：候选只用训练池（candidate_seed），冻结后用
独立复核池（roundtrip.validation_seed）一次性评估，不得据复核结果回头改参。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from level2_joint_transfer.level2_simulation import physics_from_config
from level2_joint_transfer.waveforms import MLCubicWaveform, validate_pickup_waveform

from .level2c_config import Level2cConfig
from .survival_engine import build_roundtrip_legs, run_survival


def control_grid(num_points: int) -> np.ndarray:
    """均匀控制点归一化网格（含端点 0 与 1）。"""
    return np.linspace(0.0, 1.0, int(num_points))


def manual_anchor_controls(cfg: Level2cConfig, duration_us: float) -> tuple[np.ndarray, np.ndarray]:
    """手工轨迹在控制网格上的取值（作为 ML 候选 0 的锚点）。"""
    from level2_joint_transfer.level2_simulation import build_waveform
    physics = physics_from_config(cfg.base)
    spec = {"type": "pickup_sequential", "duration_s": duration_us * 1e-6,
            "ramp_fraction": cfg.pickup.ramp_fraction, "pickup_direction": "pickup",
            "direction": "pickup"}
    wf = build_waveform(spec, physics)
    s = control_grid(cfg.pickup.ml_control_points)
    return (np.asarray(wf.center(s * wf.duration_s), dtype=float),
            np.asarray(wf.depth(s * wf.duration_s), dtype=float))


def build_ml_waveform(cfg: Level2cConfig, center_values_m, depth_values_j, duration_us=None) -> MLCubicWaveform:
    """由控制点构造 ML 三次插值波形（端点固定 0→d、0→D0）。"""
    physics = physics_from_config(cfg.base)
    duration_s = (duration_us if duration_us is not None else cfg.pickup.ml_duration_us) * 1e-6
    return MLCubicWaveform(duration_s, control_grid(cfg.pickup.ml_control_points),
                           np.asarray(center_values_m, dtype=float),
                           np.asarray(depth_values_j, dtype=float),
                           physics["aod_initial_center_m"], physics["aod_initial_depth_j"])


def evaluate_ml_candidate(cfg: Level2cConfig, center_values_m, depth_values_j,
                          pool: dict, heating) -> dict:
    """评估单个 ML 候选：固定次数连续转移的存活率（主）+ 激发/光滑（次）。

    注意：legs 必须按候选自身的波形重新构建——曾因按时长缓存 legs 导致
    所有 LHS 候选被评估在锚点波形上（训练成绩虚高、冻结候选复核崩塌），
    该 bug 已修复并加入回归测试。
    """
    physics = physics_from_config(cfg.base)
    ml = cfg.ml_optimization
    waveform = build_ml_waveform(cfg, center_values_m, depth_values_j)
    valid, reasons = validate_pickup_waveform(
        waveform, physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
        direction="pickup")
    row = {"status": "ok" if valid else "invalid", "reasons": "; ".join(reasons)}
    if not valid:
        return row
    legs = build_roundtrip_legs(waveform, cfg.roundtrip.wait_s, ml.dt_us * 1e-6)
    out = run_survival(pool["x0_m"], pool["v0_m_per_s"], physics["mass_kg"], legs,
                       ml.dt_us * 1e-6, ml.transfers_for_objective,
                       slm_depth_j=physics["slm_depth_j"], slm_waist_m=physics["slm_waist_m"],
                       aod_waist_m=physics["aod_waist_m"], heating=heating)
    n_obj = ml.transfers_for_objective
    survival = float(out["survival_fraction"][n_obj])
    # 次级目标：末次判定仍存活 shot 的平均归一化激发（越低越好；全丢失时为空）
    last_energy = [e[n_obj - 1] for e in out["judgement_energy_over_depth"] if len(e) >= n_obj]
    mean_excitation = float(np.mean(last_energy)) if last_energy else None
    smooth = waveform.smoothness_metrics(801)
    row.update({
        "survival_at_objective": survival,
        "survived_count": int(out["survival_counts"][n_obj]),
        "mean_final_excitation_alive": mean_excitation,
        "max_abs_center_velocity_m_per_s": smooth["max_abs_center_velocity_m_per_s"],
        "max_abs_center_acceleration_m_per_s2": smooth["max_abs_center_acceleration_m_per_s2"],
        "objective_total": survival
            - 0.05 * (mean_excitation if mean_excitation is not None else 1.0)
            - 1e-6 * smooth["max_abs_center_acceleration_m_per_s2"],
    })
    return row


# LHS 增量参数化的相对散布：inc = 1 + SPREAD·(2u−1)，inc>0 保持大尺度单调、
# 局部形状自由；5% 过冲校验下有效率 ~90%。非单调候选由精修阶段的直接
# 控制点扰动覆盖（见 run_ml_optimization）。直接对控制点值做 LHS 会因
# 三次样条大幅过冲几乎全军覆没（实测 63/63 invalid），故用增量参数化。
ML_LHS_SPREAD = 0.9


def lhs_ml_candidates(cfg: Level2cConfig) -> list[dict]:
    """LHS 生成候选：候选 0 为手工锚点，其余按增量参数化采样。

    中心/深度各由 13 段归一化增量生成（自由度 = 内部 12 点 × 2 = 24，
    与论文 14+14 控制点（端点固定）一致），端点严格 0→d 与 0→D0。
    """
    physics = physics_from_config(cfg.base)
    ml = cfg.ml_optimization
    n_ctrl = cfg.pickup.ml_control_points
    n_seg = n_ctrl - 1
    d = physics["aod_initial_center_m"]
    d0 = physics["aod_initial_depth_j"]
    anchor_c, anchor_d = manual_anchor_controls(cfg, cfg.pickup.ml_duration_us)
    candidates = [{"center_values_m": anchor_c.tolist(), "depth_values_j": anchor_d.tolist(),
                   "origin": "manual_anchor"}]
    sampler = qmc.LatinHypercube(d=2 * n_seg, seed=ml.candidate_seed)
    unit = sampler.random(ml.candidates - 1)
    for row in unit:
        u_c, u_d = row[:n_seg], row[n_seg:]
        inc_c = 1.0 + ML_LHS_SPREAD * (2.0 * u_c - 1.0)
        inc_d = 1.0 + ML_LHS_SPREAD * (2.0 * u_d - 1.0)
        cum_c = np.concatenate([[0.0], np.cumsum(inc_c)])
        cum_d = np.concatenate([[0.0], np.cumsum(inc_d)])
        center = cum_c / cum_c[-1] * d
        depth = cum_d / cum_d[-1] * d0
        candidates.append({"center_values_m": center.tolist(), "depth_values_j": depth.tolist(),
                           "origin": "lhs_increments"})
    return candidates


def run_ml_optimization(cfg: Level2cConfig, pool: dict, heating, output_dir: Path,
                        checkpoint_path: Path) -> list[dict]:
    """执行 LHS 评估 + 局部精修；支持检查点续跑。返回全部候选行。"""
    ml = cfg.ml_optimization
    done = {}
    if checkpoint_path.exists():
        try:
            saved = json.loads(checkpoint_path.read_text())
            done = {row["candidate_id"]: row for row in saved.get("rows", [])}
        except json.JSONDecodeError:
            done = {}
    rows = []
    candidates = lhs_ml_candidates(cfg)
    for idx, cand in enumerate(candidates):
        cid = f"ml_{idx:04d}"
        if cid in done:
            rows.append(done[cid])
            continue
        row = {"candidate_id": cid, "origin": cand["origin"],
               "center_values_m": cand["center_values_m"], "depth_values_j": cand["depth_values_j"]}
        try:
            row.update(evaluate_ml_candidate(cfg, cand["center_values_m"], cand["depth_values_j"],
                                             pool, heating))
        except Exception as exc:  # 单个候选失败记为 failed，不中断整轮搜索
            row.update({"status": "failed", "reasons": repr(exc)})
        rows.append(row)
        checkpoint_path.write_text(json.dumps({"rows": rows}, ensure_ascii=False))
        print(f"  [ml] {cid} {row['status']} survival={row.get('survival_at_objective')}")

    # 局部精修：取前 finalists 个 ok 候选做小扰动
    ok_rows = [r for r in rows if r.get("status") == "ok"]
    ok_rows.sort(key=lambda r: (-r["survival_at_objective"],
                                r.get("mean_final_excitation_alive") or 1e9))
    rng = np.random.default_rng(ml.candidate_seed + 1)
    refine_rows = []
    for rank, base in enumerate(ok_rows[:ml.finalists]):
        for step in range(ml.refine_steps):
            cid = f"ml_refine_{rank}_{step}"
            if cid in done:
                refine_rows.append(done[cid])
                continue
            center = np.array(base["center_values_m"], dtype=float)
            depth = np.array(base["depth_values_j"], dtype=float)
            physics = physics_from_config(cfg.base)
            sigma_c = ml.perturbation_sigma * physics["aod_initial_center_m"]
            sigma_d = ml.perturbation_sigma * physics["aod_initial_depth_j"]
            center[1:-1] += rng.normal(0, sigma_c, size=center.size - 2)
            depth[1:-1] += rng.normal(0, sigma_d, size=depth.size - 2)
            row = {"candidate_id": cid, "origin": f"refine_of_{base['candidate_id']}",
                   "center_values_m": center.tolist(), "depth_values_j": depth.tolist()}
            try:
                row.update(evaluate_ml_candidate(cfg, center, depth, pool, heating))
            except Exception as exc:
                row.update({"status": "failed", "reasons": repr(exc)})
            refine_rows.append(row)
            checkpoint_path.write_text(json.dumps({"rows": rows + refine_rows},
                                                  ensure_ascii=False))
    return rows + refine_rows


def select_ml_winner(cfg: Level2cConfig, rows: list[dict]) -> dict:
    """按（存活率, 激发, 光滑）字典序选唯一获胜候选并冻结。"""
    ok = [r for r in rows if r.get("status") == "ok"]
    if not ok:
        raise RuntimeError("ML 优化无有效候选")
    ok.sort(key=lambda r: (-r["survival_at_objective"],
                           r.get("mean_final_excitation_alive") or 1e9,
                           r["max_abs_center_acceleration_m_per_s2"]))
    winner = ok[0]
    return {
        "type": "ml_cubic",
        "duration_us": cfg.pickup.ml_duration_us,
        "control_points": cfg.pickup.ml_control_points,
        "center_values_m": winner["center_values_m"],
        "depth_values_j": winner["depth_values_j"],
        "provenance": {
            "candidate_id": winner["candidate_id"],
            "origin": winner["origin"],
            "training_survival_at_objective": winner["survival_at_objective"],
            "training_pool_seed": cfg.ml_optimization.candidate_seed,
            "transfers_for_objective": cfg.ml_optimization.transfers_for_objective,
            "objective": "survival_first_then_excitation_then_smoothness",
            "optimizer": "LHS + local perturbation refine（替代论文在线 ML 优化器，等价黑盒搜索）",
        },
    }
