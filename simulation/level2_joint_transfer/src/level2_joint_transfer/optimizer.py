"""分阶段波形优化：LHS 探索、局部精修、单调样条精修与冻结。

阶段 1 用 optimization_pool（默认 128 shots, dt=0.10 us）评估 LHS 候选；
阶段 2 用 selection_pool（默认 512 shots, dt=0.05 us）精修并选择；
最终波形冻结为 best_waveform.yaml，validation 池不参与本模块任何计算。
"""
from __future__ import annotations

import csv
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from .config import Level2Config
from .level2_simulation import build_waveform, evaluate_waveform_on_states
from .objectives import score_candidate
from .waveforms import validate_waveform

NPROC = max(1, (os.cpu_count() or 2) - 1)


def _evaluate_job(args):
    """子进程入口：评估一个候选（含约束检查），异常被记录为失败候选。"""
    (candidate_id, spec, states, physics, dt_s, hold_s, constraint_dict) = args
    row = {"candidate_id": candidate_id, **{k: spec.get(k) for k in (
        "move_start_fraction", "move_end_fraction", "ramp_start_fraction",
        "ramp_end_fraction", "ramp_power")}}
    if spec["type"] == "spline":
        row.update({"control_s": spec.get("control_s"),
                    "center_values_m": spec.get("center_values_m"),
                    "depth_values_j": spec.get("depth_values_j")})
    row.update({"waveform_type": spec["type"], "pool": spec.get("pool", ""),
                "dt_us": dt_s * 1e6, "status": "ok", "constraint_ok": True, "failure_reason": "",
                "is_finalist": False, "is_selected": False})
    try:
        waveform = build_waveform(spec, physics)
        if spec["type"] in ("overlapped", "spline"):
            valid, reasons = validate_waveform(waveform, constraint_dict)
            if not valid:
                row.update({"status": "invalid", "constraint_ok": False,
                            "failure_reason": "; ".join(reasons)})
                return row
        evaluation = evaluate_waveform_on_states(spec, states, physics, dt_s, hold_s)
        smoothness = waveform.smoothness_metrics()
        score = score_candidate(evaluation["summary"], smoothness, _opt_cfg_holder["cfg"])
        row.update({
            "capture_rate": score["capture_rate"], "captured_count": score["captured_count"],
            "shots": score["shots"], "margin_quantile": score["margin_quantile"],
            "excitation_term": score["excitation_term"], "smooth_penalty": score["smooth_penalty"],
            "objective_total": score["objective_total"],
            "mean_final_excitation_normalized_captured":
                evaluation["summary"]["mean_final_excitation_normalized_captured"],
            "max_abs_center_velocity_m_per_s": smoothness["max_abs_center_velocity_m_per_s"],
            "max_abs_center_acceleration_m_per_s2": smoothness["max_abs_center_acceleration_m_per_s2"],
            "max_abs_center_jerk_m_per_s3": smoothness["max_abs_center_jerk_m_per_s3"],
            "max_abs_depth_rate_J_per_s": smoothness["max_abs_depth_rate_J_per_s"],
        })
    except Exception as exc:  # 单个候选失败不得让整轮搜索丢失
        row.update({"status": "failed", "constraint_ok": False, "failure_reason": repr(exc)[:300]})
    return row


class _OptCfgHolder(dict):
    """子进程内缓存优化权重（由 Pool initializer 注入）。"""


_opt_cfg_holder = {"cfg": None}


def _init_worker(opt_cfg):
    _opt_cfg_holder["cfg"] = opt_cfg


def _config_signature(cfg: Level2Config) -> str:
    import hashlib
    payload = json.dumps({
        "opt": vars(cfg.optimization), "ens": vars(cfg.initial_ensemble),
        "wf": {k: v for k, v in vars(cfg.waveforms).items()},
        "integ": vars(cfg.integration),
    }, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def lhs_parameter_table(cfg: Level2Config) -> list[dict]:
    """Latin hypercube 生成受约束参数候选（窗口宽度与次序在映射中保证）。"""
    opt = cfg.optimization
    sampler = qmc.LatinHypercube(d=5, seed=opt.lhs_seed)
    unit = sampler.random(opt.max_candidates)
    min_width = cfg.waveforms.min_segment_fraction
    p_lo, p_hi = cfg.waveforms.ramp_power_bounds
    table = []
    for row in unit:
        u_ax, u_wx, u_ad, u_wd, u_p = row
        wx = min_width + u_wx * (1.0 - min_width)
        ax = u_ax * (1.0 - wx)
        wd = min_width + u_wd * (1.0 - min_width)
        ad = u_ad * (1.0 - wd)
        table.append({
            "move_start_fraction": ax, "move_end_fraction": ax + wx,
            "ramp_start_fraction": ad, "ramp_end_fraction": ad + wd,
            "ramp_power": p_lo + u_p * (p_hi - p_lo),
        })
    return table


def _row_sort_key(row):
    if row.get("status") != "ok":
        return (-1, ())
    return (0, (row["captured_count"], row["objective_total"]))


def run_stage1(cfg: Level2Config, states: dict, physics: dict, output_dir: Path,
               checkpoint_path: Path, history_rows: list) -> list[dict]:
    """阶段 1：LHS 候选探索（optimization_pool, dt=0.10 us），支持断点续跑。"""
    constraint = cfg.waveforms.as_constraint_dict()
    signature = _config_signature(cfg)
    done_ids = set()
    rows: list[dict] = []
    if checkpoint_path.exists():
        payload = json.loads(checkpoint_path.read_text())
        if payload.get("signature") == signature and "stage1" in payload.get("stages", {}):
            rows = payload["stages"]["stage1"]["rows"]
            done_ids = {r["candidate_id"] for r in rows}
    table = lhs_parameter_table(cfg)
    jobs = []
    for index, params in enumerate(table):
        candidate_id = f"lhs_{index:04d}"
        if candidate_id in done_ids:
            continue
        spec = dict(params)
        spec.update({"type": "overlapped", "duration_s": cfg.waveforms.optimized_duration_s,
                     "pool": "optimization"})
        jobs.append((candidate_id, spec, states, physics, cfg.integration.optimization_dt_s,
                     cfg.integration.post_transfer_hold_s, constraint))
    if jobs:
        print(f"[stage1] evaluating {len(jobs)} candidates on {NPROC} procs ...")
        started = time.time()
        with Pool(NPROC, initializer=_init_worker, initargs=(cfg.optimization,)) as pool:
            new_rows = pool.map(_evaluate_job, jobs)
        rows.extend(new_rows)
        rows.sort(key=lambda r: r["candidate_id"])
        _write_checkpoint(checkpoint_path, signature, "stage1", rows, history_rows)
        print(f"[stage1] done in {time.time() - started:.1f}s")
    return rows


def _write_checkpoint(path: Path, signature: str, stage: str, rows: list, history_rows: list):
    payload = {"signature": signature, "stages": {}}
    if path.exists():
        try:
            payload = json.loads(path.read_text())
            if payload.get("signature") != signature:
                payload = {"signature": signature, "stages": {}}
        except json.JSONDecodeError:
            payload = {"signature": signature, "stages": {}}
    payload["stages"][stage] = {"rows": rows, "count": len(rows)}
    payload["history"] = history_rows
    path.write_text(json.dumps(payload, indent=1))


def _perturb(params: dict, sigma: float, rng) -> dict:
    """对窗口分数与幂做高斯扰动并裁剪回合法范围。"""
    move_width = params["move_end_fraction"] - params["move_start_fraction"]
    ramp_width = params["ramp_end_fraction"] - params["ramp_start_fraction"]
    move_start = float(np.clip(params["move_start_fraction"] + rng.normal(0, sigma), 0.0, 1.0 - move_width))
    ramp_start = float(np.clip(params["ramp_start_fraction"] + rng.normal(0, sigma), 0.0, 1.0 - ramp_width))
    return {
        "move_start_fraction": move_start,
        "move_end_fraction": move_start + move_width,
        "ramp_start_fraction": ramp_start,
        "ramp_end_fraction": ramp_start + ramp_width,
        "ramp_power": float(np.clip(params["ramp_power"] + rng.normal(0, max(0.1, sigma * 5)), 0.5, 4.0)),
    }


def run_stage2(cfg: Level2Config, stage1_rows: list[dict], states: dict, physics: dict,
               output_dir: Path, checkpoint_path: Path, history_rows: list) -> list[dict]:
    """阶段 2：入围者在 selection_pool 上精修（dt=0.05 us）。"""
    opt = cfg.optimization
    constraint = cfg.waveforms.as_constraint_dict()
    ok_rows = [r for r in stage1_rows if r.get("status") == "ok"]
    ok_rows.sort(key=_row_sort_key, reverse=True)
    finalists = ok_rows[:opt.finalists]
    for row in finalists:
        row["is_finalist"] = True
    rng = np.random.default_rng(opt.lhs_seed + 1)
    jobs = []
    for row in finalists:
        base = {k: row[k] for k in ("move_start_fraction", "move_end_fraction",
                                    "ramp_start_fraction", "ramp_end_fraction", "ramp_power")}
        jobs.append((f"sel_{row['candidate_id']}", dict(base), "selection"))
        for step in range(opt.refine_steps_per_finalist):
            jobs.append((f"ref_{row['candidate_id']}_{step}", _perturb(base, opt.refine_perturbation_sigma, rng),
                         "selection"))
    args = []
    for candidate_id, params, pool_name in jobs:
        spec = dict(params)
        spec.update({"type": "overlapped", "duration_s": cfg.waveforms.optimized_duration_s,
                     "pool": pool_name})
        args.append((candidate_id, spec, states, physics, cfg.integration.validation_dt_s,
                     cfg.integration.post_transfer_hold_s, constraint))
    print(f"[stage2] evaluating {len(args)} refined candidates ...")
    started = time.time()
    with Pool(NPROC, initializer=_init_worker, initargs=(opt,)) as pool:
        rows = pool.map(_evaluate_job, args)
    for row in rows:
        history_rows.append({"candidate_id": row["candidate_id"], "status": row["status"],
                             "objective_total": row.get("objective_total"),
                             "captured_count": row.get("captured_count")})
    rows.sort(key=_row_sort_key, reverse=True)
    best = rows[0]
    for row in rows:
        row["is_selected"] = row["candidate_id"] == best["candidate_id"]
    _write_checkpoint(checkpoint_path, _config_signature(cfg), "stage2", rows, history_rows)
    print(f"[stage2] done in {time.time() - started:.1f}s, best={best['candidate_id']} "
          f"(capture {best.get('capture_rate', float('nan')):.4f})")
    return rows


def run_spline_refinement(cfg: Level2Config, best_row: dict, states: dict, physics: dict,
                          checkpoint_path: Path, history_rows: list) -> list[dict]:
    """由低维最优波形初始化 8 控制点 PCHIP 单调样条，在 selection_pool 上精修。"""
    opt = cfg.optimization
    if not opt.enable_monotone_spline_refinement:
        return []
    base_spec = {"type": "overlapped", "duration_s": cfg.waveforms.optimized_duration_s,
                 "move_start_fraction": best_row["move_start_fraction"],
                 "move_end_fraction": best_row["move_end_fraction"],
                 "ramp_start_fraction": best_row["ramp_start_fraction"],
                 "ramp_end_fraction": best_row["ramp_end_fraction"],
                 "ramp_power": best_row["ramp_power"]}
    base_waveform = build_waveform(base_spec, physics)
    k = opt.spline_control_points
    control_s = np.linspace(0.0, 1.0, k)
    base_center = np.asarray(base_waveform.center(control_s * cfg.waveforms.optimized_duration_s))
    base_depth = np.asarray(base_waveform.depth(control_s * cfg.waveforms.optimized_duration_s))
    rng = np.random.default_rng(opt.lhs_seed + 2)
    constraint = cfg.waveforms.as_constraint_dict()

    def make_spec(center_values, depth_values):
        return {"type": "spline", "duration_s": cfg.waveforms.optimized_duration_s,
                "control_s": control_s.tolist(),
                "center_values_m": center_values.tolist(),
                "depth_values_j": depth_values.tolist(), "pool": "selection"}

    specs = [make_spec(base_center, base_depth)]
    scale_c = physics["aod_initial_center_m"] * 0.04
    scale_d = physics["aod_initial_depth_j"] * 0.04
    for index in range(opt.spline_candidates):
        center = base_center + rng.normal(0, scale_c, k)
        depth = base_depth + rng.normal(0, scale_d, k)
        center = np.minimum.accumulate(np.clip(center, 0.0, physics["aod_initial_center_m"]))
        depth = np.minimum.accumulate(np.clip(depth, 0.0, physics["aod_initial_depth_j"]))
        center[0], center[-1] = physics["aod_initial_center_m"], 0.0
        depth[0], depth[-1] = physics["aod_initial_depth_j"], 0.0
        center = np.maximum.accumulate(center[::-1])[::-1]  # 端点重设后再保证单调
        depth = np.maximum.accumulate(depth[::-1])[::-1]
        center[0], center[-1] = physics["aod_initial_center_m"], 0.0
        depth[0], depth[-1] = physics["aod_initial_depth_j"], 0.0
        specs.append(make_spec(center, depth))

    args = [(f"spl_{i:03d}", spec, states, physics, cfg.integration.validation_dt_s,
             cfg.integration.post_transfer_hold_s, constraint) for i, spec in enumerate(specs)]
    print(f"[spline] evaluating {len(args)} spline candidates ...")
    with Pool(NPROC, initializer=_init_worker, initargs=(opt,)) as pool:
        rows = pool.map(_evaluate_job, args)
    for row in rows:
        history_rows.append({"candidate_id": row["candidate_id"], "status": row["status"],
                             "objective_total": row.get("objective_total"),
                             "captured_count": row.get("captured_count")})
    rows.sort(key=_row_sort_key, reverse=True)
    _write_checkpoint(checkpoint_path, _config_signature(cfg), "spline", rows, history_rows)
    print(f"[spline] best spline capture "
          f"{rows[0].get('capture_rate', float('nan')):.4f}" if rows and rows[0].get("status") == "ok" else "[spline] no valid spline")
    return rows


def select_and_freeze(cfg: Level2Config, stage1_rows, stage2_rows, spline_rows, output_dir: Path) -> dict:
    """从阶段 2（含样条）中选择唯一最终波形并冻结为 best_waveform.yaml。"""
    candidates = [r for r in list(stage2_rows) + list(spline_rows) if r.get("status") == "ok"]
    if not candidates:
        raise RuntimeError("没有有效候选可以冻结")
    candidates.sort(key=_row_sort_key, reverse=True)
    best = candidates[0]
    for row in candidates:
        row["is_selected"] = row["candidate_id"] == best["candidate_id"]
    if best["waveform_type"] == "overlapped":
        frozen = {
            "type": "overlapped",
            "name": "overlapped_optimized_400us",
            "duration_us": float(cfg.waveforms.optimized_duration_us),
            "move_start_fraction": float(best["move_start_fraction"]),
            "move_end_fraction": float(best["move_end_fraction"]),
            "ramp_start_fraction": float(best["ramp_start_fraction"]),
            "ramp_end_fraction": float(best["ramp_end_fraction"]),
            "ramp_power": float(best["ramp_power"]),
            "provenance": {"candidate_id": best["candidate_id"], "pool": best.get("pool"),
                           "objective_total": None if best.get("objective_total") is None
                           else float(best["objective_total"]),
                           "selection_capture_rate": None if best.get("capture_rate") is None
                           else float(best["capture_rate"])},
        }
    else:
        frozen = {
            "type": "spline",
            "name": "overlapped_optimized_400us",
            "duration_us": float(cfg.waveforms.optimized_duration_us),
            "control_s": [float(v) for v in best.get("control_s")],
            "center_values_m": [float(v) for v in best.get("center_values_m")],
            "depth_values_j": [float(v) for v in best.get("depth_values_j")],
            "provenance": {"candidate_id": best["candidate_id"], "pool": best.get("pool"),
                           "objective_total": None if best.get("objective_total") is None
                           else float(best["objective_total"]),
                           "selection_capture_rate": None if best.get("capture_rate") is None
                           else float(best["capture_rate"])},
        }
    import yaml
    (output_dir / "best_waveform.yaml").write_text(
        yaml.safe_dump(frozen, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return frozen


def frozen_spec(frozen: dict, cfg: Level2Config) -> dict:
    """把冻结 YAML 转回评估 spec（供 validate/robustness 使用）。"""
    if frozen["type"] == "overlapped":
        return {"type": "overlapped", "name": frozen.get("name", "optimized"),
                "duration_s": frozen["duration_us"] * 1e-6,
                "move_start_fraction": frozen["move_start_fraction"],
                "move_end_fraction": frozen["move_end_fraction"],
                "ramp_start_fraction": frozen["ramp_start_fraction"],
                "ramp_end_fraction": frozen["ramp_end_fraction"],
                "ramp_power": frozen["ramp_power"]}
    return {"type": "spline", "name": frozen.get("name", "optimized"),
            "duration_s": frozen["duration_us"] * 1e-6,
            "control_s": frozen["control_s"],
            "center_values_m": frozen["center_values_m"],
            "depth_values_j": frozen["depth_values_j"]}


def write_candidates_csv(path: Path, rows: list[dict]):
    """保存全部候选（含无效与失败），不得只保存最优。"""
    fieldnames = ["candidate_id", "waveform_type", "move_start_fraction", "move_end_fraction",
                  "ramp_start_fraction", "ramp_end_fraction", "ramp_power", "pool", "dt_us",
                  "status", "constraint_ok", "failure_reason", "capture_rate", "captured_count",
                  "shots", "margin_quantile", "excitation_term", "smooth_penalty",
                  "objective_total", "mean_final_excitation_normalized_captured",
                  "max_abs_center_velocity_m_per_s", "max_abs_center_acceleration_m_per_s2",
                  "max_abs_center_jerk_m_per_s3", "max_abs_depth_rate_J_per_s",
                  "is_finalist", "is_selected"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})
