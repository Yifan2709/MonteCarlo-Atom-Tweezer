"""Level 2B 实验 3：失效边界附近恒速 vs 变速剖面的双口径配对比较。

两种公平口径分别报告，缺一不可：
- same_average_speed：所有剖面 v_avg 相同（梯形峰值更低）；
- same_peak_speed：所有剖面 v_peak 相同（梯形平均速度更低）。
每个点先在扫描样本上比较，再用从未使用过的 boundary 种子独立复核。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from level2_joint_transfer.statistics import (paired_bootstrap_difference,
                                              paired_counts)
from .speed_engine import evaluate_profile_parallel
from .speed_scan import (FROZEN_MOVE_WINDOW, FROZEN_RAMP_POWER,
                         FROZEN_RAMP_WINDOW, _append_csv_row, load_flags)
from .speed_waveforms import build_speed_waveform, max_step_for, snap_dt, speed_calibers

COMPARISON_COLUMNS = [
    "dataset", "caliber", "temperature_uK", "v_target_mm_per_s", "profile",
    "duration_us", "dt_ns", "v_avg_mm_per_s", "v_peak_mm_per_s", "v_rms_mm_per_s",
    "shots", "captured", "capture_rate",
    "mean_excitation_captured", "margin_q10", "follow_loss_rate",
    "pair_vs_reference_difference", "pair_ci_low", "pair_ci_high",
    "both_captured", "only_profile", "only_reference", "both_failed",
    "seed", "status", "error",
]


def select_reference_curve(critical_table: list[dict]) -> dict:
    """选取确定比较速度的参考曲线：优先边界已覆盖的最高温度曲线。

    critical 表中 point_estimate 与 v_avg_* 的单位约定：
    point_estimate 已是 m/s；v_avg_max_mm_per_s 为 mm/s。
    """
    covered = [e for e in critical_table if e.get("boundary_covered")]
    if covered:
        chosen = max(covered, key=lambda e: e["temperature_uK"])
        boot = chosen.get("v95_bootstrap", {})
        point = boot.get("point_estimate")
        method = "bootstrap_interpolation"
        if point is None:
            point = chosen.get("v95_logistic", {}).get("point_estimate")
            method = "logistic_fit"
        if point is not None:
            return {"structure": chosen["structure"],
                    "temperature_uK": chosen["temperature_uK"],
                    "v95_m_per_s": float(point), "method": method,
                    "boundary_covered": True}
        fallback = chosen["v_avg_max_mm_per_s"]
        return {"structure": chosen["structure"],
                "temperature_uK": chosen["temperature_uK"],
                "v95_m_per_s": float(fallback) * 1e-3, "method": "fallback_vmax",
                "boundary_covered": True,
                "note": "v95 估计失败，退化为曲线最大速度"}
    fastest = max(critical_table, key=lambda e: e.get("v_avg_max_mm_per_s", 0.0)) \
        if critical_table else None
    if fastest is None:
        return {"v95_m_per_s": 0.12, "method": "default",
                "boundary_covered": False,
                "note": "无扫描曲线，使用配置默认参考速度"}
    return {"structure": fastest["structure"],
            "temperature_uK": fastest["temperature_uK"],
            "v95_m_per_s": float(fastest["v_avg_max_mm_per_s"]) * 1e-3,
            "method": "fallback_vmax_uncovered",
            "boundary_covered": False,
            "note": "全部曲线未覆盖边界，用最大扫描速度作参考并在报告中声明"}


def _profile_spec(profile: str, v_target_m_per_s: float, caliber: str,
                  distance_m: float, accel_fraction: float) -> dict:
    """按口径构造剖面 spec；返回 (spec, T_total)。

    same_average：T_move = d/v_target（所有剖面 v_avg = v_target）；
    same_peak：梯形 v_cruise = v_target；smootherstep v_peak = 1.875·v_avg
    = v_target → T_move = 1.875·d/v_target。
    """
    if profile == "trapezoid":
        if caliber == "same_average_speed":
            move_s = distance_m / v_target_m_per_s  # v_avg = v_target
        else:
            move_s = distance_m / ((1.0 - accel_fraction) * v_target_m_per_s)  # 巡航 = v_target
        duration = move_s / (FROZEN_MOVE_WINDOW[1] - FROZEN_MOVE_WINDOW[0])
        return {
            "type": "speed", "duration_s": duration,
            "move_start_fraction": FROZEN_MOVE_WINDOW[0],
            "move_end_fraction": FROZEN_MOVE_WINDOW[1],
            "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
            "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
            "ramp_power": FROZEN_RAMP_POWER,
            "move_profile": "trapezoid", "accel_fraction": accel_fraction,
            "depth_mode": "ramp_window",
        }
    if profile == "trapezoid_fast_acc":
        f = max(0.10, accel_fraction / 2.0)
        if caliber == "same_peak_speed":
            move_s = distance_m / ((1.0 - f) * v_target_m_per_s)
        else:
            move_s = distance_m / v_target_m_per_s
        duration = move_s / (FROZEN_MOVE_WINDOW[1] - FROZEN_MOVE_WINDOW[0])
        return {
            "type": "speed", "duration_s": duration,
            "move_start_fraction": FROZEN_MOVE_WINDOW[0],
            "move_end_fraction": FROZEN_MOVE_WINDOW[1],
            "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
            "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
            "ramp_power": FROZEN_RAMP_POWER,
            "move_profile": "trapezoid", "accel_fraction": f,
            "depth_mode": "ramp_window",
        }
    if profile == "frozen_smootherstep":
        if caliber == "same_peak_speed":
            move_s = 1.875 * distance_m / v_target_m_per_s
        else:
            move_s = distance_m / v_target_m_per_s
        duration = move_s / (FROZEN_MOVE_WINDOW[1] - FROZEN_MOVE_WINDOW[0])
        return {
            "type": "overlapped", "duration_s": duration,
            "move_start_fraction": FROZEN_MOVE_WINDOW[0],
            "move_end_fraction": FROZEN_MOVE_WINDOW[1],
            "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
            "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
            "ramp_power": FROZEN_RAMP_POWER,
        }
    raise ValueError(f"未知比较剖面 {profile!r}")


def _dt_for(spec: dict, physics: dict, base_dt_s: float, bcfg: dict) -> float:
    """按剖面实际峰值速度与移动时长取自适应步长。"""
    waveform = build_speed_waveform(spec, physics)
    cal = speed_calibers(waveform)
    dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                          physics["aod_waist_m"], base_dt_s,
                          min_steps_per_move=bcfg.get("min_steps_per_move", 2000),
                          max_step_displacement_over_waist=bcfg.get(
                              "max_step_displacement_over_waist", 0.01))
    return snap_dt(waveform.duration_s, dt_max)


def run_comparison(cfg2, bcfg: dict, physics: dict, critical_table: list[dict],
                   out_csv: Path, num_procs=None) -> dict:
    """执行实验 3：双口径 × 多速度点 × 三剖面的配对比较与独立复核。"""
    reference = select_reference_curve(critical_table)
    multipliers = [float(m) for m in bcfg.get("points_near_boundary",
                                              [0.5, 0.75, 1.0, 1.5, 2.0])]
    profiles = list(bcfg.get("profiles", ["trapezoid", "frozen_smootherstep",
                                          "trapezoid_fast_acc"]))
    shots = int(bcfg.get("comparison_shots", 400))
    temperature = float(reference.get("temperature_uK", bcfg["scan_temperatures_uK"][-1]))
    from .speed_engine import sample_states
    scan_states = sample_states(cfg2, int(bcfg["scan_seed"]) + 777, shots,
                                temperature_uK=temperature)
    fresh_states = sample_states(cfg2, int(bcfg["boundary_seed"]), shots,
                                 temperature_uK=temperature)
    distance = physics["aod_initial_center_m"]
    accel = float(bcfg["trapezoid_accel_fraction"])
    v_ref = float(reference["v95_m_per_s"])
    # 防御性夹取：参考速度必须落在物理合理的扫描量级内，防止估计异常
    # 生成荒谬时长（此前 1e-5 m/s 曾导致 0.29 s 转移与 5.8e6 步网格）。
    v_ref = min(max(v_ref, 2.0e-3), 0.5)
    reference = dict(reference, v95_m_per_s=v_ref)
    follow_cfg = bcfg["follow_loss"]
    bootstrap_cfg = int(bcfg.get("bootstrap_samples", 5000))

    for dataset, states, seed in (("scan_pool", scan_states, int(bcfg["scan_seed"]) + 777),
                                  ("boundary_reshot", fresh_states, int(bcfg["boundary_seed"]))):
        for caliber in ("same_average_speed", "same_peak_speed"):
            for multiplier in multipliers:
                v_target = v_ref * multiplier
                results = {}
                for profile in profiles:
                    try:
                        spec = _profile_spec(profile, v_target, caliber, distance, accel)
                        dt = _dt_for(spec, physics, bcfg["base_dt_us"] * 1e-6, bcfg)
                        out = evaluate_profile_parallel(
                            spec, states, physics, dt, hold_s=0.0,
                            follow_loss_sep_w=follow_cfg["max_separation_over_waist"],
                            follow_loss_excitation=follow_cfg["max_excitation_over_depth"],
                            work_energy_shots=2, num_procs=num_procs)
                        cal = speed_calibers(build_speed_waveform(spec, physics))
                        exc_cap = [r["final_excitation_over_depth"] for r in out["records"]
                                   if r["captured"]]
                        results[profile] = {
                            "spec": spec, "out": out, "cal": cal,
                            "mean_exc": float(np.mean(exc_cap)) if exc_cap else None,
                            "error": None,
                        }
                    except Exception as exc:
                        results[profile] = {"error": f"{type(exc).__name__}: {exc}"}
                reference_profile = "trapezoid"
                ref_ok = results.get(reference_profile, {}).get("out") is not None
                for profile in profiles:
                    entry = results[profile]
                    row = {"dataset": dataset, "caliber": caliber,
                           "temperature_uK": temperature,
                           "v_target_mm_per_s": v_target * 1e3, "profile": profile,
                           "seed": seed,
                           "status": "failed" if entry.get("error") else "ok",
                           "error": entry.get("error") or ""}
                    if entry.get("out") is not None:
                        s, cal = entry["out"]["summary"], entry["cal"]
                        row.update({
                            "duration_us": entry["spec"]["duration_s"] * 1e6,
                            "dt_ns": _dt_for(entry["spec"], physics,
                                             bcfg["base_dt_us"] * 1e-6, bcfg) * 1e9,
                            "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                            "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                            "v_rms_mm_per_s": cal["v_rms_m_per_s"] * 1e3,
                            "shots": s["shots"], "captured": s["captured"],
                            "capture_rate": s["capture_rate"],
                            "mean_excitation_captured": entry["mean_exc"] or "",
                            "margin_q10": s["margin_quantile"],
                            "follow_loss_rate": s["follow_loss_rate"],
                        })
                        if ref_ok and profile != reference_profile:
                            other = results[reference_profile]["out"]["captured_flags"]
                            mine = entry["out"]["captured_flags"]
                            counts = paired_counts(mine, other)
                            profile_seed = sum(ord(ch) for ch in profile)
                            boot = paired_bootstrap_difference(
                                mine, other, bootstrap_cfg,
                                seed + int(multiplier * 1000) + profile_seed)
                            row.update({
                                "pair_vs_reference_difference": boot["difference"],
                                "pair_ci_low": boot["ci_low"],
                                "pair_ci_high": boot["ci_high"],
                                "both_captured": counts["both_captured"],
                                "only_profile": counts["only_a_captured"],
                                "only_reference": counts["only_b_captured"],
                                "both_failed": counts["both_failed"],
                            })
                    _append_csv_row(out_csv, row, fieldnames=COMPARISON_COLUMNS)
    return {"reference": reference, "profiles": profiles,
            "multipliers": multipliers, "shots": shots,
            "temperature_uK": temperature}
