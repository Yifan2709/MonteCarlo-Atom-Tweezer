"""Level 2B 实验 5：默认噪声下的最高速度测算与同时长移动方案对比。

5a 硬件条件扫描：初始噪声模型（×1 默认幅度）开启下，改变 AOD 深度
   （激光功率）与 AOD 束腰（光斑大小）两类硬件条件，测端到端 v95
   （frozen_shape_scaled 结构），并换算 η95 = v_peak/v_ad。
5b 同时长剖面对比：在存活率陡降区取若干总时长，位移与"AOD 开始移动
   到最终停下"的时长完全相同、仅速度形状不同的四种移动方案
   （smootherstep 连续渐变 / 梯形 0.30 / 快加速梯形 0.15 / 三角 0.50
   匀加速-匀减速），比较俘获率与相对标准梯形的配对差。

约定与实验 4 一致：噪声为系综平均功率的确定性注入（无随机性），
幅度参数 assumed_sensitivity_only；初态统一用标称硬件采样（CRN，
温度远小于阱深时初态分布对硬件缩放不敏感）。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from level1_transfer_1d.analysis import wilson_interval
from level0_static_trap.potentials import trap_frequency_analytic
from level2_joint_transfer.statistics import paired_bootstrap_difference, paired_counts
from .critical_speed import (boundary_covered, bootstrap_interpolation_speeds,
                             logistic_fit_speeds)
from .noise_speed_validation import _dt_for_spec, _move_fraction, _spec_for, build_base_heating
from .speed_engine import evaluate_profile_parallel, sample_states
from .speed_scan import _append_csv_row, _save_flags

HARDWARE_COLUMNS = [
    "experiment", "condition", "aod_depth_uK", "aod_waist_um",
    "temperature_uK", "duration_us", "v_avg_mm_per_s", "v_peak_mm_per_s",
    "eta_peak", "shots", "captured", "capture_rate", "wilson_low", "wilson_high",
    "margin_q10", "margin_mean", "noise_total_injected_mean",
    "status", "error",
]

PROFILE_COLUMNS = [
    "experiment", "duration_us", "profile", "accel_fraction",
    "v_avg_mm_per_s", "v_peak_mm_per_s", "v_peak_over_avg", "eta_peak",
    "shots", "captured", "capture_rate", "wilson_low", "wilson_high",
    "margin_q10", "margin_mean", "follow_loss_rate",
    "noise_total_injected_mean",
    "pair_vs_reference_difference", "pair_ci_low", "pair_ci_high",
    "both_captured", "only_profile", "only_reference", "both_failed",
    "status", "error",
]


def condition_physics(base: dict, depth_scale: float | None = None,
                      waist_um: float | None = None) -> dict:
    """克隆物理字典并覆盖硬件条件（AOD 深度倍数 / AOD 束腰绝对值）。"""
    physics = dict(base)
    if depth_scale is not None:
        physics["aod_initial_depth_j"] = base["aod_initial_depth_j"] * float(depth_scale)
    if waist_um is not None:
        physics["aod_waist_m"] = float(waist_um) * 1e-6
    return physics


def adiabatic_speed_m_per_s(physics: dict) -> float:
    """本条件的绝热速度尺度 v_ad = w_AOD·ω_AOD。"""
    omega, _ = trap_frequency_analytic(physics["aod_initial_depth_j"],
                                       physics["aod_waist_m"], physics["mass_kg"])
    return physics["aod_waist_m"] * omega


def profile_spec_set(duration_s: float, accel_by_profile: dict) -> dict:
    """同一 duration 下四种移动方案的 spec（共享 frozen 移动/斜坡窗口）。

    各方案位移与移动时长完全一致，仅速度形状不同：
    - frozen_smootherstep：C² 连续渐变（v_peak = 1.875·v_avg）；
    - trapezoid f：加速-匀速-减速（v_peak = v_avg/(1−f)）；
      f=0.30 标准梯形、f=0.15 快加速梯形、f=0.50 三角（匀加速→匀减速）。
    """
    from .speed_scan import FROZEN_MOVE_WINDOW, FROZEN_RAMP_POWER, FROZEN_RAMP_WINDOW
    specs = {}
    for name, accel in accel_by_profile.items():
        if accel is None:  # frozen_smootherstep
            specs[name] = {
                "type": "overlapped", "duration_s": duration_s,
                "move_start_fraction": FROZEN_MOVE_WINDOW[0],
                "move_end_fraction": FROZEN_MOVE_WINDOW[1],
                "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
                "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
                "ramp_power": FROZEN_RAMP_POWER,
            }
        else:
            specs[name] = {
                "type": "speed", "duration_s": duration_s,
                "move_start_fraction": FROZEN_MOVE_WINDOW[0],
                "move_end_fraction": FROZEN_MOVE_WINDOW[1],
                "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
                "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
                "ramp_power": FROZEN_RAMP_POWER,
                "move_profile": "trapezoid", "accel_fraction": float(accel),
                "depth_mode": "ramp_window",
            }
    return specs


def run_max_speed_stage(cfg2, bcfg: dict, physics: dict, out_dir: Path,
                        num_procs=None) -> dict:
    """执行实验 5（5a 硬件条件 v95 + 5b 同时长剖面对比）并写落盘产物。"""
    probe = bcfg["max_speed_probe"]
    integration = bcfg["integration"]
    follow_cfg = bcfg["follow_loss"]
    base_dt = float(integration["base_dt_us"]) * 1e-6
    accel = float(bcfg["speed_scan"]["trapezoid_accel_fraction"])
    shots = int(probe.get("shots_per_point", 300))
    temperature = float(probe["temperature_uK"])
    structure = probe.get("v95_structure", "frozen_shape_scaled")
    durations_us = [float(t) for t in probe.get("v95_durations_us",
                                                [300, 220, 160, 120, 90, 65, 45, 32])]
    seed = int(bcfg["initial_ensemble"]["scan_seed"]) + 555

    hardware_csv = out_dir / "max_speed_hardware.csv"
    profiles_csv = out_dir / "max_speed_profiles.csv"
    new_flags = {}
    states = sample_states(cfg2, seed, shots, temperature_uK=temperature)

    def evaluate(spec, phys, heating):
        dt, waveform, cal = _dt_for_spec(spec, phys, base_dt, integration)
        result = evaluate_profile_parallel(
            spec, states, phys, dt, hold_s=0.0,
            follow_loss_sep_w=follow_cfg["max_separation_over_waist"],
            follow_loss_excitation=follow_cfg["max_excitation_over_depth"],
            work_energy_shots=integration.get("work_energy_shots", 4),
            num_procs=num_procs, heating=heating)
        return result, dt, waveform, cal

    # ------------------------------------------------ 5a 硬件条件 v95（×1 噪声）
    conditions = [{"name": "nominal", "aod_depth_scale": None, "aod_waist_um": None}]
    conditions += [dict(c) for c in probe.get("hardware_conditions", [])]
    move_frac = _move_fraction(structure)
    v95_table = []
    curves = {}
    for cond in conditions:
        phys = condition_physics(physics, cond.get("aod_depth_scale"),
                                 cond.get("aod_waist_um"))
        heating = build_base_heating(bcfg, phys)
        v_ad = adiabatic_speed_m_per_s(phys)
        curve = []
        for t_us in durations_us:
            try:
                spec = _spec_for(structure, t_us * 1e-6, None, accel, phys)
                result, dt, waveform, cal = evaluate(spec, phys, heating)
                s = result["summary"]
                captured = s["captured"]
                low, high = wilson_interval(captured, shots)
                margins = [r["capture_margin"] for r in result["records"]]
                row = {
                    "experiment": "hardware_v95", "condition": cond["name"],
                    "aod_depth_uK": phys["aod_initial_depth_j"] / 1.380649e-29,
                    "aod_waist_um": phys["aod_waist_m"] * 1e6,
                    "temperature_uK": temperature, "duration_us": t_us,
                    "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                    "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                    "eta_peak": cal["v_peak_m_per_s"] / v_ad,
                    "shots": shots, "captured": captured,
                    "capture_rate": s["capture_rate"],
                    "wilson_low": float(low), "wilson_high": float(high),
                    "margin_q10": s["margin_quantile"],
                    "margin_mean": float(np.mean(margins)),
                    "noise_total_injected_mean": s.get("mean_noise_total_energy_uK", ""),
                    "status": "ok", "error": "",
                }
                key = f"maxspeed_hw|{cond['name']}|{t_us:g}"
                new_flags[key] = result["captured_flags"]
                curve.append(row)
            except Exception as exc:
                row = {"experiment": "hardware_v95", "condition": cond["name"],
                       "temperature_uK": temperature, "duration_us": t_us,
                       "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
            _append_csv_row(hardware_csv, row, fieldnames=HARDWARE_COLUMNS)
        curves[cond["name"]] = curve
        speeds = np.array([2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac)
                           for r in curve])
        flags = [np.asarray(new_flags[f"maxspeed_hw|{cond['name']}|{float(r['duration_us']):g}"],
                            dtype=bool) for r in curve]
        rates = [float(r["capture_rate"]) for r in curve]
        entry = {
            "condition": cond["name"],
            "aod_depth_uK": phys["aod_initial_depth_j"] / 1.380649e-29,
            "aod_waist_um": phys["aod_waist_m"] * 1e6,
            "adiabatic_speed_mm_per_s": v_ad * 1e3,
            "boundary_covered": boundary_covered(rates, 0.95),
        }
        if len(curve) >= 3 and all(f.size for f in flags):
            entry["v95_bootstrap"] = bootstrap_interpolation_speeds(
                speeds, flags, 0.95, num_bootstrap=2000, seed=seed)
            entry["v95_logistic"] = logistic_fit_speeds(
                speeds, flags, 0.95, num_bootstrap=200, seed=seed + 1)
            point = entry["v95_bootstrap"].get("point_estimate")
            if point is not None:
                entry["eta95_peak"] = point * 1.875 / v_ad  # frozen 形状 v_peak=1.875·v_avg
        v95_table.append(entry)

    # ------------------------------------------------ 5b 同时长剖面对比（×1 噪声）
    heating = build_base_heating(bcfg, physics)
    v_ad = adiabatic_speed_m_per_s(physics)
    accel_by_profile = {
        "smootherstep": None,
        "trapezoid_30": 0.30,
        "trapezoid_15": 0.15,
        "triangular_50": 0.50,
    }
    reference = "trapezoid_30"
    bootstrap_cfg = int(bcfg["comparison"]["bootstrap_samples"])
    profiles_table = []
    for t_us in [float(t) for t in probe.get("profile_durations_us", [44, 39, 35])]:
        specs = profile_spec_set(t_us * 1e-6, accel_by_profile)
        results = {}
        for name, spec in specs.items():
            try:
                result, dt, waveform, cal = evaluate(spec, physics, heating)
                results[name] = {"result": result, "cal": cal}
            except Exception as exc:
                results[name] = {"error": f"{type(exc).__name__}: {exc}"}
        ref = results.get(reference, {}).get("result")
        for name in accel_by_profile:
            entry = results[name]
            error = entry.get("error")
            row = {"experiment": "same_duration_profiles", "duration_us": t_us,
                   "profile": name,
                   "accel_fraction": accel_by_profile[name] if accel_by_profile[name] is not None else "",
                   "status": "failed" if error else "ok", "error": error or ""}
            summary_item = None
            if not error:
                result, cal = entry["result"], entry["cal"]
                s = result["summary"]
                low, high = wilson_interval(s["captured"], s["shots"])
                margins = [r["capture_margin"] for r in result["records"]]
                row.update({
                    "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                    "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                    "v_peak_over_avg": cal["v_peak_m_per_s"] / cal["v_avg_m_per_s"],
                    "eta_peak": cal["v_peak_m_per_s"] / v_ad,
                    "shots": s["shots"], "captured": s["captured"],
                    "capture_rate": s["capture_rate"],
                    "wilson_low": float(low), "wilson_high": float(high),
                    "margin_q10": s["margin_quantile"],
                    "margin_mean": float(np.mean(margins)),
                    "follow_loss_rate": s["follow_loss_rate"],
                    "noise_total_injected_mean": s.get("mean_noise_total_energy_uK", ""),
                })
                if ref is not None and name != reference:
                    counts = paired_counts(result["captured_flags"], ref["captured_flags"])
                    boot = paired_bootstrap_difference(
                        result["captured_flags"], ref["captured_flags"],
                        bootstrap_cfg, seed + int(t_us * 10) + sum(ord(c) for c in name))
                    row.update({
                        "pair_vs_reference_difference": boot["difference"],
                        "pair_ci_low": boot["ci_low"], "pair_ci_high": boot["ci_high"],
                        "both_captured": counts["both_captured"],
                        "only_profile": counts["only_a_captured"],
                        "only_reference": counts["only_b_captured"],
                        "both_failed": counts["both_failed"],
                    })
                summary_item = {"profile": name,
                                "accel_fraction": accel_by_profile[name],
                                "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                                "v_peak_over_avg": cal["v_peak_m_per_s"] / cal["v_avg_m_per_s"],
                                "capture_rate": s["capture_rate"],
                                "wilson_low": float(low), "wilson_high": float(high),
                                "margin_mean": row["margin_mean"],
                                "follow_loss_rate": s["follow_loss_rate"],
                                "pair_vs_reference_difference":
                                    row.get("pair_vs_reference_difference"),
                                "pair_ci_low": row.get("pair_ci_low"),
                                "pair_ci_high": row.get("pair_ci_high")}
            _append_csv_row(profiles_csv, row, fieldnames=PROFILE_COLUMNS)
            if summary_item is not None:
                summary_item["duration_us"] = t_us
                profiles_table.append(summary_item)

    if new_flags:
        _save_flags(out_dir, new_flags)
    return {"config": {"shots_per_point": shots, "temperature_uK": temperature,
                       "v95_structure": structure,
                       "noise": "default_x1_linearized"},
            "v95_table": v95_table,
            "profiles": profiles_table}
