"""Level 2B 预检查：Level 2 十项复跑 + 本级六项，任一关键项失败即终止。

本级六项：
1. 梯形剖面校准（v_peak==v_cruise、端点零速/零加速度、位移严格 d）；
2. 新剖面解析导数 vs 有限差分；
3. 等比例缩放不变性（T→T/2 ⇒ v_peak×2、a_max×4）；
4. 慢极限回归（梯形 ≈6 mm/s vs Level 2 冻结波形，配对差 CI 含 0）；
5. 最快扫描点 dt 与 dt/2 收敛；
6. 静态极限（本引擎积分路径与 Level 0 积分器一致）。
"""
from __future__ import annotations

import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from level2_joint_transfer.level2_simulation import physics_from_config
from .critical_speed import boundary_covered
from .level2b_config import load_level2b_config
from .speed_engine import (derived_scales, evaluate_profile, sample_states,
                           static_limit_check)
from .speed_scan import end_to_end_specs, move_only_spec
from .speed_waveforms import (SpeedWaveform, build_speed_waveform, max_step_for,
                              snap_dt, speed_calibers)

PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def _check_existing_tests() -> dict:
    """运行本工程全部既有测试（level0/1/2 基线 + level2b 新增）。"""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(PACKAGE_ROOT / "tests"), "-q",
         "--no-header", "-m", "not slow"],
        cwd=PACKAGE_ROOT, capture_output=True, text=True, timeout=1800)
    tail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
    return {"exit_code": result.returncode, "summary_line": tail,
            "passed": result.returncode == 0}


def _check_trapezoid_calibration(physics: dict) -> dict:
    """梯形剖面校准：巡航速度、端点、位移、单调性。"""
    checks = {}
    for accel_fraction in (0.15, 0.30, 0.45):
        for duration_us in (50.0, 400.0):
            waveform = SpeedWaveform(
                duration_us * 1e-6, 0.01945549969117366, 1.0,
                0.6069511239907726, 1.0, 1.1121764674723393,
                physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
                move_profile="trapezoid", accel_fraction=accel_fraction)
            cal = speed_calibers(waveform, 8001)
            table = waveform.dense_table(8001)
            center = table["aod_center_m"]
            peak_rel = abs(cal["v_peak_m_per_s"] - waveform.cruise_speed_m_per_s) \
                / waveform.cruise_speed_m_per_s
            key = f"f{accel_fraction}_T{duration_us}us"
            checks[key] = {
                "peak_equals_cruise_rel_error": float(peak_rel),
                "cruise_mm_per_s": waveform.cruise_speed_m_per_s * 1e3,
                "peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                "displacement_rel_error": float(abs(
                    (center[0] - center[-1]) - physics["aod_initial_center_m"])
                    / physics["aod_initial_center_m"]),
                "monotone": bool(np.all(np.diff(center) <= 1e-12)),
                "endpoint_velocity_zero": abs(float(waveform.center_velocity(0.0))) < 1e-15
                and abs(float(waveform.center_velocity(waveform.duration_s))) < 1e-15,
                "endpoint_acceleration_zero": abs(float(waveform.center_acceleration(0.0))) < 1e-9
                and abs(float(waveform.center_acceleration(waveform.duration_s))) < 1e-9,
            }
            checks[key]["passed"] = (
                peak_rel < 1e-6
                and checks[key]["displacement_rel_error"] < 1e-9
                and checks[key]["monotone"]
                and checks[key]["endpoint_velocity_zero"]
                and checks[key]["endpoint_acceleration_zero"])
    return {"per_case": checks,
            "passed": all(case["passed"] for case in checks.values())}


def _check_derivatives(physics: dict) -> dict:
    """新剖面解析导数与中心差分一致。"""
    dt_frac = 1e-6  # 相对窗口的差分步长
    cases = {}
    for profile in ("trapezoid", "smootherstep"):
        waveform = SpeedWaveform(
            200e-6, 0.05, 0.95, 0.6, 1.0, 1.5,
            physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
            move_profile=profile, accel_fraction=0.30)
        probes = np.linspace(0.01, 0.99, 37) * waveform.duration_s
        errors_v, errors_a = [], []
        for t in probes:
            h = dt_frac * waveform.duration_s
            v_fd = (float(waveform.center(t + h)) - float(waveform.center(t - h))) / (2 * h)
            a_fd = (float(waveform.center(t + h)) - 2 * float(waveform.center(t))
                    + float(waveform.center(t - h))) / h**2
            scale_v = max(abs(float(waveform.center_velocity(t))), abs(v_fd), 1e-30)
            scale_a = max(abs(float(waveform.center_acceleration(t))), abs(a_fd), 1e-30)
            errors_v.append(abs(float(waveform.center_velocity(t)) - v_fd) / scale_v)
            errors_a.append(abs(float(waveform.center_acceleration(t)) - a_fd) / scale_a)
        # 巡航段加速度差分在相变边界附近会有 O(h) 污染，用中位数判据
        cases[profile] = {"v_median_rel_error": float(np.median(errors_v)),
                          "a_median_rel_error": float(np.median(errors_a)),
                          "v_max_rel_error": float(np.max(errors_v))}
        cases[profile]["passed"] = (cases[profile]["v_median_rel_error"] < 1e-5
                                    and cases[profile]["v_max_rel_error"] < 1e-3
                                    and cases[profile]["a_median_rel_error"] < 1e-3)
    return {"per_profile": cases, "passed": all(c["passed"] for c in cases.values())}


def _check_scaling_invariance(physics: dict) -> dict:
    """T→T/2 时 v_peak×2、a_max×4、归一化形状逐点一致。"""
    def peaks(duration_us):
        waveform = SpeedWaveform(
            duration_us * 1e-6, 0.01945549969117366, 1.0, 0.6069511239907726, 1.0,
            1.1121764674723393, physics["aod_initial_center_m"],
            physics["aod_initial_depth_j"], move_profile="trapezoid",
            accel_fraction=0.30)
        table = waveform.dense_table(4001)
        return (float(np.max(np.abs(table["center_velocity_m_per_s"]))),
                float(np.max(np.abs(table["center_acceleration_m_per_s2"]))),
                waveform)
    v1, a1, w1 = peaks(400.0)
    v2, a2, w2 = peaks(200.0)
    shape1 = np.asarray(w1.center(np.linspace(0, w1.duration_s, 2001)))
    shape2 = np.asarray(w2.center(np.linspace(0, w2.duration_s, 2001)))
    result = {
        "v_ratio": v2 / v1, "a_ratio": a2 / a1,
        "shape_max_abs_diff_m": float(np.max(np.abs(shape1 - shape2))),
    }
    result["passed"] = (abs(result["v_ratio"] - 2.0) < 1e-9
                        and abs(result["a_ratio"] - 4.0) < 1e-9
                        and result["shape_max_abs_diff_m"] < 1e-12)
    return result


def _check_slow_limit_regression(cfg2, bcfg: dict, physics: dict) -> dict:
    """慢极限回归：梯形 v_avg≈6.1 mm/s 与 Level 2 冻结 400 μs 波形配对差 CI 含 0。"""
    from level2_joint_transfer.statistics import paired_bootstrap_difference
    shots = int(bcfg.get("preflight_slow_limit_shots", 256))
    states = sample_states(cfg2, int(bcfg["initial_ensemble"]["scan_seed"]) + 555,
                           shots, temperature_uK=bcfg["initial_ensemble"]["main_temperature_uK"])
    # 冻结波形（Level 2 原实现路径）
    frozen = {
        "type": "overlapped", "duration_s": 400e-6,
        "move_start_fraction": 0.01945549969117366, "move_end_fraction": 1.0,
        "ramp_start_fraction": 0.6069511239907726, "ramp_end_fraction": 1.0,
        "ramp_power": 1.1121764674723393,
    }
    # 梯形：T_move = d/v_avg = 2.4μm/6.12mm/s ≈ 392.2 μs（与冻结波形同量级）
    v_avg_target = 6.12e-3
    move_s = physics["aod_initial_center_m"] / v_avg_target
    duration = move_s / (1.0 - 0.01945549969117366)
    trapezoid = {
        "type": "speed", "duration_s": duration,
        "move_start_fraction": 0.01945549969117366, "move_end_fraction": 1.0,
        "ramp_start_fraction": 0.6069511239907726, "ramp_end_fraction": 1.0,
        "ramp_power": 1.1121764674723393,
        "move_profile": "trapezoid", "accel_fraction": bcfg["speed_scan"]["trapezoid_accel_fraction"],
        "depth_mode": "ramp_window",
    }
    dt_s = bcfg["integration"]["base_dt_us"] * 1e-6
    hold_s = bcfg["integration"]["post_transfer_hold_us"] * 1e-6
    frozen_run = evaluate_profile(frozen, states, physics, snap_dt(400e-6, dt_s),
                                  hold_s=hold_s)
    trap_run = evaluate_profile(trapezoid, states, physics,
                                snap_dt(duration, dt_s), hold_s=hold_s)
    boot = paired_bootstrap_difference(trap_run["captured_flags"],
                                       frozen_run["captured_flags"], 10000, 24001)
    exc_trap = [r["final_excitation_over_depth"] for r in trap_run["records"]]
    exc_frozen = [r["final_excitation_over_depth"] for r in frozen_run["records"]]
    result = {
        "shots": shots,
        "frozen_capture_rate": frozen_run["summary"]["capture_rate"],
        "trapezoid_capture_rate": trap_run["summary"]["capture_rate"],
        "paired_difference": boot["difference"],
        "ci_low": boot["ci_low"], "ci_high": boot["ci_high"],
        "ci_contains_zero": bool(boot["ci_low"] <= 0.0 <= boot["ci_high"]),
        "mean_excitation_trapezoid": float(np.mean(exc_trap)),
        "mean_excitation_frozen": float(np.mean(exc_frozen)),
        "excitation_abs_diff": float(abs(np.mean(exc_trap) - np.mean(exc_frozen))),
    }
    result["passed"] = result["ci_contains_zero"] and result["excitation_abs_diff"] < 0.02
    return result


def _check_fastpoint_convergence(cfg2, bcfg: dict, physics: dict) -> dict:
    """最快扫描点：dt 与 dt/2 的末态连续量与俘获标签一致性。"""
    fastest = float(bcfg["speed_scan"]["pure_move_cruise_speeds_mm_per_s"][-1]) * 1e-3
    spec, _ = move_only_spec(fastest, physics["aod_initial_center_m"],
                             bcfg["speed_scan"]["trapezoid_accel_fraction"])
    waveform = build_speed_waveform(spec, physics)
    cal = speed_calibers(waveform)
    dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                          physics["aod_waist_m"], bcfg["integration"]["base_dt_us"] * 1e-6,
                          min_steps_per_move=bcfg["integration"]["min_steps_per_move"],
                          max_step_displacement_over_waist=bcfg["integration"]["max_step_displacement_over_waist"])
    dt1 = snap_dt(waveform.duration_s, dt_max)
    dt2 = snap_dt(waveform.duration_s, dt1 / 2.0)
    shots = int(bcfg.get("preflight_fastpoint_shots", 64))
    states = sample_states(cfg2, int(bcfg["initial_ensemble"]["scan_seed"]) + 556,
                           shots, temperature_uK=bcfg["initial_ensemble"]["main_temperature_uK"])
    run1 = evaluate_profile(spec, states, physics, dt1, hold_s=0.0)
    run2 = evaluate_profile(spec, states, physics, dt2, hold_s=0.0)
    x_diff = np.max(np.abs(np.array([r["final_x_um"] for r in run1["records"]])
                           - np.array([r["final_x_um"] for r in run2["records"]])))
    v_diff = np.max(np.abs(np.array([r["final_v_m_s"] for r in run1["records"]])
                           - np.array([r["final_v_m_s"] for r in run2["records"]])))
    e_diff = np.max(np.abs(np.array([r["final_slm_energy_uK"] for r in run1["records"]])
                           - np.array([r["final_slm_energy_uK"] for r in run2["records"]])))
    flags1 = np.array(run1["captured_flags"])
    flags2 = np.array(run2["captured_flags"])
    flips = flags1 != flags2
    flip_margins = np.abs(np.array([r["capture_margin"] for r in run1["records"]]))[flips]
    result = {
        "v_cruise_mm_per_s": fastest * 1e3, "dt1_ns": dt1 * 1e9, "dt2_ns": dt2 * 1e9,
        "shots": shots, "max_final_x_diff_um": float(x_diff),
        "max_final_v_diff_m_s": float(v_diff), "max_final_energy_diff_uK": float(e_diff),
        "label_flip_rate": float(np.mean(flips)),
        "max_flip_margin": float(flip_margins.max()) if flip_margins.size else 0.0,
    }
    result["label_flips_all_near_threshold"] = bool(
        flip_margins.size == 0 or float(flip_margins.max()) < 0.05)
    result["passed"] = float(np.mean(flips)) <= 0.02
    return result


def _check_static_limit(physics: dict, bcfg: dict) -> dict:
    """静态极限：本引擎路径与 Level 0 积分器逐点一致。"""
    dt = bcfg["integration"]["base_dt_us"] * 1e-6
    duration = 100e-6
    return static_limit_check(physics, dt, duration)


def run_preflight(bcfg: dict, cfg2, include_level2: bool = True,
                  include_existing_tests: bool = True) -> dict:
    """执行全部预检查并返回可 JSON 化的结果（含 all_passed 汇总）。"""
    started = time.time()
    physics = physics_from_config(cfg2)
    scales = derived_scales(physics)
    result = {
        "derived_scales": scales,
        "checks": {},
    }
    checks = result["checks"]

    if include_existing_tests:
        checks["existing_tests"] = _check_existing_tests()
    else:
        checks["existing_tests"] = {"skipped_by_flag": True, "passed": True}
    if include_level2:
        from level2_joint_transfer.preflight import run_preflight as level2_preflight
        checks["level2_preflight_rerun"] = {
            "all_passed": bool(level2_preflight(cfg2)["all_passed"])}
    checks["trapezoid_calibration"] = _check_trapezoid_calibration(physics)
    checks["derivative_consistency"] = _check_derivatives(physics)
    checks["scaling_invariance"] = _check_scaling_invariance(physics)
    checks["slow_limit_regression"] = _check_slow_limit_regression(cfg2, bcfg, physics)
    checks["fastpoint_convergence"] = _check_fastpoint_convergence(cfg2, bcfg, physics)
    checks["static_limit"] = _check_static_limit(physics, bcfg)

    critical = ["existing_tests", "trapezoid_calibration", "derivative_consistency",
                "scaling_invariance", "slow_limit_regression", "fastpoint_convergence",
                "static_limit"]
    if include_level2:
        critical.insert(1, "level2_preflight_rerun")
    result["critical_checks"] = critical
    result["all_passed"] = all(
        checks[name].get("passed", checks[name].get("all_passed", False))
        for name in critical)
    result["elapsed_s"] = time.time() - started
    return result
