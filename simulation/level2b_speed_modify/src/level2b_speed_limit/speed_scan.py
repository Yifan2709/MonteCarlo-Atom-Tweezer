"""Level 2B 实验 1/2 编排：纯移动跟随扫描与端到端速度扫描。

每个扫描点完成即追加写入 CSV（断点续跑跳过已完成点）；曲线级早停
（俘获率低于阈值后跳过更快点）只作用于该曲线，不影响其他曲线。
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from level1_transfer_1d.analysis import wilson_interval
from level2_joint_transfer.config import load_config as load_level2_config
from .critical_speed import (boundary_covered, bootstrap_interpolation_speeds,
                             logistic_fit_speeds)
from .speed_engine import (derived_scales, evaluate_profile_parallel,
                           physics_from_config, sample_states)
from .speed_waveforms import (build_speed_waveform, max_step_for, snap_dt,
                              speed_calibers)

# Level 2 冻结波形参数（原项目 outputs/level2_joint_transfer/demo/best_waveform.yaml，
# candidate ref_lhs_0150_3；本目录不修改 level2 基线，故以字面量冻结并标注来源）。
FROZEN_MOVE_WINDOW = (0.01945549969117366, 1.0)
FROZEN_RAMP_WINDOW = (0.6069511239907726, 1.0)
FROZEN_RAMP_POWER = 1.1121764674723393

SCAN_COLUMNS = [
    "experiment", "structure", "temperature_uK", "target",
    "duration_us", "dt_ns", "v_avg_mm_per_s", "v_peak_mm_per_s", "v_rms_mm_per_s",
    "eta_peak", "step_displacement_frac", "move_duration_us",
    "shots", "captured", "capture_rate", "wilson_low", "wilson_high",
    "follow_loss_rate", "exc_move_end_median", "exc_move_end_q90", "margin_q10",
    "failure_follow_loss", "failure_ramp_excitation", "failure_handoff",
    "work_residual_median_over_depth", "status", "error",
]


def end_to_end_specs(duration_s: float, accel_fraction: float) -> dict:
    """实验 2 的三套窗口结构（共享同一 duration）。"""
    return {
        "sequential_scaled": {
            "type": "sequential", "duration_s": duration_s, "move_fraction": 0.52},
        "frozen_shape_scaled": {
            "type": "overlapped", "duration_s": duration_s,
            "move_start_fraction": FROZEN_MOVE_WINDOW[0],
            "move_end_fraction": FROZEN_MOVE_WINDOW[1],
            "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
            "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
            "ramp_power": FROZEN_RAMP_POWER},
        "trapezoid_frozen_ramp": {
            "type": "speed", "duration_s": duration_s,
            "move_start_fraction": FROZEN_MOVE_WINDOW[0],
            "move_end_fraction": FROZEN_MOVE_WINDOW[1],
            "ramp_start_fraction": FROZEN_RAMP_WINDOW[0],
            "ramp_end_fraction": FROZEN_RAMP_WINDOW[1],
            "ramp_power": FROZEN_RAMP_POWER,
            "move_profile": "trapezoid", "accel_fraction": accel_fraction,
            "depth_mode": "ramp_window"},
    }


def move_only_spec(v_cruise_m_per_s: float, distance_m: float, accel_fraction: float,
                   profile: str = "trapezoid") -> tuple[dict, float]:
    """实验 1 的纯移动剖面：深度恒定，目标速度按峰值速度参数化。

    trapezoid：v_cruise 即巡航=峰值速度，T_move = d/[(1−f)·v_cruise]；
    smootherstep：v_peak = 1.875·v_avg，T_move = d/(v_target/1.875)。
    """
    if profile == "trapezoid":
        move_duration_s = distance_m / ((1.0 - accel_fraction) * v_cruise_m_per_s)
    else:
        move_duration_s = distance_m * 1.875 / v_cruise_m_per_s
    spec = {
        "type": "speed", "duration_s": move_duration_s,
        "move_start_fraction": 0.0, "move_end_fraction": 1.0,
        "ramp_start_fraction": 0.90, "ramp_end_fraction": 1.0,  # constant 模式下不使用
        "ramp_power": FROZEN_RAMP_POWER,
        "move_profile": profile, "accel_fraction": accel_fraction,
        "depth_mode": "constant",
    }
    return spec, move_duration_s


def _append_csv_row(path: Path, row: dict, fieldnames: list[str] | None = None):
    """追加一行到扫描/比较 CSV（不存在则先写表头，字段名随既有表头）。"""
    import csv
    if fieldnames is None:
        if path.is_file():
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader, None)
            fieldnames = header if header else list(row)
        else:
            fieldnames = list(row)
    exists = path.is_file()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in fieldnames})


def _completed_keys(path: Path) -> set[tuple]:
    """已写入且状态非 failed 的扫描点键集合（用于断点续跑）。"""
    import csv
    if not path.is_file():
        return set()
    keys = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") in ("ok", "skipped_early_stop"):
                keys.add((row["experiment"], row["structure"],
                          float(row["temperature_uK"]), row["target"]))
    return keys


def _evaluate_point(spec: dict, states: dict, physics: dict, dt_s: float,
                    hold_s: float, scales: dict, follow_cfg: dict,
                    work_energy_shots: int, hold_dt_s: float = 5.0e-8,
                    num_procs=None) -> dict:
    """评估单个扫描点并返回汇总（不含 CSV 列名相关的胶水）。"""
    waveform = build_speed_waveform(spec, physics)
    cal = speed_calibers(waveform)
    result = evaluate_profile_parallel(
        spec, states, physics, dt_s, hold_s=hold_s, hold_dt_s=hold_dt_s,
        follow_loss_sep_w=follow_cfg["max_separation_over_waist"],
        follow_loss_excitation=follow_cfg["max_excitation_over_depth"],
        work_energy_shots=work_energy_shots, num_procs=num_procs)
    summary = result["summary"]
    captured, shots = summary["captured"], summary["shots"]
    low, high = wilson_interval(captured, shots)
    eta = cal["v_peak_m_per_s"] / scales["adiabatic_speed_m_per_s"]
    step_frac = dt_s * cal["v_peak_m_per_s"] / physics["aod_waist_m"]
    return {
        "summary": summary, "captured_flags": result["captured_flags"],
        "calibers": cal, "wilson": (float(low), float(high)),
        "eta_peak": eta, "step_displacement_frac": step_frac,
    }


def _save_flags(out_dir: Path, flags: dict):
    """把逐 shot 俘获旗标并入 scan_flags.npz（断点续跑后仍可做 bootstrap）。"""
    path = out_dir / "scan_flags.npz"
    merged = {}
    if path.is_file():
        with np.load(path) as archive:
            merged = {key: archive[key] for key in archive.files}
    merged.update({key: np.asarray(value, dtype=bool) for key, value in flags.items()})
    np.savez_compressed(path, **merged)


def load_flags(out_dir: Path) -> dict:
    """读取历史俘获旗标档案。"""
    path = out_dir / "scan_flags.npz"
    if not path.is_file():
        return {}
    with np.load(path) as archive:
        return {key: archive[key] for key in archive.files}


def _flag_key(experiment: str, structure: str, temperature: float, target) -> str:
    return f"{experiment}|{structure}|{temperature}|{target}"


def run_experiment1(cfg2, bcfg: dict, physics: dict, scales: dict,
                    out_csv: Path, out_dir: Path, num_procs=None) -> list[dict]:
    """实验 1：纯移动跟随上限（梯形巡航，深度恒定 280 μK，5 μK 初态）。

    转移段按自适应 dt 规则积分，保持段（100 μs）用固定 50 ns 继续，
    两段均为同一种 velocity-Verlet。
    """
    speeds = [float(v) * 1e-3 for v in bcfg["pure_move_cruise_speeds_mm_per_s"]]
    shots = int(bcfg.get("pure_move_shots", bcfg["scan_shots_per_point"]))
    states = sample_states(cfg2, int(bcfg["scan_seed"]), shots,
                           temperature_uK=bcfg["main_temperature_uK"])
    distance = physics["aod_initial_center_m"]
    follow_cfg = bcfg["follow_loss"]
    done = _completed_keys(out_csv)
    rows = []
    new_flags = {}
    for v_cruise in speeds:
        target = f"{v_cruise * 1e3:.6g}"
        key = ("pure_move", "trapezoid_cruise", bcfg["main_temperature_uK"], target)
        if key in done:
            continue
        spec, _ = move_only_spec(v_cruise, distance, bcfg["trapezoid_accel_fraction"])
        waveform = build_speed_waveform(spec, physics)
        cal = speed_calibers(waveform)
        dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                              physics["aod_waist_m"], bcfg["base_dt_us"] * 1e-6,
                              min_steps_per_move=bcfg.get("min_steps_per_move", 2000),
                              max_step_displacement_over_waist=bcfg.get(
                                  "max_step_displacement_over_waist", 0.01))
        dt = snap_dt(waveform.duration_s, dt_max)
        try:
            point = _evaluate_point(spec, states, physics, dt,
                                    bcfg["post_transfer_hold_us"] * 1e-6, scales,
                                    follow_cfg, bcfg.get("work_energy_shots", 4),
                                    num_procs=num_procs)
            s, cal, wilson = point["summary"], point["calibers"], point["wilson"]
            row = {
                "experiment": "pure_move", "structure": "trapezoid_cruise",
                "temperature_uK": bcfg["main_temperature_uK"],
                "target": target,
                "duration_us": waveform.duration_s * 1e6, "dt_ns": dt * 1e9,
                "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                "v_rms_mm_per_s": cal["v_rms_m_per_s"] * 1e3,
                "eta_peak": point["eta_peak"],
                "step_displacement_frac": point["step_displacement_frac"],
                "move_duration_us": cal["move_duration_s"] * 1e6,
                "shots": s["shots"], "captured": s["captured"],
                "capture_rate": s["capture_rate"],
                "wilson_low": wilson[0], "wilson_high": wilson[1],
                "follow_loss_rate": s["follow_loss_rate"],
                "exc_move_end_median": s["excitation_rel_move_end_median"],
                "exc_move_end_q90": s["excitation_rel_move_end_q90"],
                "margin_q10": s["margin_quantile"],
                "failure_follow_loss": s["failure_modes"]["follow_loss"],
                "failure_ramp_excitation": s["failure_modes"]["ramp_excitation"],
                "failure_handoff": s["failure_modes"]["handoff_failure"],
                "work_residual_median_over_depth": s.get("work_residual_median_over_depth", ""),
                "status": "ok", "error": "",
            }
            new_flags[_flag_key("pure_move", "trapezoid_cruise",
                                bcfg["main_temperature_uK"], target)] = \
                point["captured_flags"]
        except Exception as exc:  # 单点异常记录为失败点，不中断整轮
            row = {"experiment": "pure_move", "structure": "trapezoid_cruise",
                   "temperature_uK": bcfg["main_temperature_uK"],
                   "target": target, "status": "failed",
                   "error": f"{type(exc).__name__}: {exc}"}
        _append_csv_row(out_csv, row, fieldnames=SCAN_COLUMNS)
        rows.append(row)
    if new_flags:
        _save_flags(out_dir, new_flags)
    return rows


def run_experiment2(cfg2, bcfg: dict, physics: dict, scales: dict,
                    out_csv: Path, out_dir: Path, num_procs=None) -> list[dict]:
    """实验 2：端到端速度上限（三结构 × 三温度 × 时长网格，CRN，曲线级早停）。"""
    durations = [float(t) * 1e-6 for t in bcfg["end_to_end_durations_us"]]
    temperatures = [float(t) for t in bcfg["scan_temperatures_uK"]]
    shots = int(bcfg["scan_shots_per_point"])
    stop_below = float(bcfg.get("stop_when_capture_below", 0.90))
    follow_cfg = bcfg["follow_loss"]
    done = _completed_keys(out_csv)
    rows = []
    new_flags = {}
    for temperature in temperatures:
        states = sample_states(cfg2, int(bcfg["scan_seed"]) + int(temperature * 1000),
                               shots, temperature_uK=temperature)
        for structure in ("sequential_scaled", "frozen_shape_scaled", "trapezoid_frozen_ramp"):
            early_stop = False
            for duration in durations:
                target = f"{duration * 1e6:.6g}"
                key = ("end_to_end", structure, temperature, target)
                if key in done:
                    continue
                if early_stop:
                    row = {"experiment": "end_to_end", "structure": structure,
                           "temperature_uK": temperature, "target": target,
                           "status": "skipped_early_stop", "error": ""}
                    _append_csv_row(out_csv, row, fieldnames=SCAN_COLUMNS)
                    rows.append(row)
                    continue
                spec = end_to_end_specs(duration, bcfg["trapezoid_accel_fraction"])[structure]
                waveform = build_speed_waveform(spec, physics)
                cal = speed_calibers(waveform)
                dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                                      physics["aod_waist_m"], bcfg["base_dt_us"] * 1e-6,
                                      min_steps_per_move=bcfg.get("min_steps_per_move", 2000),
                                      max_step_displacement_over_waist=bcfg.get(
                                          "max_step_displacement_over_waist", 0.01))
                dt = snap_dt(waveform.duration_s, dt_max)
                try:
                    point = _evaluate_point(spec, states, physics, dt, 0.0, scales,
                                            follow_cfg,
                                            bcfg.get("work_energy_shots", 4),
                                            num_procs=num_procs)
                    s, cal, wilson = point["summary"], point["calibers"], point["wilson"]
                    row = {
                        "experiment": "end_to_end", "structure": structure,
                        "temperature_uK": temperature, "target": target,
                        "duration_us": waveform.duration_s * 1e6, "dt_ns": dt * 1e9,
                        "v_avg_mm_per_s": cal["v_avg_m_per_s"] * 1e3,
                        "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3,
                        "v_rms_mm_per_s": cal["v_rms_m_per_s"] * 1e3,
                        "eta_peak": point["eta_peak"],
                        "step_displacement_frac": point["step_displacement_frac"],
                        "move_duration_us": cal["move_duration_s"] * 1e6,
                        "shots": s["shots"], "captured": s["captured"],
                        "capture_rate": s["capture_rate"],
                        "wilson_low": wilson[0], "wilson_high": wilson[1],
                        "follow_loss_rate": s["follow_loss_rate"],
                        "exc_move_end_median": s["excitation_rel_move_end_median"],
                        "exc_move_end_q90": s["excitation_rel_move_end_q90"],
                        "margin_q10": s["margin_quantile"],
                        "failure_follow_loss": s["failure_modes"]["follow_loss"],
                        "failure_ramp_excitation": s["failure_modes"]["ramp_excitation"],
                        "failure_handoff": s["failure_modes"]["handoff_failure"],
                        "work_residual_median_over_depth": s.get("work_residual_median_over_depth", ""),
                        "status": "ok", "error": "",
                    }
                    new_flags[_flag_key("end_to_end", structure, temperature, target)] = \
                        point["captured_flags"]
                    if s["capture_rate"] < stop_below:
                        early_stop = True
                except Exception as exc:
                    row = {"experiment": "end_to_end", "structure": structure,
                           "temperature_uK": temperature, "target": target,
                           "status": "failed",
                           "error": f"{type(exc).__name__}: {exc}"}
                _append_csv_row(out_csv, row, fieldnames=SCAN_COLUMNS)
                rows.append(row)
    if new_flags:
        _save_flags(out_dir, new_flags)
    return rows


def load_scan_results(path: Path) -> list[dict]:
    """读取扫描 CSV（含断点续跑前的历史行）。"""
    import csv
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def critical_speed_table(rows: list[dict], flags: dict, bcfg: dict,
                         num_bootstrap: int = 2000, seed: int = 23000) -> list[dict]:
    """对每条 end_to_end 曲线（结构×温度）给出两种方法的临界速度与覆盖判定。"""
    curves = {}
    for row in rows:
        if row.get("experiment") != "end_to_end" or row.get("status") != "ok":
            continue
        key = (row["structure"], float(row["temperature_uK"]))
        curves.setdefault(key, []).append(row)
    levels = [float(x) for x in bcfg.get("capture_target_levels", [0.95, 0.50])]
    table = []
    for curve_index, ((structure, temperature), points) in enumerate(sorted(curves.items())):
        points.sort(key=lambda r: float(r["v_avg_mm_per_s"]))
        speeds = np.array([float(r["v_avg_mm_per_s"]) for r in points]) * 1e-3
        rates = [float(r["capture_rate"]) for r in points]
        flag_arrays = [flags.get(_flag_key("end_to_end", structure, temperature, r["target"]))
                       for r in points]
        entry = {"structure": structure, "temperature_uK": temperature,
                 "n_points": len(points),
                 "v_avg_min_mm_per_s": float(speeds.min() * 1e3),
                 "v_avg_max_mm_per_s": float(speeds.max() * 1e3),
                 "min_capture_rate": float(min(rates)),
                 "max_capture_rate": float(max(rates)),
                 "boundary_covered": boundary_covered(rates, levels[0])}
        if all(f is not None for f in flag_arrays):
            for level_index, level in enumerate(levels):
                boot = bootstrap_interpolation_speeds(
                    speeds, flag_arrays, level, num_bootstrap=num_bootstrap,
                    seed=seed + curve_index * 17 + level_index)
                logi = logistic_fit_speeds(
                    speeds, flag_arrays, level,
                    num_bootstrap=max(100, num_bootstrap // 4),
                    seed=seed + 500 + curve_index * 17 + level_index)
                entry[f"v{int(level * 100)}_bootstrap"] = boot
                entry[f"v{int(level * 100)}_logistic"] = logi
        else:
            entry["flags_missing"] = True
        table.append(entry)
    return table


def fit_excitation_scaling(rows: list[dict], num_bootstrap: int = 1000,
                           seed: int = 0) -> dict:
    """实验 1 标度律：log-log 拟合 ε_med ∝ v^α（bootstrap CI）。"""
    points = [r for r in rows if r.get("experiment") == "pure_move"
              and r.get("status") == "ok" and str(r.get("exc_move_end_median", "")) != ""]
    if len(points) < 3:
        return {"fitted": False, "reason": "points_lt_3"}
    v = np.array([float(r["v_peak_mm_per_s"]) for r in points])
    eps = np.array([float(r["exc_move_end_median"]) for r in points])
    mask = eps > 0
    if mask.sum() < 3:
        return {"fitted": False, "reason": "nonpositive_median"}
    log_v, log_e = np.log10(v[mask]), np.log10(eps[mask])
    slope, intercept = np.polyfit(log_v, log_e, 1)
    residuals = log_e - (slope * log_v + intercept)
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((log_e - log_e.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rng = np.random.default_rng(seed)
    slopes = []
    for _ in range(int(num_bootstrap)):
        idx = rng.integers(0, log_v.size, size=log_v.size)
        if np.unique(idx).size < 2:
            continue
        b, a = np.polyfit(log_v[idx], log_e[idx], 1)
        slopes.append(b)
    return {
        "fitted": True, "alpha": float(slope), "intercept": float(intercept),
        "r_squared": float(r2),
        "alpha_ci_low": float(np.percentile(slopes, 2.5)) if slopes else None,
        "alpha_ci_high": float(np.percentile(slopes, 97.5)) if slopes else None,
        "n_points": int(mask.sum()),
        "v_peak_range_mm_per_s": [float(v[mask].min()), float(v[mask].max())],
    }


def estimate_budget(bcfg: dict) -> dict:
    """运行前估算总 shot×步数与最坏单点步数（写入配置与报告）。"""
    n_exp1 = len(bcfg["pure_move_cruise_speeds_mm_per_s"])
    n_exp2 = len(bcfg["end_to_end_durations_us"]) * len(bcfg["scan_temperatures_uK"]) * 3
    shots = int(bcfg["scan_shots_per_point"])
    worst_steps = 0
    for t_us in bcfg["end_to_end_durations_us"]:
        move_s = (FROZEN_MOVE_WINDOW[1] - FROZEN_MOVE_WINDOW[0]) * t_us * 1e-6
        worst_steps = max(worst_steps, int(t_us * 1e-6 / (bcfg["base_dt_us"] * 1e-6)) + 1,
                          int(move_s / (move_s / bcfg.get("min_steps_per_move", 2000))))
    return {
        "exp1_points": n_exp1, "exp2_curve_points_max": n_exp2,
        "shots_per_point": shots,
        "total_shots_upper_bound": (n_exp1 + n_exp2) * shots,
        "worst_steps_per_shot_approx": worst_steps,
    }
