"""Level 2B 实验 4：平均强度噪声加热 × 移动速度验证。

复用 level2_joint_transfer.noise_heating 的三通道确定性平均加热模型
（光子反冲 / 强度噪声-参量 / 指向噪声），回答四个子问题：

4a 边界移动：开噪（×1、×100）后端到端 v95 曲线是否移动？
4b 开关配对：边界附近若干峰值速度、两个温度下，开/关噪声的俘获率
   配对差与裕量退化（同初态 CRN，噪声无随机性故为确定性配对）。
4c 幅度敏感性：噪声整体放大 ×1→×1000 时俘获率/裕量的退化起点。
4d 通道分解：×100 下单通道（仅反冲/仅参量/仅指向）各自造成多少退化。

约定与 Level 2 完全一致：噪声为系综平均功率的确定性注入（无随机性），
幅度参数均为 assumed_sensitivity_only；scan 阶段保持无噪声（冻结基线），
本阶段只做敏感性评估、不回调任何波形参数。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from level1_transfer_1d.analysis import wilson_interval
from level0_static_trap.constants import microkelvin_to_joule
from level2_joint_transfer.noise_heating import (MeanHeating, build_mean_heating,
                                                 parametric_rate_from_rin,
                                                 pointing_power_from_psd,
                                                 recoil_energy_j,
                                                 recoil_scattering_rate_cs,
                                                 scaled_heating)
from level2_joint_transfer.statistics import paired_bootstrap_difference, paired_counts
from level0_static_trap.potentials import trap_frequency_analytic
from .critical_speed import (boundary_covered, bootstrap_interpolation_speeds,
                             logistic_fit_speeds)
from .profile_compare import _profile_spec
from .speed_engine import evaluate_profile_parallel, sample_states
from .speed_scan import (FROZEN_MOVE_WINDOW, _append_csv_row, _save_flags,
                         load_flags)
from .speed_waveforms import build_speed_waveform, max_step_for, snap_dt, speed_calibers

NOISE_COLUMNS = [
    "experiment", "sub_experiment", "noise_mode", "structure", "temperature_uK",
    "condition", "duration_us", "dt_ns", "v_peak_mm_per_s", "eta_peak",
    "step_displacement_frac", "shots", "captured", "capture_rate",
    "wilson_low", "wilson_high", "margin_q10", "margin_mean",
    "exc_move_end_median", "follow_loss_rate",
    "noise_recoil_injected_mean", "noise_parametric_injected_mean",
    "noise_pointing_injected_mean", "noise_total_injected_mean",
    "work_residual_median_over_depth",
    "pair_difference", "pair_ci_low", "pair_ci_high",
    "both_captured", "only_on", "only_off", "both_failed",
    "status", "error",
]


def build_base_heating(bcfg: dict, physics: dict) -> MeanHeating:
    """按本级配置的默认幅度构造三通道基础模型（与 Level 2 默认相同量级）。

    反冲：2·E_r·Γ_sc（Cs D2 两能级，标称 AOD 深度）；
    参量：P = Γ_par·E_ref（线性化常数功率），其中
    Γ_par = (π²/4)·ν₀²·S_RIN(2ν₀)（ν₀ 取 SLM 解析阱频）、
    E_ref = k_B·T_ref（T_ref 取噪声实验运行温度，默认 20 μK）；
    指向：P = m·ω₀⁴·S_x(ν₀)/8。全部参数 assumed_sensitivity_only。
    """
    noise = bcfg["noise_speed_validation"]
    mass = physics["mass_kg"]
    wavelength = float(noise["trap_wavelength_nm"]) * 1e-9
    scattering = recoil_scattering_rate_cs(wavelength, physics["aod_initial_depth_j"])
    recoil_power = 2.0 * recoil_energy_j(wavelength, mass) * scattering
    nu0 = trap_frequency_analytic(physics["slm_depth_j"], physics["slm_waist_m"],
                                  mass)[1]
    parametric_rate = parametric_rate_from_rin(nu0, float(noise["parametric_rin_psd_per_hz"]))
    reference_temperature_uK = float(noise.get("parametric_reference_temperature_uK", 20.0))
    parametric_power = parametric_rate * microkelvin_to_joule(reference_temperature_uK)
    pointing_power = pointing_power_from_psd(mass, nu0,
                                             float(noise["pointing_position_psd_m2_per_hz"]))
    heating = build_mean_heating(recoil_power_w=recoil_power,
                                 parametric_power_w=parametric_power,
                                 pointing_power_w=pointing_power)
    heating.labels["parametric"].update({
        "rate_per_s": float(parametric_rate),
        "reference_temperature_uK": reference_temperature_uK,
    })
    heating.labels["derived_from"] = {
        "trap_wavelength_nm": float(noise["trap_wavelength_nm"]),
        "recoil_scattering_rate_s": float(scattering),
        "recoil_energy_kB_nK": float(recoil_energy_j(wavelength, mass)
                                     / 1.380649e-29 * 1e9),
        "reference_nu0_hz": float(nu0),
        "parametric_rin_psd_per_hz": float(noise["parametric_rin_psd_per_hz"]),
        "parametric_rate_per_s": float(parametric_rate),
        "parametric_reference_temperature_uK": reference_temperature_uK,
        "pointing_position_psd_m2_per_hz": float(noise["pointing_position_psd_m2_per_hz"]),
        "recoil_power_uK_per_s": float(recoil_power / 1.380649e-29 / 1e6),
        "parametric_power_uK_per_s": float(parametric_power / 1.380649e-29 / 1e6),
        "pointing_power_uK_per_s": float(pointing_power / 1.380649e-29 / 1e6),
    }
    return heating


def single_channel_heating(base: MeanHeating, channel: str,
                           scale: float) -> MeanHeating:
    """从基础模型提取单通道并整体放大 scale 倍（其余通道关闭）。"""
    scaled = scaled_heating(base, scale)
    kwargs = {channel: getattr(scaled, channel)}
    heating = MeanHeating(recoil=kwargs.get("recoil"),
                          parametric=kwargs.get("parametric"),
                          pointing=kwargs.get("pointing"),
                          labels=dict(scaled.labels))
    heating.labels["single_channel"] = channel
    heating.labels["scale_factor"] = float(scale)
    return heating


def _dt_for_spec(spec: dict, physics: dict, base_dt_s: float, integration: dict):
    """按剖面实际峰值速度与移动时长取自适应步长。"""
    waveform = build_speed_waveform(spec, physics)
    cal = speed_calibers(waveform)
    dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                          physics["aod_waist_m"], base_dt_s,
                          min_steps_per_move=integration.get("min_steps_per_move", 2000),
                          max_step_displacement_over_waist=integration.get(
                              "max_step_displacement_over_waist", 0.01))
    return snap_dt(waveform.duration_s, dt_max), waveform, cal


def _row_from_result(sub, mode, structure, temperature, condition, spec, waveform,
                     cal, dt, physics, result, scales=None):
    """把一次评估结果整理成 CSV 行。"""
    s = result["summary"]
    captured, shots = s["captured"], s["shots"]
    low, high = wilson_interval(captured, shots)
    margins = np.array([r["capture_margin"] for r in result["records"]])
    eta = cal["v_peak_m_per_s"] / scales["adiabatic_speed_m_per_s"] if scales else ""
    return {
        "experiment": "noise", "sub_experiment": sub, "noise_mode": mode,
        "structure": structure, "temperature_uK": temperature,
        "condition": condition,
        "duration_us": waveform.duration_s * 1e6, "dt_ns": dt * 1e9,
        "v_peak_mm_per_s": cal["v_peak_m_per_s"] * 1e3, "eta_peak": eta,
        "step_displacement_frac": dt * cal["v_peak_m_per_s"] / physics["aod_waist_m"],
        "shots": shots, "captured": captured, "capture_rate": s["capture_rate"],
        "wilson_low": float(low), "wilson_high": float(high),
        "margin_q10": s["margin_quantile"], "margin_mean": float(np.mean(margins)),
        "exc_move_end_median": s["excitation_rel_move_end_median"],
        "follow_loss_rate": s["follow_loss_rate"],
        "noise_recoil_injected_mean": s.get("mean_noise_recoil_energy_uK", ""),
        "noise_parametric_injected_mean": s.get("mean_noise_parametric_energy_uK", ""),
        "noise_pointing_injected_mean": s.get("mean_noise_pointing_energy_uK", ""),
        "noise_total_injected_mean": s.get("mean_noise_total_energy_uK", ""),
        "work_residual_median_over_depth": s.get("work_residual_median_over_depth", ""),
        "status": "ok", "error": "",
    }


def _spec_for(structure: str, duration_s: float | None, v_peak_m_per_s: float | None,
              accel: float, physics: dict) -> dict:
    """按结构与时长（4a）或峰值速度（4b/4c/4d）构造 spec。"""
    from .speed_scan import end_to_end_specs, FROZEN_RAMP_POWER, FROZEN_RAMP_WINDOW
    if duration_s is not None:
        return end_to_end_specs(duration_s, accel)[structure]
    # 按峰值速度参数化：same_peak_speed 口径下各剖面的 v_peak 都等于目标
    profile = {"trapezoid_frozen_ramp": "trapezoid",
               "frozen_shape_scaled": "frozen_smootherstep"}.get(structure, "trapezoid")
    return _profile_spec(profile, v_peak_m_per_s, "same_peak_speed",
                         physics["aod_initial_center_m"], accel)


def run_noise_stage(cfg2, bcfg: dict, physics: dict, scales: dict,
                    out_csv: Path, out_dir: Path, num_procs=None) -> dict:
    """执行实验 4 的四个子实验并返回汇总（含 v95 表与头部结论）。"""
    noise_cfg = bcfg["noise_speed_validation"]
    integration = bcfg["integration"]
    follow_cfg = bcfg["follow_loss"]
    base_dt = float(integration["base_dt_us"]) * 1e-6
    accel = float(bcfg["speed_scan"]["trapezoid_accel_fraction"])
    shots = int(noise_cfg.get("shots_per_point", 300))
    base = build_base_heating(bcfg, physics)
    seed = int(bcfg["initial_ensemble"]["scan_seed"]) + 321
    new_flags = {}

    def evaluate(spec, states, heating):
        dt, waveform, cal = _dt_for_spec(spec, physics, base_dt, integration)
        return evaluate_profile_parallel(
            spec, states, physics, dt, hold_s=0.0,
            follow_loss_sep_w=follow_cfg["max_separation_over_waist"],
            follow_loss_excitation=follow_cfg["max_excitation_over_depth"],
            work_energy_shots=integration.get("work_energy_shots", 4),
            num_procs=num_procs, heating=heating), dt, waveform, cal

    # ------------------------------------------------ 4a 边界移动（v95 曲线）
    v95_cfg = noise_cfg["v95_shift"]
    v95_curves = {}
    temperature = float(v95_cfg["temperature_uK"])
    states = sample_states(cfg2, seed, shots, temperature_uK=temperature)
    modes = [("off", None)] + [(f"x{int(f)}", scaled_heating(base, float(f)))
                               for f in v95_cfg["scale_factors"]]
    for mode, heating in modes:
        curve = []
        for t_us in v95_cfg["durations_us"]:
            condition = f"{float(t_us):g}"
            try:
                spec = _spec_for(v95_cfg["structure"], float(t_us) * 1e-6,
                                 None, accel, physics)
                result, dt, waveform, cal = evaluate(spec, states, heating)
                row = _row_from_result("v95_shift", mode, v95_cfg["structure"],
                                       temperature, condition, spec, waveform,
                                       cal, dt, physics, result, scales)
                key = f"noise_v95|{mode}|{temperature}|{condition}"
                new_flags[key] = result["captured_flags"]
                row["_flag_key"] = key
            except Exception as exc:
                row = {"experiment": "noise", "sub_experiment": "v95_shift",
                       "noise_mode": mode, "structure": v95_cfg["structure"],
                       "temperature_uK": temperature, "condition": condition,
                       "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
            _append_csv_row(out_csv, row, fieldnames=NOISE_COLUMNS)
            if row.get("status") == "ok":
                curve.append(row)
        v95_curves[mode] = curve

    # v95 估计（用本轮旗标做 bootstrap + logistic 双方法）
    move_frac = _move_fraction(v95_cfg["structure"])
    v95_table = []
    for mode, curve in v95_curves.items():
        if not curve:
            continue
        points = sorted(curve, key=lambda r: float(r["duration_us"]))
        speeds = np.array([2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac)
                           for r in points])
        flags = [new_flags.get(r["_flag_key"]) for r in points]
        entry = {"noise_mode": mode}
        entry.update(_v95_from_curve(points, flags, speeds, seed))
        v95_table.append(entry)

    # ------------------------------------------------ 4b 开/关配对比较
    onoff_cfg = noise_cfg["onoff_paired"]
    onoff_summary = []
    for temperature in [float(t) for t in onoff_cfg["temperatures_uK"]]:
        states = sample_states(cfg2, seed + int(temperature * 1000), shots,
                               temperature_uK=temperature)
        for v_peak in [float(v) * 1e-3 for v in onoff_cfg["speeds_v_peak_mm_per_s"]]:
            condition = f"{v_peak * 1e3:g}"
            try:
                spec = _spec_for(onoff_cfg["structure"], None, v_peak, accel, physics)
                result_off, dt, waveform, cal = evaluate(spec, states, None)
                row_off = _row_from_result("onoff_paired", "off", onoff_cfg["structure"],
                                           temperature, condition, spec, waveform, cal,
                                           dt, physics, result_off, scales)
                _append_csv_row(out_csv, row_off, fieldnames=NOISE_COLUMNS)
                result_on, dt2, waveform2, cal2 = evaluate(spec, states, base)
                row_on = _row_from_result("onoff_paired", "x1", onoff_cfg["structure"],
                                          temperature, condition, spec, waveform2, cal2,
                                          dt2, physics, result_on, scales)
                flags_off = np.array(result_off["captured_flags"])
                flags_on = np.array(result_on["captured_flags"])
                boot = paired_bootstrap_difference(flags_on, flags_off,
                                                   int(bcfg["comparison"]["bootstrap_samples"]),
                                                   seed + int(v_peak * 1e3) + 91)
                counts = paired_counts(flags_on, flags_off)
                row_on.update({
                    "pair_difference": boot["difference"],
                    "pair_ci_low": boot["ci_low"], "pair_ci_high": boot["ci_high"],
                    "both_captured": counts["both_captured"],
                    "only_on": counts["only_a_captured"],
                    "only_off": counts["only_b_captured"],
                    "both_failed": counts["both_failed"],
                })
                _append_csv_row(out_csv, row_on, fieldnames=NOISE_COLUMNS)
                onoff_summary.append({
                    "temperature_uK": temperature,
                    "v_peak_mm_per_s": v_peak * 1e3,
                    "capture_off": result_off["summary"]["capture_rate"],
                    "capture_on": result_on["summary"]["capture_rate"],
                    "difference": boot["difference"],
                    "ci_low": boot["ci_low"], "ci_high": boot["ci_high"],
                    "ci_excludes_zero": boot["ci_excludes_zero"],
                    "margin_shift": row_on["margin_mean"] - row_off["margin_mean"],
                    "noise_injected_total": row_on["noise_total_injected_mean"],
                })
            except Exception as exc:
                row = {"experiment": "noise", "sub_experiment": "onoff_paired",
                       "noise_mode": "x1", "structure": onoff_cfg["structure"],
                       "temperature_uK": temperature, "condition": condition,
                       "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
                _append_csv_row(out_csv, row, fieldnames=NOISE_COLUMNS)

    # ------------------------------------------------ 4c 幅度敏感性
    sens_cfg = noise_cfg["amplitude_sensitivity"]
    sens_states = sample_states(cfg2, seed + 5, shots,
                                temperature_uK=float(sens_cfg["temperature_uK"]))
    sensitivity = []
    v_sens = float(sens_cfg["v_peak_mm_per_s"]) * 1e-3
    sens_spec = _spec_for(sens_cfg["structure"], None, v_sens, accel, physics)
    for scale in [float(s) for s in sens_cfg["scale_factors"]]:
        try:
            result, dt, waveform, cal = evaluate(sens_spec, sens_states,
                                                 None if scale == 0.0 else scaled_heating(base, scale))
            row = _row_from_result("amplitude_sensitivity", f"x{int(scale)}",
                                   sens_cfg["structure"],
                                   float(sens_cfg["temperature_uK"]),
                                   f"{v_sens * 1e3:g}", sens_spec, waveform, cal,
                                   dt, physics, result, scales)
            _append_csv_row(out_csv, row, fieldnames=NOISE_COLUMNS)
            sensitivity.append({"scale": scale, "capture_rate": row["capture_rate"],
                                "margin_mean": row["margin_mean"],
                                "wilson_low": row["wilson_low"],
                                "wilson_high": row["wilson_high"],
                                "noise_injected_total": row["noise_total_injected_mean"]})
        except Exception as exc:
            _append_csv_row(out_csv, {"experiment": "noise",
                                      "sub_experiment": "amplitude_sensitivity",
                                      "noise_mode": f"x{int(scale)}",
                                      "structure": sens_cfg["structure"],
                                      "temperature_uK": float(sens_cfg["temperature_uK"]),
                                      "condition": f"{v_sens * 1e3:g}", "status": "failed",
                                      "error": f"{type(exc).__name__}: {exc}"},
                            fieldnames=NOISE_COLUMNS)

    # ------------------------------------------------ 4d 通道分解
    chan_cfg = noise_cfg["channel_breakdown"]
    chan_states = sample_states(cfg2, seed + 6, shots,
                                temperature_uK=float(chan_cfg["temperature_uK"]))
    v_chan = float(chan_cfg["v_peak_mm_per_s"]) * 1e-3
    chan_spec = _spec_for(chan_cfg["structure"], None, v_chan, accel, physics)
    chan_scale = float(chan_cfg["scale_factor"])
    channels = []
    for mode, heating in [("off", None),
                          ("all_channels", scaled_heating(base, chan_scale)),
                          ("recoil_only", single_channel_heating(base, "recoil", chan_scale)),
                          ("parametric_only", single_channel_heating(base, "parametric", chan_scale)),
                          ("pointing_only", single_channel_heating(base, "pointing", chan_scale))]:
        try:
            result, dt, waveform, cal = evaluate(chan_spec, chan_states, heating)
            row = _row_from_result("channel_breakdown", mode, chan_cfg["structure"],
                                   float(chan_cfg["temperature_uK"]),
                                   f"{v_chan * 1e3:g}", chan_spec, waveform, cal,
                                   dt, physics, result, scales)
            _append_csv_row(out_csv, row, fieldnames=NOISE_COLUMNS)
            channels.append({"mode": mode, "capture_rate": row["capture_rate"],
                             "margin_mean": row["margin_mean"],
                             "recoil_injected": row["noise_recoil_injected_mean"],
                             "parametric_injected": row["noise_parametric_injected_mean"],
                             "pointing_injected": row["noise_pointing_injected_mean"],
                             "total_injected": row["noise_total_injected_mean"]})
        except Exception as exc:
            _append_csv_row(out_csv, {"experiment": "noise",
                                      "sub_experiment": "channel_breakdown",
                                      "noise_mode": mode, "structure": chan_cfg["structure"],
                                      "temperature_uK": float(chan_cfg["temperature_uK"]),
                                      "condition": f"{v_chan * 1e3:g}", "status": "failed",
                                      "error": f"{type(exc).__name__}: {exc}"},
                            fieldnames=NOISE_COLUMNS)

    if new_flags:
        _save_flags(out_dir, new_flags)
    return {"defaults": base.labels, "v95_table": v95_table,
            "onoff": onoff_summary, "sensitivity": sensitivity,
            "channels": channels,
            "config": {"shots_per_point": shots,
                       "v95_structure": v95_cfg["structure"],
                       "v95_temperature_uK": float(v95_cfg["temperature_uK"])}}


def _move_fraction(structure: str) -> float:
    """结构对应的移动窗口占时长比例（用于 duration → v_avg 换算）。"""
    if structure == "sequential_scaled":
        return 0.52
    return FROZEN_MOVE_WINDOW[1] - FROZEN_MOVE_WINDOW[0]


def _v95_from_curve(points: list[dict], flag_arrays, speeds: np.ndarray,
                    seed: int) -> dict:
    """对一条 4a 曲线做 v95 双方法估计（旗标齐全时）。"""
    rates = [float(r["capture_rate"]) for r in points]
    entry = {"n_points": len(points),
             "boundary_covered": boundary_covered(rates, 0.95),
             "min_capture_rate": float(min(rates)),
             "v_avg_max_mm_per_s": float(speeds.max() * 1e3)}
    # 统一转 ndarray：内存路径的旗标是 list，档案路径是 ndarray
    flag_arrays = [None if f is None else np.asarray(f, dtype=bool)
                   for f in flag_arrays]
    if all(f is not None for f in flag_arrays) and len(points) >= 3:
        entry["v95_bootstrap"] = bootstrap_interpolation_speeds(
            speeds, flag_arrays, 0.95, num_bootstrap=2000, seed=seed)
        entry["v95_logistic"] = logistic_fit_speeds(
            speeds, flag_arrays, 0.95, num_bootstrap=200, seed=seed + 1)
    return entry


def compute_noise_v95(out_dir: Path, seed: int = 24000) -> list[dict]:
    """从落盘的 4a 行 + 旗标档案重算 v95（含 CI），供 CLI 汇总与图使用。"""
    import csv
    out_csv = out_dir / "noise_speed_results.csv"
    flags_archive = load_flags(out_dir)
    curves = {}
    structures = {}
    if not out_csv.is_file():
        return []
    for row in csv.DictReader(out_csv.open(encoding="utf-8")):
        if row.get("sub_experiment") != "v95_shift" or row.get("status") != "ok":
            continue
        curves.setdefault(row["noise_mode"], []).append(row)
        structures[row["noise_mode"]] = row["structure"]
    table = []
    for mode, points in curves.items():
        points = sorted(points, key=lambda r: float(r["duration_us"]))
        move_frac = _move_fraction(structures[mode])
        speeds = np.array([2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac)
                           for r in points])
        flags = [flags_archive.get(
            f"noise_v95|{mode}|{points[0]['temperature_uK']}|{r['condition']}")
            for r in points]
        entry = {"noise_mode": mode}
        entry.update(_v95_from_curve(points, flags, speeds, seed))
        table.append(entry)
    return table
