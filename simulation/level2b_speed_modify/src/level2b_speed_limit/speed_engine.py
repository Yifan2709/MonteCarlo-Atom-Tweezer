"""Level 2B 评估引擎：复用 level0/1/2 的势、积分器、采样器与统计。

新增的诊断是速度研究特有的：跟随丢失（|x−c|/w 与相对激发）、失效模式
三分类（跟丢/降深段激发/交接失败）与两段式积分（转移段按自适应 dt、
保持段按基准 dt，均为同一种 velocity-Verlet）。
"""
from __future__ import annotations

import numpy as np

from level0_static_trap.potentials import (gaussian_force, gaussian_potential,
                                           trap_frequency_analytic)
from level0_static_trap.integrators import velocity_verlet
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level2_joint_transfer.level2_simulation import (is_captured,
                                                     physics_from_config,
                                                     sample_initial_states)
from level2_joint_transfer.noise_heating import MeanHeating, \
    velocity_verlet_time_dependent_heating
from level2_joint_transfer.work_energy import (hold_segment_energy_error,
                                               work_energy_along_trajectory)
from .speed_waveforms import build_speed_waveform

# 判俘获在 t=T；AOD 关闭（深度 0）时 E_f<0 判 captured。纯移动实验中
# capture 字段无物理意义（AOD 深度恒为 D0），只用于完整转移剖面。


def derived_scales(physics: dict) -> dict:
    """AOD/SLM 解析阱频、绝热速度尺度与周期。"""
    mass = physics["mass_kg"]
    omega_aod, nu_aod = trap_frequency_analytic(physics["aod_initial_depth_j"],
                                                physics["aod_waist_m"], mass)
    omega_slm, nu_slm = trap_frequency_analytic(physics["slm_depth_j"],
                                                physics["slm_waist_m"], mass)
    v_ad = physics["aod_waist_m"] * omega_aod
    return {
        "omega_aod_rad_per_s": omega_aod, "nu_aod_hz": nu_aod,
        "omega_slm_rad_per_s": omega_slm, "nu_slm_hz": nu_slm,
        "period_aod_us": 1.0e6 / nu_aod, "period_slm_us": 1.0e6 / nu_slm,
        "adiabatic_speed_m_per_s": v_ad, "adiabatic_speed_mm_per_s": v_ad * 1e3,
    }


def sample_states(cfg2, seed: int, shots: int, temperature_uK=None) -> dict:
    """复用 Level 2 的 Level 1 兼容采样器生成固定初态数组。"""
    return sample_initial_states(cfg2, seed, shots, temperature_uK=temperature_uK)


def _relative_aod_excitation(x, v, t, waveform, physics) -> np.ndarray:
    """ε_AOD(t)/D0：相对 AOD 阱底的激发除以初始阱深（时间数组输入）。"""
    mass = physics["mass_kg"]
    depth = np.asarray(waveform.depth(t), dtype=float)
    center = np.asarray(waveform.center(t), dtype=float)
    u_aod = -depth * np.exp(-2.0 * (x - center) ** 2 / physics["aod_waist_m"] ** 2)
    energy = 0.5 * mass * v**2 + u_aod
    return (energy + physics["aod_initial_depth_j"]) / physics["aod_initial_depth_j"]


def make_e_osc_fn(waveform, physics: dict):
    """参量噪声通道使用的振荡能量函数（约定与 Level 2 完全一致）。

    E_osc = E_tot 相对双阱总势阱底的近似：阱底取两个阱心处总势的较深者
    （分离阱 ≈ 真实阱底，重合阱严格给出 −(D_SLM+D_AOD)），避免用单一
    阱深作底导致束缚原子 ε_osc 被钳到 0。
    """
    slm_depth, slm_waist, slm_center = (physics["slm_depth_j"], physics["slm_waist_m"],
                                        physics["slm_center_m"])
    aod_waist = physics["aod_waist_m"]
    mass = physics["mass_kg"]

    def e_osc_fn(x, v, t):
        depth = float(waveform.depth(t))
        center = float(waveform.center(t))
        u_at_aod_center = (-depth + float(gaussian_potential(center, slm_depth,
                                                             slm_waist, slm_center)))
        u_at_slm_center = (-slm_depth + float(gaussian_potential(slm_center, depth,
                                                                 aod_waist, center)))
        floor = min(u_at_aod_center, u_at_slm_center)
        energy_total = (0.5 * mass * v * v
                        + float(gaussian_potential(x, slm_depth, slm_waist, slm_center))
                        + float(gaussian_potential(x, depth, aod_waist, center)))
        return max(0.0, energy_total - floor)

    return e_osc_fn


def evaluate_profile(spec: dict, states: dict, physics: dict, dt_s: float,
                     hold_s: float = 0.0, hold_dt_s: float | None = None,
                     follow_loss_sep_w: float = 2.0,
                     follow_loss_excitation: float = 0.5,
                     work_energy_shots: int = 0,
                     keep_trajectory_shots: int = 0,
                     heating: MeanHeating | None = None):
    """在固定初态数组上评估一个速度剖面；返回逐 shot 记录与汇总。

    dt_s 为转移段步长（须整除 T）；hold_s>0 时保持段用 hold_dt_s
    （默认 50 ns）继续积分，两段均为同一种 velocity-Verlet。
    跟随丢失阈值 follow_loss_* 一旦由调用方（配置）给定即固定使用。
    heating 给定时按平均强度确定性注入三通道噪声（无随机性），
    转移与保持段全程生效；注入能量计入功-能账本
    （恒等式 ΔE = W_ext + W_noise + R_W）。
    """
    waveform = build_speed_waveform(spec, physics)
    duration_s = float(waveform.duration_s)
    total_s = duration_s + max(0.0, hold_s)
    if not np.isclose(duration_s / dt_s, round(duration_s / dt_s), rtol=0, atol=1e-9 * max(1.0, duration_s / dt_s)):
        raise ValueError(f"T/dt 必须为整数：{duration_s / dt_s:.12g}")
    steps = int(round(duration_s / dt_s))
    times = np.arange(steps + 1, dtype=float) * dt_s
    end_index = steps
    mass = physics["mass_kg"]
    slm = (physics["slm_depth_j"], physics["slm_waist_m"], physics["slm_center_m"])
    waist = physics["aod_waist_m"]

    def total_force(x, t):
        center = float(waveform.center(t))
        depth = float(waveform.depth(t))
        return (gaussian_force(x, slm[0], slm[1], slm[2])
                + gaussian_force(x, depth, waist, center))

    hold_dt = hold_dt_s if hold_dt_s is not None else 5.0e-8
    hold_steps = int(round(hold_s / hold_dt)) if hold_s > 0 else 0
    hold_times = (duration_s + np.arange(1, hold_steps + 1, dtype=float) * hold_dt
                  if hold_steps > 0 else np.empty(0))
    e_osc_fn = make_e_osc_fn(waveform, physics) if heating is not None else None

    x0 = states["x0_m"]
    v0 = states["v0_m_per_s"]
    records = []
    trajectories = {}
    we_ids = set(range(0, x0.size, max(1, x0.size // max(1, work_energy_shots)))) \
        if work_energy_shots > 0 else set()
    keep_ids = set(range(0, x0.size, max(1, x0.size // max(1, keep_trajectory_shots)))) \
        if keep_trajectory_shots > 0 else set()

    for shot in range(x0.size):
        if heating is None:
            time, x, v, _a = velocity_verlet_time_dependent(
                total_force, mass, float(x0[shot]), float(v0[shot]), times)
            noise_work = None
        else:
            time, x, v, _a, noise_work = velocity_verlet_time_dependent_heating(
                total_force, mass, float(x0[shot]), float(v0[shot]), times,
                heating, e_osc_fn)
        # 跟随诊断：转移段内原子-阱心分离与相对初始阱深的激发
        centers = np.asarray(waveform.center(time), dtype=float)
        depths = np.asarray(waveform.depth(time), dtype=float)
        sep_over_w = np.abs(x - centers) / waist
        exc_rel = _relative_aod_excitation(x, v, time, waveform, physics)
        move_t0 = getattr(waveform, "move_t0", 0.0)
        if hasattr(waveform, "move_duration_s"):
            move_t1 = move_t0 + float(waveform.move_duration_s)
        else:
            move_t1 = duration_s
        move_mask = (time >= move_t0 - dt_s) & (time <= move_t1 + dt_s)
        max_sep_move = float(np.max(sep_over_w[move_mask]))
        max_exc_move_end = float(exc_rel[end_index])
        # 失效三分类（深度相位区分）：AOD 深度 ≥50% D0 的承载段内逃出
        # ±2w → 跟丢；深度 <50% D0 的降深段内逃出 → 降深段激发；
        # 全程未逃出但 E_f≥0 → 交接失败。纯移动（深度恒定）只有跟丢。
        carry_mask = depths >= 0.5 * physics["aod_initial_depth_j"]
        escaped_carry = bool(np.any(sep_over_w[carry_mask] > follow_loss_sep_w))
        escaped_ramp = bool(np.any(sep_over_w[~carry_mask] > follow_loss_sep_w))
        follow_lost = bool(max_sep_move > follow_loss_sep_w
                           or max_exc_move_end > follow_loss_excitation)
        # 保持段（可选）：静态极限下继续同一种积分器（开噪时保持段继续注入）
        if hold_steps > 0:
            if heating is None:
                t2, x2, v2, _ = velocity_verlet_time_dependent(
                    total_force, mass, float(x[-1]), float(v[-1]), hold_times)
                hold_noise = None
            else:
                t2, x2, v2, _, hold_noise = velocity_verlet_time_dependent_heating(
                    total_force, mass, float(x[-1]), float(v[-1]), hold_times,
                    heating, e_osc_fn)
            time_full = np.concatenate([time, t2])
            x_full = np.concatenate([x, x2])
            v_full = np.concatenate([v, v2])
            sep_hold = float(np.max(np.abs(x_full - np.asarray(
                waveform.center(time_full), dtype=float)) / waist))
            hold = hold_segment_energy_error(time_full, x_full, v_full, mass,
                                             slm[0], slm[1], slm[2], duration_s)
            hold_pp = hold["peak_to_peak_energy_error_over_slm_depth"]
            follow_lost = follow_lost or sep_hold > follow_loss_sep_w
            if noise_work is not None and hold_noise is not None:
                # 注入能量并到保持段末（与 Level 2 的"转移+保持全程生效"一致）
                for key in noise_work:
                    noise_work[key] = np.concatenate(
                        [noise_work[key], hold_noise[key]])
        else:
            hold_pp = None
        x_f, v_f = float(x[end_index]), float(v[end_index])
        final_energy = 0.5 * mass * v_f**2 + float(gaussian_potential(x_f, *slm))
        captured = is_captured(final_energy, physics.get("capture_threshold", 0.0))
        if captured:
            failure_mode = None
        elif escaped_carry:
            failure_mode = "follow_loss"
        elif escaped_ramp:
            failure_mode = "ramp_excitation"
        else:
            failure_mode = "handoff_failure"
        record = {
            "shot_id": shot,
            "x0_um": float(x0[shot]) * 1e6, "v0_m_s": float(v0[shot]),
            "final_x_um": x_f * 1e6, "final_v_m_s": v_f,
            "final_slm_energy_uK": final_energy / 1.380649e-29 / 1e6,
            "captured": captured,
            "capture_margin": -final_energy / slm[0],
            "final_excitation_over_depth": final_energy / slm[0] + 1.0,
            "max_sep_over_waist_move": max_sep_move,
            "aod_excitation_rel_at_move_end": max_exc_move_end,
            "follow_lost": follow_lost,
            "failure_mode": failure_mode,
            "hold_peak_to_peak_energy_error_over_depth": hold_pp,
            "end_time_index": end_index,
        }
        if noise_work is not None:
            # 沿用 Level 2 既有列名与单位约定：*_uK 字段数值单位实为开尔文
            # （J/kB），与 level2 validation_shots.csv 同约定以便跨列运算。
            record.update({
                "noise_recoil_energy_uK": float(noise_work["recoil_J"][end_index])
                / 1.380649e-29 / 1e6,
                "noise_parametric_energy_uK": float(noise_work["parametric_J"][end_index])
                / 1.380649e-29 / 1e6,
                "noise_pointing_energy_uK": float(noise_work["pointing_J"][end_index])
                / 1.380649e-29 / 1e6,
                "noise_total_energy_uK": float(noise_work["total_J"][end_index])
                / 1.380649e-29 / 1e6,
            })
        if shot in we_ids:
            we = work_energy_along_trajectory(time[:end_index + 1], x[:end_index + 1],
                                              v[:end_index + 1], waveform, mass,
                                              slm[0], slm[1], slm[2], waist,
                                              noise_work_j=None if noise_work is None
                                              else noise_work["total_J"][:end_index + 1])
            record["max_abs_work_energy_residual_over_depth"] = \
                we["max_abs_residual_J"] / slm[0]
            record["final_total_work_uK"] = we["final_work_total_J"] / 1.380649e-29 / 1e6
            record["final_move_work_uK"] = we["final_work_move_J"] / 1.380649e-29 / 1e6
            record["final_depth_work_uK"] = we["final_work_depth_J"] / 1.380649e-29 / 1e6
            if noise_work is not None:
                record["final_noise_work_uK"] = record["noise_total_energy_uK"]
        records.append(record)
        if shot in keep_ids:
            trajectories[shot] = {
                "time_s": time, "position_m": x, "velocity_m_per_s": v,
                "aod_center_m": centers, "aod_depth_j": np.asarray(
                    waveform.depth(time), dtype=float),
                "sep_over_waist": sep_over_w, "exc_rel_aod": exc_rel,
                "end_index": end_index,
            }

    captured_flags = [r["captured"] for r in records]
    margins = np.array([r["capture_margin"] for r in records])
    exc_move_end = np.array([r["aod_excitation_rel_at_move_end"] for r in records])
    summary = {
        "shots": len(records),
        "captured": int(sum(captured_flags)),
        "capture_rate": float(np.mean(captured_flags)),
        "margin_quantile": float(np.quantile(margins, 0.10)),
        "follow_loss_rate": float(np.mean([r["follow_lost"] for r in records])),
        "excitation_rel_move_end_median": float(np.median(exc_move_end)),
        "excitation_rel_move_end_q90": float(np.quantile(exc_move_end, 0.90)),
        "failure_modes": _failure_mode_counts(records),
    }
    residual = [r["max_abs_work_energy_residual_over_depth"] for r in records
                if "max_abs_work_energy_residual_over_depth" in r]
    if residual:
        summary["work_residual_median_over_depth"] = float(np.median(residual))
        summary["work_residual_max_over_depth"] = float(np.max(residual))
    noise_total = [r["noise_total_energy_uK"] for r in records
                   if "noise_total_energy_uK" in r]
    if noise_total:
        summary["mean_noise_total_energy_uK"] = float(np.mean(noise_total))
        summary["mean_noise_recoil_energy_uK"] = float(np.mean(
            [r["noise_recoil_energy_uK"] for r in records]))
        summary["mean_noise_parametric_energy_uK"] = float(np.mean(
            [r["noise_parametric_energy_uK"] for r in records]))
        summary["mean_noise_pointing_energy_uK"] = float(np.mean(
            [r["noise_pointing_energy_uK"] for r in records]))
    return {"spec": spec, "records": records, "summary": summary,
            "captured_flags": captured_flags, "trajectories": trajectories}


def _failure_mode_counts(records) -> dict:
    """失效模式三分类计数（captured 的 shot 不计入）。"""
    counts = {"follow_loss": 0, "ramp_excitation": 0, "handoff_failure": 0}
    for r in records:
        if r["failure_mode"] in counts:
            counts[r["failure_mode"]] += 1
    return counts


def _evaluate_chunk(args):
    """子进程入口：在初态数组的一个分片上评估剖面。"""
    (spec, states, physics, dt_s, hold_s, hold_dt_s, sep_w, exc, we_shots,
     heating, offset, size) = args
    subset = {key: value[offset:offset + size] if isinstance(value, np.ndarray) else value
              for key, value in states.items()}
    result = evaluate_profile(spec, subset, physics, dt_s, hold_s=hold_s,
                              hold_dt_s=hold_dt_s, follow_loss_sep_w=sep_w,
                              follow_loss_excitation=exc,
                              work_energy_shots=we_shots, heating=heating)
    for record in result["records"]:
        record["shot_id"] += offset
    return result


def evaluate_profile_parallel(spec: dict, states: dict, physics: dict, dt_s: float,
                              hold_s: float = 0.0, hold_dt_s: float | None = None,
                              follow_loss_sep_w: float = 2.0,
                              follow_loss_excitation: float = 0.5,
                              work_energy_shots: int = 0, num_procs: int | None = None,
                              heating: MeanHeating | None = None):
    """按 shot 分片并行评估一个剖面；初态数组本身不变（保持 CRN）。

    heating 需可 pickle（build_mean_heating / scaled_heating 构造的
    常数模型满足）；自定义 lambda 仅限 evaluate_profile 单进程路径。
    """
    import os
    from multiprocessing import Pool
    procs = num_procs or max(1, (os.cpu_count() or 2) - 1)
    total = int(states["shots"])
    chunk = max(1, (total + procs - 1) // procs)
    jobs = [(spec, states, physics, dt_s, hold_s, hold_dt_s, follow_loss_sep_w,
             follow_loss_excitation, work_energy_shots, heating, offset,
             min(chunk, total - offset)) for offset in range(0, total, chunk)]
    with Pool(len(jobs)) as pool:
        parts = pool.map(_evaluate_chunk, jobs)
    records = []
    for part in parts:
        records.extend(part["records"])
    records.sort(key=lambda r: r["shot_id"])
    captured_flags = [r["captured"] for r in records]
    margins = np.array([r["capture_margin"] for r in records])
    exc_move_end = np.array([r["aod_excitation_rel_at_move_end"] for r in records])
    summary = {
        "shots": len(records),
        "captured": int(sum(captured_flags)),
        "capture_rate": float(np.mean(captured_flags)),
        "margin_quantile": float(np.quantile(margins, 0.10)),
        "follow_loss_rate": float(np.mean([r["follow_lost"] for r in records])),
        "excitation_rel_move_end_median": float(np.median(exc_move_end)),
        "excitation_rel_move_end_q90": float(np.quantile(exc_move_end, 0.90)),
        "failure_modes": _failure_mode_counts(records),
    }
    residual = [r["max_abs_work_energy_residual_over_depth"] for r in records
                if "max_abs_work_energy_residual_over_depth" in r]
    if residual:
        summary["work_residual_median_over_depth"] = float(np.median(residual))
        summary["work_residual_max_over_depth"] = float(np.max(residual))
    noise_total = [r["noise_total_energy_uK"] for r in records
                   if "noise_total_energy_uK" in r]
    if noise_total:
        summary["mean_noise_total_energy_uK"] = float(np.mean(noise_total))
        summary["mean_noise_recoil_energy_uK"] = float(np.mean(
            [r["noise_recoil_energy_uK"] for r in records]))
        summary["mean_noise_parametric_energy_uK"] = float(np.mean(
            [r["noise_parametric_energy_uK"] for r in records]))
        summary["mean_noise_pointing_energy_uK"] = float(np.mean(
            [r["noise_pointing_energy_uK"] for r in records]))
    return {"spec": spec, "records": records, "summary": summary,
            "captured_flags": captured_flags}


def static_limit_check(physics: dict, dt_s: float, duration_s: float) -> dict:
    """静态极限：AOD/SLM 均静止时本引擎积分路径与 Level 0 积分器逐点一致。"""
    slm = (physics["slm_depth_j"], physics["slm_waist_m"], physics["slm_center_m"])
    steps = int(round(duration_s / dt_s))
    times = np.arange(steps + 1, dtype=float) * dt_s
    x_init = float(physics["aod_initial_center_m"])

    def force_static(x):
        return (gaussian_force(x, slm[0], slm[1], slm[2])
                + gaussian_force(x, physics["aod_initial_depth_j"],
                                 physics["aod_waist_m"], physics["aod_initial_center_m"]))

    t0, x0, v0, _ = velocity_verlet(force_static, physics["mass_kg"], x_init, 0.0, times)

    def force_td(x, t):
        return force_static(x)

    t1, x1, v1, _ = velocity_verlet_time_dependent(force_td, physics["mass_kg"],
                                                   x_init, 0.0, times)
    return {
        "max_position_diff": float(np.max(np.abs(x0 - x1))),
        "max_velocity_diff": float(np.max(np.abs(v0 - v1))),
        "passed": bool(np.allclose(x0, x1, rtol=0, atol=1e-15)
                       and np.allclose(v0, v1, rtol=0, atol=1e-15)),
    }
