"""Level 2 核心引擎：公共初态采样、波形评估、预检查与独立验证。

所有波形对同一样本池使用完全相同的初态数组（common random numbers）。
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

from level0_static_trap.potentials import gaussian_force, gaussian_potential
from level0_static_trap.integrators import velocity_verlet
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level1_transfer_1d.thermal import sample_thermal_initial_state
from .config import Level2Config, level1_view
from .noise_heating import MeanHeating, velocity_verlet_time_dependent_heating
from .statistics import paired_bootstrap_difference, paired_counts, summarize_capture
from .waveforms import OverlappedWaveform, SequentialWaveform, SplineWaveform, validate_waveform
from .work_energy import hold_segment_energy_error, work_energy_along_trajectory


def physics_from_config(cfg: Level2Config) -> dict:
    """提取子进程可序列化的物理常量字典。"""
    return {
        "mass_kg": cfg.atom.mass_kg,
        "slm_depth_j": cfg.slm.depth_j,
        "slm_waist_m": cfg.slm.waist_m,
        "slm_center_m": cfg.slm.center_m,
        "aod_waist_m": cfg.aod.waist_m,
        "aod_initial_depth_j": cfg.aod.initial_depth_j,
        "aod_initial_center_m": cfg.aod.initial_center_m,
        "near_threshold_fraction": cfg.classification.near_threshold_fraction_of_slm_depth,
        "capture_threshold": cfg.classification.capture_energy_threshold_J,
    }


def build_waveform(spec: dict, physics: dict):
    """由可序列化 spec 构造波形对象。"""
    kind = spec["type"]
    if kind == "sequential":
        return SequentialWaveform(spec["duration_s"], spec["move_fraction"],
                                  physics["aod_initial_center_m"], physics["aod_initial_depth_j"])
    if kind == "overlapped":
        return OverlappedWaveform(spec["duration_s"], spec["move_start_fraction"],
                                  spec["move_end_fraction"], spec["ramp_start_fraction"],
                                  spec["ramp_end_fraction"], spec["ramp_power"],
                                  physics["aod_initial_center_m"], physics["aod_initial_depth_j"])
    if kind == "spline":
        return SplineWaveform(spec["duration_s"], spec["control_s"], spec["center_values_m"],
                              spec["depth_values_j"], physics["aod_initial_center_m"],
                              physics["aod_initial_depth_j"])
    raise ValueError(f"未知波形类型 {kind!r}")


def baseline_specs(cfg: Level2Config) -> dict:
    """三个主波形中两个固定基线的 spec（第三个由优化产生）。"""
    w = cfg.waveforms
    return {
        "sequential_600us": {
            "type": "sequential", "duration_s": w.sequential_baseline_duration_s,
            "move_fraction": w.level1_move_fraction,
        },
        "sequential_400us": {
            "type": "sequential", "duration_s": w.compressed_baseline_duration_s,
            "move_fraction": w.level1_move_fraction,
        },
    }


def sample_initial_states(cfg: Level2Config, seed: int, shots: int,
                          temperature_uK=None, aod_initial_center_m=None):
    """复用 Level 1 简谐提议 + 束缚拒绝采样器，生成固定初态数组。"""
    view = level1_view(cfg, temperature_uK=temperature_uK, aod_initial_center_m=aod_initial_center_m)
    rng = np.random.default_rng(seed)
    states = []
    accepted = 0
    attempts = 0
    while accepted < shots:
        states.append(sample_thermal_initial_state(view, rng))
        accepted += 1
        attempts += 1
        if attempts > shots * 1000:
            raise RuntimeError("初态采样接受率过低")
    x0 = np.array([s.x0_m for s in states])
    v0 = np.array([s.v0_m_s for s in states])
    initial_energy = np.array([s.aod_total_energy_J for s in states])
    initial_excitation = np.array([s.aod_excitation_energy_J for s in states])
    return {
        "seed": int(seed), "shots": int(shots), "attempts": attempts,
        "acceptance_rate": shots / attempts,
        "x0_m": x0, "v0_m_per_s": v0,
        "initial_aod_energy_J": initial_energy,
        "initial_excitation_J": initial_excitation,
        "mean_x0_m": float(x0.mean()), "std_x0_m": float(x0.std(ddof=1)),
        "mean_v0_m_s": float(v0.mean()), "std_v0_m_s": float(v0.std(ddof=1)),
    }


def is_captured(final_energy_j: float, threshold_j: float = 0.0) -> bool:
    """严格俘获判据：E_f < 阈值判俘获，E_f = 0 不俘获。"""
    return bool(final_energy_j < threshold_j)


def _evaluate_chunk(args):
    """子进程入口：在初态数组的一个分片上评估波形。"""
    spec, states, physics, dt_s, hold_s, alignment_offset_m, work_energy, offset, size, heating = args
    subset = {key: value[offset:offset + size] if isinstance(value, np.ndarray) else value
              for key, value in states.items()}
    result = evaluate_waveform_on_states(spec, subset, physics, dt_s, hold_s,
                                         alignment_offset_m=alignment_offset_m,
                                         work_energy=work_energy, heating=heating)
    for record in result["records"]:
        record["shot_id"] += offset
    result["captured_flags"] = [r["captured"] for r in result["records"]]
    return result


def evaluate_waveform_parallel(spec: dict, states: dict, physics: dict, dt_s: float, hold_s: float,
                               alignment_offset_m: float = 0.0, work_energy: bool = False,
                               num_procs: int | None = None, heating: MeanHeating | None = None):
    """按 shot 分片并行评估一个波形；初态数组本身不变（保持 CRN）。

    heating 需可 pickle（build_mean_heating 构造的常数模型满足）；
    lambda 自定义函数仅限 evaluate_waveform_on_states 单进程路径。
    """
    import os
    from multiprocessing import Pool
    procs = num_procs or max(1, (os.cpu_count() or 2) - 1)
    total = int(states["shots"])
    chunk = max(1, (total + procs - 1) // procs)
    jobs = [(spec, states, physics, dt_s, hold_s, alignment_offset_m, work_energy, offset,
             min(chunk, total - offset), heating) for offset in range(0, total, chunk)]
    with Pool(len(jobs)) as pool:
        parts = pool.map(_evaluate_chunk, jobs)
    records = []
    for part in parts:
        records.extend(part["records"])
    captured_flags = [r["captured"] for r in records]
    margins = np.array([r["capture_margin"] for r in records])
    captured_excitation = [r["final_excitation_over_depth"] for r in records if r["captured"]]
    summary = {
        "shots": len(records),
        "captured": int(sum(captured_flags)),
        "capture_rate": float(np.mean(captured_flags)),
        "margin_quantile": float(np.quantile(margins, 0.10)),
        "mean_final_excitation_normalized_captured":
            float(np.mean(captured_excitation)) if captured_excitation else None,
        "max_final_excitation_normalized_captured":
            float(np.max(captured_excitation)) if captured_excitation else None,
        "near_threshold_count": int(sum(r["near_threshold"] for r in records)),
    }
    return {"spec": spec, "records": records, "summary": summary,
            "captured_flags": captured_flags}


def evaluate_waveform_on_states(spec: dict, states: dict, physics: dict, dt_s: float,
                                 hold_s: float, alignment_offset_m: float = 0.0,
                                 work_energy: bool = False, keep_trajectory_ids=(),
                                 heating: MeanHeating | None = None):
    """在固定初态数组上评估一个波形；返回逐 shot 记录与汇总。

    俘获判据在转移结束时刻 t=T（AOD 深度严格为 0）用 SLM 单阱能量判定。
    heating 给定时，按平均强度确定性注入三通道噪声加热（无随机性）；
    注入能量计入功-能账本（恒等式 ΔE = W_ext + W_noise + R_W），
    转移与保持段全程生效，注入能量随时间累积使系统越发不稳定。
    """
    waveform = build_waveform(spec, physics)
    duration_s = waveform.duration_s
    total_s = duration_s + hold_s
    if not math.isclose(total_s / dt_s, round(total_s / dt_s), abs_tol=1e-9):
        raise ValueError(f"(T+hold)/dt 必须为整数：{total_s/dt_s:.12g}")
    steps = int(round(total_s / dt_s))
    times = np.arange(steps + 1, dtype=float) * dt_s
    end_index = int(round(duration_s / dt_s))
    mass = physics["mass_kg"]
    slm = (physics["slm_depth_j"], physics["slm_waist_m"], physics["slm_center_m"])

    def total_force(x, t):
        center = float(waveform.center(t)) + alignment_offset_m
        depth = float(waveform.depth(t))
        return (gaussian_force(x, slm[0], slm[1], slm[2])
                + gaussian_force(x, depth, physics["aod_waist_m"], center))

    def e_osc_fn(x, v, t):
        """参量通道使用的振荡能量：E_tot 相对双阱总势阱底的近似。

        阱底取两个阱心处总势的较深者（分离阱 ≈ 真实阱底，重合阱严格给出
        -(D_SLM+D_AOD)）；避免用单一阱深作底导致束缚原子 ε_osc 被钳到 0。
        """
        depth = float(waveform.depth(t))
        center = float(waveform.center(t)) + alignment_offset_m
        u_at_aod_center = (-depth
                           + float(gaussian_potential(center, slm[0], slm[1], slm[2])))
        u_at_slm_center = (-slm[0]
                           + float(gaussian_potential(slm[2], depth, physics["aod_waist_m"], center)))
        floor = min(u_at_aod_center, u_at_slm_center)
        energy_total = (0.5 * mass * v * v
                        + float(gaussian_potential(x, slm[0], slm[1], slm[2]))
                        + float(gaussian_potential(x, depth, physics["aod_waist_m"], center)))
        return max(0.0, energy_total - floor)

    x0 = states["x0_m"]
    v0 = states["v0_m_per_s"]
    near_threshold = physics["near_threshold_fraction"]
    records = []
    trajectories = {}
    for shot in range(x0.size):
        if heating is None:
            time, x, v, _a = velocity_verlet_time_dependent(
                total_force, mass, float(x0[shot]), float(v0[shot]), times)
            noise_work = None
        else:
            time, x, v, _a, noise_work = velocity_verlet_time_dependent_heating(
                total_force, mass, float(x0[shot]), float(v0[shot]), times, heating, e_osc_fn)
        x_f, v_f = float(x[end_index]), float(v[end_index])
        final_energy = 0.5 * mass * v_f**2 + float(gaussian_potential(x_f, *slm))
        captured = is_captured(final_energy, physics.get("capture_threshold", 0.0))
        margin = -final_energy / slm[0]
        record = {
            "shot_id": shot,
            "x0_um": float(x0[shot]) * 1e6, "v0_m_s": float(v0[shot]),
            "initial_aod_energy_uK": float(states["initial_aod_energy_J"][shot]) / 1.380649e-29 / 1e6,
            "initial_excitation_uK": float(states["initial_excitation_J"][shot]) / 1.380649e-29 / 1e6,
            "final_x_um": x_f * 1e6, "final_v_m_s": v_f,
            "final_slm_energy_uK": final_energy / 1.380649e-29 / 1e6,
            "captured": captured,
            "near_threshold": bool(abs(final_energy) / slm[0] < near_threshold),
            "capture_margin": margin,
            "final_excitation_over_depth": final_energy / slm[0] + 1.0,
            "excitation_change_uK": (final_energy + slm[0] - float(states["initial_excitation_J"][shot]))
                                    / 1.380649e-29 / 1e6,
            "end_time_index": end_index,
        }
        if noise_work is not None:
            record.update({
                "noise_recoil_energy_uK": float(noise_work["recoil_J"][end_index]) / 1.380649e-29 / 1e6,
                "noise_parametric_energy_uK": float(noise_work["parametric_J"][end_index]) / 1.380649e-29 / 1e6,
                "noise_pointing_energy_uK": float(noise_work["pointing_J"][end_index]) / 1.380649e-29 / 1e6,
                "noise_total_energy_uK": float(noise_work["total_J"][end_index]) / 1.380649e-29 / 1e6,
            })
        if work_energy:
            we = work_energy_along_trajectory(time[:end_index + 1], x[:end_index + 1], v[:end_index + 1],
                                              waveform, mass, slm[0], slm[1], slm[2], physics["aod_waist_m"],
                                              noise_work_j=None if noise_work is None
                                              else noise_work["total_J"][:end_index + 1])
            record.update({
                "final_move_work_uK": we["final_work_move_J"] / 1.380649e-29 / 1e6,
                "final_depth_work_uK": we["final_work_depth_J"] / 1.380649e-29 / 1e6,
                "final_total_work_uK": we["final_work_total_J"] / 1.380649e-29 / 1e6,
                "max_abs_work_energy_residual_over_depth": we["max_abs_residual_J"] / slm[0],
            })
        hold = hold_segment_energy_error(time, x, v, mass, slm[0], slm[1], slm[2], duration_s)
        record["hold_peak_to_peak_energy_error_over_depth"] = hold["peak_to_peak_energy_error_over_slm_depth"]
        records.append(record)
        if shot in keep_trajectory_ids:
            trajectory = {
                "time_s": time, "position_m": x, "velocity_m_per_s": v,
                "aod_center_m": np.asarray(waveform.center(time), dtype=float) + alignment_offset_m,
                "aod_depth_j": np.asarray(waveform.depth(time), dtype=float),
                "end_index": end_index,
            }
            if noise_work is not None:
                trajectory["noise_work"] = noise_work
            trajectories[shot] = trajectory

    captured_flags = [r["captured"] for r in records]
    margins = np.array([r["capture_margin"] for r in records])
    captured_excitation = [r["final_excitation_over_depth"] for r in records if r["captured"]]
    summary = {
        "shots": len(records),
        "captured": int(sum(captured_flags)),
        "capture_rate": float(np.mean(captured_flags)),
        "margin_quantile": float(np.quantile(margins, 0.10)),
        "mean_final_excitation_normalized_captured":
            float(np.mean(captured_excitation)) if captured_excitation else None,
        "max_final_excitation_normalized_captured":
            float(np.max(captured_excitation)) if captured_excitation else None,
        "near_threshold_count": int(sum(r["near_threshold"] for r in records)),
    }
    if heating is not None:
        summary["mean_noise_total_energy_uK"] = float(
            np.mean([r["noise_total_energy_uK"] for r in records]))
    return {"spec": spec, "records": records, "summary": summary, "trajectories": trajectories,
            "captured_flags": captured_flags}


def ideal_run(spec: dict, physics: dict, dt_s: float, hold_s: float, alignment_offset_m=0.0,
              heating: MeanHeating | None = None):
    """理想初态（AOD 阱底静止）单条轨迹及功-能核算。"""
    states = {"x0_m": np.array([physics["aod_initial_center_m"]]),
              "v0_m_per_s": np.array([0.0]),
              "initial_aod_energy_J": np.array([-physics["aod_initial_depth_j"]]),
              "initial_excitation_J": np.array([0.0])}
    result = evaluate_waveform_on_states(spec, states, physics, dt_s, hold_s,
                                         alignment_offset_m=alignment_offset_m, work_energy=True,
                                         keep_trajectory_ids=(0,), heating=heating)
    waveform = build_waveform(spec, physics)
    trajectory = result["trajectories"][0]
    noise_work_j = None
    if "noise_work" in trajectory:
        noise_work_j = trajectory["noise_work"]["total_J"]
    we = work_energy_along_trajectory(trajectory["time_s"], trajectory["position_m"],
                                      trajectory["velocity_m_per_s"], waveform, physics["mass_kg"],
                                      physics["slm_depth_j"], physics["slm_waist_m"],
                                      physics["slm_center_m"], physics["aod_waist_m"],
                                      noise_work_j=noise_work_j)
    result["work_energy"] = we
    result["waveform_name"] = spec.get("name", spec["type"])
    return result


def run_validation(cfg: Level2Config, specs: dict, states: dict, physics: dict,
                   heating: MeanHeating | None = None):
    """对多个波形在同一 validation 初态数组上评估（CRN），并做配对统计。"""
    results = {}
    for name, spec in specs.items():
        spec = dict(spec)
        spec["name"] = name
        results[name] = evaluate_waveform_parallel(
            spec, states, physics, cfg.integration.validation_dt_s,
            cfg.integration.post_transfer_hold_s, work_energy=True, heating=heating)
    names = list(specs)
    comparisons = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            flags_a = results[a]["captured_flags"]
            flags_b = results[b]["captured_flags"]
            comparisons[f"{a}_vs_{b}"] = {
                "paired_counts": paired_counts(flags_a, flags_b),
                "capture_rate_difference_bootstrap": paired_bootstrap_difference(
                    flags_a, flags_b, cfg.robustness.bootstrap_samples, cfg.robustness.bootstrap_seed + i * 7 + j),
            }
    return results, comparisons
