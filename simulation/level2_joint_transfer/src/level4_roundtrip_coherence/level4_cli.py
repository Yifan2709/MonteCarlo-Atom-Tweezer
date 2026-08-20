"""Level 4 命令行入口：preflight / single_round / repeated_rounds / ablations / sensitivity / all。"""
from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import numpy as np
import yaml

from level0_static_trap.io_utils import (ensure_output_directory, save_json,
                                         save_markdown, save_yaml, validate_pngs)
from level3_3d_transport_lensing.gaussian_3d import KB
from . import level4_visualization as viz
from . import protocol_segments as ps
from . import upstream_artifacts
from .dd_sequences import toggling_function
from .dephasing import conditional_contrast, contrast_bootstrap_ci
from .energy_ledger import segment_ledger_rows, segment_ledger_summary
from .level4_config import Level4Config, config_as_dict, load_config
from .noise_models import IntensityNoise
from .preflight4 import Preflight4
from .roundtrip_simulation import simulate_roundtrip
from .roundtrip_statistics import (evaluate_outcomes, heating_slope,
                                   paired_bootstrap_difference, paired_counts,
                                   summarize_success, survival_fits)
from .sampling4 import sample_slm_states

US = 1e-6
HBAR = 1.054571817e-34
DD_MODES = ("none", "spin_echo", "xy4_per_long_move")

SINGLE_ROUND_CONDITIONS = [
    {"name": "protocol_a_nominal_lensing", "protocol": "protocol_a",
     "lensing": True},
    {"name": "protocol_a_lensing_off", "protocol": "protocol_a",
     "lensing": False},
    {"name": "protocol_b_frozen", "protocol": "protocol_b", "lensing": True},
    {"name": "static_idle_same_duration", "protocol": "static_idle",
     "lensing": False},
]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Level 4 端到端往返协议与半经典相干性")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("preflight", "single_round",
                                            "repeated_rounds", "ablations",
                                            "sensitivity", "all"),
                        default="all")
    parser.add_argument("--output-dir", help="覆盖输出目录")
    parser.add_argument("--apply-visual-review", metavar="JSON",
                        help="只合并视觉审阅结论，不重跑任何阶段")
    return parser


def _write_rows_csv(path: Path, rows: list[dict]):
    """把 dict 行列表写成 CSV（空行保护）。"""
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({k for r in rows for k in r}, key=str)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k))
                             for k in fieldnames})


def build_protocols(cfg: Level4Config, resolved: dict) -> dict:
    """构建全部协议时间线（nominal lensing + 关闭变体）。"""
    vs = resolved["nominal_vs_m_per_s"]
    return {
        "protocol_a": ps.build_protocol_a(cfg, vs),
        "protocol_a_no_lensing": ps.build_protocol_a(cfg, None),
        "protocol_b": ps.build_protocol_b(cfg, resolved["level2_spec"], vs),
        "static_idle": ps.build_static_idle(cfg),
    }


def _pulses_for(timeline, dt, modes=DD_MODES, extra_shift_s=0.0):
    """生成各 DD 模式的脉冲集合（自动避开网格与边界）。"""
    pulses = {}
    for mode in modes:
        if mode == "none":
            pulses[mode] = np.empty(0)
            continue
        base = timeline.dd_pulse_times(mode, dt, None)
        if extra_shift_s:
            base = base + extra_shift_s
            base = np.array([ps._snap_off_grid(p, dt) for p in base])
        pulses[mode] = base
    return pulses


def _outcome_rows(res, n_shots):
    """从 RoundtripResult 提取逐 shot 结果与 checkpoint 匹配。"""
    matches = []
    names = []
    for rec in res.checkpoints:
        cl = rec["classification"]
        matches.append(np.array([b == rec["checkpoint"]["expected_basin"]
                                 for b in cl["basin"]]))
        names.append(rec["checkpoint"]["name"])
    if not matches:
        matches = [np.ones(n_shots, dtype=bool)]
        names = ["final_hold_end"]
    final_retained = matches[-1]
    out = evaluate_outcomes(matches, final_retained, checkpoint_names=names)
    return out, matches, names


def run_condition(timeline, dt, pos0, vel0, cfg, condition_name,
                  noise=None, keep_ids=(), record_stride=0, n_repeats=1,
                  phase_snap_steps=()):
    """运行一个条件的整批初态，返回结果与逐 shot 摘要。"""
    pulses = _pulses_for(timeline, dt)
    res = simulate_roundtrip(
        timeline, dt, pos0, vel0, timeline.physics, pulses, noise=noise,
        keep_ids=keep_ids, record_stride=record_stride,
        ambiguous_tol=cfg.ambiguous_basin_tolerance, n_repeats=n_repeats,
        phase_snap_steps=phase_snap_steps)
    out, matches, names = _outcome_rows(res, len(pos0))
    phases = res.phase.phases(cfg.eta_slm, cfg.eta_aod)
    return {"result": res, "outcomes": out, "matches": matches,
            "checkpoint_names": names, "phases": phases, "pulses": pulses,
            "condition": condition_name}


# ---------------------------------------------------------------- preflight
def stage_preflight(cfg, resolved, protocols, output_dir) -> dict:
    """20 项预检查；失败则抛异常停止大样本。"""
    manifest, _ = upstream_artifacts.build_manifest(cfg)
    save_json(output_dir / "upstream_manifest.json", manifest)
    started = time.time()
    pf = Preflight4(cfg, manifest, resolved, protocols)
    checks = pf.run()
    checks["elapsed_s_total"] = time.time() - started
    save_json(output_dir / "preflight.json", checks)
    if not checks["all_passed"]:
        failed = [n for n in checks["critical_checks"]
                  if not checks[n].get("passed", False)]
        raise RuntimeError(f"preflight 未通过，停止主 validation：{failed}")
    print(f"[preflight] {len(checks['critical_checks'])} 项检查全部通过"
          f"（{checks['elapsed_s_total']:.0f} s）")
    return checks


def _frozen_validation_yaml(cfg, resolved):
    """冻结的主 validation 条件（在打开 validation 前写入）。"""
    return {
        "frozen_at_stage": "before_single_round",
        "upstream_nominal_vs_m_per_s": resolved["nominal_vs_m_per_s"],
        "single_round_conditions": [dict(c) for c in SINGLE_ROUND_CONDITIONS],
        "dd_modes": list(DD_MODES),
        "eta_slm": cfg.eta_slm, "eta_aod": cfg.eta_aod,
        "validation_dt_us": cfg.validation_dt_us,
        "single_round_seed": cfg.single_round_seed,
        "single_round_shots": cfg.single_round_validation_shots,
        "repeated_round": {
            "protocols": ["protocol_a", "protocol_b"],
            "dd_modes": ["none", "xy4_per_long_move"],
            "shots": cfg.repeated_round_shots, "max_rounds": cfg.max_rounds,
            "dt_us": cfg.repeated_dt_us, "seed": cfg.repeated_round_seed},
        "long_tail": {"shots": cfg.long_tail_shots,
                      "max_rounds": cfg.long_tail_max_rounds,
                      "dt_us": cfg.repeated_dt_us, "seed": cfg.long_tail_seed},
        "note": "validation 数据只读；开发/敏感性种子与 validation 种子隔离，"
                "不得依据 validation 结果调参后重跑",
    }


# ------------------------------------------------------------ single round
def stage_single_round(cfg, resolved, protocols, output_dir, plot_records):
    """Stage B：冻结条件的 2000-shot 独立单轮 validation。"""
    dt = cfg.validation_dt_us * US
    states = sample_slm_states(cfg, cfg.single_round_seed,
                               cfg.single_round_validation_shots)
    timeline_map = {"protocol_a": protocols["protocol_a"],
                    "protocol_b": protocols["protocol_b"],
                    "static_idle": protocols["static_idle"]}
    shot_rows = []
    ckpt_rows = []
    phase_rows = []
    ledger_rows = []
    summaries = {}
    stage_counts = {}
    phase_integrals = {}
    checkpoint_stats = {}
    nominal_ckpt_records = None
    nominal_phase_breakdown = None
    for cond in SINGLE_ROUND_CONDITIONS:
        timeline = (protocols["protocol_a"] if cond["protocol"] == "protocol_a"
                    and cond["lensing"] else
                    protocols["protocol_a_no_lensing"]
                    if cond["protocol"] == "protocol_a" else
                    timeline_map[cond["protocol"]])
        started = time.time()
        run = run_condition(timeline, dt, states["pos_m"],
                            states["vel_m_per_s"], cfg, cond["name"])
        res = run["result"]
        out = run["outcomes"]
        n = len(states["pos_m"])
        summary = summarize_success(out["roundtrip_success"], cond["name"])
        summary.update({
            "final_retained": int(out["final_retained"].sum()),
            "final_retained_fraction": float(out["final_retained"].mean()),
            "recaptured": int(out["recaptured"].sum()),
            "absorbing_failure": int(out["absorbing_failure"].sum()),
            "first_failure_counts": {
                s: int((out["first_failure_stage"] == s).sum())
                for s in np.unique(out["first_failure_stage"])},
            "max_abs_residual_uK": float(np.max(np.abs(res.residual))
                                         / KB * 1e6),
        })
        summaries[cond["name"]] = summary
        stage_counts[cond["name"]] = summary["first_failure_counts"]
        # checkpoint 逐步生存
        prev = np.ones(n, dtype=bool)
        stats = []
        from level1_transfer_1d.analysis import wilson_interval
        for k, (match, name) in enumerate(zip(run["matches"],
                                              run["checkpoint_names"])):
            cond_prob = float(match[prev].mean()) if prev.any() else 0.0
            cum = float(match.mean())
            low, high = wilson_interval(int(match.sum()), n)
            stats.append({"checkpoint": name, "conditional": cond_prob,
                          "cumulative": cum,
                          "cumulative_wilson_low": float(low),
                          "cumulative_wilson_high": float(high),
                          "n_at_risk": int(prev.sum()),
                          "n_fail_here": int((prev & ~match).sum())})
            prev = prev & match
        checkpoint_stats[cond["name"]] = stats
        # 逐 shot 行
        grid = timeline.build_grid(dt)
        exc_final = (res.checkpoints[-1]["e_total"] + grid["d_slm"][-1]) \
            / KB * 1e6 if res.checkpoints else np.zeros(n)
        for i in range(n):
            row = {
                "condition": cond["name"], "shot_id": i,
                "x0_um": states["pos_m"][i, 0] * 1e6,
                "y0_um": states["pos_m"][i, 1] * 1e6,
                "z0_um": states["pos_m"][i, 2] * 1e6,
                "vx0_m_s": states["vel_m_per_s"][i, 0],
                "vy0_m_s": states["vel_m_per_s"][i, 1],
                "vz0_m_s": states["vel_m_per_s"][i, 2],
                "e0_uK": res.e0[i] / KB * 1e6,
                "xf_um": res.final_pos[i, 0] * 1e6,
                "yf_um": res.final_pos[i, 1] * 1e6,
                "zf_um": res.final_pos[i, 2] * 1e6,
                "vxf_m_s": res.final_vel[i, 0],
                "vyf_m_s": res.final_vel[i, 1],
                "vzf_m_s": res.final_vel[i, 2],
                "e_final_uK": res.e_final[i] / KB * 1e6,
                "excitation_final_uK": float(exc_final[i]),
                "final_retained": bool(out["final_retained"][i]),
                "roundtrip_success": bool(out["roundtrip_success"][i]),
                "recaptured": bool(out["recaptured"][i]),
                "first_failure_stage": str(out["first_failure_stage"][i]),
                "w_ext_total_uK": res.w_total[i] / KB * 1e6,
                "w_residual_uK": res.residual[i] / KB * 1e6,
                "max_radial_um": res.max_radial[i] * 1e6,
                "max_axial_um": res.max_axial[i] * 1e6,
            }
            for mode in DD_MODES:
                row[f"phase_{mode}_rad"] = float(run["phases"][mode][i])
            shot_rows.append(row)
        # checkpoint 状态行
        for rec in res.checkpoints:
            cl = rec["classification"]
            ck = rec["checkpoint"]
            for i in range(n):
                ckpt_rows.append({
                    "condition": cond["name"], "shot_id": i,
                    "checkpoint": ck["name"], "expected_basin":
                        ck["expected_basin"], "basin": cl["basin"][i],
                    "x_um": rec["pos"][i, 0] * 1e6,
                    "y_um": rec["pos"][i, 1] * 1e6,
                    "z_um": rec["pos"][i, 2] * 1e6,
                    "vx_m_s": rec["vel"][i, 0], "vy_m_s": rec["vel"][i, 1],
                    "vz_m_s": rec["vel"][i, 2],
                    "e_slm_uK": cl["e_slm"][i] / KB * 1e6,
                    "e_aod_uK": cl["e_aod"][i] / KB * 1e6,
                    "e_total_uK": rec["e_total"][i] / KB * 1e6,
                    "near_threshold": bool(cl["near_slm_threshold"][i]
                                           or cl["near_aod_threshold"][i]),
                    "ambiguous": bool(cl["ambiguous"][i]),
                    "cum_work_uK": rec["cum_work"][i] / KB * 1e6,
                    "separation_um": cl["separation_m"] * 1e6,
                })
        # 相位记录（每 shot × DD × 分段 + 总相位）
        seg_phases = res.phase.segment_phases(cfg.eta_slm, cfg.eta_aod)
        seg_names = grid["segment_names"]
        if cond["name"] == "protocol_a_nominal_lensing":
            # 供 basin/能量演化/分段相位三张诊断图复用（仅 nominal 主条件）
            nominal_ckpt_records = res.checkpoints
            nominal_phase_breakdown = {
                "seg_phases": seg_phases, "seg_names": seg_names,
                "success": out["roundtrip_success"]}
        for mode in DD_MODES:
            seg_mat = seg_phases[mode]
            for i in range(n):
                for s in range(seg_mat.shape[0]):
                    phase_rows.append({
                        "condition": cond["name"], "shot_id": i, "round": 1,
                        "dd_mode": mode, "segment": seg_names[s],
                        "phase_rad": float(seg_mat[s, i]),
                        "phase_total_rad": float(run["phases"][mode][i]),
                    })
        # 功-能账本
        rows, seg_works = segment_ledger_rows(res, grid, cond["name"])
        ledger_rows.extend(rows)
        phase_integrals[cond["name"]] = {
            "a_slm": res.phase.state[DD_MODES[0]]["a_slm"].copy(),
            "a_aod": res.phase.state[DD_MODES[0]]["a_aod"].copy(),
        }
        for mode in DD_MODES:
            phase_integrals[cond["name"]][f"a_slm_{mode}"] = \
                res.phase.state[mode]["a_slm"].copy()
            phase_integrals[cond["name"]][f"a_aod_{mode}"] = \
                res.phase.state[mode]["a_aod"].copy()
        print(f"[single_round] {cond['name']}: "
              f"{summary['successes']}/{summary['shots']} "
              f"({time.time() - started:.0f}s)")
    _write_rows_csv(output_dir / "single_round_shots.csv", shot_rows)
    _write_rows_csv(output_dir / "checkpoint_states.csv", ckpt_rows)
    _write_rows_csv(output_dir / "phase_records.csv", phase_rows)
    _write_rows_csv(output_dir / "energy_ledger.csv", ledger_rows)

    # 配对比较
    def flags(name, key="roundtrip_success"):
        return np.array([r[key] for r in shot_rows if r["condition"] == name])

    comparisons = {}
    pairs = [
        ("lensing_off_vs_nominal", "protocol_a_lensing_off",
         "protocol_a_nominal_lensing"),
        ("protocol_b_vs_a", "protocol_b_frozen", "protocol_a_nominal_lensing"),
        ("static_idle_vs_a", "static_idle_same_duration",
         "protocol_a_nominal_lensing"),
    ]
    for name, a, b in pairs:
        fa, fb = flags(a), flags(b)
        comparisons[name] = {
            "pair": (a, b),
            "paired_counts": paired_counts(fa, fb),
            "bootstrap": paired_bootstrap_difference(fa, fb,
                                                     cfg.bootstrap_samples,
                                                     cfg.bootstrap_seed),
        }
    # 对比度统计
    contrast = {}
    for cond in SINGLE_ROUND_CONDITIONS:
        mask = flags(cond["name"])
        for mode in DD_MODES:
            phi = np.array([r[f"phase_{mode}_rad"] for r in shot_rows
                            if r["condition"] == cond["name"]])
            stats = conditional_contrast(phi, mask)
            stats["mode"] = mode
            stats["condition"] = cond["name"]
            stats["bootstrap_ci"] = contrast_bootstrap_ci(
                phi, mask, n_boot=2000, seed=cfg.bootstrap_seed)
            contrast[f"{cond['name']}|{mode}"] = stats
    # 压力池：采集失败阶段代表轨迹（诊断种子，非 validation）
    stress = _collect_stress_representatives(cfg, protocols, output_dir)
    # 图
    plot_records.append(viz.plot_checkpoint_survival(
        checkpoint_stats, output_dir / "checkpoint_survival.png"))
    plot_records.append(viz.plot_failure_stage_breakdown(
        stage_counts, output_dir / "failure_stage_breakdown.png"))
    plot_records.append(viz.plot_paired_comparisons(
        comparisons, output_dir / "paired_comparisons.png"))
    if nominal_ckpt_records:
        plot_records.append(viz.plot_basin_classification(
            nominal_ckpt_records, output_dir / "basin_classification.png"))
        plot_records.append(viz.plot_checkpoint_energy_evolution(
            nominal_ckpt_records,
            output_dir / "checkpoint_energy_evolution.png"))
    if nominal_phase_breakdown:
        plot_records.append(viz.plot_segment_phase_breakdown(
            nominal_phase_breakdown["seg_phases"],
            nominal_phase_breakdown["seg_names"],
            nominal_phase_breakdown["success"],
            output_dir / "segment_phase_breakdown.png"))
    if stress.get("trajs"):
        plot_records.append(viz.plot_representative_roundtrips(
            stress["trajs"], output_dir / "representative_roundtrips.png"))
    # DD toggling 与代表相位轨迹
    toggling_data = _representative_toggling(cfg, protocols["protocol_a"])
    if toggling_data:
        plot_records.append(viz.plot_dd_toggling(
            toggling_data, output_dir / "dd_toggling_and_phase.png"))
    return {"summaries": summaries, "comparisons": comparisons,
            "contrast": contrast, "checkpoint_stats": checkpoint_stats,
            "stage_counts": stage_counts,
            "initial_states_meta": {k: v for k, v in states.items()
                                    if not isinstance(v, np.ndarray)},
            "phase_integrals": phase_integrals,
            "stress": {k: v for k, v in stress.items() if k != "trajs"}}


def _representative_toggling(cfg, timeline, shot_seed=None):
    """单个代表轨迹的三种 DD 模式 y(t)、δω 的 SLM/AOD 分解与累计 φ。"""
    dt = 0.10 * US
    states = sample_slm_states(cfg, cfg.development_seed + 400, 1)
    pulses = _pulses_for(timeline, dt)
    res = simulate_roundtrip(timeline, dt, states["pos_m"],
                             states["vel_m_per_s"], timeline.physics,
                             pulses, keep_ids=[0], record_stride=20)
    if not res.trajectory:
        return None
    from scipy.integrate import cumulative_trapezoid
    t = np.asarray(res.trajectory["time_s"])
    pos = np.asarray(res.trajectory["position_m"])[:, 0, :]
    grid = timeline.build_grid(dt)
    idx = np.searchsorted(grid["times"], t)
    d_slm = grid["d_slm"][idx]
    d_aod = grid["d_aod"][idx]
    from level3_3d_transport_lensing.gaussian_3d import static_potential
    from level3_3d_transport_lensing.aod_lensing import \
        lensing_potential_and_force
    u_slm = static_potential(pos[:, 0], pos[:, 1], pos[:, 2], d_slm,
                             timeline.physics.slm_waist_m,
                             timeline.physics.slm_zr_m)
    u_aod = lensing_potential_and_force(
        pos[:, 0] - grid["c1"][idx], pos[:, 1] - grid["c2"][idx], pos[:, 2],
        d_aod, timeline.physics.aod_waist_m, timeline.physics.aod_zr_m,
        grid["zs1"][idx], grid["zs2"][idx])[0]
    delta_omega = (cfg.eta_slm * u_slm + cfg.eta_aod * u_aod) / HBAR
    modes = {}
    for mode in DD_MODES:
        y = toggling_function(t, pulses[mode])
        modes[mode] = {"y": y,
                       "cumulative_phi": cumulative_trapezoid(
                           y * delta_omega, t, initial=0.0),
                       "pulses_s": pulses[mode]}
    return {"time_s": t, "modes": modes,
            "delta_omega": delta_omega,
            "delta_omega_slm": cfg.eta_slm * u_slm / HBAR,
            "delta_omega_aod": cfg.eta_aod * u_aod / HBAR}


def _collect_stress_representatives(cfg, protocols, output_dir):
    """高温压力池：为 representative_roundtrips 提供失败/recaptured 代表。"""
    dt = cfg.repeated_dt_us * US
    states = sample_slm_states(cfg, cfg.development_seed + 300, 512,
                               temperature_uK=50.0)
    timeline = protocols["protocol_a"]
    keep = list(range(6))
    run = run_condition(timeline, dt, states["pos_m"], states["vel_m_per_s"],
                        cfg, "stress_diagnostic_50uK", keep_ids=keep,
                        record_stride=max(1, int(round(0.5 * US / dt))))
    out = run["outcomes"]
    res = run["result"]
    rep_rows = []
    picked = {}
    for i in range(len(out["roundtrip_success"])):
        stage = str(out["first_failure_stage"][i])
        key = ("success" if out["roundtrip_success"][i] else
               f"fail_{stage}" if not out["recaptured"][i] else "recaptured")
        picked.setdefault(key, [])
        if len(picked[key]) < cfg.representative_shots_per_failure_stage:
            picked[key].append(i)
    for key, ids in picked.items():
        for rank, i in enumerate(ids):
            rep_rows.append({"class": key, "shot_id": i, "rank": rank,
                             "condition": "stress_diagnostic_50uK"})
    # 需要绘制轨迹的类别 → 选其在批次内的实际 shot id（保证失败/重俘获
    # 轨迹一定入图，而非固定画前 6 个成功 shot）
    wanted = [i for ids in picked.values() for i in ids]
    kept = [i for i in wanted if i in keep]
    for i in wanted:
        if len(kept) >= len(keep):
            break
        if i not in kept:
            kept.append(i)
    if kept != keep:
        # 所选轨迹 id 与前 6 个不同 → 用同一初态再跑一次小批，只记录所选
        run = run_condition(timeline, dt, states["pos_m"],
                            states["vel_m_per_s"], cfg,
                            "stress_diagnostic_50uK", keep_ids=kept,
                            record_stride=max(1, int(round(0.5 * US / dt))))
        res = run["result"]
    # 轨迹时序转储（成功 + 失败代表）
    traj_rows = []
    if res.trajectory:
        t = np.asarray(res.trajectory["time_s"])
        pos = np.asarray(res.trajectory["position_m"])
        vel = np.asarray(res.trajectory["velocity_m_per_s"])
        w = np.asarray(res.trajectory["w_total_j"])
        e = np.asarray(res.trajectory["energy_j"])
        for slot, i in enumerate(kept):
            label = None
            for key, ids in picked.items():
                if i in ids:
                    label = key
            label = label or "success"
            for k in range(len(t)):
                traj_rows.append({
                    "condition": "stress_diagnostic_50uK", "class": label,
                    "shot_id": i, "time_us": t[k] * 1e6,
                    "x_um": pos[k, slot, 0] * 1e6,
                    "y_um": pos[k, slot, 1] * 1e6, "z_um": pos[k, slot, 2] * 1e6,
                    "vx_m_s": vel[k, slot, 0], "vy_m_s": vel[k, slot, 1],
                    "vz_m_s": vel[k, slot, 2],
                    "w_total_uK": w[k, slot] / KB * 1e6,
                    "energy_uK": e[k, slot] / KB * 1e6})
    _write_rows_csv(output_dir / "representative_trajectories.csv",
                    traj_rows or rep_rows)
    return {"outcome_counts": {k: len(v) for k, v in picked.items()},
            "success_rate": float(np.mean(out["roundtrip_success"])),
            "temperature_uK": 50.0,
            "trajs": _stress_traj_map(res, kept, picked),
            "note": "高温压力池仅用于采集失败阶段代表轨迹（诊断种子池，"
                    "非 validation 数据）"}


def _stress_traj_map(res, keep, picked):
    """把 keep_ids 的轨迹整理成 {标签: time/pos/energy/work} 映射。"""
    if not res.trajectory:
        return {}
    t = np.asarray(res.trajectory["time_s"])
    pos = np.asarray(res.trajectory["position_m"])
    energy = np.asarray(res.trajectory["energy_j"])
    work = np.asarray(res.trajectory["w_total_j"])
    out = {}
    for slot, i in enumerate(keep):
        label = None
        for key, ids in picked.items():
            if i in ids:
                label = key
        label = label or "success"
        out[f"shot{i}({label})"] = {
            "time_s": t,
            "position_m": pos[:, slot:slot + 1, :],
            "energy_uK": energy[:, slot] / KB * 1e6,
            "work_uK": work[:, slot] / KB * 1e6}
    return out


# --------------------------------------------------------- repeated rounds
def stage_repeated_rounds(cfg, resolved, protocols, output_dir, plot_records):
    """Stage C：500-shot/10-round 重复往返 + 128-shot/40-round 长尾。"""
    dt = cfg.repeated_dt_us * US
    rounds = cfg.max_rounds
    metrics_rows = []
    failure_rows = []
    heating = {}
    contrast_curve = {}
    survival = {}
    fits = {}
    long_tail = {}
    for pname in ("protocol_a", "protocol_b"):
        timeline = protocols[pname]
        states = sample_slm_states(cfg, cfg.repeated_round_seed,
                                   cfg.repeated_round_shots)
        n = cfg.repeated_round_shots
        n_steps_one = int(round(timeline.total_duration_s / dt))
        snaps = {r * n_steps_one: r for r in range(1, rounds + 1)}
        started = time.time()
        run = run_condition(timeline, dt, states["pos_m"], states["vel_m_per_s"],
                            cfg, f"{pname}|repeated", n_repeats=rounds,
                            phase_snap_steps=snaps)
        res = run["result"]
        # 按轮提取 checkpoint
        ck_by_round = {}
        for rec in res.checkpoints:
            ck_by_round.setdefault(rec["checkpoint"]["round"], []).append(rec)
        ever_ok = np.ones(n, dtype=bool)
        absorb_curve, absorb_low, absorb_high = [], [], []
        nonabs_curve = []
        heat_med, heat_q10, heat_q90, heat_q95 = [], [], [], []
        contrast_pts = {}
        for mode in DD_MODES:
            contrast_pts[mode] = []
        d_slm_bottom = -cfg.slm_initial_depth_uK * 1e-6 * KB
        stage_cum = {}
        for r in range(1, rounds + 1):
            recs = ck_by_round.get(r, [])
            matches = [np.array([b == rec["checkpoint"]["expected_basin"]
                                 for b in rec["classification"]["basin"]])
                       for rec in recs]
            names = [rec["checkpoint"]["name"] for rec in recs]
            out = evaluate_outcomes(matches, matches[-1],
                                    checkpoint_names=names)
            ever_ok = ever_ok & out["roundtrip_success"]
            # 本轮新增失败阶段的累计（吸收语义：每 shot 只记一次）
            for stage in np.unique(out["first_failure_stage"]):
                if stage == "none":
                    continue
                stage_cum[str(stage)] = stage_cum.get(str(stage), 0) + int(
                    (out["first_failure_stage"] == stage).sum())
            # 吸收成功 = 截至该轮无任何失败
            absorb = ever_ok
            s = summarize_success(absorb, f"{pname}|round{r}")
            absorb_curve.append(s["probability"])
            absorb_low.append(s["wilson_low"])
            absorb_high.append(s["wilson_high"])
            final_ret = matches[-1]
            nonabs_curve.append(float(final_ret.mean()))
            exc = (recs[-1]["e_total"] - d_slm_bottom) / KB * 1e6
            exc_retained = exc[final_ret]
            heat_med.append(float(np.median(exc_retained))
                            if exc_retained.size else None)
            heat_q10.append(float(np.quantile(exc_retained, 0.10))
                            if exc_retained.size else None)
            heat_q90.append(float(np.quantile(exc_retained, 0.90))
                            if exc_retained.size else None)
            heat_q95.append(float(np.quantile(exc_retained, 0.95))
                            if exc_retained.size else None)
            # 相位快照 → contrast
            snap = res.phase_snaps.get(r)
            phi_mode = {}
            for mode in DD_MODES:
                phi = ((cfg.eta_slm * snap[mode]["a_slm"]
                        + cfg.eta_aod * snap[mode]["a_aod"]) / HBAR)
                phi_mode[mode] = phi
                stats = conditional_contrast(phi, absorb)
                contrast_pts[mode].append(stats["contrast"])
            for i in np.where(~out["roundtrip_success"])[0]:
                failure_rows.append({
                    "protocol": pname, "round": r, "shot_id": int(i),
                    "checkpoint": str(out["first_failure_checkpoint"][i]),
                    "stage": str(out["first_failure_stage"][i]),
                    "recaptured_later": bool(final_ret[i]),
                })
            near = np.array([rec["classification"]["ambiguous"][i]
                             for rec in recs for i in range(n)]).mean()
            metrics_rows.append({
                "protocol": pname, "round": r,
                "risk_set": int(ever_ok.sum()),
                "absorbing_success": int(absorb.sum()),
                "absorbing_fraction": s["probability"],
                "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
                "non_absorbing_retained": int(final_ret.sum()),
                "non_absorbing_fraction": float(final_ret.mean()),
                "first_failure_counts": json.dumps(stage_cum, ensure_ascii=False),
                "excitation_median_uK": heat_med[-1],
                "excitation_q90_uK": heat_q90[-1],
                "excitation_q95_uK": heat_q95[-1],
                "contrast_none": contrast_pts["none"][-1],
                "contrast_spin_echo": contrast_pts["spin_echo"][-1],
                "contrast_xy4": contrast_pts["xy4_per_long_move"][-1],
                "usable_coherent_fraction_xy4":
                    (contrast_pts["xy4_per_long_move"][-1] or 0.0)
                    * s["probability"],
                "ambiguous_fraction_mean": float(near),
            })
        rlist = list(range(1, rounds + 1))
        survival[f"{pname}|none"] = {"rounds": rlist,
                                     "absorbing": absorb_curve,
                                     "absorbing_low": absorb_low,
                                     "absorbing_high": absorb_high,
                                     "non_absorbing": nonabs_curve}
        heating[f"{pname}"] = {"rounds": rlist, "median": heat_med,
                               "q10": heat_q10, "q90": heat_q90,
                               "q95": heat_q95}
        for mode in DD_MODES:
            contrast_curve[f"{pname}|{mode}"] = {"rounds": rlist,
                                                 "contrast": contrast_pts[mode]}
        fits[f"{pname}|none"] = survival_fits(rlist, absorb_curve)
        print(f"[repeated] {pname}: absorbing {absorb_curve[-1]:.4f}"
              f" non-absorbing {nonabs_curve[-1]:.4f}"
              f"（{time.time() - started:.0f}s）")
    # 长尾
    timeline = protocols["protocol_a"]
    states = sample_slm_states(cfg, cfg.long_tail_seed, cfg.long_tail_shots)
    n = cfg.long_tail_shots
    n_steps_one = int(round(timeline.total_duration_s / dt))
    snaps = {r * n_steps_one: r
             for r in cfg.checkpoint_rounds if r <= cfg.long_tail_max_rounds}
    started = time.time()
    run = run_condition(timeline, dt, states["pos_m"], states["vel_m_per_s"],
                        cfg, "long_tail", n_repeats=cfg.long_tail_max_rounds,
                        phase_snap_steps=snaps)
    res = run["result"]
    ck_by_round = {}
    for rec in res.checkpoints:
        ck_by_round.setdefault(rec["checkpoint"]["round"], []).append(rec)
    ever_ok = np.ones(n, dtype=bool)
    tail_curve = []
    d_slm_bottom = -cfg.slm_initial_depth_uK * 1e-6 * KB
    for r in sorted(ck_by_round):
        recs = ck_by_round[r]
        matches = [np.array([b == rec["checkpoint"]["expected_basin"]
                             for b in rec["classification"]["basin"]])
                   for rec in recs]
        names = [rec["checkpoint"]["name"] for rec in recs]
        out = evaluate_outcomes(matches, matches[-1], checkpoint_names=names)
        ever_ok = ever_ok & out["roundtrip_success"]
        exc = (recs[-1]["e_total"] - d_slm_bottom) / KB * 1e6
        exc_ret = exc[matches[-1]]
        snap = res.phase_snaps.get(r)
        row = {
            "protocol": "protocol_a_long_tail", "round": r,
            "risk_set": int(ever_ok.sum()),
            "absorbing_success": int(ever_ok.sum()),
            "absorbing_fraction": float(ever_ok.mean()),
            "non_absorbing_retained": int(matches[-1].sum()),
            "non_absorbing_fraction": float(matches[-1].mean()),
            "excitation_median_uK": float(np.median(exc_ret))
            if exc_ret.size else None,
            "excitation_q90_uK": float(np.quantile(exc_ret, 0.9))
            if exc_ret.size else None,
        }
        if snap is not None:
            phi_xy4 = ((cfg.eta_slm * snap["xy4_per_long_move"]["a_slm"]
                        + cfg.eta_aod * snap["xy4_per_long_move"]["a_aod"])
                       / HBAR)
            stats = conditional_contrast(phi_xy4, ever_ok)
            row["contrast_xy4"] = stats["contrast"]
        metrics_rows.append(row)
        tail_curve.append(float(ever_ok.mean()))
    long_tail = {"rounds": sorted(ck_by_round), "survival": tail_curve,
                 "shots": n, "max_rounds": cfg.long_tail_max_rounds,
                 "dt_us": cfg.repeated_dt_us}
    print(f"[long_tail] 40-round absorbing survival "
          f"{tail_curve[-1]:.4f}（{time.time() - started:.0f}s）")
    _write_rows_csv(output_dir / "repeated_round_metrics.csv", metrics_rows)
    _write_rows_csv(output_dir / "failure_events.csv", failure_rows)
    # 加热斜率
    slopes = {}
    for pname in heating:
        finite = [(r, m) for r, m in zip(heating[pname]["rounds"],
                                          heating[pname]["median"]) if m]
        slopes[pname] = heating_slope([r for r, _ in finite],
                                      [m for _, m in finite])
    plot_records.append(viz.plot_survival_vs_round(
        survival, fits, output_dir / "survival_vs_round.png"))
    plot_records.append(viz.plot_heating_vs_round(
        heating, output_dir / "heating_vs_round.png"))
    if failure_rows:
        plot_records.append(viz.plot_failure_events_by_round(
            failure_rows, ("protocol_a", "protocol_b"),
            output_dir / "failure_events_by_round.png"))
    plot_records.append(viz.plot_long_tail_survival(
        long_tail, output_dir / "long_tail_survival.png"))
    if any(v is not None for data in contrast_curve.values()
           for v in data["contrast"]):
        plot_records.append(viz.plot_usable_coherent_fraction(
            contrast_curve, survival,
            output_dir / "usable_coherent_fraction.png"))
    return {"metrics_rows": metrics_rows, "survival": survival,
            "heating": heating, "heating_slopes": slopes, "fits": fits,
            "contrast_curve": contrast_curve, "long_tail": long_tail}


# -------------------------------------------------------------- ablations
def stage_ablations(cfg, resolved, protocols, output_dir, plot_records,
                    single_round=None):
    """Stage D：消融（归因用，不重选主协议）。"""
    dt = cfg.validation_dt_us * US
    states = sample_slm_states(cfg, cfg.sensitivity_seed, cfg.sensitivity_shots)
    pos0, vel0 = states["pos_m"], states["vel_m_per_s"]
    n = cfg.sensitivity_shots
    vs = resolved["nominal_vs_m_per_s"]
    conditions = [
        ("transfer_only", ps.build_transfer_only(cfg, vs)),
        ("transport_only", ps.build_transport_only(cfg, vs)),
        ("protocol_a_no_remote_hold", ps.build_protocol_a_no_remote_hold(cfg, vs)),
        ("protocol_a_constant_jerk_long_move",
         ps.build_protocol_a(cfg, vs, long_profile="constant_jerk")),
    ]
    rows = []
    for name, timeline in conditions:
        started = time.time()
        run = run_condition(timeline, dt, pos0, vel0, cfg, name)
        out = run["outcomes"]
        s = summarize_success(out["roundtrip_success"], name)
        exc = (run["result"].checkpoints[-1]["e_total"]
               + max(cfg.slm_initial_depth_uK, 0) * 1e-6 * KB) / KB * 1e6
        rows.append({
            "condition": name, "shots": s["shots"],
            "successes": s["successes"], "success_fraction": s["probability"],
            "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
            "final_retained_fraction": float(out["final_retained"].mean()),
            "recaptured": int(out["recaptured"].sum()),
            "median_excitation_uK": float(np.median(exc)),
            "contrast_none": conditional_contrast(
                run["phases"]["none"], out["roundtrip_success"])["contrast"],
            "contrast_xy4": conditional_contrast(
                run["phases"]["xy4_per_long_move"],
                out["roundtrip_success"])["contrast"],
            "max_abs_residual_uK": float(np.max(np.abs(
                run["result"].residual)) / KB * 1e6),
        })
        print(f"[ablation] {name}: {s['successes']}/{s['shots']}"
              f"（{time.time() - started:.0f}s）")
    # 无透镜 / static idle 复用单轮 validation 结果
    if single_round:
        for src, label in (("protocol_a_lensing_off", "no_lensing"),
                           ("static_idle_same_duration", "static_idle")):
            s = single_round["summaries"][src]
            rows.append({
                "condition": label, "shots": s["shots"],
                "successes": s["successes"],
                "success_fraction": s["probability"],
                "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
                "final_retained_fraction": s["final_retained_fraction"],
                "recaptured": s["recaptured"],
                "note": "复用 single_round validation 结果"})
        # no_differential_light_shift：经典结果同 nominal，相位为 0
        s = single_round["summaries"]["protocol_a_nominal_lensing"]
        rows.append({
            "condition": "no_differential_light_shift", "shots": s["shots"],
            "successes": s["successes"], "success_fraction": s["probability"],
            "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
            "final_retained_fraction": s["final_retained_fraction"],
            "recaptured": s["recaptured"], "contrast_none": 1.0,
            "contrast_xy4": 1.0,
            "note": "η=0 时相位恒为 0；经典动力学与 nominal 完全相同"})
        s = single_round["summaries"]["protocol_a_nominal_lensing"]
        rows.append({
            "condition": "no_intensity_noise", "shots": s["shots"],
            "successes": s["successes"], "success_fraction": s["probability"],
            "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
            "final_retained_fraction": s["final_retained_fraction"],
            "recaptured": s["recaptured"],
            "note": "nominal 运行默认不加入随机技术噪声"})
    else:
        timeline = protocols["protocol_a_no_lensing"]
        run = run_condition(timeline, dt, pos0, vel0, cfg, "no_lensing")
        s = summarize_success(run["outcomes"]["roundtrip_success"], "no_lensing")
        rows.append({"condition": "no_lensing", "shots": s["shots"],
                     "successes": s["successes"],
                     "success_fraction": s["probability"],
                     "wilson_low": s["wilson_low"],
                     "wilson_high": s["wilson_high"]})
    _write_rows_csv(output_dir / "ablations.csv", rows)
    plot_records.append(viz.plot_ablation_comparison(
        rows, output_dir / "ablation_comparison.png"))
    return {"rows": rows}


# ------------------------------------------------------------ sensitivity
def _noise_demo_traces(cfg, n_shots=256, n_show=6, duration_us=2000.0):
    """为噪声诊断图生成 quasistatic 与 OU 示例轨迹（assumed 参数）。"""
    rms_grid = [r for r in cfg.noise_rms_grid if r and r > 0]
    rms = float(np.median(rms_grid)) if rms_grid else 0.02
    tau_grid = [t for t in cfg.noise_tau_us_grid if t]
    tau_us = float(np.median(tau_grid)) if tau_grid else 50.0
    dt = 1.0 * US
    times = np.arange(0.0, duration_us * US + 1e-12, dt)
    ou = IntensityNoise("ornstein_uhlenbeck", n_shots, rms,
                        tau_s=tau_us * US,
                        cross_correlation=cfg.noise_cross_correlation,
                        seed=cfg.noise_seed + 77)
    eps_slm = [ou.eps_slm[:n_show].copy()]
    eps_aod = [ou.eps_aod[:n_show].copy()]
    for _ in range(len(times) - 1):
        ou.advance(dt)
        eps_slm.append(ou.eps_slm[:n_show].copy())
        eps_aod.append(ou.eps_aod[:n_show].copy())
    quasi = IntensityNoise("quasistatic", n_shots, rms, tau_s=None,
                           cross_correlation=cfg.noise_cross_correlation,
                           seed=cfg.noise_seed + 78)
    return {
        "ou": {"times_s": times, "eps_slm": np.stack(eps_slm),
               "eps_aod": np.stack(eps_aod), "rms": rms,
               "tau_s": tau_us * US},
        "quasi": {"eps_slm": quasi.eps_slm, "eps_aod": quasi.eps_aod,
                  "rms": rms},
    }


def stage_sensitivity(cfg, resolved, protocols, output_dir, plot_records,
                      single_round=None):
    """Stage E：one-factor-at-a-time 敏感性（冻结主协议后运行）。"""
    dt = cfg.repeated_dt_us * US
    rows = []
    eta_free = []
    plot_records.append(viz.plot_noise_diagnostic(
        _noise_demo_traces(cfg), output_dir / "noise_model_diagnostic.png"))

    def add_row(factor, value, run, note="", contrast=None):
        out = run["outcomes"]
        s = summarize_success(out["roundtrip_success"], f"{factor}={value}")
        rows.append({
            "factor": factor, "value": value, "shots": s["shots"],
            "successes": s["successes"], "success_fraction": s["probability"],
            "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
            "recaptured": int(out["recaptured"].sum()),
            "contrast_none": conditional_contrast(
                run["phases"]["none"], out["roundtrip_success"])["contrast"],
            "contrast_xy4": conditional_contrast(
                run["phases"]["xy4_per_long_move"],
                out["roundtrip_success"])["contrast"],
            "note": note,
        })

    vs = resolved["nominal_vs_m_per_s"]
    base_states = sample_slm_states(cfg, cfg.sensitivity_seed,
                                    cfg.sensitivity_shots)
    # 温度
    for temp in (3.0, 5.0, 10.0):
        states = sample_slm_states(cfg, cfg.sensitivity_seed,
                                   cfg.sensitivity_shots, temperature_uK=temp)
        run = run_condition(protocols["protocol_a"], dt, states["pos_m"],
                            states["vel_m_per_s"], cfg, f"temp={temp}")
        add_row("temperature_uK", temp, run)
        print(f"[sensitivity] T={temp}μK: "
              f"{run['outcomes']['roundtrip_success'].mean():.4f}")
    # 深度
    for factor, values, kwargs in (
            ("slm_initial_depth_uK", (162.0, 198.0), {}),
            ("slm_transport_depth_uK", (54.0, 66.0), {}),
            ("aod_depth_uK", (252.0, 308.0), {})):
        for v in values:
            variant = replace(cfg)
            if factor == "slm_initial_depth_uK":
                variant = replace(cfg, slm_initial_depth_uK=v)
            elif factor == "slm_transport_depth_uK":
                variant = replace(cfg, slm_transport_depth_uK=v)
            else:
                variant = replace(cfg, aod_transport_depth_uK=v)
            timeline = ps.build_protocol_a(variant, vs)
            run = run_condition(timeline, dt, base_states["pos_m"],
                                base_states["vel_m_per_s"], cfg, f"{factor}={v}")
            add_row(factor, v, run)
            print(f"[sensitivity] {factor}={v}: "
                  f"{run['outcomes']['roundtrip_success'].mean():.4f}")
    # 对准（SLM 中心横向偏移）
    for axis, offsets in ((0, (-0.10, -0.05, 0.05, 0.10)),
                          (1, (-0.10, -0.05, 0.05, 0.10))):
        for off in offsets:
            center = [0.0, 0.0]
            center[axis] = off
            variant = replace(cfg, slm_center_um=tuple(center))
            timeline = ps.build_protocol_a(variant, vs)
            run = run_condition(timeline, dt, base_states["pos_m"],
                                base_states["vel_m_per_s"], cfg,
                                f"align_{off}")
            add_row(f"alignment_axis{axis + 1}_um", off, run)
            print(f"[sensitivity] axis{axis + 1} offset {off}μm: "
                  f"{run['outcomes']['roundtrip_success'].mean():.4f}")
    # lensing 严重度（Level 3 冻结相邻场景）
    for sev, v_s in zip(("0.75", "1.25"), cfg.sensitivity_vs_m_per_s):
        timeline = ps.build_protocol_a(cfg, v_s)
        run = run_condition(timeline, dt, base_states["pos_m"],
                            base_states["vel_m_per_s"], cfg, f"vs={v_s:.4f}")
        add_row("lensing_vs_m_per_s", f"sev{sev}", run,
                note=f"Level 3 冻结场景 severity {sev}")
        print(f"[sensitivity] v_s(sev{sev})={v_s:.4f}: "
              f"{run['outcomes']['roundtrip_success'].mean():.4f}")
    # η（由 nominal validation 的 A_S/A_A 线性重算，免重积分）
    if single_round and "protocol_a_nominal_lensing" in \
            single_round["phase_integrals"]:
        integ = single_round["phase_integrals"]["protocol_a_nominal_lensing"]
        for eta in (1.2e-4, 1.3e-4, 1.5e-4, 1.6e-4):
            phi_none = (eta * integ["a_slm_none"] + eta * integ["a_aod_none"]) \
                / HBAR
            phi_xy4 = (eta * integ["a_slm_xy4_per_long_move"]
                       + eta * integ["a_aod_xy4_per_long_move"]) / HBAR
            c_none = conditional_contrast(phi_none)["contrast"]
            c_xy4 = conditional_contrast(phi_xy4)["contrast"]
            rows.append({
                "factor": "eta", "value": eta, "shots": len(phi_none),
                "successes": None, "success_fraction": None,
                "wilson_low": None, "wilson_high": None,
                "contrast_none": c_none, "contrast_xy4": c_xy4,
                "note": "相位线性重算（全体 shots；经典成功率不随 η 变化，"
                        "此处 C 未按 success 条件化，与 validation 的 "
                        "C_cond 定义不同，仅作 η 标度检查）"})
            eta_free.append({"eta": eta, "contrast_none": c_none,
                             "contrast_xy4": c_xy4})
    # 技术噪声网格（assumed）
    idx = 0
    for rms in cfg.noise_rms_grid:
        for tau_us in cfg.noise_tau_us_grid:
            idx += 1
            noise = IntensityNoise(
                cfg.noise_model, cfg.sensitivity_shots, rms,
                tau_s=None if tau_us is None else tau_us * US,
                cross_correlation=cfg.noise_cross_correlation,
                seed=cfg.noise_seed + idx)
            run = run_condition(protocols["protocol_a"], dt,
                                base_states["pos_m"], base_states["vel_m_per_s"],
                                cfg, f"noise={rms}@{tau_us}", noise=noise)
            label = "quasistatic" if tau_us is None else f"tau={tau_us}us"
            add_row("intensity_noise", f"rms={rms},{label}", run,
                    note="assumed_sensitivity_only；rms=0.05 为夸张化演示")
            print(f"[sensitivity] noise rms={rms} τ={tau_us}: "
                  f"contrast_none="
                  f"{conditional_contrast(run['phases']['none'], run['outcomes']['roundtrip_success'])['contrast']:.4f}")
    # DD pulse timing offset
    for shift_um in (-0.5, 0.5):
        pulses = _pulses_for(protocols["protocol_a"], dt,
                             modes=("none", "spin_echo", "xy4_per_long_move"),
                             extra_shift_s=shift_um * 1e-6)
        res = simulate_roundtrip(protocols["protocol_a"], dt,
                                 base_states["pos_m"], base_states["vel_m_per_s"],
                                 protocols["protocol_a"].physics, pulses,
                                 ambiguous_tol=cfg.ambiguous_basin_tolerance)
        matches = [np.array([b == rec["checkpoint"]["expected_basin"]
                             for b in rec["classification"]["basin"]])
                   for rec in res.checkpoints]
        names = [rec["checkpoint"]["name"] for rec in res.checkpoints]
        out = evaluate_outcomes(matches, matches[-1], checkpoint_names=names)
        phases = res.phase.phases(cfg.eta_slm, cfg.eta_aod)
        s = summarize_success(out["roundtrip_success"], f"pulse={shift_um}")
        rows.append({
            "factor": "dd_pulse_offset_um", "value": shift_um,
            "shots": s["shots"], "successes": s["successes"],
            "success_fraction": s["probability"], "wilson_low": s["wilson_low"],
            "wilson_high": s["wilson_high"],
            "contrast_none": conditional_contrast(
                phases["none"], out["roundtrip_success"])["contrast"],
            "contrast_xy4": conditional_contrast(
                phases["xy4_per_long_move"],
                out["roundtrip_success"])["contrast"],
            "note": "XY4 全体脉冲整体平移（理想瞬时脉冲模型）"})
        print(f"[sensitivity] pulse offset {shift_um}μm: "
              f"{s['probability']:.4f}")
    _write_rows_csv(output_dir / "sensitivity.csv", rows)
    plot_records.append(viz.plot_sensitivity_panels(
        rows, output_dir / "sensitivity_panels.png"))
    return {"rows": rows, "eta_free": eta_free}


# ------------------------------------------------------------------ 输出
def _write_timeline_csv(cfg, protocols, output_dir):
    """protocol_timeline.csv：时间、分段、深度、中心/速度、lensing、DD 状态。"""
    rows = []
    dt = 0.5 * US
    for name in ("protocol_a", "protocol_b"):
        timeline = protocols[name]
        grid = timeline.build_grid(0.05 * US)
        pulses = _pulses_for(timeline, 0.05 * US)
        times = grid["times"][::10]  # 0.5 μs 采样
        pulse_set = np.sort(np.concatenate([pulses[m] for m in DD_MODES
                                            if pulses[m].size]))
        seg_names = grid["segment_names"]
        for mode in DD_MODES:
            y = toggling_function(grid["times"], pulses[mode])
            for i in range(0, len(grid["times"]), 10):
                near = any(abs(grid["times"][i] - p) < 0.3 * US
                           for p in pulse_set)
                rows.append({
                    "protocol": name, "time_us": grid["times"][i] * 1e6,
                    "segment": seg_names[grid["segment_index"][i]],
                    "d_slm_uK": grid["d_slm"][i] / KB * 1e6,
                    "d_aod_uK": grid["d_aod"][i] / KB * 1e6,
                    "c1_um": grid["c1"][i] * 1e6, "c2_um": grid["c2"][i] * 1e6,
                    "v1_m_s": grid["c1d"][i], "v2_m_s": grid["c2d"][i],
                    "a1_m_s2": grid["c1dd"][i], "a2_m_s2": grid["c2dd"][i],
                    "zs1_um": grid["zs1"][i] * 1e6,
                    "zs2_um": grid["zs2"][i] * 1e6,
                    "dd_mode": mode, "toggling_y": float(y[i]),
                    "pulse_marker": bool(near),
                })
    _write_rows_csv(output_dir / "protocol_timeline.csv", rows)
    return rows


def _ledger_summary_rows(output_dir, protocols, cfg, resolved):
    """重算 ledger 汇总行（均值/最大残差）用于画图与 metrics。"""
    timeline = protocols["protocol_a"]
    states = sample_slm_states(cfg, cfg.development_seed + 500, 64)
    dt = 0.10 * US
    run = run_condition(timeline, dt, states["pos_m"], states["vel_m_per_s"],
                        cfg, "ledger_summary")
    rows, _seg = segment_ledger_rows(run["result"],
                                     timeline.build_grid(dt), "protocol_a",
                                     per_shot=False)
    return segment_ledger_summary(_seg)


def _build_summary(cfg, resolved, preflight, single, repeated, ablations,
                   sensitivity, ledger_summary, output_dir=None) -> str:
    lines = ["# Level 4 端到端往返协议与半经典相干性", "",
             "组合 Level 2 冻结 transfer 波形与 Level 3 冻结三维 AOD 运输，"
             "完成 pick-up—split—长运输—返回—merge—drop-off 往返协议的"
             "经典蒙特卡洛验证与半经典差分光移相位分析。", "",
             f"- Protocol A 总时长 "
             f"{ps.protocol_a_total_us(cfg):.0f} μs；Protocol B 总时长 "
             f"{ps.protocol_b_total_us(cfg, resolved['level2_duration_us']):.0f} μs。",
             f"- nominal v_s = {resolved['nominal_vs_m_per_s']:.4f} m/s"
             f"（Level 3 severity 1.0 冻结场景）。", ""]
    if preflight:
        lines.append(f"- preflight：{len(preflight['critical_checks'])} 项检查"
                     f"{'全部通过' if preflight['all_passed'] else '存在失败'}"
                     f"（{preflight.get('elapsed_s_total', 0):.0f} s）。")
    if single:
        lines += ["", "## 单轮 validation（validation_pool，2000 shots，"
                  "dt=0.05 μs）", "",
                  "| 条件 | roundtrip success | Wilson 95% CI | final retained |"
                  " recaptured |", "|---|---:|---|---:|---:|"]
        for name, s in single["summaries"].items():
            lines.append(
                f"| {name} | {s['successes']}/{s['shots']} "
                f"({s['probability']:.4f}) | "
                f"[{s['wilson_low']:.4f}, {s['wilson_high']:.4f}] | "
                f"{s['final_retained_fraction']:.4f} | {s['recaptured']} |")
        lines += ["", "### 配对比较（同初态）"]
        for name, comp in single["comparisons"].items():
            boot = comp["bootstrap"]
            lines.append(f"- {name}: 差 {boot['difference']:+.4f}，95% CI "
                         f"[{boot['ci_low']:+.4f}, {boot['ci_high']:+.4f}]"
                         f"{'（显著）' if boot['ci_excludes_zero'] else '（含 0）'}。")
        lines += ["", "### 相干对比度（条件于 roundtrip success）"]
        for key, c in single["contrast"].items():
            if c["contrast"] is None:
                continue
            lines.append(f"- {key}: C={c['contrast']:.4f}，"
                         f"mean direction={c['mean_direction']:+.3f} rad，"
                         f"circular std={c['circular_std']:.4f}。")
    if repeated:
        lines += ["", "## 重复往返（500 shots × 10 rounds，dt=0.10 μs）", ""]
        for pname in ("protocol_a", "protocol_b"):
            key = f"{pname}|none"
            surv = repeated["survival"][key]
            lines.append(f"- {pname}: 10 轮 absorbing survival "
                         f"{surv['absorbing'][-1]:.4f}，非吸收 final retention "
                         f"{surv['non_absorbing'][-1]:.4f}。")
            slope = repeated["heating_slopes"].get(pname, {})
            if slope.get("slope_uK_per_round") is not None:
                lines.append(f"- {pname} 加热斜率 "
                             f"{slope['slope_uK_per_round']:+.4f} μK/轮"
                             f"（R²={slope['r_squared']:.3f}）。")
        lt = repeated["long_tail"]
        lines.append(f"- 长尾（128 shots × {lt['max_rounds']} rounds）："
                     f"末轮 absorbing survival {lt['survival'][-1]:.4f}。")
    if ablations:
        lines += ["", "## 消融", ""]
        for row in ablations["rows"]:
            frac = row.get("success_fraction")
            frac_str = f"{frac:.4f}" if frac is not None else "n/a"
            lines.append(f"- {row['condition']}: {frac_str}"
                         f"（{row.get('note', '')}）")
    if sensitivity:
        lines += ["", "## 敏感性（one-factor-at-a-time）", ""]
        for row in sensitivity["rows"]:
            frac = row.get("success_fraction")
            frac_str = f"{frac:.4f}" if frac is not None else "—"
            extra = ""
            if row.get("contrast_none") is not None:
                extra = f"，C(none)={row['contrast_none']:.4f}"
            lines.append(f"- {row['factor']} = {row['value']}: {frac_str}{extra}")
    if ledger_summary:
        lines += ["", "## 功-能账本（Protocol A 代表批）", ""]
        for row in ledger_summary:
            lines.append(
                f"- {row['segment']}: ΔE 均值 {row['mean_delta_e_uK']:+.3f} μK，"
                f"W_ext 均值 {row['mean_w_ext_uK']:+.3f} μK，"
                f"max|R| {row['max_abs_residual_uK']:.2e} μK。")
    lines += ["", "## 诊断图清单", ""]
    for name, desc in viz.FIGURE_CATALOG.items():
        if output_dir is not None:
            exists = (Path(output_dir) / name).is_file()
            mark = "" if exists else "（本次未生成："
            if exists:
                lines.append(f"- `{name}`：{desc}")
            else:
                reason = ("全程无 first-failure 事件，无可绘数据）"
                          if name == "failure_events_by_round.png"
                          else "该 stage 未运行）")
                lines.append(f"- `{name}`：{desc} {mark}{reason}")
        else:
            lines.append(f"- `{name}`：{desc}")
    lines += ["", "## 与论文的关系与边界", "",
              "- 本级输出的是经典协议成功概率（checkpoint basin 归属 + 末态"
              "E_SLM<0）与半经典轨迹相位对比度，**不得**称为论文的 transport "
              "fidelity 或 IRB fidelity：未包含 Clifford twirling、随机基准拟合、"
              "Raman 散射、真空碰撞、脉冲误差与成像误差。",
              "- usable_coherent_fraction = P(roundtrip_success)×C_cond 只是"
              "工程代理，不是量子过程保真度。",
              "- remote 54 μs 仅作为相位积累 hold，不输出任何 gate fidelity。",
              "- 深度 ramp 形状（smootherstep）与 η=1.3e-4 为模型假设；"
              "噪声参数标记 assumed。", "",
              "## Level 5 建议", "",
              "- 有限时长/含脉冲误差的 Rabi 驱动与完整自旋哈密顿量；",
              "- 真实噪声谱（磁场、强度 OU 谱）与 pulse-area error；",
              "- Clifford/IRB 随机基准协议拟合；",
              "- 多原子并行操作的相互作用与成像误差。", "",
              "图片完成程序化检查；未进行人工视觉审阅。", ""]
    return "\n".join(lines)


def _apply_visual_review(output_dir: Path, findings_path: Path) -> int:
    """把视觉审阅结论合并进既有产物，不重跑任何阶段。"""
    findings = json.loads(Path(findings_path).read_text(encoding="utf-8"))
    per_image = findings.get("images", {})
    issues = {k: v.get("issues", "") for k, v in per_image.items()
              if not v.get("layout_ok", True)}
    patch = {
        "visual_review_performed": True,
        "visual_review_by": findings.get("reviewer", "unspecified"),
        "visual_review_findings": per_image,
        "visual_review_all_ok": not issues,
    }
    for name in ("image_validation.json", "metrics.json"):
        path = output_dir / name
        if not path.is_file():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        data.update(patch)
        if name == "metrics.json":
            data["visual_review_summary"] = {
                "reviewed_images": len(per_image), "issues": issues,
                "note": findings.get("note", "")}
        save_json(path, data)
    summary_path = output_dir / "summary.md"
    if summary_path.is_file():
        text = summary_path.read_text(encoding="utf-8")
        verdict = ("全部通过" if not issues else
                   "发现 %d 处问题：%s" % (len(issues), "; ".join(issues.values())))
        text = text.replace(
            "图片完成程序化检查；未进行人工视觉审阅。",
            f"图片完成程序化检查，并经逐图视觉审阅（{len(per_image)} 张，"
            f"{findings.get('reviewer', 'unspecified')}）：{verdict}。")
        summary_path.write_text(text, encoding="utf-8")
    print(f"[visual_review] merged into {output_dir}（issues={list(issues)}）")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    cfg = load_config(Path(args.config))
    output_dir = ensure_output_directory(
        Path(args.output_dir) if args.output_dir else
        cfg.project_root / cfg.output_directory, cfg.project_root)
    if args.apply_visual_review:
        return _apply_visual_review(output_dir, Path(args.apply_visual_review))
    manifest, resolved = upstream_artifacts.build_manifest(cfg)
    protocols = build_protocols(cfg, resolved)
    save_yaml(output_dir / "config_used.yaml", {
        **config_as_dict(cfg),
        "resolved_upstream": {
            "nominal_vs_m_per_s": resolved["nominal_vs_m_per_s"],
            "level2_waveform": manifest["level2_waveform"],
            "level3_lensing": manifest["level3_lensing"],
        }})
    save_json(output_dir / "upstream_manifest.json", manifest)
    stage = args.stage
    results = {"preflight": None, "single_round": None, "repeated_rounds": None,
               "ablations": None, "sensitivity": None}
    plot_records = []
    if stage in ("preflight", "all"):
        results["preflight"] = stage_preflight(cfg, resolved, protocols,
                                               output_dir)
    timeline_rows = _write_timeline_csv(cfg, protocols, output_dir)
    checkpoints_map = {
        "protocol_a": protocols["protocol_a"].build_grid(
            0.05 * US)["checkpoints"],
        "protocol_b": protocols["protocol_b"].build_grid(
            0.05 * US)["checkpoints"],
    }
    plot_records.append(viz.plot_upstream_provenance(
        manifest, output_dir / "upstream_provenance.png"))
    init_states = sample_slm_states(cfg, cfg.single_round_seed,
                                    cfg.single_round_validation_shots)
    plot_records.append(viz.plot_initial_ensemble(
        init_states, output_dir / "initial_ensemble.png"))
    if results["preflight"]:
        plot_records.append(viz.plot_preflight_checks(
            results["preflight"], output_dir / "preflight_checks.png"))
    plot_records.append(viz.plot_protocol_timeline(
        {"protocol_a": protocols["protocol_a"],
         "protocol_b": protocols["protocol_b"]},
        checkpoints_map, output_dir / "protocol_timeline.png",
        dd_pulse_sets={
            "protocol_a": _pulses_for(protocols["protocol_a"], 0.05 * US),
            "protocol_b": _pulses_for(protocols["protocol_b"], 0.05 * US)}))
    grid_a = protocols["protocol_a"].build_grid(0.05 * US)
    grid_b = protocols["protocol_b"].build_grid(0.05 * US)
    plot_records.append(viz.plot_trap_path_and_depths(
        grid_a, output_dir / "trap_path_and_depths.png", grid_b=grid_b))
    if stage in ("single_round", "all", "repeated_rounds", "ablations",
                 "sensitivity"):
        frozen = _frozen_validation_yaml(cfg, resolved)
        frozen_path = output_dir / "frozen_validation_protocols.yaml"
        if not frozen_path.is_file():
            save_yaml(frozen_path, frozen)
    if stage in ("single_round", "all", "repeated_rounds", "ablations",
                 "sensitivity"):
        results["single_round"] = stage_single_round(
            cfg, resolved, protocols, output_dir, plot_records)
    if stage in ("repeated_rounds", "all"):
        results["repeated_rounds"] = stage_repeated_rounds(
            cfg, resolved, protocols, output_dir, plot_records)
    if stage in ("ablations", "all"):
        results["ablations"] = stage_ablations(
            cfg, resolved, protocols, output_dir, plot_records,
            single_round=results["single_round"])
    if stage in ("sensitivity", "all"):
        results["sensitivity"] = stage_sensitivity(
            cfg, resolved, protocols, output_dir, plot_records,
            single_round=results["single_round"])
    # 通用图
    ledger_summary = _ledger_summary_rows(output_dir, protocols, cfg, resolved)
    plot_records.append(viz.plot_energy_ledger(
        ledger_summary, output_dir / "energy_ledger.png"))
    phi_data = {}
    contrast_data = {}
    if results["single_round"]:
        integ = results["single_round"]["phase_integrals"][
            "protocol_a_nominal_lensing"]
        for mode in DD_MODES:
            phi_data[mode] = ((cfg.eta_slm * integ[f"a_slm_{mode}"]
                               + cfg.eta_aod * integ[f"a_aod_{mode}"]) / HBAR)
    if results["repeated_rounds"]:
        contrast_data = results["repeated_rounds"]["contrast_curve"]
    if phi_data or contrast_data:
        plot_records.append(viz.plot_phase_and_contrast(
            phi_data, contrast_data, output_dir / "phase_and_contrast.png"))
    if results["preflight"] and results["preflight"]["dt_convergence"].get(
            "rows"):
        plot_records.append(viz.plot_timestep_convergence(
            results["preflight"]["dt_convergence"]["rows"],
            output_dir / "timestep_convergence.png"))
    image_validation = validate_pngs([r for r in plot_records if r])
    save_json(output_dir / "image_validation.json", image_validation)
    metrics = {
        "level4": {"name": cfg.name, "stage": stage},
        "upstream_manifest": {"files": manifest["files"],
                              "level3_lensing": manifest["level3_lensing"],
                              "level2_waveform": manifest["level2_waveform"],
                              "reconstructed_from_config":
                                  manifest["reconstructed_from_config"]},
        "preflight": results["preflight"],
        "frozen_validation_protocols": _frozen_validation_yaml(cfg, resolved),
        "single_round": None if not results["single_round"] else {
            "summaries": results["single_round"]["summaries"],
            "comparisons": results["single_round"]["comparisons"],
            "contrast": results["single_round"]["contrast"],
            "checkpoint_stats": results["single_round"]["checkpoint_stats"],
            "first_failure_distribution":
                results["single_round"]["stage_counts"],
            "initial_states_meta":
                results["single_round"]["initial_states_meta"],
            "stress_pool": results["single_round"]["stress"],
        },
        "repeated_rounds": None if not results["repeated_rounds"] else {
            "survival": results["repeated_rounds"]["survival"],
            "heating": results["repeated_rounds"]["heating"],
            "heating_slopes": results["repeated_rounds"]["heating_slopes"],
            "survival_fits": results["repeated_rounds"]["fits"],
            "long_tail": results["repeated_rounds"]["long_tail"],
        },
        "ablations": None if not results["ablations"] else {
            "rows": results["ablations"]["rows"],
            "note": "消融用于归因，未据此重选主协议"},
        "sensitivity": None if not results["sensitivity"] else {
            "rows": results["sensitivity"]["rows"],
            "eta_linear_recalc": results["sensitivity"]["eta_free"],
        },
        "work_energy": {"ledger_summary": ledger_summary},
        "comparability": {
            "comparable_to_paper": [
                "分段时长锚点（Extended Data Fig. 10e）、2.4 μm 分离、"
                "375 μm 对角运输距离、深度量级"],
            "not_comparable": [
                "经典成功概率 ≠ 论文 IRB/transport fidelity",
                "C_cond 与 usable_coherent_fraction ≠ 量子过程保真度",
                "lensing severity 为 Level 3 无量纲场景，非论文测量值",
                "η、噪声参数为假设值"],
        },
        "image_validation_all_passed": image_validation["all_passed"],
        "visual_review_performed": False,
    }
    save_json(output_dir / "metrics.json", metrics)
    save_markdown(output_dir / "summary.md",
                  _build_summary(cfg, resolved, results["preflight"],
                                 results["single_round"],
                                 results["repeated_rounds"],
                                 results["ablations"],
                                 results["sensitivity"], ledger_summary,
                                 output_dir=output_dir))
    print(f"Level 4 stage '{stage}' completed: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
