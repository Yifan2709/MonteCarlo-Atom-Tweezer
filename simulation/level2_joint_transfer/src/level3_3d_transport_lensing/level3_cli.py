"""Level 3 命令行入口：preflight / potential_scan / coarse_mc / validate / sensitivity / all。"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import yaml

from level0_static_trap.io_utils import (ensure_output_directory, save_json,
                                         save_markdown, save_yaml, validate_pngs)
from . import level3_visualization as viz
from .gaussian_3d import KB
from .integrators_3d import run_trajectory
from .level3_config import Level3Config, config_as_dict, load_config
from .level3_simulation import (compare_paired, deterministic_scan,
                                make_control, physics_from_config, resolve_vs,
                                run_condition, sample_states)
from .preflight3 import run_preflight
from .transport_statistics import summarize_retention

FROZEN_CONDITIONS = [
    {"name": "paper_anchor_diagonal", "geometry": "diagonal",
     "profile": "adiabatic_sine", "distance_um": 610.0, "duration_us": 1600.0},
    {"name": "paper_anchor_constant_jerk", "geometry": "diagonal",
     "profile": "constant_jerk", "distance_um": 610.0, "duration_us": 1600.0},
    {"name": "straight_control", "geometry": "straight",
     "profile": "adiabatic_sine", "distance_um": 510.0, "duration_us": 1600.0},
    {"name": "straight_control_constant_jerk", "geometry": "straight",
     "profile": "constant_jerk", "distance_um": 510.0, "duration_us": 1600.0},
]
VALIDATION_SEVERITIES = (0.0, 1.0)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Level 3 三维 AOD 长距离运输与柱面透镜")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("preflight", "potential_scan", "coarse_mc",
                                            "validate", "sensitivity", "all"),
                        default="all")
    parser.add_argument("--output-dir", help="覆盖输出目录")
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


def stage_preflight(cfg: Level3Config, output_dir: Path) -> dict:
    """16 项预检查，失败抛异常。"""
    started = time.time()
    checks = run_preflight(cfg)
    checks["elapsed_s"] = time.time() - started
    save_json(output_dir / "preflight.json", checks)
    if not checks["all_passed"]:
        failed = [k for k in checks["critical_checks"]
                  if not (checks[k].get("passed", False)
                          if isinstance(checks[k], dict) else False)]
        raise RuntimeError(f"preflight 未通过，停止大规模扫描：{failed}")
    print(f"[preflight] {len(checks['critical_checks'])} 项检查全部通过"
          f"（{checks['elapsed_s']:.0f} s）")
    return checks


def stage_potential_scan(cfg: Level3Config, output_dir: Path) -> dict:
    """阶段 A：确定性势阱与轨迹扫描。"""
    physics = physics_from_config(cfg)
    minima_rows, traj_rows = deterministic_scan(cfg, physics)
    _write_rows_csv(output_dir / "potential_minima_scan.csv", minima_rows)
    _write_rows_csv(output_dir / "trajectory_profiles.csv", traj_rows)
    print(f"[potential_scan] {len(minima_rows)} 条势阱诊断 / {len(traj_rows)} 条轨迹点")
    return {"minima_rows": minima_rows, "traj_rows": traj_rows, "physics": physics}


def stage_coarse_mc(cfg: Level3Config, output_dir: Path) -> dict:
    """阶段 B：粗蒙特卡洛相图（128 shots, dt=0.10 μs, 公共 scan_pool）。"""
    physics = physics_from_config(cfg)
    states = sample_states(cfg, physics, cfg.scan_seed, cfg.scan_shots)
    rows = []
    grid = []
    for geometry, distance in (("straight", 510.0), ("diagonal", 610.0)):
        for profile in ("adiabatic_sine", "constant_jerk"):
            for duration in (600.0, 800.0, 1000.0, 1200.0, 1600.0, 2000.0):
                for severity in (0.0, 0.75, 1.00, 1.25):
                    grid.append((geometry, distance, profile, duration, severity))
    for index, (geometry, distance, profile, duration, severity) in enumerate(grid):
        v_s = resolve_vs(cfg, severity)
        started = time.time()
        try:
            out = run_condition(physics, states, profile, geometry, distance,
                                duration, severity, v_s, cfg.axis_signs,
                                cfg.scan_dt_us, cfg.post_transport_hold_us, cfg)
            summary = out["summary"]
            rows.append({
                "geometry": geometry, "distance_um": distance, "profile": profile,
                "duration_us": duration, "severity": severity,
                "shots": summary["shots"], "retained": summary["retained"],
                "retention_fraction": summary["retention_fraction"],
                "wilson_low": summary["wilson_low"],
                "wilson_high": summary["wilson_high"],
                "near_threshold": sum(r["near_threshold"] for r in out["records"]),
                "median_delta_eps_uK_retained":
                    summary["median_delta_eps_uK_retained"],
                "median_max_axial_um": float(np.median(
                    [r["max_axial_um"] for r in out["records"]])),
                "split_fraction": summary["split_fraction"],
                "max_abs_w_residual_over_depth":
                    summary["max_abs_w_residual_over_depth"],
                "status": "ok",
            })
        except Exception as exc:
            rows.append({"geometry": geometry, "distance_um": distance,
                         "profile": profile, "duration_us": duration,
                         "severity": severity, "status": f"failed:{exc!r}"[:200]})
        print(f"[coarse_mc] {index + 1}/{len(grid)} {geometry} {profile} "
              f"{duration:.0f}us sev{severity} -> "
              f"{rows[-1].get('retention_fraction', float('nan')):.3f} "
              f"({time.time() - started:.1f}s)")
    _write_rows_csv(output_dir / "coarse_scan.csv", rows)
    return {"rows": rows}


def stage_validate(cfg: Level3Config, output_dir: Path, plot_records: list) -> dict:
    """阶段 C：冻结主条件的 2000-shot 独立验证（一次性 validation_pool）。"""
    physics = physics_from_config(cfg)
    states = sample_states(cfg, physics, cfg.validation_seed, cfg.validation_shots)
    all_records = []
    summaries = {}
    keep_ids = tuple(range(5))
    conditions = []
    for cond in FROZEN_CONDITIONS:
        for severity in VALIDATION_SEVERITIES:
            conditions.append(dict(cond, severity=severity))
    for cond in conditions:
        v_s = resolve_vs(cfg, cond["severity"])
        started = time.time()
        out = run_condition(physics, states, cond["profile"], cond["geometry"],
                            cond["distance_um"], cond["duration_us"],
                            cond["severity"], v_s, cfg.axis_signs,
                            cfg.validation_dt_us, cfg.post_transport_hold_us, cfg,
                            keep_trajectory_ids=keep_ids
                            if cond["name"] == "paper_anchor_diagonal" else ())
        key = f"{cond['name']}|sev{cond['severity']}"
        summaries[key] = out["summary"]
        for r in out["records"]:
            r["condition"] = key
        all_records.extend(out["records"])
        if cond["name"] == "paper_anchor_diagonal":
            trajectories = out.get("trajectories")
        print(f"[validate] {key}: {out['summary']['retained']}/"
              f"{out['summary']['shots']} ({time.time() - started:.0f}s)")
    _write_rows_csv(output_dir / "validation_shots.csv", all_records)

    # 0.025 μs 子集复算（convergence_pool 种子）
    conv_states = sample_states(cfg, physics, cfg.convergence_seed, 256)
    convergence_rows = []
    base_labels = None
    for dt_us in (0.10, 0.05, 0.025):
        out = run_condition(physics, conv_states, "adiabatic_sine", "diagonal", 610.0,
                            1600.0, 1.0, resolve_vs(cfg, 1.0), cfg.axis_signs, dt_us,
                            cfg.post_transport_hold_us, cfg)
        labels = np.array([r["retained"] for r in out["records"]])
        energies = np.array([r["e_final_uK"] for r in out["records"]])
        if base_labels is None:
            base_labels = labels
            base_energy = energies
        convergence_rows.append({
            "dt_us": dt_us,
            "max_final_energy_diff_uK": float(np.max(np.abs(energies - base_energy))),
            "max_abs_w_residual_over_depth":
                out["summary"]["max_abs_w_residual_over_depth"],
            "label_flip_rate": float(np.mean(labels != base_labels)),
            "retention_fraction": out["summary"]["retention_fraction"],
        })
    # 代表轨迹（成功 / 丢失 / 双阱经历）保存
    rep_rows = []
    anchor = [r for r in all_records
              if r["condition"] == "paper_anchor_diagonal|sev1.0"]
    for label, pred in (("retained", lambda r: r["retained"]),
                        ("lost", lambda r: not r["retained"]),
                        ("split_exposed", lambda r: r["split_exposed"])):
        picks = [r for r in anchor if pred(r)][:cfg.sensitivity_shots // 100 + 1]
        for r in picks:
            rep_rows.append({"class": label, "shot_id": r["shot_id"],
                             "branch": r["branch"], "mean_z_um": r["mean_z_um"],
                             "e_final_uK": r["e_final_uK"]})
    _write_rows_csv(output_dir / "representative_trajectories.csv", rep_rows)

    # 配对比较
    comparisons = {}
    by_cond = {}
    for r in all_records:
        by_cond.setdefault(r["condition"], []).append(r)
    pairs = [
        ("sine_vs_constant_jerk_diagonal_sev1.0",
         "paper_anchor_diagonal|sev1.0", "paper_anchor_constant_jerk|sev1.0"),
        ("sine_vs_constant_jerk_straight_sev1.0",
         "straight_control|sev1.0", "straight_control_constant_jerk|sev1.0"),
        ("straight_vs_diagonal_sine_sev1.0",
         "straight_control|sev1.0", "paper_anchor_diagonal|sev1.0"),
        ("lensing_off_vs_sev1.0_paper_anchor",
         "paper_anchor_diagonal|sev0.0", "paper_anchor_diagonal|sev1.0"),
    ]
    for name, key_a, key_b in pairs:
        comparisons[name] = compare_paired(by_cond[key_a], by_cond[key_b],
                                           cfg.validation_seed)

    # 功-能代表曲线：从 validation 记录重建 shot 0 的 ΔE / W_ext / R_W 汇总
    work_cases = {}
    for label, key in (("ideal_diagonal", "paper_anchor_diagonal|sev0.0"),
                       ("lensing_diagonal", "paper_anchor_diagonal|sev1.0"),
                       ("lensing_straight", "straight_control|sev1.0")):
        recs = by_cond.get(key, [])
        if recs:
            work_cases[label] = {
                "time_us": np.array([0.0, 1600.0]),
                "delta_E": np.array([0.0, recs[0]["e_final_uK"] - recs[0]["e0_uK"]]),
                "W_ext": np.array([0.0, recs[0]["w_total_uK"]]),
                "R_W": np.array([0.0, recs[0]["w_residual_uK"]]) * 100,
            }
    plot_records.append(viz.plot_heating(all_records, output_dir / "heating_comparison.png"))
    plot_records.append(viz.plot_final_energy(all_records,
                                              output_dir / "final_energy_distributions.png"))
    if convergence_rows:
        plot_records.append(viz.plot_timestep_convergence(
            convergence_rows, output_dir / "timestep_convergence.png"))
    if work_cases and any(w.get("delta_E") is not None for w in work_cases.values()):
        usable = {k: v for k, v in work_cases.items() if v.get("delta_E") is not None}
        plot_records.append(viz.plot_work_energy(usable,
                                                 output_dir / "work_energy_balance.png"))
    axial_trajs = {}
    if trajectories:
        pos = trajectories["position_m"]  # (n_rec, n_ids, 3)
        anchor_records = by_cond.get("paper_anchor_diagonal|sev1.0", [])
        for i in range(pos.shape[1]):
            rec = anchor_records[i] if i < len(anchor_records) else None
            label = f"shot {rec['shot_id'] if rec else i}" \
                    f"({'retained' if rec and rec['retained'] else 'lost'})"
            axial_trajs[label] = {
                "time_s": trajectories["time_s"],
                "position_m": pos[:, i:i + 1, :],
                "velocity_m_per_s": trajectories["velocity_m_per_s"][:, i:i + 1, :],
            }
    if axial_trajs:
        plot_records.append(viz.plot_axial_dynamics(
            axial_trajs, output_dir / "axial_dynamics.png"))
    return {"summaries": summaries, "comparisons": comparisons,
            "convergence_rows": convergence_rows,
            "validation_states_meta": {k: v for k, v in states.items()
                                       if not isinstance(v, np.ndarray)},
            "conditions": [dict(c) for c in conditions]}


def stage_sensitivity(cfg: Level3Config, output_dir: Path, plot_records: list) -> dict:
    """阶段 D：温度/深度/束腰/severity/符号敏感性（冻结模型后 one-factor-at-a-time）。"""
    base = {"profile": "adiabatic_sine", "geometry": "diagonal", "distance_um": 610.0,
            "duration_us": 1600.0, "severity": 1.0}
    rows = []
    scans = (
        ("temperature", [3.0, 5.0, 10.0], "temperature_uK"),
        ("depth", [160.0, 280.0, 920.0], "depth_uK"),
        ("waist", [1.11, 1.17, 1.23], "waist_um"),
        ("severity", [0.50, 0.75, 1.00, 1.25], "severity"),
    )
    jobs = []
    for scan, values, field in scans:
        for value in values:
            jobs.append((scan, value, field))
    jobs.append(("axis_signs_opposite", 1.0, "opposite_signs"))
    for index, (scan, value, field) in enumerate(jobs):
        physics = physics_from_config(
            cfg, depth_uK=value if field == "depth_uK" else None,
            waist_um=value if field == "waist_um" else None)
        states = sample_states(
            cfg, physics, cfg.scan_seed + 7000 + index, cfg.sensitivity_shots,
            temperature_uK=value if field == "temperature_uK" else None)
        severity = value if field == "severity" else base["severity"]
        signs = (-1, 1) if field == "opposite_signs" else cfg.axis_signs
        v_s = resolve_vs(cfg, severity)
        out = run_condition(physics, states, base["profile"], base["geometry"],
                            base["distance_um"], base["duration_us"], severity, v_s,
                            signs, cfg.validation_dt_us, cfg.post_transport_hold_us, cfg)
        s = out["summary"]
        rows.append({"scan": scan, "value": value, "shots": s["shots"],
                     "retained": s["retained"],
                     "retention_fraction": s["retention_fraction"],
                     "wilson_low": s["wilson_low"], "wilson_high": s["wilson_high"],
                     "median_delta_eps_uK_retained": s["median_delta_eps_uK_retained"],
                     "split_fraction": s["split_fraction"]})
        print(f"[sensitivity] {scan}={value}: {s['retained']}/{s['shots']}")
    _write_rows_csv(output_dir / "sensitivity.csv", rows)
    plot_records.append(viz.plot_sensitivity(rows, output_dir / "sensitivity_panels.png"))
    return {"rows": rows}


def _build_summary(cfg, preflight, pot_scan, coarse, validation, sensitivity) -> str:
    lines = ["# Level 3 三维 AOD 长距离运输与柱面透镜效应", "",
             "隔离研究 AOD 单阱长距离搬运：三维经典模型 + 论文 SI 启发的 crossed-AOD "
             "柱面透镜势；留阱率仅指经典末态仍被束缚，不代表量子保真度。", ""]
    physics = physics_from_config(cfg)
    lines.append(f"- z_R = {physics['z_r_m']*1e6:.3f} μm，"
                 f"ω_r/2π = {physics['omega_r']/2/np.pi/1e3:.2f} kHz，"
                 f"ω_z/2π = {physics['omega_z']/2/np.pi/1e3:.2f} kHz。")
    if preflight:
        lines.append(f"- preflight：{len(preflight['critical_checks'])} 项检查"
                     f"{'全部通过' if preflight['all_passed'] else '存在失败'}"
                     f"（{preflight.get('elapsed_s', 0):.0f} s）。")
    if validation:
        lines += ["", "## 冻结主条件独立验证（validation_pool，一次性）", "",
                  "| 条件 | 留阱 | 留阱率 | Wilson 95% CI |",
                  "|---|---:|---:|---|"]
        for key, s in validation["summaries"].items():
            lines.append(f"| {key} | {s['retained']}/{s['shots']} | "
                         f"{s['retention_fraction']:.4f} | "
                         f"[{s['wilson_low']:.4f}, {s['wilson_high']:.4f}] |")
        lines += ["", "### 配对比较"]
        for name, comp in validation["comparisons"].items():
            boot = comp["bootstrap"]
            lines.append(f"- {name}: 差 {boot['difference']:+.4f}，"
                         f"95% CI [{boot['ci_low']:+.4f}, {boot['ci_high']:+.4f}]"
                         f"{'（显著）' if boot['ci_excludes_zero'] else '（含 0，不显著）'}。")
    if sensitivity:
        lines += ["", "## 敏感性（one-factor-at-a-time，每点 "
                 f"{cfg.sensitivity_shots} shots）", ""]
        for row in sensitivity["rows"]:
            lines.append(f"- {row['scan']} = {row['value']}: "
                         f"{row['retention_fraction']:.4f}"
                         f"（Wilson [{row['wilson_low']:.4f}, {row['wilson_high']:.4f}]）")
    lines += ["", "## 边界", "",
              "经典留阱率不可等同于论文的 99.953(2)% 量子通道保真度；lensing severity "
              "为无量纲扫描场景而非论文测量值。未包含：真空损失、光子散射、退相干、"
              "pick-up/drop-off 组合、重复往返。", "",
              "图片完成程序化检查；未进行人工视觉审阅。", ""]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    cfg = load_config(Path(args.config))
    output_dir = ensure_output_directory(
        Path(args.output_dir) if args.output_dir else
        cfg.project_root / cfg.output_directory, cfg.project_root)
    save_yaml(output_dir / "config_used.yaml", config_as_dict(cfg))

    plot_records = []
    stage = args.stage
    results = {"preflight": None, "potential_scan": None, "coarse_mc": None,
               "validation": None, "sensitivity": None}
    physics = physics_from_config(cfg)
    # 通用诊断图（无 MC 依赖）
    plot_records.append(viz.plot_static_slices(physics,
                                               output_dir / "static_3d_trap_slices.png"))
    if stage in ("preflight", "all"):
        results["preflight"] = stage_preflight(cfg, output_dir)
    if stage in ("potential_scan", "all"):
        results["potential_scan"] = stage_potential_scan(cfg, output_dir)
        if results["potential_scan"]["minima_rows"]:
            plot_records.append(viz.plot_bifurcation(
                physics, output_dir / "lensing_potential_bifurcation.png"))
            plot_records.append(viz.plot_profiles(
                output_dir / "trajectory_profile_comparison.png"))
            plot_records.append(viz.plot_straight_vs_diagonal(
                results["potential_scan"]["minima_rows"], physics,
                output_dir / "straight_vs_diagonal.png"))
    if stage in ("coarse_mc", "all"):
        results["coarse_mc"] = stage_coarse_mc(cfg, output_dir)
        plot_records.append(viz.plot_survival_phase(
            results["coarse_mc"]["rows"], output_dir / "survival_phase_diagram.png"))
    if stage in ("validate", "all"):
        results["validation"] = stage_validate(cfg, output_dir, plot_records)
    if stage in ("sensitivity", "all"):
        results["sensitivity"] = stage_sensitivity(cfg, output_dir, plot_records)

    image_validation = validate_pngs([r for r in plot_records if r])
    save_json(output_dir / "image_validation.json", image_validation)
    metrics = {
        "level3": {"name": cfg.name, "stage": stage},
        "physics": {"z_r_m": physics["z_r_m"], "omega_r": physics["omega_r"],
                    "omega_z": physics["omega_z"], "depth_uK": cfg.depth_uK,
                    "waist_um": cfg.waist_um,
                    "waist_source": cfg.waist_source},
        "seeds": {"scan": cfg.scan_seed, "validation": cfg.validation_seed,
                  "convergence": cfg.convergence_seed},
        "lensing": {
            "calibration_mode": cfg.calibration_mode,
            "fixed_vs_m_per_s": {
                str(sev): resolve_vs(cfg, sev)
                for sev in cfg.reference_axis_vmax_over_vs if sev > 0},
            "scenario_note": "normalized_reference_scan 的 v_s 为无量纲场景参数，"
                             "非论文测量值"},
        "frozen_validation_conditions": [dict(c) for c in FROZEN_CONDITIONS],
        "preflight": results["preflight"],
        "coarse_scan": None if not results["coarse_mc"] else {
            "n_grid": len(results["coarse_mc"]["rows"]),
            "failed": sum(1 for r in results["coarse_mc"]["rows"]
                          if r.get("status", "ok") != "ok")},
        "validation": None if not results["validation"] else {
            "summaries": results["validation"]["summaries"],
            "comparisons": results["validation"]["comparisons"],
            "convergence_rows": results["validation"]["convergence_rows"],
            "initial_state_meta": results["validation"]["validation_states_meta"]},
        "sensitivity": None if not results["sensitivity"] else {
            "rows": results["sensitivity"]["rows"]},
        "image_validation_all_passed": image_validation["all_passed"],
        "visual_review_performed": False,
    }
    save_json(output_dir / "metrics.json", metrics)
    save_markdown(output_dir / "summary.md",
                  _build_summary(cfg, results["preflight"], results["potential_scan"],
                                 results["coarse_mc"], results["validation"],
                                 results["sensitivity"]))
    print(f"Level 3 stage '{stage}' completed: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
