"""Level 2 命令行入口：preflight / optimize / validate / robustness / all。"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import yaml

from level0_static_trap.io_utils import (ensure_output_directory, save_csv, save_json,
                                         save_markdown, save_yaml, validate_pngs)
from . import level2_visualization as viz
from .config import Level2Config, config_as_dict, load_config, with_output_override
from .level2_simulation import (baseline_specs, build_waveform, evaluate_waveform_on_states,
                                ideal_run, physics_from_config, run_validation,
                                sample_initial_states)
from .noise_heating import heating_from_config
from .optimizer import (frozen_spec, restore_history, run_spline_refinement, run_stage1,
                        run_stage2, select_and_freeze, write_candidates_csv)
from .preflight import run_preflight
from .robustness import alignment_scan, temperature_scan, timestep_scan
from .statistics import summarize_capture


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Level 2 同步波形优化与鲁棒验证")
    parser.add_argument("--config", required=True, help="YAML 配置路径")
    parser.add_argument("--stage", choices=("preflight", "optimize", "validate", "robustness", "all"),
                        default="all")
    parser.add_argument("--output-dir", help="覆盖输出目录；相对路径按项目根解析")
    return parser


def _frozen_or_none(output_dir: Path):
    path = output_dir / "best_waveform.yaml"
    if not path.exists():
        return None
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _main_specs(cfg: Level2Config, frozen) -> dict:
    specs = dict(baseline_specs(cfg))
    specs["overlapped_optimized_400us"] = frozen_spec(frozen, cfg)
    return specs


def stage_preflight(cfg: Level2Config, output_dir: Path) -> dict:
    """运行十项预检查并保存 preflight.json；失败时抛异常。"""
    started = time.time()
    checks = run_preflight(cfg)
    checks["elapsed_s"] = time.time() - started
    save_json(output_dir / "preflight.json", checks)
    if not checks["all_passed"]:
        def _passed(name):
            value = checks[name]
            return bool(value) if isinstance(value, bool) else bool(value.get("passed", False))
        failed = [name for name in checks["critical_checks"] if not _passed(name)]
        raise RuntimeError(f"preflight 未通过，禁止启动优化：{failed}")
    print(f"[preflight] all {len(checks['critical_checks'])} critical checks passed")
    return checks


def stage_optimize(cfg: Level2Config, output_dir: Path, plot_records: list) -> dict:
    """三阶段优化并冻结最终波形。"""
    physics = physics_from_config(cfg)
    ens = cfg.initial_ensemble
    opt_pool = sample_initial_states(cfg, ens.optimization_seed, ens.optimization_shots)
    sel_pool = sample_initial_states(cfg, ens.selection_seed, ens.selection_shots)
    checkpoint = output_dir / "optimization_checkpoint.json"
    # 续跑时先恢复检查点里的评估历史；阶段缓存命中不再追加，靠行数据兜底补齐
    history: list[dict] = restore_history(checkpoint)

    stage1_rows = run_stage1(cfg, opt_pool, physics, output_dir, checkpoint, history)
    known_ids = {entry.get("candidate_id") for entry in history}
    missing_stage1 = [{"candidate_id": r["candidate_id"], "status": r["status"],
                       "objective_total": r.get("objective_total"),
                       "captured_count": r.get("captured_count")}
                      for r in stage1_rows if r["candidate_id"] not in known_ids]
    # 阶段 1 按候选 id 顺序前置，保持整体评估先后次序
    history[:0] = missing_stage1
    stage2_rows = run_stage2(cfg, stage1_rows, sel_pool, physics, output_dir, checkpoint, history)
    best_stage2 = [r for r in stage2_rows if r.get("status") == "ok"][0]
    spline_rows = run_spline_refinement(cfg, best_stage2, sel_pool, physics, checkpoint, history)
    frozen = select_and_freeze(cfg, stage1_rows, stage2_rows, spline_rows, output_dir)
    known_ids = {entry.get("candidate_id") for entry in history}
    for rows in (stage2_rows, spline_rows):
        history.extend([{"candidate_id": r["candidate_id"], "status": r["status"],
                         "objective_total": r.get("objective_total"),
                         "captured_count": r.get("captured_count")}
                        for r in rows if r["candidate_id"] not in known_ids])

    all_rows = stage1_rows + stage2_rows + spline_rows
    if cfg.output.save_candidate_table:
        write_candidates_csv(output_dir / "candidates.csv", all_rows)
    with (output_dir / "optimization_history.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "status", "objective_total",
                                                    "captured_count"])
        writer.writeheader()
        for row in history:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in
                             ["candidate_id", "status", "objective_total", "captured_count"]})
    if history:
        plot_records.append(viz.plot_optimization_history(history,
                                                          output_dir / "optimization_history.png"))
    valid = sum(1 for r in all_rows if r.get("status") == "ok")
    invalid = sum(1 for r in all_rows if r.get("status") == "invalid")
    failed = sum(1 for r in all_rows if r.get("status") == "failed")
    print(f"[optimize] candidates total={len(all_rows)} ok={valid} invalid={invalid} failed={failed}")
    print(f"[optimize] frozen -> {frozen['type']} provenance={frozen['provenance']}")
    return {"rows": all_rows, "frozen": frozen,
            "counts": {"total": len(all_rows), "ok": valid, "invalid": invalid, "failed": failed}}


def _write_waveforms_csv(cfg: Level2Config, specs: dict, physics: dict, output_dir: Path):
    """三个主波形在公共时间轴上的密集波形表及导数。"""
    max_duration = max(build_waveform(dict(spec, name=n), physics).duration_s
                       for n, spec in specs.items())
    times = np.linspace(0.0, max_duration, 3001)
    columns = {
        "time_s": times, "time_us": times * 1e6,
    }
    for name, spec in specs.items():
        waveform = build_waveform(dict(spec, name=name), physics)
        table = waveform.dense_table(len(times))
        prefix = name
        columns[f"{prefix}_aod_center_m"] = table["aod_center_m"]
        columns[f"{prefix}_aod_center_um"] = table["aod_center_m"] * 1e6
        columns[f"{prefix}_aod_depth_J"] = table["aod_depth_J"]
        columns[f"{prefix}_aod_depth_uK"] = table["aod_depth_J"] / 1.380649e-29 / 1e6
        columns[f"{prefix}_center_velocity_m_per_s"] = table["center_velocity_m_per_s"]
        columns[f"{prefix}_center_acceleration_m_per_s2"] = table["center_acceleration_m_per_s2"]
        columns[f"{prefix}_depth_rate_J_per_s"] = table["depth_rate_J_per_s"]
    save_csv(output_dir / "waveforms.csv", columns)


def _validation_records_rows(results: dict) -> list[dict]:
    rows = []
    for name, result in results.items():
        for record in result["records"]:
            row = {"waveform": name}
            row.update(record)
            row.pop("end_time_index", None)
            rows.append(row)
    return rows


def _write_validation_csv(output_dir: Path, rows: list[dict]):
    if not rows:
        return
    fieldnames = ["waveform", "shot_id", "x0_um", "v0_m_s", "initial_aod_energy_uK",
                  "initial_excitation_uK", "final_x_um", "final_v_m_s", "final_slm_energy_uK",
                  "captured", "near_threshold", "capture_margin", "final_excitation_over_depth",
                  "excitation_change_uK", "final_move_work_uK", "final_depth_work_uK",
                  "final_total_work_uK", "max_abs_work_energy_residual_over_depth",
                  "noise_recoil_energy_uK", "noise_parametric_energy_uK",
                  "noise_pointing_energy_uK", "noise_total_energy_uK",
                  "hold_peak_to_peak_energy_error_over_depth"]
    with (output_dir / "validation_shots.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})


def stage_validate(cfg: Level2Config, output_dir: Path, plot_records: list) -> dict:
    """独立验证：三个主波形在同一 validation 初态数组上一次性评估。"""
    frozen = _frozen_or_none(output_dir)
    if frozen is None:
        raise RuntimeError("缺少 best_waveform.yaml；请先运行 --stage optimize")
    physics = physics_from_config(cfg)
    heating = heating_from_config(cfg)
    specs = _main_specs(cfg, frozen)
    states = sample_initial_states(cfg, cfg.initial_ensemble.validation_seed,
                                   cfg.initial_ensemble.validation_shots)
    started = time.time()
    results, comparisons = run_validation(cfg, specs, states, physics, heating=heating)
    print(f"[validate] {len(specs)} waveforms x {states['shots']} shots in {time.time()-started:.1f}s"
          + ("（含平均噪声加热）" if heating is not None else ""))

    _write_waveforms_csv(cfg, specs, physics, output_dir)
    _write_validation_csv(output_dir, _validation_records_rows(results))

    # 代表轨迹（成功 / 阈值附近 / 失败）重放
    optimized = results["overlapped_optimized_400us"]
    groups = {"captured": [], "near_threshold": [], "not_captured": []}
    for record in optimized["records"]:
        key = "not_captured" if not record["captured"] else (
            "near_threshold" if record["near_threshold"] else "captured")
        groups[key].append(record["shot_id"])
    keep_ids = set()
    for key, ids in groups.items():
        keep_ids.update(ids[:cfg.output.representative_shots_per_outcome])
    representative = evaluate_waveform_on_states(
        dict(specs["overlapped_optimized_400us"]), states, physics,
        cfg.integration.validation_dt_s, cfg.integration.post_transfer_hold_s,
        keep_trajectory_ids=tuple(sorted(keep_ids)), heating=heating)

    # 理想轨迹功-能 + 一个普通 MC 轨迹功-能
    ideal = ideal_run(dict(specs["overlapped_optimized_400us"]), physics,
                      cfg.integration.validation_dt_s, cfg.integration.post_transfer_hold_s,
                      heating=heating)
    mc_shot = 0
    mc_states = {key: value[mc_shot:mc_shot + 1] for key, value in states.items()
                 if isinstance(value, np.ndarray)}
    mc_eval = evaluate_waveform_on_states(
        dict(specs["overlapped_optimized_400us"]), mc_states, physics,
        cfg.integration.validation_dt_s, cfg.integration.post_transfer_hold_s,
        work_energy=True, keep_trajectory_ids=(0,), heating=heating)
    from .work_energy import work_energy_along_trajectory
    mc_trajectory = mc_eval["trajectories"][0]
    mc_noise_work = mc_trajectory.get("noise_work", {}).get("total_J")
    mc_work = work_energy_along_trajectory(mc_trajectory["time_s"], mc_trajectory["position_m"],
                                           mc_trajectory["velocity_m_per_s"],
                                           build_waveform(specs["overlapped_optimized_400us"], physics),
                                           physics["mass_kg"], physics["slm_depth_j"],
                                           physics["slm_waist_m"], physics["slm_center_m"],
                                           physics["aod_waist_m"], noise_work_j=mc_noise_work)

    # 图
    tables = {name: build_waveform(dict(spec, name=name), physics).dense_table(2001)
              for name, spec in specs.items()}
    plot_records.append(viz.plot_waveform_comparison(tables, output_dir / "waveform_comparison.png"))
    plot_records.append(viz.plot_potential_evolution(specs["overlapped_optimized_400us"], physics,
                                                     output_dir / "potential_evolution.png"))
    labeled = {}
    for shot in sorted(keep_ids):
        record = optimized["records"][shot]
        label = ("near-threshold" if record["near_threshold"] and record["captured"]
                 else "captured" if record["captured"] else "failed")
        labeled[f"shot {shot} ({label})"] = representative["trajectories"][shot]
    if labeled:
        plot_records.append(viz.plot_representative_trajectories(labeled,
                                                                 output_dir / "representative_trajectories.png"))
    plot_records.append(viz.plot_final_phase_space(optimized["records"], physics,
                                                   output_dir / "final_phase_space.png"))
    training_summaries = {}
    frozen_data = yaml.safe_load((output_dir / "best_waveform.yaml").read_text(encoding="utf-8"))
    selection_rate = frozen_data["provenance"].get("selection_capture_rate")
    selection_shots = cfg.initial_ensemble.selection_shots
    selection_captured = round(float(selection_rate) * selection_shots)
    from level1_transfer_1d.analysis import wilson_interval
    selection_low, selection_high = wilson_interval(selection_captured, selection_shots)
    training_summaries["overlapped_optimized_400us"] = {
        "capture_fraction": float(selection_rate),
        "wilson_low": float(selection_low), "wilson_high": float(selection_high),
    }
    validation_summaries = {name: summarize_capture(result["captured_flags"], name)
                            for name, result in results.items()}
    plot_records.append(viz.plot_capture_rate_comparison(validation_summaries, training_summaries,
                                                         output_dir / "capture_rate_comparison.png"))
    plot_records.append(viz.plot_excitation_and_margin(
        {name: result["records"] for name, result in results.items()},
        output_dir / "excitation_and_margin.png"))
    plot_records.append(viz.plot_work_energy_balance(
        {"ideal": ideal["work_energy"], "MC shot 0": mc_work},
        output_dir / "work_energy_balance.png"))

    summaries = {name: summarize_capture(result["captured_flags"], name)
                 for name, result in results.items()}
    excitation_stats = {}
    for name, result in results.items():
        captured_records = [r for r in result["records"] if r["captured"]]
        excitation_stats[name] = {
            "mean_final_excitation_over_depth_captured":
                float(np.mean([r["final_excitation_over_depth"] for r in captured_records]))
                if captured_records else None,
            "mean_margin_captured":
                float(np.mean([r["capture_margin"] for r in captured_records])) if captured_records else None,
            "margin_q10_all": float(np.quantile([r["capture_margin"] for r in result["records"]], 0.10)),
            "mean_excitation_change_uK_all":
                float(np.mean([r["excitation_change_uK"] for r in result["records"]])),
            "near_threshold_count": sum(r["near_threshold"] for r in result["records"]),
            "work_residual_max_over_depth":
                max(r["max_abs_work_energy_residual_over_depth"] for r in result["records"]),
            "work_residual_median_over_depth":
                float(np.median([r["max_abs_work_energy_residual_over_depth"] for r in result["records"]])),
            "hold_peak_to_peak_max_over_depth":
                max(r["hold_peak_to_peak_energy_error_over_depth"] for r in result["records"]),
        }
    excluded_keys = {"x0_m", "v0_m_per_s", "initial_aod_energy_J", "initial_excitation_J"}
    return {"states_meta": {k: v for k, v in states.items() if k not in excluded_keys},
            "results": results, "comparisons": comparisons, "summaries": summaries,
            "excitation_stats": excitation_stats, "validation_states": states,
            "representative_groups": {k: len(v) for k, v in groups.items()}}


def stage_robustness(cfg: Level2Config, output_dir: Path, plot_records: list) -> dict:
    """冻结参数后的温度、对准与步长敏感性检查。"""
    frozen = _frozen_or_none(output_dir)
    if frozen is None:
        raise RuntimeError("缺少 best_waveform.yaml；请先运行 --stage optimize")
    physics = physics_from_config(cfg)
    heating = heating_from_config(cfg)
    specs = _main_specs(cfg, frozen)
    rows = temperature_scan(cfg, specs, heating=heating) + alignment_scan(cfg, specs, heating=heating)
    validation_states = sample_initial_states(cfg, cfg.initial_ensemble.validation_seed,
                                              cfg.initial_ensemble.validation_shots)
    rows += timestep_scan(cfg, specs, validation_states, heating=heating)
    fieldnames = ["scan", "condition", "waveform", "shots", "captured", "capture_fraction",
                  "wilson_low", "wilson_high", "mean_capture_margin", "label_flips_vs_reference",
                  "final_energy_max_abs_diff_uK"]
    with (output_dir / "robustness.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})
    plot_records.append(viz.plot_robustness_heatmaps(rows, output_dir / "robustness_heatmaps.png"))
    plot_records.append(viz.plot_timestep_convergence(
        [r for r in rows if r["scan"] == "timestep"], output_dir / "timestep_convergence.png"))
    return {"rows": rows}


def _build_summary(cfg: Level2Config, preflight, optimize, validate, robustness,
                   image_validation) -> str:
    lines = ["# Level 2 同步波形优化与鲁棒蒙特卡洛验证", "",
             "在 Level 0/1 的一维经典模型内，比较 600 μs 顺序基线、400 μs 压缩顺序基线与 "
             "400 μs 优化同步波形；优化只用 optimization/selection 池，validation 池只报告一次。", ""]
    if preflight:
        # critical_checks 只含 9 个关键项；对外报告与图8/PPT 一致的检查总数
        meta_keys = {"waveform_constraints_all_valid", "critical_checks", "all_passed", "elapsed_s"}
        check_count = len([k for k in preflight if k not in meta_keys])
        lines.append(f"- 预检查：{check_count} 项检查全部通过"
                     f"（{preflight['elapsed_s']:.0f} s，含现有 Level 0/1 测试）。")
    if optimize:
        counts = optimize["counts"]
        lines.append(f"- 优化：共评估 {counts['total']} 个候选（有效 {counts['ok']}，"
                     f"无效 {counts['invalid']}，失败 {counts['failed']}），最终波形已冻结至 best_waveform.yaml。")
        frozen = optimize["frozen"]
        if frozen["type"] == "overlapped":
            lines.append(f"- 冻结参数：move [{frozen['move_start_fraction']:.3f}, {frozen['move_end_fraction']:.3f}]，"
                         f"ramp [{frozen['ramp_start_fraction']:.3f}, {frozen['ramp_end_fraction']:.3f}]，"
                         f"ramp_power = {frozen['ramp_power']:.3f}。")
    if validate:
        noise = cfg.noise_mean_heating
        if noise.enabled:
            lines.append(f"- 平均噪声加热开启（确定性常数，assumed）：注入能量随时间累积，"
                         f"详见 metrics.json 的 noise_mean_heating 段。")
        lines += ["", "## 独立验证（validation 池，一次性）", "",
                  "| 波形 | 俘获 | 俘获率 | Wilson 95% CI | 平均激发(已俘获) |",
                  "|---|---:|---:|---|---:|"]
        for name, summary in validate["summaries"].items():
            stats = validate["excitation_stats"][name]
            exc = stats["mean_final_excitation_over_depth_captured"]
            lines.append(f"| {name} | {summary['captured']}/{summary['shots']} | "
                         f"{summary['capture_fraction']:.4f} | "
                         f"[{summary['wilson_low']:.4f}, {summary['wilson_high']:.4f}] | "
                         f"{exc:.4f} |" if exc is not None else
                         f"| {name} | {summary['captured']}/{summary['shots']} | "
                         f"{summary['capture_fraction']:.4f} | "
                         f"[{summary['wilson_low']:.4f}, {summary['wilson_high']:.4f}] | n/a |")
        lines += ["", "### 配对比较"]
        for key, comparison in validate["comparisons"].items():
            counts = comparison["paired_counts"]
            boot = comparison["capture_rate_difference_bootstrap"]
            lines.append(
                f"- {key}：都成功 {counts['both_captured']}，仅前者 {counts['only_a_captured']}，"
                f"仅后者 {counts['only_b_captured']}，都失败 {counts['both_failed']}；"
                f"俘获率差 {boot['difference']:+.4f}，95% CI [{boot['ci_low']:+.4f}, {boot['ci_high']:+.4f}]"
                f"{'（不含 0，统计显著）' if boot['ci_excludes_zero'] else '（含 0，未达显著）'}。")
    if robustness:
        temperature_rows = [r for r in robustness["rows"] if r["scan"] == "temperature"]
        lines += ["", "## 鲁棒性", "",
                  f"- 温度扫描：{len(temperature_rows)} 个组合（3–10 μK × 三波形，每点 "
                  f"{cfg.robustness.shots_per_condition} shots）。",
                  f"- 末端对准 ±0.10 μm 与步长 0.10/0.05/0.025 μs 敏感性已写入 robustness.csv。"]
    lines += ["", "## 边界", "",
              "本 Level 2 仍是一维经典模型：不含二维/轴向运动、AOD 柱面透镜效应、光子散射、技术噪声、量子内态、"
              "相干性、加热/损失随机过程、多原子相互作用与长距离运输。经典俘获率不可对标实验转移保真度。"]
    images = image_validation.get("images", {}) if image_validation else {}
    if images:
        failed = sorted(name for name, record in images.items() if not record.get("passed"))
        if failed:
            lines.append(f"- 程序化图片检查：{len(images) - len(failed)}/{len(images)} 张通过；"
                         f"未通过：{', '.join(failed)}。")
        else:
            lines.append(f"- 程序化图片检查：{len(images)}/{len(images)} 张全部通过。")
    lines += ["- 未进行人工视觉审阅；程序化检查不能证明标签无重叠或版面美观。", ""]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    cfg = load_config(Path(args.config))
    if args.output_dir:
        cfg = with_output_override(cfg, args.output_dir)
    output_dir = ensure_output_directory(cfg.output.resolved_directory, cfg.project_root)
    save_yaml(output_dir / "config_used.yaml", config_as_dict(cfg))

    plot_records = []
    preflight_data = optimize_data = validate_data = robustness_data = None
    stage = args.stage
    if stage in ("preflight", "all"):
        preflight_data = stage_preflight(cfg, output_dir)
    if stage in ("optimize", "all"):
        optimize_data = stage_optimize(cfg, output_dir, plot_records)
    if stage in ("validate", "all"):
        validate_data = stage_validate(cfg, output_dir, plot_records)
    if stage in ("robustness", "all"):
        robustness_data = stage_robustness(cfg, output_dir, plot_records)

    image_validation = validate_pngs([r for r in plot_records if r])
    # 单独运行某个 stage 时与既有产物合并，而不是覆盖丢失其他阶段的结果
    previous_metrics_path = output_dir / "metrics.json"
    previous_metrics = {}
    if previous_metrics_path.exists():
        try:
            previous_metrics = json.loads(previous_metrics_path.read_text())
        except (json.JSONDecodeError, OSError):
            previous_metrics = {}
    previous_images = {}
    previous_image_path = output_dir / "image_validation.json"
    if previous_image_path.exists():
        try:
            previous_images = json.loads(previous_image_path.read_text()).get("images", {})
        except (json.JSONDecodeError, OSError):
            previous_images = {}
    merged_images = {**previous_images, **image_validation["images"]}
    merged_validation = {"all_passed": all(item["passed"] for item in merged_images.values()),
                         "visual_review_performed": False,
                         "visual_review_findings": [],
                         "images": merged_images}

    if preflight_data is None and (output_dir / "preflight.json").exists():
        preflight_data = json.loads((output_dir / "preflight.json").read_text())
    if optimize_data is None and (output_dir / "best_waveform.yaml").exists():
        counts = {"total": 0, "ok": 0, "invalid": 0, "failed": 0}
        candidates_path = output_dir / "candidates.csv"
        if candidates_path.exists():
            with candidates_path.open(newline="", encoding="utf-8") as handle:
                statuses = [row["status"] for row in csv.DictReader(handle)]
            counts.update(total=len(statuses), ok=statuses.count("ok"),
                          invalid=statuses.count("invalid"), failed=statuses.count("failed"))
        optimize_data = {"counts": counts,
                         "frozen": yaml.safe_load((output_dir / "best_waveform.yaml").read_text())}
    if robustness_data is None and (output_dir / "robustness.csv").exists():
        with (output_dir / "robustness.csv").open(newline="", encoding="utf-8") as handle:
            robustness_data = {"rows": [dict(row) for row in csv.DictReader(handle)]}
    validation_section = None if not validate_data else {
        "summaries": validate_data["summaries"],
        "excitation_stats": validate_data["excitation_stats"],
        "comparisons": validate_data["comparisons"],
        "initial_state_meta": validate_data["states_meta"],
        "representative_groups": validate_data["representative_groups"],
    }

    save_json(output_dir / "image_validation.json", merged_validation)
    save_json(output_dir / "metrics.json", {
        "level2": {"name": cfg.level2.name, "stage": stage},
        "seeds": {"optimization": cfg.initial_ensemble.optimization_seed,
                  "selection": cfg.initial_ensemble.selection_seed,
                  "validation": cfg.initial_ensemble.validation_seed},
        "noise_mean_heating": (heating_from_config(cfg).describe()
                               if cfg.noise_mean_heating.enabled
                               else {"enabled": False, "model": "mean_rate_deterministic"}),
        "preflight": preflight_data,
        "optimization": (optimize_data or {}).get("counts"),
        "frozen_waveform": (optimize_data or {}).get("frozen"),
        "validation": validation_section,
        "robustness": None if not robustness_data else {"rows": robustness_data["rows"]},
        "image_validation_all_passed": merged_validation["all_passed"],
        "visual_review_performed": False,
    })
    save_markdown(output_dir / "summary.md",
                  _build_summary(cfg, preflight_data, optimize_data, validate_data,
                                 robustness_data, merged_validation))
    print(f"Level 2 stage '{stage}' completed: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
