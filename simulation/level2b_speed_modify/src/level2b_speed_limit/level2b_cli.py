"""Level 2B CLI：preflight / scan / compare / all 四个阶段。

用法：
  python -m level2b_speed_limit.level2b_cli \
      --config configs/level2b_speed_limit.yaml --stage all \
      --output-dir outputs/level2b_speed_limit/demo
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Sequence

from level0_static_trap.io_utils import (ensure_output_directory, save_csv,
                                         save_json, save_markdown, save_yaml,
                                         validate_pngs)
from level2_joint_transfer.level2_simulation import physics_from_config
from .level2b_config import config_with_output, load_level2b_config
from .speed_engine import derived_scales
from .speed_scan import (critical_speed_table, estimate_budget,
                         fit_excitation_scaling, load_flags, load_scan_results,
                         run_experiment1, run_experiment2)

PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Level 2B 移动速度边界与变速剖面验证")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", default="all",
                        choices=["preflight", "scan", "compare", "noise",
                                 "maxspeed", "all"])
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--skip-level2-preflight", action="store_true",
                        help="跳过 Level 2 十项复跑（调试用；正式运行不得使用）")
    parser.add_argument("--skip-existing-tests", action="store_true",
                        help="跳过既有测试运行（调试/冒烟用；正式运行不得使用）")
    return parser


def _no_nan(value):
    """递归把 NaN/Infinity 替换为 None，保证严格 JSON。"""
    if isinstance(value, float):
        return None if (math.isnan(value) or math.isinf(value)) else value
    if isinstance(value, dict):
        return {k: _no_nan(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_no_nan(v) for v in value]
    return value


def _critical_speeds_csv(table: list[dict]) -> dict:
    """把临界速度表展平为 CSV 列。"""
    columns = {"structure": [], "temperature_uK": [], "boundary_covered": [],
               "v_avg_max_mm_per_s": [], "min_capture_rate": []}
    for level_key in ("v95", "v50"):
        for method in ("bootstrap", "logistic"):
            columns[f"{level_key}_{method}_mm_per_s"] = []
            columns[f"{level_key}_{method}_ci_low"] = []
            columns[f"{level_key}_{method}_ci_high"] = []
            columns[f"{level_key}_{method}_status"] = []
    for entry in table:
        columns["structure"].append(entry["structure"])
        columns["temperature_uK"].append(entry["temperature_uK"])
        columns["boundary_covered"].append(entry["boundary_covered"])
        columns["v_avg_max_mm_per_s"].append(entry["v_avg_max_mm_per_s"])
        columns["min_capture_rate"].append(entry["min_capture_rate"])
        for level_key in ("v95", "v50"):
            for method in ("bootstrap", "logistic"):
                item = entry.get(f"{level_key}_{method}", {})
                point = item.get("point_estimate")
                columns[f"{level_key}_{method}_mm_per_s"].append(
                    point * 1e3 if point is not None else None)
                columns[f"{level_key}_{method}_ci_low"].append(
                    item.get("ci_low", None) if point is not None else None)
                columns[f"{level_key}_{method}_ci_high"].append(
                    item.get("ci_high", None) if point is not None else None)
                columns[f"{level_key}_{method}_status"].append(item.get("status", "missing"))
    lengths = {len(v) for v in columns.values()}
    assert len(lengths) == 1, f"临界速度表列长不一致：{lengths}"
    return {k: ["" if x is None else x for x in v] for k, v in columns.items()}


def stage_preflight(bcfg: dict, cfg2, out_dir: Path, include_level2=True,
                    include_tests=True) -> dict:
    """预检查阶段。"""
    from .preflight2b import run_preflight
    result = run_preflight(bcfg, cfg2, include_level2=include_level2,
                           include_existing_tests=include_tests)
    save_json(out_dir / "preflight.json", _no_nan(result))
    if not result["all_passed"]:
        failed = [name for name in result["critical_checks"]
                  if not result["checks"].get(name, {}).get(
                      "passed", result["checks"].get(name, {}).get("all_passed", False))]
        print(f"[preflight] 失败项：{failed}", file=sys.stderr)
    return result


def stage_scan(bcfg: dict, cfg2, out_dir: Path, num_procs=None) -> dict:
    """扫描阶段：实验 1 + 实验 2 + 临界速度表 + 标度律。"""
    from .speed_scan import (run_experiment1, run_experiment2)
    physics = physics_from_config(cfg2)
    scales = derived_scales(physics)
    scan_cfg = bcfg["speed_scan"]
    merged = {**bcfg, **bcfg["initial_ensemble"], "speed_scan": scan_cfg}
    merged = dict(bcfg)
    merged.update({
        "pure_move_cruise_speeds_mm_per_s": scan_cfg["pure_move_cruise_speeds_mm_per_s"],
        "end_to_end_durations_us": scan_cfg["end_to_end_durations_us"],
        "scan_temperatures_uK": bcfg["initial_ensemble"]["scan_temperatures_uK"],
        "main_temperature_uK": bcfg["initial_ensemble"]["main_temperature_uK"],
        "scan_shots_per_point": bcfg["initial_ensemble"]["scan_shots_per_point"],
        "scan_seed": bcfg["initial_ensemble"]["scan_seed"],
        "boundary_seed": bcfg["initial_ensemble"]["boundary_seed"],
        "trapezoid_accel_fraction": scan_cfg["trapezoid_accel_fraction"],
        "stop_when_capture_below": scan_cfg["stop_when_capture_below"],
        "base_dt_us": bcfg["integration"]["base_dt_us"],
        "post_transfer_hold_us": bcfg["integration"]["post_transfer_hold_us"],
        "min_steps_per_move": bcfg["integration"]["min_steps_per_move"],
        "max_step_displacement_over_waist": bcfg["integration"]["max_step_displacement_over_waist"],
        "work_energy_shots": bcfg["integration"].get("work_energy_shots", 4),
        "follow_loss": bcfg["follow_loss"],
        "capture_target_levels": scan_cfg["capture_target_levels"],
        "bootstrap_samples": bcfg["comparison"]["bootstrap_samples"],
    })
    out_csv = out_dir / "speed_scan_results.csv"
    t0 = time.time()
    rows1 = run_experiment1(cfg2, merged, physics, scales, out_csv, out_dir,
                            num_procs=num_procs)
    rows2 = run_experiment2(cfg2, merged, physics, scales, out_csv, out_dir,
                            num_procs=num_procs)
    elapsed = time.time() - t0
    rows = load_scan_results(out_csv)
    flags = load_flags(out_dir)
    table = critical_speed_table(rows, flags, merged)
    fit = fit_excitation_scaling(rows)
    budget = estimate_budget(merged)
    if table:
        save_csv(out_dir / "critical_speeds.csv", _critical_speeds_csv(table))
    save_json(out_dir / "scaling_fit.json", _no_nan(fit))
    ok = sum(1 for r in rows if r.get("status") == "ok")
    skipped = sum(1 for r in rows if r.get("status") == "skipped_early_stop")
    failed = sum(1 for r in rows if r.get("status") == "failed")
    return {"rows": rows, "critical_table": table, "scaling": fit,
            "budget": budget, "elapsed_s": elapsed,
            "counts": {"ok": ok, "skipped_early_stop": skipped, "failed": failed}}


def stage_compare(bcfg: dict, cfg2, out_dir: Path, critical_table: list[dict],
                  num_procs=None) -> dict:
    """比较阶段：实验 3 双口径配对比较与独立复核。"""
    from .profile_compare import COMPARISON_COLUMNS, run_comparison
    merged = dict(bcfg)
    merged.update({
        "scan_seed": bcfg["initial_ensemble"]["scan_seed"],
        "boundary_seed": bcfg["initial_ensemble"]["boundary_seed"],
        "scan_temperatures_uK": bcfg["initial_ensemble"]["scan_temperatures_uK"],
        "trapezoid_accel_fraction": bcfg["speed_scan"]["trapezoid_accel_fraction"],
        "base_dt_us": bcfg["integration"]["base_dt_us"],
        "min_steps_per_move": bcfg["integration"]["min_steps_per_move"],
        "max_step_displacement_over_waist": bcfg["integration"]["max_step_displacement_over_waist"],
        "follow_loss": bcfg["follow_loss"],
    })
    out_csv = out_dir / "comparison_results.csv"
    if not out_csv.is_file():
        import csv
        with out_csv.open("w", encoding="utf-8", newline="") as handle:
            csv.DictWriter(handle, fieldnames=COMPARISON_COLUMNS).writeheader()
    result = run_comparison(cfg2, merged, physics_from_config(cfg2),
                            critical_table, out_csv, num_procs=num_procs)
    from .speed_scan import load_scan_results
    comparison_rows = load_scan_results(out_csv)
    result["rows"] = comparison_rows
    return result


def stage_noise(bcfg: dict, cfg2, out_dir: Path, num_procs=None) -> dict:
    """噪声阶段（实验 4）：v95 偏移 / 开关配对 / 幅度敏感性 / 通道分解。

    结果增量写入 noise_speed_results.csv（断点续跑跳过已完成的 4a 点），
    v95 双方法估计用旗标档案重算后写入 noise_v95.json。
    """
    import csv
    from .noise_speed_validation import (NOISE_COLUMNS, compute_noise_v95,
                                         run_noise_stage)
    from .speed_engine import derived_scales, physics_from_config as _pfc
    physics = _pfc(cfg2)
    scales = derived_scales(physics)
    out_csv = out_dir / "noise_speed_results.csv"
    if not out_csv.is_file():
        with out_csv.open("w", encoding="utf-8", newline="") as handle:
            csv.DictWriter(handle, fieldnames=NOISE_COLUMNS).writeheader()
    summary = run_noise_stage(cfg2, bcfg, physics, scales, out_csv, out_dir,
                              num_procs=num_procs)
    # 以落盘数据 + 旗标档案重算 v95（与 run 内估计独立，防内存/磁盘不一致）
    summary["v95_table"] = compute_noise_v95(out_dir)
    save_json(out_dir / "noise_v95.json", _no_nan(summary))
    return summary


def stage_maxspeed(bcfg: dict, cfg2, out_dir: Path, num_procs=None) -> dict:
    """最高速度测算（实验 5）：默认噪声下硬件条件 v95 + 同时长剖面对比。

    结果增量写入 max_speed_hardware.csv / max_speed_profiles.csv，
    汇总（v95 表 + 剖面对比）写入 max_speed_summary.json。
    """
    from .max_speed_probe import run_max_speed_stage
    from .speed_engine import physics_from_config as _pfc
    physics = _pfc(cfg2)
    summary = run_max_speed_stage(cfg2, bcfg, physics, out_dir, num_procs=num_procs)
    save_json(out_dir / "max_speed_summary.json", _no_nan(summary))
    return summary


def _reshot_agreement(comparison_rows: list[dict]) -> list[dict]:
    """扫描池 vs 独立复核的俘获率一致性（同口径同点同剖面）。"""
    from level1_transfer_1d.analysis import wilson_interval
    scan = {(r["caliber"], r["v_target_mm_per_s"], r["profile"]): r
            for r in comparison_rows if r["dataset"] == "scan_pool" and r["status"] == "ok"}
    out = []
    for r in comparison_rows:
        if r["dataset"] != "boundary_reshot" or r["status"] != "ok":
            continue
        key = (r["caliber"], r["v_target_mm_per_s"], r["profile"])
        base = scan.get(key)
        if base is None:
            continue
        s1, n1 = int(float(base["captured"])), int(float(base["shots"]))
        s2, n2 = int(float(r["captured"])), int(float(r["shots"]))
        lo1, hi1 = wilson_interval(s1, n1)
        lo2, hi2 = wilson_interval(s2, n2)
        out.append({"caliber": r["caliber"], "v_target_mm_per_s": float(r["v_target_mm_per_s"]),
                    "profile": r["profile"],
                    "scan_rate": s1 / n1, "reshot_rate": s2 / n2,
                    "scan_ci": [float(lo1), float(hi1)],
                    "reshot_ci": [float(lo2), float(hi2)],
                    "ci_overlap": bool(not (hi1 < lo2 or hi2 < lo1))})
    return out


def _build_summary(bcfg, preflight, scan, compare, image_validation, scales,
                   noise=None, maxspeed=None) -> str:
    lines = ["# Level 2B AOD 移动速度边界与变速剖面验证",
             "",
             "在 Level 0/1/2 一维经典模型内，把移动速度作为显式扫描变量：实验 1 纯移动跟随上限、实验 2 端到端速度上限（三结构 × 三温度）、实验 3 恒速 vs 变速双口径配对比较。", ""]
    if preflight is not None:
        lines.append(f"- 预检查：{'全部通过' if preflight['all_passed'] else '存在失败项'}"
                     f"（{preflight['elapsed_s']:.0f} s）。")
    if scan is not None:
        c = scan["counts"]
        lines.append(f"- 扫描：{c['ok']} 点完成、{c['skipped_early_stop']} 点曲线级早停跳过、"
                     f"{c['failed']} 点失败（耗时 {scan['elapsed_s']:.0f} s）。")
        covered = [t for t in scan["critical_table"] if t["boundary_covered"]]
        if covered:
            for t in covered:
                v95 = t.get("v95_bootstrap", {})
                point = v95.get("point_estimate")
                text = f"{point*1e3:.1f} mm/s" if point is not None else "估计失败"
                lines.append(f"- {t['structure']} @ {t['temperature_uK']:g} μK：v95 ≈ {text}，边界已覆盖。")
        else:
            lines.append("- 所有曲线俘获率均未跌破 95%：扫描区间未覆盖失效边界，"
                         "结论只能写 v95 > 最大扫描速度，不得宣称已测得最大速度。")
        if scan["scaling"].get("fitted"):
            s = scan["scaling"]
            lines.append(f"- 实验 1 标度律：ε_med ∝ v^{s['alpha']:.2f}"
                         f"（95% CI [{s['alpha_ci_low']:.2f}, {s['alpha_ci_high']:.2f}]，R²={s['r_squared']:.3f}）。")
    if compare is not None:
        ref = compare["reference"]
        lines.append(f"- 实验 3 参考曲线：{ref.get('structure', '默认')} @ {ref.get('temperature_uK', '?'):g} μK"
                     f"（{ref['method']}，v_ref = {ref['v95_m_per_s']*1e3:.1f} mm/s）。")
        agree = compare.get("reshot_agreement", [])
        if agree:
            overlap = sum(1 for a in agree if a["ci_overlap"])
            lines.append(f"- 独立复核：{overlap}/{len(agree)} 个 (口径, 速度, 剖面) 组合的 95% CI 与扫描池重叠。")
    if noise is not None:
        lines.append("")
        lines.append("## 实验 4：平均强度噪声加热 × 移动速度")
        lines.append("")
        defaults = noise.get("defaults", {}).get("derived_from", {})
        if defaults:
            lines.append(f"- 默认幅度（assumed_sensitivity_only）：反冲 {defaults.get('recoil_power_uK_per_s', 0):.2g} μK/s、"
                         f"参量 Γ_par = {defaults.get('parametric_rate_per_s', 0):.1f} s⁻¹、"
                         f"指向 {defaults.get('pointing_power_uK_per_s', 0):.2g} μK/s。")
        for entry in noise.get("v95_table", []):
            boot = entry.get("v95_bootstrap", {})
            point = boot.get("point_estimate")
            text = f"{point*1e3:.1f} mm/s" if point is not None else "未测得/未覆盖"
            lines.append(f"- 4a v95[{entry['noise_mode']}] = {text}"
                         f"（covered={entry['boundary_covered']}）。")
        sig = [a for a in noise.get("onoff", []) if a.get("ci_excludes_zero")]
        if sig:
            for a in sig:
                lines.append(f"- 4b 显著配对差 @ {a['v_peak_mm_per_s']:g} mm/s, {a['temperature_uK']:g} μK："
                             f"{a['difference']:+.4f} CI[{a['ci_low']:+.4f},{a['ci_high']:+.4f}]。")
        else:
            lines.append("- 4b 开/关配对差全部 CI 含 0：默认幅度噪声在单次转移内不可分辨。")
        sens = noise.get("sensitivity", [])
        if sens:
            first_bad = next((s for s in sens if s["capture_rate"] < 0.99), None)
            if first_bad:
                lines.append(f"- 4c 俘获率退化始于 ×{first_bad['scale']:g} 档"
                             f"（rate={first_bad['capture_rate']:.3f}）。")
            else:
                lines.append("- 4c 全部放大档俘获率仍 ≥ 0.99。")
        chans = {c["mode"]: c for c in noise.get("channels", [])}
        if "parametric_only" in chans and "off" in chans:
            chan_scale = bcfg["noise_speed_validation"]["channel_breakdown"]["scale_factor"]
            lines.append(f"- 4d 通道分解（×{chan_scale:g}）：仅参量 rate={chans['parametric_only']['capture_rate']:.3f}、"
                         f"仅反冲 {chans.get('recoil_only', {}).get('capture_rate', float('nan')):.3f}、"
                         f"仅指向 {chans.get('pointing_only', {}).get('capture_rate', float('nan')):.3f}、"
                         f"关 {chans['off']['capture_rate']:.3f}。")
    if maxspeed is not None:
        for entry in maxspeed.get("v95_table", []):
            point = entry.get("v95_bootstrap", {}).get("point_estimate")
            eta = entry.get("eta95_peak")
            text = (f"{point * 1e3:.1f} mm/s（η95 ≈ {eta:.2f}）" if point is not None
                    else "未测得/未覆盖")
            lines.append(f"- 5a v95[{entry['condition']}]（默认噪声 ×1）= {text}。")
        by_duration = {}
        for p in maxspeed.get("profiles", []):
            by_duration.setdefault(p["duration_us"], []).append(p)
        for t_us, items in sorted(by_duration.items()):
            rates = "、".join(f"{p['profile']} {p['capture_rate']:.3f}" for p in items)
            lines.append(f"- 5b 同时长剖面对比 T={t_us:g} μs：{rates}。")
    lines.append(f"- 图片程序化检查：{'全部通过' if image_validation.get('all_passed') else '存在失败'}；"
                 "未进行人工视觉审阅（visual_review_performed: false）。")
    lines += ["", "## 边界", "",
              "本 Level 仍是一维经典模型且不含 AOD 器件带宽/RF 扫描速率限制：测得的是该模型下的动力学速度上限，不是器件规格。", ""]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    bcfg_raw, cfg2 = load_level2b_config(args.config)
    out_dir = Path(args.output_dir) if args.output_dir \
        else Path(bcfg_raw["output"]["directory"])
    if not out_dir.is_absolute():
        out_dir = PACKAGE_ROOT / out_dir
    resolved = out_dir.resolve()
    if "sources" in resolved.parts:
        print("输出目录不得位于 sources/ 内", file=sys.stderr)
        return 2
    ensure_output_directory(resolved, PACKAGE_ROOT)
    bcfg = config_with_output(bcfg_raw, resolved)
    save_yaml(resolved / "config_used.yaml", bcfg)

    preflight = scan = compare = None
    image_validation = {"all_passed": False}

    if args.stage in ("preflight", "all"):
        preflight = stage_preflight(bcfg, cfg2, resolved,
                                    include_level2=not args.skip_level2_preflight,
                                    include_tests=not args.skip_existing_tests)
        if args.stage == "all" and not preflight["all_passed"]:
            print("[abort] 预检查未全部通过，按规程终止后续阶段", file=sys.stderr)
            return 1
        if args.stage == "preflight":
            return 0 if preflight["all_passed"] else 1

    if args.stage in ("scan", "all"):
        if args.stage == "scan" and not (resolved / "preflight.json").is_file():
            print("[scan] 未找到 preflight.json，请先运行 preflight", file=sys.stderr)
            return 2
        scan = stage_scan(bcfg, cfg2, resolved)

    if args.stage in ("compare", "all"):
        if scan is None:
            rows = load_scan_results(resolved / "speed_scan_results.csv")
            flags = load_flags(resolved)
            merged_cfg = {**bcfg}
            scan = {"critical_table": critical_speed_table(rows, flags, merged_cfg)}
        compare = stage_compare(bcfg, cfg2, resolved, scan["critical_table"])
        compare["reshot_agreement"] = _reshot_agreement(compare["rows"])

    noise = None
    if args.stage in ("noise", "all"):
        noise = stage_noise(bcfg, cfg2, resolved)

    maxspeed = None
    if args.stage in ("maxspeed", "all"):
        maxspeed = stage_maxspeed(bcfg, cfg2, resolved)

    if args.stage == "all":
        physics = physics_from_config(cfg2)
        scales = derived_scales(physics)
        v_ref = compare["reference"]["v95_m_per_s"] if compare else 0.05
        from .level2b_visualization import generate_all
        plot_records = generate_all(cfg2, bcfg, physics, resolved,
                                    scan["critical_table"], compare["rows"], v_ref,
                                    noise_summary=noise)
        image_validation = validate_pngs(plot_records)
        save_json(resolved / "image_validation.json", _no_nan(image_validation))

        metrics = {
            "level2b": bcfg["level2b"], "stage": args.stage,
            "seeds": {k: bcfg["initial_ensemble"][k]
                      for k in ("scan_seed", "boundary_seed")},
            "derived_scales": scales,
            "preflight_all_passed": preflight["all_passed"] if preflight else None,
            "scan": {"counts": scan["counts"], "budget": scan["budget"],
                     "elapsed_s": scan["elapsed_s"]},
            "critical_speeds": _no_nan(scan["critical_table"]),
            "scaling_fit": _no_nan(scan["scaling"]),
            "comparison": {"reference": _no_nan(compare["reference"]),
                           "reshot_agreement": _no_nan(compare.get("reshot_agreement", []))},
            "noise": _no_nan(noise) if noise else None,
            "max_speed": _no_nan(maxspeed) if maxspeed else None,
            "boundary_covered_any": any(t["boundary_covered"] for t in scan["critical_table"]),
            "image_validation_all_passed": image_validation.get("all_passed", False),
            "visual_review_performed": False,
        }
        save_json(resolved / "metrics.json", _no_nan(metrics))
        summary = _build_summary(bcfg, preflight, scan, compare, image_validation,
                                 scales, noise, maxspeed=maxspeed)
        save_markdown(resolved / "summary.md", summary)
        print(f"[done] 输出目录：{resolved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
