"""Level 2 预检查：不过关不得启动优化。全部结果写入 preflight.json。"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from level0_static_trap.integrators import velocity_verlet
from level0_static_trap.potentials import gaussian_force
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from .config import Level2Config
from .level2_simulation import (baseline_specs, evaluate_waveform_on_states, ideal_run,
                                physics_from_config, sample_initial_states)
from .waveforms import validate_waveform

IDEAL_RESIDUAL_TOLERANCE = 5.0e-4
MC_RESIDUAL_MEDIAN_TOLERANCE = 1.0e-3
LABEL_FLIP_TOLERANCE = 0.02
HOLD_PEAK_TO_PEAK_TOLERANCE = 1.0e-4

# 测试目录与 Level 1 配置属于包本身，按包位置解析，避免依赖 config 文件所在目录
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
LEVEL1_CONFIG = PACKAGE_ROOT / "configs" / "transfer_1d.yaml"


def _check_static_limit(physics: dict, dt_s: float) -> dict:
    """时变积分器退化为静态力时应与 Level 0 积分器逐点一致。"""
    times = np.arange(1001) * dt_s
    x0 = physics["aod_initial_center_m"] + 0.1 * physics["aod_waist_m"]
    force = lambda x: gaussian_force(x, physics["aod_initial_depth_j"],
                                     physics["aod_waist_m"], physics["aod_initial_center_m"])
    old = velocity_verlet(force, physics["mass_kg"], x0, 0.0, times)
    new = velocity_verlet_time_dependent(lambda x, t: force(x), physics["mass_kg"], x0, 0.0, times)
    position_ok = bool(np.allclose(old[1], new[1], rtol=1e-13, atol=1e-18))
    velocity_ok = bool(np.allclose(old[2], new[2], rtol=1e-13, atol=1e-15))
    return {"position_match": position_ok, "velocity_match": velocity_ok,
            "passed": position_ok and velocity_ok}


def _check_existing_tests(project_root: Path) -> dict:
    """运行现有 Level 0/Level 1 测试并记录通过数。

    排除 slow 标记的集成测试：它们内部会调用本 CLI 的 preflight/all stage，
    不排除会在 CLI 冒烟测试中造成递归 pytest。
    """
    result = subprocess.run([sys.executable, "-m", "pytest", str(PACKAGE_ROOT / "tests"), "-q",
                             "--no-header", "-m", "not slow"],
                            cwd=project_root, capture_output=True, text=True, timeout=1800)
    tail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
    return {"exit_code": result.returncode, "summary_line": tail,
            "passed": result.returncode == 0}


def _check_level1_reproduction(cfg: Level2Config) -> dict:
    """用固定种子复现 Level 1 sequential 600 us 的蒙特卡洛结果。"""
    from level1_transfer_1d.config import load_config as load_level1
    from level1_transfer_1d.simulation import monte_carlo as level1_monte_carlo
    level1_cfg = load_level1(LEVEL1_CONFIG)
    records = level1_monte_carlo(level1_cfg)
    captured = sum(r["captured_in_slm"] for r in records)
    fraction = captured / len(records)
    historical_path = cfg.project_root.parent / "level1_transfer_1d" / "outputs" / \
        "level1_transfer_1d" / "demo" / "metrics.json"
    comparison = {"reproduced_shots": len(records), "reproduced_captured": captured,
                  "reproduced_fraction": fraction, "seed": level1_cfg.monte_carlo.seed}
    if historical_path.exists():
        historical = json.loads(historical_path.read_text())
        historical_fraction = historical["monte_carlo"]["classical_capture_fraction"]
        comparison.update({
            "historical_fraction": historical_fraction,
            "matches_historical": bool(fraction == historical_fraction),
            "passed": bool(fraction == historical_fraction),
        })
    else:
        comparison.update({"historical_available": False,
                           "passed": bool(fraction == 1.0 or fraction >= 0.9)})
    return comparison


def _convergence_order(errors: list[float]) -> float:
    """由相邻误差估计收敛阶 log2(e_coarse/e_fine)。"""
    if errors[1] <= 0 or errors[0] <= 0:
        return float("nan")
    return float(np.log2(errors[0] / errors[1]))


def run_preflight(cfg: Level2Config) -> dict:
    """按顺序执行十项预检查；关键项失败时抛出异常阻止优化。"""
    physics = physics_from_config(cfg)
    hold_s = cfg.integration.post_transfer_hold_s
    constraint = cfg.waveforms.as_constraint_dict()
    nominal_overlapped = {
        "type": "overlapped", "name": "nominal_overlapped", "duration_s": cfg.waveforms.optimized_duration_s,
        "move_start_fraction": 0.0, "move_end_fraction": cfg.waveforms.level1_move_fraction,
        "ramp_start_fraction": cfg.waveforms.level1_move_fraction, "ramp_end_fraction": 1.0,
        "ramp_power": 2.0,
    }
    specs = dict(baseline_specs(cfg))
    specs["nominal_overlapped_400us"] = nominal_overlapped

    checks = {}
    # 1. 现有测试
    checks["existing_tests"] = _check_existing_tests(cfg.project_root)
    # 2. 静态极限
    checks["static_limit"] = _check_static_limit(physics, cfg.integration.validation_dt_s)
    # 3. Level 1 复现
    checks["level1_reproduction"] = _check_level1_reproduction(cfg)
    # 4. 波形端点/范围/单调性
    waveform_checks = {}
    from .level2_simulation import build_waveform
    for name, spec in specs.items():
        waveform = build_waveform(spec, physics)
        valid, reasons = validate_waveform(waveform, constraint)
        waveform_checks[name] = {"valid": valid, "reasons": reasons,
                                 "smoothness": waveform.smoothness_metrics(2001)}
    checks["waveform_constraints"] = waveform_checks
    # 5-6. 理想初态三波形 + 功-能平衡
    ideal_results, residual_checks = {}, {}
    depth_scale = max(physics["slm_depth_j"], physics["aod_initial_depth_j"])
    for name, spec in specs.items():
        run = ideal_run(spec, physics, cfg.integration.validation_dt_s, hold_s)
        ideal_results[name] = {
            "final_slm_energy_over_depth": run["records"][0]["final_slm_energy_uK"] * 1.380649e-29 * 1e6 / physics["slm_depth_j"],
            "captured": run["records"][0]["captured"],
            "max_abs_residual_over_depth": run["work_energy"]["max_abs_residual_J"] / depth_scale,
        }
        residual_checks[name] = ideal_results[name]["max_abs_residual_over_depth"] < IDEAL_RESIDUAL_TOLERANCE
    checks["ideal_runs"] = ideal_results
    checks["ideal_work_energy"] = {"per_waveform_passed": residual_checks,
                                   "tolerance": IDEAL_RESIDUAL_TOLERANCE,
                                   "passed": all(residual_checks.values())}
    # 7. 理想轨迹步长收敛（0.10 / 0.05 / 0.025）
    convergence = {}
    for name, spec in specs.items():
        fine = ideal_run(spec, physics, cfg.integration.convergence_dt_s, hold_s)
        fine_state = (fine["trajectories"][0]["position_m"][-1], fine["trajectories"][0]["velocity_m_per_s"][-1])
        errors = []
        for dt in (cfg.integration.optimization_dt_s, cfg.integration.validation_dt_s):
            coarse = ideal_run(spec, physics, dt, hold_s)
            trajectory = coarse["trajectories"][0]
            dx = abs(float(trajectory["position_m"][-1]) - float(fine_state[0]))
            dv = abs(float(trajectory["velocity_m_per_s"][-1]) - float(fine_state[1]))
            errors.append((dx, dv))
        convergence[name] = {
            "position_error_coarse_m": errors[0][0], "position_error_fine_m": errors[1][0],
            "velocity_error_coarse_m_per_s": errors[0][1], "velocity_error_fine_m_per_s": errors[1][1],
            "position_convergence_order": _convergence_order([errors[0][0], errors[1][0]]),
            "velocity_convergence_order": _convergence_order([errors[0][1], errors[1][1]]),
        }
    position_orders = [c["position_convergence_order"] for c in convergence.values()]
    checks["ideal_timestep_convergence"] = {
        "per_waveform": convergence,
        "min_position_order": float(np.nanmin(position_orders)),
        "passed": bool(np.nanmin(position_orders) > 1.0),
    }
    # 8. MC 子集配对步长检查（128 shots，0.05 vs 0.025）
    pool = sample_initial_states(cfg, cfg.initial_ensemble.optimization_seed,
                                 min(128, cfg.initial_ensemble.optimization_shots))
    dt_label_checks = {}
    reference_waveforms = dict(baseline_specs(cfg))
    reference_waveforms["nominal_overlapped_400us"] = nominal_overlapped
    for name, spec in reference_waveforms.items():
        coarse = evaluate_waveform_on_states(spec, pool, physics, cfg.integration.validation_dt_s, hold_s)
        fine = evaluate_waveform_on_states(spec, pool, physics, cfg.integration.convergence_dt_s, hold_s)
        flags_c = np.array(coarse["captured_flags"])
        flags_f = np.array(fine["captured_flags"])
        flips = flags_c != flags_f
        energies_c = np.array([r["final_slm_energy_uK"] for r in coarse["records"]])
        margins_c = np.array([r["capture_margin"] for r in coarse["records"]])
        dt_label_checks[name] = {
            "label_flip_rate": float(flips.mean()),
            "flip_count": int(flips.sum()),
            "median_abs_margin_of_flipped": float(np.median(np.abs(margins_c[flips]))) if flips.any() else None,
            "max_final_energy_abs_diff_uK": float(np.max(np.abs(
                energies_c - np.array([r["final_slm_energy_uK"] for r in fine["records"]])))),
        }
    flip_rates = [v["label_flip_rate"] for v in dt_label_checks.values()]
    residual_medians = {}
    for name in reference_waveforms:
        evaluation = evaluate_waveform_on_states(
            {**reference_waveforms[name], "name": name}, pool, physics,
            cfg.integration.validation_dt_s, hold_s, work_energy=True)
        residuals = [r["max_abs_work_energy_residual_over_depth"] for r in evaluation["records"]]
        residual_medians[name] = {
            "median": float(np.median(residuals)),
            "p95": float(np.percentile(residuals, 95)),
            "max": float(np.max(residuals)),
        }
    checks["mc_timestep_paired"] = {
        "per_waveform": dt_label_checks, "shots": pool["shots"],
        "max_label_flip_rate": float(max(flip_rates)),
        "passed": bool(max(flip_rates) <= LABEL_FLIP_TOLERANCE),
    }
    checks["mc_work_energy_residual"] = {
        "per_waveform": residual_medians,
        "tolerance": MC_RESIDUAL_MEDIAN_TOLERANCE,
        "passed": bool(all(v["median"] < MC_RESIDUAL_MEDIAN_TOLERANCE for v in residual_medians.values())),
    }
    # 9. 保持段能量检查（理想轨迹）
    hold_checks = {}
    for name, spec in reference_waveforms.items():
        run = ideal_run(spec, physics, cfg.integration.validation_dt_s, hold_s)
        record = run["records"][0]
        hold_checks[name] = record["hold_peak_to_peak_energy_error_over_depth"]
    checks["hold_segment"] = {
        "per_waveform_peak_to_peak_over_depth": hold_checks,
        "tolerance": HOLD_PEAK_TO_PEAK_TOLERANCE,
        "passed": bool(all(v is not None and v < HOLD_PEAK_TO_PEAK_TOLERANCE for v in hold_checks.values())),
    }
    # 10. 关键项汇总
    critical = ["existing_tests", "static_limit", "level1_reproduction", "ideal_work_energy",
                "ideal_timestep_convergence", "mc_timestep_paired", "mc_work_energy_residual",
                "hold_segment"]
    checks["waveform_constraints_all_valid"] = bool(
        all(v["valid"] for v in waveform_checks.values()))
    critical.append("waveform_constraints_all_valid")
    checks["critical_checks"] = critical

    def _passed(name):
        value = checks[name]
        return bool(value) if isinstance(value, bool) else bool(value.get("passed", False))

    checks["all_passed"] = bool(all(_passed(name) for name in critical))
    return checks
