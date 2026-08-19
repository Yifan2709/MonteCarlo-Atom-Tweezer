"""Level 2 测试：波形、功-能、统计隔离、判据与 CLI 冒烟。"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from level2_joint_transfer.config import load_config
from level2_joint_transfer.level2_simulation import (
    baseline_specs, build_waveform, evaluate_waveform_on_states, ideal_run, is_captured,
    physics_from_config, sample_initial_states)
from level2_joint_transfer.optimizer import _evaluate_job, frozen_spec
from level2_joint_transfer.preflight import _check_static_limit, _convergence_order
from level2_joint_transfer.statistics import paired_bootstrap_difference, paired_counts, summarize_capture
from level2_joint_transfer.waveforms import (OverlappedWaveform, SequentialWaveform,
                                             SplineWaveform, validate_waveform)
from level2_joint_transfer.work_energy import (hold_segment_energy_error,
                                               work_energy_along_trajectory)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "configs" / "level2_joint_transfer.yaml"
CONSTRAINT = {"min_segment_fraction": 0.10, "ramp_power_bounds": [0.5, 4.0]}


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG)


@pytest.fixture(scope="module")
def physics(cfg):
    return physics_from_config(cfg)


def nominal_overlapped(cfg):
    return {"type": "overlapped", "duration_s": cfg.waveforms.optimized_duration_s,
            "move_start_fraction": 0.05, "move_end_fraction": 0.60,
            "ramp_start_fraction": 0.40, "ramp_end_fraction": 1.0, "ramp_power": 2.0}


# ---------------------------------------------------------------- 1. 端点
def test_waveform_endpoints(cfg, physics):
    specs = dict(baseline_specs(cfg))
    specs["overlapped"] = nominal_overlapped(cfg)
    for name, spec in specs.items():
        waveform = build_waveform(spec, physics)
        t0, t1 = 0.0, waveform.duration_s
        assert math_close(float(waveform.center(t0)), physics["aod_initial_center_m"])
        assert abs(float(waveform.center(t1))) < 1e-13
        assert math_close(float(waveform.depth(t0)), physics["aod_initial_depth_j"])
        assert float(waveform.depth(t1)) < 1e-24


def math_close(a, b, rel=1e-12):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


# ---------------------------------------------------------------- 2. 单调与界限
def test_overlapped_monotone_and_bounded(cfg, physics):
    waveform = build_waveform(nominal_overlapped(cfg), physics)
    table = waveform.dense_table(4001)
    assert np.all(np.diff(table["aod_center_m"]) <= 1e-15)
    assert np.all(np.diff(table["aod_depth_J"]) <= 1e-21)
    assert np.all(table["aod_center_m"] >= -1e-15)
    assert np.all(table["aod_center_m"] <= physics["aod_initial_center_m"] + 1e-15)
    assert np.all(table["aod_depth_J"] >= -1e-24)


# ---------------------------------------------------------------- 3. 非法参数被拒绝
@pytest.mark.parametrize("params,fragment", [
    ({"move_start_fraction": 0.5, "move_end_fraction": 0.4, "ramp_start_fraction": 0.0,
      "ramp_end_fraction": 0.9, "ramp_power": 2.0}, "移动窗口"),
    ({"move_start_fraction": 0.0, "move_end_fraction": 0.9, "ramp_start_fraction": 0.0,
      "ramp_end_fraction": 0.04, "ramp_power": 2.0}, "降深窗口"),
    ({"move_start_fraction": 0.0, "move_end_fraction": 0.9, "ramp_start_fraction": 0.0,
      "ramp_end_fraction": 0.9, "ramp_power": 0.6}, "奇异"),
    ({"move_start_fraction": 0.0, "move_end_fraction": 0.9, "ramp_start_fraction": 0.0,
      "ramp_end_fraction": 0.9, "ramp_power": 9.0}, "ramp_power"),
])
def test_invalid_parameters_rejected(cfg, physics, params, fragment):
    spec = dict(nominal_overlapped(cfg), **params)
    valid, reasons = validate_waveform(build_waveform(spec, physics), CONSTRAINT)
    assert not valid
    assert any(fragment in r for r in reasons)


# ---------------------------------------------------------------- 4. 解析导数 vs 有限差分
@pytest.mark.parametrize("attr,delta", [("center_velocity", 1e-9), ("depth_rate", 1e-9)])
@pytest.mark.parametrize("spec_name", ["overlapped", "sequential_600us", "sequential_400us"])
def test_analytic_derivatives(cfg, physics, attr, delta, spec_name):
    if spec_name == "overlapped":
        spec = nominal_overlapped(cfg)
    else:
        spec = baseline_specs(cfg)[spec_name]
    waveform = build_waveform(spec, physics)
    times = np.linspace(0.02, 0.98, 41) * waveform.duration_s
    function = getattr(waveform, attr)
    base_function = waveform.center if attr == "center_velocity" else waveform.depth
    analytic = np.asarray(function(times), dtype=float)
    numeric = (np.asarray(base_function(times + delta), dtype=float)
               - np.asarray(base_function(times - delta), dtype=float)) / (2 * delta)
    mask = np.abs(analytic) > 1e-15
    assert np.allclose(analytic[mask], numeric[mask], rtol=5e-3), attr


def test_second_derivative_matches_finite_difference(cfg, physics):
    waveform = build_waveform(nominal_overlapped(cfg), physics)
    times = np.linspace(0.05, 0.95, 31) * waveform.duration_s
    delta = 1e-8
    center = lambda t: np.asarray(waveform.center(t), dtype=float)
    numeric = (center(times + delta) - 2 * center(times) + center(times - delta)) / delta**2
    analytic = np.asarray(waveform.center_acceleration(times), dtype=float)
    mask = np.abs(analytic) > 1e-9
    assert np.allclose(analytic[mask], numeric[mask], rtol=5e-2)


def test_third_derivative_matches_finite_difference(cfg, physics):
    waveform = build_waveform(nominal_overlapped(cfg), physics)
    # 取移动窗口内部（jerk 非零区域）采样，避免端点零区放大相对误差
    times = np.linspace(0.15, 0.55, 21) * waveform.duration_s
    delta = 1e-7
    accel = lambda t: np.asarray(waveform.center_acceleration(t), dtype=float)
    numeric = (accel(times + delta) - accel(times - delta)) / (2 * delta)
    analytic = np.asarray(waveform.center_jerk(times), dtype=float)
    mask = np.abs(analytic) > 1e-6
    assert mask.any()
    assert np.allclose(analytic[mask], numeric[mask], rtol=1e-1)


def test_ramp_power_just_below_one_singular_rejected(cfg, physics):
    spec = dict(nominal_overlapped(cfg), ramp_power=0.95)
    valid, reasons = validate_waveform(build_waveform(spec, physics), CONSTRAINT)
    assert not valid
    assert any("奇异" in reason for reason in reasons)


# ---------------------------------------------------------------- 5. 与 Level 1 一致
def test_sequential_600us_matches_level1(cfg, physics):
    from level1_transfer_1d.config import load_config as load_level1
    from level1_transfer_1d.simulation import propagate as level1_propagate
    level1_cfg = load_level1(PROJECT_ROOT / "configs" / "transfer_1d.yaml")
    t, x, v, _ = level1_propagate(level1_cfg, level1_cfg.aod.initial_center_m, 0.0)
    spec = baseline_specs(cfg)["sequential_600us"]
    states = {"x0_m": np.array([level1_cfg.aod.initial_center_m]),
              "v0_m_per_s": np.array([0.0]),
              "initial_aod_energy_J": np.array([-level1_cfg.aod.depth_j]),
              "initial_excitation_J": np.array([0.0])}
    result = evaluate_waveform_on_states(spec, states, physics, level1_cfg.trajectory.dt_s, 0.0,
                                         keep_trajectory_ids=(0,))
    end = result["records"][0]["end_time_index"]
    trajectory = result["trajectories"][0]
    assert trajectory["position_m"][end] == pytest.approx(x[-1], rel=1e-13, abs=1e-18)
    assert trajectory["velocity_m_per_s"][end] == pytest.approx(v[-1], rel=1e-13, abs=1e-15)


# ---------------------------------------------------------------- 6. 静态极限
def test_static_limit(cfg, physics):
    assert _check_static_limit(physics, cfg.integration.validation_dt_s)["passed"]


# ---------------------------------------------------------------- 7. 功-能平衡与收敛
def test_work_energy_balance_and_dt_convergence(cfg, physics):
    spec = nominal_overlapped(cfg)
    coarse = ideal_run(spec, physics, cfg.integration.validation_dt_s, 1e-5)
    fine = ideal_run(spec, physics, cfg.integration.convergence_dt_s, 1e-5)
    depth_scale = max(physics["slm_depth_j"], physics["aod_initial_depth_j"])
    assert coarse["work_energy"]["max_abs_residual_J"] / depth_scale < 5e-4
    assert (fine["work_energy"]["max_abs_residual_J"]
            < coarse["work_energy"]["max_abs_residual_J"] * 0.6)


def test_convergence_order_helper():
    assert _convergence_order([4e-6, 1e-6]) == pytest.approx(2.0)


# ---------------------------------------------------------------- 8. 保持段
def test_hold_segment_bounded(cfg, physics):
    run = ideal_run(nominal_overlapped(cfg), physics, cfg.integration.validation_dt_s,
                    cfg.integration.post_transfer_hold_s)
    error = run["records"][0]["hold_peak_to_peak_energy_error_over_depth"]
    assert error is not None and error < 1e-4


# ---------------------------------------------------------------- 9. CRN 初态逐行相同
def test_common_random_numbers(cfg, physics):
    states = sample_initial_states(cfg, cfg.initial_ensemble.optimization_seed, 24)
    spec_a, spec_b = nominal_overlapped(cfg), dict(nominal_overlapped(cfg), ramp_power=3.0)
    result_a = evaluate_waveform_on_states(spec_a, states, physics, 2e-8, 1e-5)
    result_b = evaluate_waveform_on_states(spec_b, states, physics, 2e-8, 1e-5)
    assert [r["x0_um"] for r in result_a["records"]] == [r["x0_um"] for r in result_b["records"]]
    assert [r["v0_m_s"] for r in result_a["records"]] == [r["v0_m_s"] for r in result_b["records"]]


# ---------------------------------------------------------------- 10. 种子互不相同
def test_seeds_must_differ(tmp_path):
    data = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    data["initial_ensemble"]["selection_seed"] = data["initial_ensemble"]["optimization_seed"]
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ValueError, match="种子"):
        load_config(path)


# ---------------------------------------------------------------- 11. validation 池隔离
def test_optimizer_never_touches_validation_pool():
    source = (PROJECT_ROOT / "src" / "level2_joint_transfer" / "optimizer.py").read_text(encoding="utf-8")
    assert "validation_seed" not in source
    assert "validation_shots" not in source


class _SequentialPool:
    """单进程替代 multiprocessing.Pool：执行 initializer 后顺序 map。"""

    def __init__(self, *args, initializer=None, initargs=(), **kwargs):
        if initializer is not None:
            initializer(*initargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False

    def map(self, function, jobs):
        return [function(job) for job in jobs]


def test_stage_optimize_never_samples_validation_pool(tmp_path, monkeypatch):
    """接口分层 + spy：优化全流程只允许抽取 optimization/selection 池。"""
    from level2_joint_transfer import cli as level2_cli
    from level2_joint_transfer import optimizer as optimizer_module
    requested_seeds = []
    real_sampler = level2_cli.sample_initial_states

    def spy_sampler(config, seed, shots):
        requested_seeds.append(seed)
        return real_sampler(config, seed, shots)

    monkeypatch.setattr(level2_cli, "sample_initial_states", spy_sampler)
    monkeypatch.setattr(optimizer_module, "Pool", _SequentialPool)
    tiny_cfg = load_config(_tiny_config(tmp_path))
    output = tmp_path / "outputs" / "level2_joint_transfer" / "tiny_test"
    output.mkdir(parents=True, exist_ok=True)
    level2_cli.stage_optimize(tiny_cfg, output)
    assert requested_seeds == [tiny_cfg.initial_ensemble.optimization_seed,
                               tiny_cfg.initial_ensemble.selection_seed]


# ---------------------------------------------------------------- 12. Wilson 区间
def test_wilson_interval_known_case():
    summary = summarize_capture([True] * 3 + [False] * 7)
    assert summary["captured"] == 3 and summary["shots"] == 10
    assert summary["wilson_low"] == pytest.approx(0.1078, abs=5e-4)
    assert summary["wilson_high"] == pytest.approx(0.6032, abs=5e-4)


# ---------------------------------------------------------------- 13. 配对统计
def test_paired_statistics_constructed():
    a = np.array([1, 1, 1, 0, 1, 1, 0, 1])
    b = np.array([1, 0, 1, 0, 1, 0, 0, 1])
    counts = paired_counts(a.astype(bool), b.astype(bool))
    assert counts["both_captured"] == 4 and counts["only_a_captured"] == 2
    assert counts["only_b_captured"] == 0 and counts["both_failed"] == 2
    boot = paired_bootstrap_difference(a.astype(bool), b.astype(bool), 4000, 7)
    assert boot["difference"] == pytest.approx(2 / 8)
    assert boot["ci_low"] <= boot["difference"] <= boot["ci_high"]


# ---------------------------------------------------------------- 14-15. 判据与近阈值样本
def test_capture_criterion_strict():
    assert is_captured(-1e-30)
    assert not is_captured(0.0)
    assert not is_captured(1e-30)


def test_near_threshold_samples_kept(cfg, physics):
    states = sample_initial_states(cfg, cfg.initial_ensemble.optimization_seed, 32)
    result = evaluate_waveform_on_states(nominal_overlapped(cfg), states, physics,
                                         cfg.integration.validation_dt_s, 1e-5)
    assert len(result["records"]) == 32
    for record in result["records"]:
        assert record["captured"] == is_captured(record["final_slm_energy_uK"] * 1.380649e-29 * 1e6)


# ---------------------------------------------------------------- 16. 步长减半收敛
def test_dt_halving_continuous_convergence(cfg, physics):
    states = sample_initial_states(cfg, cfg.initial_ensemble.optimization_seed, 8)
    spec = nominal_overlapped(cfg)
    energies = {}
    for label, dt in (("coarse", cfg.integration.validation_dt_s), ("fine", cfg.integration.convergence_dt_s)):
        evaluation = evaluate_waveform_on_states(spec, states, physics, dt, 1e-5)
        energies[label] = np.array([r["final_slm_energy_uK"] for r in evaluation["records"]])
    assert np.max(np.abs(energies["coarse"] - energies["fine"])) < 1.0


# ---------------------------------------------------------------- 17. 失败候选被记录
def test_failed_candidate_recorded(cfg, physics):
    states = sample_initial_states(cfg, cfg.initial_ensemble.optimization_seed, 4)
    bad_spline = {"type": "spline", "duration_s": cfg.waveforms.optimized_duration_s,
                  "control_s": [0.0, 0.3, 0.6, 1.0],
                  "center_values_m": [2.4e-6, 3e-6, 1e-6, 0.0],  # 非单调
                  "depth_values_j": [1e-27, 8e-28, 3e-28, 0.0]}
    row = _evaluate_job(("bad_0", bad_spline, states, physics, 1e-7, 1e-5, CONSTRAINT))
    assert row["status"] in ("invalid", "failed")
    assert row["failure_reason"]


# ---------------------------------------------------------------- 18-19. CLI 冒烟 + 输出 schema
TINY_OVERRIDE = {
    "initial_ensemble": {"optimization_shots": 6, "selection_shots": 12, "validation_shots": 16,
                         "optimization_seed": 31001, "selection_seed": 31002,
                         "validation_seed": 31003},
    "integration": {"post_transfer_hold_us": 10.0, "optimization_dt_us": 0.02,
                    "validation_dt_us": 0.01, "convergence_dt_us": 0.005},
    "waveforms": {"sequential_baseline_duration_us": 60.0, "compressed_baseline_duration_us": 40.0,
                  "optimized_duration_us": 40.0},
    "optimization": {"max_candidates": 3, "finalists": 1, "refine_steps_per_finalist": 1,
                     "enable_monotone_spline_refinement": True, "spline_candidates": 2},
    "robustness": {"temperatures_uK": [4.0, 6.0], "final_alignment_offsets_um": [0.0, 0.05],
                   "shots_per_condition": 4, "timestep_subset_shots": 4,
                   "bootstrap_samples": 200},
}


def _tiny_config(tmp_path):
    data = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    for section, values in TINY_OVERRIDE.items():
        data[section].update(values)
    data["output"]["directory"] = "outputs/level2_joint_transfer/tiny_test"
    path = tmp_path / "tiny.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


@pytest.mark.slow
def test_cli_stages_smoke(tmp_path):
    from level2_joint_transfer import cli as level2_cli
    config_path = _tiny_config(tmp_path)
    for stage in ("optimize", "validate", "robustness"):
        assert level2_cli.main(["--config", str(config_path), "--stage", stage]) == 0
    output = tmp_path / "outputs" / "level2_joint_transfer" / "tiny_test"
    required = ["config_used.yaml", "candidates.csv", "best_waveform.yaml", "waveforms.csv",
                "validation_shots.csv", "robustness.csv", "metrics.json",
                "image_validation.json", "summary.md"]
    for name in required:
        assert (output / name).is_file(), name
    with (output / "validation_shots.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        assert {"waveform", "shot_id", "captured", "near_threshold", "capture_margin"} <= set(rows[0])
        shot_sets = {}
        for row in rows:
            shot_sets.setdefault(row["waveform"], []).append(row["shot_id"])
        assert shot_sets["sequential_600us"] == shot_sets["overlapped_optimized_400us"]

    def reject_constant(value):
        raise ValueError(f"non-strict JSON constant {value}")

    metrics = json.loads((output / "metrics.json").read_text(), parse_constant=reject_constant)
    assert set(metrics["seeds"]) == {"optimization", "selection", "validation"}
    image_validation = json.loads((output / "image_validation.json").read_text(),
                                  parse_constant=reject_constant)
    assert image_validation["all_passed"]
    frozen = yaml.safe_load((output / "best_waveform.yaml").read_text(encoding="utf-8"))
    assert frozen_spec(frozen, load_config(config_path))["duration_s"] == pytest.approx(40e-6)
    with (output / "candidates.csv").open(newline="", encoding="utf-8") as handle:
        candidate_ids = [row["candidate_id"] for row in csv.DictReader(handle)]
    assert any(candidate_id.startswith("spl_") for candidate_id in candidate_ids)


@pytest.mark.slow
def test_cli_preflight_and_all_stages(tmp_path):
    """preflight 与 all 两个 stage 的冒烟测试。"""
    from level2_joint_transfer import cli as level2_cli
    config_path = _tiny_config(tmp_path)
    assert level2_cli.main(["--config", str(config_path), "--stage", "preflight"]) == 0
    output = tmp_path / "outputs" / "level2_joint_transfer" / "tiny_test"
    preflight = json.loads((output / "preflight.json").read_text())
    assert preflight["all_passed"]
    assert level2_cli.main(["--config", str(config_path), "--stage", "all"]) == 0


@pytest.mark.slow
def test_optimize_checkpoint_resume(tmp_path):
    """检查点续跑：重跑 optimize 应复用三个阶段的既有结果且候选表不变。"""
    from level2_joint_transfer import cli as level2_cli
    config_path = _tiny_config(tmp_path)
    assert level2_cli.main(["--config", str(config_path), "--stage", "optimize"]) == 0
    output = tmp_path / "outputs" / "level2_joint_transfer" / "tiny_test"
    checkpoint = json.loads((output / "optimization_checkpoint.json").read_text())

    def read_candidates():
        with (output / "candidates.csv").open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    first = read_candidates()
    assert {"stage1", "stage2", "spline"} <= set(checkpoint["stages"])
    assert level2_cli.main(["--config", str(config_path), "--stage", "optimize"]) == 0
    second = read_candidates()
    assert first == second
    checkpoint2 = json.loads((output / "optimization_checkpoint.json").read_text())
    for stage in ("stage1", "stage2", "spline"):
        assert checkpoint["stages"][stage]["rows"] == checkpoint2["stages"][stage]["rows"]


# ---------------------------------------------------------------- 20. sources/ 拒绝
def test_output_into_sources_rejected(tmp_path):
    from level0_static_trap.io_utils import ensure_output_directory
    with pytest.raises(ValueError, match="sources"):
        ensure_output_directory(tmp_path / "sources" / "out", tmp_path)


def test_sampler_matches_level1_statistics(cfg):
    states = sample_initial_states(cfg, 123, 2000)
    assert states["acceptance_rate"] > 0.9
    assert np.all(states["initial_aod_energy_J"] < 0)
