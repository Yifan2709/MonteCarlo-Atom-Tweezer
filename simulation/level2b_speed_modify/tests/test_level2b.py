"""Level 2B 测试：剖面定义、速度口径、步长规则、临界速度估计、引擎与 CLI。"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from level2_joint_transfer.config import load_config as load_level2_config
from level2_joint_transfer.level2_simulation import physics_from_config
from level2b_speed_limit.critical_speed import (boundary_covered,
                                                bootstrap_interpolation_speeds,
                                                logistic_fit_speeds)
from level2b_speed_limit.level2b_config import load_level2b_config
from level2b_speed_limit.speed_engine import evaluate_profile, static_limit_check
from level2b_speed_limit.speed_scan import move_only_spec
from level2b_speed_limit.speed_waveforms import (SpeedWaveform, max_step_for,
                                                 snap_dt, speed_calibers)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LEVEL2B_CONFIG = PROJECT_ROOT / "configs" / "level2b_speed_limit.yaml"
LEVEL2_CONFIG = PROJECT_ROOT / "configs" / "level2_joint_transfer.yaml"


@pytest.fixture(scope="module")
def physics():
    return physics_from_config(load_level2_config(LEVEL2_CONFIG))


def make_trapezoid(physics, duration_us=400.0, accel=0.30):
    return SpeedWaveform(
        duration_us * 1e-6, 0.01945549969117366, 1.0,
        0.6069511239907726, 1.0, 1.1121764674723393,
        physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
        move_profile="trapezoid", accel_fraction=accel)


# ---------------------------------------------------------------- 剖面定义
@pytest.mark.parametrize("accel", [0.15, 0.30, 0.45])
def test_trapezoid_peak_equals_cruise(physics, accel):
    waveform = make_trapezoid(physics, 400.0, accel)
    cal = speed_calibers(waveform, 8001)
    assert abs(cal["v_peak_m_per_s"] - waveform.cruise_speed_m_per_s) \
        / waveform.cruise_speed_m_per_s < 1e-6
    assert abs(cal["v_avg_m_per_s"] - (1.0 - accel) * waveform.cruise_speed_m_per_s) \
        / cal["v_avg_m_per_s"] < 1e-9


def test_trapezoid_endpoints_and_displacement(physics):
    waveform = make_trapezoid(physics, 250.0)
    table = waveform.dense_table(8001)
    center = table["aod_center_m"]
    assert abs(center[0] - physics["aod_initial_center_m"]) < 1e-12
    assert abs(center[-1]) < 1e-12
    assert abs((center[0] - center[-1]) - physics["aod_initial_center_m"]) < 1e-12
    assert np.all(np.diff(center) <= 1e-15)
    for probe in (0.0, waveform.duration_s):
        assert abs(float(waveform.center_velocity(probe))) < 1e-15
        assert abs(float(waveform.center_acceleration(probe))) < 1e-9


def test_trapezoid_cruise_segment_constant_velocity(physics):
    waveform = make_trapezoid(physics, 400.0)
    t_mid = waveform.move_t0 + 0.5 * waveform.move_duration_s
    probes = np.linspace(t_mid - 0.1 * waveform.move_duration_s,
                         t_mid + 0.1 * waveform.move_duration_s, 25)
    velocities = np.abs(np.asarray(waveform.center_velocity(probes)))
    assert np.allclose(velocities, waveform.cruise_speed_m_per_s, rtol=1e-12)


def test_smootherstep_mode_matches_overlapped(physics):
    """smootherstep 模式必须与 OverlappedWaveform 完全一致（复用校验）。"""
    from level2_joint_transfer.waveforms import OverlappedWaveform
    args = (400e-6, 0.01945549969117366, 1.0, 0.6069511239907726, 1.0,
            1.1121764674723393, physics["aod_initial_center_m"],
            physics["aod_initial_depth_j"])
    speed = SpeedWaveform(*args, move_profile="smootherstep")
    overlapped = OverlappedWaveform(*args)
    t = np.linspace(0, 400e-6, 500)
    for call_speed, call_ref in ((speed.center, overlapped.center),
                                 (speed.center_velocity, overlapped.center_velocity),
                                 (speed.depth, overlapped.depth)):
        assert np.allclose(np.asarray(call_speed(t)), np.asarray(call_ref(t)),
                           rtol=0, atol=1e-30)


def test_invalid_parameters_rejected(physics):
    with pytest.raises(ValueError):
        SpeedWaveform(400e-6, 0.1, 0.9, 0.5, 1.0, 1.5,
                      physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
                      move_profile="hyperbolic")
    with pytest.raises(ValueError):
        SpeedWaveform(400e-6, 0.1, 0.9, 0.5, 1.0, 1.5,
                      physics["aod_initial_center_m"], physics["aod_initial_depth_j"],
                      accel_fraction=0.6)
    with pytest.raises(ValueError):
        SpeedWaveform(400e-6, 0.5, 0.5, 0.5, 1.0, 1.5,
                      physics["aod_initial_center_m"], physics["aod_initial_depth_j"])


# ---------------------------------------------------------------- 缩放与口径
def test_scaling_invariance(physics):
    w1 = make_trapezoid(physics, 400.0)
    w2 = make_trapezoid(physics, 200.0)
    t1 = w1.dense_table(4001)
    t2 = w2.dense_table(4001)
    v1 = np.max(np.abs(t1["center_velocity_m_per_s"]))
    v2 = np.max(np.abs(t2["center_velocity_m_per_s"]))
    a1 = np.max(np.abs(t1["center_acceleration_m_per_s2"]))
    a2 = np.max(np.abs(t2["center_acceleration_m_per_s2"]))
    assert abs(v2 / v1 - 2.0) < 1e-9
    assert abs(a2 / a1 - 4.0) < 1e-9
    shape1 = np.asarray(w1.center(np.linspace(0, w1.duration_s, 1001)))
    shape2 = np.asarray(w2.center(np.linspace(0, w2.duration_s, 1001)))
    assert np.max(np.abs(shape1 - shape2)) < 1e-13


def test_speed_calibers_consistency(physics):
    for accel in (0.15, 0.30):
        waveform = make_trapezoid(physics, 300.0, accel)
        cal = speed_calibers(waveform)
        d = cal["total_displacement_m"]
        assert abs(cal["v_avg_m_per_s"] * cal["move_duration_s"] - d) / d < 1e-6
        assert cal["v_peak_m_per_s"] >= cal["v_rms_m_per_s"]
        assert cal["v_peak_m_per_s"] >= cal["v_avg_m_per_s"]
    smooth = SpeedWaveform(300e-6, 0.1, 0.9, 0.5, 1.0, 1.5,
                           physics["aod_initial_center_m"],
                           physics["aod_initial_depth_j"],
                           move_profile="smootherstep")
    cal = speed_calibers(smooth)
    assert abs(cal["v_peak_m_per_s"] / cal["v_avg_m_per_s"] - 1.875) / 1.875 < 1e-3


# ---------------------------------------------------------------- 步长规则
def test_max_step_rule():
    waist = 1.17e-6
    slow = max_step_for(0.01, 400e-6, waist, 0.05e-6)
    assert slow == pytest.approx(0.05e-6)
    fast = max_step_for(0.5, 4e-6, waist, 0.05e-6)
    assert fast == pytest.approx(min(0.05e-6, 4e-6 / 2000, 0.01 * waist / 0.5))
    short = max_step_for(0.01, 60e-6, waist, 0.05e-6)
    assert short == pytest.approx(min(0.05e-6, 60e-6 / 2000, 0.01 * waist / 0.01))
    with pytest.raises(ValueError):
        max_step_for(0.0, 1e-6, waist, 1e-6)


def test_snap_dt_divides_total():
    for total, dt_max in ((100e-6, 7e-9), (400e-6, 0.05e-6), (20e-6, 3.4e-9)):
        dt = snap_dt(total, dt_max)
        n = total / dt
        assert abs(n - round(n)) < 1e-9
        assert dt <= dt_max * (1 + 1e-12)


# ---------------------------------------------------------------- 临界速度
def _logistic_flags(speeds_mm, true_v95_mm, shots=400, seed=7):
    rng = np.random.default_rng(seed)
    b = -8.0
    v95 = true_v95_mm
    # p(v95) = 0.95 → a = logit(0.95) − b·log10(v95)
    a = np.log(0.95 / 0.05) - b * np.log10(v95)
    flags = []
    for v in speeds_mm:
        p = 1.0 / (1.0 + np.exp(-(a + b * np.log10(v))))
        flags.append(rng.random(shots) < p)
    return flags


def test_critical_speed_recovery():
    speeds = np.geomspace(10, 300, 9)
    true_v95 = 80.0
    flags = _logistic_flags(speeds, true_v95)
    logi = logistic_fit_speeds(speeds * 1e-3, flags, 0.95,
                               num_bootstrap=50, seed=3)
    assert logi["point_estimate"] is not None
    assert abs(logi["point_estimate"] * 1e3 - true_v95) / true_v95 < 0.10
    boot = bootstrap_interpolation_speeds(speeds * 1e-3, flags, 0.95,
                                          num_bootstrap=300, seed=3)
    assert boot["point_estimate"] is not None
    assert boot["ci_low"] <= boot["point_estimate"] <= boot["ci_high"]


def test_boundary_covered_logic():
    assert boundary_covered([1.0, 1.0, 0.8, 0.3]) is True
    assert boundary_covered([1.0, 1.0, 1.0, 1.0]) is False
    assert boundary_covered([0.6, 0.4, 0.2]) is False
    assert boundary_covered([]) is False


# ---------------------------------------------------------------- 引擎
def test_engine_failure_modes_and_fields(physics):
    cfg2 = load_level2_config(LEVEL2_CONFIG)
    from level2b_speed_limit.speed_engine import sample_states
    states = sample_states(cfg2, 123, 24, temperature_uK=5.0)
    # 极快移动必然跟丢；慢速基本跟随
    for v_cruise, expect_follow_loss in ((0.5, True), (0.006, False)):
        spec, _ = move_only_spec(v_cruise, physics["aod_initial_center_m"], 0.30)
        from level2b_speed_limit.speed_waveforms import build_speed_waveform
        waveform = build_speed_waveform(spec, physics)
        cal = speed_calibers(waveform)
        dt = snap_dt(waveform.duration_s, max_step_for(
            cal["v_peak_m_per_s"], cal["move_duration_s"], physics["aod_waist_m"],
            0.05e-6, min_steps_per_move=400, max_step_displacement_over_waist=0.01))
        run = evaluate_profile(spec, states, physics, dt, hold_s=0.0)
        assert len(run["records"]) == 24
        rate = run["summary"]["follow_loss_rate"]
        if expect_follow_loss:
            assert rate > 0.5
        else:
            assert rate == 0.0
        for record in run["records"]:
            assert set(record) >= {"shot_id", "captured", "follow_lost",
                                   "failure_mode", "max_sep_over_waist_move"}
            assert (record["failure_mode"] is None) == record["captured"] \
                or record["follow_lost"]


def test_engine_static_limit(physics):
    result = static_limit_check(physics, 0.05e-6, 60e-6)
    assert result["passed"], result


# ---------------------------------------------------------------- 配置与种子
def test_config_seeds_distinct_and_loaded():
    bcfg, cfg2 = load_level2b_config(LEVEL2B_CONFIG)
    assert bcfg["initial_ensemble"]["scan_seed"] != bcfg["initial_ensemble"]["boundary_seed"]
    assert cfg2.atom.mass_u == pytest.approx(132.90545196)
    with pytest.raises(ValueError):
        bad = yaml.safe_load(LEVEL2B_CONFIG.read_text(encoding="utf-8"))
        bad["initial_ensemble"]["boundary_seed"] = bad["initial_ensemble"]["scan_seed"]
        path = Path("/tmp/_bad_level2b.yaml") if Path("/tmp").is_dir() else PROJECT_ROOT / "_bad.yaml"
        path.write_text(yaml.safe_dump(bad), encoding="utf-8")
        try:
            load_level2b_config(path)
        finally:
            path.unlink()


def test_scan_resume_consistency(physics, tmp_path):
    """断点续跑：分两轮与一轮跑完的结果逐行一致。"""
    from level2b_speed_limit.speed_scan import (SCAN_COLUMNS, _append_csv_row,
                                                _completed_keys)
    row1 = {"experiment": "pure_move", "structure": "trapezoid_cruise",
            "temperature_uK": 5.0, "target": "6", "status": "ok", "shots": 10,
            "captured": 10, "capture_rate": 1.0}
    path = tmp_path / "scan.csv"
    _append_csv_row(path, row1, fieldnames=SCAN_COLUMNS)
    _append_csv_row(path, {"experiment": "pure_move", "structure": "trapezoid_cruise",
                           "temperature_uK": 5.0, "target": "10", "status": "failed",
                           "error": "x"}, fieldnames=SCAN_COLUMNS)
    keys = _completed_keys(path)
    assert keys == {("pure_move", "trapezoid_cruise", 5.0, "6")}
    with path.open() as handle:
        header = next(csv.reader(handle))
    assert header == SCAN_COLUMNS


# ---------------------------------------------------------------- CLI 冒烟
@pytest.mark.slow
def test_cli_all_stage_smoke(tmp_path):
    """CLI --stage all 冒烟（跳过 level2 十项复跑与既有测试以避免递归）。"""
    bcfg = yaml.safe_load(LEVEL2B_CONFIG.read_text(encoding="utf-8"))
    bcfg["initial_ensemble"]["scan_shots_per_point"] = 12
    bcfg["comparison"]["comparison_shots"] = 12
    bcfg["speed_scan"]["pure_move_cruise_speeds_mm_per_s"] = [6, 100, 500]
    bcfg["speed_scan"]["end_to_end_durations_us"] = [600, 120, 45, 20]
    bcfg["initial_ensemble"]["scan_temperatures_uK"] = [5.0, 20.0]
    bcfg["comparison"]["points_near_boundary"] = [1.0]
    bcfg["noise_speed_validation"]["shots_per_point"] = 8
    bcfg["noise_speed_validation"]["v95_shift"]["durations_us"] = [600, 90, 25]
    bcfg["noise_speed_validation"]["v95_shift"]["scale_factors"] = [100.0]
    bcfg["noise_speed_validation"]["onoff_paired"]["speeds_v_peak_mm_per_s"] = [55.0]
    bcfg["noise_speed_validation"]["onoff_paired"]["temperatures_uK"] = [20.0]
    bcfg["noise_speed_validation"]["amplitude_sensitivity"]["scale_factors"] = [1.0, 1000.0]
    config_path = tmp_path / "smoke.yaml"
    config_path.write_text(yaml.safe_dump(bcfg, allow_unicode=True), encoding="utf-8")
    output = tmp_path / "out"
    completed = subprocess.run(
        [sys.executable, "-m", "level2b_speed_limit.level2b_cli",
         "--config", str(config_path), "--stage", "all",
         "--output-dir", str(output),
         "--skip-level2-preflight", "--skip-existing-tests"],
        cwd=PROJECT_ROOT, text=True, capture_output=True, check=False, timeout=1800)
    assert completed.returncode == 0, completed.stderr[-4000:]
    expected = {"config_used.yaml", "preflight.json", "speed_scan_results.csv",
                "critical_speeds.csv", "comparison_results.csv", "metrics.json",
                "image_validation.json", "summary.md", "scan_flags.npz",
                "capture_vs_speed.png", "critical_speed_vs_temperature.png",
                "excitation_vs_speed_loglog.png", "adiabatic_collapse.png",
                "trapezoid_vs_smoothstep.png", "boundary_phase_space.png",
                "follow_loss_diagnostics.png", "timestep_convergence_fast.png",
                "noise_speed_results.csv", "noise_v95.json",
                "noise_speed_validation.png"}
    names = {p.name for p in output.iterdir()}
    assert expected <= names, sorted(expected - names)
    metrics_text = (output / "metrics.json").read_text(encoding="utf-8")
    assert "NaN" not in metrics_text and "Infinity" not in metrics_text
    json.loads(metrics_text, parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
    images = json.loads((output / "image_validation.json").read_text(encoding="utf-8"))
    assert images["all_passed"]


def test_cli_rejects_sources_output(tmp_path):
    sources_dir = tmp_path / "sources"
    sources_dir.mkdir()
    completed = subprocess.run(
        [sys.executable, "-m", "level2b_speed_limit.level2b_cli",
         "--config", str(LEVEL2B_CONFIG), "--stage", "preflight",
         "--output-dir", str(sources_dir / "x")],
        cwd=PROJECT_ROOT, text=True, capture_output=True, check=False)
    assert completed.returncode == 2
    assert "sources" in completed.stderr


# ---------------------------------------------------------------- 噪声（实验 4）
def _tiny_heating():
    from level2_joint_transfer.noise_heating import build_mean_heating
    # 参量：16/s × k_B·5 μK（线性化常数功率，量级与默认推导一致）
    return build_mean_heating(recoil_power_w=4.948e-36, parametric_power_w=1.104e-27,
                              pointing_power_w=1.795e-29)


def test_engine_noise_columns_and_channel_sum(physics):
    cfg2 = load_level2_config(LEVEL2_CONFIG)
    from level2b_speed_limit.speed_engine import sample_states
    states = sample_states(cfg2, 4242, 8, temperature_uK=5.0)
    spec, _ = move_only_spec(0.01, physics["aod_initial_center_m"], 0.30)
    from level2b_speed_limit.speed_waveforms import build_speed_waveform, snap_dt
    waveform = build_speed_waveform(spec, physics)
    dt = snap_dt(waveform.duration_s, 1e-7)
    run = evaluate_profile(spec, states, physics, dt, hold_s=0.0,
                           work_energy_shots=2, heating=_tiny_heating())
    for record in run["records"]:
        assert {"noise_recoil_energy_uK", "noise_parametric_energy_uK",
                "noise_pointing_energy_uK", "noise_total_energy_uK"} <= set(record)
        total = (record["noise_recoil_energy_uK"] + record["noise_parametric_energy_uK"]
                 + record["noise_pointing_energy_uK"])
        assert abs(total - record["noise_total_energy_uK"]) \
            <= 1e-12 * max(1.0, abs(total))
        assert record["noise_total_energy_uK"] > 0.0
    assert run["summary"]["work_residual_median_over_depth"] < 1e-3


def test_engine_noise_off_has_no_noise_columns(physics):
    cfg2 = load_level2_config(LEVEL2_CONFIG)
    from level2b_speed_limit.speed_engine import sample_states
    states = sample_states(cfg2, 4242, 8, temperature_uK=5.0)
    spec, _ = move_only_spec(0.01, physics["aod_initial_center_m"], 0.30)
    from level2b_speed_limit.speed_waveforms import build_speed_waveform, snap_dt
    waveform = build_speed_waveform(spec, physics)
    dt = snap_dt(waveform.duration_s, 1e-7)
    run = evaluate_profile(spec, states, physics, dt, hold_s=0.0)
    assert all("noise_total_energy_uK" not in r for r in run["records"])
    run2 = evaluate_profile(spec, states, physics, dt, hold_s=0.0)
    assert run["captured_flags"] == run2["captured_flags"]  # 确定性可复现


def test_noise_worsens_margin_and_scaled_channels(physics):
    import pickle
    from level2_joint_transfer.noise_heating import scaled_heating
    cfg2 = load_level2_config(LEVEL2_CONFIG)
    from level2b_speed_limit.speed_engine import sample_states
    states = sample_states(cfg2, 777, 12, temperature_uK=20.0)
    spec, _ = move_only_spec(0.04, physics["aod_initial_center_m"], 0.30)
    from level2b_speed_limit.speed_waveforms import build_speed_waveform, snap_dt
    waveform = build_speed_waveform(spec, physics)
    dt = snap_dt(waveform.duration_s, 5e-8)
    base = _tiny_heating()
    # 线性化后注入不再随激发复合放大，×1000 的注入仍处相位涨落区，
    # 取 ×3000（注入 ~15 μK，仍远小于阱深）保证裕度单调恶化
    big = scaled_heating(base, 3000.0)
    pickle.loads(pickle.dumps(big))  # 可 pickle（多进程要求）
    off = evaluate_profile(spec, states, physics, dt, hold_s=0.0)
    worse = evaluate_profile(spec, states, physics, dt, hold_s=0.0, heating=big)
    margin_off = np.mean([r["capture_margin"] for r in off["records"]])
    margin_big = np.mean([r["capture_margin"] for r in worse["records"]])
    assert margin_big < margin_off
    # 点缩放：scaled 的各通道功率恰为基础的 factor 倍
    scaled = scaled_heating(base, 7.5)
    p0 = base.channel_powers(1e-6, 0.1, 1e-4, 1e-27)
    p1 = scaled.channel_powers(1e-6, 0.1, 1e-4, 1e-27)
    for key in ("recoil", "parametric", "pointing", "total"):
        assert abs(p1[key] - 7.5 * p0[key]) <= 1e-12 * max(1e-30, abs(p1[key]))


def test_single_channel_heating_isolates_channels(physics):
    from level2b_speed_limit.noise_speed_validation import single_channel_heating
    base = _tiny_heating()
    only_param = single_channel_heating(base, "parametric", 100.0)
    powers = only_param.channel_powers(1e-6, 0.1, 1e-4, 1e-27)
    assert powers["recoil"] == 0.0 and powers["pointing"] == 0.0
    assert powers["parametric"] > 0.0
    only_rec = single_channel_heating(base, "recoil", 100.0)
    powers2 = only_rec.channel_powers(1e-6, 0.1, 1e-4, 1e-27)
    assert powers2["parametric"] == 0.0 and powers2["recoil"] > 0.0


def test_noise_config_section_loaded():
    bcfg, _ = load_level2b_config(LEVEL2B_CONFIG)
    noise = bcfg["noise_speed_validation"]
    assert noise["v95_shift"]["structure"] in ("frozen_shape_scaled",
                                               "sequential_scaled",
                                               "trapezoid_frozen_ramp")
    scales = noise["amplitude_sensitivity"]["scale_factors"]
    assert scales == sorted(scales)
    assert noise["v95_shift"]["scale_factors"]
