"""Level 2C 测试：pick-up 波形、方向切换、采样阱、往返引擎、拟合与 CLI。

不测试“模拟必须复现论文数值”（物理结论不是单元测试）；测试的是波形定义正确、
判据方向切换正确、池隔离/CRN 无泄漏、引擎记账自洽、拟合器在构造数据上无偏。
"""

from __future__ import annotations

import numpy as np
import pytest

from level2_joint_transfer.level2_simulation import (build_waveform, evaluate_waveform_on_states,
                                                     ideal_run, physics_from_config,
                                                     sample_initial_states)
from level2_joint_transfer.config import load_config
from level2_joint_transfer.waveforms import (MLCubicWaveform, PickupSequentialWaveform,
                                             _constant_jerk_terms, validate_pickup_waveform)
from level2c_pickup_survival.level2c_config import load_level2c_config
from level2c_pickup_survival.ml_optimize import (build_ml_waveform, lhs_ml_candidates,
                                                 manual_anchor_controls)
from level2c_pickup_survival.reference_data import compare_to_reference, load_reference_curves
from level2c_pickup_survival.survival_engine import (DisorderDraws, build_roundtrip_legs,
                                                     draw_disorder, run_roundtrip_detailed,
                                                     run_survival, time_reverse_waveform)
from level2c_pickup_survival.survival_fit import bootstrap_fit, fit_survival_curve, survival_model

KB_UK = 1.380649e-29
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG2 = PROJECT_ROOT / "configs" / "level2_joint_transfer.yaml"
CONFIG2C = PROJECT_ROOT / "configs" / "level2c_pickup_survival.yaml"


@pytest.fixture(scope="module")
def cfg2():
    return load_config(CONFIG2)


@pytest.fixture(scope="module")
def cfg2c():
    return load_level2c_config(CONFIG2C)


@pytest.fixture(scope="module")
def physics(cfg2):
    return physics_from_config(cfg2)


# ---------------------------------------------------------------- 波形

def test_constant_jerk_profile_boundary_and_displacement():
    """constant-jerk 剖面：端点速度/加速度为零，总位移严格 span。"""
    x, v, a, j = _constant_jerk_terms(np.linspace(0, 1, 5001), 2.4e-6, 208e-6)
    assert abs(x[0]) < 1e-30 and abs(x[-1] - 2.4e-6) < 1e-20
    assert abs(v[0]) < 1e-30 and abs(v[-1]) < 1e-25
    assert abs(a[0]) < 1e-30 and abs(a[-1]) < 1e-20
    # 单调性与 jerk 分段符号
    assert np.all(np.diff(x) >= 0)
    mid = len(x) // 2
    assert j[10] > 0 and j[mid] < 0 and j[-10] > 0


def test_pickup_waveform_endpoints_monotone_and_validation():
    """pick-up：c 0→d、D 0→D0 单调不减；端点速度/加速度为零。"""
    wf = PickupSequentialWaveform(400e-6, 0.48, 2.4e-6, 280.0 * KB_UK)
    valid, reasons = validate_pickup_waveform(wf, 2.4e-6, 280.0 * KB_UK)
    assert valid, reasons
    t = np.linspace(0, 400e-6, 4001)
    assert np.all(np.diff(wf.center(t)) >= -1e-25)
    assert np.all(np.diff(wf.depth(t)) >= -1e-35)
    assert wf.center_velocity(0.0) == 0.0 and wf.center_velocity(400e-6) == 0.0


def test_pickup_dropoff_is_exact_time_reversal():
    """direction='dropoff' 必须与 pick-up 严格时间反演对称。"""
    pu = PickupSequentialWaveform(400e-6, 0.48, 2.4e-6, 280.0 * KB_UK)
    do = PickupSequentialWaveform(400e-6, 0.48, 2.4e-6, 280.0 * KB_UK, direction="dropoff")
    t = np.linspace(0, 400e-6, 3001)
    assert np.allclose(do.center(t), pu.center(400e-6 - t), atol=1e-20)
    assert np.allclose(do.depth(t), pu.depth(400e-6 - t), atol=1e-35)
    valid, reasons = validate_pickup_waveform(do, 2.4e-6, 280.0 * KB_UK, direction="dropoff")
    assert valid, reasons


def test_pickup_derivatives_match_finite_difference():
    """解析导数与有限差分一致（jerk 分段常数，跳变点邻域除外）。"""
    wf = PickupSequentialWaveform(400e-6, 0.48, 2.4e-6, 280.0 * KB_UK)
    t = np.linspace(1e-9, 400e-6 - 1e-9, 4001)
    h = 1e-9
    fd_v = (wf.center(t + h) - wf.center(t - h)) / (2 * h)
    assert np.max(np.abs(fd_v - wf.center_velocity(t))) / np.max(np.abs(wf.center_velocity(t))) < 1e-6
    # 深度变化率：避开升深段末端 kink（解析值在该处不连续）
    fd_d = (wf.depth(t + h) - wf.depth(t - h)) / (2 * h)
    an_d = wf.depth_rate(t)
    scale_d = np.max(np.abs(an_d))
    mask_d = (np.abs(an_d) > 1e-3 * scale_d) & (np.abs(t - 192e-6) > 5e-7)
    assert np.max(np.abs(fd_d[mask_d] - an_d[mask_d])) / scale_d < 1e-6
    # 加速度：避开 jerk 跳变点（208µs 处的 1/4 与 3/4 点）
    fd_a = (wf.center_velocity(t + h) - wf.center_velocity(t - h)) / (2 * h)
    jump1 = 192e-6 + 208e-6 * 0.25
    jump2 = 192e-6 + 208e-6 * 0.75
    away = (np.abs(t - jump1) > 1e-6) & (np.abs(t - jump2) > 1e-6) & (t > 192.5e-6)
    assert np.max(np.abs(fd_a[away] - wf.center_acceleration(t)[away])) \
        / np.max(np.abs(wf.center_acceleration(t))) < 1e-6


def test_time_reverse_waveform_ml_cubic():
    """ML 三次样条的时间反演：控制点反转等价于 c(T−t)。"""
    s = np.linspace(0, 1, 14)
    center = np.linspace(0, 2.4e-6, 14) ** 1.1 / (2.4e-6) ** 0.1  # 单调非线性
    depth = np.linspace(0, 280.0 * KB_UK, 14)
    wf = MLCubicWaveform(400e-6, s, center, depth, 2.4e-6, 280.0 * KB_UK)
    rev = time_reverse_waveform(wf)
    t = np.linspace(0, 400e-6, 2001)
    assert np.allclose(rev.center(t), wf.center(400e-6 - t), atol=1e-12 * 2.4e-6)
    assert np.allclose(rev.depth(t), wf.depth(400e-6 - t), atol=1e-12 * 280.0 * KB_UK)


def test_ml_cubic_overshoot_rejected():
    """ML 样条严重过冲被校验拒绝；温和摆动被允许。"""
    s = np.linspace(0, 1, 14)
    bad_center = np.linspace(0, 2.4e-6, 14)
    bad_center[5] = 4.0e-6  # 67% 过冲
    bad = MLCubicWaveform(400e-6, s, bad_center, np.linspace(0, 280.0 * KB_UK, 14),
                          2.4e-6, 280.0 * KB_UK)
    valid, reasons = validate_pickup_waveform(bad, 2.4e-6, 280.0 * KB_UK)
    assert not valid and any("过冲" in r for r in reasons)
    good_center = np.linspace(0, 2.4e-6, 14)
    good_center[5] += 0.02 * 2.4e-6  # 2% 温和摆动
    good = MLCubicWaveform(400e-6, s, good_center, np.linspace(0, 280.0 * KB_UK, 14),
                           2.4e-6, 280.0 * KB_UK)
    valid2, _ = validate_pickup_waveform(good, 2.4e-6, 280.0 * KB_UK)
    assert valid2


def test_ml_cubic_negative_depth_rejected():
    s = np.linspace(0, 1, 14)
    depth = np.linspace(0, 280.0 * KB_UK, 14)
    depth[3] = -0.2 * 280.0 * KB_UK
    bad = MLCubicWaveform(400e-6, s, np.linspace(0, 2.4e-6, 14), depth, 2.4e-6, 280.0 * KB_UK)
    valid, reasons = validate_pickup_waveform(bad, 2.4e-6, 280.0 * KB_UK)
    assert not valid


# ---------------------------------------------------------------- 采样阱与方向切换

def test_slm_well_sampling(cfg2):
    """well='slm' 采样：全部样本束缚于 SLM 阱，分布中心在 SLM 阱心。"""
    states = sample_initial_states(cfg2, 777, 128, well="slm")
    assert states["sampling_well"] == "slm"
    assert np.all(states["initial_aod_energy_J"] < 0)
    assert abs(states["mean_x0_m"]) < 3 * states["std_x0_m"] / np.sqrt(128)
    # 与 AOD 采样（280µK 深阱、中心 2.4µm）明显不同
    states_aod = sample_initial_states(cfg2, 777, 128, well="aod")
    assert states_aod["mean_x0_m"] > 2.0e-6


def test_direction_switch_criterion(cfg2, physics):
    """pick-up 判据为 AOD 单阱束缚；同一初态在 drop-off 判据下含义不同。"""
    spec = {"type": "pickup_sequential", "duration_s": 400e-6, "ramp_fraction": 0.48,
            "pickup_direction": "pickup", "direction": "pickup", "name": "pu"}
    states = {"x0_m": np.array([0.0]), "v0_m_per_s": np.array([0.0]),
              "initial_aod_energy_J": np.array([-physics["slm_depth_j"]]),
              "initial_excitation_J": np.array([0.0])}
    out = evaluate_waveform_on_states(spec, states, physics, 0.05e-6, 100e-6)
    record = out["records"][0]
    assert record["captured"] and record["direction"] == "pickup"
    assert not record["escaped_during_hold"]
    assert record["final_aod_energy_uK"] < 0
    # 理想 pick-up 末态应靠近 AOD 终点中心 2.4 µm
    assert abs(record["final_hold_x_um"] - 2.4) < 0.2


def test_pickup_escaped_atom_fails(cfg2, physics):
    """保持段内逃出 AOD 阱的 shot 必须判失败（能量远超阱深的高速初态）。

    逃逸速度：sqrt(2·D_AOD/m) ≈ 0.19 m/s；取 0.25 m/s 确保飞出。
    """
    spec = {"type": "pickup_sequential", "duration_s": 400e-6, "ramp_fraction": 0.48,
            "pickup_direction": "pickup", "direction": "pickup", "name": "pu"}
    states = {"x0_m": np.array([0.0, 0.0]), "v0_m_per_s": np.array([0.0, 0.25]),
              "initial_aod_energy_J": np.array([-physics["slm_depth_j"], 0.5 * physics["mass_kg"] * 0.25 ** 2]),
              "initial_excitation_J": np.array([0.0, 0.0])}
    out = evaluate_waveform_on_states(spec, states, physics, 0.05e-6, 100e-6)
    assert out["records"][0]["captured"]
    assert not out["records"][1]["captured"]


def test_dropoff_path_unchanged(cfg2, physics):
    """drop-off 路径保持原样：Level 1 顺序波形理想轨迹仍被 SLM 俘获。"""
    spec = {"type": "sequential", "duration_s": 600e-6, "move_fraction": 0.52, "name": "seq600"}
    out = ideal_run(spec, physics, 0.05e-6, 100e-6)
    assert out["records"][0]["captured"]
    assert out["records"][0]["direction"] == "dropoff"


# ---------------------------------------------------------------- 往返引擎

def test_roundtrip_conservative_lossless(cfg2c, physics):
    """保守模型：热 shots 连续 4 次转移零丢失，survival 恒 1。"""
    wf = build_waveform({"type": "pickup_sequential", "duration_s": 400e-6,
                         "ramp_fraction": cfg2c.pickup.ramp_fraction,
                         "pickup_direction": "pickup", "direction": "pickup"}, physics)
    pool = sample_initial_states(cfg2c.base, cfg2c.roundtrip.scan_seed, 24, well="slm")
    legs = build_roundtrip_legs(wf, cfg2c.roundtrip.wait_s, 0.1e-6)
    out = run_survival(pool["x0_m"], pool["v0_m_per_s"], physics["mass_kg"], legs, 0.1e-6, 4,
                       slm_depth_j=physics["slm_depth_j"], slm_waist_m=physics["slm_waist_m"],
                       aod_waist_m=physics["aod_waist_m"])
    assert np.all(out["survival_fraction"] == 1.0)
    assert np.all(out["transfers_survived"] == 4)


def test_time_reversal_returns_to_origin(cfg2c, physics):
    """时间反演预检查：理想 shot 一个往返回到 SLM 阱底附近（远小于热尺度）。

    物理注记：constant-jerk 移动的非绝热残余使原子在 AOD 底附近留下
    ~a_max/ω_AOD 量级的振荡（本工况 ~2 mm/s），等待段把部分位移转成速度，
    因此往返不严格等于恒等映射；正确判据是残余 ≪ 5 µK 热尺度
    （σ_v ≈ 17.7 mm/s，σ_x ≈ 111 nm），取 0.5σ 为容差。
    """
    wf = build_waveform({"type": "pickup_sequential", "duration_s": 400e-6,
                         "ramp_fraction": 0.48, "pickup_direction": "pickup",
                         "direction": "pickup"}, physics)
    legs = build_roundtrip_legs(wf, cfg2c.roundtrip.wait_s, 0.1e-6)
    detail = run_roundtrip_detailed(0.0, 0.0, physics["mass_kg"], legs, 0.1e-6,
                                    physics["slm_depth_j"], physics["slm_waist_m"],
                                    physics["aod_waist_m"], physics["aod_initial_depth_j"])
    import math
    kbt = 5.0 * KB_UK
    omega = math.sqrt(4 * physics["slm_depth_j"] / (physics["mass_kg"] * physics["slm_waist_m"] ** 2))
    sigma_v = math.sqrt(kbt / physics["mass_kg"])
    sigma_x = sigma_v / omega
    assert abs(detail["x_m"][-1]) < 0.5 * sigma_x
    assert abs(detail["v_m_s"][-1]) < 0.5 * sigma_v


def test_survival_engine_loss_bookkeeping(cfg2c, physics):
    """人为注入损耗（超快 100 µs 档 + 高温初态）：丢失记账自洽且 survival 单调不增。"""
    wf = build_waveform({"type": "pickup_sequential", "duration_s": 100e-6,
                         "ramp_fraction": 0.48, "pickup_direction": "pickup",
                         "direction": "pickup"}, physics)
    pool = sample_initial_states(cfg2c.base, 999, 48, temperature_uK=30.0, well="slm")
    legs = build_roundtrip_legs(wf, cfg2c.roundtrip.wait_s, 0.1e-6)
    out = run_survival(pool["x0_m"], pool["v0_m_per_s"], physics["mass_kg"], legs, 0.1e-6, 8,
                       slm_depth_j=physics["slm_depth_j"], slm_waist_m=physics["slm_waist_m"],
                       aod_waist_m=physics["aod_waist_m"])
    frac = out["survival_fraction"]
    assert np.all(np.diff(frac) <= 1e-15)
    lost = (out["transfers_survived"] < 8).sum()
    assert lost == 48 - out["survival_counts"][8]
    for sid in range(48):
        assert len(out["judgement_energy_over_depth"][sid]) == out["transfers_survived"][sid] + (
            1 if out["transfers_survived"][sid] < 8 else 0)


def test_disorder_draws_deterministic_and_disabled_noop(cfg2c):
    """位点差异：同种子同结果；enabled=false 严格无效应。"""
    d1 = draw_disorder(cfg2c.disorder, 16, 555)
    d2 = draw_disorder(cfg2c.disorder, 16, 555)
    assert np.all(d1.alignment_offset_m == d2.alignment_offset_m)
    assert np.all(d1.slm_depth_factor == d2.slm_depth_factor)
    disabled_cfg = type(cfg2c.disorder)(False, *tuple(getattr(cfg2c.disorder, f) for f in
                                                      ("alignment_sigma_um", "aod_depth_rel_sigma",
                                                       "slm_depth_rel_sigma", "aod_waist_rel_sigma",
                                                       "shot_jitter_alignment_sigma_um")))
    d0 = draw_disorder(disabled_cfg, 16, 555)
    assert np.all(d0.slm_depth_factor == 1.0) and np.all(d0.alignment_offset_m == 0.0)


def test_crn_same_pool_across_waveforms(cfg2c):
    """common random numbers：同一池种子生成的初态数组逐行相同。"""
    a = sample_initial_states(cfg2c.base, 4242, 32, well="slm")
    b = sample_initial_states(cfg2c.base, 4242, 32, well="slm")
    assert np.array_equal(a["x0_m"], b["x0_m"]) and np.array_equal(a["v0_m_per_s"], b["v0_m_per_s"])


# ---------------------------------------------------------------- 拟合与参考数据

def test_survival_fit_recovers_known_params():
    """在构造的已知 (p0, p, b) 曲线上回收参数（误差 < 1%）。"""
    n = np.arange(1, 61)
    truth = (0.995, 0.999, 0.05)
    curve = survival_model(n, *truth)
    fit = fit_survival_curve(n, curve)
    assert abs(fit["p0"] - truth[0]) < 0.01
    assert abs(fit["p"] - truth[1]) < 1e-3
    assert abs(fit["b"] - truth[2]) < 0.01
    assert not fit["degenerate_full_survival"]


def test_survival_fit_degenerate_full_survival():
    """全存活曲线：不强行拟合，给 Wilson 下界并标记 degenerate。"""
    n = np.arange(1, 61)
    fit = fit_survival_curve(n, np.ones_like(n, dtype=float), shot_count=256)
    assert fit["degenerate_full_survival"]
    assert fit["p_lower_wilson95"] > 0.99
    assert fit["total_transfer_trials"] == 256 * 60


def test_bootstrap_fit_runs_and_brackets():
    ts = np.full(200, 60)
    rng = np.random.default_rng(0)
    losses = rng.random(200) < 0.2
    ts[losses] = rng.integers(5, 55, size=int(losses.sum()))
    out = bootstrap_fit(ts, 60, 50, 12345)
    assert out["num_bootstrap"] == 50
    assert out["p_median"] is None or 0.9 < out["p_median"] <= 1.0


def test_reference_csv_loads_and_compares(cfg2c):
    curves = load_reference_curves(cfg2c.reference_csv)
    assert set(curves) == {"ml_400us", "manual_600us", "manual_400us", "manual_200us"}
    n_ml, s_ml = curves["ml_400us"]
    assert n_ml[0] == 5 and n_ml[-1] == 60
    assert 0.9 < s_ml[0] < 1.0 and 0.5 < s_ml[-1] < 0.8
    # 200 µs 曲线只到 n=31
    n_200, _ = curves["manual_200us"]
    assert n_200[-1] == 31
    comp = compare_to_reference(n_ml, s_ml, curves["ml_400us"], 0.02)
    assert comp["matches_reference"]  # 与自身对比必然匹配


def test_manual_anchor_controls_endpoints(cfg2c):
    center, depth = manual_anchor_controls(cfg2c, 400.0)
    assert center.size == 14 and depth.size == 14
    assert abs(center[0]) < 1e-20 and abs(center[-1] - 2.4e-6) < 1e-15
    assert abs(depth[0]) < 1e-40 and abs(depth[-1] - 280.0 * KB_UK) / (280.0 * KB_UK) < 1e-9


def test_lhs_candidates_count_and_anchor(cfg2c):
    cands = lhs_ml_candidates(cfg2c)
    assert len(cands) == cfg2c.ml_optimization.candidates
    assert cands[0]["origin"] == "manual_anchor"
    assert all(len(c["center_values_m"]) == 14 for c in cands)


def test_build_ml_waveform_via_spec(physics):
    """build_waveform 支持 ml_cubic 与 pickup_sequential 两类新 spec。"""
    s = np.linspace(0, 1, 14)
    spec = {"type": "ml_cubic", "duration_s": 400e-6, "control_s": s,
            "center_values_m": np.linspace(0, 2.4e-6, 14),
            "depth_values_j": np.linspace(0, 280.0 * KB_UK, 14)}
    wf = build_waveform(spec, physics)
    assert isinstance(wf, MLCubicWaveform)
    assert abs(wf.center(400e-6) - 2.4e-6) < 1e-12


def test_ml_candidates_evaluated_on_own_waveform(cfg2c):
    """回归测试：legs 必须按候选自身波形构建。

    曾按时长缓存 legs 导致全部 LHS 候选被评估在锚点波形上（训练成绩虚高、
    冻结候选在复核池崩塌）。不同合法候选的激发统计必须不同。
    """
    from level2_joint_transfer.noise_heating import heating_from_config
    from level2c_pickup_survival.ml_optimize import evaluate_ml_candidate
    pool = sample_initial_states(cfg2c.base, cfg2c.ml_optimization.candidate_seed, 16, well="slm")
    heating = heating_from_config(cfg2c.base)
    ac, ad = manual_anchor_controls(cfg2c, 400.0)
    r_a = evaluate_ml_candidate(cfg2c, ac, ad, pool, heating)
    wig_c = np.linspace(0, 2.4e-6, 14)
    wig_c[4:9] -= 0.04 * 2.4e-6
    r_b = evaluate_ml_candidate(cfg2c, wig_c, ad, pool, heating)
    assert r_a["status"] == r_b["status"] == "ok"
    assert r_a["mean_final_excitation_alive"] != r_b["mean_final_excitation_alive"]
    assert r_a["max_abs_center_acceleration_m_per_s2"] \
        != r_b["max_abs_center_acceleration_m_per_s2"]


def test_survival_engine_jitter_with_midrun_loss(cfg2c, physics):
    """回归测试：抖动数组尺寸必须跟随活跃池收缩。

    高温 + 超快档制造中途丢失；jitter_sigma>0 时曾因尺寸错位崩溃。
    """
    wf = build_waveform({"type": "pickup_sequential", "duration_s": 100e-6,
                        "ramp_fraction": 0.48, "pickup_direction": "pickup",
                        "direction": "pickup"}, physics)
    pool = sample_initial_states(cfg2c.base, 555, 32, temperature_uK=30.0, well="slm")
    disorder = draw_disorder(cfg2c.disorder, 32, 777)
    legs = build_roundtrip_legs(wf, cfg2c.roundtrip.wait_s, 0.1e-6)
    out = run_survival(pool["x0_m"], pool["v0_m_per_s"], physics["mass_kg"], legs, 0.1e-6, 8,
                       slm_depth_j=physics["slm_depth_j"], slm_waist_m=physics["slm_waist_m"],
                       aod_waist_m=physics["aod_waist_m"],
                       disorder=disorder, jitter_sigma_m=0.01e-6)
    frac = out["survival_fraction"]
    assert np.all(np.diff(frac) <= 1e-15)
    assert out["shots"] == 32


# ---------------------------------------------------------------- 2D 横向模型（第 12 步）

def test_2d_sampler_bound_and_symmetric(cfg2c):
    """2D 采样器：全部样本 E<0，y 方向零均值对称。"""
    from level1_transfer_1d.thermal import sample_thermal_initial_state_2d
    from level2_joint_transfer.config import level1_view
    view = level1_view(cfg2c.base)
    rng = np.random.default_rng(321)
    states = [sample_thermal_initial_state_2d(view, rng, well="slm") for _ in range(200)]
    assert all(s.total_energy_J < 0 for s in states)
    ys = np.array([s.y0_m for s in states])
    assert abs(ys.mean()) < 3 * ys.std(ddof=1) / np.sqrt(200)
    vys = np.array([s.vy0_m_s for s in states])
    assert abs(vys.mean()) < 3 * vys.std(ddof=1) / np.sqrt(200)


def test_2d_roundtrip_conservative_lossless(cfg2c, physics):
    """2D 保守模型：400 µs 档 4 次转移零丢失。"""
    from level1_transfer_1d.thermal import sample_thermal_initial_state_2d
    from level2_joint_transfer.config import level1_view
    from level2c_pickup_survival.survival_engine_2d import run_survival_2d
    wf = build_waveform({"type": "pickup_sequential", "duration_s": 400e-6,
                        "ramp_fraction": 0.48, "pickup_direction": "pickup",
                        "direction": "pickup"}, physics)
    view = level1_view(cfg2c.base)
    rng = np.random.default_rng(654)
    states = [sample_thermal_initial_state_2d(view, rng, well="slm") for _ in range(16)]
    legs = build_roundtrip_legs(wf, cfg2c.roundtrip.wait_s, 0.1e-6)
    out = run_survival_2d(
        np.array([s.x0_m for s in states]), np.array([s.y0_m for s in states]),
        np.array([s.vx0_m_s for s in states]), np.array([s.vy0_m_s for s in states]),
        physics["mass_kg"], legs, 0.1e-6, 4,
        slm_depth_j=physics["slm_depth_j"], slm_waist_m=physics["slm_waist_m"],
        aod_waist_m=physics["aod_waist_m"])
    assert np.all(out["survival_fraction"] == 1.0)
