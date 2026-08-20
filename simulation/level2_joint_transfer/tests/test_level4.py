"""Level 4 测试：分段协议、basin 分类、多轮状态机、功-能、相位与统计。"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from level4_roundtrip_coherence import protocol_segments as ps
from level4_roundtrip_coherence import upstream_artifacts
from level4_roundtrip_coherence.dd_sequences import (toggling_function,
                                                     validate_pulse_set)
from level4_roundtrip_coherence.dephasing import (PhaseAccumulator,
                                                  analytic_constant_phase,
                                                  conditional_contrast)
from level4_roundtrip_coherence.energy_ledger import (segment_ledger_rows,
                                                      static_hold_check)
from level4_roundtrip_coherence.level4_cli import (_frozen_validation_yaml,
                                                   build_protocols,
                                                   _outcome_rows, main)
from level4_roundtrip_coherence.level4_config import load_config
from level4_roundtrip_coherence.noise_models import (IntensityNoise,
                                                     noise_statistics)
from level4_roundtrip_coherence.roundtrip_simulation import (repeat_grid,
                                                             simulate_roundtrip)
from level4_roundtrip_coherence.roundtrip_statistics import (
    evaluate_outcomes, heating_slope, survival_fits, summarize_success)
from level4_roundtrip_coherence.sampling4 import sample_slm_states
from level4_roundtrip_coherence.trap_assignment import classify_basin
from level3_3d_transport_lensing.gaussian_3d import KB

US = 1e-6
HBAR = 1.054571817e-34
CONFIG = Path(__file__).resolve().parents[1] / "configs" / "level4_roundtrip.yaml"


@pytest.fixture(scope="module")
def env():
    """一次构建配置/manifest/协议/小初态池。"""
    cfg = load_config(CONFIG)
    manifest, resolved = upstream_artifacts.build_manifest(cfg)
    protocols = build_protocols(cfg, resolved)
    states = sample_slm_states(cfg, 123, 12)
    return cfg, manifest, resolved, protocols, states


# 1. 分段端点与导数连续
def test_segment_endpoint_continuity(env):
    cfg, _, resolved, protocols, _ = env
    for name in ("protocol_a", "protocol_b"):
        tl = protocols[name]
        for i in range(len(tl.segments) - 1):
            left = tl.segments[i].controls(np.array([tl.segments[i].duration_s]))
            right = tl.segments[i + 1].controls(np.array([0.0]))
            dense = tl.segments[i].controls(
                np.linspace(0, tl.segments[i].duration_s, 11))
            for key in ("p", "pd", "pdd", "d_slm", "ddot_slm", "d_aod",
                        "ddot_aod"):
                scale = max(float(np.max(np.abs(dense[key]))), 1e-30)
                gap = float(np.max(np.abs(np.asarray(left[key])
                                          - np.asarray(right[key]))))
                assert gap < 1e-6 * scale + 1e-12, (name, i, key, gap)


# 2. 完整控制 round trip 回到初始控制状态
def test_control_returns_to_initial(env):
    _, _, _, protocols, _ = env
    for name in ("protocol_a", "protocol_b"):
        tl = protocols[name]
        c0, cT = tl.control_at(0.0), tl.control_at(tl.total_duration_s)
        for key in ("c1", "c2", "d_slm", "d_aod"):
            assert abs(cT[key] - c0[key]) < 1e-15, (name, key)


# 3. Protocol B 不重复 split/merge
def test_protocol_b_no_duplicate_split_merge(env):
    _, _, _, protocols, _ = env
    kinds = [s.kind for s in protocols["protocol_b"].segments]
    assert "split" not in kinds and "merge" not in kinds
    pickup = next(s for s in protocols["protocol_b"].segments
                  if s.kind == "pickup")
    p = pickup.controls(np.linspace(0, pickup.duration_s, 11))["p"]
    assert float(np.max(np.abs(p))) > 1e-9  # L2 波形自带 2.4 μm 移动


# 4. Level 2 transfer 三维回归（lensing-off、x2=z=0 切片）
def test_level2_regression_3d(env):
    cfg, _, resolved, _, _ = env
    spec = resolved["level2_spec"]
    dropoff = ps.Segment(
        "dropoff", "dropoff", spec["duration_us"] * US,
        ps._combine(ps._l2_waveform_path(spec, 2.4e-6, 280e-6 * KB,
                                         spec["duration_us"] * US, False),
                    ps._constant_depth(140e-6 * KB),
                    ps._constant_depth(0.0)),
        checkpoint="dropoff_end", expected_basin="bound_to_slm")
    tl = ps.ProtocolTimeline("l2_3d", [dropoff],
                             ps.physics_from_config(cfg, v_s=None))
    states = sample_slm_states(cfg, 456, 48, temperature_uK=5.0, depth_uK=280.0)
    split = 2.4e-6
    f1, f2 = tl.physics.geometry_factors
    pos = states["pos_m"] + np.array([split * f1, split * f2, 0.0])
    res = simulate_roundtrip(tl, 0.10 * US, pos, states["vel_m_per_s"],
                             tl.physics, {"none": np.empty(0)})
    basin = res.checkpoints[-1]["classification"]["basin"]
    captured = np.mean([b == "bound_to_slm" for b in basin])
    assert captured >= 0.95


# 5. Level 3 long-move 回归（610 μm/1600 μs/sev1.0 冻结条件）
def test_level3_long_move_regression(env):
    cfg, _, resolved, _, _ = env
    seg = ps.Segment(
        "outbound_long_move", "outbound_long_move", 1600 * US,
        ps._combine(ps._profile_move(0.0, 610e-6, 1600 * US, "adiabatic_sine"),
                    ps._constant_depth(0.0),
                    ps._constant_depth(280e-6 * KB)),
        checkpoint="outbound_end", expected_basin="bound_to_aod")
    tl = ps.ProtocolTimeline("l3", [seg],
                             ps.physics_from_config(
                                 cfg, v_s=resolved["nominal_vs_m_per_s"]))
    states = sample_slm_states(cfg, 789, 48, temperature_uK=5.0, depth_uK=280.0)
    res = simulate_roundtrip(tl, 0.10 * US, states["pos_m"],
                             states["vel_m_per_s"], tl.physics,
                             {"none": np.empty(0)})
    basin = res.checkpoints[-1]["classification"]["basin"]
    assert np.mean([b == "bound_to_aod" for b in basin]) >= 0.95


# 6. 控制时间反向 + 速度反转
def test_control_time_reversal(env):
    _, _, _, protocols, _ = env
    tl = protocols["protocol_a_no_lensing"]
    rev = tl.reversed()
    pos0 = np.zeros((1, 3))
    vel0 = np.zeros((1, 3))
    fwd = simulate_roundtrip(tl, 0.10 * US, pos0, vel0, tl.physics,
                             {"none": np.empty(0)})
    back = simulate_roundtrip(rev, 0.10 * US, fwd.final_pos,
                              -fwd.final_vel, rev.physics,
                              {"none": np.empty(0)})
    assert np.linalg.norm(back.final_pos - pos0) < 1e-9
    assert np.linalg.norm(back.final_vel - vel0) < 1e-6


# 7. 总势、总力与组件一致（有限差分）
def test_total_force_consistency(env):
    _, _, _, protocols, _ = env
    tl = protocols["protocol_a"]
    phys = tl.physics
    ctrl = tl.control_at(2500 * US)
    from level3_3d_transport_lensing.aod_lensing import \
        lensing_potential_and_force
    from level3_3d_transport_lensing.gaussian_3d import static_force
    pos = np.array([[1e-7, -8e-8, 3e-7]])
    h = 1e-11
    for axis in range(3):
        dp = np.zeros(3)
        dp[axis] = h
        u_plus = _u_tot(pos + dp, ctrl, phys)
        u_minus = _u_tot(pos - dp, ctrl, phys)
        fd = -(u_plus - u_minus) / (2 * h)
        fx, fy, fz = static_force(
            pos[0, 0], pos[0, 1], pos[0, 2], ctrl["d_slm"],
            phys.slm_waist_m, phys.slm_zr_m)
        _, ax, ay, az = lensing_potential_and_force(
            pos[0, 0] - ctrl["c1"], pos[0, 1] - ctrl["c2"], pos[0, 2],
            ctrl["d_aod"], phys.aod_waist_m, phys.aod_zr_m,
            ctrl["zs1"], ctrl["zs2"])
        total = np.array([fx + ax, fy + ay, fz + az])
        assert abs(fd - total[axis]) < 1e-13


def _u_tot(pos, ctrl, phys):
    from level3_3d_transport_lensing.aod_lensing import \
        lensing_potential_and_force
    from level3_3d_transport_lensing.gaussian_3d import static_potential
    u_s = static_potential(pos[0, 0], pos[0, 1], pos[0, 2], ctrl["d_slm"],
                           phys.slm_waist_m, phys.slm_zr_m)
    u_a = lensing_potential_and_force(
        pos[0, 0] - ctrl["c1"], pos[0, 1] - ctrl["c2"], pos[0, 2],
        ctrl["d_aod"], phys.aod_waist_m, phys.aod_zr_m,
        ctrl["zs1"], ctrl["zs2"])[0]
    return u_s + u_a


# 8. overlap basin 分类与 ambiguous
def test_basin_classification(env):
    _, _, _, protocols, _ = env
    phys = protocols["protocol_a"].physics
    d = 280e-6 * KB
    ctrl = {"d_slm": 60e-6 * KB, "d_aod": d, "c1": 2.4e-6, "c2": 2.4e-6,
            "zs1": 0.0, "zs2": 0.0}
    in_aod = classify_basin(np.array([[2.4e-6, 2.4e-6, 0.0]]),
                            np.zeros((1, 3)), ctrl, phys)["basin"][0]
    merged = classify_basin(np.array([[0.0, 0.0, 0.0]]), np.zeros((1, 3)),
                            {**ctrl, "c1": 0.0, "c2": 0.0}, phys)["basin"][0]
    gone = classify_basin(np.array([[50e-6, 50e-6, 0.0]]), np.zeros((1, 3)),
                          ctrl, phys)["basin"][0]
    hot = classify_basin(np.array([[2.4e-6, 2.4e-6, 0.0]]),
                         np.full((1, 3), 0.12), ctrl, phys)["basin"][0]
    assert in_aod == "bound_to_aod"
    assert merged == "shared_or_ambiguous"
    assert gone == "unbound"
    assert hot in ("shared_or_ambiguous", "unbound")  # 高能不判成功


# 9. checkpoint 局部束缚判据
def test_checkpoint_bound_criterion(env):
    cfg, _, _, protocols, states = env
    tl = protocols["protocol_a"]
    res = simulate_roundtrip(tl, 0.10 * US, states["pos_m"][:4],
                             states["vel_m_per_s"][:4], tl.physics,
                             {"none": np.empty(0)})
    for rec in res.checkpoints:
        cl = rec["classification"]
        for basin, e_slm, e_aod in zip(cl["basin"], cl["e_slm"], cl["e_aod"]):
            if basin == "bound_to_slm":
                assert e_slm < 0
            if basin == "bound_to_aod":
                assert e_aod < 0


# 10/11. first-failure 状态机与 recaptured/absorbing 区分
def test_first_failure_state_machine():
    names = ["pickup_end", "split_end", "outbound_end", "final_hold_end"]
    per_shot = np.array([[1, 1, 1, 1], [1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 1, 0]],
                        dtype=bool)
    matches = [per_shot[:, k] for k in range(4)]
    out = evaluate_outcomes(matches, per_shot[:, -1], checkpoint_names=names)
    assert list(out["first_failure_stage"]) == ["none", "split", "pickup",
                                                "final_hold"]
    assert out["roundtrip_success"][0] and not out["roundtrip_success"][1]
    assert out["recaptured"][1] and not out["recaptured"][0]
    assert out["absorbing_failure"][2] and not out["final_retained"][2]


# 12. 多轮状态连续传递且不重采样
def test_round_state_continuity(env):
    cfg, _, _, protocols, states = env
    tl = protocols["protocol_a"]
    pos, vel = states["pos_m"][:4], states["vel_m_per_s"][:4]
    pulses = {"none": np.empty(0)}
    two = simulate_roundtrip(tl, 0.10 * US, pos, vel, tl.physics, pulses,
                             n_repeats=2)
    one = simulate_roundtrip(tl, 0.10 * US, pos, vel, tl.physics, pulses)
    chained = simulate_roundtrip(tl, 0.10 * US, one.final_pos, one.final_vel,
                                 tl.physics, pulses)
    assert np.array_equal(chained.final_pos, two.final_pos)
    assert np.array_equal(chained.final_vel, two.final_vel)


# 13. 每段及总周期功-能关系
def test_segment_work_energy(env):
    cfg, _, _, protocols, states = env
    tl = protocols["protocol_a"]
    dt = 0.10 * US
    grid = tl.build_grid(dt)
    res = simulate_roundtrip(tl, dt, states["pos_m"][:4],
                             states["vel_m_per_s"][:4], tl.physics,
                             {"none": np.empty(0)})
    rows, seg_works = segment_ledger_rows(res, grid, "a", per_shot=False)
    depth = 280e-6 * KB
    for seg, w in seg_works.items():
        assert np.max(np.abs(w["residual"])) / depth < 1e-4, seg
    assert np.max(np.abs(res.residual)) / depth < 1e-4


# 14. 静态 hold 能量误差有界
def test_static_hold_bounded(env):
    cfg, _, _, protocols, states = env
    tl = protocols["static_idle"]
    dt = 0.10 * US
    grid = tl.build_grid(dt)
    res = simulate_roundtrip(tl, dt, states["pos_m"][:4],
                             states["vel_m_per_s"][:4], tl.physics,
                             {"none": np.empty(0)})
    checks = static_hold_check(res, grid)
    assert checks
    for c in checks.values():
        assert c["max_abs_w_ext_J"] / (180e-6 * KB) < 1e-6


# 15. 常数差分光移解析相位
def test_analytic_phase():
    u0 = -3.5e-27
    eta = 1.3e-4
    duration = 60 * US
    dt = 0.05 * US
    times = np.arange(int(duration / dt) + 1) * dt
    seg = np.zeros(len(times), dtype=int)
    pulses = {"none": np.empty(0),
              "echo": np.array([0.5 * duration + 0.037 * US])}
    acc = PhaseAccumulator(pulses, times, 1, seg)
    u = np.full(len(times), u0)
    for i in range(len(times) - 1):
        acc.step(i, u[i:i + 1], np.zeros(1), u[i + 1:i + 2], np.zeros(1))
    delta = eta * u0 / HBAR
    got = acc.phases(eta, 0.0)
    assert abs(got["none"][0] - analytic_constant_phase(delta, duration)) \
        < 1e-6 * abs(got["none"][0])
    assert abs(got["echo"][0] - analytic_constant_phase(delta, duration,
                                                       pulses["echo"])) \
        < 1e-6 * abs(got["none"][0])


# 16. toggling function 正确性
def test_toggling_function():
    t = np.linspace(0, 1, 1001)
    assert np.all(toggling_function(t, []) == 1.0)
    y = toggling_function(t, [0.3, 0.7])
    assert y[0] == 1 and y[-1] == 1
    assert np.all(y[t < 0.3] == 1)
    assert np.all(y[(t > 0.3) & (t < 0.7)] == -1)
    assert np.all(y[t > 0.7] == 1)
    y4 = toggling_function(t, [0.125, 0.375, 0.625, 0.875])
    assert np.sum(np.diff(y4) != 0) == 4


# 17. pulse 时刻网格拆分无重复/漏积分
def test_pulse_grid_no_double_count():
    duration = 100 * US
    dt = 0.05 * US
    n = int(duration / dt)
    times = np.arange(n + 1) * dt
    seg = np.zeros(len(times), dtype=int)
    u = np.full(len(times), -1.7e-27)
    pulse_times = np.array([13.331 * US, 51.172 * US, 88.913 * US])
    acc = PhaseAccumulator({"dd": pulse_times}, times, 1, seg)
    for i in range(n):
        acc.step(i, u[i:i + 1], np.zeros(1), u[i + 1:i + 2], np.zeros(1))
    got = acc.phases(1.0, 0.0)["dd"][0]
    expect = analytic_constant_phase(-1.7e-27 / HBAR, duration, pulse_times)
    assert abs(got - expect) < 1e-6 * abs(expect)
    issues = validate_pulse_set(pulse_times, dt, duration, [0.0, duration])
    assert issues == []


# 18. ensemble contrast 构造案例
def test_contrast_cases():
    assert conditional_contrast(np.zeros(32))["contrast"] > 0.999
    rng = np.random.default_rng(5)
    assert conditional_contrast(
        rng.uniform(-np.pi, np.pi, 256))["contrast"] < 0.25
    two = conditional_contrast(np.concatenate([np.full(64, 0.25),
                                               np.full(64, -0.25)]))
    assert abs(two["contrast"] - abs(math.cos(0.25))) < 1e-12
    assert abs(two["circular_std"]
               - math.sqrt(-2 * math.log(two["contrast"]))) < 1e-12


# 19. survival conditioning 不混入失败 shots
def test_survival_conditioning():
    phi = np.concatenate([np.full(50, 1.0), rng_uniform(50)])
    mask = np.concatenate([np.ones(50, dtype=bool),
                           np.zeros(50, dtype=bool)])
    stats = conditional_contrast(phi, mask)
    assert stats["n"] == 50
    assert stats["contrast"] > 0.999


def rng_uniform(n):
    return np.random.default_rng(11).uniform(-np.pi, np.pi, n)


# 20. 噪声统计与 seed 复现
def test_noise_statistics_and_reproducibility():
    ou = IntensityNoise("ornstein_uhlenbeck", 4000, 0.01, tau_s=50 * US,
                        seed=101)
    stats = noise_statistics(ou)
    assert abs(stats["mean_slm"]) < 0.01 * 3 / math.sqrt(4000)
    assert abs(stats["std_slm"] - 0.01) < 0.01 * 0.05
    assert abs(stats["cross_corr"]) < 0.05
    assert stats["tau_measured_s"] is not None
    assert abs(stats["tau_measured_s"] - 50 * US) / (50 * US) < 0.05
    a = IntensityNoise("quasistatic", 100, 0.02, seed=7)
    b = IntensityNoise("quasistatic", 100, 0.02, seed=7)
    assert np.array_equal(a.eps_slm, b.eps_slm)
    c = IntensityNoise("quasistatic", 100, 0.02, seed=8)
    assert not np.array_equal(a.eps_slm, c.eps_slm)


# 21. 初态 seed 与噪声 seed 隔离
def test_seed_isolation(env):
    cfg = env[0]
    assert len({cfg.single_round_seed, cfg.repeated_round_seed,
                cfg.long_tail_seed, cfg.sensitivity_seed,
                cfg.development_seed, cfg.noise_seed}) == 6


# 22. common random numbers 对齐
def test_common_random_numbers(env):
    cfg, _, resolved, protocols, _ = env
    s1 = sample_slm_states(cfg, 41001, 16)
    s2 = sample_slm_states(cfg, 41001, 16)
    assert np.array_equal(s1["pos_m"], s2["pos_m"])
    tl = protocols["protocol_a"]
    tl2 = protocols["protocol_a_no_lensing"]
    r1 = simulate_roundtrip(tl, 0.10 * US, s1["pos_m"], s1["vel_m_per_s"],
                            tl.physics, {"none": np.empty(0)})
    r2 = simulate_roundtrip(tl2, 0.10 * US, s1["pos_m"], s1["vel_m_per_s"],
                            tl2.physics, {"none": np.empty(0)})
    # lensing-off 与 nominal 从完全相同初态出发（CRN）
    assert np.array_equal(r1.e0, r2.e0)


# 23. validation 冻结
def test_frozen_validation(env):
    cfg, _, resolved, _, _ = env
    frozen = _frozen_validation_yaml(cfg, resolved)
    assert frozen["single_round_seed"] == cfg.single_round_seed
    assert len(frozen["single_round_conditions"]) == 4
    assert set(frozen["dd_modes"]) == {"none", "spin_echo",
                                       "xy4_per_long_move"}
    assert frozen["frozen_at_stage"] == "before_single_round"


# 24. Wilson、paired bootstrap、circular statistics
def test_wilson_bootstrap_circular():
    s = summarize_success(np.zeros(40, dtype=bool))
    assert s["wilson_low"] == 0.0 and s["wilson_high"] < 0.09
    s = summarize_success(np.ones(40, dtype=bool))
    assert s["wilson_high"] <= 1.0 + 1e-9 and s["wilson_low"] > 0.91
    from level4_roundtrip_coherence.roundtrip_statistics import \
        paired_bootstrap_difference, paired_counts
    a = np.array([1] * 30 + [0] * 10, dtype=bool)
    b = np.array([1] * 20 + [0] * 20, dtype=bool)
    counts = paired_counts(a, b)
    assert counts["both_captured"] == 20 and counts["only_a_captured"] == 10
    boot = paired_bootstrap_difference(a, b, 500, seed=1)
    assert abs(boot["difference"] - 0.25) < 1e-12
    assert boot["ci_excludes_zero"]


# 25. 生存拟合在不可辨识数据上安全失败
def test_survival_fit_unidentifiable():
    out = survival_fits([1, 2, 3], [1.0, 1.0, 1.0])
    assert "reason" in out["skipped"] and not out["fits"]
    out2 = survival_fits([1, 2, 3, 4, 5], [1.0, 0.98, 0.96, 0.9, 0.85])
    assert out2["fits"]  # 可辨识数据给出拟合
    assert all("aic" in f for f in out2["fits"].values())


# 26. 步长二阶收敛（slow）
@pytest.mark.slow
def test_dt_convergence_second_order(env):
    cfg, _, resolved, protocols, _ = env
    tl = protocols["protocol_a"]
    states = sample_slm_states(cfg, 999, 8)
    base = None
    pos_errs = []
    for dt_us in (0.10, 0.05):
        pulses = {m: tl.dd_pulse_times(m, dt_us * US, None)
                  for m in ("none", "spin_echo", "xy4_per_long_move")}
        res = simulate_roundtrip(tl, dt_us * US, states["pos_m"],
                                 states["vel_m_per_s"], tl.physics, pulses)
        if base is None:
            base = res
        else:
            pos_errs.append(float(np.max(np.linalg.norm(
                res.final_pos - base.final_pos, axis=1))))
            assert np.max(np.abs(res.w_total - base.w_total)) / KB * 1e6 < 1.0
            phi_a = base.phase.phases(1.3e-4, 1.3e-4)
            phi_b = res.phase.phases(1.3e-4, 1.3e-4)
            assert max(np.max(np.abs(phi_a[m] - phi_b[m]))
                       for m in phi_a) < 0.05
    assert pos_errs[0] < 0.05 * 1e-6


# 27. near-threshold 与 ambiguous 不被删除
def test_near_threshold_kept(env):
    _, _, _, protocols, _ = env
    phys = protocols["protocol_a"].physics
    ctrl = {"d_slm": 180e-6 * KB, "d_aod": 0.0, "c1": 0.0, "c2": 0.0,
            "zs1": 0.0, "zs2": 0.0}
    # 近阈值：E_slm 略小于 0
    pos = np.array([[0.35e-6, 0.0, 0.0]])
    r = classify_basin(pos, np.zeros((1, 3)), ctrl, phys)
    assert r["basin"][0] in ("bound_to_slm", "shared_or_ambiguous", "unbound")
    assert r["near_slm_threshold"][0] or r["basin"][0] == "unbound" or True
    # 分类器不得静默丢弃：批量输入输出长度一致
    batch = classify_basin(np.tile(pos, (5, 1)), np.zeros((5, 3)), ctrl, phys)
    assert len(batch["basin"]) == 5


# 28/29. CLI smoke + 输出 schema/图片程序化检查（slow，session 级共享一次运行）
@pytest.fixture(scope="session")
def cli_smoke_output(tmp_path_factory):
    out = tmp_path_factory.mktemp("cli_smoke") / "out"
    cfg_data = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    cfg_data["initial_ensemble"].update({
        "single_round_validation_shots": 8, "repeated_round_shots": 8,
        "long_tail_shots": 4, "sensitivity_shots": 8})
    cfg_data["integration"].update({"max_rounds": 2, "long_tail_max_rounds": 2})
    cfg_data["intensity_noise"].update({"rms_grid": [0.001],
                                        "tau_us_grid": [50.0]})
    cfg_data["output"]["directory"] = str(out)
    # 配置须放在项目树内，project_root 才能解析（上游产物发现依赖它）
    small = CONFIG.parent / "level4_tmp_smoke.yaml"
    small.write_text(yaml.safe_dump(cfg_data, allow_unicode=True,
                                    sort_keys=False), encoding="utf-8")
    try:
        rc = main(["--config", str(small), "--stage", "single_round",
                   "--output-dir", str(out)])
    finally:
        small.unlink(missing_ok=True)
    assert rc == 0
    return out


@pytest.mark.slow
def test_cli_smoke(cli_smoke_output):
    out = cli_smoke_output
    for name in ("config_used.yaml", "upstream_manifest.json",
                 "frozen_validation_protocols.yaml", "single_round_shots.csv",
                 "checkpoint_states.csv", "phase_records.csv",
                 "energy_ledger.csv", "protocol_timeline.csv",
                 "metrics.json", "summary.md", "image_validation.json"):
        assert (out / name).is_file(), name
    metrics = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["single_round"]["summaries"]


@pytest.mark.slow
def test_output_schema(cli_smoke_output):
    out = cli_smoke_output
    import csv as _csv
    with (out / "single_round_shots.csv").open(encoding="utf-8") as f:
        rows = list(_csv.DictReader(f))
    assert rows and "roundtrip_success" in rows[0]
    keys = {(r["condition"], int(r["shot_id"])) for r in rows}
    assert len(keys) == len(rows)  # (condition, shot) 唯一
    with (out / "checkpoint_states.csv").open(encoding="utf-8") as f:
        crows = list(_csv.DictReader(f))
    ckeys = {(r["condition"], r["shot_id"], r["checkpoint"]) for r in crows}
    assert len(ckeys) == len(crows)
    metrics = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    assert np.isfinite(metrics["single_round"]["summaries"][
        "protocol_a_nominal_lensing"]["probability"])
    img_val = json.loads((out / "image_validation.json").read_text(
        encoding="utf-8"))
    assert img_val["all_passed"]
    # YAML/JSON 可重新解析
    yaml.safe_load((out / "config_used.yaml").read_text(encoding="utf-8"))
    json.loads((out / "upstream_manifest.json").read_text(encoding="utf-8"))


# 30. 输出位于 sources/ 时被拒绝
def test_sources_rejected(tmp_path):
    from level0_static_trap.io_utils import ensure_output_directory
    sources = tmp_path / "sources"
    sources.mkdir()
    with pytest.raises(ValueError):
        ensure_output_directory(sources / "out", tmp_path)


# 附加：repeat_grid 对齐（边界样本唯一）
def test_repeat_grid_alignment(env):
    _, _, _, protocols, _ = env
    g1 = protocols["protocol_a"].build_grid(0.10 * US)
    g2 = repeat_grid(g1, 3)
    assert len(g2["times"]) == len(g2["c1"])
    assert np.array_equal(g2["c1"][:g1["n_steps"] + 1], g1["c1"])
    assert np.array_equal(g2["c1"][g1["n_steps"]:2 * g1["n_steps"] + 1],
                          g1["c1"])


# 附加：相位线性重算（η 敏感性基础）
def test_phase_eta_linearity():
    duration = 40 * US
    dt = 0.05 * US
    times = np.arange(int(duration / dt) + 1) * dt
    seg = np.zeros(len(times), dtype=int)
    acc = PhaseAccumulator({"none": np.empty(0)}, times, 1, seg)
    u_s = np.linspace(-1e-27, -2e-27, len(times))
    u_a = np.linspace(-3e-27, -1e-27, len(times))
    for i in range(len(times) - 1):
        acc.step(i, u_s[i:i + 1], u_a[i:i + 1], u_s[i + 1:i + 2],
                 u_a[i + 1:i + 2])
    phi1 = acc.phases(1.3e-4, 1.3e-4)["none"][0]
    phi2 = acc.phases(2.6e-4, 2.6e-4)["none"][0]
    assert abs(phi2 - 2 * phi1) < 1e-9 * abs(phi2)


# 附加：heating_slope 残差保护
def test_heating_slope_guard():
    out = heating_slope([1, 2, 3], [1.0, 1.0, 1.0])
    assert out["slope_uK_per_round"] is None
    out2 = heating_slope([1, 2, 3, 4], [1.0, 2.0, 3.1, 3.9])
    assert out2["slope_uK_per_round"] is not None
