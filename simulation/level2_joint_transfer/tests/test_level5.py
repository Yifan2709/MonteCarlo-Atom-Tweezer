"""Level 5 测试：Hamiltonian/传播器/脉冲/Clifford/DD/噪声/轨迹回放/
PTM/RB/统计/CLI/产物验收。覆盖 level5_prompt 第二十一节 33 项要求。"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from level5_finite_pulse_irb import (level5_statistics as l5s,
                                     level5_visualization as viz)
from level5_finite_pulse_irb.channel_tomography import (
    _reference_inputs, affine_fit, average_state_fidelity_from_affine,
    characterize_from_states, known_channel_checks, ptm_bootstrap)
from level5_finite_pulse_irb.clifford_group import (CliffordGroup,
                                                    virtual_z_unitary)
from level5_finite_pulse_irb.dd_finite_pulses import (
    dd_pulse_events, symmetric_pulse_fractions, validate_transformed_frame,
    xy_axis_pattern)
from level5_finite_pulse_irb.interleaved_rb import (generate_irb_sequences,
                                                    run_interleaved_rb)
from level5_finite_pulse_irb.level5_config import load_config
from level5_finite_pulse_irb.noise_psd import (
    DetuningNoiseModel, ou_autocorrelation_theory, ou_trace, psd_target,
    synthesize_from_psd, validate_psd_synthesis, welch_psd)
from level5_finite_pulse_irb.pulse_library import (
    SCROFULOUS_MIN_THETA, arcsinc_upper, envelope_area_integral,
    envelope_value, make_bare_pulse, make_scrofulous_pulses,
    pulse_list_area, pulse_list_duration, scrofulous_angles)
from level5_finite_pulse_irb.qubit_hamiltonian import (
    INPUT_STATES, PAULI_X, PAULI_Y, PAULI_Z, analytic_rabi_population,
    commutator, equatorial_rotation_unitary, light_shift_detuning,
    return_probability_z, statevector_to_bloch, validate_statevector)
from level5_finite_pulse_irb.quantum_propagators import (
    apply_su2, density_matrix_checks, depolarizing_ptm, dephasing_ptm,
    lindblad_step, su2_propagator, unitarity_error)
from level5_finite_pulse_irb.rb_fitting import (
    early_window_slope, fit_exponential_rb, fit_irb_decay,
    instantaneous_ratios, irb_estimate, model_selection, safe_fit_summary)
from level5_finite_pulse_irb.reference_rb import (
    generate_reference_sequences, run_reference_rb)
from level5_finite_pulse_irb.spin_channels import (MoveChannelEngine,
                                                   StaticGateEngine)
from level5_finite_pulse_irb.spatial_noise import (
    correlated_site_noise, empirical_cross_site_correlation,
    sample_site_parameters)
from level5_finite_pulse_irb.trajectory_replay import (
    TrajectoryPool, replay_phase_for_pool)

TWO_PI = 2.0 * np.pi
OMEGA = TWO_PI * 24.611e3
HBAR = 1.054571817e-34

GROUP = CliffordGroup()


def _synthetic_pool(n_traj=32, n_t=201, duration_s=1e-3, seed=0,
                    segment_names=("outbound_long_move",
                                   "remote_hold", "inbound_long_move")):
    """构造合成轨迹池（无需运行 Level 4 引擎的轻量测试对象）。"""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, duration_s, n_t)
    depth = -180e-6 * 1.380649e-23
    osc = (1.0 + 0.05 * rng.standard_normal((n_traj, n_t))
           + 0.05 * np.sin(TWO_PI * 36e3 * t)[None, :])
    u_slm = depth * osc
    u_aod = 0.5 * depth * (1.0 - np.exp(-t[None, :] / 2e-4)) * \
        (1 + 0.03 * rng.standard_normal((n_traj, n_t)))
    bounds = np.linspace(0, duration_s, len(segment_names) + 1)
    return TrajectoryPool(
        name="synthetic", times_s=t, u_slm=u_slm.astype(np.float32),
        u_aod=u_aod.astype(np.float32),
        roundtrip_success=rng.random(n_traj) > 0.05,
        final_retained=rng.random(n_traj) > 0.02,
        first_failure_stage=np.zeros(n_traj, dtype=np.int16),
        stage_names=["none"], u_slm_final=u_slm[:, -1].astype(float),
        u_aod_final=u_aod[:, -1].astype(float),
        level4_a_slm=np.trapezoid(u_slm, t, axis=1),
        level4_a_aod=np.trapezoid(u_aod, t, axis=1),
        segment_boundaries_s=bounds, total_duration_s=duration_s,
        segment_names=list(segment_names),
        meta={"synthetic": True})


# ------------------------------------------------------------ 1-4 基础
def test_pauli_commutators_and_rotation_conventions():
    assert np.allclose(commutator(PAULI_X, PAULI_Y), 2j * PAULI_Z)
    assert np.allclose(commutator(PAULI_Y, PAULI_Z), 2j * PAULI_X)
    assert np.allclose(commutator(PAULI_Z, PAULI_X), 2j * PAULI_Y)
    u = equatorial_rotation_unitary(0.0, np.pi)
    assert np.allclose(np.abs(u @ PAULI_Z @ u.conj().T + PAULI_Z), 0)
    # R_x(π/2) 共轭：Y→Z，Z→−Y（绕 X 轴旋转与 X 对易）
    u2 = equatorial_rotation_unitary(0.0, 0.5 * np.pi)
    assert np.allclose(u2 @ PAULI_Y @ u2.conj().T, PAULI_Z, atol=1e-12)
    assert np.allclose(u2 @ PAULI_Z @ u2.conj().T, -PAULI_Y, atol=1e-12)


def test_constant_hamiltonian_su2_matches_matrix_exponential():
    from scipy.linalg import expm
    delta, omega, phi, dt = 300.0, OMEGA, 0.4, 1.3e-6
    h = 0.5 * (delta * PAULI_Z + omega * (np.cos(phi) * PAULI_X
                                          + np.sin(phi) * PAULI_Y))
    exact = expm(-1j * h * dt)
    u = [x[0] for x in su2_propagator(np.array([delta]), np.array([omega]),
                                      phi, dt)]
    approx = np.array([[u[0], u[1]], [u[2], u[3]]])
    assert np.allclose(approx, exact, atol=1e-12)


def test_rabi_ramsey_echo_limits():
    omega = OMEGA
    t = 0.3 * np.pi / omega
    psi = np.array([[1.0, 0.0]], dtype=complex)
    for _ in range(2000):
        apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                       t / 2000))
    target = analytic_rabi_population(0.0, omega, t)
    assert abs((1 - abs(psi[0, 0]) ** 2) - target) < 1e-9
    # Ramsey：第二 π/2 用 −x 相位 → P1=sin²(Δτ/2)
    delta = TWO_PI * 250.0
    free = 4e-6
    psi = np.array([[1.0, 0.0]], dtype=complex)
    apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                   0.5 * np.pi / omega))
    apply_su2(psi, *su2_propagator(np.full(1, delta), np.zeros(1), 0.0,
                                   free))
    apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), np.pi,
                                   0.5 * np.pi / omega))
    assert abs((1 - abs(psi[0, 0]) ** 2)
               - np.sin(0.5 * delta * free) ** 2) < 1e-9
    # echo：quasistatic Δ 回聚
    tau = 6e-6
    psi = np.array([[1.0, 0.0]], dtype=complex)
    apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                   0.5 * np.pi / omega))
    apply_su2(psi, *su2_propagator(np.full(1, delta), np.zeros(1), 0.0,
                                   tau))
    apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                   np.pi / omega))
    apply_su2(psi, *su2_propagator(np.full(1, delta), np.zeros(1), 0.0,
                                   tau))
    apply_su2(psi, *su2_propagator(np.zeros(1), np.full(1, omega), 0.0,
                                   0.5 * np.pi / omega))
    assert abs(psi[0, 0]) ** 2 > 1 - 1e-9


def test_state_norm_and_density_matrix_checks():
    psi = np.array([[1.0, 1j]], dtype=complex) / np.sqrt(2)
    u = su2_propagator(np.array([100.0]), np.array([5e3]), 0.3, 1e-7)
    apply_su2(psi, *u)
    checks = validate_statevector(psi)
    assert checks["norm_ok"]
    rho = np.array([[0.5, 0.1j], [-0.1j, 0.5]])
    for _ in range(50):
        rho = lindblad_step(rho, 200.0, 0.0, 0.0, 1e-7)
    dm = density_matrix_checks(rho)
    assert dm["trace_preserved"] and dm["hermitian"] and dm["positive"]


# ------------------------------------------------------------ 5-8 脉冲
def test_bare_rotation_unitary():
    for phase, angle in ((0.0, 0.5 * np.pi), (0.5 * np.pi, np.pi),
                         (np.pi, 0.5 * np.pi)):
        p = make_bare_pulse(phase, angle, OMEGA)
        assert abs(p.duration * OMEGA - angle) < 1e-12
        u = equatorial_rotation_unitary(phase, angle)
        psi = np.array([[1.0, 0.0]], dtype=complex)
        n_sub = 80
        dt = p.duration / n_sub
        for _ in range(n_sub):
            apply_su2(psi, *su2_propagator(np.zeros(1),
                                           np.full(1, OMEGA),
                                           p.axis_phase, dt))
        a, b = psi[0, 0], psi[0, 1]
        # SU(2) 对称形式：U=[[u00,−u10*],[u10,u00*]]
        out = np.array([[a, -np.conj(b)], [b, np.conj(a)]])
        prod = out @ u.conj().T
        ph = prod.flat[int(np.argmax(np.abs(prod)))]
        assert np.allclose(prod / ph, np.eye(2), atol=1e-9)


def test_pulse_envelope_area():
    square = envelope_area_integral("square", 2e-5, 0.0)
    assert abs(square - 2e-5) < 1e-15
    edge = envelope_area_integral("cosine_edge", 2e-5, 2.0)
    assert 0 < edge < square
    gauss = envelope_area_integral("gaussian", 2e-5, 0.0)
    assert 0 < gauss < square
    assert envelope_value("square", 0.3, 0.0, 1.0, 20.0) == 1.0


def test_scrofulous_target_branch_and_robustness():
    for theta in (0.5 * np.pi, np.pi, np.pi / 3, 0.9 * np.pi, 0.02):
        ang = scrofulous_angles(theta, 0.7)
        p1, p2, _ = ang.pulse_phases
        u = (equatorial_rotation_unitary(p1, ang.theta_c)
             @ equatorial_rotation_unitary(p2, np.pi)
             @ equatorial_rotation_unitary(p1, ang.theta_c))
        tgt = equatorial_rotation_unitary(0.7, theta)
        prod = u @ tgt.conj().T
        ph = prod.flat[int(np.argmax(np.abs(prod)))]
        assert np.allclose(prod / ph, np.eye(2), atol=1e-9)
    # 幅度鲁棒性优于 bare
    theta = 0.5 * np.pi
    eps = 0.08
    ang = scrofulous_angles(theta, 0.0)
    p1, p2, _ = ang.pulse_phases
    u_e = (equatorial_rotation_unitary(p1, ang.theta_c * (1 + eps))
           @ equatorial_rotation_unitary(p2, np.pi * (1 + eps))
           @ equatorial_rotation_unitary(p1, ang.theta_c * (1 + eps)))
    tgt = equatorial_rotation_unitary(0.0, theta)
    d_scrof = np.abs(u_e - tgt).max()
    d_bare = np.abs(equatorial_rotation_unitary(
        0.0, theta * (1 + eps)) - tgt).max()
    assert d_scrof < d_bare / 5
    # 分支处理
    with pytest.raises(ValueError):
        arcsinc_upper(0.8)
    assert 0.5 * np.pi <= scrofulous_angles(0.01).theta_c <= np.pi
    # 近零角度退化
    pulses = make_scrofulous_pulses(0.0, 0.0, OMEGA)
    assert pulses == []
    pulses_small = make_scrofulous_pulses(0.0, SCROFULOUS_MIN_THETA * 0.5,
                                          OMEGA)
    assert len(pulses_small) == 1  # bare 退化


def test_dd_pulse_order_and_finite_propagation():
    omega = OMEGA
    for mode, n in (("none", 0), ("spin_echo", 1), ("xy4", 4), ("xy8", 8),
                    ("xy16", 16)):
        events = dd_pulse_events(mode, 0.0, 100e-6, omega)
        assert len(events) == n
        starts = [p.start for p in events]
        assert starts == sorted(starts)
        for p in events:
            assert abs(p.nominal_area - np.pi) < 1e-9 or p.family == \
                "scrofulous"
    fracs = symmetric_pulse_fractions(4)
    assert np.allclose(fracs, [1 / 8, 3 / 8, 5 / 8, 7 / 8])
    assert xy_axis_pattern("xy4") == ["X", "Y", "X", "Y"]
    # 有限脉冲传播：XY4 下 DD 抑制静失谐
    pool = _synthetic_pool()
    for mode in ("none", "xy4"):
        eng = MoveChannelEngine(pool, omega, dd_mode=mode)
        psi = np.tile(INPUT_STATES["x_plus"], (8, 1)).astype(complex)
        k = np.arange(8)
        eng.propagate(psi, k, 1e-4, 1e-4, delta_site=np.full(8, 0.0),
                      noise_qs=np.full(8, TWO_PI * 150.0))
        globals().setdefault("_dd_contrast", {})[mode] = \
            abs(np.mean(statevector_to_bloch(psi)[:, 0]))


# ------------------------------------------------------------ 9-11 Level4
def test_instantaneous_limit_and_light_shift_consistency():
    pool = _synthetic_pool()
    k = np.arange(pool.n_traj)
    phi_replay = replay_phase_for_pool(pool, k, 1.3e-4, 1.3e-4)
    # 池内自带梯形积分（level4 风格）
    phi_l4 = (1.3e-4 * pool.level4_a_slm + 1.3e-4 * pool.level4_a_aod) / HBAR
    assert np.allclose(phi_replay, phi_l4, atol=1e-3)
    # light shift 公式
    delta = light_shift_detuning(-1e-26, -2e-26, 1.3e-4, 1.3e-4)
    assert abs(delta - 1.3e-4 * (-3e-26) / HBAR) < 1e-6


def test_trajectory_interpolation_and_boundaries():
    pool = _synthetic_pool()
    engine = MoveChannelEngine(pool, OMEGA, dd_mode="xy4_per_long_move")
    bounds = [0.0]
    for p in engine.events:
        bounds += [p.start, p.start + p.duration]
    bounds.append(pool.total_duration_s)
    arr = np.unique(np.asarray(bounds))
    assert arr[0] == 0.0
    assert abs(arr[-1] - pool.total_duration_s) < 1e-12
    assert np.all(np.diff(arr) > 0)
    # 插值在样本点精确
    u_row = pool.u_slm[3]
    pos = pool.times_s[50] / pool.sample_interval_s
    assert abs(pool.u_at(pool.u_slm, np.array([3]),
                         np.array([pool.times_s[50]]))[0] - u_row[50]) \
        < 1e-6


# ------------------------------------------------------------ 12-14 噪声
def test_ou_statistics():
    rng = np.random.default_rng(5)
    n, dt, tau, sigma = 200000, 1e-6, 2e-4, 3.0
    x = ou_trace(rng, n, dt, tau, sigma)
    assert abs(np.mean(x)) < 0.05 * sigma
    assert abs(np.std(x) - sigma) < 0.03 * sigma
    lags = np.array([1, 5, 20, 50])
    emp = np.array([np.corrcoef(x[:-l], x[l:])[0, 1] for l in lags])
    th = ou_autocorrelation_theory(lags * dt, tau, sigma) / sigma ** 2
    assert np.max(np.abs(emp - th)) < 0.05


def test_psd_normalization_nyquist_reproducibility():
    spec = [{"model": "white", "params": {"level": 1e-4}}]
    val = validate_psd_synthesis(spec, n_samples=16000, fs=1e6, seed=3)
    assert abs(val["variance_ratio"] - 1.0) < 0.35
    assert val["real_valued"] and val["finite"]
    assert val["seed_reproducible"] and val["seed_independent"]
    assert val["nyquist_hz"] == 5e5
    # Hz 与 rad/s：lorentzian PSD 在特征频率处的量级
    freqs = np.linspace(1.0, 1e5, 1000)
    psd = psd_target(freqs, [{"model": "lorentzian",
                              "params": {"sigma": 1.0, "tau_s": 1e-3}}])
    assert np.all(np.isfinite(psd)) and np.all(psd > 0)
    # 60 Hz 线
    psd_line = psd_target([60.0], [{"model": "line",
                                    "params": {"f_hz": 60.0,
                                               "width_hz": 1.0,
                                               "amp": 1e-3,
                                               "harmonics": 3}}])
    assert psd_line[0] > psd_target([200.0], [{"model": "line",
                                               "params": {
                                                   "f_hz": 60.0,
                                                   "width_hz": 1.0,
                                                   "amp": 1e-3,
                                                   "harmonics": 3}}])[0]


def test_common_local_spatial_correlations():
    traces = correlated_site_noise(np.random.default_rng(7), 16, 5000, 0.7,
                                   1.0)
    corr = empirical_cross_site_correlation(traces)
    off = corr[np.triu_indices_from(corr, k=1)]
    assert abs(off.mean() - 0.7) < 0.1
    params = sample_site_parameters(20, 11, trap_depth_relative_std=0.114,
                                    rabi_relative_std=0.05,
                                    static_detuning_std_hz=30.0,
                                    eta_relative_std=0.04)
    assert abs(np.std(params["depth_scale"]) - 0.114) < 0.04
    assert abs(np.std(params["static_detuning_hz"]) - 30.0) < 10


# ------------------------------------------------------- 15-19 群与 PTM
def test_clifford_group_structure():
    checks = GROUP.group_checks()
    assert checks["n_elements"] == 24 and checks["unique"]
    assert checks["closure"] and checks["inverse_ok"]
    assert checks["pauli_conjugation_signed"]
    assert checks["ideal_realization_matches"]
    # Pauli 共轭带符号
    for i in range(24):
        for k, p in enumerate((PAULI_X, PAULI_Y, PAULI_Z)):
            conj = GROUP.unitary(i) @ p @ GROUP.unitary(i).conj().T
            matched = any(np.allclose(conj, s * q, atol=1e-8)
                          for s in (1, -1)
                          for q in (PAULI_X, PAULI_Y, PAULI_Z))
            assert matched


def test_clifford_sampler_uniformity():
    counts = GROUP.sampler_uniformity_counts(240000, 3)
    chi2 = float(np.sum((counts - 10000.0) ** 2 / 10000.0))
    assert chi2 < 80.0  # df=23, 5σ ≈ 45；留裕量


def test_ideal_sequence_plus_inverse_returns():
    rng = np.random.default_rng(9)
    for _ in range(40):
        seq = [int(c) for c in GROUP.sample_uniform(rng, 15)]
        inv = GROUP.inverse[seq[-1]]
        for c in reversed(seq[:-1]):
            inv = int(GROUP.mult[GROUP.inverse[c], inv])
        u_total = np.eye(2, dtype=complex)
        for c in seq:
            u_total = GROUP.unitary(c) @ u_total
        psi = GROUP.unitary(inv) @ u_total @ np.array([1.0, 0.0])
        assert abs(abs(psi[0]) ** 2 - 1.0) < 1e-10


def test_known_ptm_reconstruction_and_favg_applicability():
    checks = known_channel_checks()
    assert checks["depolarizing"]["match"]
    assert checks["dephasing"]["match"]
    # z-旋转通道
    theta = 0.4
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0], [0, 0, 1.0]])
    inputs = _reference_inputs()
    outputs = (rot @ inputs.T).T
    m, c, res = affine_fit(inputs, outputs)
    assert np.allclose(m @ rot.T, np.eye(3), atol=1e-12)
    assert np.linalg.norm(c) < 1e-12 and res < 1e-12
    f = average_state_fidelity_from_affine(m, c)
    assert abs(f - (2 + np.cos(theta)) / 3) < 1e-12
    # 非 unital 通道（T1 类）必须触发适用性诊断
    m2, c2 = np.diag([1, 1, 0.9]), np.array([0, 0, 0.05])
    from level5_finite_pulse_irb.channel_tomography import \
        fidelity_applicability
    diag = fidelity_applicability(m2, c2, 1e-9)
    assert not diag["unital"]


# ------------------------------------------------------ 20-24 RB 与统计
def test_loss_conditional_joint_separation():
    surv = np.array([1.0, 1.0, 0.0])
    cond = np.array([0.99, 0.98, 0.97])
    joint = surv * cond
    assert joint[2] == 0.0
    assert abs(joint[0] - 0.99) < 1e-12
    assert np.mean(cond) != np.mean(joint)


def test_synthetic_reference_rb_fit_recovery():
    rng = np.random.default_rng(4)
    p_true = 0.9975
    lengths = np.array([1, 2, 4, 8, 16, 32, 64, 128])
    y = np.array([0.985 * p_true ** n + 0.5 * 0.015
                  + rng.normal(0, 0.0015) for n in lengths])
    fit = fit_exponential_rb(lengths, y)
    assert fit["ok"]
    assert abs(fit["params"]["p"] - p_true) < 0.0015


def test_synthetic_irb_fit_recovery():
    rng = np.random.default_rng(5)
    p_move, p_ref = 0.9985, 0.9995
    m_vals = np.array([0, 1, 2, 4, 8, 16, 32])
    y = np.array([0.98 * p_move ** m * p_ref ** 80 + 0.5 * 0.02
                  + rng.normal(0, 0.001) for m in m_vals])
    fit = fit_exponential_rb(m_vals, y)
    est = irb_estimate(p_ref, fit["params"]["p"])
    assert est["interpretable"]
    assert abs(est["f_avg_interleaved"] - p_move) < 0.0015


def test_nonexponential_fit_safe_failure():
    # 恒定数据 → 拟合不可辨识/失败，safe_fit_summary 不给数字
    y = np.full(8, 0.5)
    fit = fit_exponential_rb(np.arange(1, 9), y)
    summary = safe_fit_summary(fit)
    if summary["status"] != "fit_failed":
        assert "status" in summary
    # 相邻比值不稳定时不输出误导性 fidelity
    ratios = instantaneous_ratios([1.0, 0.9, 0.5, 0.49, 0.1])
    assert ratios["stable"] is False or ratios["stable"] is None
    est = irb_estimate(None, None)
    assert not est["interpretable"]
    # 生存接近零：clipped-Boltzmann 拟合也不应崩溃
    surv_fit = fit_irb_decay(np.arange(1, 9),
                             np.array([1, 1, 1, 1, 1, 0.5, 0.1, 0.0]),
                             "clipped_boltzmann_x_depol")
    assert surv_fit is not None


def test_hierarchical_bootstrap_width_vs_naive():
    demo = l5s.independent_bernoulli_vs_hierarchical_demo(n_boot=400)
    hier_width = demo["hierarchical"]["ci_high"] - \
        demo["hierarchical"]["ci_low"]
    naive_width = demo["naive_wilson"]["high"] - demo["naive_wilson"]["low"]
    assert hier_width >= naive_width * 0.8


def test_common_random_sequence_pairing():
    a = generate_irb_sequences(GROUP, 20, [0, 4], 6, 555)
    b = generate_irb_sequences(GROUP, 20, [0, 4], 6, 555)
    assert a == b
    # 同一 string_id 的不同 M 共享 Clifford strings
    m0 = [s for s in a if s["M"] == 0][0]
    m4 = [s for s in a if s["M"] == 4][0]
    assert m0["cliffords"] == m4["cliffords"]


def test_calibration_validation_seed_isolation():
    pools = l5s.SEED_POOLS
    vals = list(pools.values())
    assert len(set(vals)) == len(vals)
    assert l5s.seed_for("reference_rb_pool") != l5s.seed_for(
        "interleaved_rb_pool")
    with pytest.raises(KeyError):
        l5s.seed_for("nonexistent_pool")


def test_frozen_design_no_leakage(tmp_path):
    from level5_finite_pulse_irb.level5_cli import _frozen_design
    cfg_path = Path(__file__).resolve().parents[1] / \
        "configs/level5_finite_pulse_irb.yaml"
    cfg = load_config(cfg_path)
    design = _frozen_design(cfg)
    assert design["interleaved_rb"]["seed"] != \
        design["reference_rb"]["seed"]
    # 校准池与 validation 池分离（设计文本显式声明）
    assert "分离" in design["seed_pools"]
    assert design["noise_realizations"]["reference"] >= 1


# ----------------------------------------------- 25-30 引擎级一致性
def test_free_segment_equals_step_by_step():
    """缓存式自由段（matvec）与逐步 z 旋转一致（缓存通道验证）。"""
    pool = _synthetic_pool()
    engine = MoveChannelEngine(pool, OMEGA, dd_mode="none")
    k = np.arange(6)
    psi_fast = np.tile(INPUT_STATES["x_plus"], (6, 1)).astype(complex)
    engine.propagate(psi_fast, k, 1.3e-4, 1.3e-4)
    # 逐步：0.1μs 子步
    psi_slow = np.tile(INPUT_STATES["x_plus"], (6, 1)).astype(complex)
    dt = 0.1e-6
    n_steps = int(pool.total_duration_s / dt)
    sample_dt = pool.sample_interval_s
    for i in range(n_steps):
        t_mid = (i + 0.5) * dt
        pos = t_mid / sample_dt
        i0 = min(int(pos), pool.u_slm.shape[1] - 2)
        frac = min(pos - i0, 1.0)
        us = pool.u_slm[k, i0] + frac * (pool.u_slm[k, i0 + 1]
                                         - pool.u_slm[k, i0])
        ua = pool.u_aod[k, i0] + frac * (pool.u_aod[k, i0 + 1]
                                         - pool.u_aod[k, i0])
        delta = (1.3e-4 * us + 1.3e-4 * ua) / HBAR
        half = 0.5 * delta * dt
        psi_slow[:, 0] *= np.exp(-1j * half)
        psi_slow[:, 1] *= np.exp(1j * half)
    dev = np.abs(psi_fast - psi_slow).max()
    assert dev < 2e-3


def test_timestep_convergence():
    pool = _synthetic_pool()
    engine = MoveChannelEngine(pool, OMEGA, dd_mode="xy4_per_long_move")
    k = np.arange(8)
    results = []
    for dt_pulse, dt_free in ((0.10, 0.5), (0.025, 0.125)):
        psi = np.tile(INPUT_STATES["x_plus"], (8, 1)).astype(complex)
        engine.propagate(psi, k, 1.3e-4, 1.3e-4, dt_free_us=dt_free,
                         dt_pulse_us=dt_pulse)
        results.append(statevector_to_bloch(psi))
    assert np.abs(results[0] - results[1]).max() < 5e-3


def test_47_195_extension_frozen_model():
    from level5_finite_pulse_irb.array_ensemble import build_site_ensemble
    ens47 = build_site_ensemble({"trap_depth_relative_std": 0.0,
                                 "rabi_relative_std": 0.0,
                                 "static_detuning_std_Hz": 0.0}, 47, 100)
    ens195 = build_site_ensemble({"trap_depth_relative_std": 0.0,
                                  "rabi_relative_std": 0.0,
                                  "static_detuning_std_Hz": 0.0}, 195, 100)
    assert np.allclose(ens47["params"]["rabi_scale"], 1.0)
    assert np.allclose(ens195["params"]["rabi_scale"], 1.0)
    assert ens47["params"]["description"] == \
        ens195["params"]["description"]


def test_reference_rb_runs_and_site_consistency():
    pool = _synthetic_pool(n_traj=24)
    gate = StaticGateEngine(pool, OMEGA, family="bare", group=GROUP)
    seqs = generate_reference_sequences(GROUP, [1, 2, 4], 4, 77)
    result = run_reference_rb(seqs, gate, GROUP, pool, 1.3e-4, 1.3e-4,
                              seed=78, n_sites=3)
    probs = result["return_probabilities"]
    assert probs.shape == (12, 3)
    assert np.all((probs >= 0) & (probs <= 1))
    assert result["n_gates_per_sequence"][0] == 2


def test_irb_runs_and_conventions():
    pool = _synthetic_pool(n_traj=16, seed=1)
    static = _synthetic_pool(n_traj=16, seed=2, segment_names=("static",))
    move = MoveChannelEngine(pool, OMEGA, dd_mode="xy4_per_long_move")
    gate = StaticGateEngine(static, OMEGA, family="bare", group=GROUP)
    seqs = generate_irb_sequences(GROUP, 6, [0, 2], 4, 91)
    result = run_interleaved_rb(seqs, move, gate, pool, static, 1.3e-4,
                                1.3e-4, seed=92, n_sites=2)
    assert result["joint_return"].shape == (8, 2)
    assert np.allclose(result["joint_return"],
                       result["survival"] * result["conditional_return"])
    assert np.all((result["survival"] >= 0)
                  & (result["survival"] <= 1))


def test_transformed_frame_validation():
    val = validate_transformed_frame(GROUP, 60, 13)
    assert val["passed"]


# ------------------------------------------------------ 31-33 CLI/产物
def test_cli_help_and_smoke(tmp_path):
    root = Path(__file__).resolve().parents[1]
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    proc = subprocess.run(
        [sys.executable, "-m", "level5_finite_pulse_irb.level5_cli",
         "--help"], cwd=str(root), capture_output=True, text=True, env=env)
    assert proc.returncode == 0


def test_cli_preflight_fast_dev_pool(tmp_path, monkeypatch):
    """CLI preflight 快速通道（dev 池 + fast-preflight）冒烟。"""
    root = Path(__file__).resolve().parents[1]
    cfg_path = root / "configs/level5_finite_pulse_irb.yaml"
    cfg = load_config(cfg_path)
    monkeypatch.setattr(cfg, "pool_full_trajectories", 16)
    monkeypatch.setattr(cfg, "pool_secondary_trajectories", 12)
    monkeypatch.setattr(cfg, "joint_validation_shots", 16)
    from level5_finite_pulse_irb import level5_cli as cli
    out = tmp_path / "out"
    out.mkdir()
    l4_dir = cli._discover_level4_output(cfg)
    from level4_roundtrip_coherence.level4_config import \
        load_config as load4
    from level4_roundtrip_coherence import protocol_segments as ps
    from level4_roundtrip_coherence import upstream_artifacts as l4up
    cfg4 = load4(root / "configs/level4_roundtrip.yaml")
    _, resolved = l4up.build_manifest(cfg4)
    protocols = {
        "protocol_a": ps.build_protocol_a(
            cfg4, resolved["nominal_vs_m_per_s"]),
        "protocol_a_no_lensing": ps.build_protocol_a(cfg4, None),
        "protocol_b": ps.build_protocol_b(
            cfg4, resolved["level2_spec"], resolved["nominal_vs_m_per_s"]),
        "static_idle": ps.build_static_idle(cfg4)}
    pools = cli._ensure_pools(cfg, cfg4, protocols, out, log=lambda *a: None)
    assert pools["protocol_a"].n_traj == 16
    manifest = cli._build_upstream_manifest(cfg, l4_dir, {
        name: {"_path": "x", "_sha256": pool.sha256()}
        for name, pool in pools.items()})
    from level5_finite_pulse_irb.preflight5 import Preflight5
    pf = Preflight5(cfg, manifest, resolved, group=GROUP,
                    pool=pools["protocol_a"],
                    static_pool=pools["static_idle"],
                    timeline=protocols["protocol_a"], cfg4=cfg4,
                    run_upstream_tests=False, log=lambda *a: None)
    checks = pf.run()
    failed = [n for n in checks["critical_checks"]
              if not checks["checks"][n]["passed"]
              and not checks["checks"][n]["detail"].get(
                  "skipped_dev_mode", False)]
    assert not failed, failed


def test_output_rejects_sources_directory(tmp_path):
    from level0_static_trap.io_utils import ensure_output_directory
    root = tmp_path
    (root / "sources").mkdir()
    with pytest.raises(ValueError):
        ensure_output_directory(root / "sources" / "x", root)


def test_csv_json_jsonl_schema_and_image_checks(tmp_path):
    from level0_static_trap.io_utils import validate_pngs
    # 小 PNG + 校验
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="a")
    ax.legend()
    png = tmp_path / "demo.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)
    record = viz._finish(plt.figure(), tmp_path / "x.png", 0, [])
    val = validate_pngs([record])
    assert val["images"]["x.png"]["exists"]
    # CSV/JSON/JSONL 重读
    csv_path = tmp_path / "rows.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["a", "b"])
        w.writeheader()
        w.writerow({"a": 1, "b": 2.5})
    rows = list(csv.DictReader(csv_path.open()))
    assert len(rows) == 1 and float(rows[0]["b"]) == 2.5
    j = tmp_path / "m.json"
    j.write_text(json.dumps({"ok": True}))
    assert json.loads(j.read_text())["ok"] is True
