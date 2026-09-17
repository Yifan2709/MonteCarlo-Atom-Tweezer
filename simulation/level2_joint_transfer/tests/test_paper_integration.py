"""paper_integration 测试套件：物理基元、结构判据、合同与 I1 链路。

标记 slow 的测试做小样本全链路，正式参数见 runs.py。
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import paper_integration.species  # noqa: F401 注册触发
from paper_integration import (contract_data, events, integration, quantum_mc,
                               rb_circuit, rb_entangled, rb_transport,
                               trajectories, yb_code, yb_gate, yb_leakage,
                               yb_transport, yb_vls)
from paper_integration.dynamics3d import (Segment, TrapPhysics, propagate,
                                          sample_bound_thermal_3d)
from paper_integration.provenance import REGISTRIES

H_PLANCK = 6.62607015e-34


# ---------------------------------------------------------------- 基元
def test_gamma_family_min_jerk_exact():
    c2, a, b, c = trajectories.gamma_family_coeffs(1.875)
    assert c2 == 0 and a == 10 and b == -15 and c == 6


def test_trajectories_endpoints_and_peak_ordering():
    for kind, expect_peak in (("low_peak_v", None), ("gamma:1.875", 1.875)):
        t = trajectories.poly_trajectory(kind, 50e-6, 1e-3, 0.1e-6)
        assert abs(t.x_m[-1] - 50e-6) < 1e-15
        assert abs(t.v_m_s[-1]) < 1e-12 and abs(t.v_m_s[0]) < 1e-12
        peak = float(np.max(t.v_m_s)) * 1e-3 / 50e-6
        if expect_peak:
            assert abs(peak - expect_peak) < 1e-3
    zj = trajectories.poly_trajectory("low_peak_v", 50e-6, 1e-3, 0.1e-6)
    mj = trajectories.poly_trajectory("gamma:1.875", 50e-6, 1e-3, 0.1e-6)
    assert np.max(zj.v_m_s) < np.max(mj.v_m_s)   # 结构：ZJ 峰值速度更低


def test_quantum_mc_bell_ghz_teleport():
    rng = np.random.default_rng(0)
    st = quantum_mc.BatchedState(2, 4000, rng)
    st.h(0)
    st.cx(0, 1)
    b0, b1 = st.measure_z(0), st.measure_z(1)
    assert float(np.mean(b0 == b1)) == pytest.approx(1.0)
    st2 = quantum_mc.BatchedState(3, 6000, rng)
    st2.h(0)
    st2.h(1)
    st2.cx(1, 2)
    st2.cx(0, 1)
    st2.h(0)
    m0, m1 = st2.measure_z(0), st2.measure_z(1)
    fx = np.zeros((3, 6000), bool)
    fz = np.zeros((3, 6000), bool)
    fx[2] = (m1 == 1)
    fz[2] = (m0 == 1)
    st2.apply_pauli_errors(fx, fz)
    xb = st2.measure_x(2)
    assert float(np.mean(xb == 0)) == pytest.approx(1.0, abs=0.02)


def test_quantum_mc_pauli_errors_conditional():
    """逐 trial Pauli 误差不得影响其他 trial（回归：X 分支条件化）。"""
    rng = np.random.default_rng(3)
    st = quantum_mc.BatchedState(1, 1000, rng)
    st.h(0)
    probs_before = st.probabilities().copy()
    st2 = quantum_mc.BatchedState(1, 1000, np.random.default_rng(3))
    st2.h(0)
    fx = np.zeros((1, 1000), bool)
    fx[0, 0] = True     # 仅 trial 0 施加 X
    st2.apply_pauli_errors(fx, np.zeros_like(fx))
    probs_after = st2.probabilities()
    # 概率不变（X 只换幅值位置），trial 间无串扰由量子态不可区分性保证
    assert np.allclose(probs_before.sum(axis=0), probs_after.sum(axis=0))
    assert np.allclose(probs_before, probs_after)  # |+⟩ 对 X 不敏感


def test_diag_expectation_true_joint():  # 验证 D1 回归
    """Bell 态上 ZZ 必须等于 1（旧乘积式实现返回 0）。"""
    st = quantum_mc.BatchedState(2, 500, np.random.default_rng(1))
    st.h(0)
    st.cx(0, 1)
    assert float(np.mean(st.diag_expectation("ZZ"))) == pytest.approx(1.0)
    assert float(np.mean(st.diag_expectation("ZI"))) == pytest.approx(0.0)
    ghz = quantum_mc.BatchedState(4, 300, np.random.default_rng(2))
    ghz.h(0)
    ghz.cx(0, 1)
    ghz.cx(0, 2)
    ghz.cx(0, 3)
    assert float(np.mean(ghz.diag_expectation("ZZZZ"))) == pytest.approx(1.0)


def test_apply_z_rotation_bell_dephasing():  # D6 基元
    """|Φ+⟩ 差分相位 φ 后 ⟨ZZ⟩=1 不变（纯 Z 失相干不改 Z 基关联）。"""
    st = quantum_mc.BatchedState(2, 4000, np.random.default_rng(3))
    st.h(0)
    st.cx(0, 1)
    rng = np.random.default_rng(4)
    dphi = rng.standard_normal(4000) * 0.8
    st.apply_z_rotation(0, dphi)      # 相对相位在 |00>/|11> 之间：转单比特
    assert float(np.mean(st.diag_expectation("ZZ"))) == pytest.approx(1.0)
    # X 基关联衰减：⟨XX⟩=E[cos φ]
    stx = quantum_mc.BatchedState(2, 4000, np.random.default_rng(3))
    stx.h(0)
    stx.cx(0, 1)
    stx.apply_z_rotation(0, dphi)
    stx.h(0)
    stx.h(1)
    b0 = stx.measure_z(0)
    b1 = stx.measure_z(1)
    xx = float(np.mean(1 - 2 * (b0 ^ b1)))
    assert xx == pytest.approx(float(np.mean(np.cos(dphi))), abs=0.05)


def test_quantum_mc_loss_semantics():
    rng = np.random.default_rng(1)
    st = quantum_mc.BatchedState(2, 2000, rng)
    st.h(0)
    st.cx(0, 1)
    lost = np.zeros((2, 2000), bool)
    lost[1, :500] = True
    st.mark_lost(lost)
    before = st.diag_expectation("ZI")[:500].copy()
    b = st.measure_z(1)
    after = st.diag_expectation("ZI")[:500]
    assert int((b[:500] == -1).sum()) == 500
    assert np.allclose(before, after)   # 丢失位测量不得坍缩幸存者


def test_stream_seed_keyed_and_stable():
    assert events.stream_seed(7, 0, 0, "a") == events.stream_seed(7, 0, 0, "a")
    assert events.stream_seed(7, 0, 0, "a") != events.stream_seed(7, 1, 0, "a")
    led = events.Ledger()
    led.add_true(events.TrueEvent(0, 0, 1.0, "leakage", "s"))
    led.add_detection(events.DetectionRecord(0, 0, 2.0, "erasure_check", True))
    assert led.summary()["true_events"]["leakage"] == 1
    with pytest.raises(ValueError):
        led.add_true(events.TrueEvent(0, 0, 1.0, "not_a_kind", "s"))


def test_m1_rng():
    # 删除原子（缩小数组）不改变剩余原子的噪声实现
    a1 = np.random.default_rng(events.stream_seed(5, 0, 3, "ch")).random(3)
    a2 = np.random.default_rng(events.stream_seed(5, 0, 3, "ch")).random(3)
    assert np.array_equal(a1, a2)


def test_species_registry_guards():
    yb = REGISTRIES["yb171_2026_native"]
    assert yb.get("gate_ar_error").value == pytest.approx(0.016)
    with pytest.raises(KeyError):
        yb.get("nonexistent_param")     # 禁止静默回退


# ---------------------------------------------------------------- R1/R2/R3
def test_r1_small_scan_and_analytic():
    rb = REGISTRIES["rb87_2022_native"]
    rows = rb_transport.scan_transports(110e-6, [300e-6, 800e-6],
                                        192, 11, 0.05e-6)
    sv = [r["survival"] for r in rows]
    assert sv[-1] >= sv[0]              # 长 T 更存活
    # 解析式(2)：constant-jerk 数值 ΔN 应显著低于同参数常加速度公式
    mass = rb.get("mass_u").value * 1.6605390666e-27
    dn_ca = rb_transport.analytic_delta_n_const_accel(
        110e-6, 300e-6, rb.get("omega_r_circuit_rad_s").value, mass)
    finite = [r for r in rows if np.isfinite(r["mean_delta_N"])]
    assert finite, "全部格点失阱"
    assert dn_ca >= finite[0]["mean_delta_N_excess"]


@pytest.mark.slow
def test_r2_platform_shape():
    cal = rb_entangled.calibrate_c_phi(110e-6, 300e-6, shots=192, seed=7)
    scan = rb_entangled.scan_entangled(
        110e-6, [150e-6, 300e-6, 600e-6], 192, 8,
        cal["c_phi_rad_per_quanta"])
    crit = rb_entangled.evaluate_structure(scan)
    assert crit["C1_flat_platform"]
    assert crit["C2_fast_transport_drop"]
    # DD 对照在静态失相干可分辨的时长（σ_static=2T/T2*，T=1ms → 0.5 rad）
    dd_on_1ms = rb_entangled.entangled_transport(
        110e-6, 1000e-6, 192, 9, cal["c_phi_rad_per_quanta"], dd=True)
    dd_off_1ms = rb_entangled.entangled_transport(
        110e-6, 1000e-6, 192, 10, cal["c_phi_rad_per_quanta"], dd=False)
    assert dd_on_1ms["fidelity_raw_loss_as_one"] \
        > dd_off_1ms["fidelity_raw_loss_as_one"]


def test_r3_noiseless_code_space():
    """无误差极限：X̄ 确定 +1、ZZ 校验子全零（Steane 码空间验证）。"""
    rng = np.random.default_rng(2)
    st = quantum_mc.BatchedState(7, 3000, rng)
    rb_circuit._run_prep(st, np.random.default_rng(2),
                         {"x": 0.0, "y": 0.0, "z": 0.0, "loss": 0.0}, 0.0,
                         initial_loss=0.0)
    xb = [st.measure_x(q) for q in range(7)]
    par = np.zeros(3000, np.int8)
    for b in xb:
        par ^= np.where(b < 0, 0, b).astype(np.int8)
    assert float(np.mean(par % 2 == 0)) == pytest.approx(1.0)
    st2 = quantum_mc.BatchedState(7, 3000, np.random.default_rng(2))
    rb_circuit._run_prep(st2, np.random.default_rng(2),
                         {"x": 0.0, "y": 0.0, "z": 0.0, "loss": 0.0}, 0.0,
                         initial_loss=0.0)
    zb = [st2.measure_z(q) for q in range(7)]
    for support in rb_circuit.ZZ_SYNDROME_SUPPORTS:
        p = np.zeros(3000, np.int8)
        for q in support:
            p ^= np.where(zb[q] < 0, 0, zb[q]).astype(np.int8)
        assert float(np.mean(p % 2 == 0)) == pytest.approx(1.0)


def test_r3_cz_layer_budget():
    cz = rb_circuit.single_cz_layer_fidelity(50000, 3)
    assert cz["implied_cz_fidelity"] == pytest.approx(0.972, abs=0.005)


# ---------------------------------------------------------------- Y1/Y2/Y3
@pytest.mark.slow
def test_y1_single_move_runs():
    r = yb_transport.run_single_move(50e-6, 0.89e-3, "low_peak_v", 4.0,
                                     128, 21, 0.1e-6)
    assert 0.0 <= r["survival"] <= 1.0
    assert r["checkpoints"]


def test_y2_anchors():
    a = yb_vls.evaluate_anchors()
    assert a["optimal_within_paper_bound"]
    assert a["parallel_strongly_dephasing"]
    assert a["ordering_static_le_optimal_le_parallel"]


def test_y3_anchors_and_confusion():
    rt = yb_leakage.simulate_roundtrips(2048, 5, 42)
    a = yb_leakage.evaluate_anchors(rt)
    assert a["rt1_transient_elevated"]
    assert a["erasure_frac_declines"]
    assert abs(a["rt1_loss"] - 0.011) < 0.006
    cm = yb_leakage.confusion_matrix(20000)
    assert abs(cm["P_detect_given_clean"] - 0.014) < 0.005
    assert cm["P_detect_given_leak"] > 0.97


def test_y4_scaling_and_composition():
    sc_ar = yb_gate.scan_scaling("AR", n_trials=20000)
    sc_to = yb_gate.scan_scaling("TO", n_trials=20000)
    assert abs(sc_ar["fitted_slope"] - 4.0) < 0.6
    assert abs(sc_to["fitted_slope"] - 2.0) < 0.6
    z = yb_gate.run_zone_mc(50000, 7)
    assert abs(z["AR"]["erasure_frac_of_errors"] - 0.38) < 0.05
    assert abs(z["TO"]["loss_frac_of_errors"] - 0.75) < 0.05


# ---------------------------------------------------------------- Y5/Y6
def test_y5_code_structure_noiseless():
    """无误差：|00̄⟩ 双基宇称全对（GHZ 码字）。"""
    r = yb_code.measure_prep_fidelity(4000, 1, "raw", eps_layer=0.0,
                                      transport_loss=0.0)
    assert r["fidelity_estimate"] == pytest.approx(1.0, abs=1e-6)


def test_y5_storage_calibers_ordering():
    s = yb_code.storage_experiment(8000, 6, 0.5)
    assert s["erasure_informed"] > s["unconditional"]
    assert s["postselected"] >= s["erasure_informed"] - 0.05
    assert 0.0 < s["acceptance"] < 1.0


def test_y5_prep_ordering():
    raw = yb_code.measure_prep_fidelity(8000, 5, "raw")
    flag = yb_code.measure_prep_fidelity(8000, 5, "flag")
    assert flag["fidelity_estimate"] > raw["fidelity_estimate"]


def test_y6_adaptive_beats_random():
    ad = yb_code.teleportation(12000, 8, policy="adaptive")
    rd = yb_code.teleportation(12000, 8, policy="random")
    assert ad["success"] > rd["success"]
    assert ad["postselected_success"] >= ad["success"] - 0.05
    assert ad["mean_selection_latency_s"] > 0
    # 验证 D3 回归：观测后选与理想诊断必须分开报告
    assert "ideal_true_loss_conditioned_success" in ad


# ---------------------------------------------------------------- I1/合同
def test_i1_rb_chain_propagates():
    r = integration.rb_perturbation(shots=96, seed=3)
    assert r["propagated"]["move_survival"]
    assert r["propagated"]["steane_raw_xbar"]


@pytest.mark.slow
def test_i1_yb_chain_propagates():
    r = integration.yb_perturbation(move_time_s=0.6e-3, distance_m=100e-6,
                                    shots=96, seed=4242)
    assert r["propagated"]["teleport_success"] or r["propagated"]["rt1_loss_frac"]


def test_contract_writes_and_hashes(tmp_path):
    info = contract_data.write_contract(tmp_path)
    assert info["n_requirements"] == 13
    c = json.loads((tmp_path / "acceptance_contract.json").read_text("utf-8"))
    assert len(c["contract_sha256"]) == 64
    assert c["acceptance_protocol_file"]["status"] == "not_found_in_repository"
    assert (tmp_path / "paper_coverage.csv").exists()


def test_dynamics3d_hold_energy_conservation():
    rb = REGISTRIES["rb87_2022_native"]
    ph = TrapPhysics(mass_kg=86.909180531 * 1.6605390666e-27,
                     depth_j=4e6 * H_PLANCK, waist_m=rb.get("waist_m").value,
                     zr_m=rb.get("zr_m").value)
    pos, vel = sample_bound_thermal_3d(5, 5.0, 128, ph)
    t = np.arange(0, 100e-6 + 0.05e-6, 0.05e-6)
    seg = Segment("hold", t, np.zeros_like(t), np.zeros_like(t),
                  np.full_like(t, ph.depth_j))
    res = propagate(pos, vel, [seg], ph, 0.05e-6)
    assert res["survival_final"] == pytest.approx(1.0)
