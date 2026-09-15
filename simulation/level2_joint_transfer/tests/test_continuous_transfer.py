from pathlib import Path
import numpy as np
import pytest
from continuous_transfer.engine import field, spin_step, run_continuous
from continuous_transfer.protocol import build_protocol, waveform_for, ideal_checkpoint_unitaries
from continuous_transfer.cli import ROOT, case_inputs, quantum_summary
from level2c_pickup_survival.survival_engine import DisorderDraws, build_roundtrip_legs, run_survival
from level3_3d_transport_lensing.aod_lensing import lensing_potential_and_force
from level5_finite_pulse_irb.quantum_propagators import su2_propagator
import yaml


@pytest.fixture
def inputs():
    settings = yaml.safe_load((ROOT/"configs/continuous_level2c_345.yaml").read_text())
    return settings, case_inputs(settings,40,16,23002,"baseline_1d")


def test_accelerated_force_is_level3_force():
    rng = np.random.default_rng(8)
    for _ in range(50):
        x,y,z = rng.normal(size=3)*2e-6
        d,w,zr = 280*1.380649e-29,1.17e-6,4.07e-6
        s1,s2 = rng.normal(size=2)*zr
        actual = np.array(field(x,y,z,d,w,zr,s1,s2))
        expected = np.array(lensing_potential_and_force(x,y,z,d,w,zr,s1,s2))
        np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=1e-40)
        h = 1e-11
        fd = -(field(x+h,y,z,d,w,zr,s1,s2)[0]-field(x-h,y,z,d,w,zr,s1,s2)[0])/(2*h)
        np.testing.assert_allclose(actual[1],fd,rtol=2e-8,atol=1e-35)


def test_spin_kernel_matches_level5_and_preserves_unitarity():
    rng = np.random.default_rng(3)
    a,b = 1+0j,0+0j
    expected = np.eye(2,dtype=complex)
    for _ in range(100):
        delta,omega,phi,dt = rng.normal()*2000,24611*2*np.pi,rng.uniform(0,6.3),1e-7
        a,b = spin_step(a,b,delta,omega,phi,dt)
        expected = np.asarray(su2_propagator(delta,omega,phi,dt)).reshape(2,2) @ expected
    np.testing.assert_allclose([a,b],expected[0],atol=1e-13)
    assert abs(abs(a)**2+abs(b)**2-1) < 1e-13


@pytest.mark.parametrize("arm",["matched","extended"])
def test_protocol_switches_count_and_exact_finite_pulse_edges(inputs,arm):
    settings,(cfg,*_) = inputs
    wave = waveform_for("manual_400us",cfg,ROOT/settings["ml_waveform"])
    protocol = build_protocol(wave,100e-6,.05e-6,arm=arm)
    assert protocol.duration_s == pytest.approx((900 if arm=="matched" else 3854)*1e-6)
    assert protocol.distance_per_round_m == pytest.approx((4.8 if arm=="matched" else 754.8)*1e-6)
    assert len(protocol.pulses) == (4 if arm=="matched" else 8)
    grid = protocol.controls
    for pulse in protocol.pulses:
        assert np.min(abs(grid[:,0]-pulse.start)) < 1e-18
        assert np.min(abs(grid[:,0]-(pulse.start+pulse.duration))) < 1e-18
    for i,(start,end,_,_) in enumerate(protocol.segments):
        assert np.all(np.diff(grid[start:end+1,0])>0)
        assert grid[start,6] == (1 if protocol.names[i] in ("pickup","dropoff") else 0)
    ideal = ideal_checkpoint_unitaries(protocol,30)
    np.testing.assert_allclose(abs(np.trace(ideal[-1])),2,atol=1e-12)


def test_1d_limit_matches_original_level2c(inputs):
    settings,(cfg,physics,pos,vel,_,_,_) = inputs
    wave = waveform_for("manual_400us",cfg,ROOT/settings["ml_waveform"])
    protocol = build_protocol(wave,100e-6,.05e-6,dd=False)
    draws = DisorderDraws.disabled(len(pos))
    new = run_continuous(pos,vel,protocol,rounds=3,physics=physics,draws=draws,
                         jitter_sigma_m=0,power_w=0,vs_m_s=0,eta=0)
    old = run_survival(pos[:,0],vel[:,0],physics["mass_kg"],
        build_roundtrip_legs(wave,100e-6,.05e-6),.05e-6,6,
        physics["slm_depth_j"],physics["slm_waist_m"],physics["aod_waist_m"],disorder=draws,
        record_transfers=(6,))
    np.testing.assert_array_equal(new["alive"].sum(axis=0),old["survival_counts"])
    snapshot = old["snapshots"][6]
    ids = snapshot["shot_ids"]
    kept = new["alive"][ids,6]
    np.testing.assert_allclose(new["states"][ids[kept],6,0],snapshot["x_m"][kept],rtol=1e-5,atol=1e-12)
    np.testing.assert_allclose(new["states"][ids[kept],6,3],snapshot["v_m_s"][kept],rtol=1e-5,atol=1e-8)
    # An independently continued second round equals one continuous two-round run.
    first = run_continuous(pos,vel,protocol,rounds=1,physics=physics,draws=draws,
                           jitter_sigma_m=0,power_w=0,vs_m_s=0,eta=0)
    keep = first["alive"][:,-1]
    second = run_continuous(first["final_states"][keep,:3],first["final_states"][keep,3:],
        protocol,rounds=1,physics=physics,draws=DisorderDraws.disabled(int(keep.sum())),
        jitter_sigma_m=0,power_w=0,vs_m_s=0,eta=0)
    new2 = run_continuous(pos,vel,protocol,rounds=2,physics=physics,draws=draws,
                          jitter_sigma_m=0,power_w=0,vs_m_s=0,eta=0)
    np.testing.assert_allclose(second["final_states"],new2["final_states"][keep],rtol=1e-10,atol=1e-13)


def test_absorbing_loss_heat_accounting_and_no_spin_feedback(inputs):
    settings,(cfg,physics,pos,vel,draws,_,_) = inputs
    pos = np.zeros_like(pos)
    vel = np.zeros_like(vel)
    vel[0,0] = 10.0  # deliberately fail the first pickup
    wave = waveform_for("manual_400us",cfg,ROOT/settings["ml_waveform"])
    protocol = build_protocol(wave,100e-6,.05e-6)
    kwargs = dict(rounds=2,physics=physics,draws=draws,jitter_sigma_m=1e-8,power_w=1e-27,vs_m_s=.539)
    result = run_continuous(pos,vel,protocol,**kwargs)
    zero = run_continuous(pos,vel,protocol,eta=0,**kwargs)
    assert result["alive"][0].tolist()==[True,False,False,False,False]
    assert np.all(np.diff(result["alive"].astype(int),axis=1)<=0)
    np.testing.assert_allclose(result["injected_j"],result["elapsed_s"]*kwargs["power_w"],rtol=1e-10)
    assert result["elapsed_s"][0] < result["elapsed_s"][1]
    np.testing.assert_array_equal(result["alive"],zero["alive"])
    np.testing.assert_array_equal(result["final_states"],zero["final_states"])
    summary = quantum_summary(zero,protocol,2)
    for row in summary["checkpoints"]:
        if row["survivors"]:
            assert row["favg_dd"] == pytest.approx(1,abs=1e-10)
    np.testing.assert_allclose(summary["final_conditional_ptm_relative_to_ideal"],np.eye(4),atol=1e-10)


def test_continuous_spin_matches_direct_level5_replay(inputs):
    settings,(cfg,physics,_,_,_,_,_) = inputs
    wave = waveform_for("manual_400us",cfg,ROOT/settings["ml_waveform"])
    protocol = build_protocol(wave,100e-6,.05e-6)
    # A static overlapping trap makes delta constant while finite pulses retain
    # their actual schedule. Direct Level 5 propagation checks every substep.
    protocol.controls[:,1:5] = 0
    protocol.controls[:,5] = physics["slm_depth_j"]
    protocol.controls[:,6] = 1
    eta = 1.3e-4
    result = run_continuous(np.zeros((1,3)),np.zeros((1,3)),protocol,rounds=2,
        physics=physics,draws=DisorderDraws.disabled(1),jitter_sigma_m=0,power_w=0,vs_m_s=0,eta=eta)
    expected = np.eye(2,dtype=complex)
    delta = -2*physics["slm_depth_j"]*eta/1.054571817e-34
    for _ in range(2):
        for start,end,_,_ in protocol.segments:
            for i in range(start,end):
                c = protocol.controls[i]
                dt = protocol.controls[i+1,0]-c[0]
                expected = np.asarray(su2_propagator(delta,c[7],c[8],dt)).reshape(2,2) @ expected
    np.testing.assert_allclose(result["spin_ab"][0,-1],expected[0],atol=1e-10)


def test_ml_is_unmodified_stage1_digitized_replay(inputs):
    settings,(cfg,*_) = inputs
    import json
    raw = json.loads((ROOT/settings["ml_waveform"]).read_text())
    wave = waveform_for("ml_400us",cfg,ROOT/settings["ml_waveform"])
    t = np.array(raw["grid_us"])*1e-6
    np.testing.assert_allclose(wave.center(t),np.array(raw["pos_um"])*1e-6,atol=1e-20)
    np.testing.assert_allclose(wave.depth(t),np.array(raw["depth_mK"])*1.380649e-26,atol=1e-40)
