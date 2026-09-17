"""Independent mechanical/noise checks for the corrected Fig6d entry."""
import numpy as np
import pytest
from continuous_transfer.engine import field,run_continuous
from continuous_transfer.noisy_survival import beam_noise,run_noisy_survival
from continuous_transfer.paper_waveforms import PaperManualWaveform
from continuous_transfer.protocol import build_protocol
from level2c_pickup_survival.survival_engine import DisorderDraws,time_reverse_waveform


def test_pointing_force_derivative_including_axial_coupling():
    rng=np.random.default_rng(26)
    for _ in range(30):
        xyz=rng.normal(size=3)*1e-6
        args=(280*1.380649e-29,1.17e-6,4.1e-6,.2e-6,-.3e-6)
        analytic=np.array(beam_noise(*xyz,*args)[4:])
        h=1e-11
        plus=xyz.copy();plus[0]+=h
        minus=xyz.copy();minus[0]-=h
        fd=(np.array(field(*plus,*args)[1:])-np.array(field(*minus,*args)[1:]))/(2*h)
        np.testing.assert_allclose(analytic,fd,rtol=2e-7,atol=1e-32)


def test_paper_cubic_constant_jerk_and_exact_reverse():
    w=PaperManualWaveform(400e-6,.48,2.4e-6,280*1.380649e-29)
    t=np.linspace(0,w.duration_s,201)
    reverse=time_reverse_waveform(w)
    np.testing.assert_allclose(reverse.center(t),w.center(w.duration_s-t),atol=2e-21)
    np.testing.assert_allclose(reverse.depth(t),w.depth(w.duration_s-t),atol=1e-41)
    u=np.linspace(.01,.99,101)
    np.testing.assert_allclose(w.center(w.ramp_duration_s+u*w.move_duration_s)/2.4e-6,3*u*u-2*u**3)
    assert w.center_velocity(0)==0 and w.center_velocity(w.duration_s)==0
    assert np.ptp(w.center_jerk(w.ramp_duration_s+u*w.move_duration_s))==0


def test_zero_diffusion_matches_continuous_mechanics_and_stage_bookkeeping():
    physics=dict(mass_kg=132.90545196*1.66053906660e-27,slm_depth_j=140*1.380649e-29,
                 slm_waist_m=1.17e-6,aod_waist_m=1.17e-6)
    p=np.random.default_rng(6).normal(size=(12,3))*1e-7
    v=np.zeros_like(p);draws=DisorderDraws.disabled(len(p))
    wave=PaperManualWaveform(400e-6,.48,2.4e-6,280*1.380649e-29)
    protocol=build_protocol(wave,100e-6,.025e-6,dd=False)
    kw=dict(rounds=2,physics=physics,draws=draws,jitter_sigma_m=0,vs_m_s=0)
    a=run_continuous(p,v,protocol,power_w=0,eta=0,**kw)
    b=run_noisy_survival(p,v,protocol,mode='diffusion',rin=0,pointing=0,**kw)
    np.testing.assert_array_equal(a['alive'],b['alive'])
    np.testing.assert_allclose(a['states'],b['states'],rtol=1e-10,atol=1e-14,equal_nan=True)
    assert np.all(b['noise_ledger'][:,:2]==0)
    # SLM-off adds positive potential work, SLM-on removes it; no state reset.
    assert np.all(b['switch_work_j'][:,:,1]>=0)
    assert np.all(b['switch_work_j'][:,:,2]<=0)
    assert np.isfinite(b['stage_energy_j'][:,:,0]).all()


def test_noise_streams_are_reproducible_across_thread_counts():
    import numba
    physics=dict(mass_kg=132.90545196*1.66053906660e-27,slm_depth_j=140*1.380649e-29,
                 slm_waist_m=1.17e-6,aod_waist_m=1.17e-6)
    p=np.random.default_rng(6).normal(size=(24,3))*1e-7
    wave=PaperManualWaveform(400e-6,.48,2.4e-6,280*1.380649e-29)
    protocol=build_protocol(wave,100e-6,.05e-6,dd=False)
    kw=dict(rounds=2,physics=physics,draws=DisorderDraws.disabled(len(p)),
            jitter_sigma_m=0,vs_m_s=0,rin=1e-8,pointing=1e-24,slm_rin=0,slm_pointing=0,seed=918)
    previous=numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        a=run_noisy_survival(p,np.zeros_like(p),protocol,**kw)
        numba.set_num_threads(min(4,previous))
        b=run_noisy_survival(p,np.zeros_like(p),protocol,**kw)
    finally:
        numba.set_num_threads(previous)
    for key in a:np.testing.assert_array_equal(a[key],b[key])
