"""Regression guards for confirmed physical errors, not curve-shape criteria."""
import numpy as np
from paper_integration.dynamics3d import TrapPhysics,sample_bound_thermal_3d,Segment,propagate
from paper_integration.trajectories import poly_trajectory
from level3_3d_transport_lensing.gaussian_3d import static_potential
from long_distance_transport.model import Config,run,KB

def test_actual_site_all_initially_bound_and_no_round_reset():
    c=Config(shots=64,protocol='B',distance_um=0,repeats=3,disorder=.05,temperature_uK=30)
    _,r=run(c)
    assert (r['energy_j'][:,0]<0).all()
    assert (r['attempts']>=1).all()
    # A perfect rest is a continuous orbit, not a re-drawn cloud.
    c=Config(shots=8,distance_um=0,duration_us=200,hold_us=100,repeats=2)
    _,a=run(c)
    b=Config(shots=8,distance_um=0,duration_us=500,hold_us=100,repeats=1)
    _,b=run(b)
    np.testing.assert_allclose(a['states'][:,-1],b['states'][:,-1],atol=2e-10)

def test_legacy_sampler_retries_bound_states_without_boolean_shape_error():
    ph=TrapPhysics(133*1.66054e-27,50e-6*KB,1e-6,3e-6)
    p,v=sample_bound_thermal_3d(7,70,127,ph)
    e=.5*ph.mass_kg*np.sum(v*v,axis=1)+static_potential(*p.T,ph.depth_j,ph.waist_m,ph.zr_m)
    assert np.all(e<0)

def test_legacy_uniform_translation_uses_trap_frame():
    ph=TrapPhysics(133*1.66054e-27,280e-6*KB,1e-6,3e-6)
    t=np.linspace(0,100e-6,1001);V=.7
    s=Segment('uniform',t,t*V,np.zeros_like(t),np.full_like(t,ph.depth_j),np.full_like(t,V))
    r=propagate(np.zeros((2,3)),np.tile([V,0,0],(2,1)),[s],ph,999.)
    assert r['alive'].all()
    np.testing.assert_allclose(r['final_energy_j'],-ph.depth_j,rtol=1e-10)

def test_exact_trajectory_time_and_connection_semantics():
    t=poly_trajectory('zero_jerk',2,3,.7)
    assert t.t_s[-1]==3 and abs(t.x_m[-1]-2)<1e-12
    assert abs(t.j_m_s3[0])<1e-12 and abs(t.j_m_s3[-1])<1e-12
    assert t.a_m_s2[0]>0 and t.a_m_s2[-1]<0
    t=poly_trajectory('constant_accel',1,1,.0001)
    assert np.max(np.diff(t.x_m))<.001
