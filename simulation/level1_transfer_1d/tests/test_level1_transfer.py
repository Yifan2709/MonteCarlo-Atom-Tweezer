from pathlib import Path
import numpy as np
from level0_static_trap.integrators import velocity_verlet
from level0_static_trap.potentials import gaussian_force
from level1_transfer_1d.analysis import classify_final_state
from level1_transfer_1d.config import load_config
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level1_transfer_1d.simulation import monte_carlo
from level1_transfer_1d.waveforms import aod_center, aod_depth

CONFIG=Path(__file__).resolve().parents[1]/"configs"/"transfer_1d.yaml"
def test_time_dependent_integrator_static_limit():
    c=load_config(CONFIG); times=np.arange(1001)*c.trajectory.dt_s; force=lambda x: gaussian_force(x,c.aod.depth_j,c.aod.waist_m,c.aod.initial_center_m)
    old=velocity_verlet(force,c.atom.mass_kg,c.aod.initial_center_m+.1*c.aod.waist_m,0,times); new=velocity_verlet_time_dependent(lambda x,t:force(x),c.atom.mass_kg,c.aod.initial_center_m+.1*c.aod.waist_m,0,times)
    assert np.allclose(old[1],new[1],rtol=1e-13,atol=1e-18); assert np.allclose(old[2],new[2],rtol=1e-13,atol=1e-15)
def test_waveform_endpoints():
    c=load_config(CONFIG); m=c.transfer.move_duration_s; T=c.transfer.duration_s
    assert np.allclose([aod_center(0,c),aod_center(m,c),aod_center(T,c)],[c.aod.initial_center_m,0,0],atol=1e-18)
    assert np.allclose([aod_depth(0,c),aod_depth(m,c),aod_depth(T,c)],[c.aod.depth_j,c.aod.depth_j,0],rtol=1e-14,atol=0)
def test_final_capture_classifier():
    c=load_config(CONFIG); assert classify_final_state(c.slm.center_m,0,c)["captured_in_slm"]
    assert not classify_final_state(c.slm.center_m,10,c)["captured_in_slm"]
def test_five_shot_monte_carlo_reproducibility():
    c=load_config(CONFIG); a=monte_carlo(c,5); b=monte_carlo(c,5)
    for left,right in zip(a,b):
        for key in ("x0_um","v0_m_s","final_slm_energy_uK","captured_in_slm"): assert left[key]==right[key]
