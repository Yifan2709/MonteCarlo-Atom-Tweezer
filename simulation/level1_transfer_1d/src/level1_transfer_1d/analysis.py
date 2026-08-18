"""Capture and excitation diagnostics."""
from level0_static_trap.potentials import gaussian_potential
def classify_final_state(x_m, v_m_s, config):
    energy = 0.5*config.atom.mass_kg*v_m_s**2 + float(gaussian_potential(x_m, config.slm.depth_j, config.slm.waist_m, config.slm.center_m)); captured = energy < 0
    return {"final_slm_energy_J": energy, "final_excitation_energy_J": energy+config.slm.depth_j, "final_energy_over_slm_depth": energy/config.slm.depth_j, "captured_in_slm": bool(captured), "capture_margin": -energy/config.slm.depth_j if captured else None}
def wilson_interval(successes, total, z=1.959963984540054):
    p=successes/total; d=1+z*z/total; c=(p+z*z/(2*total))/d; r=z*((p*(1-p)/total+z*z/(4*total*total))**0.5)/d
    return c-r, c+r
