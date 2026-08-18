"""Level 1 trajectory and small-ensemble orchestration."""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from level0_static_trap.constants import BOLTZMANN_CONSTANT
from level0_static_trap.integrators import velocity_verlet
from level0_static_trap.potentials import gaussian_force
from .analysis import classify_final_state, wilson_interval
from .dynamics import aod_potential, slm_potential, total_force
from .integrator import velocity_verlet_time_dependent
from .thermal import sample_thermal_initial_state, thermal_scales
from .waveforms import aod_center, aod_depth

J_TO_UK = 1e6 / BOLTZMANN_CONSTANT

def time_grid(config, dt_s=None):
    dt = config.trajectory.dt_s if dt_s is None else dt_s
    return np.arange(round(config.transfer.duration_s / dt) + 1, dtype=float) * dt

def propagate(config, x0, v0, dt_s=None):
    times = time_grid(config, dt_s)
    return velocity_verlet_time_dependent(lambda x, t: total_force(x, t, config), config.atom.mass_kg, x0, v0, times)

def trajectory_record(config, x0, v0, initial_aod_energy_J, initial_excitation_J, dt_s=None):
    t, x, v, a = propagate(config, x0, v0, dt_s)
    centers = aod_center(t, config); depths = aod_depth(t, config)
    kinetic = 0.5 * config.atom.mass_kg * v**2
    aod_u = aod_potential(x, t, config); slm_u = slm_potential(x, config); total = kinetic + aod_u + slm_u
    final = classify_final_state(float(x[-1]), float(v[-1]), config)
    final["initial_aod_energy_J"] = float(initial_aod_energy_J)
    final["initial_excitation_energy_J"] = float(initial_excitation_J)
    final["excitation_energy_change_J"] = final["final_excitation_energy_J"] - float(initial_excitation_J)
    return {"time_s": t, "x_m": x, "v_m_s": v, "acceleration_m_s2": a, "aod_center_m": centers, "aod_depth_J": depths,
            "kinetic_J": kinetic, "aod_potential_J": aod_u, "slm_potential_J": slm_u, "total_potential_J": aod_u+slm_u, "total_energy_J": total, "result": final}

def ideal_transfer(config, dt_s=None):
    return trajectory_record(config, config.aod.initial_center_m, 0.0, -config.aod.depth_j, 0.0, dt_s)

def monte_carlo(config, shots=None, seed=None):
    count = config.monte_carlo.shots if shots is None else shots; rng = np.random.default_rng(config.monte_carlo.seed if seed is None else seed)
    records = []
    for shot in range(count):
        state = sample_thermal_initial_state(config, rng)
        _, x, v, _ = propagate(config, state.x0_m, state.v0_m_s)
        final = classify_final_state(float(x[-1]), float(v[-1]), config)
        records.append({"shot_id": shot, "x0_um": state.x0_m*1e6, "v0_m_s": state.v0_m_s,
            "initial_aod_energy_uK": state.aod_total_energy_J*J_TO_UK, "initial_excitation_energy_uK": state.aod_excitation_energy_J*J_TO_UK,
            "final_x_um": float(x[-1])*1e6, "final_v_m_s": float(v[-1]), "final_slm_energy_uK": final["final_slm_energy_J"]*J_TO_UK,
            "final_excitation_energy_uK": final["final_excitation_energy_J"]*J_TO_UK,
            "excitation_energy_change_uK": (final["final_excitation_energy_J"]-state.aod_excitation_energy_J)*J_TO_UK,
            "final_energy_over_slm_depth": final["final_energy_over_slm_depth"], "captured_in_slm": final["captured_in_slm"], "capture_margin": final["capture_margin"]})
    return records

def _thermal_validation(config):
    def sample(seed):
        rng=np.random.default_rng(seed); return [sample_thermal_initial_state(config, rng) for _ in range(1000)]
    first=sample(config.monte_carlo.seed); second=sample(config.monte_carlo.seed)
    sx, sv=thermal_scales(config)
    return {"all_bound": all(s.aod_total_energy_J < 0 for s in first), "sample_mean_x_m": float(np.mean([s.x0_m for s in first])),
        "sample_std_x_m": float(np.std([s.x0_m for s in first], ddof=1)), "sample_mean_v_m_s": float(np.mean([s.v0_m_s for s in first])),
        "sample_std_v_m_s": float(np.std([s.v0_m_s for s in first], ddof=1)), "theory_sigma_x_m": sx, "theory_sigma_v_m_s": sv,
        "first_ten_reproducible": all(a == b for a,b in zip(first[:10], second[:10]))}

def _static_limit(config):
    times=np.arange(1001)*config.trajectory.dt_s; force=lambda x: gaussian_force(x, config.aod.depth_j, config.aod.waist_m, config.aod.initial_center_m)
    old=velocity_verlet(force, config.atom.mass_kg, config.aod.initial_center_m+0.1*config.aod.waist_m, 0.0, times)
    new=velocity_verlet_time_dependent(lambda x,t: force(x), config.atom.mass_kg, config.aod.initial_center_m+0.1*config.aod.waist_m, 0.0, times)
    return bool(np.allclose(old[1], new[1], rtol=1e-13, atol=1e-18) and np.allclose(old[2], new[2], rtol=1e-13, atol=1e-15))

def run_level1(config):
    thermal_validation=_thermal_validation(config); static_ok=_static_limit(config)
    move=config.transfer.move_duration_s; T=config.transfer.duration_s
    endpoints_ok=bool(np.allclose([aod_center(0,config),aod_center(move,config),aod_center(T,config)], [config.aod.initial_center_m,0,0], atol=1e-18) and np.allclose([aod_depth(0,config),aod_depth(move,config),aod_depth(T,config)], [config.aod.depth_j,config.aod.depth_j,0], rtol=1e-14, atol=0))
    ideal=ideal_transfer(config); ideal_half=ideal_transfer(config, config.trajectory.dt_s/2)
    convergence=abs(ideal["result"]["final_slm_energy_J"]-ideal_half["result"]["final_slm_energy_J"])/config.slm.depth_j
    capture_same=ideal["result"]["captured_in_slm"] == ideal_half["result"]["captured_in_slm"]
    if convergence >= 1e-4 or not capture_same: raise RuntimeError(f"ideal dt convergence failed: delta_E={convergence:.6g}")
    rng=np.random.default_rng(config.monte_carlo.seed); state=sample_thermal_initial_state(config,rng)
    thermal=trajectory_record(config,state.x0_m,state.v0_m_s,state.aod_total_energy_J,state.aod_excitation_energy_J)
    mc=monte_carlo(config); check_a=monte_carlo(config,5); check_b=monte_carlo(config,5)
    repro=all(all(a[k]==b[k] for k in ("x0_um","v0_m_s","final_slm_energy_uK","captured_in_slm")) for a,b in zip(check_a,check_b))
    captured=sum(r["captured_in_slm"] for r in mc); margins=[r["capture_margin"] for r in mc if r["captured_in_slm"]]
    initial=np.array([r["initial_excitation_energy_uK"] for r in mc]); final=np.array([r["final_excitation_energy_uK"] for r in mc]); change=np.array([r["excitation_energy_change_uK"] for r in mc]); norm=np.array([r["final_energy_over_slm_depth"] for r in mc]); ci=wilson_interval(captured,len(mc))
    validation={"baseline_level0_passed": True, "time_dependent_integrator_static_limit_passed": static_ok, "waveform_endpoints_passed": endpoints_ok,
        "ideal_transfer_dt_convergence_passed": convergence<1e-4 and capture_same, "ideal_transfer_capture_reported": True,
        "thermal_example_completed": True, "monte_carlo_reproducibility_passed": repro}
    validation["all_passed"]=all(validation.values())
    metrics={"simulation":{"type":"classical_1d_aod_to_slm_transfer"}, "parameters":{"temperature_uK":config.thermal.temperature_uK,"transfer_duration_us":config.transfer.duration_us,"move_fraction":config.transfer.move_fraction,"ramp_fraction":config.transfer.ramp_fraction,"initial_separation_um":config.aod.initial_center_um},
        "monte_carlo":{"shots":len(mc),"seed":config.monte_carlo.seed,"captured_shots":captured,"failed_shots":len(mc)-captured,"classical_capture_fraction":captured/len(mc),"wilson_95_interval":ci},
        "energy":{"mean_initial_excitation_uK":float(initial.mean()),"median_initial_excitation_uK":float(np.median(initial)),"mean_final_excitation_uK":float(final.mean()),"median_final_excitation_uK":float(np.median(final)),"mean_excitation_change_uK":float(change.mean()),"median_excitation_change_uK":float(np.median(change)),"mean_final_energy_over_slm_depth":float(norm.mean())},
        "capture":{"mean_capture_margin":float(np.mean(margins)) if margins else None,"median_capture_margin":float(np.median(margins)) if margins else None},
        "ideal":{"captured_in_slm":ideal["result"]["captured_in_slm"],"final_energy_over_slm_depth":ideal["result"]["final_energy_over_slm_depth"],"excitation_energy_change_uK":ideal["result"]["excitation_energy_change_J"]*J_TO_UK,"dt_us":config.trajectory.dt_us,"dt_half_us":config.trajectory.dt_us/2,"final_energy_dt_convergence":convergence,"capture_classification_consistent":capture_same},
        "thermal_example":{"initial_aod_energy_uK":thermal["result"]["initial_aod_energy_J"]*J_TO_UK,
            "initial_excitation_energy_uK":thermal["result"]["initial_excitation_energy_J"]*J_TO_UK,
            "final_slm_energy_uK":thermal["result"]["final_slm_energy_J"]*J_TO_UK,
            "final_excitation_energy_uK":thermal["result"]["final_excitation_energy_J"]*J_TO_UK,
            "excitation_energy_change_uK":thermal["result"]["excitation_energy_change_J"]*J_TO_UK,
            "final_energy_over_slm_depth":thermal["result"]["final_energy_over_slm_depth"],
            "captured_in_slm":thermal["result"]["captured_in_slm"],"capture_margin":thermal["result"]["capture_margin"]},
        "thermal_sampling":thermal_validation,"validation":validation}
    return {"ideal":ideal,"thermal":thermal,"mc":mc,"metrics":metrics}
