"""Command-line entry point and artifact writer."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from level0_static_trap.io_utils import ensure_output_directory, save_csv, save_json, save_markdown, save_yaml
from .config import config_as_dict, load_config, with_output_override
from .simulation import J_TO_UK, run_level1
from .visualization import plot_distribution, plot_energy, plot_phase, plot_position, plot_potential_snapshot, plot_waveform
from .waveforms import aod_center, aod_depth

def _trajectory_columns(r,c):
    return {"time_us":r["time_s"]*1e6,"x_um":r["x_m"]*1e6,"v_m_s":r["v_m_s"],"aod_center_um":r["aod_center_m"]*1e6,
        "aod_depth_uK":r["aod_depth_J"]*J_TO_UK,"slm_depth_uK":np.full_like(r["time_s"],c.slm.depth_uK),"kinetic_uK":r["kinetic_J"]*J_TO_UK,
        "aod_potential_uK":r["aod_potential_J"]*J_TO_UK,"slm_potential_uK":r["slm_potential_J"]*J_TO_UK,"total_potential_uK":r["total_potential_J"]*J_TO_UK,"total_energy_uK":r["total_energy_J"]*J_TO_UK}
def _summary(c,m):
    ideal=m["ideal"]; thermal=m["thermal_example"]; mc=m["monte_carlo"]; cap=m["capture"]; energy=m["energy"]; val=m["validation"]
    status=lambda x:"PASS" if x else "FAIL"
    return f"""# Level 1: 1D classical AOD -> SLM transfer

This is a finite-temperature classical mechanical model, not experimental survival, transfer fidelity, or quantum fidelity. It is not directly comparable to the literature 99.81% benchmark.

## Q1 - Time-dependent integrator static limit

{status(val['time_dependent_integrator_static_limit_passed'])}

## Q2 - Waveform

- Initial separation: {c.aod.initial_center_um:.3f} um
- Move duration: {c.transfer.duration_us*c.transfer.move_fraction:.3f} us
- Ramp duration: {c.transfer.duration_us*c.transfer.ramp_fraction:.3f} us
- AOD depth: {c.aod.depth_uK:.1f} -> 0.0 uK
- SLM depth: {c.slm.depth_uK:.1f} uK

## Q3 - Ideal atom

- Final energy / SLM depth: {ideal['final_energy_over_slm_depth']:.8f}
- Captured: {ideal['captured_in_slm']}

## Q4 - Timestep convergence

- dt / dt-half: {ideal['dt_us']:.3f} / {ideal['dt_half_us']:.3f} us
- Normalized final-energy difference: {ideal['final_energy_dt_convergence']:.3e}
- Capture classification consistent: {ideal['capture_classification_consistent']}

## Q5 - Thermal example

- Initial AOD energy: {thermal['initial_aod_energy_uK']:.6f} uK
- Final SLM energy: {thermal['final_slm_energy_uK']:.6f} uK
- Excitation-energy change: {thermal['excitation_energy_change_uK']:.6f} uK
- Captured: {thermal['captured_in_slm']}

## Q6 - Classical capture fraction

- Captured / total: {mc['captured_shots']} / {mc['shots']}
- classical_capture_fraction: {mc['classical_capture_fraction']:.6f}
- Wilson 95% interval: [{mc['wilson_95_interval'][0]:.6f}, {mc['wilson_95_interval'][1]:.6f}]

## Q7 - Capture margin

- Mean / median capture margin: {cap['mean_capture_margin']} / {cap['median_capture_margin']}
- See `final_slm_energy_distribution.png`.

## Q8 - Classical motional excitation proxy

- Mean excitation-energy change: {energy['mean_excitation_change_uK']:.6f} uK
- Median excitation-energy change: {energy['median_excitation_change_uK']:.6f} uK

## Validation

All passed: {val['all_passed']}
"""
def write_outputs(c,result):
    out=ensure_output_directory(c.output.resolved_directory,c.project_root); save_yaml(out/"config_used.yaml",config_as_dict(c)); ideal=result["ideal"]; thermal=result["thermal"]; mc=result["mc"]
    t=ideal["time_s"]; save_csv(out/"transfer_waveform.csv",{"time_us":t*1e6,"aod_center_um":aod_center(t,c)*1e6,"aod_depth_uK":aod_depth(t,c)*J_TO_UK,"slm_depth_uK":np.full_like(t,c.slm.depth_uK)})
    save_csv(out/"ideal_transfer_trajectory.csv",_trajectory_columns(ideal,c)); save_csv(out/"thermal_example_trajectory.csv",_trajectory_columns(thermal,c))
    keys=list(mc[0]); save_csv(out/"mc_shots.csv",{k:[r[k] if r[k] is not None else "" for r in mc] for k in keys})
    plot_waveform(t,c,out/"transfer_waveform.png")
    for time,name,title in ((0,"potential_start.png","Start"),(c.transfer.move_duration_s,"potential_merged.png","Merged"),(c.transfer.duration_s,"potential_final.png","Final SLM-only")): plot_potential_snapshot(c,time,out/name,title)
    plot_position(ideal,c,out/"ideal_transfer_position.png","Ideal transfer position"); plot_energy(ideal,out/"ideal_transfer_energy.png","Ideal mechanical energy (time dependent)")
    plot_position(thermal,c,out/"thermal_example_position.png","Thermal example position"); plot_energy(thermal,out/"thermal_example_energy.png","Thermal example mechanical energy")
    plot_phase(mc,out/"initial_phase_space.png",True); plot_phase(mc,out/"final_phase_space.png",False)
    plot_distribution([r["final_energy_over_slm_depth"] for r in mc],out/"final_slm_energy_distribution.png","Final SLM energy / depth",0)
    plot_distribution([r["excitation_energy_change_uK"] for r in mc],out/"excitation_energy_change_distribution.png","Excitation-energy change (uK)")
    save_json(out/"metrics.json",result["metrics"]); save_markdown(out/"summary.md",_summary(c,result["metrics"])); return out
def main(argv=None):
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); p.add_argument("--output-dir"); args=p.parse_args(argv)
    c=load_config(Path(args.config)); c=with_output_override(c,args.output_dir) if args.output_dir else c; result=run_level1(c); out=write_outputs(c,result)
    print(f"Level 1 completed: {out}"); print(f"classical_capture_fraction={result['metrics']['monte_carlo']['classical_capture_fraction']:.6f}"); return 0
if __name__ == "__main__": raise SystemExit(main())
