"""Reproducible Level 2C→3/4/5 continuous simulation and reference acceptance."""
from __future__ import annotations
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import yaml
from level2c_pickup_survival.level2c_config import load_level2c_config
from level2c_pickup_survival.survival_engine import draw_disorder
from level2c_pickup_survival.reference_data import load_reference_curves, compare_to_reference
from level2_joint_transfer.level2_simulation import physics_from_config, sample_initial_states
from level2_joint_transfer.noise_heating import heating_from_config
from level3_3d_transport_lensing.thermal_sampling_3d import sample_initial_states_3d
from .engine import ACCELERATED, run_continuous
from .protocol import waveform_for, build_protocol, ideal_checkpoint_unitaries
from .initialization import sample_site_bound_harmonic, slm_energies

ROOT = Path(__file__).resolve().parents[2]
ARMS = ("baseline_1d", "matched_lensing_off", "matched", "extended")


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path, value):
    temporary = path.with_suffix(".partial.json")
    temporary.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def wilson(success, n):
    p = np.asarray(success) / n
    z = 1.959963984540054
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return np.maximum(0,center-half), np.minimum(1,center+half)


def quantum_summary(result, protocol, rounds):
    targets = ideal_checkpoint_unitaries(protocol, rounds)
    rows = []
    # PTM is the conditional, trace-preserving channel on surviving shots.
    pauli = np.array([np.eye(2), [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], dtype=complex)
    final_ptm = None
    for n, target in enumerate(targets):
        mask = result["alive"][:,n]
        if not mask.any():
            rows.append({"n":n,"survivors":0,"contrast_none":None,"contrast_dd":None,
                         "favg_none":None,"favg_dd":None})
            continue
        a,b = result["spin_ab"][mask,n].T
        u = np.stack((a,b,-b.conj(),a.conj()),axis=-1).reshape(-1,2,2)
        relative = np.einsum("ij,bjk->bik",target.conj().T,u)
        trace = np.trace(relative,axis1=1,axis2=2)
        plus = np.einsum("bij,j->bi",relative,np.ones(2)/np.sqrt(2))
        coherence = 2*np.mean(plus[:,0]*plus[:,1].conj())
        phases = result["phase_rad"][mask,n]
        rows.append({"n":n,"survivors":int(mask.sum()),
                     "contrast_none":float(abs(np.mean(np.exp(1j*phases)))),
                     "contrast_dd":float(abs(coherence)),
                     "favg_none":float(np.mean((4*np.cos(phases/2)**2+2)/6)),
                     "favg_dd":float(np.mean((np.abs(trace)**2+2)/6)),
                     "unitarity_error":float(np.max(abs(abs(a)**2+abs(b)**2-1)))})
        if n == 2*rounds:
            final_ptm = np.array([[np.mean(np.trace(
                pauli[i] @ relative @ pauli[j] @ relative.conj().transpose(0,2,1),
                axis1=1,axis2=2)).real/2 for j in range(4)] for i in range(4)])
    return {"checkpoints":rows,"final_conditional_ptm_relative_to_ideal":final_ptm,
            "definition":"Favg=(mean |Tr(Uideal† U)|²+2)/6; survivors only; no RB ratio fit",
            "feedback":"No spin-dependent force, so none/DD have identical classical survival"}


def case_inputs(settings, temperature, shots, seed, arm):
    cfg = load_level2c_config(ROOT / settings["level2c_config"])
    cfg = replace(cfg, base=replace(cfg.base, initial_ensemble=replace(
        cfg.base.initial_ensemble, temperature_uK=temperature)))
    if settings.get("parametric_reference_temperature_uK") is not None:
        tref = float(settings["parametric_reference_temperature_uK"])
        if not np.isfinite(tref) or tref < 0:
            raise ValueError("Heating reference temperature must be finite and nonnegative")
        cfg = replace(cfg, base=replace(cfg.base, noise_mean_heating=replace(
            cfg.base.noise_mean_heating, parametric_reference_temperature_uK=tref)))
    physics = physics_from_config(cfg.base)
    if physics["slm_center_m"] != 0:
        raise ValueError("This bridge requires the Level 2C SLM center at zero")
    draws = draw_disorder(cfg.disorder,shots,seed+7)
    initial_sampling = settings.get("initial_sampling", "site_bound_harmonic")
    if initial_sampling not in ("site_bound_harmonic", "legacy_nominal"):
        raise ValueError(f"Unknown initial sampling policy {initial_sampling}")
    if initial_sampling == "site_bound_harmonic":
        pool = sample_site_bound_harmonic(seed, temperature, physics, draws.slm_depth_factor,
                                          dimensions=1 if arm == "baseline_1d" else 3)
        pos,vel = pool["pos_m"],pool["vel_m_per_s"]
        sampling = pool["sampler"]
    elif arm == "baseline_1d":
        pool = sample_initial_states(cfg.base, seed, shots, well="slm")
        pos,vel = np.zeros((shots,3)),np.zeros((shots,3))
        pos[:,0],vel[:,0] = pool["x0_m"],pool["v0_m_per_s"]
        sampling = "original Level 2C 1D harmonic rejection"
    else:
        m,d,w = physics["mass_kg"],physics["slm_depth_j"],physics["slm_waist_m"]
        zr = np.pi*w*w/1061e-9
        pool = sample_initial_states_3d(seed,shots,temperature,m,d,w,zr,
                                        np.sqrt(4*d/(m*w*w)),np.sqrt(2*d/(m*zr*zr)))
        pos,vel = pool["pos_m"],pool["vel_m_per_s"]
        sampling = pool["sampler"]
    heating = heating_from_config(cfg.base)
    powers = heating.channel_powers(0,0,0,0) if heating else {"total":0.0}
    return cfg,physics,pos,vel,draws,powers,sampling


def run_case(settings, *, stage, temperature, key, arm, shots, seed, dt_us, out):
    cfg,physics,pos,vel,draws,powers,sampling = case_inputs(settings,temperature,shots,seed,arm)
    wave = waveform_for(key,cfg,ROOT/settings["ml_waveform"])
    protocol = build_protocol(wave,cfg.roundtrip.wait_s,dt_us*1e-6,
        arm="extended" if arm == "extended" else "matched",
        long_distance_m=settings["long_distance_um"]*1e-6,
        long_duration_s=settings["long_duration_us"]*1e-6,
        remote_hold_s=settings["remote_hold_us"]*1e-6,
        omega=2*np.pi*settings["rabi_frequency_hz"])
    rounds = settings["max_one_way_transfers"]//2
    vs = settings["vs_m_s"] if arm in ("matched","extended") else 0.0
    t0 = time.perf_counter()
    result = run_continuous(pos,vel,protocol,rounds=rounds,physics=physics,draws=draws,
        jitter_sigma_m=cfg.disorder.shot_jitter_alignment_sigma_um*1e-6 if cfg.disorder.enabled else 0,
        power_w=powers["total"],vs_m_s=vs,eta=settings["eta"])
    runtime = time.perf_counter()-t0
    success = result["alive"].sum(axis=0)
    survival = success/shots
    refs = load_reference_curves(cfg.reference_csv)
    n = np.arange(2*rounds+1)
    nr,sr = refs[key]
    keep = nr <= n[-1]
    comparison = compare_to_reference(n,survival,(nr[keep],sr[keep]),cfg.acceptance_tolerance_band)
    comparison["rmse"] = float(np.sqrt(np.mean(np.asarray(comparison["deviation"])**2)))
    # The previous report smoothed only fixed anchors; preserve as secondary view.
    anchors = [5,8,12,16,20,30,31] if key == "manual_200us" else [5,8,12,16,20,30,40,50,60]
    anchors = np.array([v for v in anchors if nr.min() <= v <= min(nr.max(),n[-1])])
    anchor_comparison = None
    if len(anchors):
        reference = np.minimum.accumulate(np.interp(anchors,nr,sr))
        anchor_comparison = compare_to_reference(n,survival,(anchors,reference),cfg.acceptance_tolerance_band)
    low,high = wilson(success,shots)
    case_id = f"{stage}_T{temperature:g}_{arm}_{key}_N{shots}_dt{dt_us:g}"
    if settings.get("initial_sampling", "site_bound_harmonic") == "site_bound_harmonic":
        case_id += "_initSite"
    if settings.get("parametric_reference_temperature_uK") is not None:
        case_id += f"_Tref{settings['parametric_reference_temperature_uK']:g}"
    source_paths = [ROOT/settings["level2c_config"],ROOT/settings["ml_waveform"],cfg.reference_csv]
    source_paths += sorted(Path(__file__).parent.glob("*.py"))
    metadata = {"case_id":case_id,"stage":stage,"temperature_uK":temperature,"waveform":key,
        "arm":arm,"shots":shots,"seed":seed,"dt_us":dt_us,"runtime_s":runtime,
        "rounds":rounds,"one_way_transfers":2*rounds,"sampling":sampling,
        "duration_per_round_us":protocol.duration_s*1e6,
        "distance_per_round_um":protocol.distance_per_round_m*1e6,
        "dd_mode":protocol.dd_mode,"pulses_per_round":len(protocol.pulses),
        "vs_m_s":vs,"vs_provenance":"relative Level 3 severity scenario, uncalibrated",
        "eta":settings["eta"],"rabi_frequency_hz":settings["rabi_frequency_hz"],
        "heating_power_uK_per_s":{k:v/1.380649e-29 for k,v in powers.items()},
        "parametric_Tref_uK":(temperature if cfg.base.noise_mean_heating.parametric_reference_temperature_uK is None
                               else cfg.base.noise_mean_heating.parametric_reference_temperature_uK),
        "initial_unbound_actual_site":int(np.sum(slm_energies(pos,vel,physics,draws.slm_depth_factor) >= 0)),
        "heating_policy":"Level 2C constant total power; vector velocity rescale, no multiplication by dimension",
        "jitter_policy":"independent round/leg/original-shot streams; shared across waveforms and protocol arms",
        "sources":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "settings":settings,"n":n,"survival_counts":success,"survival":survival,
        "wilson95_low":low,"wilson95_high":high,"comparison_raw":comparison,
        "comparison_v8_anchors":anchor_comparison,
        "paper_protocol_comparable":arm != "extended",
        "losses_by_stage":{name:int(np.sum(result["lost_stage"]==i)) for i,name in enumerate(protocol.names)},
        "quantum":quantum_summary(result,protocol,rounds),
        "numeric_backend":"numba" if ACCELERATED else "python",
        "runtime_versions":{"python":platform.python_version(),"numpy":np.__version__}}
    np.savez_compressed(out/f"{case_id}.npz",**result,initial_pos=pos,initial_vel=vel,
                        initial_slm_depth_factor=draws.slm_depth_factor,
                        initial_energy_j=slm_energies(pos,vel,physics,draws.slm_depth_factor))
    write_json(out/f"{case_id}.json",metadata)
    print(f"{case_id}: S({comparison['n'][-1]})={comparison['simulated'][-1]:.4f} "
          f"maxdev={comparison['max_abs_deviation']:.4f} pass={comparison['matches_reference']} "
          f"({runtime:.1f}s)",flush=True)
    return metadata


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path,default=ROOT/"configs/continuous_level2c_345.yaml")
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--stage",choices=("preflight","scan","validation","convergence","all"),default="scan")
    p.add_argument("--arms",nargs="+",choices=ARMS)
    p.add_argument("--temperatures",nargs="+",type=float)
    p.add_argument("--waveforms",nargs="+")
    p.add_argument("--shots",type=int)
    p.add_argument("--transfers",type=int)
    p.add_argument("--threads",type=int)
    args = p.parse_args(argv)
    settings = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if args.transfers:
        settings["max_one_way_transfers"] = args.transfers
    if settings["max_one_way_transfers"] < 2 or settings["max_one_way_transfers"] % 2:
        p.error("Continuous roundtrip count requires a positive even transfer count")
    if ACCELERATED:
        from numba import set_num_threads
        set_num_threads(args.threads or settings["threads"])
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.stage == "preflight":
        import subprocess
        return subprocess.call([sys.executable,"-m","pytest",str(ROOT/"tests/test_continuous_transfer.py"),"-q"])
    stages = ("scan","validation","convergence") if args.stage == "all" else (args.stage,)
    for stage in stages:
        shots = args.shots or settings[f"{stage}_shots"]
        seed = settings["scan_seed" if stage == "scan" else "validation_seed"]
        dt_us = settings["convergence_dt_us" if stage == "convergence" else "dt_us"]
        for temperature in args.temperatures or settings["temperatures_uK"]:
            for arm in args.arms or settings["arms"]:
                for key in args.waveforms or settings["waveforms"]:
                    run_case(settings,stage=stage,temperature=temperature,key=key,arm=arm,
                             shots=shots,seed=seed,dt_us=dt_us,out=args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
