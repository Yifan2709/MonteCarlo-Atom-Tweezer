"""用引擎自己的手工波形（constant-jerk）复跑 2×2 交叉，消除剖面混杂。"""
import json, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, "src")
from scipy.interpolate import CubicSpline
from level2c_pickup_survival.level2c_config import load_level2c_config
from level2c_pickup_survival.survival_engine import RoundtripLegs, draw_disorder, run_survival
from level2_joint_transfer.level2_simulation import physics_from_config, build_waveform, sample_initial_states
from level2_joint_transfer.noise_heating import heating_from_config

UK_J = 1.380649e-29
wf = json.load(open("/tmp/scan_stage1/paper_ml_waveform.json"))
grid_us = np.array(wf["grid_us"]); depth_mK = np.array(wf["depth_mK"]); pos_um = np.array(wf["pos_um"])
ml_depth = CubicSpline(grid_us, depth_mK); ml_pos = CubicSpline(grid_us, pos_um)

cfg = load_level2c_config(Path("configs/scan_stage1/scan_T40.yaml"))
physics = physics_from_config(cfg.base)
manual = build_waveform({"type":"pickup_sequential","duration_s":400e-6,
                         "ramp_fraction":cfg.pickup.ramp_fraction,
                         "pickup_direction":"pickup","direction":"pickup"}, physics)
man_depth = lambda t_us: manual.depth(np.asarray(t_us)*1e-6)/UK_J/1000.0   # mK
man_pos   = lambda t_us: manual.center(np.asarray(t_us)*1e-6)*1e6          # µm

def build_legs(depth_fn, center_fn, duration_s, wait_s, dt_s):
    n_pu=int(round(duration_s/dt_s)); n_w=int(round(wait_s/dt_s))
    t1=np.arange(n_pu+1)*dt_s*1e6
    cp=np.asarray(center_fn(t1))*1e-6; dp=np.asarray(depth_fn(t1))*1000.0*UK_J
    cw=np.full(n_w+1,cp[-1]); dw=np.full(n_w+1,dp[-1])
    slm=np.concatenate([np.ones(n_pu+1),np.zeros(n_w+1),np.ones(n_pu+1)])
    return RoundtripLegs(cp,dp,cw,dw,cp[::-1].copy(),dp[::-1].copy(),slm,n_pu,n_w,n_pu)

heating = heating_from_config(cfg.base)
pool = sample_initial_states(cfg.base, cfg.roundtrip.validation_seed, cfg.roundtrip.validation_shots, well="slm")
disorder = draw_disorder(cfg.disorder, pool["shots"], cfg.roundtrip.validation_seed+7)
jitter = cfg.disorder.shot_jitter_alignment_sigma_um*1e-6

def run(df,cf,tag):
    legs=build_legs(df,cf,400e-6,cfg.roundtrip.wait_s,cfg.roundtrip.dt_s)
    out=run_survival(pool["x0_m"],pool["v0_m_per_s"],physics["mass_kg"],legs,
        cfg.roundtrip.dt_s,cfg.roundtrip.max_one_way_transfers,
        slm_depth_j=physics["slm_depth_j"],slm_waist_m=physics["slm_waist_m"],
        aod_waist_m=physics["aod_waist_m"],heating=heating,disorder=disorder,jitter_sigma_m=jitter)
    s60=float(out["survival_fraction"][-1])
    legs_lost={}
    for leg in out["lost_leg"]:
        if leg: legs_lost[leg]=legs_lost.get(leg,0)+1
    print(f"[{tag:8s}] S(60)={s60:.4f}  几何等效每程损失={(1-s60**(1/60))*100:.3f}%  丢失阶段={legs_lost}")
    return s60

print("引擎剖面版 2×2（T40 复核池 512 shots）：")
a=run(man_depth,man_pos,"M,M")
b=run(ml_depth,man_pos,"ML,M")
c=run(man_depth,ml_pos,"M,ML")
d=run(ml_depth,ml_pos,"ML,ML")
print(f"\n位置换成ML的增益: {c-a:+.3f} / {d-b:+.3f}")
print(f"深度换成ML的增益: {b-a:+.3f} / {d-c:+.3f}")
