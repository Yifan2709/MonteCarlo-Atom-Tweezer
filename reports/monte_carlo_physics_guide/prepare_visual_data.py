"""Numerical values for the revised figure-led presentation, not bitmap images."""
from pathlib import Path
import json
import importlib.util
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
B = HERE / '.build' / 'visual-v2'
B.mkdir(parents=True, exist_ok=True)
D = json.loads((HERE / '.build' / 'data.json').read_text())
kb = 1.380649e-23
mass = 132.90545196 * 1.66053906660e-27

# Physical snapshots at the nominal final waist, with no invented atom trajectory.
x = np.linspace(-1.4, 4, 109)
D['transfer'] = {'x': x.tolist(), 'panels': []}
for title, aod_depth, center, slm_depth in [
    ('起点', 0, 0, 140), ('两阱重合', 280, 0, 140),
    ('移动阱离开', 280, 2.4, 140), ('等待', 280, 2.4, 0)]:
    u1 = -slm_depth * np.exp(-2*x*x/1.11**2)
    u2 = -aod_depth * np.exp(-2*(x-center)**2/1.11**2)
    D['transfer']['panels'].append({'title': title, 'slm': u1.tolist(), 'aod': u2.tolist(), 'total': (u1+u2).tolist()})

# Explicit 400 + 100 + 400 us short-transfer control schedule.
t = np.linspace(0,900,451)
local = np.where(t<=400,t,np.where(t<500,400,900-t))
q = np.clip((local-192)/208,0,1)
slm = np.where((t>400)&(t<500),0.,140.)
D['timeline'] = {'t':t.tolist(),'pos':(2.4*(3*q*q-2*q*q*q)).tolist(),
    'aod':(280*np.minimum(local/192,1)**2).tolist(),'slm':slm.tolist()}

# Actual-site initialization at final nominal T and waist, for illustrations only.
src = ROOT / 'simulation/level2_joint_transfer/src/continuous_transfer/initialization.py'
spec = importlib.util.spec_from_file_location('init_for_visual',src)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
factors = 1+np.random.default_rng(2709).normal(0,.114,512)
phys = {'slm_depth_j':140*kb*1e-6,'slm_waist_m':1.11e-6,'mass_kg':mass,'slm_center_m':0.}
pool = mod.sample_site_bound_harmonic(27091,19,phys,factors)
D['initial'] = {'x':(pool['pos_m'][:180,0]*1e6).tolist(),'z':(pool['pos_m'][:180,2]*1e6).tolist(),
    'v':(pool['vel_m_per_s'][:180,0]*100).tolist()}

# Weak-noise harmonic energy illustration: dE = Gamma E dt + sqrt(Gamma) E dW.
# Exact solution, same mean as deterministic drift. Not a production calibration.
rng = np.random.default_rng(27092)
ts = np.linspace(0,.03,121); gamma = 10.; e0 = 40.
W = np.column_stack([np.zeros(8),np.cumsum(rng.normal(0,np.sqrt(ts[1]),(8,120)),axis=1)])
E = e0*np.exp(.5*gamma*ts[None,:]+np.sqrt(gamma)*W)
D['diffusion'] = {'t_ms':(ts*1000).tolist(),'mean':(e0*np.exp(gamma*ts)).tolist(),'paths':E.tolist()}

# Representative checkpoint paths from that same explicit pedagogical model.
check = E[:,::12]-140
D['checkpoints'] = {'n':list(range(11)), 'paths':check.tolist()}

def clean(v):
    if isinstance(v, dict): return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list): return [clean(x) for x in v]
    return v
(B/'data.json').write_text(json.dumps(clean(D),ensure_ascii=False))
print('Prepared teaching data and retained the original archived comparisons.')
