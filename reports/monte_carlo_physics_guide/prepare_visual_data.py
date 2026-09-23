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

# Five method pages: analytical force and explicitly labelled 1D teaching paths.
# These are numerical chart data, not illustrations rendered to bitmaps.
w = 1.11e-6
depth = 140*kb*1e-6
xx = np.linspace(-1.5, 3., 181)*1e-6
u_static = -140*np.exp(-2*(xx/w)**2)
u_moving = -280*np.exp(-2*((xx-1.2e-6)/w)**2)
fx = np.linspace(-1.5, 1.5, 151)*1e-6
def u1(q, d=depth):
    return -d*np.exp(-2*(q/w)**2)
def force1(q, d=depth):
    return u1(q,d)*4*q/w**2

dt = .0125e-6
def path1(q0, v0, duration_us, *, rin=0., seed=0, d=depth):
    q, v = float(q0), float(v0)
    rg = np.random.default_rng(seed)
    times, qs, vs, energies = [], [], [], []
    for k in range(round(duration_us*1e-6/dt)+1):
        if k % 20 == 0:
            times.append(k*dt*1e6); qs.append(q*1e6); vs.append(v*100)
            energies.append((.5*mass*v*v+u1(q,d))/(kb*1e-6))
        if k == round(duration_us*1e-6/dt): break
        v += .5*dt*force1(q,d)/mass
        q += dt*v
        v += .5*dt*force1(q,d)/mass
        if rin and (k+1) % 4 == 0:
            v += force1(q,d)/mass*np.sqrt(rin*4*dt/2)*rg.normal()
    return {'t':times, 'x':qs, 'v':vs, 'energy':energies}

motion = path1(.4e-6,0,32)
teach_phys = {**phys, 'slm_depth_j':2*depth}
teach_pool = mod.sample_site_bound_harmonic(270930,19,teach_phys,np.ones(3))
ensemble = [path1(teach_pool['pos_m'][i,0],teach_pool['vel_m_per_s'][i,0],100,
    rin=5.75e-9,seed=270931+i,d=2*depth) for i in range(3)]
at = .4e-6
u_at = u1(at)/(kb*1e-6)
loss = {'x':at*1e6, 'u':u_at, 'k':[40.,140.],
        'e':[u_at+40,u_at+140],
        'v':[float(np.sqrt(2*k*kb*1e-6/mass)*100) for k in [40.,140.]]}
D['method'] = {
    'potential':{'x':(xx*1e6).tolist(),'slm':u_static.tolist(),
        'aod':u_moving.tolist(),'total':(u_static+u_moving).tolist()},
    'force':{'x':(fx*1e6).tolist(),'u':(u1(fx)/(kb*1e-6)).tolist(),
        'f':(force1(fx)/1e-21).tolist(), 'point_u':float(u_at),
        'point_f':float(force1(at)/1e-21)},
    'motion':motion,'ensemble':ensemble,'loss':loss}
# Internal consistency checks for the educational quantities.
assert force1(-at)>0 and force1(at)<0 and force1(0)==0
assert loss['e'][0]<0<loss['e'][1]
assert max(motion['energy'])-min(motion['energy']) < .001
assert all(np.isfinite(v) for r in ensemble for v in r['x'])

def clean(v):
    if isinstance(v, dict): return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list): return [clean(x) for x in v]
    return v
(B/'data.json').write_text(json.dumps(clean(D),ensure_ascii=False))
print('Prepared teaching data and retained the original archived comparisons.')
