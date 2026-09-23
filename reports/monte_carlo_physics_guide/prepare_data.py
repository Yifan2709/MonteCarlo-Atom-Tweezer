"""Prepare editable chart values. No images are drawn and no full run is repeated."""
from pathlib import Path
import csv, json, hashlib, importlib.util
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = HERE / '.build'
BUILD.mkdir(exist_ok=True)
CAL = ROOT / 'reports/movement_effects_v1_calibration_20260916'
def read_csv(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def clean(v):
    if isinstance(v, dict): return {k: clean(x) for k,x in v.items()}
    if isinstance(v, (list,tuple)): return [clean(x) for x in v]
    if isinstance(v, np.ndarray): return clean(v.tolist())
    if isinstance(v, np.generic): return v.item()
    return v

data = {}
for key, fn in [('curves','survival_curves.csv'),('points','point_residuals.csv'),('metrics','metrics.csv')]:
    data[key] = read_csv(CAL/fn)

# Formula-based teaching charts use nominal 1.17 um, 140/280 uK, T=15 uK.
kb, mass, w = 1.380649e-23, 132.90545196*1.66053906660e-27, 1.17e-6
zr = np.pi*w*w/1061e-9
x = np.linspace(-2.5,2.5,101)
data['well'] = {'x':x, 'slm':-140*np.exp(-2*x*x/1.17**2),
                'aod':-280*np.exp(-2*x*x/1.17**2),
                'harmonic':-140+2*140*x*x/1.17**2,
                'force':-4*140*kb*1e-6*(x*1e-6)/(w*w)*np.exp(-2*x*x/1.17**2)/1e-21}
z = np.linspace(-8,8,129)
data['axial'] = {'x':z, 'radial':-140*np.exp(-2*z*z/1.17**2),
                 'axial':-140/(1+(z*1e-6/zr)**2)}
data['widths'] = []
for T,D,ww,wl in [(5,280,1.17,1055),(15,140,1.17,1061),(19,140,1.11,1061)]:
    zrr=np.pi*(ww*1e-6)**2/(wl*1e-9)
    data['widths'].append({'T':T,'D':D,'w':ww,'sr':np.sqrt(T*ww**2/(4*D)),
      'sz':np.sqrt(T*(zrr*1e6)**2/(2*D)), 'sv':np.sqrt(kb*T*1e-6/mass)*100,
      'fr':np.sqrt(4*D*kb*1e-6/(mass*(ww*1e-6)**2))/(2*np.pi)*1e-3,
      'fz':np.sqrt(2*D*kb*1e-6/(mass*zrr*zrr))/(2*np.pi)*1e-3})

# Directly use the repository's actual-site, bound initial-state preparation.
src = ROOT/'simulation/level2_joint_transfer/src/continuous_transfer/initialization.py'
spec=importlib.util.spec_from_file_location('initialization',src)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
rng=np.random.default_rng(20260922)
factors=1+rng.normal(0,.114,4096)
phys={'slm_depth_j':140*kb*1e-6,'slm_waist_m':w,'mass_kg':mass,'slm_center_m':0.}
pool=mod.sample_site_bound_harmonic(20260923,15,phys,factors)
data['sample']={'x':pool['pos_m'][:180,0]*1e6,'v':pool['vel_m_per_s'][:180,0]*100,
                'minE':pool['initial_energy_j'].min()/kb*1e6,
                'maxE':pool['initial_energy_j'].max()/kb*1e6,
                'retries':int(np.sum(pool['attempts_per_site']>1))}
counts,edges=np.histogram(140*factors,bins=np.arange(80,205,5))
data['depthHist']={'x':(edges[:-1]+edges[1:])/2,'y':counts}
gx=np.linspace(-.8,.8,81)
data['gaussian']={'x':gx,'y':np.exp(-.5*(gx/data['widths'][1]['sr'])**2)}

# Manual controls and the archived, PDF-derived ML controls. Not-a-knot spline.
def spline(x,y,q):
    n=len(x); h=np.diff(x); A=np.zeros((n,n)); b=np.zeros(n)
    A[0,0],A[0,1],A[0,2]=-h[1],h[0]+h[1],-h[0]
    A[-1,-3],A[-1,-2],A[-1,-1]=-h[-1],h[-1]+h[-2],-h[-2]
    for i in range(1,n-1):
        A[i,i-1:i+2]=[h[i-1],2*(h[i-1]+h[i]),h[i]]
        b[i]=6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1])
    M=np.linalg.solve(A,b)
    k=np.clip(np.searchsorted(x,q)-1,0,n-2); hh=x[k+1]-x[k]
    a=(x[k+1]-q)/hh; bb=(q-x[k])/hh
    return a*y[k]+bb*y[k+1]+((a**3-a)*M[k]+(bb**3-bb)*M[k+1])*hh**2/6
t=np.linspace(0,400,161); q=np.clip((t-192)/208,0,1)
nodes=json.loads((CAL/'ml_vector_nodes.json').read_text())
data['wave']={'t':t,'manual_pos':2.4*(3*q*q-2*q*q*q),
              'manual_depth':280*np.minimum(t/192,1)**2,
              'ml_pos':spline(np.array(nodes['grid_us']),np.array(nodes['pos_um']),t),
              'ml_depth':spline(np.array(nodes['grid_us']),1000*np.array(nodes['depth_mK']),t)}
xx=np.linspace(-1.7,4.1,117)
data['snapshots']={'x':xx,'curves':[]}
for d,c in [(0,0),(280,0),(280,1.2),(280,2.4)]:
    data['snapshots']['curves'].append(-140*np.exp(-2*xx*xx/1.17**2)-d*np.exp(-2*(xx-c)**2/1.17**2))
trajpath=ROOT/'simulation/level1_transfer_1d/outputs/level1_transfer_1d/demo/thermal_example_trajectory.csv'
traj=read_csv(trajpath)[::40]
data['trajectory']={k:[float(r[k]) for r in traj] for k in ['time_us','x_um','aod_center_um','kinetic_uK','total_energy_uK','total_potential_uK']}

# Quantitative illustrative contrasts, not fitted or observed distributions.
data['toyEnergy']={'x':[20,40,60,80,100,120,140,160],
                  'narrow':[0,0,20,60,20,0,0,0],
                  'wide':[20,10,20,10,10,10,10,10]}
# Equal mean 80 uK for both (wide corrected at 80/160 bins).
data['toyEnergy']['wide']=[20,10,20,10,10,10,10,10]
# Wide mean = 80 uK, narrow mean = 80 uK.
assert sum(a*b for a,b in zip(data['toyEnergy']['x'],data['toyEnergy']['wide'])) == 8000
s=np.linspace(0,3,61)
data['coherence']={'x':s,'y':np.exp(-s*s/2)}
ts=np.linspace(0,1,101)
data['echo']={'x':ts,'free':ts,'echo':np.where(ts<=.5,ts,1-ts)}
zz=np.linspace(-12,12,121); zrA=np.pi*1.17**2/1.055
data['lensing']={'x':zz,'curves':[]}
for delta in [0,zrA,2*zrA]:
    aa=1+((zz-delta)/zrA)**2;bb=1+((zz+delta)/zrA)**2
    data['lensing']['curves'].append(-280/np.sqrt(aa*bb))
sources=[src,trajpath,CAL/'survival_curves.csv',CAL/'point_residuals.csv',CAL/'metrics.csv',CAL/'ml_vector_nodes.json',ROOT/'docs/monte_carlo_plain_guide.md']
data['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
(BUILD/'data.json').write_text(json.dumps(clean(data),ensure_ascii=False),encoding='utf-8')
print(json.dumps({'widths':clean(data['widths']),'initial_resampled_sites':data['sample']['retries'],'groups':sorted(set(r['group'] for r in data['curves']))},ensure_ascii=False))
