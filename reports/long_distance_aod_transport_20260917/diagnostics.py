from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from continuous_transfer.engine import field
from long_distance_transport.model import Config,KB,metrics
ROOT=Path(__file__).resolve().parent

def main():
    c=Config();D=c.depth_uK*KB*1e-6;w=c.waist;zr=np.pi*w*w/c.wavelength
    ac=metrics(c)['maximum_static_force_acceleration']
    rows=[]
    for ratio in np.linspace(.01,.999,100):
        a=ratio*ac
        f=lambda x:float(field(x,0.,0.,D,w,zr,0.,0.)[1]-c.mass*a)
        xm=brentq(f,-w/2,0,xtol=1e-18);xs=brentq(f,-8*w,-w/2,xtol=1e-18)
        ue=lambda x:-D*np.exp(-2*x*x/w**2)+c.mass*a*x
        rows.append(dict(a_over_acrit=ratio,acceleration=a,minimum_x_um=xm*1e6,saddle_x_um=xs*1e6,barrier_uK=(ue(xs)-ue(xm))/KB*1e6))
    pd.DataFrame(rows).to_csv(ROOT/'inertial_barrier.csv',index=False)
    lens=[]
    for geom in ['x','diagonal','opposite']:
        for v in np.linspace(0,1.2,121):
            vx,vy=(v,0) if geom=='x' else (v/np.sqrt(2),v/np.sqrt(2))
            shifts=np.array([vx,vy])*5e-6/w*zr
            if geom=='opposite':shifts[1]*=-1
            zs=np.linspace(min(0,min(shifts))-5*zr,max(0,max(shifts))+5*zr,10001)
            a=1+((zs-shifts[0])/zr)**2;b=1+((zs-shifts[1])/zr)**2
            u=-D/np.sqrt(a*b);ix=int(np.argmin(u));dz=zs[1]-zs[0]
            curvature=(u[ix+1]-2*u[ix]+u[ix-1])/dz**2
            lens.append(dict(geometry=geom,path_speed_m_s=v,focus1_um=shifts[0]*1e6,focus2_um=shifts[1]*1e6,
                axial_minimum_um=zs[ix]*1e6,actual_depth_uK=-u[ix]/KB*1e6,
                axial_frequency_Hz=np.sqrt(max(curvature,0)/c.mass)/(2*np.pi)))
    tab=pd.DataFrame(lens);tab.to_csv(ROOT/'lensing_geometry.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    b=pd.DataFrame(rows);axs[0].plot(b.a_over_acrit,b.barrier_uK);axs[0].set(xlabel='a / acrit',ylabel='Tilted barrier (uK)',title=f'acrit={ac:.0f} m/s²')
    for geom,g in tab.groupby('geometry'):
        axs[1].plot(g.path_speed_m_s,g.actual_depth_uK,label=geom)
        axs[2].plot(g.path_speed_m_s,g.axial_frequency_Hz/1000,label=geom)
    axs[1].set(xlabel='Path speed (m/s)',ylabel='Actual depth (uK)',title='tau=5 us');axs[2].set(xlabel='Path speed (m/s)',ylabel='Axial frequency (kHz)')
    for ax in axs:ax.grid(alpha=.2)
    axs[2].legend();fig.tight_layout();fig.savefig(ROOT/'force_lens_diagnostics.png');plt.close(fig)
    # Stage-by-stage diagnostics from fixed controls and fresh-seed runs.
    ledger=[];residual=[]
    for path in (ROOT/'runs/confirm').glob('*.json'):
        r=json.loads(path.read_text());cfg=r['config']
        if cfg['repeats']!=50:continue
        raw=np.load(path.with_suffix('.npz'));st=raw['stage_energy_j']/KB*1e6;alive=raw['alive'];e=raw['energy_j']/KB*1e6
        fac=raw['site_factors'][:,1 if cfg['protocol']=='B' else 0]
        depth=cfg['slm_depth_uK'] if cfg['protocol']=='B' else cfg['depth_uK']
        for n in [1,10,50]:
            mask=alive[:,n];energy=e[mask,n]+depth*fac[mask]
            if len(energy):
                delta=energy-(e[mask,0]+depth*fac[mask])
                q=raw['states'][mask,n,:3].copy();v=raw['states'][mask,n,3:]
                q[:,0]-=cfg['distance_um']*1e-6 if n%2 else 0.
                if cfg['protocol']=='A':q[:,0]-=raw['offsets'][mask]
                moments={f'mean_{axis}_um':float(q[:,i].mean()*1e6) for i,axis in enumerate('xyz')}
                moments.update({f'mean_v{axis}_m_s':float(v[:,i].mean()) for i,axis in enumerate('xyz')})
                moments.update({f'rms_{axis}_um':float(np.sqrt(np.mean(q[:,i]**2))*1e6) for i,axis in enumerate('xyz')})
                residual.append(dict(protocol=cfg['protocol'],duration_us=cfg['duration_us'],rin=cfg['rin'],seed=cfg['seed'],repeat=n,
                    retained=len(energy),mean_above_actual_bottom_uK=float(energy.mean()),p95_above_actual_bottom_uK=float(np.quantile(energy,.95)),
                    paired_survivor_mean_energy_change_uK=float(delta.mean()),**moments))
        lost=alive[:,:-1]&~alive[:,1:]
        for j,label in enumerate(r['stage_labels']):
            finite=np.isfinite(st[:,:,j])&alive[:,:-1]
            ledger.append(dict(file=path.name,protocol=cfg['protocol'],duration_us=cfg['duration_us'],rin=cfg['rin'],seed=cfg['seed'],stage=label,
                mean_relative_energy_uK=float(np.mean(st[:,:,j][finite])),
                positive_energy_before_later_failure=int(np.sum(lost&(st[:,:,j]>=0))),
                total_endpoint_failures=int(lost.sum()),
                note='first positive diagnostic is not a proven irreversible escape time'))
    pd.DataFrame(ledger).to_csv(ROOT/'loss_timing_diagnostics.csv',index=False)
    pd.DataFrame(residual).to_csv(ROOT/'residual_excitation.csv',index=False)

if __name__=='__main__':main()
