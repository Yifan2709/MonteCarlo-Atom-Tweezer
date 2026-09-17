"""Analysis from raw results only; no survival curve serves as a dynamics input."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from long_distance_transport.model import interval
ROOT=Path(__file__).resolve().parent
TR=['rb_cubic','zero_jerk','min_jerk']
FORMAL={'confirm','operations_confirm','repeat_confirm','final_single'}
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'figure.dpi':130})

def save(fig,name):
    fig.tight_layout();fig.savefig(ROOT/(name+'.png'));fig.savefig(ROOT/(name+'.svg'));plt.close(fig)

def tables():
    rows=[]
    for path in sorted((ROOT/'runs').glob('*/*.json')):
        r=json.loads(path.read_text());c=r['config']
        previous=c['shots']
        for o in r['observations']:
            rows.append({**c,**r['kinematics'],**o,'group':path.parent.name,'file':str(path.relative_to(ROOT)),
                'previous_alive':previous,'conditional_loss':1-o['alive']/previous if previous else None,
                'elapsed_total_us':r['kinematics']['complete_oneway_us']*o['repeat']})
            previous=o['alive']
    df=pd.DataFrame(rows);df.to_csv(ROOT/'all_points.csv',index=False)
    formal=df[df.group.isin(FORMAL)].copy()
    keys=list(json.loads(next((ROOT/'runs/coarse').glob('*.json')).read_text())['config'])
    keys=[k for k in keys if k not in ['shots','seed','repeats']]+['repeat']
    # Fresh validation seeds at the same production timestep are additional
    # evidence, including unfavorable draws. Do not discard them from a final
    # condition already selected for independent confirmation. Other timesteps
    # remain numerical controls, never extra independent atom trials.
    extra=df[df.group.isin(['convergence','boundary_steps','repeat_steps'])&(df.dt_us==.05)]
    extra=extra.merge(formal[keys].drop_duplicates(),on=keys,how='inner')
    formal=pd.concat([formal,extra],ignore_index=True)
    formal=formal.sort_values('shots').drop_duplicates(keys+['seed'],keep='last')
    pooled=[];energy_cache={}
    for values,g in formal.groupby(keys,dropna=False):
        row=g.iloc[0].to_dict();n=int(g.shots.sum());k=int(g.alive.sum());lo,hi=interval(k,n)
        row['representative_file']=row.pop('file');row['representative_seed']=row.pop('seed')
        distributions=[]
        for source in g.itertuples():
            if source.file not in energy_cache:
                with np.load((ROOT/source.file).with_suffix('.npz')) as raw:
                    energy_cache[source.file]=(raw['energy_j'],raw['alive'])
            en,retained=energy_cache[source.file];r=int(source.repeat)
            distributions.append(en[retained[:,r],r]/1.380649e-29)
        energy=np.concatenate(distributions);previous=int(g.previous_alive.sum())
        row.update(shots=n,alive=k,survival=k/n,lower=lo,upper=hi,seeds=','.join(map(str,g.seed)),seed_count=len(g),
            source_files=';'.join(g.file),
            previous_alive=previous,conditional_loss=1-k/previous if previous else None,
            median_energy_uK=float(np.median(energy)) if len(energy) else None,
            energy_p95_uK=float(np.quantile(energy,.95)) if len(energy) else None,
            status99='supported' if lo>=.99 else 'failed' if hi<.99 else 'uncertain',
            status999='supported' if lo>=.999 else 'failed' if hi<.999 else 'uncertain')
        pooled.append(row)
    pool=pd.DataFrame(pooled);pool.to_csv(ROOT/'independent_confirmation.csv',index=False)
    boundaries=[]
    for key,g in pool[pool.repeat.isin([1,10,50])].groupby(['distance_um','trajectory','protocol','tau_us','rin','disorder','repeat']):
        for target,status in [(.99,'status99'),(.999,'status999')]:
            good=g[g[status]=='supported'].sort_values('duration_us')
            boundaries.append(dict(zip(['distance_um','trajectory','protocol','tau_us','rin','disorder','repeat'],key))|
                dict(target=target,supported_times_us=';'.join(str(x) for x in good.duration_us.unique()),
                    fastest_confirmed_time_us=float(good.duration_us.min()) if len(good) else None,
                    maximum_confirmed_average_m_s=float(good.average_m_s.max()) if len(good) else None,
                    maximum_confirmed_peak_m_s=float(good.peak_m_s.max()) if len(good) else None,
                    minimum_complete_us=float(good.complete_oneway_us.min()) if len(good) else None,
                    note='statistical support in discrete points only; read timestep and parameter uncertainty separately; intervals NOT certified; no global optimum'))
    pd.DataFrame(boundaries).to_csv(ROOT/'boundaries.csv',index=False)
    return df,pool

def convergence(df):
    sub=df[(df.group.isin(['convergence','boundary_steps','repeat_steps']))&df.repeat.isin([1,10,50])]
    keys=['distance_um','duration_us','trajectory','protocol','tau_us','rin','repeat','seed','noise_us']
    rows=[]
    for key,g in sub.groupby(keys):
        g=g.sort_values('dt_us',ascending=False)
        for (_,a),(_,b) in zip(g.iloc[:-1].iterrows(),g.iloc[1:].iterrows()):
            if a.dt_us==b.dt_us:continue
            ra=np.load((ROOT/a.file).with_suffix('.npz'));rb=np.load((ROOT/b.file).with_suffix('.npz'));n=int(a['repeat'])
            aa,bb=ra['alive'][:,n],rb['alive'][:,n]
            ea,eb=ra['energy_j'][aa,n],rb['energy_j'][bb,n]
            diff=aa.astype(int)-bb.astype(int)
            rows.append(dict(zip(keys,key))|dict(dt_coarse=a.dt_us,dt_fine=b.dt_us,
                survival_difference=float(diff.mean()),paired_SE=float(diff.std(ddof=1)/np.sqrt(len(diff))),
                label_disagreement=float(np.mean(aa!=bb)),energy_KS=float(ks_2samp(ea,eb).statistic) if len(ea) and len(eb) else None,
                median_energy_difference_uK=(float(np.median(ea)-np.median(eb))/1.380649e-29) if len(ea) and len(eb) else None))
    pd.DataFrame(rows).to_csv(ROOT/'mechanical_convergence.csv',index=False)
    sub[sub.seed==26091732].to_csv(ROOT/'noise_step_convergence.csv',index=False)

def plots(df,pool):
    base=df[(df.group=='coarse')&(df.depth_uK==280)&(df.distance_um==510)]
    fig,axs=plt.subplots(1,2,figsize=(12,4))
    for tr in TR:
        for tau in [0.,5.]:
            g=base[(base.trajectory==tr)&(base.tau_us==tau)].sort_values('duration_us')
            for ax,col in zip(axs,['duration_us','average_m_s']):ax.plot(g[col],g.survival,'o-',label=f'{tr}, tau={tau:g} us')
    axs[0].set(xlabel='Move duration (us)',ylabel='Mechanical survival',xscale='log');axs[1].set(xlabel='Average focus speed (m/s)',xscale='log')
    axs[1].legend(fontsize=7);save(fig,'survival_time_speed')
    fig,axs=plt.subplots(2,3,figsize=(13,7))
    for i,tau in enumerate([0.,5.]):
        for j,tr in enumerate(TR):
            g=df[(df.group=='coarse')&(df.depth_uK==280)&(df.trajectory==tr)&(df.tau_us==tau)]
            sc=axs[i,j].scatter(g.duration_us,g.distance_um,c=g.survival,vmin=0,vmax=1,cmap='viridis',s=90,marker='s')
            axs[i,j].set(xscale='log',xlabel='Move duration (us)',ylabel='Distance (um)',title=f'{tr}; tau={tau:g} us')
    fig.subplots_adjust(right=.87,wspace=.35,hspace=.4)
    cax=fig.add_axes([.91,.16,.014,.67])
    fig.colorbar(sc,cax=cax,label='Screening survival')
    fig.savefig(ROOT/'distance_time_region.png');fig.savefig(ROOT/'distance_time_region.svg');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    b=pool[(pool.protocol=='A')&(pool['repeat']==1)&(pool.disorder==0)&(pool.status99=='supported')]
    for ax,tau in zip(axs,[0.,2.,5.]):
        for tr in TR:
            g=b[(b.tau_us==tau)&(b.trajectory==tr)].groupby('distance_um').average_m_s.max()
            ax.plot(g.index,g.values,'o-',label=tr)
        ax.set(xlabel='Distance (um)',ylabel='Max confirmed average speed (m/s)',title=f'99% mechanical, tau={tau:g} us')
    axs[-1].legend(fontsize=8);save(fig,'boundary_comparison')
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for ax,tr in zip(axs,TR):
        for proto,rin,dis,label in [('A',0.,0.,'A ideal'),('A',5.75e-9,0.,'A noise'),('B',0.,0.,'B ideal'),('B',0.,.02,'B disorder'),('B',5.75e-9,.02,'B disorder+noise')]:
            g=df[(df.group.isin(['coarse','mechanisms']))&(df.trajectory==tr)&(df.protocol==proto)&(df.rin==rin)&(df.disorder==dis)&(df.tau_us==5)&(df.depth_uK==280)&(df.distance_um==510)&(df.geometry=='x')].sort_values('duration_us')
            ax.plot(g.duration_us,g.survival,'o-',label=label)
        ax.set(xlabel='Move duration (us)',ylabel='Survival',title=tr,xscale='log')
    axs[-1].legend(fontsize=7);save(fig,'mechanisms')
    fig,axs=plt.subplots(2,3,figsize=(13,7))
    for i,proto in enumerate(['A','B']):
        for j,tr in enumerate(TR):
            for n in [1,10,50]:
                for noise in [False,True]:
                    g=df[(df.group=='repeat')&(df.protocol==proto)&(df.trajectory==tr)&(df['repeat']==n)&((df.rin>0)==noise)].sort_values('duration_us')
                    axs[i,j].plot(g.duration_us,g.survival,'o--' if noise else 'o-',label=f'{n} moves'+(' +noise' if noise else ''))
            axs[i,j].set(title=f'{proto}; {tr}',xscale='log',xlabel='Move duration (us)',ylabel='Cumulative survival')
    axs[-1,-1].legend(fontsize=7);save(fig,'repeated_transport')
    fig,axs=plt.subplots(1,3,figsize=(14,4))
    for noise in [False,True]:
        g=df[(df.group=='repeat')&(df.protocol=='B')&(df.trajectory=='zero_jerk')&(df.duration_us==1600)&((df.rin>0)==noise)].sort_values('repeat')
        if len(g):
            label='noise' if noise else 'no noise';raw=np.load((ROOT/g.iloc[0].file).with_suffix('.npz'))
            excess=raw['energy_j']/1.380649e-29+140*raw['site_factors'][:,1,None]
            for n in [0,1,10,50]:
                mask=raw['alive'][:,n]
                if mask.any():axs[0].hist(excess[mask,n],bins=35,histtype='step',density=True,label=f'{label}, n={n}')
            axs[1].plot(g['repeat'],g.conditional_loss,label=label)
            med=[np.median(excess[raw['alive'][:,n],n]) if raw['alive'][:,n].any() else np.nan for n in range(1,51)]
            axs[2].plot(range(1,51),med,label=label)
    axs[0].set(xlabel='Survivor energy above final trap bottom (uK)',ylabel='Density');axs[1].set(xlabel='Move number',ylabel='Conditional loss');axs[2].set(xlabel='Move number',ylabel='Survivor median energy (uK)')
    for ax in axs:ax.legend(fontsize=7)
    save(fig,'energy_conditional_loss')
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    sens=df[(df.group=='sensitivity')&(df.trajectory=='zero_jerk')&(df['repeat']==10)]
    for ax,col,default in zip(axs,['temperature_uK','waist','tau_us'],[19.,1.11e-6,5.]):
        g=sens.copy()
        for other,val in [('temperature_uK',19.),('waist',1.11e-6),('tau_us',5.),('ramp_us',200.)]:
            if other!=col:g=g[g[other]==val]
        for val,h in g.groupby(col):
            h=h.sort_values('duration_us');ax.plot(h.duration_us,h.survival,'o-',label=f'{val*1e6:g} um' if col=='waist' else f'{val:g}')
        ax.set(xlabel='Move duration (us)',ylabel='10-move B survival',title=col);ax.legend(fontsize=8)
    save(fig,'sensitivity')
    fig,axs=plt.subplots(1,2,figsize=(12,4))
    for dt,g in df[(df.group=='convergence')&(df.seed==26091731)&(df['repeat']==1)&(df.protocol=='A')&(df.disorder==0)].groupby('dt_us'):
        axs[0].scatter(g.duration_us,g.survival,label=f'dt={dt:g} us',alpha=.7)
    g=pool[(pool.distance_um==510)&(pool['repeat']==1)&(pool.protocol=='A')&(pool.disorder==0)]
    for tr,h in g.groupby('trajectory'):
        axs[1].errorbar(h.average_m_s,h.survival,yerr=[h.survival-h.lower,h.upper-h.survival],fmt='o',label=tr)
    axs[0].set(xlabel='Move duration (us)',ylabel='Survival',title='Timestep controls');axs[1].set(xlabel='Average speed (m/s)',ylabel='Survival',title='Independent seeds, exact 95% binomial CI',ylim=(.98,1.001))
    for ax in axs:ax.legend(fontsize=7)
    save(fig,'sampling_timestep')

def paper_plots(df):
    yb=pd.read_csv(ROOT/'zhang_source/dataset/extended_data_fig2b.csv')
    fig,axs=plt.subplots(1,2,figsize=(12,4));res=[]
    for ax,tr in zip(axs,TR[1:]):
        tag=tr if tr=='zero_jerk' else 'min_jerk'
        t=yb[f'exp_moving_time_us_{tag}'];p=yb[f'exp_survival_{tag}'];e=yb[f'exp_survival_err_{tag}']
        ax.errorbar(t,p,e,fmt='o',color='black',ms=3,label='Author experimental source')
        for on in ['w','wo']:
            ax.plot(yb[f'sim_{on}_lensing_moving_time_us_{tag}'],yb[f'sim_{on}_lensing_survival_{tag}'],label=f'Author simulation {on} lens')
        for tau in [0.,2.,5.]:
            g=df[(df.group=='native')&(df.species=='Yb')&(df.trajectory==tr)&(df.tau_us==tau)].sort_values('duration_us')
            ax.plot(g.duration_us,g.survival,'--',label=f'Our conditional tau={tau:g} us')
            exp=pd.DataFrame({'duration_us':t,'exp':p,'exp_error':e}).dropna();merged=g.merge(exp,on='duration_us')
            merged['residual']=merged.survival-merged.exp
            merged.to_csv(ROOT/f'yb_residual_{tr}_tau{tau:g}.csv',index=False)
            res.append(dict(species='Yb',trajectory=tr,tau_us=tau,points=len(merged),RMSE=float(np.sqrt(np.mean(merged.residual**2))),max_error=float(abs(merged.residual).max()),status='conditional comparison; apparatus input incomplete'))
        ax.set(xlabel='Move duration (us)',ylabel='Survival',title=tr);ax.legend(fontsize=6)
    save(fig,'yb_paper_comparison')
    rb=pd.read_csv(ROOT/'rb_ed2a_digitized.csv')
    g=df[(df.group=='native')&(df.species=='Rb')].sort_values('duration_us')
    fig,axs=plt.subplots(1,2,figsize=(12,4))
    axs[0].plot(rb.separation_speed_m_s,rb.survival,'o',label='Paper ED2a vector markers')
    axs[0].plot(2*g.average_m_s,g.survival,'o-',label='Our conditional Gaussian MC (10 uK)')
    # Interpolation is analysis of comparison points, never a loss model input.
    pred=np.interp(rb.duration_us,g.duration_us,g.survival,left=np.nan,right=np.nan)
    rb['our_survival_interpolated_for_comparison']=pred;rb['residual']=pred-rb.survival
    rb.to_csv(ROOT/'rb_comparison_residuals.csv',index=False)
    res.append(dict(species='Rb',trajectory='rb_cubic',tau_us=0,points=int(np.isfinite(pred).sum()),RMSE=float(np.sqrt(np.nanmean(rb.residual**2))),max_error=float(np.nanmax(abs(rb.residual))),status='initial temperature and lens not independently calibrated; plot error bars unavailable'))
    h=pd.read_csv(ROOT/'rb_harmonic_averaging.csv');axs[1].loglog(h.duration_us,h.exact_mean_N,label='Exact Fourier + frequency average');axs[1].loglog(h.duration_us,h.asym_mean_N,'--',label='Paper large-omegaT Eq2')
    axs[0].set(xlabel='Pair separation speed 2D/T (m/s)',ylabel='Atom retention');axs[1].set(xlabel='Move duration (us)',ylabel='Added oscillator quanta')
    for ax in axs:ax.legend(fontsize=8)
    save(fig,'rb_paper_comparison');pd.DataFrame(res).to_csv(ROOT/'paper_comparison_metrics.csv',index=False)

if __name__=='__main__':
    df,pool=tables();convergence(df);plots(df,pool);paper_plots(df)
    print('analyzed',len(df),'observations;',len(pool),'pooled independent conditions')
