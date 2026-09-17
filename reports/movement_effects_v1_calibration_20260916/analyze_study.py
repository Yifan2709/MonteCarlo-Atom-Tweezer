"""Recount frozen outputs, quantify residuals, and draw publication-readable plots."""
from run_study import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

LABELS=['Manual 200 us','Manual 400 us','Manual 600 us','ML 400 us']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                    'savefig.facecolor':'white'})

def csvwrite(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load_group(ids,wave):
    meta=[json.loads((HERE/'runs'/i/(wave+'.json')).read_text(encoding='utf-8')) for i in ids]
    raw=[np.load(HERE/'runs'/i/(wave+'.npz')) for i in ids]
    for m,z in zip(meta,raw):
        assert np.array_equal(z['alive'].sum(axis=0),m['counts'])
        assert np.all(z['initial_energy_j']<0)
        assert np.all(np.diff(z['alive'].astype(int),axis=1)<=0)
    keys=set.intersection(*[set(z.files) for z in raw])
    result={k:np.concatenate([z[k] for z in raw],axis=0) for k in keys}
    counts=result['alive'].sum(axis=0);shots=len(result['alive']);s=counts/shots;lo,hi=wilson(counts,shots)
    return dict(meta=meta,raw=result,counts=counts,shots=shots,s=s,lo=lo,hi=hi)

def main():
    frozen=json.loads((HERE/'final_lock.json').read_text(encoding='utf-8'))
    groups={'Baseline':frozen['baseline_ids'],'Corrected formulas / controls':frozen['corrected_ids'],
            'Joint model':frozen['validation_ids']}
    colors=['#7B8794','#D18A29','#136FAD']
    fits=np.genfromtxt(ROOT/'reference/fig6d_vector_v2/published_fit_reference.csv',delimiter=',',names=True)
    metrics=[];csvrows=[];resrows=[];hazards=[];energyrows=[];checks=[];stagestats=[];convergence=[]
    overview,axes=plt.subplots(2,2,figsize=(12.6,9.2))
    hfig,hax=plt.subplots(2,2,figsize=(12,8))
    efig,eax=plt.subplots(2,2,figsize=(12,8))
    joint=[]
    for wi,(wave,label) in enumerate(zip(WAVES,LABELS)):
        data={k:load_group(v,wave) for k,v in groups.items()}
        exp=[p for p in POINTS if p['waveform']==wave]
        nx=np.array([int(p['n']) for p in exp]);sy=np.array([float(p['survival']) for p in exp])
        el=np.array([float(p['error_lower']) for p in exp]);eh=np.array([float(p['error_upper']) for p in exp])
        valid=np.isfinite(fits['survival_'+wave]);nf=fits['n_one_way_transfers'][valid].astype(int);yf=fits['survival_'+wave][valid]
        for group,d in data.items():
            s=d['s'];r=s[nx]-sy;rf=s[nf]-yf
            metrics.append(dict(group=group,waveform=wave,shots=d['shots'],rmse_pp=100*np.sqrt(np.mean(r*r)),
                max_abs_pp=100*max(abs(r)),early_bias_pp=100*np.mean(r[nx<=20]),late_bias_pp=100*np.mean(r[nx>20]),
                early_rmse_pp=100*np.sqrt(np.mean(r[nx<=20]**2)),late_rmse_pp=100*np.sqrt(np.mean(r[nx>20]**2)),
                fit_rmse_pp=100*np.sqrt(np.mean(rf*rf)),fit_max_pp=100*max(abs(rf)),
                S2=s[2],S10=s[10],S30=s[30],S60=s[60]))
            for i in range(61):
                csvrows.append(dict(group=group,waveform=wave,n=i,shots=d['shots'],survivors=int(d['counts'][i]),
                    survival=s[i],wilson95_low=d['lo'][i],wilson95_high=d['hi'][i],
                    published_fit=float(yf[list(nf).index(i)]) if i in nf else ''))
            for n,y,delta,l,u in zip(nx,sy,r,el,eh):
                resrows.append(dict(group=group,waveform=wave,n=n,experiment=y,mc=s[n],residual_pp=100*delta,
                                    paper_error_lower=l,paper_error_upper=u))
            m=d['meta'][0]
            checks.append(dict(group=group,waveform=wave,initial_unbound=0,recount=True,monotone=True,
                initial_hash=m['initial_pos_sha256'],velocity_hash=m['initial_vel_sha256']))
        def draw(ax):
            for (group,d),color in zip(data.items(),colors):
                n=np.arange(61);ax.fill_between(n,d['lo'],d['hi'],color=color,alpha=.14,lw=0)
                ax.plot(n,d['s'],color=color,lw=1.8,label=group)
            ax.plot(nf,yf,color='#252525',ls='--',lw=1,label='Published fitted line')
            ax.errorbar(nx,sy,yerr=np.maximum(0,np.vstack([sy-el,eh-sy])),fmt='D',ms=4,
                color='black',ecolor='black',capsize=2,label='Published markers')
            ax.set(title=label,xlabel='Number of one-way transfers',ylabel='Survival',
                   xlim=(0,32 if wi==0 else 61),ylim=(0,1.035))
            ax.yaxis.set_major_formatter(PercentFormatter(1));ax.grid(alpha=.14)
        draw(axes.flat[wi]);fig,ax=plt.subplots(figsize=(8.2,5.7));draw(ax)
        ax.legend(loc='lower left',fontsize=8,frameon=False)
        fig.text(.08,.016,'Shading: pointwise Wilson 95% Monte Carlo intervals; excludes model and timestep error.',fontsize=8)
        fig.tight_layout(rect=(0,.045,1,1));fig.savefig(HERE/f'fig6d_{wave}.png',dpi=190);fig.savefig(HERE/f'fig6d_{wave}.svg');plt.close(fig)
        d=data['Joint model'];s=d['s'];z=d['raw']
        hazard=np.divide(d['counts'][:-1]-d['counts'][1:],d['counts'][:-1],out=np.full(60,np.nan),where=d['counts'][:-1]>0)
        ax=hax.flat[wi];ax.plot(np.arange(1,61),hazard*100,color=colors[-1],alpha=.3,lw=.8,label='Each observation')
        for (group,gd),color in zip(data.items(),colors):
            # Fixed two-transfer conditional loss: actual experimental observation cadence.
            hs=1-np.sqrt(gd['s'][2::2]/gd['s'][:-2:2])
            ax.plot(np.arange(2,61,2),100*hs,'o-',ms=2,color=color,label=group)
        erate=1-(sy[1:]/sy[:-1])**(1/np.diff(nx))
        ax.scatter(nx[1:],100*erate,label='Paper interval geometric rate',marker='x',color='black')
        ax.set(title=label,xlabel='One-way transfer n',ylabel='Conditional loss (%)',xlim=(0,32 if wi==0 else 61));ax.grid(alpha=.15)
        for i,h in enumerate(hazard,1):hazards.append(dict(waveform=wave,n=i,conditional_loss=h,at_risk=int(d['counts'][i-1])))
        # Compare excitation in the same SLM after each completed round, survivors only.
        e=(z['energy_j']+140*1.380649e-29*z['initial_slm_depth_factor'][:,None])/1.380649e-29
        for n in range(0,61,2):
            es=e[z['alive'][:,n],n]
            qs=np.quantile(es,[.1,.5,.9]) if len(es) else [np.nan]*3
            energyrows.append(dict(waveform=wave,n=n,survivors=len(es),q10_uK=qs[0],q50_uK=qs[1],q90_uK=qs[2],
                                   mean_uK=np.mean(es) if len(es) else np.nan))
        er=[r for r in energyrows if r['waveform']==wave];nn=[r['n'] for r in er]
        ax=eax.flat[wi];ax.fill_between(nn,[r['q10_uK'] for r in er],[r['q90_uK'] for r in er],alpha=.2,label='10-90% of surviving atoms')
        ax.plot(nn,[r['q50_uK'] for r in er],label='Median among survivors')
        ax.plot(nn,[r['mean_uK'] for r in er],'--',label='Mean among survivors')
        ax.set(title=label,xlabel='One-way transfers (SLM observations)',ylabel='Excitation E/kB (uK)');ax.grid(alpha=.15)
        if 'stage_energy_j' in z:
            stage=z['stage_energy_j'];switch=z['switch_work_j'];lost=z['lost_round_stage']
            waitloss=lost[:,1]==1;ids=np.where(waitloss)[0];rr=lost[ids,0]
            at_start=stage[ids,rr,0]+switch[ids,rr,1]
            stagestats.append(dict(waveform=wave,wait_end_recorded_loss=int(waitloss.sum()),
                already_positive_after_switch=int(np.sum(at_start>=0)),
                fraction_positive_before_wait=float(np.mean(at_start>=0)) if len(ids) else None,
                drop_end_recorded_loss=int(np.sum(lost[:,1]==2)),
                mean_realized_noise_work_uK=float(z['noise_ledger'][:,0].mean()/1.380649e-29),
                mean_expected_noise_work_uK=float(z['noise_ledger'][:,1].mean()/1.380649e-29)))
        if frozen.get('refinement_ids'):
            fine=load_group(frozen['refinement_ids'],wave)
            coarse=load_group([frozen['validation_ids'][0]],wave)
            dif=fine['raw']['alive'].astype(float)-coarse['raw']['alive'].astype(float)
            delta=dif.mean(axis=0);worst=np.argmax(abs(delta))
            convergence.append(dict(waveform=wave,kind='mechanical_dt_same_noise_increments',
                coarse_dt_us=coarse['meta'][0]['case']['dt_us'],fine_dt_us=fine['meta'][0]['case']['dt_us'],
                max_abs_delta_pp=100*max(abs(delta)),worst_n=int(worst),
                paired_se_pp=100*dif[:,worst].std(ddof=1)/np.sqrt(len(dif)),
                final_label_disagreement=float(np.mean(dif[:,-1]!=0))))
    # Verify actual arrays, not merely recorded hashes, across four waveforms in every formal run.
    for ids in groups.values():
        for caseid in ids:
            z0=np.load(HERE/'runs'/caseid/(WAVES[0]+'.npz'))
            for wave in WAVES[1:]:
                z=np.load(HERE/'runs'/caseid/(wave+'.npz'))
                for key in ('initial_pos','initial_vel','initial_energy_j','initial_slm_depth_factor'):
                    assert np.array_equal(z[key],z0[key]),(caseid,wave,key)
    handles,labels=axes[0,0].get_legend_handles_labels()
    overview.legend(handles,labels,loc='lower center',ncol=3,fontsize=9,frameon=False)
    overview.text(.5,.075,'Joint model N=8192; controls N=4096. Shading: pointwise Wilson 95% Monte Carlo intervals.',
                  ha='center',fontsize=9)
    overview.suptitle('Fig. 6d | Shared preparation and shared apparatus parameters',fontsize=16,y=.99)
    overview.tight_layout(rect=(0,.11,1,.97));overview.savefig(HERE/'fig6d_overview.png',dpi=190);overview.savefig(HERE/'fig6d_overview.svg')
    for fig,axs,name,title in [(hfig,hax,'conditional_loss','Loss conditional on survival to the previous observation'),
                               (efig,eax,'energy_evolution','SLM excitation distribution; selection by survival is included')]:
        axs[0,0].legend(fontsize=7);fig.suptitle(title);fig.tight_layout(rect=(0,0,1,.96));fig.savefig(HERE/(name+'.png'),dpi=170)
    if (HERE/'extra_timestep_check.json').exists():
        convergence.append(json.loads((HERE/'extra_timestep_check.json').read_text(encoding='utf-8')))
    for filename,rows in [('metrics.csv',metrics),('survival_curves.csv',csvrows),('point_residuals.csv',resrows),
                         ('conditional_loss.csv',hazards),('energy_evolution.csv',energyrows),('timestep_check.csv',convergence),
                         ('loss_recording_diagnostics.csv',stagestats)]:
        if rows:csvwrite(HERE/filename,rows)
    for group in groups:
        vals=[r for r in metrics if r['group']==group]
        joint.append(dict(group=group,joint_rmse_pp=float(np.sqrt(np.mean([r['rmse_pp']**2 for r in vals])))))
    scans=[]
    for p in (HERE/'runs').glob('*/summary.json'):
        r=json.loads(p.read_text(encoding='utf-8'));scans.append(dict(**r['case'],joint_rmse_pp=r['joint_rmse_pp'],curve_rmse_pp=r['curve_rmse_pp']))
    write(HERE/'all_cases.json',scans);write(HERE/'integrity.json',checks)
    write(HERE/'analysis_summary.json',dict(joint=joint,metrics=metrics,stages=stagestats,convergence=convergence))
    print(json.dumps(joint),flush=True)

if __name__=='__main__':main()
