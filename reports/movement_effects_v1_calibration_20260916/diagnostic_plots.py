"""Display discrete screening results without treating them as confidence regions."""
from analyze_study import *

def main():
    fig,axs=plt.subplots(1,2,figsize=(11,5.3),layout='constrained')
    for ax,queue,title in zip(axs,['scan_queue.json','aod_queue.json'],
        ['Equal SLM/AOD PSD (excluded by static hold)','AOD PSD; quiet SLM candidate']):
        cases=[json.loads((RUNS/c['id']/'summary.json').read_text(encoding='utf-8'))
               for c in json.loads((HERE/queue).read_text())]
        ts=sorted({r['case']['temperature'] for r in cases});rs=sorted({r['case']['rin'] for r in cases})
        z=np.array([[next(r['joint_rmse_pp'] for r in cases if r['case']['temperature']==t and r['case']['rin']==s)
                     for s in rs] for t in ts])
        im=ax.imshow(z,origin='lower',cmap='viridis_r',vmin=0,vmax=30,aspect='auto')
        for i,t in enumerate(ts):
            for j,s in enumerate(rs):ax.text(j,i,f'{z[i,j]:.1f}',ha='center',va='center',fontsize=9,
                                           color='white' if z[i,j]>14 else 'black')
        ax.set(xticks=np.arange(len(rs)),xticklabels=[f'{r/1e-9:g}' for r in rs],
               yticks=np.arange(len(ts)),yticklabels=ts,xlabel='RIN PSD (1e-9 /Hz)',ylabel='Common T0 (uK)',title=title)
    fig.colorbar(im,ax=axs,label='Joint marker RMSE (percentage points)',shrink=.8)
    fig.suptitle('Discrete initial screens | N=256; waist=1.17 um; lensing off\nNot a confidence map; source waveforms differ as recorded',fontsize=12)
    fig.savefig(HERE/'parameter_screen.png',dpi=180)
    rows=[]
    for q in ['diagnostic_queue.json','aod_queue.json']:
        for c in json.loads((HERE/q).read_text()):
            for wave in WAVES:
                m=json.loads((RUNS/c['id']/(wave+'.json')).read_text())
                rows.append(dict(case=c['id'],waveform=wave,shots=m['case']['shots'],temperature=m['case']['temperature'],
                    rin=m['case']['rin'],rmse_pp=m['rmse_pp'],early_bias_pp=m['early_bias_pp'],late_bias_pp=m['late_bias_pp'],
                    s2=m['survival'][2],s10=m['survival'][10],s30=m['survival'][30],s60=m['survival'][60]))
    csvwrite(HERE/'diagnostic_residuals.csv',rows)

if __name__=='__main__':main()
