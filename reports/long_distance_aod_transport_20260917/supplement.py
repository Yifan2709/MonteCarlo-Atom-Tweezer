"""Additional evidence tables. Reads simulations; never feeds curves into dynamics."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from long_distance_transport.model import interval
ROOT=Path(__file__).resolve().parent

def md(df):
    # No optional tabulate dependency.
    return '|'+ '|'.join(map(str,df.columns))+'|\n|'+ '|'.join(['---']*len(df.columns))+'|\n'+''.join('|'+ '|'.join(str(x) for x in row)+'|\n' for row in df.itertuples(index=False,name=None))

def main():
    df=pd.read_csv(ROOT/'all_points.csv');pool=pd.read_csv(ROOT/'independent_confirmation.csv')
    b=pd.read_csv(ROOT/'boundaries.csv')
    pure=b[(b.protocol=='A')&(b['repeat']==1)&(b.disorder==0)]
    pure.to_csv(ROOT/'single_boundaries.csv',index=False)
    view=pure[['distance_um','trajectory','tau_us','target','fastest_confirmed_time_us','maximum_confirmed_average_m_s','maximum_confirmed_peak_m_s','minimum_complete_us']].round(6)
    text='# 全部条件性边界\n\n逐点95% Clopper–Pearson下界认证。只表示已计算时间点，不认证点间连续区间。参数未知与数值误差未纳入二项区间。\n\n## A 单程，无技术噪声与位点散布\n\n'+md(view)
    text+='\n## B 与重复 A/B\n\n条件含2%位点深度散布、50nm固定对准误差。重复单位为单程。完整累计操作时间=每程完整时间×重复次数。空值表示未找到统计支持点，不等于证明不存在。\n\n'
    ops=b[~((b.protocol=='A')&(b['repeat']==1)&(b.disorder==0))].copy()
    text+=md(ops[['distance_um','trajectory','protocol','rin','repeat','target','fastest_confirmed_time_us','maximum_confirmed_average_m_s','minimum_complete_us']].round(6).fillna('未支持'))
    text+='\n全部支持时间点列表见 [boundaries.csv](boundaries.csv)，逐点区间与未确定/失败点见 [independent_confirmation.csv](independent_confirmation.csv)。\n'
    (ROOT/'BOUNDARIES.md').write_text(text,encoding='utf-8')
    # Independent seeds are presented separately, not hidden by pooling.
    df[df.group.isin(['confirm','operations_confirm','repeat_confirm','final_single'])&df['repeat'].isin([1,10,50])].to_csv(ROOT/'seedwise_confirmation.csv',index=False)
    # Native paper controls use identical initial pools at all three timesteps.
    rows=[]
    native=df[df.group=='native_steps']
    for key,g in native.groupby(['species','trajectory','duration_us','tau_us']):
        g=g.sort_values('dt_us',ascending=False)
        for (_,a),(_,c) in zip(g.iloc[:-1].iterrows(),g.iloc[1:].iterrows()):
            x=np.load((ROOT/a.file).with_suffix('.npz'));y=np.load((ROOT/c.file).with_suffix('.npz'))
            aa=x['alive'][:,1];bb=y['alive'][:,1]
            rows.append(dict(zip(['species','trajectory','duration_us','tau_us'],key))|dict(dt_coarse=a.dt_us,dt_fine=c.dt_us,
                survival_difference=float(aa.mean()-bb.mean()),label_disagreement=float(np.mean(aa!=bb)),
                energy_KS=float(ks_2samp(x['energy_j'][aa,1],y['energy_j'][bb,1]).statistic)))
    pd.DataFrame(rows).to_csv(ROOT/'native_convergence.csv',index=False)
    nonmono=[]
    screening=df[df.group.isin(['coarse','refine','repeat','repeat_distances'])]
    keys=['distance_um','trajectory','protocol','tau_us','depth_uK','rin','disorder','repeat']
    for key,g in screening.groupby(keys):
        if key[-1] not in [1,10,50]:continue
        g=g.sort_values(['duration_us','shots']).drop_duplicates('duration_us',keep='last')
        for (_,a),(_,c) in zip(g.iloc[:-1].iterrows(),g.iloc[1:].iterrows()):
            if c.upper<a.lower:
                nonmono.append(dict(zip(keys,key))|dict(time_short_us=a.duration_us,time_long_us=c.duration_us,
                    survival_short=a.survival,survival_long=c.survival,short_lower=a.lower,long_upper=c.upper,
                    interpretation='screening evidence only; no familywise correction'))
    pd.DataFrame(nonmono).to_csv(ROOT/'nonmonotonic_screening.csv',index=False)
    # Robustness: fastest single-move 99% point, continuously repeated ten times.
    edge=df[(df.group=='edge_robustness')&df['repeat'].isin([1,10])].copy()
    edge.to_csv(ROOT/'edge_robustness.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for ax,L in zip(axs,[270,510,610]):
        g=edge[(edge.distance_um==L)&(edge['repeat']==10)]
        for tr,h in g.groupby('trajectory'):
            ax.scatter(h.average_m_s,h.survival,label=tr,alpha=.7)
        ax.set(xlabel='Average speed (m/s)',ylabel='10-move survival',title=f'{L} um; edge perturbations',ylim=(0,1.02));ax.grid(alpha=.2)
    axs[-1].legend(fontsize=8);fig.tight_layout();fig.savefig(ROOT/'edge_robustness.png');plt.close(fig)
    hc=df[(df.group=='handoff_controls')&df['repeat'].isin([1,10,50])]
    hc.to_csv(ROOT/'handoff_controls.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for ax,tr in zip(axs,['rb_cubic','zero_jerk','min_jerk']):
        for geom in ['x','diagonal','opposite']:
            g=df[(df.group.isin(['coarse','mechanisms']))&(df.geometry==geom)&(df.trajectory==tr)&(df.depth_uK==280)&(df.distance_um==510)&(df.protocol=='A')&(df.rin==0)&(df.tau_us==5)].sort_values('duration_us')
            ax.plot(g.duration_us,g.survival,'o-',label=geom)
        ax.set(xlabel='Move duration (us)',ylabel='Single-move survival',title=tr,xscale='log');ax.grid(alpha=.2)
    axs[-1].legend();fig.tight_layout();fig.savefig(ROOT/'geometry_survival.png');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for ax,tr in zip(axs,['rb_cubic','zero_jerk','min_jerk']):
        for depth in [140,280,560]:
            for tau in [0,5]:
                g=df[(df.group=='coarse')&(df.depth_uK==depth)&(df.distance_um==510)&(df.trajectory==tr)&(df.tau_us==tau)].sort_values('duration_us')
                ax.plot(g.peak_accel_m_s2/g.maximum_static_force_acceleration,g.survival,'o--' if tau else 'o-',label=f'D={depth}, tau={tau}')
        ax.set(xlabel='Peak a / static acrit',ylabel='Single-move survival',title=tr,xscale='log');ax.grid(alpha=.2)
    axs[-1].legend(fontsize=7);fig.tight_layout();fig.savefig(ROOT/'depth_acceleration.png');plt.close(fig)
    queue=[]
    for q in ROOT.glob('*_queue.json'):
        stage=q.name.removesuffix('_queue.json');n=len(json.loads(q.read_text()));finished=len(list((ROOT/'runs'/stage).glob('*.json')))
        queue.append(dict(stage=stage,requested=n,completed=finished,complete=finished==n))
    pd.DataFrame(queue).to_csv(ROOT/'queue_completion.csv',index=False)
    print('Supplementary evidence written. Queue completion:',all(r['complete'] for r in queue))

if __name__=='__main__':main()
