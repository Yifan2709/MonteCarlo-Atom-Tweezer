"""Compact, reproducible console summary used to check report statements."""
from pathlib import Path
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parent
d=pd.read_csv(R/'all_points.csv');p=pd.read_csv(R/'independent_confirmation.csv');b=pd.read_csv(R/'boundaries.csv')
print('PURE SINGLE 510')
print(b[(b.distance_um==510)&(b.protocol=='A')&(b['repeat']==1)&(b.disorder==0)][['trajectory','tau_us','target','fastest_confirmed_time_us','maximum_confirmed_average_m_s']].to_string(index=False))
print('B SINGLE 510')
print(b[(b.distance_um==510)&(b.protocol=='B')&(b['repeat']==1)][['trajectory','tau_us','rin','target','fastest_confirmed_time_us']].to_string(index=False))
print('ZERO JERK TAU5 N50 510 CONFIRMED')
print(p[(p.distance_um==510)&(p.trajectory=='zero_jerk')&(p['repeat']==50)][['protocol','duration_us','rin','shots','survival','lower','upper','status99','status999']].sort_values(['protocol','rin','duration_us']).to_string(index=False))
print('MOVING VS NO-MOVE HANDOFF')
print(d[(d.group=='handoff_controls')&(d['repeat']==50)][['distance_um','duration_us','alignment_nm','ramp_us','survival','lower','upper']].sort_values(['duration_us','ramp_us','alignment_nm','distance_um']).to_string(index=False))
print('MECHANISMS MINJERK 510')
print(d[(d.group=='coarse')&(d.distance_um==510)&(d.depth_uK==280)&(d.trajectory=='min_jerk')&d.duration_us.isin([400,600,800,1000,1200])][['duration_us','tau_us','survival','peak_accel_m_s2']].sort_values(['duration_us','tau_us']).to_string(index=False))
print('N10 BASELINE EDGE ROBUSTNESS')
edge=d[(d.group=='edge_robustness')&(d['repeat']==10)&(d.temperature_uK==19)&(d.waist==1.11e-6)&(d.disorder==0)]
print(edge.survival.describe().to_string())
print('CONVERGENCE MAX')
c=pd.read_csv(R/'mechanical_convergence.csv')
for n in [1,10,50]:
    g=c[(c['repeat']==n)&(c.seed==26091731)]
    if len(g):print(n,g[['survival_difference','label_disagreement','energy_KS']].abs().max().to_dict())
g=c[c.seed==26091771]
print('single final boundaries',g[['survival_difference','label_disagreement','energy_KS']].abs().max().to_dict())
