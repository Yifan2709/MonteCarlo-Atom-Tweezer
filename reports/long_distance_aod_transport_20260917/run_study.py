"""Resume-safe study queues. Each result retains full per-atom checkpoints."""
from pathlib import Path
from dataclasses import replace,asdict
import argparse,json,hashlib,time,sys,platform,subprocess
import numpy as np
import pandas as pd
from long_distance_transport.model import Config,run,KB,U
ROOT=Path(__file__).resolve().parent
TR=['rb_cubic','zero_jerk','min_jerk']
TIMES=[100,150,200,300,400,600,800,1000,1200,1600,2200,3000,4500,6000]

def cases(stage):
    base=Config()
    if stage=='native':
        # Rb measured frequency and depth determine waist. D=55 um per atom.
        mass=86.909180531*U;dep=6.62607015e-34*4e6
        rb=replace(base,species='Rb',mass=mass,wavelength=828e-9,
            waist=np.sqrt(4*dep/(mass*(2*np.pi*40e3)**2)),depth_uK=dep/KB*1e6,
            temperature_uK=10,distance_um=55,shots=2048,dt_us=.05)
        for T in np.arange(90,331,10):yield replace(rb,duration_us=float(T))
        # Yb native species and wavelength; waist/temp/tau remain explicit hypotheses.
        # Full-circuit source waveform has 2 MHz AOD move depth; not asserted ED2 depth.
        yb=replace(base,species='Yb',mass=170.9363316*U,wavelength=487e-9,
            waist=.77e-6,depth_uK=2e6*6.62607015e-34/KB*1e6,
            temperature_uK=5,distance_um=50,shots=1024,dt_us=.05)
        tab=pd.read_csv(ROOT/'zhang_source/dataset/extended_data_fig2b.csv')
        ts=tab.exp_moving_time_us_zero_jerk.dropna().values
        for tr in TR[1:]:
            for tau in [0.,2.,5.]:
                for T in ts:yield replace(yb,trajectory=tr,tau_us=tau,duration_us=float(T))
    elif stage=='native_steps':
        # Independent step refinement at points on the native transition curves.
        native=list(cases('native'))
        for species,tr,target in [('Rb','rb_cubic',120.),('Rb','rb_cubic',160.),
                                  ('Yb','zero_jerk',230.),('Yb','min_jerk',270.)]:
            choices=[c for c in native if c.species==species and c.trajectory==tr
                     and c.tau_us==(0. if species=='Rb' else 2.)]
            c=min(choices,key=lambda c:abs(c.duration_us-target))
            for dt in [.1,.05,.025]:
                yield replace(c,dt_us=dt,shots=4096,seed=26091761)
    elif stage=='coarse':
        for L in [270.,510.,610.]:
            for tr in TR:
                for tau in [0.,2.,5.]:
                    for T in TIMES:yield replace(base,distance_um=L,trajectory=tr,tau_us=tau,duration_us=T)
        for dep in [140.,560.]:
            for tau in [0.,5.]:
                for tr in TR:
                    for T in TIMES:yield replace(base,depth_uK=dep,tau_us=tau,trajectory=tr,duration_us=T)
    elif stage=='mechanisms':
        # No simultaneous fitting: paired changes in defined groups.
        for tr in TR:
            for tau in [0.,5.]:
                for T in [400,600,800,1000,1200,1600,2200,3000,4500,6000]:
                    for proto,noise,disorder in [('A',True,False),('B',False,False),('B',False,True),('B',True,True)]:
                        yield replace(base,trajectory=tr,tau_us=tau,duration_us=T,protocol=proto,
                            rin=5.75e-9 if noise else 0.,pointing=1e-24 if noise else 0.,
                            disorder=.02 if disorder else 0.,alignment_nm=50 if disorder else 0.)
        for geom in ['diagonal','opposite']:
            for tr in TR:
                for T in [400,600,800,1000,1200,1600,2200,3000]:
                    yield replace(base,geometry=geom,trajectory=tr,tau_us=5,duration_us=T)
    elif stage=='repeat':
        for tr in TR:
            for proto in ['A','B']:
                for noise in [False,True]:
                    for T in [600,800,1000,1200,1600,2200,3000,4500,6000]:
                        yield replace(base,trajectory=tr,protocol=proto,repeats=50,shots=512,
                            tau_us=5,duration_us=T,rin=5.75e-9 if noise else 0.,
                            pointing=1e-24 if noise else 0.,disorder=.02,alignment_nm=50)
    elif stage=='repeat_distances':
        for c in cases('repeat'):
            for L in [270.,610.]:
                yield replace(c,distance_um=L)
    elif stage=='repeat_steps':
        # Additional repeated-operation controls for the other two trajectory
        # families, beyond the full zero-jerk time/noise refinement grid.
        for tr,ts in [('rb_cubic',[2200.,3000.]),('min_jerk',[1600.,3000.])]:
            for T in ts:
                for dt in [.1,.05,.025]:
                    yield replace(base,trajectory=tr,duration_us=T,tau_us=5,
                        repeats=50,shots=2048,seed=26091773,dt_us=dt,
                        disorder=.02,alignment_nm=50)
    elif stage=='edge_robustness':
        df=pd.read_csv(ROOT/'independent_confirmation.csv')
        df=df[(df.protocol=='A')&(df['repeat']==1)&(df.disorder==0)&(df.status99=='supported')]
        for _,g in df.groupby(['distance_um','trajectory','tau_us']):
            r=g.sort_values('duration_us').iloc[0]
            c=replace(base,distance_um=float(r.distance_um),trajectory=r.trajectory,tau_us=float(r.tau_us),
                duration_us=float(r.duration_us),shots=1024,seed=26091791,repeats=10)
            yield c
            yield replace(c,disorder=.02)
            yield replace(c,temperature_uK=30)
            yield replace(c,waist=1.17e-6)
            yield replace(c,duration_us=c.duration_us*.99)
            yield replace(c,duration_us=c.duration_us*1.01)
            if c.tau_us:yield replace(c,tau_us=c.tau_us*1.2)
    elif stage=='handoff_controls':
        for duration in [1600.,6000.]:
            for offset in [0.,50.]:
                for L in [0.,510.]:
                    yield replace(base,protocol='B',trajectory='zero_jerk',tau_us=5,
                        duration_us=duration,distance_um=L,repeats=50,shots=2048,seed=26091781,
                        disorder=.02,alignment_nm=offset,dt_us=.1)
        for ramp in [400.,800.]:
            yield replace(base,protocol='B',trajectory='zero_jerk',tau_us=5,
                duration_us=1600,distance_um=0,repeats=50,shots=2048,seed=26091781,
                disorder=.02,alignment_nm=50,dt_us=.1,ramp_us=ramp)
    elif stage=='boundary_steps':
        df=pd.read_csv(ROOT/'independent_confirmation.csv')
        df=df[(df.protocol=='A')&(df['repeat']==1)&(df.disorder==0)&(df.status99=='supported')]
        for _,g in df.groupby(['distance_um','trajectory','tau_us']):
            r=g.sort_values('duration_us').iloc[0]
            for dt in [.1,.05,.025]:
                yield replace(base,distance_um=float(r.distance_um),trajectory=r.trajectory,tau_us=float(r.tau_us),duration_us=float(r.duration_us),dt_us=dt,shots=4096,seed=26091771)
    elif stage=='final_single':
        # Stronger confirmation: include the 99% edge, a margin beyond zero-loss
        # screening, and a slow point. Seeds 41/42 never used for physics selection.
        for L in [270.,510.,610.]:
            for tr in TR:
                for tau in [0.,2.,5.]:
                    rs=[]
                    for group in ['coarse','refine']:
                        for path in (ROOT/'runs'/group).glob('*.json'):
                            r=json.loads(path.read_text());c=r['config']
                            if c['distance_um']==L and c['trajectory']==tr and c['tau_us']==tau and c['depth_uK']==280:rs.append(r)
                    rs.sort(key=lambda r:r['config']['duration_us'])
                    edge=next(r for r in rs if r['observations'][0]['survival']>=.99)
                    zero=next(r for r in rs if r['observations'][0]['survival']==1.)
                    ts={edge['config']['duration_us'],float(round(zero['config']['duration_us']*1.2/5)*5),3000.,6000.}
                    for t in sorted(ts):
                        for seed in [26091741,26091742]:yield replace(base,distance_um=L,trajectory=tr,tau_us=tau,duration_us=t,shots=4096,seed=seed,dt_us=.05)
    elif stage=='operations_confirm':
        selected=[]
        # Final single-operation confirmations for every trajectory in B.
        for tr in TR:
            for tau in [0.,5.]:
                for noise in [False,True]:
                    rs=[]
                    for path in (ROOT/'runs/mechanisms').glob('*.json'):
                        r=json.loads(path.read_text());c=r['config']
                        if c['trajectory']==tr and c['tau_us']==tau and c['protocol']=='B' and c['disorder']==.02 and bool(c['rin'])==noise:
                            rs.append(r)
                    rs.sort(key=lambda r:r['config']['duration_us'])
                    good=[r for r in rs if r['observations'][0]['survival']>=.999]
                    if good:
                        ix=rs.index(good[0]);selected.extend(rs[max(0,ix-1):ix+1])
                    else:
                        selected.extend(sorted(rs,key=lambda r:r['observations'][0]['survival'],reverse=True)[:2])
        for r in selected:
            for seed in [26091721,26091722]:yield replace(Config(**r['config']),shots=4096,seed=seed,dt_us=.05)
    elif stage=='repeat_confirm':
        selected=[]
        for L in [270.,510.,610.]:
            for tr in TR:
                for proto in ['A','B']:
                    rs=[]
                    for group in ['repeat','repeat_distances']:
                        for path in (ROOT/'runs'/group).glob('*.json'):
                            r=json.loads(path.read_text());c=r['config']
                            if c['distance_um']==L and c['trajectory']==tr and c['protocol']==proto and c['rin']==0:rs.append(r)
                    rs.sort(key=lambda r:r['config']['duration_us'])
                    indices=set()
                    for n in [10,50]:
                        good=[r for r in rs if r['observations'][n-1]['survival']>=.999]
                        if good:indices.add(rs.index(good[0]))
                    if not indices and rs:indices.add(int(np.argmax([r['observations'][-1]['survival'] for r in rs])))
                    for ix in indices:selected.append(rs[ix])
        for r in selected:
            for seed in [26091721,26091722]:yield replace(Config(**r['config']),shots=4096,seed=seed,dt_us=.05)
    elif stage=='sensitivity':
        for tr in TR:
            for T in [600,1000,1600,2200,3000]:
                nominal=replace(base,trajectory=tr,tau_us=5,duration_us=T,protocol='B',repeats=10,shots=512)
                for temp in [10.,19.,30.]:yield replace(nominal,temperature_uK=temp)
                for waist in [1.17e-6,1.23e-6]:yield replace(nominal,waist=waist)
                for tau in [2.,3.5,6.5]:yield replace(nominal,tau_us=tau)
                for ramp in [100.,400.]:yield replace(nominal,ramp_us=ramp)
    elif stage=='refine':
        # All transitions near 99% including nonmonotonic edges, not binary search.
        for L in [270.,510.,610.]:
            for tr in TR:
                for tau in [0.,2.,5.]:
                    rows=[]
                    for path in (ROOT/'runs/coarse').glob('*.json'):
                        r=json.loads(path.read_text());c=r['config']
                        if c['distance_um']==L and c['trajectory']==tr and c['tau_us']==tau and c['depth_uK']==280:
                            rows.append((c['duration_us'],r['observations'][-1]['survival']))
                    rows.sort()
                    for (t0,p0),(t1,p1) in zip(rows[:-1],rows[1:]):
                        if (p0-.99)*(p1-.99)<=0 or (.9<min(p0,p1)<.999):
                            for t in np.linspace(t0,t1,5)[1:-1]:
                                yield replace(base,distance_um=L,trajectory=tr,tau_us=tau,duration_us=float(t),shots=2048)
    elif stage in ('confirm','convergence'):
        # Queue is built once from existing screening results; parameters not fitted.
        queue=ROOT/'confirmation_selection.json'
        if not queue.exists():
            selected=[]
            for L in [270.,510.,610.]:
                for tr in TR:
                    for tau in [0.,2.,5.]:
                        matches=[]
                        for group in ['coarse','refine']:
                            for path in (ROOT/'runs'/group).glob('*.json'):
                                r=json.loads(path.read_text());c=r['config']
                                if c['distance_um']==L and c['trajectory']==tr and c['tau_us']==tau and c['depth_uK']==280:
                                    matches.append(r)
                        matches.sort(key=lambda r:r['config']['duration_us'])
                        good=[r for r in matches if r['observations'][-1]['survival']>=.999]
                        if good:
                            r=good[0];idx=matches.index(r)
                            selected.append(r['config'])
                            if idx:selected.append(matches[idx-1]['config'])
            # Confirm multiple operational protocols at common representative times.
            for proto in ['A','B']:
                for noise in [False,True]:
                    for T in [1600.,3000.,6000.]:
                        selected.append(asdict(replace(base,protocol=proto,repeats=50,tau_us=5,
                            duration_us=T,trajectory='zero_jerk',rin=5.75e-9 if noise else 0.,
                            pointing=1e-24 if noise else 0.,disorder=.02,alignment_nm=50)))
            queue.write_text(json.dumps(selected,indent=2))
        selected=json.loads(queue.read_text())
        if stage=='confirm':
            for val in selected:
                for seed in [26091721,26091722]:yield replace(Config(**val),shots=4096,seed=seed,dt_us=.05)
        else:
            # At least one near-boundary condition per trajectory and lens scenario.
            filtered=[v for v in selected if v['distance_um']==510 and v['repeats']==1]
            filtered+=selected[-12:]
            for val in filtered:
                for dt in [.1,.05,.025]:yield replace(Config(**val),shots=2048,seed=26091731,dt_us=dt)
            for val in selected[-12:]:
                if val['rin']:
                    for noise_dt in [.2,.1,.05]:yield replace(Config(**val),shots=2048,seed=26091732,dt_us=.025,noise_us=noise_dt)
    else:raise ValueError(stage)

def execute(stage,start=0,stop=None):
    checks=json.loads((ROOT/'benchmarks.json').read_text())
    if not checks['passed']:raise RuntimeError('benchmarks required')
    out=ROOT/'runs'/stage;out.mkdir(parents=True,exist_ok=True)
    configs=list(cases(stage))
    physics_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path('simulation/level2_joint_transfer/src/long_distance_transport/model.py'),Path('simulation/level2_joint_transfer/src/paper_integration/exact_trajectories.py')]}
    (ROOT/(stage+'_queue.json')).write_text(json.dumps([asdict(c) for c in configs],indent=2))
    for i,c in enumerate(configs):
        if i<start or (stop is not None and i>=stop):continue
        cfg=asdict(c);key=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()[:18]
        dest=out/(key+'.json')
        if dest.exists():
            previous=json.loads(dest.read_text())
            if any(previous['source_sha256'].get(k)!=v for k,v in physics_hashes.items()):
                raise RuntimeError(f'Physics changed; preserve old output and use a fresh study directory: {dest}')
            continue
        t=time.perf_counter();result,raw=run(c)
        np.savez_compressed(out/(key+'.npz'),**raw)
        result['runtime_s']=time.perf_counter()-t
        result['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path('simulation/level2_joint_transfer/src/long_distance_transport/model.py'),Path('simulation/level2_joint_transfer/src/paper_integration/exact_trajectories.py')]}
        dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
        last=result['observations'][-1]
        print(stage,i+1,len(configs),key,c.species,c.trajectory,c.protocol,c.distance_um,c.duration_us,c.tau_us,c.repeats,last['survival'],round(result['runtime_s'],2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stages',nargs='+')
    parser.add_argument('--start',type=int,default=0);parser.add_argument('--stop',type=int)
    args=parser.parse_args()
    for stage in args.stages:execute(stage,args.start,args.stop)
