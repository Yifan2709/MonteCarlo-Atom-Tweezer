"""Explicit Fig6d entry; all controls shared except published waveform/duration."""
from pathlib import Path
import sys,os,json,time,hashlib,argparse,csv,platform
from dataclasses import replace,asdict
HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]
ROOT=WORK/'simulation/level2_joint_transfer'
RUNS=HERE/'runs'
sys.path[:0]=[str(ROOT/'src'),str(WORK/'tmp/audit_deps'),str(WORK/'tmp/continuous_deps')]
os.environ.setdefault('NUMBA_NUM_THREADS','12')
import numpy as np
import yaml,numba
from continuous_transfer.cli import case_inputs,clean,wilson
from continuous_transfer.protocol import waveform_for,build_protocol
from continuous_transfer.engine import run_continuous
from continuous_transfer.noisy_survival import run_noisy_survival
from continuous_transfer.initialization import sample_site_bound_harmonic,slm_energies
from level2_joint_transfer.noise_heating import recoil_energy_j,recoil_scattering_rate_cs
WAVES=['manual_200us','manual_400us','manual_600us','ml_400us']
SETTINGS=yaml.safe_load((ROOT/'configs/continuous_fig6d_unified_T15.yaml').read_text(encoding='utf-8'))
SOURCE_HASHES={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'src').rglob('*.py')}
with (ROOT/'reference/fig6d_vector_v2/experimental_marker_centers.csv').open(encoding='utf-8-sig') as f:
    POINTS=list(csv.DictReader(f))

def write(path,obj):
    path.write_text(json.dumps(clean(obj),ensure_ascii=False,indent=2,default=str),encoding='utf-8')

def digest(a):return hashlib.sha256(np.asarray(a).tobytes()).hexdigest()

def evaluate(case):
    allowed={'id','temperature','rin','pointing','vs','jitter_nm','alignment_nm','manual','model',
             'shots','seed','dt_us','noise_us','waist_um','aod_waist_um','no_aod_disorder',
             'ml_endpoints','ml_vector','recoil','slm_rin','slm_pointing','absorb_at_checkpoints'}
    if set(case)-allowed:raise ValueError(f'Unknown case keys: {set(case)-allowed}')
    case={**dict(temperature=15.,rin=1e-8,pointing=1e-24,vs=SETTINGS['vs_m_s'],
        jitter_nm=10.,alignment_nm=50.,manual='paper_cubic',model='diffusion',
        shots=256,seed=26091611,dt_us=.025,noise_us=.05),**case}
    out=RUNS/case['id'];out.mkdir(parents=True,exist_ok=True)
    existing=out/'summary.json'
    if existing.exists():
        saved=json.loads(existing.read_text(encoding='utf-8'))
        if saved['case'] != case: raise RuntimeError('Case id reused with different configuration')
        return saved
    cfg,physics,pos,vel,draws,powers,sampling=case_inputs(SETTINGS,case['temperature'],case['shots'],case['seed'],'matched')
    # Rebuild the same per-site initialization if (and only if) common waist changes.
    if 'waist_um' in case:
        physics=dict(physics,slm_waist_m=case['waist_um']*1e-6,aod_waist_m=case['waist_um']*1e-6)
        pool=sample_site_bound_harmonic(case['seed'],case['temperature'],physics,draws.slm_depth_factor)
        pos,vel=pool['pos_m'],pool['vel_m_per_s']
    if 'aod_waist_um' in case:physics=dict(physics,aod_waist_m=case['aod_waist_um']*1e-6)
    draws=replace(draws,alignment_offset_m=draws.alignment_offset_m*(case['alignment_nm']/50))
    if case.get('no_aod_disorder',False):
        draws=replace(draws,aod_depth_factor=np.ones(case['shots']),aod_waist_factor=np.ones(case['shots']))
    e0=slm_energies(pos,vel,physics,draws.slm_depth_factor)
    assert np.all(e0<0)
    rows=[]
    for wave in WAVES:
        dest=out/(wave+'.json')
        if dest.exists():
            rows.append(json.loads(dest.read_text(encoding='utf-8')));continue
        wf=waveform_for(wave,cfg,ROOT/SETTINGS['ml_waveform'],manual_profile=case['manual'])
        if wave=='ml_400us' and case.get('ml_vector',False):
            from level2_joint_transfer.waveforms import MLCubicWaveform
            raw=json.loads((HERE/'ml_vector_nodes.json').read_text(encoding='utf-8'))
            t=np.asarray(raw['grid_us'])*1e-6
            wf=MLCubicWaveform(t[-1],t/t[-1],np.asarray(raw['pos_um'])*1e-6,
                np.asarray(raw['depth_mK'])*1.380649e-26,2.4e-6,280*1.380649e-29)
        if wave=='ml_400us' and case.get('ml_endpoints',False):
            from scipy.interpolate import CubicSpline
            wf.center_values_m=wf.center_values_m.copy();wf.depth_values_j=wf.depth_values_j.copy()
            wf.center_values_m[[0,-1]]=[0,2.4e-6];wf.depth_values_j[[0,-1]]=[0,280*1.380649e-29]
            wf._center_interp=CubicSpline(wf.control_s,wf.center_values_m,bc_type=((1,0.),(1,0.)))
            wf._depth_interp=CubicSpline(wf.control_s,wf.depth_values_j)
        constant=case['model'] in ('legacy','corrected_constant','off')
        protocol=build_protocol(wf,100e-6,case['dt_us']*1e-6,dd=constant)
        started=time.perf_counter()
        power=powers['total']
        if case['model']=='legacy': power=powers['parametric']/4+powers['pointing']/2+powers['recoil']
        if case['model']=='off':power=0
        if constant:
            result=run_continuous(pos,vel,protocol,rounds=30,physics=physics,draws=draws,
                jitter_sigma_m=case['jitter_nm']*1e-9,power_w=power,vs_m_s=case['vs'],eta=0)
        else:
            # Small recoil term; local intensity weighting. Two-level Cs D2 approximation.
            rec=2*recoil_energy_j(1055e-9,physics['mass_kg'])*recoil_scattering_rate_cs(1055e-9,1.)
            result=run_noisy_survival(pos,vel,protocol,rounds=30,physics=physics,draws=draws,
                jitter_sigma_m=case['jitter_nm']*1e-9,vs_m_s=case['vs'],mode=case['model'],
                rin=case['rin'],pointing=case['pointing'],noise_step_s=case['noise_us']*1e-6,
                seed=case['seed']+1000,recoil_per_joule=rec if case.get('recoil',True) else 0.,
                slm_rin=case.get('slm_rin'),slm_pointing=case.get('slm_pointing'),
                absorb_at_checkpoints=case.get('absorb_at_checkpoints',True))
        elapsed=time.perf_counter()-started
        counts=result['alive'].sum(axis=0);s=counts/case['shots'];low,high=wilson(counts,case['shots'])
        points=[p for p in POINTS if p['waveform']==wave]
        nx=np.array([int(p['n']) for p in points]);yy=np.array([float(p['survival']) for p in points])
        residual=s[nx]-yy
        row=dict(waveform=wave,case=case,runtime_s=elapsed,n=list(range(61)),survival=s,
            counts=counts,wilson_low=low,wilson_high=high,point_n=nx,point_residual=residual,
            mse=float(np.mean(residual**2)),rmse_pp=100*np.sqrt(np.mean(residual**2)),
            early_bias_pp=100*np.mean(residual[nx<=20]),late_bias_pp=100*np.mean(residual[nx>20]),
            initial_unbound=0,initial_pos_sha256=digest(pos),initial_vel_sha256=digest(vel),
            effective_initial_temperature_uK=case['temperature'],
            Tref_used=constant,reference_spectrum_convention='one-sided per Hz',
            slm_sites_sha256=digest(draws.slm_depth_factor),physics=physics,
            full_base_config=asdict(cfg),full_settings=SETTINGS,
            fixed_power_uK_s=power/1.380649e-29 if constant else None,
            versions=dict(python=platform.python_version(),numpy=np.__version__,numba=numba.__version__),
            source_sha256=SOURCE_HASHES)
        np.savez_compressed(out/(wave+'.npz'),**result,initial_pos=pos,initial_vel=vel,
            initial_energy_j=e0,initial_slm_depth_factor=draws.slm_depth_factor,
            aod_depth_factor=draws.aod_depth_factor,aod_waist_factor=draws.aod_waist_factor,
            alignment_offset_m=draws.alignment_offset_m)
        write(dest,row);rows.append(clean(row))
        print(case['id'],wave,'rmse_pp',round(row['rmse_pp'],3),'S60',round(s[-1],4),
              'seconds',round(elapsed,1),flush=True)
    summary=dict(case=case,joint_rmse_pp=100*np.sqrt(np.mean([r['mse'] for r in rows])),
                 curve_rmse_pp=[r['rmse_pp'] for r in rows],s60=[r['survival'][-1] for r in rows])
    write(existing,summary);print('SUMMARY',json.dumps(summary),flush=True)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('queue');p.add_argument('--out-dir',type=Path);args=p.parse_args()
    if args.out_dir is not None:RUNS=args.out_dir.resolve()
    numba.set_num_threads(12)
    queue_path=Path(args.queue)
    content=queue_path.read_text(encoding='utf-8-sig')
    cases=json.loads(content) if queue_path.suffix.lower()=='.json' else yaml.safe_load(content)
    if isinstance(cases,dict):cases=[cases]
    for case in cases:evaluate(case)
