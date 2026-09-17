"""Freeze selection before independent validation; never select on validation data."""
from run_study import *
from datetime import datetime, timezone
import subprocess

if __name__ == '__main__':
    destination=HERE/'final_lock.json'
    if destination.exists():raise RuntimeError('Selection is already locked')
    queue=json.loads((HERE/'selection_queue.json').read_text(encoding='utf-8'))
    candidates=[json.loads((RUNS/c['id']/'summary.json').read_text(encoding='utf-8')) for c in queue]
    best=min(candidates,key=lambda c:c['joint_rmse_pp'])
    parameters={k:v for k,v in best['case'].items() if k not in ('id','shots','seed','dt_us','noise_us')}
    assert not any((RUNS/i/'summary.json').exists() for i in ['validation_seed21','validation_seed22'])
    final=[dict(id='formal_baseline',model='legacy',manual='legacy_piecewise',shots=4096,
                seed=26091621,dt_us=.0125),
           dict(id='formal_corrected',model='corrected_constant',manual='paper_cubic',
                ml_vector=True,shots=4096,seed=26091621,dt_us=.0125)]
    final.extend(dict(parameters,id=f'validation_seed{s}',shots=4096,seed=26091600+s,
                      dt_us=.0125,noise_us=.05) for s in (21,22))
    checks=[dict(parameters,id='refinement_seed21',shots=4096,seed=26091621,dt_us=.00625,noise_us=.05),
            dict(parameters,id='noise_half_seed23',shots=4096,seed=26091623,dt_us=.0125,noise_us=.025),
            dict(parameters,id='noise_reference_seed23',shots=4096,seed=26091623,dt_us=.0125,noise_us=.05)]
    changes=[('reference',{}),('jitter_off',dict(jitter_nm=0)),('alignment_off',dict(alignment_nm=0)),
             ('aod_disorder_off',dict(no_aod_disorder=True)),('pointing_off',dict(pointing=0)),
             ('recoil_off',dict(recoil=False)),('slm_small',dict(slm_rin=1e-12,slm_pointing=1e-26)),
             ('ml_original',dict(ml_vector=False)),('ml_clamped',dict(ml_endpoints=True)),
             ('lens_old',dict(vs=SETTINGS['vs_m_s'])),
             ('lens_filled',dict(vs=.412425*parameters.get('waist_um',1.17)/1.17)),
             ('mean_only',dict(model='local_mean')),('rin_off',dict(rin=0))]
    sensitivity=[dict(parameters,**dict(delta,id='sensitivity_'+name,shots=1024,seed=26091612,
                                       dt_us=.0125,noise_us=.05)) for name,delta in changes]
    write(destination,dict(locked_utc=datetime.now(timezone.utc).isoformat(),parameters=parameters,
        selected=best,selection_candidates=candidates,selection_rule='Minimum mean per-curve marker MSE; equal four curves',
        baseline_ids=['formal_baseline'],corrected_ids=['formal_corrected'],
        validation_ids=['validation_seed21','validation_seed22'],refinement_ids=['refinement_seed21'],
        mechanical_dt_us=.0125,noise_step_us=.05,validation_seed_policy='Never reselect physical parameters after these runs',
        noise_step_check='Compare h=.05 and .025 us; fresh seed23, 4096 atoms, same initial pool, distinct noise paths',
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=WORK,text=True).strip(),
        source_sha256=SOURCE_HASHES))
    write(HERE/'final_queue.json',final)
    write(HERE/'validation_checks_queue.json',checks)
    write(HERE/'sensitivity_queue.json',sensitivity)
    (HERE/'final_model.yaml').write_text(yaml.safe_dump(final[2],sort_keys=False),encoding='utf-8')
    print(json.dumps(best),flush=True)
