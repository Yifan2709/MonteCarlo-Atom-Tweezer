"""Check reference immutability, shared configuration and deliverable completeness."""
from run_study import *
import subprocess

if __name__=='__main__':
    lock=json.loads((HERE/'final_lock.json').read_text(encoding='utf-8'))
    base=(HERE/'base_commit.txt').read_text().strip()
    references=[ROOT/'reference/fig6d_vector_v2/experimental_marker_centers.csv',
                ROOT/'reference/fig6d_vector_v2/published_fit_reference.csv',ROOT/SETTINGS['ml_waveform'],
                WORK/'2403.12021v4.pdf']
    for p in references:
        old=subprocess.check_output(['git','show',base+':'+p.relative_to(WORK).as_posix()],cwd=WORK)
        now=p.read_bytes()
        if p.suffix!='.pdf':old=old.replace(b'\r\n',b'\n');now=now.replace(b'\r\n',b'\n')
        assert old==now,str(p)
    for key in ('baseline_ids','corrected_ids','validation_ids','refinement_ids'):
        for caseid in lock[key]:
            rows=[json.loads((RUNS/caseid/(w+'.json')).read_text()) for w in WAVES]
            assert all(r['case']==rows[0]['case'] for r in rows)
            if key=='validation_ids':
                for r in rows:
                    assert all(r['case'][k]==v for k,v in lock['parameters'].items())
                    assert r['Tref_used'] is False and r['fixed_power_uK_s'] is None
    for queue in ['final_queue.json','validation_checks_queue.json','sensitivity_queue.json']:
        for c in json.loads((HERE/queue).read_text()):
            summary=RUNS/c['id']/'summary.json'
            assert summary.exists(),str(summary)
            actual=json.loads(summary.read_text())['case']
            assert all(actual[k]==v for k,v in c.items())
    names=['REPORT.md','PARAMETERS.md','PHYSICS.md','PROTOCOL.md','CONFIGURATION.md','MECHANISMS.md','VALIDATION.md',
           'final_lock.json','final_model.yaml','source_manifest.json','source_final.zip','raw_inventory.csv',
           'survival_curves.csv','point_residuals.csv','metrics.csv','scan_history.csv','extra_timestep_check.json',
           'noise_step_check.csv','independent_seed_check.csv','mc_bootstrap.json','sensitivity_results.csv',
           'conditional_loss.csv','energy_evolution.csv','timestep_check.csv','reproduce.ps1','requirements_exact.txt',
           'fig6d_overview.png','fig6d_overview.svg']
    names += [f'fig6d_{w}.{ext}' for w in WAVES for ext in ['png','svg']]
    for name in names:assert (HERE/name).stat().st_size>0,name
    checked=dict(reference_files_unchanged=True,shared_parameters=True,all_queued_runs_complete=True,
                 deliverable_count=len(names),source_commit=lock['source_commit'])
    write(HERE/'delivery_check.json',checked);print(json.dumps(checked))
