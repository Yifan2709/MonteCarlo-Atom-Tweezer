"""Audit every raw record and replay a frozen production condition.

This checks archive integrity and reproducibility, not experimental validity.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from long_distance_transport.model import Config,run
ROOT=Path(__file__).resolve().parent

def main():
    failures=[];count=0;atoms=0;versions={};physical_versions={}
    for queue in ROOT.glob('*_queue.json'):
        stage=queue.name.removesuffix('_queue.json')
        for cfg in json.loads(queue.read_text()):
            key=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()[:18]
            path=ROOT/'runs'/stage/(key+'.json')
            if not path.exists():failures.append((str(path),'incomplete queue'))
    for path in sorted((ROOT/'runs').glob('*/*.json')):
        result=json.loads(path.read_text());c=result['config'];count+=1;atoms+=c['shots']
        key=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()[:18]
        if key!=path.stem:failures.append((str(path),'configuration hash mismatch'))
        for name,digest in result['source_sha256'].items():
            versions.setdefault(name,set()).add(digest)
        with np.load(path.with_suffix('.npz')) as raw:
            a=raw['alive'];e=raw['energy_j'];s=raw['states']
            if a.shape!=(c['shots'],c['repeats']+1):failures.append((str(path),'shape'))
            if not a[:,0].all() or not (e[:,0]<0).all():failures.append((str(path),'initially unbound'))
            if (a[:,1:]&~a[:,:-1]).any():failures.append((str(path),'resurrected labels'))
            if not np.isfinite(s[a]).all():failures.append((str(path),'invalid bound state'))
            if not (e[a]<0).all():failures.append((str(path),'unbound success'))
            for o in result['observations']:
                if int(a[:,o['repeat']].sum())!=o['alive']:failures.append((str(path),'count mismatch'))
    for name in ['model.py','exact_trajectories.py']:
        src=Path('simulation/level2_joint_transfer/src')/('long_distance_transport' if name=='model.py' else 'paper_integration')/name
        current=hashlib.sha256(src.read_bytes()).hexdigest()
        physical_versions[name]=dict(current_sha256=current,recorded_sha256=sorted(versions[name]))
        if versions[name]!={current}:failures.append((str(src),'mixed physics versions'))
    path=sorted((ROOT/'runs/final_single').glob('*.json'))[0]
    record=json.loads(path.read_text());_,again=run(Config(**record['config']))
    with np.load(path.with_suffix('.npz')) as saved:
        replay={}
        for key in ['states','energy_j','alive','initial_pos','initial_vel','site_factors','offsets','attempts']:
            replay[key]=bool(np.array_equal(saved[key],again[key],equal_nan=True))
            if not replay[key]:failures.append((str(path),'replay mismatch '+key))
    out=dict(records_checked=count,atom_trials=atoms,failures=failures,
        physical_versions=physical_versions,runner_versions=sorted(versions['run_study.py']),
        runner_note='Queue-only extensions during execution; all final configurations retained in individual JSON and queue files.',
        replay_file=str(path.relative_to(ROOT)),replay_bitwise_equal=replay,
        scope='All raw records checked against metadata; bound-state labels and no resurrection; independent physics tests are separate.')
    (ROOT/'artifact_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))
    if failures:raise RuntimeError('Artifact verification failed')

if __name__=='__main__':main()
