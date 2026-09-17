"""Target the largest dt discrepancy without reselecting model parameters."""
import run_study as s
import numpy as np
import json
import argparse

if __name__=='__main__':
    lock=json.loads((s.HERE/'final_lock.json').read_text(encoding='utf-8'))
    parser=argparse.ArgumentParser();parser.add_argument('--runs-dir',type=s.Path);args=parser.parse_args()
    base_runs=(args.runs_dir.resolve() if args.runs_dir else s.HERE/'runs')
    metadata=(base_runs if args.runs_dir else s.HERE)
    s.numba.set_num_threads(12)
    s.WAVES=['manual_600us'];s.RUNS=metadata/'timestep_extra'
    case=dict(lock['parameters'],id='fine600',shots=4096,seed=26091621,dt_us=.003125,noise_us=.05)
    s.write(metadata/'extra_timestep_config.json',case)
    s.evaluate(case)
    a=np.load(base_runs/'refinement_seed21/manual_600us.npz')['alive']
    b=np.load(s.RUNS/'fine600/manual_600us.npz')['alive']
    dif=b.astype(float)-a.astype(float);delta=dif.mean(axis=0);worst=np.argmax(abs(delta))
    s.write(metadata/'extra_timestep_check.json',dict(waveform='manual_600us',
        kind='additional_mechanical_dt_same_noise_increments',coarse_dt_us=.00625,fine_dt_us=.003125,
        max_abs_delta_pp=100*max(abs(delta)),worst_n=int(worst),
        paired_se_pp=100*dif[:,worst].std(ddof=1)/np.sqrt(len(dif)),
        final_label_disagreement=float(np.mean(dif[:,-1]!=0))))
