"""Numerical checks and ablations; never feed these outputs back into selection."""
from analyze_study import *

def mdtable(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows)

def compare(a,b,wave,kind):
    aa=load_group([a],wave);bb=load_group([b],wave)
    delta=bb['raw']['alive'].astype(float)-aa['raw']['alive'].astype(float)
    mean=delta.mean(axis=0);worst=np.argmax(abs(mean))
    return dict(kind=kind,waveform=wave,reference=a,comparison=b,max_delta_pp=100*max(abs(mean)),
        worst_n=int(worst),delta_at_worst_pp=100*mean[worst],
        paired_se_at_worst_pp=100*delta[:,worst].std(ddof=1)/np.sqrt(len(delta)),
        final_delta_pp=100*mean[-1],final_label_disagreement=float(np.mean(delta[:,-1]!=0)))

def bootstrap(ids,repetitions=500):
    # Resample original atoms, preserving all n correlations and cross-curve common draws.
    zz=[[load_group([i],wave)['raw']['alive'] for wave in WAVES] for i in ids]
    exp=[[p for p in POINTS if p['waveform']==wave] for wave in WAVES]
    ns=[np.array([int(p['n']) for p in points]) for points in exp]
    ys=[np.array([float(p['survival']) for p in points]) for points in exp]
    rng=np.random.default_rng(26091641);out=[]
    for _ in range(repetitions):
        ss=np.zeros((4,61));total=0
        for cohort in zz:
            count=len(cohort[0]);idx=rng.integers(0,count,count);total+=count
            for j,alive in enumerate(cohort):ss[j]+=alive[idx].sum(axis=0)
        rms=np.array([100*np.sqrt(np.mean((s[n]/total-y)**2)) for s,n,y in zip(ss,ns,ys)])
        out.append([*rms,np.sqrt(np.mean(rms*rms))])
    out=np.array(out)
    return dict(repetitions=repetitions,seed=26091641,
        intervals={name:list(np.quantile(out[:,j],[.025,.5,.975])) for j,name in enumerate([*WAVES,'joint'])},
        meaning='Monte Carlo resampling of atoms; not parameter confidence or experimental likelihood')

def main():
    frozen=json.loads((HERE/'final_lock.json').read_text(encoding='utf-8'))
    noise=[compare('noise_reference_seed23','noise_half_seed23',w,'noise_h_halved') for w in WAVES]
    independent=[compare('validation_seed21','validation_seed22',w,'independent_seed') for w in WAVES]
    csvwrite(HERE/'noise_step_check.csv',noise);csvwrite(HERE/'independent_seed_check.csv',independent)
    cases=json.loads((HERE/'sensitivity_queue.json').read_text(encoding='utf-8'))
    sensitivity=[]
    for c in cases[1:]:
        for wave in WAVES:
            r=compare('sensitivity_reference',c['id'],wave,'sensitivity')
            m=json.loads((RUNS/c['id']/(wave+'.json')).read_text(encoding='utf-8'))
            r.update(rmse_pp=m['rmse_pp'],early_bias_pp=m['early_bias_pp'],late_bias_pp=m['late_bias_pp'])
            sensitivity.append(r)
    csvwrite(HERE/'sensitivity_results.csv',sensitivity)
    preparation=[]
    for caseid in frozen['validation_ids']:
        d=load_group([caseid],WAVES[0]);z=d['raw'];mass=d['meta'][0]['physics']['mass_kg']
        excitation=z['initial_energy_j']/1.380649e-29+140*z['initial_slm_depth_factor']
        preparation.append(dict(id=caseid,shots=d['shots'],unbound=int(np.sum(z['initial_energy_j']>=0)),
            proposed_T_uK=d['meta'][0]['case']['temperature'],
            kinetic_T_uK=float(mass*np.mean(np.sum(z['initial_vel']**2,axis=1))/(3*1.380649e-29)),
            excitation_quantiles_uK=list(np.quantile(excitation,[.1,.5,.9,.99]))))
    boot=bootstrap(frozen['validation_ids']);write(HERE/'mc_bootstrap.json',boot)
    write(HERE/'validation_summary.json',dict(noise=noise,independent=independent,preparation=preparation,bootstrap=boot))
    def comparisons(rows):return mdtable(['曲线','最大存活率差/pp','该点n','差值配对SE/pp','末点差/pp'],[
        [r['waveform'],f"{r['max_delta_pp']:.3f}",r['worst_n'],f"{r['paired_se_at_worst_pp']:.3f}",f"{r['final_delta_pp']:+.3f}"] for r in rows])
    srows=[]
    for c in cases:
        m=json.loads((RUNS/c['id']/'summary.json').read_text(encoding='utf-8'))
        srows.append([c['id'],*[f'{x:.2f}' for x in m['curve_rmse_pp']],f"{m['joint_rmse_pp']:.2f}"])
    text='# 锁定后数值与敏感性验证\n\n所有本页结果均未用于重新选择物理参数。\n\n'
    text+='## 噪声步长\n\nh=0.05→0.025 us，机械dt固定0.0125 us，N=4096，独立于选参的seed26091623。初态相同，但噪声路径没有作Brownian bridge耦合，因此这是分布比较；个体轨迹不同是预期结果。配对SE按每个初始样本的差值计算，最大差值选择本身有多重比较，不把其SE当作预先指定点的显著性检验。\n\n'+comparisons(noise)
    text+='\n## 两个新种子\n\n各4096原子，所有点均保留；差异以MC误差评估，不选择表现较好的种子。\n\n'+comparisons(independent)
    text+='\n最终合并结果逐点Wilson95%区间最大半宽约1.08 pp（N=8192、p≈0.5）。区间只描述本计算中独立原子/噪声实现的抽样误差；真实实验同一光路的同时原子可共享噪声，且位点/重复测量可能相关，故这些区间不是实验样本误差。\n'
    text+='\nMC bootstrap保留同一原子的全部n相关性及四曲线共同初态，500次重抽；只估计有限MC样本引起的指标波动：\n\n'
    text+=mdtable(['曲线/指标','2.5%','50%','97.5%'],[[k,*[f'{x:.3f}' for x in v]] for k,v in boot['intervals'].items()])
    text+='\n## 初态\n\n下表Tkin=m〈v²〉/(3kB)；T0是谐振提议分布温度。同一位点E<0拒绝采样会改变其实际动能矩，四组使用完全相同的最终数组，没有分别截断或筛选。\n\n'
    text+=mdtable(['种子/批次','N','未束缚数','T0/uK','Tkin/uK','激发能10/50/90/99%分位/uK'],[[r['id'],r['shots'],r['unbound'],r['proposed_T_uK'],f"{r['kinetic_T_uK']:.3f}",'/'.join(f'{x:.2f}' for x in r['excitation_quantiles_uK'])] for r in preparation])
    text+='\n## 单项消融\n\n以下N=1024、seed26091612、dt=0.0125 us固定；每次仅改变所标机制。它们是敏感性检查，不能把其偶然更小RMSE当作实验机制证据或据此改锁定参数。完整逐曲线差值、早晚偏差见sensitivity_results.csv。\n\n'
    text+=mdtable(['对照','200 RMSE/pp','400 RMSE/pp','600 RMSE/pp','ML RMSE/pp','联合RMSE/pp'],srows)
    (HERE/'VALIDATION.md').write_text(text,encoding='utf-8')

if __name__=='__main__':main()
