"""Render a factual Chinese report from locked, independently validated results."""
from analyze_study import *
from datetime import datetime,timezone
import subprocess,zipfile

def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows)

def main():
    lock=json.loads((HERE/'final_lock.json').read_text(encoding='utf-8'))
    summary=json.loads((HERE/'analysis_summary.json').read_text(encoding='utf-8'))
    c=lock['parameters'];w=c.get('waist_um',1.17)
    rows=[]
    def param(name,old,new,unit,source,adjust,domain,effect,degeneracy):
        rows.append(dict(parameter=name,before=old,manual_200us=new,manual_400us=new,manual_600us=new,ml_400us=new,
            unit=unit,shared=True,source=source,adjustment=adjust,range_basis=domain,effect=effect,degeneracy=degeneracy))
    param('T0 harmonic proposal',15,c['temperature'],'uK','物理工作域内未测标定量','联合调整','4.3–28；微波阶段温度锚点与T<=U/5工作域，非实验置信区间','全部四组初态','热尾/实际制备、RIN、腰')
    param('SLM/AOD nominal waist',1.17,w,'um','SLM论文测量1.17(6)；AOD匹配入瞳假设','共同束腰联合调整','1.11–1.23；采用报告括号不确定范围，未拟合独立AOD腰','阱频、机械激发、噪声耦合','T0、RIN、像差')
    param('AOD RIN PSD',1e-8,c['rin'],'1/Hz','弱白噪声有效模型，尚未装置测量','联合调整','0–1e-8；0–100kHz积分rms<=3.2%弱扰动工作域','能量相关扩散及长时累积损失','谱形、T0、腰')
    param('SLM RIN PSD',1e-8,c.get('slm_rin',c['rin']),'1/Hz','论文132s静态寿命约束；短时安静近似','结构修正，非拟合','固定0；另1e-12敏感性','移除与独立静态寿命矛盾的强加热','噪声是否仅操作期间出现')
    param('Tref',15,'不使用','uK','旧线性化单模参考能量','在物理扩散模型中删除','只在历史/系数修正对照保留15','避免独立加热旋钮','旧固定功率与初温退化')
    param('Heating reference frequency',36007,'不使用固定频率','Hz','原配置显式采用标称AOD径向解析频率','局部力自动确定响应','旧固定功率仍用36007；新模型不再整段套同一个谐振频率','方向、深度、阶段加热差异','RIN与束腰')
    param('Constant total heating',485.631567,'不使用固定功率','kB*uK/s','旧单模参量+指向+固定反冲功率','由局部力扩散替代','历史值保留，修正系数对照的精确值见运行JSON','累计平均能量和逃逸尾','Tref')
    param('Recoil scattering',4.23,'随局部强度计算','1/s','旧SLM+AOD两能级估计；非直接测量','改为局部反冲扩散','新模型采用Cs D2两能级近似，另关闭对照','小幅加性注能','两能级近似、散射各向异性')
    param('SLM depth',140,140,'uK','论文Fig6 Methods直接报告','固定','未报告可信校准误差，未擅自拟合','初阱及放下束缚能','T0、腰')
    param('AOD target depth',280,280,'uK','论文Fig6 Methods直接报告','固定','ML矢量回放终点按论文280；原回放280.5167另存','拾取及等待束缚能','腰、噪声')
    param('Transport distance',2.4,2.4,'um','论文直接报告','固定','ML原始4nm偏置修正来自PDF控制点，非存活拟合','机械激发','波形端点')
    param('Wait',100,100,'us','论文直接报告','固定','100 us；具体SLM关断边沿未报告','暴露时间及残留清除','切换边沿')
    param('SLM wavelength',1061,1061,'nm','论文直接报告','固定','Gaussian Rayleigh range','轴向束缚','光学像差')
    param('AOD wavelength',1055,1055,'nm','论文直接报告','固定','Gaussian Rayleigh range','轴向束缚','光学像差')
    param('Cs mass',132.90545196,132.90545196,'u','原配置Cs-133质量','固定','不移植Rb/Yb','动力学','无')
    param('SLM depth relative sigma',.114,.114,'fraction','论文全阵列实测；子阵列未单独给出','固定','未逐组调散布','位点频率/阱深分布','初温')
    param('AOD depth relative sigma',.02,.02,'fraction','无装置测量假设','不拟合','零散布对照；无可用测量范围','静态位点差异','标定、腰')
    param('AOD waist relative sigma',.02,.02,'fraction','无装置测量假设','不拟合','零散布对照','静态位点差异','标定、腰')
    param('Alignment sigma',50,c.get('alignment_nm',50),'nm','无装置测量假设','不拟合','固定及关闭敏感性','交接误差','位置标定')
    param('Segment alignment jitter sigma',10,c.get('jitter_nm',10),'nm','无带宽依据的原模型假设','不拟合','固定及关闭敏感性；不等同连续PSD','势阱段边沿跳变','指向PSD')
    param('AOD pointing PSD',1e-24,c.get('pointing',1e-24),'m²/Hz','原安静光路敏感性假设','不拟合','固定与关闭对照；只沿运输x轴','加性力扩散','段间抖动')
    param('SLM pointing PSD',1e-24,c.get('slm_pointing',1e-24),'m²/Hz','短时安静近似','结构修正，非拟合','0；另1e-26敏感性','静态加热','静态寿命')
    param('Lensing effective vs',SETTINGS['vs_m_s'],c.get('vs',0),'m/s','旧severity情景；论文SI与器件参数支持额外0.412425情景','模型情景比较，不连续拟合','0在代码中表示关闭，物理上vs趋于无穷；非零须测光束填充率','快速运动的轴向激发','几何、束腰')
    param('Parametric coefficient','pi²/4','pi²','dimensionless','Savard Eq12及本次独立推导','确认修正','PSD仍为单边每Hz，定义未更改','平均加热乘4','PSD幅度')
    param('Pointing coefficient','1/8','1/4','dimensionless','Savard Eq15谱转换及本次独立推导','确认修正','同一PSD定义','平均加热乘2','PSD幅度')
    param('Heating implementation','P=Gamma kB Tref','Gaussian force diffusion','model','Ito力噪声推导，谐振与实际核基准通过','物理结构修正','没有额外运动加热；SLM/AOD噪声分别规定但四组共用','能量分布方差/逃逸尾部','谱平坦性假设')
    param('Preparation','site-bound harmonic','site-bound harmonic','model','已修复的原初始化，保持','不调整','每个实际位点E<0；没有逐组截断/热尾/预筛选','共同初态池','非平衡制备未知')
    param('Observation','E>=0 absorbing','E>=0 absorbing','model','实验末次SLM成像的近似','未拟合阈值','等待末AOD，放下末SLM；另允许再捕获对照','丢失记录时点','成像前等待/关断波形')
    param('Manual timing','48/52','48/52','percent','论文直接报告','固定','二次升深/单段三次移动','不同速度下机械激发','控制校准')
    param('ML interpolation','raster CubicSpline','PDF-marker CubicSpline','model','独立PDF几何，未使用存活率','波形来源修正','未夹紧导数、未单调平滑、位置/深度同钟','回放忠实度','缺原AWG数据')
    csvwrite(HERE/'parameters_before_after.csv',rows)
    parameter_md='# 四组实际参数和物理约束\n\n'+table(['参数','修正前','200 us','400 us','600 us','ML400 us','单位','来源','操作 / 范围'],
        [[r['parameter'],r['before'],r['manual_200us'],r['manual_400us'],r['manual_600us'],r['ml_400us'],r['unit'],r['source'],r['adjustment']+'；'+r['range_basis']] for r in rows])
    parameter_md+='\n操作时长分别为200/400/600/400 us，每30往返总时长为15/27/39/27 ms；这两行是实验允许的差异。主要影响与退化关系详见CSV的effect和degeneracy列。所有数值均为实际传给积分器的共同配置，Tref字段仅在历史对照使用。\n'
    (HERE/'PARAMETERS.md').write_text(parameter_md,encoding='utf-8')

    metrics=summary['metrics'];joint=summary['joint']
    numerical=table(['曲线','基线RMSE/pp','确认修正后/pp','最终独立验证/pp','最终早期偏差/pp','最终后期偏差/pp','最终最大散点偏差/pp'],[
        [label,*[f"{next(r for r in metrics if r['waveform']==wave and r['group']==g['group'])['rmse_pp']:.2f}" for g in joint],
         *[f"{next(r for r in metrics if r['waveform']==wave and r['group']=='Joint model')[key]:+.2f}" for key in ['early_bias_pp','late_bias_pp']],
         f"{next(r for r in metrics if r['waveform']==wave and r['group']=='Joint model')['max_abs_pp']:.2f}"]
        for wave,label in zip(WAVES,LABELS)])
    conver=table(['曲线','dt变化/us','最大S差/pp','该点配对SE/pp','末点个体标签不一致'],[
        [r['waveform'],f"{r['coarse_dt_us']}→{r['fine_dt_us']}",f"{r['max_abs_delta_pp']:.3f}",f"{r['paired_se_pp']:.3f}",f"{100*r['final_label_disagreement']:.2f}%"] for r in summary['convergence']])
    stages=table(['曲线','等待末记录损失','等待开始已正能','其中占比','放下末记录损失','平均噪声功/kB(uK)'],[
        [r['waveform'],r['wait_end_recorded_loss'],r['already_positive_after_switch'],f"{100*r['fraction_positive_before_wait']:.1f}%",r['drop_end_recorded_loss'],f"{r['mean_realized_noise_work_uK']:.2f}"] for r in summary['stages']])
    allcases=json.loads((HERE/'all_cases.json').read_text(encoding='utf-8'))
    fieldnames=sorted(set().union(*(r.keys() for r in allcases)))
    csvwrite(HERE/'scan_history.csv',[{k:r.get(k,'') for k in fieldnames} for r in allcases])
    opt=sorted([r for r in allcases if r['id'].startswith(('refine','selection'))],key=lambda r:r['joint_rmse_pp'])
    opttable=table(['候选','N','T/uK','AOD RIN/Hz^-1','腰/um','vs/m s^-1','联合RMSE/pp'],[
        [r['id'],r['shots'],r['temperature'],r['rin'],r.get('waist_um',1.17),r['vs'],f"{r['joint_rmse_pp']:.3f}"] for r in opt[:12]])
    obs=[]
    for wave in WAVES:
        a=np.load(HERE/'runs/observation_absorbing'/(wave+'.npz'))['alive']
        b=np.load(HERE/'runs/observation_recapture'/(wave+'.npz'))['alive']
        obs.append(dict(waveform=wave,max_difference_pp=float(100*np.max(abs(a.mean(axis=0)-b.mean(axis=0)))),
                        final_difference_pp=float(100*(b[:,-1].mean()-a[:,-1].mean()))))
    write(HERE/'observation_comparison.json',obs)
    ablation=[]
    for case in json.loads((HERE/'diagnostic_queue.json').read_text(encoding='utf-8')):
        m=json.loads((RUNS/case['id']/'summary.json').read_text(encoding='utf-8'))
        ablation.append([case['id'],m['case']['shots'],*[f'{x:.2f}' for x in m['curve_rmse_pp']],f"{m['joint_rmse_pp']:.2f}"])
    ablationtable=table(['机制对照','N','200 RMSE/pp','400 RMSE/pp','600 RMSE/pp','ML RMSE/pp','联合/pp'],ablation)
    report=f'''# Fig. 6d 物理诊断、共同模型修正与联合标定

当前目标分支为 `movement_effects_v1`，开始工作时通过 GitHub 官方API和实际git fetch核实最新提交仍为 **c4775da561be0bf746448c52c5d2a5228937fbdd**。保留旧审计目录和所有历史参考数据，未更改论文散点、归一化或评价区间。本目录共保存{len(allcases)}组完整四曲线计算，详见 `scan_history.csv`。源码精确版本以 `source_manifest.json`、逐运行源文件哈希和源码归档为准。

修正后的生产源码提交为 `{lock['source_commit']}`（本地提交，未推送远端）。研究入口、配置、队列和报告也纳入最终源码归档。

**结论：共同模型显著改善了失配，但尚不能把四条曲线全部解释到实验误差或±2个百分点。** 最终结果是已探索物理工作域中的共同参数模型，不能当作装置参数已被唯一测定，也不宣称找到全局最优或实验损失机制已被证实。旧基线的200 us损失偏大与长时曲线损失不足确实同时存在；它们对应机械操作与随机能量演化两个不同问题，统一放大固定功率不足以解决。

## 最终实际共同模型

四组共同 `T0={c['temperature']} uK`，SLM/AOD标称共同腰 `{w} um`，AOD单边RIN `{c['rin']:.4g} /Hz`；SLM在短时操作窗口采用安静近似。`Tref`从最终扩散模型中删除。透镜情景 `vs={c.get('vs',0)}`（0表示关闭透镜的短途近似，不表示真实声速为零）。原50 nm静态偏差、10 nm段间抖动、2% AOD静态散布没有参与拟合。四组仅保留实验操作时长与波形的差异。

完整四列参数、单位、来源、范围依据、调整理由和退化关系见 [PARAMETERS.md](PARAMETERS.md) 与 [parameters_before_after.csv](parameters_before_after.csv)。实际调用链与配置覆盖顺序见 [CONFIGURATION.md](CONFIGURATION.md)。温度和PSD范围来自明确披露的模型工作域，**不是实验测量的置信区间**。共同腰位于测量允许区间内，但AOD等腰仍依赖入瞳匹配假设。最终腰已触及1.11 um边界，不能据此进一步放宽范围或称其为精确测量。

## 三阶段定量对比

评价使用原图37个实验散点：200 us共7点，其余每条10点。每条先平均点误差，再对四条等权平均；这个规则在选参前冻结。早期固定n<=20，后期n>20。下表残差为模拟减实验，正号表示存活率过高。

{numerical}

联合RMSE：{' → '.join(f"{r['joint_rmse_pp']:.2f} pp" for r in joint)}（基线→确认问题修正→独立验证）。确认修正阶段仍保留固定功率作为逐项对照，不能把它当作最终完整噪声模型。该阶段修正PSD系数、手工波形和ML控制点来源，保持T0=Tref=15 uK。

必须单独指出：**400 us没有实现整体改善**（RMSE约6.43→6.54 pp，差异本身不应过度解读）；其早期平均偏差下降，但后期仍高估存活约10 pp。200/600/ML显著改善，不能掩盖这一点。最终200 us早期存活略低、n=30略高；600 us与ML的早晚存活均偏低，说明共同噪声对不同操作的相对响应及初期损失形状仍不完全符合实验。

最终参数在全部验证前锁定，见 [final_lock.json](final_lock.json)。主结果使用两个新种子26091621/26091622，每种子4096原子，合计8192；机械dt=0.0125 us、噪声h=0.05 us。独立验证结果没有反馈用于重选物理参数。基线和确认修正阶段各使用4096原子。误差条仅来自论文图形，置信口径、每点样本量和相关性未明确，故没有构造精确实验似然或统计“通过”。±2 pp只作工程参照。

![四图总览](fig6d_overview.png)

- [手工200 us](fig6d_manual_200us.png) · [SVG](fig6d_manual_200us.svg)
- [手工400 us](fig6d_manual_400us.png) · [SVG](fig6d_manual_400us.svg)
- [手工600 us](fig6d_manual_600us.png) · [SVG](fig6d_manual_600us.svg)
- [ML400 us](fig6d_ml_400us.png) · [SVG](fig6d_ml_400us.svg)

阴影是逐点Wilson 95% Monte Carlo抽样区间，不含模型误差、时间步误差或真实实验位点相关性。200 us辅助拟合线只到n=31，未外推到60评价；该图横轴限制在32。其余参考拟合线到60。主目标始终是散点，不将密集拟合线视为独立测量。[逐点结果](survival_curves.csv)、[散点残差](point_residuals.csv)、[分曲线指标](metrics.csv)。

## 机制解释与证据强弱

**代码/理论确认。** 单边每Hz PSD下应为Γ=π²ν²S(2ν)、P_x=mω⁴S_x/4，旧系数分别小4倍和2倍。谱定义未改变。旧 `P=ΓkBTref` 是小Γt、单模参考能量的线性化近似；基线整段Γt已不普遍很小，也没有三方向的能量与频率反馈。三维热平衡谐振能量为3kBT，不应把单模kBT误当总能量或再乘错维度因子。

手工图线支持单段 `3s²−2s³`，而旧代码为+J/−J/+J分段S曲线。在15 uK诊断中，单独改回图线支持的函数明显改善200 us，证明这一实现差别足以改变快速操作损失；它同时使慢操作更稳定，反而暴露了后期损失缺失，不能仅靠该修正解释四条实验曲线。

**物理合理但未获装置实验证实。** 新模型对完整时变高斯力施加由PSD推导的随机增量，自动包含能量、方向、深度、阶段与非谐性依赖。平均和方差通过独立谐振子及实际生产核基准；未另加运输加热。原PSD幅值不是测量值，代入新公式直接使用会使长时曲线损失过大。保留同一局部平均漂移、去掉扩散的对照给出不同尾部与存活率，说明平均功率不能代替分布演化；但不证明实验噪声一定是白噪声。

独立静态保持预测排除了“把拟合量级白噪声无条件施加到SLM和AOD两光路”的装置解释：仅SLM强度噪声就使200 ms存活约4.1%，与论文约132 s的无冷却寿命严重不符。更受约束的候选是AOD光路的附加技术噪声、SLM短时近似安静；各光路参数在四组间完全共用。若噪声仅在特定电子操作期间出现，静态寿命约束需要相应重审，当前没有据此再引入自由阶段倍率。

声透镜的vs是有效速度，绝非声速。厂家典型705 m/s和论文SI关系，在填满15 mm孔径、c=7.5条件下给出约0.412425 m/s；该情景与关闭透镜近似都进行了共同参数对照，见扫描记录。未知实际光束填充率使它尚不可唯一确定，不能将旧severity=1值包装为装置标定。

ML位置、深度使用同一时间轴。新的PDF控制点回放对位置图线RMSE约0.000047 um，旧栅格回放约0.00436 um；不夹紧导数、不作额外平滑。旧回放和零速度端点处理另存敏感性对照，均未用于逐曲线补偿。

**剩余猜测。** 非白技术谱、真实AOD束腰/焦平面/非高斯像差、RF功率到阱深的时变标定，以及实际SLM开关边沿可能改变400 us的累积激发。这些目前没有独立证据，未增加自由函数或任意每次技术丢失概率去吸收残差。初态依然是实际位点内的谐振提议拒绝采样，未宣称有限深高斯阱达到严格热平衡。较大样本不能消除这些结构性不确定性。

一个直接的取舍证据是选参池中T=19 uK、腰1.11 um、同N=2048及同seed的RIN对照：从5.75e-9升到6.5e-9/Hz，400 us RMSE由6.86改善至5.39 pp，但600 us由2.71恶化至5.83 pp，ML由3.37恶化至4.37 pp。因此继续统一增大白噪声幅值不能消除目前的相反残差；也不能据此释放400 us独立噪声倍率。这个判断由局部共同参数对照支持，不等于证明所有有物理依据的更完整模型都失败。

特别是n=60，实验400 us与600 us分别约41.8%和57.3%，最终模型却均约52%。这项对操作时长的相对响应是最明确的未解释结构性残差，而不是只缺一个所有曲线共同的温度或加热常数。

详细推导、适用条件和原始文献见 [PHYSICS.md](PHYSICS.md)。每个候选机制的预期可观察结果、实际对照和证据等级见 [MECHANISMS.md](MECHANISMS.md)；[离散初筛图](parameter_screen.png)和[早晚残差原表](diagnostic_residuals.csv)保留未选中的模型。

下面为15 uK、同一seed26091611的机制对照（N见表，初期dt=0.025 us）；系数、手工轨迹、透镜、抖动、扩散均逐项比较。初期较小样本用于辨认较大机制差异，不将百分点以内差异当成显著结论。扩散和平均漂移两行仍是随后被静态寿命排除的“两光路同谱”诊断，不能与最终AOD附加噪声模型混淆。

{ablationtable}

## 能量演化、条件损失与记录时点

[条件损失图](conditional_loss.png)与[逐次CSV](conditional_loss.csv)显示奇偶交接差异；图中不同实验间隔转换为每单程等效几何损失率，不冒充实测逐单次危险率。[能量演化图](energy_evolution.png)与[CSV](energy_evolution.csv)比较每个完整往返回到SLM的幸存者分布；高能原子优先丢失会改变剩余分布，所以幸存者平均能量的下降不表示发生冷却。

{stages}

最后一列是每个原始原子在失去束缚前累计的平均随机噪声功，不能当成三维温升。等待末记录的损失中，已有一部分在拾取结束、关闭SLM后就为正能；因此绝不能用pickup分类为零证明拾取没有激发。新数组保存分段末能量与切换功，而位置和速度连续。

实际比例为约98.6%–99.9%：绝大多数等待末记录的损失在等待开始时已经正能。相同AOD谱产生的平均累计噪声功却分别约15.3/45.0/65.8/45.7 kB·uK，差异自然来自波形、暴露时间、局部力及提前损失；没有逐曲线加热倍率。

在T=18 uK、RIN=7e-9/Hz、腰1.17 um的1024原子对照中，允许正能原子继续传播并可能再捕获，四组S60分别变化{' / '.join(f"{r['final_difference_pp']:.3f}" for r in obs)} pp，说明此观测近似在所测条件下不是主要残差来源；这不是对所有参数的严格定理。最终仍沿用统一E=0束缚边界，未调整阈值。[观测对照](observation_comparison.json)。

## 数值验证与可辨识性

{conver}

上表数值收敛比较覆盖0–60全部检查点；这不是对200 us论文拟合线作区间外评价。
600 us第一次减半的最大差约2.1 pp，因此追加dt=0.003125 us的针对性复核，保持同一初态与噪声增量；这个检查没有用于重选温度、PSD、束腰或主报告步长。

追加复核的最大差降到0.98 pp。噪声h=0.05→0.025 us的四曲线最大分布差为约1.12/1.05/1.44/0.71 pp，末点差均不超过0.37 pp。两个新种子的合并联合RMSE为4.26 pp，按原子重抽样的95% MC区间约4.05–4.51 pp；该区间完全不包括装置参数和实验不确定性。数值验证支持总体趋势与主要残差的判断，**不支持逐原子逃逸轨迹已收敛或亚百分点数值精度**。

机械dt减半时保持同一个噪声h和每原子随机序列。存活率接近不表示个体逃逸轨迹已严格收敛，末点标签不一致比例如上；较大的个体差异应作为混沌逃逸边界的数值敏感性报告。另有新种子的噪声h减半、两独立种子差异、保留时间相关性的MC重抽样以及逐项消融，数值表和解释见 [VALIDATION.md](VALIDATION.md)；不按哪一步长更贴实验选步长。MC区间针对独立模拟噪声实现，不能替代真实实验同光路原子间共享噪声的误差模型。

`noise_benchmarks.json`的参量/指向平均增长与能量方差在三档步长均符合推导；`gaussian_benchmark.json`进一步验证实际全高斯力核的径向/轴向增长率和注能账本。`tests_rng.xml`验证噪声在不同线程数下可重复。原实际位点束缚初始化和连续传播测试均通过。整套非slow回归245项直接通过，另两项安装环境问题在独立venv中修复后通过；新物理/线程检查另外验证。未对本轮无关的10项slow测试宣称执行。

以下较大样本候选显示温度、RIN与共同束腰的退化，差异不能作为精确参数置信区间。选择依据是固定共同目标；没有四条曲线分别选温再拼接。

{opttable}

全部已探索范围：T=4.3–28 uK，RIN=0–1e-8/Hz，共同腰1.11–1.23 um；透镜关闭、旧severity、器件几何情景；原/矢量ML、端点边界、静态误差与抖动、平均漂移/随机扩散、吸收/再捕获均有有限对照。范围边缘和未测噪声谱限制最终解释，不能以有限扫描证明“任何共同模型均不可能”。

## 代码修改与复现

- `noise_heating.py`：修正两个PSD系数与说明，去掉纯音Floquet指数能锁定随机谱系数的错误说法。
- `paper_waveforms.py`、`protocol.py`、`cli.py`：单段三次手工轨迹，保留显式legacy选项；波形与时序不随拟合曲线改变。
- `survival_engine.py`：反向操作保留波形子类，避免新单段三次在drop-off路径被错误退回旧S曲线。
- `noisy_survival.py`：完整高斯力扩散，分光路但全曲线共用谱；分开机械/噪声步长，增加切换功、分段能量与再捕获诊断。
- 本目录 `run_study.py`：显式共同配置、固定散点评分、逐原子NPZ/逐运行完整配置/源码哈希；`waveform_audit.py`从PDF独立读取控制点。

新目录下的历史输出不覆盖。每组 `.npz`含原始alive、状态、能量、实际位点与初态；`.json`含所有实际case参数、底层配置、种子、N、dt/h、源码哈希和指标。参考PDF和参考CSV哈希记录于 `source_manifest.json`。`source_diffusion_v1.zip/v2.zip`保留早期模型，`source_final.zip`冻结最终实现。

当前目录已准备独立venv，下面命令在本机可直接使用。迁移环境时安装 `pip install -e "simulation/level2_joint_transfer[continuous,dev]"` 及 `pdfplumber`，准确版本见 `requirements_exact.txt`；可用脚本的 `-Python` 指定解释器。运行示例：

```powershell
.\\reports\\movement_effects_v1_calibration_20260916\\reproduce.ps1 -Phase final -NewRunDirectory D:\\work\\fig6d_reproduction_new
.\\reports\\movement_effects_v1_calibration_20260916\\reproduce.ps1 -Phase benchmarks
.\\reports\\movement_effects_v1_calibration_20260916\\reproduce.ps1 -Phase analyze
```

精确队列见 `diagnostic_queue.json`、`scan_queue.json`、`aod_queue.json`、`waist_queue.json`、`lens_queue.json`、`opt3_queue.json`、`refine*_queue.json`、`selection_queue.json`、`final_queue.json`、`validation_checks_queue.json`、`sensitivity_queue.json`。已存在同名同配置结果会复用，要求实际重算时必须用新的 `-NewRunDirectory`。分析当前归档用 `-Phase analyze`。不要直接运行仓库旧默认连续协议并把它当成本报告配置。

## 最需要的独立实验输入

1. 转移前、140 uK SLM中的同一制备温度及能量分布；最好同步做释放再捕获与旁带测温，以区分初温和初态非平衡。
2. AOD和SLM各自原子平面的RIN与指向PSD，覆盖轴向/径向共振及两倍频，同时测两光路相关性；单独AOD静态保持随阱深/时间的损失可直接检验拟合噪声候选。
3. 实际14个ML控制量、AWG时钟/端点约束、RF到位置/深度校准、SLM关断与重开边沿、成像前保持时间。
4. 195个参与位点的AOD腰/深度/焦平面分布与光束在AOD内的直径；该信息可决定透镜情景与400 us剩余机械激发。

这些测量能区分候选机制；继续增加未经约束的拟合自由度不能替代它们。
'''
    (HERE/'REPORT.md').write_text(report,encoding='utf-8')

if __name__=='__main__':main()
