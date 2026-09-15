import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
import {execFileSync} from 'node:child_process';
const here=path.dirname(fileURLToPath(import.meta.url)),root=path.resolve(here,'../..');
const runtime='C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=path.join(runtime,'node/node_modules');
process.env.RUNTIME_NODE=path.join(runtime,'node/bin/node.exe');
process.env.RUNTIME_PYTHON=path.join(runtime,'python/python.exe');
const skill='C:/Users/Administrator/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations';
const require=createRequire(path.join(runtime,'node/node_modules/_resolver.cjs'));
const {Presentation,PresentationFile,FileBlob}=await import(pathToFileURL(require.resolve('@oai/artifact-tool')).href);
const {applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const ePath=path.join(here,'output/Level5_vector_reference_v2.evidence.json');
const e=JSON.parse(await fs.readFile(ePath,'utf8'));
const build=path.join(here,'.build'),outDir=path.join(here,'output');
await fs.mkdir(build,{recursive:true});await fs.mkdir(outDir,{recursive:true});
const rev=process.env.REPORT_REV||'vector-v2-r1';
const output=process.env.REPORT_PPTX||path.join(outDir,'Level5_MonteCarlo_论文对照与全流程_2026-09-15_参考数据修正版.pptx');
const p=Presentation.create({slideSize:{width:1280,height:720}});
const C={bg:'#FFFFFF',navy:'#17324D',text:'#344054',blue:'#1467B3',teal:'#007F78',red:'#BB3E35',gray:'#687787',light:'#F1F5F8',line:'#D9E2EA',white:'#FFFFFF'};
const FONT='Microsoft YaHei',titles=[],textBySlide=[],tableOwners=[],chartOwners=[];
const paper=path.join(root,'2403.12021v4.pdf');
const codeDir=path.join(root,'simulation/level2_joint_transfer');
const labels=['手工 200 μs','手工 400 μs','手工 600 μs','论文 ML 400 μs'];
function tx(s,t,x,y,w,h,size=25,color=C.text,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=t;sh.text.style={typeface:FONT,fontSize:size,color,bold,autoFit:'none'};
 textBySlide[p.slides.items.indexOf(s)]?.push(t);return sh;
}
function slide(title,sub,ref=''){
 const s=p.slides.add();s.background.fill=C.bg;titles.push(title);textBySlide.push([]);
 tx(s,title,58,32,1164,62,40,C.navy,true);
 tx(s,sub,60,105,1160,47,21,C.gray);
 tx(s,'Level 5 论文对照与全流程',60,682,680,20,13,C.gray);
 tx(s,String(p.slides.items.length).padStart(2,'0'),1160,680,60,24,15,C.gray);
 s.speakerNotes.textFrame.setText(`Primary source: ${paper}\n${ref}\nCorrected evidence: ${ePath}\nReference: ${codeDir}/reference/fig6d_vector_v2/fig6d_vectors.json\nFit paths, 37 marker centers and original error-bar endpoints are recovered independently from PDF vector objects. No smoothing or refitting. Error-bar geometry does not establish a confidence convention. Dense fit samples are not independent measurements.\nArchived 2026-09-14 simulations are re-compared, not rerun. The previously frozen temperatures are retained.\nSimulation: ${codeDir}/src/continuous_transfer/{protocol,engine,cli}.py\nConfig: ${codeDir}/configs/continuous_level2c_345.yaml\nBaseline physical assumptions: ${codeDir}/configs/level2c_pickup_survival.yaml\nLegacy result status: ${codeDir}/outputs/level5_finite_pulse_irb/demo/RESULT_STATUS.md\nScope: this deck is a standalone scientific comparison. Continuous mechanical + finite-pulse results and legacy independent RB evidence have different validation status.`);
 return s;
}
function table(s,rows,x=60,y=178,w=1160,h=365,widths=[1.7,4.7],size=24){
 const t=s.tables.add({rows:rows.length,columns:rows[0].length,left:x,top:y,width:w,height:h,columnTracks:widths.map(value=>({mode:'fr',value})),values:rows});
 t.borders.assign({style:'solid',fill:C.line,width:1});
 for(let r=0;r<rows.length;r++)for(let c=0;c<rows[0].length;c++){
  const cell=t.getCell(r,c);cell.fill=r===0?C.navy:r%2?C.white:C.light;
  cell.text.style={typeface:FONT,fontSize:size,bold:r===0,color:r===0?C.white:C.text,verticalAlignment:'middle',autoFit:'none'};
 }
 tableOwners.push(p.slides.items.indexOf(s)+1);textBySlide[p.slides.items.indexOf(s)].push(rows.map(r=>r.join(' / ')).join('\n'));return t;
}
function note(s,t,y=606,color=C.blue,size=24){tx(s,t,60,y,1160,64,size,color,true);}
async function img(s,file,x,y,w,h){s.images.add({blob:new Uint8Array(await fs.readFile(path.join(build,'paper',file))),contentType:'image/png',alt:`论文原图 ${file}`,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
function chart(s,categories,series,x,y,w,h,ytitle='存活概率 S(n)'){
 const c=s.charts.add('line',{position:{left:x,top:y,width:w,height:h},categories,series:series.map(z=>({...z,values:z.values.map(v=>Number(v.toPrecision(12)))})),lineOptions:{smooth:false},hasLegend:true,chartFill:C.white,plotAreaFill:C.white});
 applyPresentationChartFont(c,{fontFamily:FONT});
 c.legend={position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:16,fill:C.gray}};
 c.xAxis={visible:true,textStyle:{typeface:FONT,fontSize:13,fill:C.gray},line:{fill:C.line,width:1},majorGridlines:null};
 c.yAxis={visible:true,title:ytitle,min:0,max:1,numberFormatCode:'0.0',textStyle:{typeface:FONT,fontSize:16,fill:C.gray},majorGridlines:{fill:C.line,width:1}};
 applyPresentationChartFont(c,{fontFamily:FONT});
 chartOwners.push(p.slides.items.indexOf(s)+1);return c;
}
const f=x=>Number(x).toFixed(4),pp=x=>(100*x).toFixed(2),qcp=r=>r.quantum.checkpoints.at(-1);
function comparisonChart(s,row,x,y,w,h){
 const ref=e.reference_vectors.series[row.waveform], end=row.comparison_paper_fit.n.at(-1);
 const simn=Array.from({length:end+1},(_,n)=>n);
 const c=s.charts.add('scatter',{position:{left:x,top:y,width:w,height:h},
  scatterOptions:{style:'lineWithMarkers'},hasLegend:true,chartFill:C.white,plotAreaFill:C.white,
  series:[
   {name:'论文拟合线',xValues:ref.display_n,values:ref.display_fit,marker:{symbol:'none'},line:{fill:C.navy,width:2.5}},
   {name:'论文散点',xValues:ref.markers.map(m=>m.n),values:ref.markers.map(m=>m.survival),fill:C.navy,marker:{symbol:'diamond',size:4},line:{fill:'none',width:0}},
   {name:'连续仿真',xValues:simn,values:simn.map(n=>row.survival[n]),marker:{symbol:'none'},line:{fill:C.red,width:2.5}}
  ].map(r=>({...r,values:r.values.map(v=>Number(v.toPrecision(12))),xValues:r.xValues.map(v=>Number(v.toPrecision(12)))}))});
 c.legend={position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:15,fill:C.gray}};
 c.xAxis={visible:true,min:0,max:end,majorUnit:end<40?5:10,numberFormatCode:'0',textStyle:{typeface:FONT,fontSize:14,fill:C.gray},line:{fill:C.line,width:1},majorGridlines:null};
 c.yAxis={visible:true,title:'存活概率 S(n)',min:0,max:1.02,majorUnit:.2,numberFormatCode:'0.0',textStyle:{typeface:FONT,fontSize:15,fill:C.gray},majorGridlines:{fill:C.line,width:1}};
 applyPresentationChartFont(c,{fontFamily:FONT});chartOwners.push(p.slides.items.indexOf(s)+1);return c;
}
// 01
{
 const s=slide('Level 5 Monte Carlo','论文对照、差距原因与完整流程');
 tx(s,'与论文匹配到哪里',60,246,1160,84,60,C.navy,true);
 tx(s,'从原子热运动到内部态演化，再到实验可比观测量',64,365,1150,67,31,C.text);
 tx(s,'依据 2026-09-14 连续结果，2026-09-15 修正论文参考并重新对照\n论文：A tweezer array with 6100 highly coherent atomic qubits（v4）',64,540,1150,85,22,C.gray);
}
// 02
{
 const s=slide('与论文的重合程度','按实验过程、协议条件和观测量分别判断，当前没有可定义的“总体复现百分比”','Fig. 4–6, pp. 5–7; Extended Data Fig. 10, p. 19.');
 table(s,[['论文内容','过程覆盖','定量匹配状态'],['Fig. 6 短转移存活','较高：拾取、2.4 μm 分离、等待、放回','可直接检验的四条曲线均未通过'],['Fig. 4 单比特控制与相干性','部分：两能级、有限脉冲、已有 RB 框架','T₂*、T₂、Clifford 保真度未同条件验收'],['Fig. 5 AOD 长运输','部分：三维运动、透镜效应、DD','最新连续结果未复现 610 μm IRB 协议'],['Ext. Data Fig. 10 组合操作','部分：拾取、长距离去回、放回','375 μm 重合，完整时序与背景阱不同']],60,174,1160,385,[2.1,3.5,3.6],23);
 note(s,'当前证据支持“部分物理过程已覆盖”，尚不支持“Level 5 已定量复现论文”。',591,C.red,25);
}
// 03
{
 const s=slide('Fig. 6：短转移过程的主要重合','论文将一次 pick-up & split 或一次 merge & drop-off 计为一次单向转移','Fig. 6a–d, p. 7; Methods, Large-scale coherent atom transfer.');
 await img(s,'fig6abc.png',60,182,570,356);
 tx(s,'主要条件已沿用',674,180,520,43,29,C.blue,true);
 tx(s,'SLM 140 μK，AOD 280 μK\n分离 2.4 μm，等待 100 μs\n等待期间关闭 SLM\n手工 200 / 400 / 600 μs\n回放论文 ML 400 μs 波形',674,239,535,256,26);
 tx(s,'仍有差异',674,518,500,40,27,C.red,true);
 tx(s,'实际波形响应、转移温度和 AOD 标定尚未确定。ML 回放也未复现论文的实验优化过程。',674,569,535,94,23);
 tx(s,'论文原图 Fig. 6a–c',60,562,560,30,18,C.gray);
 tx(s,'早期单次转移保真度 99.81(3)%\n对应约前 12 次单向转移',60,607,570,65,24,C.navy,true);
}
// 04
{
 const s=slide('移动距离与“60 次”的计数','当前配置的 n=60 是 60 次单向短转移，即 30 轮拾取与放回','Fig. 5, p. 6; Fig. 6, p. 7; Extended Data Fig. 10, p. 19; continuous protocol settings.');
 table(s,[['项目','一次 / 一轮','当前 n=60 累计'],['论文 Fig. 5 AOD-only','对角单向 610 μm，1.6 ms','与短转移 n 的操作定义不同'],['论文 Fig. 6 短转移','单向 2.4 μm，往返 4.8 μm','名义路程 144 μm'],['当前短转移对照','400 + 100 + 400 μs / 轮','30 轮，27 ms，144 μm'],['当前组合运输','每轮另加 375 μm 去程与回程','30 轮，115.62 ms，22.644 mm']],60,174,1160,365,[2.4,3.9,3.3],23);
 note(s,'组合运输包含长距离去程与回程。若确指 60 轮往返，路程为 45.288 mm。',567,C.blue,25);
 tx(s,'路程按手工波形名义位移计算。ML 回放存在细小折返。当前组合时序并非论文组合协议的逐项复现。',60,635,1160,37,19,C.gray);
}
// 05–06
for(const pair of [[0,1],[2,3]]){
 const s=slide(pair[0]===0?'Fig. 6：快速手工波形仍有偏差':'Fig. 6：慢波形与 ML 末段存活仍偏高','沿用此前冻结参数与 4096 样本结果；论文拟合线、实验散点和误差条分别读取','Fig. 6d top panel; independent vector fit paths, marker centers and error bars; evidence.selected; dt = 0.0125 μs, seed = 23003.');
 pair.forEach((idx,col)=>{
  const r=e.selected[idx],cp=r.comparison_paper_fit,x=60+col*598;
  tx(s,`${labels[idx]}，T = ${r.temperature_uK} μK`,x,166,562,42,27,C.navy,true);
  comparisonChart(s,r,x,224,560,327);
  tx(s,`n=${cp.n.at(-1)}：论文拟合 ${f(cp.reference.at(-1))}，仿真 ${f(cp.simulated.at(-1))}\n对拟合线的最大偏差 ${pp(cp.max_abs_deviation)} 个百分点`,x,561,560,72,22,C.text);
 });
 tx(s,'横轴为单向转移数。比较区间仍为 n≥5；四条均未通过 ±2 个百分点容差，该容差与论文误差条分别处理。',60,641,1160,35,18,C.red,true);
}
// 07
{
 const s=slide('一组共同参数仍不能同时解释四条曲线','共同扫描较优格点 T = Tref = 15 μK；下表为独立 4096 样本、精细步长复核','evidence.temperature_scan, precision[15]; equal weight of four curves. Tref changes mean heating together with initial T.');
 table(s,[['波形 / 比较末点','论文拟合 S','仿真 S','对拟合线最大偏差'],...e.precision['15'].map((r,i)=>[`${labels[i]} / n=${r.comparison_paper_fit.n.at(-1)}`,f(r.comparison_paper_fit.reference.at(-1)),f(r.comparison_paper_fit.simulated.at(-1)),`${pp(r.comparison_paper_fit.max_abs_deviation)} 个百分点`])],60,184,1160,337,[3.4,1.5,1.3,2.9],24);
 note(s,'200 μs 损失过强，600 μs 与 ML 损失不足，调同一个温度无法兼顾。',555,C.red,25);
 tx(s,'Tref 随 T 改变，所以扫描同时改变初始热分布和平均加热强度。该最优格点不能解释为测得的转移温度。',60,625,1160,46,22);
}
// 08
{
 const s=slide('增加三维运动与透镜效应后的变化','同为 T = 5 μK 的诊断对照；512 样本，最大步长 0.05 μs','evidence.base[5]; paper Methods says cylindrical lensing is not detrimental at the short transfer speeds.');
 table(s,[['波形末点 S','论文拟合','原 1D','3D 无透镜','3D 含透镜'],...labels.map((label,i)=>[label,f(e.base['5'].matched[i].comparison_paper_fit.reference.at(-1)),...['baseline_1d','matched_lensing_off','matched'].map(a=>f(e.base['5'][a][i].comparison_paper_fit.simulated.at(-1)))])],60,179,1160,321,[2.3,1.3,1.3,1.6,1.6],24);
 tx(s,'模型中的透镜效应主要压低最快的 200 μs 曲线，慢波形仍缺少衰减。',60,537,1160,51,26,C.blue,true);
 tx(s,'论文称其短转移速度下柱面透镜不构成明显限制。当前透镜参数仅是未标定场景，不能据此认定实验短转移损失来自透镜效应。',60,600,1160,68,23,C.red);
}
// 09
{
 const s=slide('Fig. 4：相干性和单比特门尚未完成同条件复现','数值常数进入模型，并不等于对应的实验衰减曲线已得到验证','Fig. 4, p. 5; microwave/RB Methods; configs/level5_finite_pulse_irb.yaml.');
 table(s,[['论文观测','论文结果 / 条件','当前证据'],['Rabi 振荡','Ω/2π = 24.611(1) kHz','已采用频率，未验收空间平均衰减'],['Ramsey 相干','系综 T₂* = 14.0(1) ms','只作参考锚点，未按同协议拟合'],['XY16 相干','T₂ = 12.6(1) s，π 间隔 12.5 ms','存在 DD 模型，未验证这条长时曲线'],['随机 Clifford 门','平均保真度 99.9834(2)%','已有独立 RB 框架，暂无同条件有效结果']],60,176,1160,359,[2.1,3.7,4.1],23);
 note(s,'主要缺口是实际噪声、空间分布与测量协议的约束。',569,C.blue,27);
 tx(s,'论文微波阶段约 4.3 μK 的温度，不能直接代替拾取/运输时的温度。',60,632,1160,39,23);
}
// 10
{
 const s=slide('Fig. 5：610 μm 长距离相干运输','论文分别测量存活、IRB 返回概率，再提取早期单次瞬时保真度','Fig. 5a–d, p. 6; Randomized benchmarking of coherent transport Methods.');
 await img(s,'fig5.png',60,194,566,387);
 tx(s,'论文的比较目标',678,176,525,43,29,C.blue,true);
 tx(s,'10 个 AOD 阱，深度 280 μK\n610 μm 对角移动，1.6 ms\n每次移动配合 XY4\n固定 80 个 Clifford，改变移动数\n早期单次保真度 99.953(2)%',678,232,538,229,25);
 tx(s,'当前重合与差异',678,478,538,42,28,C.red,true);
 tx(s,'已覆盖运动、透镜与有限脉冲。最新连续结果采用 375 μm 组合段，尚无相同 Clifford 与测量协议的 610 μm 结果。',678,535,538,132,24);
 tx(s,'论文原图 Fig. 5',60,603,566,34,18,C.gray);
}
// 11
{
 const s=slide('组合运输：375 μm 重合，完整协议仍不同','论文 Extended Data Fig. 10e–f 的一次单向操作含门、拾取/分离和长距离移动','Extended Data Fig. 10e–f, p. 19; combined transfer and long moves Methods, pp. 28–29.');
 await img(s,'ext10e.png',60,171,566,351);
 tx(s,'论文：每轮约 5.508 ms',676,170,545,41,28,C.blue,true);
 tx(s,'拾取 850 μs + 分离 400 μs\n长移动 1.45 ms，去回各一次\n两端各有 54 μs 单比特门\nSLM 180 降至 60 μK\n47 个 AOD，保留静态阵列通道',676,222,544,224,24);
 tx(s,'当前：每轮 3.854 ms',676,471,544,43,28,C.red,true);
 tx(s,'短转移每程合计 400 μs\n含 100 μs 等待和远端 54 μs 保持\n长移动期间 SLM 关闭，无静态阵列通道',676,525,544,118,23);
 tx(s,'论文早期单次组合保真度 99.87(1)%\n当前结果尚不能与该数值直接相减。',60,566,568,88,24,C.navy,true);
}
// 12
{
 const s=slide('当前 Level 5 内部态结果的含义','手工 400 μs、T = 5 μK；相干对比度按存活原子统计，相对于理想脉冲演化','evidence.base[5].matched[1] and extended[1]; quantum_summary. Conditional coherence C = 2|mean rho01| of ideal-relative + state.');
 const r=e.base['5'].matched[1],z=e.base['5'].extended[1],ns=r.quantum.checkpoints.map(a=>a.n);
 chart(s,ns.map(String),[{name:'短转移 无 DD',values:r.quantum.checkpoints.map(a=>a.contrast_none),line:{fill:C.gray,width:2}},{name:'短转移 XY4',values:r.quantum.checkpoints.map(a=>a.contrast_dd),line:{fill:C.blue,width:3}}],60,181,718,387,'条件相干对比度 C(n)');
 table(s,[['n=60 / 30 轮','短转移','组合运输'],['存活 S',f(r.survival[60]),f(z.survival[60])],['无 DD 的 C',f(qcp(r).contrast_none),f(qcp(z).contrast_none)],['有 DD 的 C',f(qcp(r).contrast_dd),f(qcp(z).contrast_dd)]],814,204,406,304,[1.7,1.1,1.1],22);
 tx(s,'本参数场景中，DD 提升短转移的累计相干对比度，对组合运输的改善不明显。',60,587,1160,48,25,C.blue,true);
 tx(s,'这些是连续模型的累计条件量，尚未转换成论文的 IRB 返回曲线或单次瞬时保真度。',60,641,1160,34,21,C.red);
}
// 13
{
 const s=slide('论文指标与当前输出的对应关系','存活、条件相干性、累计通道保真度和单次操作保真度的分母与含义不同','Fig. 5d/6d and corresponding Methods; quantum_summary definition and feedback fields.');
 table(s,[['指标','分母 / 含义','当前对照状态'],['存活 S(n)','存活数 / 初始原子数','Fig. 6d 顶图可直接对照'],['相干对比度 C(n)','仅存活原子，某个叠加态的相干幅度','用于诊断，不能替代通道保真度'],['累计条件 Favg(n)','存活通道相对理想脉冲的平均保真度','涵盖此前所有操作，尚非 IRB 结果'],['IRB 返回概率 R(n)','完整随机门序列后回到初态的概率','最新连续流程尚未生成'],['单次瞬时保真度 Fₙ','联合 S、R 的拟合与损失处理得到','尚无同协议的可比仿真值']],60,175,1160,420,[1.9,4.2,3.7],23);
 tx(s,'例如 30 轮后的累计 Favg，不能与论文 99.953% 的早期单次长移动保真度比较。',60,625,1160,48,23,C.red,true);
}
// 14: explicitly requested process diagram, native editable objects.
{
 const s=slide('Monte Carlo 全流程与状态传递','每个样本表示一个原子的一整段实验历史，运动和内部态共同保留到下一轮','continuous_transfer engine/protocol/cli; no state-dependent force; no resampling across rounds.');
 const node=(label,x,y,w,h)=>{const n=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:C.light,line:{fill:C.blue,width:1.4}});n.text=label;n.text.style={typeface:FONT,fontSize:24,color:C.navy,bold:true,alignment:'center',verticalAlignment:'middle',autoFit:'none'};return n;};
 const a=node('① 初始系综\nr₀、v₀、阱参数',60,198,245,112),b=node('② 控制时序\n位置、阱深、脉冲',365,198,245,112),c=node('③ 三维动力学\nr(t)、v(t)、能量',670,198,245,112),d=node('④ 内部态演化\n光移、相位、U(t)',670,446,245,112),g=node('⑤ 检查点与下一轮\n保留状态，失阱终止',975,198,245,112),fnode=node('⑥ 系综观测\nS、C、累计 Favg',975,446,245,112);
 const conn=(x,y,from='right',to='left')=>s.shapes.connect(x,y,{kind:'straight',fromSide:from,toSide:to,line:{fill:C.blue,width:2},tail:{type:'arrow',width:'med',length:'med'}});
 conn(a,b);conn(b,c);conn(c,g);conn(c,d,'bottom','top');conn(d,fnode);conn(g,fnode,'bottom','top');
 tx(s,'轨迹势能\nV(t)',636,350,144,67,21,C.gray);
 tx(s,'存活标记\n运动状态',982,355,106,60,19,C.gray);
 tx(s,'持续传递的数据',60,412,538,40,29,C.blue,true);
 tx(s,'原子编号、位置 r、速度 v、内部态 U\n存活标记、累计时间、累计加热\n下一段 / 下一轮沿用上一步的输出',60,472,550,123,25);
 tx(s,'控制表驱动物理演化。内部态目前接受运动产生的光移，但不反向改变机械运动。',60,631,1160,43,23,C.red);
}
// 15
{
 const s=slide('过程 1：生成热初态与原子间差异','Monte Carlo 随机性来自初始热分布、静态参数散布和逐段对准抖动','thermal_sampling_3d.py; draw_disorder; case_inputs. Sampling is a harmonic proposal with exact E<0 rejection, not the full finite-depth canonical ensemble.');
 table(s,[['环节','参数与状态'],['输入','Cs-133 质量 m，温度 T（μK），SLM 深度（μK）、束腰（μm）、波长（nm），样本数 N 与随机种子'],['处理','抽样三维位置与速度；用完整高斯势的总能量 E<0 接受束缚初态。为每个原子抽取阱深、束腰和对准偏差'],['输出','原子编号 i，r₀=(x,y,z)（m），v₀（m/s），初始能量 E₀（J），各原子的静态参数因子'],['传给下一过程','控制时序对同一批原子生效。r₀、v₀ 和静态差异进入第一段动力学']],60,175,1160,374,[1.3,7.4],24);
 note(s,'热初态采用简谐近似抽样；转移时的真实温度与高能尾部分布仍待约束。',596,C.blue,25);
}
// 16
{
 const s=slide('过程 2：把实验步骤变成随时间变化的控制量','短转移对照与组合运输采用独立的时序定义','continuous_transfer/protocol.py, build_protocol; manual pickup_spec and digitized ML CubicSpline.');
 table(s,[['环节','参数与状态'],['输入','波形类型，单程时长（200 / 400 / 600 μs），短位移 2.4 μm，等待 100 μs；组合额外输入 375 μm、1.45 ms 与远端保持 54 μs'],['处理','顺序安排拾取、等待、可选长距离去回、放回；安排 SLM 开关、AOD 升降深和有限宽度脉冲'],['输出','控制表 cₓ(t)、cᵧ(t)（m），dc/dt（m/s），AOD 深度 D(t)（J），SLM 开关，Ω(t)（rad/s）与脉冲相位 φ(t)（rad）'],['传给下一过程','位置、深度和 SLM 开关决定势场；AOD 速度决定透镜修正；Ω、φ 同步驱动内部态']],60,175,1160,403,[1.3,7.4],23);
 tx(s,'手工波形先升深再移动。ML 同步升深与移动，当前使用论文曲线的数字化回放。',60,613,1160,49,25,C.blue,true);
}
// 17
{
 const s=slide('过程 3：运动、平均加热与瞬时光阱','位置和速度贯穿拾取、移动、等待与放回，不在阶段交界处重新抽样','engine.total_field/integrate; level2_joint_transfer/noise_heating.py; heating_policy.');
 table(s,[['环节','参数与状态'],['输入','当前 r（m）、v（m/s），两束光的深度、束腰与中心，AOD 扫描速度，透镜参数 vₛ（m/s），总平均加热功率 P（J/s）'],['处理','叠加 SLM 与 AOD 的三维高斯势；势能梯度产生力。每个时间步推进 r、v，并向总动能加入 P·dt'],['输出','r(t)、v(t)、动能与势能（J），累计注入能量（J），沿轨迹的总势能 V(t)'],['传给下一过程','r、v 和能量进入失阱判断；V(t) 同时进入光移与内部态传播；下一步继续沿用运动状态']],60,175,1160,389,[1.3,7.4],23);
 tx(s,'当前加热表示平均功率，未生成逐时刻随机噪声的完整谱与稀有强踢。',60,605,1160,52,25,C.red,true);
}
// Mean-heating processes, part of stage 3.
{
 const s=slide('平均加热的三个过程','各过程先给出平均功率，再共同进入三维动力学的每个时间步','level2_joint_transfer/noise_heating.py; level2c_pickup_survival.yaml; current continuous metadata.heating_power_uK_per_s. All noise amplitudes are scenario assumptions.');
 table(s,[['过程','输入参数','输出与衔接'],['光子反冲','散射率 Γsc（s⁻¹）\n波长 λ（nm）与原子质量 m','散射带来的平均加热功率\nPrec（J/s）'],['强度噪声 / 参量加热','强度噪声谱 S_RIN（1/Hz）\n参考阱频 ν₀（Hz）、Tref（μK）','固定参考能量下的平均功率\nPpar（J/s）'],['指向噪声','位置噪声谱 Sx（m²/Hz）\n参考阱频 ν₀ 与质量 m','阱中心抖动对应的平均功率\nPpoint（J/s）']],60,175,1160,350,[2.1,4.4,3.7],24);
 note(s,'总功率 P = Prec + Ppar + Ppoint，每步注入动能 P·dt。',564,C.blue,27);
 tx(s,'三项目前都是恒定平均功率。转移中实际噪声谱、随机强踢和 Raman 内部态跃迁尚未在此流程中展开。',60,631,1160,42,22,C.red);
}
// 18
{
 const s=slide('过程 4：轨迹光移与有限脉冲内部态演化','两能级对应 Cs 钟态；同一原子的运动轨迹决定其经历的失谐历史','engine.spin_step; quantum_summary and ideal_checkpoint_unitaries. Nominal δ = η V/ħ, sign according to potential convention.');
 table(s,[['环节','参数与状态'],['输入','沿轨迹的势能 V(t)（J），差分光移系数 η=1.3×10⁻⁴，Ω/2π=24.611 kHz，脉冲相位 φ(t)，前一时刻的 U'],['处理','由 Δ(t)=ηV(t)/ℏ 得到失谐（rad/s），积分两能级演化。裸 π 脉冲约 20.3 μs，脉冲期间运动继续'],['输出','累计相位（rad）与每原子的内部态演化矩阵 U(t)；据此可求任意输入量子态的输出'],['传给下一过程','U(t) 与原子编号、r、v 一起保留。只有存活原子的内部态进入条件相干性和通道统计']],60,175,1160,391,[1.3,7.4],23);
 tx(s,'当前使用 XY4。论文的随机 Clifford 序列与随之变换的 DD 旋转轴，尚未接入这条连续流程。',60,604,1160,62,24,C.red,true);
}
// 19
{
 const s=slide('过程 5：失阱判断与跨轮连续传递','论文重复操作会积累热运动和相位，连续模型保留这两种记忆','engine.integrate judge fields; protocol.segments: pickup_wait AOD, long segment AOD extra checkpoints, dropoff SLM.');
 table(s,[['环节','参数与状态'],['输入','当前 r、v、U，接收阱的深度与中心，原子编号，以及当前阶段和轮次'],['检查点','等待结束检查 AOD 束缚，放回结束检查 SLM 束缚。组合运输额外在去程、远端保持和回程末端检查 AOD'],['判据与输出','接收阱中的 E = mv²/2 + V < 0 才存活。输出 alive、失阱阶段、失阱轮次和检查点处的完整状态'],['传给下一轮','存活原子沿用 r、v、U、时间与累计加热。失阱原子终止，不补入新原子，也不重置温度']],60,175,1160,387,[1.3,7.4],23);
 tx(s,'损失在指定检查点判定。旧独立 RB 路径逐次抽取单轮轨迹，其跨轮记忆与当前连续模型不同。',60,605,1160,62,24,C.blue,true);
}
// 20
{
 const s=slide('过程 6：系综统计与实验读出之间的接口','当前连续计算终止于 S、条件 C 和累计 Favg；复现论文 IRB 还需要完整读出协议','quantum_summary; paper Randomized benchmarking of coherent transport and Large-scale coherent atom transfer Methods.');
 table(s,[['环节','输入与输出'],['当前已完成的汇总','输入所有原子的 alive、U 与理想脉冲目标。输出 S(n)、条件 C(n)、累计条件 Favg(n)，并分别与 Fig. 6d 拟合线和实验散点比较'],['论文需要的序列','固定总 Clifford 数，改变插入操作数。移动与静态参考使用相应的同一随机门 / DD 序列，最后施加逆门'],['论文需要的观测','测得存活概率 S(n) 与返回概率 R(n)，保留失阱的处理约定、初态制备与读出条件'],['论文需要的提取','按论文相同模型联合分析 S、R，提取随操作历史变化的单次瞬时保真度 Fₙ 与置信区间']],60,175,1160,405,[2,6.7],23);
 tx(s,'历史 IRB 和 47/195 站点结果受实现与统计口径问题影响。代码已修复，结果尚未同协议重跑。',60,611,1160,57,23,C.red,true);
}
// 21
{
 const s=slide('关键输入的来源与不确定性','论文实测值、近似模型和场景假设分别影响不同的比较结论','Fig. 4 and Methods; Ext. Data Fig. 2d; level2c and continuous configs; legacy level5 config.');
 table(s,[['输入类别','当前值 / 来源','对解释的限制'],['论文约束较强','短转移 2.4 μm、100 μs\n140 / 280 μK，η 与 Rabi 频率','有助于固定基本尺度，仍需匹配实验阶段'],['初态与束腰','转移 T 扫描 5–40 μK\nAOD 束腰假设等于 1.17 μm','温度、阱频和高能尾会改变逃逸概率'],['位点与对准差异','SLM 深度 σ=11.4% 取自论文\nAOD 深度/束腰 2%、对准 0.05 μm 为假设','采用论文散布不等于匹配其站点分布'],['噪声与透镜','RIN、指向噪声为灵敏度假设\nvₛ≈0.539 m/s 为相对场景','尚不能据此反推装置的真实损失机制']],60,175,1160,392,[1.7,4.3,3.8],23);
 tx(s,'逐段对准抖动 σ=0.01 μm 也是假设。扫描默认使 Tref=T，需要分开温度与加热的效应。',60,612,1160,54,23,C.blue,true);
}
// 22
{
 const s=slide('匹配不上的可能原因与证据强弱','已观察到的协议差异可以确认，具体实验损失机制仍需数据约束','Sources: direct protocol comparison; precision/convergence evidence; paper Methods short-transfer loss includes likely experimental imperfections; coherence noise discussion.');
 table(s,[['证据层次','原因','目前可以下的判断'],['已确认','时间、SLM 背景、Clifford/DD 与指标定义不同','阻止长运输/内部态结果直接验收'],['模型内证据','快波形易损失、慢波形过于稳定\n温度与平均加热耦合','单参数不能解释四条存活曲线'],['待实测检验','真实位置/功率响应、对准与像差\n温度尾部、强噪声事件、磁噪声','可能改变逃逸或相干性，贡献尚未定量'],['不足以解释全部','样本量与时间步误差','精细步长最大变化约 1.29 个百分点\n仍小于冻结参数下 6.78–9.76 点偏差']],60,175,1160,390,[1.6,4.2,4.2],23);
 tx(s,'论文也仅将部分短转移指数损失归因于“可能的实验不完美”。现有数据不足以进一步认定具体机制。',60,610,1160,57,23,C.red,true);
}
// 23
{
 const s=slide('达到论文匹配所需的证据','下一步验收应同时固定实验条件、整条观测曲线与误差口径','Paper Fig. 4–6 and Ext. Data Fig. 10; current numerical evidence and protocol audit.');
 table(s,[['比较对象','可检验的达标条件'],['Fig. 6 存活','同一套物理参数解释手工 200/400/600 μs 和 ML 400 μs 的完整曲线，而非逐波形换温度'],['Fig. 4 相干与门','按相同阱深、脉冲时序与空间平均方式再现 Rabi、Ramsey、XY16 和参考 RB'],['Fig. 5 / 组合 IRB','分别采用论文 610 μm 和组合 375 μm 协议，生成同条件的 S、R 与单次瞬时保真度'],['差距原因','用实测温度、AOD 位置/深度响应、噪声与对准标定约束模型，再判断残差由哪类过程造成']],60,175,1160,379,[2.1,7.6],24);
 note(s,'当前结论：流程覆盖已有进展，短转移定量匹配仍未通过，论文高保真指标尚待同协议验证。',595,C.navy,25);
}
const rawCandidate=path.join(build,`candidate-${rev}.pptx`),candidate=path.join(build,`candidate-${rev}-errorbars.pptx`);
await (await PresentationFile.exportPptx(p)).save(rawCandidate);
execFileSync(path.join(runtime,'python/python.exe'),[path.join(here,'add_paper_error_bars.py'),rawCandidate,candidate,path.join(codeDir,'reference/fig6d_vector_v2/fig6d_vectors.json')],{stdio:'inherit'});
const uniq=a=>[...new Set(a)];
const result=await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:output,pythonExecutable:path.join(runtime,'python/python.exe'),
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 requiredNativeTableOwnerSlides:uniq(tableOwners),requiredNativeChartOwnerSlides:uniq(chartOwners),materializeLiteralChartWorkbooks:true,
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...uniq(tableOwners).flatMap(n=>['--require-native-table-slide',String(n)])],
 fontPolicy:{basis:'design',families:[FONT]},verifyArtifactToolImport:true,receiptPath:path.join(build,`validation-${rev}.json`)});
console.log(JSON.stringify({output,slides:p.slides.items.length,receipt:result.receiptPath}));
await fs.writeFile(path.join(build,`slide-content-${rev}.md`),titles.map((t,i)=>`## ${i+1}. ${t}\n\n${textBySlide[i].join('\n\n')}`).join('\n\n'),'utf8');
const final=await PresentationFile.importPptx(await FileBlob.load(output));
const previews=path.join(build,`previews-${rev}`);await fs.mkdir(previews,{recursive:true});
for(let i=0;i<final.slides.items.length;i++){
 const png=await final.export({slide:final.slides.getItem(i),format:'png',scale:1});
 await fs.writeFile(path.join(previews,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
 console.log(`Rendered ${i+1}/${final.slides.items.length}`);
}
await fs.writeFile(path.join(build,`delivery-${rev}.json`),JSON.stringify({output,titles,tableOwners:uniq(tableOwners),chartOwners:uniq(chartOwners),source:ePath},null,2));
