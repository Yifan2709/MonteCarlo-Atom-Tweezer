// Extend v9 with actual continuous Monte Carlo evidence, preserving source charts.
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
import {execFileSync} from 'node:child_process';
const here=path.dirname(fileURLToPath(import.meta.url)), root=path.resolve(here,'../..');
const modules=process.env.RUNTIME_NODE_MODULES;
const require=createRequire(path.join(modules,'_resolver.cjs'));
const {FileBlob,PresentationFile}=await import(pathToFileURL(require.resolve('@oai/artifact-tool')).href);
const skill=process.env.PRESENTATIONS_SKILL_DIR, python=process.env.RUNTIME_PYTHON;
const {applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const stem='AOD_SLM_continuous_level2c_345_2026-09-14_v10';
const source=path.join(here,'output/AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9_reviewed.pptx');
const output=process.env.REPORT_PPTX || path.join(here,'output',stem+'_reviewed.pptx');
const build=path.join(here,'.build/v10');
await fs.mkdir(build,{recursive:true});
const e=JSON.parse(await fs.readFile(path.join(here,'output',stem+'.evidence.json'),'utf8'));
const p=await PresentationFile.importPptx(await FileBlob.load(source));
const snap=await p.inspect({kind:'slide,textbox,shape,table,chart',maxChars:300000});
const records=snap.ndjson.split('\n').filter(Boolean).map(x=>JSON.parse(x));
await fs.writeFile(path.join(build,'before.ndjson'),snap.ndjson);
function replace(page,oldText,newText){
 const item=records.find(r=>r.kind==='textbox'&&r.slide===page&&r.text===oldText);
 if(!item)throw new Error(`Missing text on page ${page}: ${oldText}`);
 p.resolve(item.id).text.replace(oldText,newText);
}
replace(1,'短转移、三维长运输、完整往返与有限脉冲','连续全流程重算与 Fig. 6d 存活曲线验收');
replace(1,'2026-09-14  ·  v9  ·  5e55bb0 + 工作区更新','2026-09-14  ·  v10  ·  连续仿真与独立复核');
replace(19,'Level 2C 最新波形尚未接入 Level 4/5。','新增 full 入口已接入 Level 2C 连续全流程。');
replace(19,'本版补齐现有实现的入口、结果和边界；新的连续物理协议需要另外定义并验证。','上表保留独立模块关系。第 33–43 页报告新的连续协议、重算结果及验收。');
replace(31,'本轮代码修复与验证','v9 代码修复与验证记录');
replace(31,'尚未重跑 2000-shot 正式验证与完整 47/195 站点 RB。','此页为 v9 记录。新增 full 计算与 4096-shot 复核见后续页。');
replace(32,'当前覆盖范围与后续物理验证','独立模块范围与连续扩展');
replace(32,'Level 2C、3、4、5 均纳入仓库入口和本版报告','第 1–32 页保留既有模型和历史证据，后续页给出连续全流程的新结果');
const pending32=records.find(r=>r.kind==='textbox'&&r.slide===32&&r.text.includes('修复后的 Level 5 全量重跑'));
if(!pending32)throw new Error('Missing slide 32 scope text');
const scope32=p.resolve(pending32.id);
scope32.text.replace('修复后的 Level 5 全量重跑与目标门核查','原有 RB 独立协议仍需单独重跑与目标门核查');
scope32.text.replace('论文波形定义、绝对 AOD 标定及噪声系数核对','论文波形定义、绝对 AOD 标定及噪声系数仍待核对');
scope32.text.replace('Level 2C 最新波形接入三维全协议后的连续性验证','最新波形连续接入、重算及数值复核已在后续页完成');

const C={bg:'#F7F9FC',navy:'#102A43',text:'#344054',blue:'#1570EF',teal:'#0E9384',orange:'#F79009',red:'#D92D20',muted:'#667085',line:'#EAECF0',white:'#FFFFFF'};
const FONT='Microsoft YaHei';
const names=['手工 200 μs','手工 400 μs','手工 600 μs','论文 ML 400 μs'];
function text(s,t,x,y,w,h,size=22,color=C.text,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=t;sh.text.style={typeface:FONT,fontSize:size,color,bold,autoFit:'none'};return sh;
}
function slide(title,subtitle,files=[]){
 const s=p.slides.add();s.background.fill=C.bg;
 text(s,title,64,34,1152,54,31,C.navy,true);
 text(s,subtitle,64,98,1152,44,18,C.muted);
 text(s,'连续全流程重算  ·  参数与数值来源见配套报告 v10',64,685,1060,20,11,C.muted);
 text(s,String(p.slides.items.length).padStart(2,'0'),1160,684,56,22,12,C.muted);
 s.speakerNotes.textFrame.setText(['Sources:',path.join(here,'output',stem+'.md'),path.join(here,'output',stem+'.evidence.json'),...files.map(f=>path.join(root,f))].join('\n'));
 return s;
}
function table(s,values,x,y,w,h,widths,size=21){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnTracks:widths.map(value=>({mode:'fr',value})),values});
 t.borders.assign({style:'solid',fill:C.line,width:1});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){
  const cell=t.getCell(r,c);cell.fill=r===0?C.navy:(r%2?C.white:'#F2F4F7');
  cell.text.style={typeface:FONT,fontSize:size,bold:r===0,color:r===0?C.white:C.text,verticalAlignment:'middle',autoFit:'none'};
 }
 return t;
}
function chart(s,categories,series){
 const c=s.charts.add('line',{position:{left:64,top:165,width:820,height:420},categories,
   series:series.map(r=>({...r,values:r.values.map(v=>Number(Number(v).toPrecision(12)))})),
   lineOptions:{smooth:false},hasLegend:true,chartFill:C.white,plotAreaFill:C.white});
 applyPresentationChartFont(c,{fontFamily:FONT});
 c.legend={position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:14,fill:C.muted}};
 c.xAxis={visible:true,textStyle:{typeface:FONT,fontSize:13,fill:C.muted},line:{fill:C.line,width:1},majorGridlines:null};
 c.yAxis={visible:true,title:'Survival S(n)',min:0,max:1,numberFormatCode:'0.0',textStyle:{typeface:FONT,fontSize:15,fill:C.muted},majorGridlines:{fill:C.line,width:1}};
 return c;
}
const end=r=>r.survival[r.comparison_raw.n.at(-1)], ref=r=>r.comparison_raw.reference.at(-1), f=v=>Number(v).toFixed(4);
// 33
{
const s=slide('连续全流程接入后的验收结果','原始四条波形保持不变，使用独立样本检查整条存活曲线');
text(s,'本次测试未找到同时拟合四条曲线的共同参数',64,162,1152,48,29,C.red,true);
table(s,[['波形','逐曲线扫描选择 T','4096-shot 末点','论文末点','整曲线通过 ±0.02'],
 ...e.selected.map((r,i)=>[names[i],r.temperature_uK+' μK',f(end(r)),f(ref(r)),'否'])],64,235,1152,287,[1.7,1.8,1.6,1.3,2],21);
text(s,'表中不同 T 仅用于诊断，不能视为同一实验条件。',64,549,1152,42,25,C.navy,true);
text(s,'200 μs 比较 n=31，其余比较 n=60。独立 seed=23003，dt=0.0125 μs。',64,608,1152,50,21);
}
// 34
{
const s=slide('连续协议与距离口径','matched 沿用 Fig. 6d，extended 另加长距离去程与回程');
table(s,[['项目','matched 短转移对照','extended 组合运输'],
['每轮阶段','拾取、100 μs 等待、放回','另插入去程、远端保持、回程'],
['新增长距离','无','375 μm + 375 μm，xy 对角'],
['400 μs 波形单轮时长','900 μs','3854 μs'],
['手工波形单轮名义路程','4.8 μm','754.8 μm'],
['n=60 的含义','30 轮，共约 144 μm','30 轮，共约 22.644 mm'],
['对应论文验收','沿用原 Level 2C 标准','作为扩展预测单列']],64,165,1152,387,[1.7,2.6,2.8],21);
text(s,'每个原子的运动状态、内部态和累计加热贯穿所有段与所有轮。',64,578,1152,44,24,C.blue,true);
text(s,'失阱后停止演化。SLM 开关改变势场，位置和速度不重置。',64,629,1152,34,20);
}
// 35–38
for(let i=0;i<4;i++){
 const r=e.selected[i],b=e.base['5'].baseline_1d[i],m=e.base['5'].matched[i],cp=r.comparison_raw,n=cp.n;
 const s=slide(names[i]+'：连续存活曲线','按原始 Fig. 6d 数字化点验收，要求所有点偏差不超过 0.02',[r.result_file,b.result_file,m.result_file]);
 chart(s,n.map(String),[
  {name:'Paper',values:cp.reference,line:{fill:C.navy,width:3}},
  {name:'1D T5',values:n.map(n=>b.survival[n]),line:{fill:C.muted,width:2}},
  {name:'3D T5',values:n.map(n=>m.survival[n]),line:{fill:C.blue,width:2}},
  {name:`3D T${r.temperature_uK} / 4096`,values:cp.simulated,line:{fill:C.red,width:3}}
 ]);
 text(s,'冻结参数复核',920,174,288,42,25,C.navy,true);
 text(s,`T = Tref = ${r.temperature_uK} μK\n末点 ${f(end(r))}\n论文 ${f(ref(r))}`,920,231,288,146,24);
 text(s,`最大偏差 ${f(cp.max_abs_deviation)}\n最差位置 n=${e.confidence[i].worst_n}\n整曲线未通过`,920,407,288,142,22,C.red,true);
 text(s,`灰/蓝线为 512-shot、dt=0.05 μs 对照。红线为 4096-shot、dt=0.0125 μs 独立复核。`,64,610,1152,52,20);
}
// 39
{
const s=slide('共同参数扫描与波形间矛盾','候选 T 与默认 Tref 同步改变，结果不能单独解释为温度标定');
table(s,[['T = Tref / μK','四曲线共同 RMSE','全部通过 ±0.02'],
 ...e.temperature_scan.map(r=>[String(r.temperature_uK),r.pooled_rmse.toFixed(5),'否'])],64,166,725,428,[1.8,2,2],22);
text(s,'共同最优格点 15 μK',836,175,366,48,26,C.blue,true);
text(s,'四条曲线等权\n256 shots，seed=23001\n冻结后用新随机池复核',836,238,366,138,22);
text(s,'逐曲线选择分别为\n5 / 15 / 25 / 25 μK\n共同参数无法兼顾',836,416,366,142,24,C.navy,true);
text(s,'扫描范围有限。未通过表示当前已测试模型和参数不足以复现，不是证明任何补充模型都无法拟合。',64,617,1152,49,20);
}
// 40
{
const s=slide('长距离运输对存活的影响','512-shot 独立复核，400 μs 波形每轮额外演化 2954 μs');
table(s,[['波形末点','论文','短协议 T5','长协议 T5','长协议 T40'],
 ...names.map((name,i)=>[name,f(ref(e.base['5'].matched[i])),f(end(e.base['5'].matched[i])),f(end(e.base['5'].extended[i])),f(end(e.base['40'].extended[i]))])],64,170,1152,287,[2,1.1,1.5,1.5,1.5],22);
text(s,'长运输降低了部分存活率，600 μs 和 ML 的低温结果仍接近 1。',64,489,1152,52,25,C.blue,true);
text(s,'extended 已改变 Fig. 6d 的实验协议，不能把增加路程后的接近作为原波形拟合通过。',64,565,1152,69,23);
}
// 41
{
const s=slide('有限脉冲接入后的内部态结果','同一条连续轨迹同时累计无 DD 相位和有限 DD 的 SU(2) 演化');
const rr=[['matched',5],['matched',40],['extended',5],['extended',40]].map(([a,t])=>({a,t,r:e.base[String(t)][a][1]}));
table(s,[['400 μs / n=60','存活 S','C 无 DD','C 有限 DD','Favg 有限 DD'],
 ...rr.map(({a,t,r})=>{let q=r.quantum.checkpoints.at(-1);return [`${a} T${t}`,f(r.survival.at(-1)),f(q.contrast_none),f(q.contrast_dd),f(q.favg_dd)];})],64,166,1152,290,[2.1,1.3,1.5,1.5,1.6],21);
text(s,'纯内部态幺正演化不会补出原子损失。',64,492,1152,46,28,C.blue,true);
text(s,'本模型没有自旋依赖机械力。DD 影响相干性，机械存活率与无 DD 完全相同。',64,550,1152,53,23);
text(s,'Favg 相对理想脉冲序列、只对存活原子取平均。此处不使用旧 RB/IRB 比值拟合。',64,617,1152,43,20);
}
// 42
{
const s=slide('步长与统计复核','相同 4096 初态、无序和逐段抖动，比较 dt=0.025 与 0.0125 μs');
table(s,[['最大 |ΔS(n)|','T5','T15','T25'],...names.map((name,i)=>[name,...[5,15,25].map(t=>f(e.convergence.find(r=>r.temperature_uK===t&&r.waveform===e.selected[i].waveform).max_survival_change))])],64,172,1152,286,[2.1,1.4,1.4,1.4],23);
const max=Math.max(...e.convergence.map(r=>r.max_survival_change));
text(s,`精细检查最大曲线变化 ${f(max)}，原 128-shot 检查最大变化 ${f(Math.max(...e.original_step_check.map(r=>r.max_survival_change)))}。`,64,493,1152,51,24,C.blue,true);
text(s,'4096 shots 在 S≈0.5 处的单点 95% 区间半宽约 0.0153。',64,561,1152,42,23);
text(s,'轨迹标签可能因混沌演化改变。分别报告逐原子标签差和系综曲线差，避免混为一谈。',64,618,1152,43,20);
}
// 43
{
const s=slide('仍未拟合的原因与模型边界','本次已完成连续接入、原标准复核及冻结参数的大样本确认');
table(s,[['观察 / 假设','对结论的影响'],
['三维热能与有限阱深','同温度有更多运动能量，T40 参数出现过度损失'],
['速度相关的柱面透镜','对快速 200 μs 波形影响大，对慢波形不足以产生所需衰减'],
['温度与加热参考能量耦合','共同扫描参数无法同时兼顾四条曲线'],
['装置与波形仍有未标定项','绝对 AOD 透镜、噪声谱、实际波形响应仍需数据约束'],
['内部态与机械态的耦合范围','现有有限脉冲模型没有自旋依赖机械力或 Raman 损失']],64,166,1152,358,[2.2,4.7],22);
text(s,`${e.result_case_count} 个实际计算格点已归档，104 项相关代码检查通过。`,64,557,1152,45,25,C.blue,true);
text(s,'full 入口可重现本流程。原独立 RB 协议及实验损失标定需要另外的验证。',64,617,1152,44,22);
}

const draft=path.join(build,'candidate.pptx'),candidate=path.join(build,'candidate.with-source-data.pptx');
await (await PresentationFile.exportPptx(p)).save(draft);
execFileSync(python,[path.join(here,'restore_chart_workbooks.py'),source,draft,candidate],{stdio:'inherit'});
const tableSlides=[2,3,4,5,6,7,9,10,16,19,20,22,23,24,26,28,29,30,31,33,34,39,40,41,42,43];
const result=await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:output,pythonExecutable:python,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 explicitTotalSlideCount:43,requiredNativeTableOwnerSlides:tableSlides,
 requiredNativeChartOwnerSlides:[6,8,9,11,12,13,14,15,17,21,25,27,35,36,37,38],
 requiredEmbeddedWorkbookChartOwnerSlides:[],materializeLiteralChartWorkbooks:true,
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableSlides.flatMap(n=>['--require-native-table-slide',String(n)])],
 fontPolicy:{basis:'reference',families:[FONT,'Calibri'],referencePath:source,referenceSha256:crypto.createHash('sha256').update(await fs.readFile(source)).digest('hex')},
 verifyArtifactToolImport:true,receiptPath:path.join(build,path.basename(output)+'.validation.json')});
console.log(JSON.stringify({output,slides:p.slides.items.length,validationReceipt:result.receiptPath}));
const final=await PresentationFile.importPptx(await FileBlob.load(output));
const previews=path.join(build,'final-previews');await fs.mkdir(previews,{recursive:true});
for(let i=0;i<final.slides.items.length;i++){
 const blob=await final.export({slide:final.slides.getItem(i),format:'png',scale:1});
 await fs.writeFile(path.join(previews,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const montage=await final.export({format:'webp',montage:true,scale:.5});
await fs.writeFile(path.join(previews,'montage.webp'),new Uint8Array(await montage.arrayBuffer()));
console.log(`Rendered ${final.slides.items.length} final slides`);
