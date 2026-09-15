// Extend the inspected v8 source deck, preserving its editable evidence.
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
import {execFileSync} from 'node:child_process';

const here=path.dirname(fileURLToPath(import.meta.url));
const root=path.resolve(here,'../..');
const modules=process.env.RUNTIME_NODE_MODULES || path.join(here,'.build/node_modules');
process.env.RUNTIME_NODE_MODULES=modules;
const require=createRequire(path.join(modules,'_resolver.cjs'));
const {FileBlob,PresentationFile}=await import(pathToFileURL(require.resolve('@oai/artifact-tool')).href);
const skill=process.env.PRESENTATIONS_SKILL_DIR;
const python=process.env.RUNTIME_PYTHON;
if(!skill||!python)throw new Error('Set PRESENTATIONS_SKILL_DIR and RUNTIME_PYTHON');
const {applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const stem='AOD_SLM_MonteCarlo_levels2c_3_4_5_2026-09-14_v9';
const source=path.join(here,'output/AOD_SLM_MonteCarlo_step_by_step_2026-09-12_v8.pptx');
const output=process.env.REPORT_PPTX || path.join(here,'output',`${stem}_reviewed.pptx`);
const build=path.join(here,'.build/v9');
await fs.mkdir(build,{recursive:true});
const evidence=JSON.parse(await fs.readFile(path.join(here,'output',`${stem}.evidence.json`),'utf8'));
const p=await PresentationFile.importPptx(await FileBlob.load(source));
const snap=await p.inspect({kind:'slide,textbox,shape,table,chart,layout',maxChars:250000});
await fs.writeFile(path.join(build,'before.ndjson'),snap.ndjson);
const records=snap.ndjson.split('\n').filter(Boolean).map(x=>JSON.parse(x));
function replaceText(page,oldText,newText){
 const r=records.find(r=>r.kind==='textbox'&&r.slide===page&&r.text===oldText);
 if(!r)throw new Error(`Missing source text on slide ${page}: ${oldText}`);
 p.resolve(r.id).text.replace(oldText,newText);
}
replaceText(1,'Step-by-step Calculation and Evidence Report','Level 2C / 3 / 4 / 5 综合计算与证据报告');
replaceText(1,'从论文数据、物理输入到存活曲线、敏感性扫描和验收','短转移、三维长运输、完整往返与有限脉冲');
replaceText(1,'2026-09-11  ·  gpt_level4_5  ·  commit 5e55bb0','2026-09-14  ·  v9  ·  5e55bb0 + 工作区更新');
replaceText(1,'Fig. 6c / 6d','Level 2C–5');
replaceText(10,'单因素敏感性扫描设计','历史参数扫描设计与复核');
replaceText(10,'每次只改变一个参数；波形、随机池、阱参数和数值步长保持一致','复核：温度扫描同时改变默认参量加热参考能量，不能视为纯初温效应');
replaceText(11,'单因素扫描结果','历史参数扫描结果');
replaceText(11,'解释：T40 和 R3e−7 能把损失推到论文量级附近；J0.10 虽然损失大，但 400/600 几乎同样差，无法解释论文的时长排序。','T40 同时改变初温与参量加热参考能量。数值保留作历史对照，不能据此单独标定初温。');
replaceText(18,'这次任务完成了什么','Level 2C 已有结果与复核边界');
replaceText(18,'用相同随机池计算 baseline、温度、RIN、对准和抖动场景，确认 T40 最接近论文的时长排序。','历史 T40 曲线较接近论文时长排序，但初温与加热参考能量同时变化，机制归因需要拆分。');

const C={bg:'#F7F9FC',navy:'#102A43',text:'#344054',blue:'#1570EF',teal:'#0E9384',orange:'#F79009',red:'#D92D20',muted:'#667085',line:'#EAECF0',white:'#FFFFFF'};
const FONT='Microsoft YaHei';
function text(s,t,x,y,w,h,size=22,color=C.text,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=t;sh.text.style={typeface:FONT,fontSize:size,color,bold,autoFit:'none'};return sh;
}
function slide(title,subtitle,sources){
 const s=p.slides.add();s.background.fill=C.bg;
 text(s,title,64,34,1152,54,31,C.navy,true);
 text(s,subtitle,64,98,1152,42,18,C.muted);
 const n=p.slides.items.length;
 text(s,`Level 2C–5 综合报告  ·  历史数值及状态详见配套报告`,64,685,1060,20,11,C.muted);
 text(s,String(n).padStart(2,'0'),1160,684,56,22,12,C.muted);
 s.speakerNotes.textFrame.setText(['Sources:',...sources.map(x=>path.join(root,x)), 'Existing numerical outputs are historical and were not rerun for this report.'].join('\n'));
 return s;
}
function table(s,values,x,y,w,h,widths,size=20){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnTracks:widths.map(value=>({mode:'fr',value})),values});
 t.borders.assign({style:'solid',fill:C.line,width:1});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){
  const cell=t.getCell(r,c);cell.fill=r===0?C.navy:(r%2?'#FFFFFF':'#F2F4F7');
  cell.text.style={typeface:FONT,fontSize:size,bold:r===0,color:r===0?C.white:C.text,verticalAlignment:'middle',autoFit:'none'};
 }
 return t;
}
function chart(s,categories,series,x=64,y=158,w=800,h=370,yTitle='Probability',yMin=0,yMax=1){
 series=series.map(r=>({...r,values:r.values.map(v=>Number(Number(v).toPrecision(12)))}));
 const c=s.charts.add('line',{position:{left:x,top:y,width:w,height:h},categories,series,lineOptions:{smooth:false},hasLegend:true,chartFill:C.white,plotAreaFill:C.white});
 applyPresentationChartFont(c,{fontFamily:FONT});
 c.legend={position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:16,fill:C.muted}};
 c.xAxis={visible:true,textStyle:{typeface:FONT,fontSize:15,fill:C.muted},line:{fill:C.line,width:1},majorGridlines:null};
 c.yAxis={visible:true,title:yTitle,min:yMin,max:yMax,numberFormatCode:'0.0',textStyle:{typeface:FONT,fontSize:15,fill:C.muted},majorGridlines:{fill:C.line,width:1}};
 return c;
}
const outbase='simulation/level2_joint_transfer/outputs/';
const srcbase='simulation/level2_joint_transfer/src/';
const cfgbase='simulation/level2_joint_transfer/configs/';
const src3=[cfgbase+'level3_3d_transport.yaml',outbase+'level3_3d_transport_lensing/demo/metrics.json'];
const src4=[cfgbase+'level4_roundtrip.yaml',outbase+'level4_end_to_end_roundtrip/demo/metrics.json'];
const src5=[cfgbase+'level5_finite_pulse_irb.yaml',outbase+'level5_finite_pulse_irb/demo/RESULT_STATUS.md'];

// 19: relationship and scope
{
const s=slide('Level 2C、3、4、5 的计算关系','统一入口管理各模块；每个层级保留明确的物理对象与统计单位',['README.md',srcbase+'simulation_workflow/cli.py']);
table(s,[['层级','计算对象','主要输出','当前衔接'],['2C','一维短距离拾取与放回','单程存活 S(n)','使用 Level 2 引擎'],['3','三维 AOD 单次长运输','留阱、机械激发、透镜效应','给 Level 4 提供运输模型'],['4','完整拾取—运输—返回—放回','连续多轮成功率与相位对比度','使用 Level 2 冻结波形与 Level 3'],['5','冻结轨迹上的有限脉冲量子态','PTM、RB/IRB、站点统计','使用 Level 4 协议与轨迹池']],64,160,1152,350,[.7,2,2,2.4],20);
text(s,'Level 2C 最新波形尚未接入 Level 4/5。',64,540,1130,42,25,C.blue,true);
text(s,'本版补齐现有实现的入口、结果和边界；新的连续物理协议需要另外定义并验证。',64,596,1130,58,21);
}
// 20: Level 3 model
{
const s=slide('Level 3：三维长运输模型','AOD 中的热初态，经一次运输和 100 μs 静止保持后判断是否仍然束缚',src3);
table(s,[['输入 / 方法','当前实现'],['原子与光阱','Cs-133；1055 nm；w = 1.17 μm；D/kB = 280 μK'],['初态与传播','T = 5 μK；三维谐振提议 + 束缚拒绝；velocity-Verlet'],['几何与距离','直线 / 对角；270 / 510 / 610 μm（路径总长度）'],['轨迹与时长','sine / 四段 constant-jerk / minimum-jerk；400–2000 μs'],['透镜效应','扫频速度决定两轴焦点偏移；有效深度下降、轴向分岔'],['判据与步长','末态 E < 0；正式 dt = 0.05 μs；0.025 μs 复核']],64,160,1152,390,[1.4,4.5],21);
text(s,'对角 610 μm：每轴位移约 431 μm。',64,578,1100,40,24,C.blue,true);
text(s,'本级不含 SLM 拾取/放回、重复往返或量子内部态。',64,626,1100,32,20);
}
// 21: native numerical evidence
{
const s=slide('Level 3：长运输的经典留阱结果','已有 96 格点粗扫描是配置范围的子集；下图为 adiabatic-sine，128 shots/点',[...src3,outbase+'level3_3d_transport_lensing/demo/coarse_scan.csv']);
const rows=evidence.level3.coarse;
const pick=(geometry,distance,severity)=>rows.filter(r=>r.profile==='adiabatic_sine'&&r.geometry===geometry&&Number(r.distance_um)===distance&&Number(r.severity)===severity).sort((a,b)=>a.duration_us-b.duration_us);
const straight=pick('straight',510,1),diag=pick('diagonal',610,1),off=pick('straight',510,0);
const cats=straight.map(r=>`${Number(r.duration_us)} μs`);
for(const arr of [diag,off])if(arr.map(r=>r.duration_us).join()!==straight.map(r=>r.duration_us).join())throw new Error('Coarse scan x mismatch');
chart(s,cats,[{name:'Straight 510 / sev1',values:straight.map(r=>+r.retention_fraction),line:{fill:C.orange,width:3}},{name:'Diagonal 610 / sev1',values:diag.map(r=>+r.retention_fraction),line:{fill:C.blue,width:3}},{name:'Straight / lensing off',values:off.map(r=>+r.retention_fraction),line:{fill:C.teal,width:2}}]);
text(s,'冻结主条件',900,174,302,42,25,C.navy,true);
text(s,'8 个条件\n每个 2000/2000 留阱\n95% 区间约\n[0.9981, 1.0000]',900,230,302,194,23);
text(s,'范围：510 μm 直线、610 μm 对角，1.6 ms，sine / constant-jerk。',900,430,302,102,19);
text(s,'主条件未观察到失阱，不能推断实验零损失；留阱率也不是量子门保真度。',64,590,1152,64,23,C.blue,true);
}
// 22: calibrated limitations and compensation
{
const s=slide('Level 3：透镜标定与深度补偿边界','默认 severity 是相对场景；器件临界速度没有独立实验标定',src3);
text(s,'severity = 1 对应 v_s ≈ 0.53917 m/s',64,164,1120,44,26,C.blue,true);
text(s,'该值由 610 μm / 1.6 ms 的参考对角轨迹推导，不能作为论文实测硬件参数。',64,214,1120,56,22);
table(s,[['510 μm 直线 / 1.6 ms','未补偿','探索性补偿'],['功率倍率上限','1','2'],['最低有效深度','118.4 μK','236.8 μK'],['留阱','500/500','500/500'],['中位激发能增量','15.833 μK','17.647 μK']],64,290,1152,240,[2.4,1.6,1.6],21);
text(s,'所需峰值倍率约 2.36，超过上限；改善阱深没有自动降低加热。',64,556,1150,43,23,C.navy,true);
text(s,'四段 constant-jerk 的 T⁻⁶ 标度仍需与论文的具体波形定义核对。',64,614,1150,40,21);
}
// 23: protocol A/B
{
const s=slide('Level 4：完整组合往返协议','短距离 split/merge 与长距离去程/回程都包含在一轮操作中',src4);
table(s,[['Protocol A 阶段','时长 / μs','AOD 路程'],['初始保持','100','0'],['pick-up','850','0'],['split','400','2.4 μm'],['长距离去程','1450','375 μm'],['远端保持窗口','54','0'],['长距离回程','1450','375 μm'],['merge','400','2.4 μm'],['drop-off','850','0'],['最终保持','100','0']],64,160,700,470,[2.5,1.3,1.5],20);
text(s,'Protocol A',820,174,360,38,26,C.navy,true);
text(s,'总时长 5654 μs\nSLM：180→60→180 μK\nAOD 峰值：280 μK',820,222,365,135,23);
text(s,'Protocol B',820,378,360,38,26,C.navy,true);
text(s,'总时长 3954 μs\nSLM 恒 140 μK\nLevel 2 冻结波形自带\n短移动，不重复 split/merge',820,426,372,151,22);
text(s,'54 μs 只表示保持窗口，没有实际两比特门。',820,590,372,64,19,C.muted);
}
// 24: classical and coherence distinct
{
const s=slide('Level 4：单轮成功与条件相干对比度','2000 shots/条件，dt = 0.05 μs；相位模型采用理想瞬时 DD',src4);
table(s,[['单轮条件','成功数','成功率'],['Protocol A / lensing','2000/2000','1.0000'],['Protocol A / lensing off','2000/2000','1.0000'],['Protocol B / frozen','2000/2000','1.0000'],['同时间静止保持','2000/2000','1.0000']],64,176,700,305,[2.4,1.3,1.1],20);
const phase=evidence.level4.contrast;
table(s,[['Protocol A DD','C_cond'],...['none','spin_echo','xy4_per_long_move'].map((k,i)=>[['None','Spin echo','每段长移动 XY4'][i],phase['protocol_a_nominal_lensing|'+k].contrast.toFixed(4)])],800,176,416,305,[2,1],20);
text(s,'成功率的 Wilson 95% 区间约 [0.9981, 1.0000]。',64,520,1148,40,23);
text(s,'相位 φ = ∫ y(t)ηU/ħ dt；C_cond = |⟨exp(iφ)⟩|。',64,574,1148,42,23,C.blue,true);
text(s,'条件对比度与 P(success)×C_cond 都不能直接替代论文 IRB 保真度。',64,626,1148,32,20);
}
// 25: continuous rounds
{
const s=slide('Level 4：连续重复往返','每轮继承原子的位置与速度；下图为 500-shot Protocol A 历史结果',[...src4,outbase+'level4_end_to_end_roundtrip/demo/repeated_round_metrics.csv']);
const rows=evidence.level4.repeated.filter(r=>r.protocol==='protocol_a').sort((a,b)=>a.round-b.round);
chart(s,rows.map(r=>r.round),[['None','contrast_none',C.orange],['Spin echo','contrast_spin_echo',C.teal],['XY4','contrast_xy4',C.blue]].map(([name,key,color])=>({name,values:rows.map(r=>+r[key]),line:{fill:color,width:3}})),64,165,806,373,'C_cond');
text(s,'第 10 轮',910,181,285,40,26,C.navy,true);
const last=rows.at(-1);
text(s,`成功 ${last.absorbing_success}/500\nC_none = ${(+last.contrast_none).toFixed(4)}\nC_echo = ${(+last.contrast_spin_echo).toFixed(4)}\nC_XY4 = ${(+last.contrast_xy4).toFixed(4)}`,910,235,292,180,22);
text(s,'40 轮长尾\n128/128 成功',910,447,285,92,24,C.navy,true);
text(s,'留阱与相干性反映不同退化通道；上述曲线依赖当前重复协议及相位统计定义。',64,580,1152,75,23);
}
// 26: finite pulse scope
{
const s=slide('Level 5：有限脉冲单比特动力学','在 Level 4 冻结势能轨迹上计算内部态，运动轨迹不受量子态反作用',src5);
text(s,'H/ħ = [Δσz + Ω(cosφ σx + sinφ σy)] / 2',64,159,1152,52,28,C.blue,true);
table(s,[['部分','内容'],['驱动与光移','Ω/(2π)=24.611 kHz；η_SLM=η_AOD=1.3×10⁻⁴'],['有限脉冲','bare / SCROFULOUS；解析 SU(2) 子步传播'],['动力学解耦','none、spin-echo、XY4/XY16、每段长移动 XY4'],['通道表征','六输入态、PTM、条件平均门保真度与 bootstrap 区间'],['随机基准测试','24 个 Clifford；reference RB；固定 N=80 的 M 扫描'],['站点与噪声','47/195 站点；共同/局域噪声与站点参数可配置']],64,236,1152,348,[1.2,4.8],21);
text(s,'默认噪声与站点散布为零；完整 transformed Clifford frame 未接入正式传播。',64,614,1152,45,21);
}
// 27: qualified channel results
{
const s=slide('Level 5：单通道 PTM 的历史结果','以下结果不经过 IRB 索引路径；目标门定义决定 F_avg 的含义',[...src5,outbase+'level5_finite_pulse_irb/demo/metrics.json']);
const keys=['none','spin_echo','xy4_per_long_move','xy16'];
const vals=keys.map(k=>Number(evidence.level5.channel[`protocol_a|${k}|nominal`].f_avg.toPrecision(12)));
const c=s.charts.add('bar',{position:{left:64,top:169,width:802,height:384},categories:['None','Spin echo','XY4 per move','XY16'],series:[{name:'Historical F_avg vs I',values:vals,valuesFormatCode:'0.0000',fill:C.blue,dataLabelOverrides:vals.map((v,idx)=>({idx,text:v.toFixed(4),position:'outEnd',textStyle:{typeface:FONT,fontSize:17,fill:C.text}}))}],barOptions:{direction:'column',grouping:'clustered'},hasLegend:false,dataLabels:{showValue:true,position:'outEnd',textStyle:{typeface:FONT,fontSize:17,fill:C.text}}});
applyPresentationChartFont(c,{fontFamily:FONT});
c.xAxis={textStyle:{typeface:FONT,fontSize:17,fill:C.muted}};c.yAxis={min:0,max:1,title:'F_avg relative to I',textStyle:{typeface:FONT,fontSize:16,fill:C.muted}};
text(s,'目标门：恒等 I',906,177,300,40,25,C.navy,true);
text(s,'spin-echo 的净 π 操作\n与恒等门不同。\n约 1/3 的读数不能\n直接解释为失相干。',906,237,300,196,23);
text(s,'各条件经典 survival = 1\n数值来自旧通道输出。',906,463,300,86,21,C.muted);
text(s,'比较 DD 性能前，需要统一目标门或移除理想净旋转；这些数值不是论文实验 IRB。',64,589,1152,66,23,C.blue,true);
}
// 28: RB design and validity
{
const s=slide('Level 5：RB/IRB 设计与结果状态','旧多站点 IRB 需要重跑，当前报告不发布其保真度或阵列性能结论',[...src5,srcbase+'level5_finite_pulse_irb/rb_fitting.py',srcbase+'level5_finite_pulse_irb/interleaved_rb.py']);
table(s,[['项目','当前设计 / 状态'],['Reference RB','Clifford 长度 1–1000；72 序列/点；47 站点'],['运输插入','固定 N=80；扫描 M=0–80；每点 72 strings'],['每个 move','默认完整组合往返；独立抽取 Level 4 单轮冻结轨迹'],['历史 IRB / 阵列','移动索引串写影响多站点结果；修复后尚未全量重跑'],['拟合解释','p(M) 与 reference p(n) 自变量不同，禁止套标准比值']],64,165,1152,330,[1.4,4.6],22);
text(s,'历史 F_avg = 1.001650 不可解释。',64,528,1152,45,28,C.red,true);
text(s,'新版保留原始数据与状态说明；报告衰减参数时也必须检验目标、设计和拟合质量。',64,590,1152,65,22);
}
// 29: counts
{
const s=slide('运输距离与“60 次”的统一口径','总路程统计实际运动段；静止保持段不增加距离',[cfgbase+'level2c_pickup_survival.yaml',...src3,...src4,...src5]);
table(s,[['场景','一次操作','次数 / 累计路程'],['Level 2C','2.4 μm 单程短转移','n=60：144 μm；30 次短往返'],['Level 3','一次长运输 270/510/610 μm','无默认重复往返'],['Level 4','2.4 + 375 + 375 + 2.4 μm','每轮 754.8 μm；默认 10/40 轮'],['Level 5 默认 M=60','60 个完整组合往返通道','45,288 μm = 45.288 mm']],64,171,1152,339,[1.3,2.8,3.1],22);
text(s,'Level 4/5 的完整组合操作包含长距离运输段。',64,543,1152,43,26,C.blue,true);
text(s,'Level 5 的路程是操作定义的累计计数；它没有连续传播同一个经典原子 60 轮。',64,607,1152,50,22);
}
// 30: unified execution
{
const s=slide('仓库统一入口与产物追溯','从仓库根目录运行，仍由各层级原 CLI 执行其物理模型',['run_simulation.py',srcbase+'simulation_workflow/cli.py','README.md']);
table(s,[['命令','用途'],['python run_simulation.py list','查看层级、操作单位和依赖'],['python run_simulation.py check --json','检查源文件 / 配置 / 历史指标，并输出 SHA-256'],['python run_simulation.py run 3 --stage preflight','Level 3 数值预检查'],['python run_simulation.py run 4 --stage single_round','Level 4 单轮阶段'],['python run_simulation.py run 5 --stage preflight','Level 5 预检查与轨迹池准备']],64,171,1152,341,[3.7,2.7],19);
text(s,'默认写入新目录；分阶段续跑需显式指定同一个 --output-dir。',64,543,1152,44,23,C.blue,true);
text(s,'新上游结果需在配置中明确路径；自动发现不保证选中最新目录。--dry-run 只显示命令。',64,608,1152,51,21);
}
// 31: actual verification
{
const s=slide('本轮代码修复与验证','单元测试证明指定行为，不能代替正式规模重跑或实验校准',[srcbase+'level5_finite_pulse_irb/interleaved_rb.py',srcbase+'level5_finite_pulse_irb/rb_fitting.py','simulation/level2_joint_transfer/tests/test_level5_irb_regression.py']);
table(s,[['修复 / 验证项','结果'],['多站点移动索引','sequence × n_sites + site，逐站点唯一更新'],['插入移动后的时间','将 move 总时长计入后续门的回放时间'],['IRB 公式适用性','固定 N 的 M 扫描不输出标准比值保真度'],['解析回归','1/3/47 站点；M=0/1/2/60；状态、吸收损失、耗时'],['模型与新增回归','89 项通过；8 项慢速 / CLI / preflight 未运行'],['统一入口','6 项通过；本轮合计 95 项测试通过']],64,169,1152,367,[2,4],21);
text(s,'尚未重跑 2000-shot 正式验证与完整 47/195 站点 RB。',64,572,1152,47,25,C.blue,true);
text(s,'旧 IRB 与阵列数值已标记失效待重跑，原始输出保留用于追溯。',64,629,1152,33,20);
}
// 32: complete scope and remaining physical work
{
const s=slide('当前覆盖范围与后续物理验证','Level 2C、3、4、5 均纳入仓库入口和本版报告',['README.md','reports/stage_review/output/'+stem+'.md']);
text(s,'已覆盖',64,166,330,46,28,C.blue,true);
text(s,'短转移存活与波形比较\n三维 AOD 长运输和透镜效应\n完整组合往返与连续运动相位\n有限脉冲、PTM 和多站点 RB/IRB 框架',64,227,1100,183,25);
text(s,'优先补齐的物理证据',64,439,1120,43,28,C.navy,true);
text(s,'修复后的 Level 5 全量重跑与目标门核查\n论文波形定义、绝对 AOD 标定及噪声系数核对\nLevel 2C 最新波形接入三维全协议后的连续性验证',64,498,1120,143,24);
}

const draft=path.join(build,'candidate.pptx');
const candidate=path.join(build,'candidate.with-source-data.pptx');
await (await PresentationFile.exportPptx(p)).save(draft);
execFileSync(python,[path.join(here,'restore_chart_workbooks.py'),source,draft,candidate],{stdio:'inherit'});
const requiredTables=[2,3,4,5,6,7,9,10,16,19,20,22,23,24,26,28,29,30,31];
const result=await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:output,pythonExecutable:python,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 explicitTotalSlideCount:32,requiredNativeTableOwnerSlides:requiredTables,
 requiredNativeChartOwnerSlides:[6,8,9,11,12,13,14,15,17,21,25,27],
 requiredEmbeddedWorkbookChartOwnerSlides:[],
 materializeLiteralChartWorkbooks:true,
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...requiredTables.flatMap(n=>['--require-native-table-slide',String(n)])],
 fontPolicy:{basis:'reference',families:[FONT,'Calibri'],referencePath:source,referenceSha256:crypto.createHash('sha256').update(await fs.readFile(source)).digest('hex')},
 verifyArtifactToolImport:true,receiptPath:path.join(build,path.basename(output)+'.validation.json')});
console.log(JSON.stringify({output,slides:p.slides.items.length,validationReceipt:result.receiptPath}));
// Render the actual final file, not just in-memory authored objects.
const final=await PresentationFile.importPptx(await FileBlob.load(output));
const previews=path.join(build,'final-previews');await fs.mkdir(previews,{recursive:true});
for(let i=0;i<final.slides.items.length;i++){
 const blob=await final.export({slide:final.slides.getItem(i),format:'png',scale:1});
 await fs.writeFile(path.join(previews,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const montage=await final.export({format:'webp',montage:true,scale:.5});
await fs.writeFile(path.join(previews,'montage.webp'),new Uint8Array(await montage.arrayBuffer()));
console.log(`Rendered ${final.slides.items.length} final slides`);
