import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
const here=path.dirname(fileURLToPath(import.meta.url)),root=path.resolve(here,'../..');
const build=path.join(here,'.build','visual-v2');
const runtime=process.env.RUNTIME_NODE_MODULES||'/Users/imok/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.RUNTIME_NODE_MODULES=runtime;
const require=createRequire(path.join(runtime,'_resolve.cjs'));
const {Presentation,PresentationFile}=await import(pathToFileURL(require.resolve('@oai/artifact-tool')).href);
const skill=process.env.PRESENTATIONS_SKILL_DIR||'/Users/imok/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations';
const {applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const python=process.env.RUNTIME_PYTHON||'/Users/imok/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const D=JSON.parse(await fs.readFile(path.join(build,'data.json'),'utf8'));
const p=Presentation.create({slideSize:{width:1280,height:720}});
const FONT='Hiragino Sans GB';
const C={text:'#182F3D',muted:'#536977',blue:'#186DA0',orange:'#D36E2B',green:'#21877A',red:'#B44147',gray:'#8395A1',pale:'#ABC4D3',grid:'#E0E7EC'};
const cols=[C.blue,C.orange,C.green,C.red];
const notes=[],chartOwners=[],tableOwners=[];
const SRC={guide:'docs/monte_carlo_plain_guide.md',init:'simulation/level2_joint_transfer/src/continuous_transfer/initialization.py',protocol:'simulation/level2_joint_transfer/src/continuous_transfer/protocol.py',engine:'simulation/level2_joint_transfer/src/continuous_transfer/engine.py',noise:'simulation/level2_joint_transfer/src/continuous_transfer/noisy_survival.py',cli:'simulation/level2_joint_transfer/src/continuous_transfer/cli.py',phase:'simulation/level2_joint_transfer/src/level4_roundtrip_coherence/dephasing.py',report:'reports/movement_effects_v1_calibration_20260916/REPORT.md',physics:'reports/movement_effects_v1_calibration_20260916/PHYSICS.md',params:'reports/movement_effects_v1_calibration_20260916/PARAMETERS.md',valid:'reports/movement_effects_v1_calibration_20260916/VALIDATION.md',data:'reports/movement_effects_v1_calibration_20260916/survival_curves.csv',points:'reports/movement_effects_v1_calibration_20260916/point_residuals.csv',metrics:'reports/movement_effects_v1_calibration_20260916/metrics.csv',trajectory:'simulation/level1_transfer_1d/outputs/level1_transfer_1d/demo/thermal_example_trajectory.csv',paper:'https://arxiv.org/html/2403.12021v4',savard:'https://doi.org/10.1103/PhysRevA.56.R1095'};
function text(s,t,x,y,w,h,size=27,color=C.text,bold=false){
 const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 sh.text=t;sh.text.style={typeface:FONT,fontSize:size,color,bold,autoFit:'none',verticalAlignment:'top'};return sh;
}
function page(title,caption,note,refs=['guide']){
 const s=p.slides.add();s.background.fill='#FFFFFF';const n=p.slides.items.length;
 text(s,title,58,35,1165,67,43,C.text,true);
 if(caption)text(s,caption,62,112,1155,48,23,C.muted);
 text(s,n>24?`附录 ${n-24}`:String(n).padStart(2,'0'),1124,676,108,28,16,C.gray);
 const sources=refs.map(v=>SRC[v]||v);
 s.speakerNotes.textFrame.setText(`讲解\n${note}\n\n来源与口径\n${sources.map(v=>v.startsWith('http')?v:`${v}\nhttps://github.com/Yifan2709/MonteCarlo-Atom-Tweezer/blob/4d5e200/${v}`).join('\n')}`);
 notes.push({n,title,caption,note,sources});return s;
}
function takeaway(s,t,color=C.blue){text(s,t,64,622,1150,48,28,color,true);}
function label(s,t,x,y,w=400,color=C.text,size=25){text(s,t,x,y,w,42,size,color,true);}
const rnd=v=>Number(Number(v).toPrecision(12));
function plot(s,series,o={}){
 const {x=57,y=178,w=1165,h=422,xlabel='',ylabel='',xmin,xmax,ymin,ymax,xstep,ystep,legend=true,percent=false,xf='0',yf='0',axisFont=21,legendFont=21,axes=true}=o;
 const first=series[0].x;
 const cc=s.charts.add('scatter',{position:{left:x,top:y,width:w,height:h},categories:first.map(v=>String(rnd(v))),
  series:series.map((v,i)=>({name:v.name||'',xValues:v.x.map(rnd),values:v.y.map(rnd),fill:v.color||cols[i%4],line:{fill:v.color||cols[i%4],width:v.points?0:(v.width||3)},marker:{symbol:v.points||v.markers?(v.symbol||'circle'):'none',size:v.size||7}})),
  scatterOptions:{style:'lineWithMarkers'},hasLegend:legend,legend:{position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:legendFont,fill:C.text}},
  xAxis:{visible:axes,title:{text:xlabel,textStyle:{typeface:FONT,fontSize:axisFont,fill:C.muted}},min:xmin,max:xmax,majorUnit:xstep,numberFormatCode:xf,position:'bottom',tickLabelPosition:axes?'low':'none',textStyle:{typeface:FONT,fontSize:axisFont,fill:C.muted},line:{fill:axes?C.gray:'none',width:1},majorGridlines:null},
  yAxis:{visible:axes,title:{text:ylabel,textStyle:{typeface:FONT,fontSize:axisFont,fill:C.muted}},min:ymin,max:ymax,majorUnit:ystep,numberFormatCode:percent?'0%':yf,position:'left',tickLabelPosition:axes?'low':'none',textStyle:{typeface:FONT,fontSize:axisFont,fill:C.muted},line:{fill:axes?C.gray:'none',width:1},majorGridlines:axes?{fill:C.grid,width:.7}:null},
  chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},plotAreaLine:{fill:'none',width:0}});
 applyPresentationChartFont(cc,{fontFamily:FONT});chartOwners.push(p.slides.items.indexOf(s)+1);return cc;
}
function bars(s,categories,series,o={}){
 const {x=64,y=185,w=1148,h=402,ylabel='',ymin=0,ymax,percent=false,legend=true,labels=true,yf='0',ystep}=o;
 const cc=s.charts.add('bar',{position:{left:x,top:y,width:w,height:h},categories,
  series:series.map((v,i)=>({name:v.name,values:v.y.map(rnd),fill:v.color||cols[i],valuesFormatCode:percent?'0.0%':'0.00'})),
  barOptions:{direction:'column',grouping:'clustered',gapWidth:100},hasLegend:legend,
  legend:{position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:22,fill:C.text}},
  dataLabels:{showValue:labels,position:'outEnd',textStyle:{typeface:FONT,fontSize:22,fill:C.text}},
  xAxis:{visible:true,position:'bottom',tickLabelPosition:'low',textStyle:{typeface:FONT,fontSize:23,fill:C.text},line:{fill:C.gray,width:1},majorGridlines:null},
  yAxis:{visible:true,title:{text:ylabel,textStyle:{typeface:FONT,fontSize:22,fill:C.muted}},min:ymin,max:ymax,majorUnit:ystep,numberFormatCode:percent?'0%':yf,position:'left',tickLabelPosition:'low',textStyle:{typeface:FONT,fontSize:21,fill:C.muted},line:{fill:C.gray,width:1},majorGridlines:{fill:C.grid,width:.7}},
  chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},plotAreaLine:{fill:'none',width:0}});
 applyPresentationChartFont(cc,{fontFamily:FONT});chartOwners.push(p.slides.items.indexOf(s)+1);return cc;
}
function table(s,values,widths,o={}){
 const {x=65,y=185,w=1150,h=388,size=27}=o;
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnTracks:widths.map(value=>({mode:'fr',value})),values});
 t.borders.assign({style:'solid',fill:C.grid,width:.7});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){const cell=t.getCell(r,c);cell.fill=r===0?'#EDF4F8':'#FFFFFF';cell.text.style={typeface:FONT,fontSize:size,color:r===0?C.blue:C.text,bold:r===0,verticalAlignment:'middle',autoFit:'none'};}
 tableOwners.push(p.slides.items.indexOf(s)+1);return t;
}
const keys=['manual_200us','manual_400us','manual_600us','ml_400us'];
const names=['手工 200 μs','手工 400 μs','手工 600 μs','优化 400 μs'];
const curve=(k,g='Joint model')=>D.curves.filter(r=>r.group===g&&r.waveform===k);
const points=k=>D.points.filter(r=>r.group==='Joint model'&&r.waveform===k);
const metric=(k,g='Joint model')=>D.metrics.find(r=>r.group===g&&r.waveform===k);
const ser=(name,x,y,color,extra={})=>({name,x,y,color,...extra});
function comparison(s,k,o={}){
 const max=k==='manual_200us'?30:60, rows=curve(k).filter(r=>+r.n<=max), ref=points(k);
 plot(s,[ser('模拟',rows.map(v=>+v.n),rows.map(v=>+v.survival),C.blue),ser('实验中心值',ref.map(v=>+v.n),ref.map(v=>+v.experiment),C.orange,{points:true})],
 {xlabel:'单程转移次数 n',ylabel:'存活比例',xmin:0,xmax:max,ymin:0,ymax:1,xstep:max===30?10:20,ystep:.2,percent:true,...o});
}

// 1: minimal cover, immediately tied to the observable.
{
 const s=page('光镊原子的 Monte Carlo','', '本汇报围绕短距离反复交接的存活实验。先提出要解释的曲线，再沿着一颗原子的实际操作顺序讲述，随后引入初态差异和随机噪声，最后回到实验结果。主报告24页，参数与推导放在6页附录。右图为归档最终共同模型的手工400微秒结果，说明模拟最后得到一条存活曲线。',['data','report']);
 text(s,'一颗颗算运动\n最后数剩下多少',65,214,530,174,53,C.blue,true);
 text(s,'图解汇报',69,446,500,50,29,C.muted);
 text(s,'Yifan2709',69,565,430,40,25,C.muted);
 const r=curve(keys[1]);plot(s,[ser('手工 400 μs 模拟示例',r.map(v=>+v.n),r.map(v=>+v.survival),C.blue)],{x:609,y:187,w:610,h:420,xlabel:'转移次数',ylabel:'存活比例',xmin:0,xmax:60,ymin:0,ymax:1,percent:true,xstep:20,ystep:.5,legend:true});
}
// 2: experimental question before any background.
{
 const s=page('反复交接以后，还能剩下多少原子','论文实验散点，不同操作得到不同的存活曲线',
 '这是整个模拟要回答的实验问题。横轴每一次计单程，接走和送回各计一次。手工200微秒的实验只到30次，其余到60次。图只显示论文提取的中心值，误差条上下限保留在来源文件中。现在先不解释为什么，只提出要由模型回答的两个差别：衰减形状为何不同，优化波形为何保留更多原子。',['points','paper']);
 plot(s,keys.map((k,i)=>{const r=points(k);return ser(names[i],r.map(v=>+v.n),r.map(v=>+v.experiment),cols[i],{points:true,symbol:['circle','square','triangle','diamond'][i],size:9});}),{xlabel:'单程转移次数 n',ylabel:'实验存活比例',xmin:0,xmax:60,ymin:0,ymax:1,percent:true,xstep:10,ystep:.2});
 takeaway(s,'模拟要解释曲线的差别，也要解释差别从哪里来');
}
// 3: the entire calculation as three quantitative views, no preliminary lecture.
{
 const s=page('一次虚拟实验的三个环节','抽一个初态，沿着势阱运动，再判断能否留下',
 '三幅图表示同一计算的逻辑顺序，不是三个独立实验。左图来自实际位点束缚初态采样的教学样本。中图是仓库历史一维600微秒放回轨迹。右图用同一个有限高斯阱表示束缚与未束缚总能量。真实计算在三维中逐颗重复这些步骤，之后将结局汇总成概率。此页用于建立总体方向，下面逐一展开。',['init','trajectory','noise']);
 label(s,'1  抽取初态',76,179);label(s,'2  计算运动',471,179);label(s,'3  判断去留',866,179);
 plot(s,[ser('',D.initial.x,D.initial.v,C.blue,{points:true,size:4})],{x:51,y:232,w:380,h:327,xlabel:'位置（μm）',ylabel:'速度（cm/s）',legend:false,xf:'0.0',axisFont:18});
 plot(s,[ser('',D.trajectory.time_us,D.trajectory.x_um,C.blue)],{x:447,y:232,w:379,h:327,xlabel:'时间（μs）',ylabel:'位置（μm）',legend:false,xmin:0,xmax:600,xstep:300,axisFont:18,yf:'0.0'});
 plot(s,[ser('',D.well.x,D.well.slm,C.blue),ser('',[-2.5,2.5],[0,0],C.gray),ser('',[-1,1],[-45,-45],C.green),ser('',[-1,1],[25,25],C.red)],{x:842,y:232,w:379,h:327,xlabel:'位置（μm）',ylabel:'能量 / kB（μK）',legend:false,ymin:-150,ymax:50,ystep:50,axisFont:18,xf:'0.0'});
 text(s,'绿线留下，红线丢失',924,565,292,35,21,C.muted);
 takeaway(s,'许多次虚拟实验的结果，组成一条存活曲线');
}
// 4: one full round.
{
 const s=page('一轮交接：接走、等待、送回','以手工 400 μs 操作为例，一轮共 900 μs',
 '上图是移动阱中心，下图是两束光的名义深度。接走400微秒，等待100微秒，送回400微秒。等待时SLM关闭，AOD独自束缚原子。放回是控制波形的时间反向，不是把原子的速度反向。30轮包含60次单程交接。深度以U/kB的正值表示阱深，非原子温度。边沿画为理想开关，未引入额外平滑。',['protocol']);
 label(s,'接走 400 μs',126,171,320);label(s,'等待 100 μs',543,171,300);label(s,'送回 400 μs',909,171,310);
 plot(s,[ser('移动阱中心',D.timeline.t,D.timeline.pos,C.blue)],{y:218,h:166,xlabel:'',ylabel:'位置（μm）',xmin:0,xmax:900,ymin:0,ymax:3,xstep:100,legend:false,axisFont:19});
 plot(s,[ser('静止阱',D.timeline.t,D.timeline.slm,C.gray),ser('移动阱',D.timeline.t,D.timeline.aod,C.orange)],{y:395,h:217,xlabel:'时间（μs）',ylabel:'阱深（μK）',xmin:0,xmax:900,ymin:0,ymax:300,xstep:100,ystep:100,axisFont:19,legendFont:19});
 takeaway(s,'30 轮连续进行，每轮都接着上一轮的末态');
}
// 5: physical potential snapshots.
{
 const s=page('接走原子时，势阱怎样变化','横向势能截面，名义束腰 1.11 μm',
 '蓝线为两束光相加后的势能。左侧先只有原点的静止阱，随后移动阱在原点加深，再带着深阱移到2.4微米。实际原子是否跟随，由当时的位置、速度和受力共同决定，势阱最低点不能直接当成原子的轨迹。四幅图是解析势阱截面，不含噪声及位点差异，也没有声称它们是实验图像。',['protocol','engine','params']);
 D.transfer.panels.forEach((v,i)=>{label(s,v.title,71+i*301,190,294,C.text,26);plot(s,[ser('总势能',D.transfer.x,v.total,C.blue)],{x:35+i*306,y:250,w:303,h:330,xlabel:'位置（μm）',ylabel:i===0?'U / kB（μK）':'',xmin:-2,xmax:4,ymin:-450,ymax:0,xstep:2,ystep:150,legend:false,axisFont:18});});
 takeaway(s,'原子能否跟上，要看它在变化的势阱里怎样运动');
}
// 6: one actual trajectory.
{
 const s=page('原子会跟随，也会在阱里振动','仓库一维放回示例：600 μs，初始温度 5 μK',
 '这幅图直接读取历史一维演示的热原子轨迹，不能当成最新三维19微开尔文标定样本。灰线是AOD中心，蓝线是原子。整体趋势相似，但原子在中心附近振动；中心停止后，振动仍继续。这幅图说明必须求原子的实际运动，不能只把光镊中心轨迹当作原子轨迹。牛顿方程的逐小步推进方法在附录用公式说明。',['trajectory','engine']);
 plot(s,[ser('光镊中心',D.trajectory.time_us,D.trajectory.aod_center_um,C.gray),ser('原子位置',D.trajectory.time_us,D.trajectory.x_um,C.blue)],{xlabel:'时间（μs）',ylabel:'位置（μm）',xmin:0,xmax:600,ymin:-.5,ymax:3,xstep:100,ystep:.5,yf:'0.0'});
 takeaway(s,'每一步：由势阱的坡度算力，再更新位置和速度');
}
// 7: initial ensemble, no lecture on MC definitions.
{
 const s=page('同样的操作，会遇到不同的初态','按仓库束缚初态方法抽样，T = 19 μK，仅作分布示意',
 '从同一温度下抽取180个可见样本。左图是横向位置与速度，右图比较横向与光轴方向的位置。沿光轴的束缚较弱，所以范围更宽。实际初态在各自的位点深度下抽取，并要求完整三维高斯势中的总能量小于零。初始化只做一次，后续轮次不重新抽取。这里用新教学种子，图中样本不是8192条定量结果的一个子集。',['init','params']);
 label(s,'位置相同，速度也可能不同',83,177,560);label(s,'沿光轴的分布更宽',701,177,510);
 plot(s,[ser('',D.initial.x,D.initial.v,C.blue,{points:true,size:5})],{x:53,y:229,w:578,h:364,xlabel:'横向位置 x（μm）',ylabel:'横向速度（cm/s）',legend:false,xf:'0.0'});
 plot(s,[ser('',D.initial.x,D.initial.z,C.orange,{points:true,size:5})],{x:658,y:229,w:565,h:364,xlabel:'横向位置 x（μm）',ylabel:'沿光轴位置 z（μm）',legend:false,xf:'0.0',yf:'0.0'});
 takeaway(s,'Monte Carlo 把这些初态差异逐颗保留下来');
}
// 8: stochastic energy histories.
{
 const s=page('随机噪声让原子的能量走向不同','谐振近似中的教学示例，不是本次标定轨迹',
 '示例用dE=ΓE dt+√Γ E dW的精确解，Γ=10/s，初始高于阱底能量40微开尔文，八个随机实现。橙线是理论平均40 exp(Γt)。教学模型对应弱强度噪声的相位平均能量扩散，用来说明平均上升可以伴随单颗能量上下变化。正式标定直接对完整高斯势中的力施加随机扰动，三维传播，不使用这条简化能量方程。',['physics','noise','savard']);
 const r=D.diffusion;plot(s,[...r.paths.map(y=>ser('单颗原子',r.t_ms,y,C.pale,{width:1.8})),ser('理论平均',r.t_ms,r.mean,C.orange,{width:4})],{xlabel:'时间（ms）',ylabel:'高于阱底的能量（μK）',xmin:0,xmax:30,ymin:0,legend:false});
 label(s,'粗橙线：平均能量',876,195,338,C.orange,24);
 takeaway(s,'光强和指向的涨落，像不断给原子小的随机推力');
}
// 9: the reason mean heating is not enough.
{
 const s=page('平均能量一样，损失也可以不同','两组构造分布：每组 100 颗，平均能量均为 80 μK',
 '这是专门构造的离散分布，不能当作生产模型的直方图。两组均有100颗、平均高于阱底能量80微开尔文。取固定阱深140微开尔文时，能量达到或超过140的原子在该静态模型中未束缚。窄分布没有，宽分布有20颗。正式多阱时变过程仍在指定目标阱检查能量，这里只展示分布尾部为什么重要。',['physics','noise']);
 bars(s,D.toyEnergy.x.map(v=>String(v)),[{name:'较集中',y:D.toyEnergy.narrow,color:C.blue},{name:'较分散',y:D.toyEnergy.wide,color:C.orange}],{ylabel:'原子数',labels:false,ymax:65});
 text(s,'高于阱底的能量（μK）',453,573,650,36,22,C.muted);
 takeaway(s,'达到 140 μK 的原子：较集中 0 颗，较分散 20 颗');
}
// 10: checkpoint meaning.
{
 const s=page('在指定时刻，判断原子是否仍被束缚','等待结束查移动阱，送回结束查静止阱',
 '这里用两口分开的高斯阱解释能量判据，线条为势能，水平短线是两个示例总能量。远处势能设为零，所以目标阱中的动能加势能小于零记为束缚。零以上记为丢失，这是该模型的检查点近似，不能说每一个微小时间步都立刻执行丢失判定。等待结束在AOD中判断，送回结束在SLM中判断。示例水平线不是计算得到的粒子轨迹。',['noise','protocol']);
 const x=D.well.x;label(s,'等待结束：移动阱',116,177,560);label(s,'送回结束：静止阱',750,177,460);
 for(let i=0;i<2;i++){plot(s,[ser('势能',x,i?D.well.slm:D.well.aod,C.blue),ser('零能量',[-2.5,2.5],[0,0],C.gray),ser('束缚',[-1.7,1.7],[-65,-65],C.green),ser('丢失',[-1.7,1.7],[35,35],C.red)],{x:55+i*610,y:234,w:562,h:347,xlabel:'相对阱中心的位置（μm）',ylabel:'能量 / kB（μK）',xmin:-2.5,xmax:2.5,ymin:-300,ymax:75,ystep:100,legend:false,axisFont:20,xf:'0.0'});}
 takeaway(s,'总能量 E < 0：留下；E ≥ 0：记为丢失');
}
// 11: continuous repeated history.
{
 const s=page('反复交接时，同一批原子一路往下走','手工 400 μs，最终共同模型，初始 8192 颗',
 '这是归档真实存活计数，不是每轮重新抽样。原子将位置和速度带入下一轮，上一轮造成的激发继续影响下一轮。丢失轨迹提前结束，不补原子，所以分母始终是8192。每一轮提供等待和放回两个检查点，即单程计数增加2。原子总能量不必单调上升，但在吸收式丢失规则下存活数只能不增。',['data','noise']);
 const r=curve(keys[1]);plot(s,[ser('仍被束缚',r.map(v=>+v.n),r.map(v=>+v.survivors),C.blue)],{xlabel:'单程转移次数 n',ylabel:'剩余原子数',xmin:0,xmax:60,ymin:0,ymax:8500,xstep:10,ystep:2000,legend:false});
 const last=r.at(-1);label(s,`第 60 次：${last.survivors} 颗`,878,230,332,C.blue,27);
 takeaway(s,'保留每颗原子的历史，丢失后不补充新原子');
}
// 12: converts count to probability at the actual end point.
{
 const r=curve(keys[1]),last=r.at(-1),s=page('计数以后，得到存活概率','手工 400 μs 的第 60 次转移',
 `初始8192颗，末点存活${last.survivors}颗，丢失${8192-(+last.survivors)}颗。相除得到${(+last.survival*100).toFixed(4)}%。这一步才把许多轨迹汇总成概率。8192来自两个独立新种子各4096颗。抽样误差在p约0.5时95%区间半宽约1.08个百分点，但不包含物理模型和时间步误差。`,['data','valid']);
 bars(s,['留下','丢失'],[{name:'原子数',y:[+last.survivors,8192-(+last.survivors)],color:C.blue}],{x:65,y:215,w:626,h:372,ylabel:'原子数',legend:false,labels:false,ymax:8192});
 text(s,`${last.survivors} / 8192`,765,248,451,80,43,C.text,true);
 text(s,`${(+last.survival*100).toFixed(1)}%`,765,353,451,110,74,C.blue,true);
 takeaway(s,'样本更多，比例更稳定；物理假设仍需要单独验证');
}
// 13: the actual competing controls.
{
 const s=page('四种操作的差别，放在同一套物理条件下','400 μs 手工与优化波形的直接比较',
 '两图由归档的手工三次位置曲线、48%二次加深阶段和PDF提取的优化控制点计算。优化节点用原来的三次样条，不额外夹紧导数或平滑。温度、束腰和噪声共用，200/400/600手工仅改变时长。避免一条曲线配一套温度或噪声。优化曲线控制的是阱中心与深度，原子轨迹仍需要另算。',['protocol','params','reports/movement_effects_v1_calibration_20260916/ml_vector_nodes.json']);
 label(s,'移动到哪里',132,178,500);label(s,'势阱有多深',764,178,460);
 plot(s,[ser('手工',D.wave.t,D.wave.manual_pos,C.blue),ser('优化',D.wave.t,D.wave.ml_pos,C.orange)],{x:55,y:230,w:567,h:358,xlabel:'时间（μs）',ylabel:'中心位置（μm）',xmin:0,xmax:400,ymin:0,ymax:2.6,xstep:100,yf:'0.0',axisFont:20});
 plot(s,[ser('手工',D.wave.t,D.wave.manual_depth,C.blue),ser('优化',D.wave.t,D.wave.ml_depth,C.orange)],{x:659,y:230,w:565,h:358,xlabel:'时间（μs）',ylabel:'阱深（μK）',xmin:0,xmax:400,ymin:0,ymax:320,xstep:100,ystep:100,axisFont:20});
 takeaway(s,'共用初温、束腰和噪声，只改变操作时长与波形');
}
// 14–16: results, with one message each.
{
 const s=page('快速交接的曲线明显接近实验','手工 200 μs，比较范围只到第 30 次',
 '最终共同模型相对旧基线的RMSE从19.11下降至2.76个百分点。但不能将全部改善归因于噪声，因为同时确认了手工控制曲线、谱系数和初态条件等问题。这里显示实验散点中心值，最大点偏差4.47个百分点，仍未全部进入±2个百分点参照带。曲线的30次以后没有实验散点，未展示或验收。',['data','points','metrics','report']);
 comparison(s,keys[0]);takeaway(s,'曲线误差从 19.11 降到 2.76 个百分点');
}
{
 const s=page('400 μs 的后期损失仍然不足','手工 400 μs，模型后期高估存活',
 '本组是剩余差距最明显的一组。n>20的平均残差约+9.99个百分点，总RMSE6.54，相比基线6.43没有整体改善。后期过多原子留在模型中，不能用另三条曲线的改善掩盖。图只显示实验中心值，不构造实验似然或统计通过判断。',['data','points','metrics','report']);
 comparison(s,keys[1]);takeaway(s,'后期平均高估约 10 个百分点',C.red);
}
{
 const s=page('600 μs 与优化波形也更接近实验','两组仍有一些存活率偏低的点',
 '左手工600微秒RMSE3.44个百分点，右优化400微秒3.19个百分点。二者的早期与后期平均残差均偏负，意味着模拟损失略多。图中是实验中心值，不显示误差条。缩放范围一致，便于比较两种操作的结果。详细每点残差见来源。',['data','points','metrics']);
 label(s,'手工 600 μs',139,174,450);label(s,'优化 400 μs',767,174,450);
 comparison(s,keys[2],{x:56,y:230,w:563,h:365,axisFont:20,legendFont:19});
 comparison(s,keys[3],{x:660,y:230,w:563,h:365,axisFont:20,legendFont:19});
 takeaway(s,'曲线误差分别为 3.44 和 3.19 个百分点');
}
// 17: evidence summary as bars, not a large table.
{
 const s=page('三条曲线改善，400 μs 仍有缺口','与全部对应实验散点比较，误差越小越接近',
 'RMSE是把各点偏差平方后平均再开方，单位是百分点。四组分别比较，不能只报联合平均。旧基线与最终模型用的参考散点及归一化相同。此处误差是模拟和实验的差距，不是模拟自身统计误差。±2个百分点仍是比较参照，而非由论文误差条推导的严格验收标准。',['metrics','report']);
 bars(s,names,[{name:'原基线',y:keys.map(k=>+metric(k,'Baseline').rmse_pp),color:C.gray},{name:'共同模型',y:keys.map(k=>+metric(k).rmse_pp),color:C.blue}],{ylabel:'曲线误差 RMSE（百分点）',ymax:22});
 takeaway(s,'共同模型有进展，但四条曲线尚未全部解释到目标范围',C.red);
}
// 18: a decisive mismatch rather than four more textual analyses.
{
 const s=page('最需要解释的是 400 与 600 μs 的差别','第 60 次转移，模型给出几乎相同的存活率',
 '直接从point_residuals.csv取末点：手工400实验约41.4%，600约57.3%。归档报告正文写41.8%与CSV有小差别，此处以原散点CSV为准。最终模拟400约52.0%，600约52.7%。这是模型对操作时长的相对反应仍不正确的具体表现，不只是整体曲线上下平移的问题。',['points','metrics','report']);
 const kk=[keys[1],keys[2]];bars(s,['手工 400 μs','手工 600 μs'],[{name:'实验中心值',y:kk.map(k=>+points(k).find(r=>+r.n===60).experiment),color:C.orange},{name:'模拟',y:kk.map(k=>+metric(k).S60),color:C.blue}],{ylabel:'存活比例',ymax:.7,percent:true});
 takeaway(s,'实验有明显差别，模型还没有把这个差别算出来');
}
// 19: a real local sensitivity test.
{
 const s=page('统一增大噪声，不能同时修正这些差距','归档中的共同参数对照，三组使用同一批初态',
 '局部选参对照为T19微开尔文、束腰1.11微米、N2048、同seed。RIN由5.75e-9增大到6.5e-9/Hz，400的RMSE6.86到5.39，600由2.71到5.83，优化400由3.37到4.37。它表明这一局部调节无法同时改善三组，不证明所有更完整物理模型都不可能解释。此页数值来自选参对照，样本池与上一页8192最终验证不同。',['report']);
 plot(s,[ser('手工 400 μs',[5.75,6.5],[6.86,5.39],C.blue,{markers:true,size:9}),ser('手工 600 μs',[5.75,6.5],[2.71,5.83],C.orange,{markers:true,size:9}),ser('优化 400 μs',[5.75,6.5],[3.37,4.37],C.green,{markers:true,size:9})],{xlabel:'AOD 强度噪声谱（10⁻⁹ /Hz）',ylabel:'曲线误差（百分点）',xmin:5.75,xmax:6.5,xstep:.75,ymin:0,ymax:8,xf:'0.00',ystep:2});
 takeaway(s,'400 μs 更接近了，另外两组却偏得更多');
}
// 20: phase is an extension, with native quantitative phase lines.
{
 const s=page('原子留下以后，还要问相位是否整齐','相位示意：轨迹差异带来不同的相位积累，橙线表示平均值',
 '一个箭头代表一次轨迹中两个内部态的相对相位。这里用单位圆上的数值线段表示8个相位。左侧全部同向，中间共同转90度，右侧在圆周等间隔分散。共同相位偏移会改变条纹位置，但不降低C；彼此分散使平均向量变短。不同原子不需要相互干涉，统计可以是同一个位点的重复实验。图为精确构造示例，非本次计算输出。',['phase','noise']);
 const circle=Array.from({length:81},(_,i)=>2*Math.PI*i/80);
 for(let j=0;j<3;j++){
  label(s,['相位对齐','一起转过 90°','相位散开'][j],88+j*407,182,366);
  const phases=Array.from({length:8},(_,i)=>j===0?0:j===1?Math.PI/2:2*Math.PI*i/8);
  const ss=[ser('',circle.map(Math.cos),circle.map(Math.sin),C.grid,{width:1}),...phases.map(a=>ser('',[0,Math.cos(a)],[0,Math.sin(a)],C.blue,{width:2.3}))];
  const mx=phases.reduce((v,a)=>v+Math.cos(a),0)/8,my=phases.reduce((v,a)=>v+Math.sin(a),0)/8;
  ss.push(ser('',[0,mx],[0,my],C.orange,{width:5}));
  plot(s,ss,{x:80+j*409,y:235,w:315,h:315,xmin:-1.1,xmax:1.1,ymin:-1.1,ymax:1.1,legend:false,axes:false,xlabel:'',ylabel:''});
  label(s,`相干度 C = ${Math.hypot(mx,my).toFixed(0)}`,122+j*407,550,330,C.blue,26);
 }
 takeaway(s,'存活率高，不代表内部状态的相位仍然一致');
}
// 21: echo as an actual function plot.
{
 const s=page('回波可以抵消一部分相位差','示例：三条线代表三种稳定的相位转速',
 '用三种固定失谐比较自由相位与理想回波后的残余相位。中间施加瞬时π脉冲，相位误差在后半程以相反符号累积，末端抵消。图取相对时间0到1，不是具体脉冲实验数据。运输过程中的光移随位置和光强改变，所以真实回波未必完全抵消。仓库有限脉冲计算还保留约20.3微秒脉冲期间的原子运动，需要传播内部量子态。',['phase','engine','protocol']);
 label(s,'不加回波',149,179,470);label(s,'中间加一次理想回波',742,179,465);
 [false,true].forEach((echo,j)=>plot(s,[.6,1,1.4].map((k,i)=>ser('',D.echo.x,(echo?D.echo.echo:D.echo.free).map(v=>k*v),[C.pale,C.blue,C.orange][i])),{x:56+j*608,y:232,w:565,h:362,xlabel:'相对时间 t / T',ylabel:'相位（rad）',xmin:0,xmax:1,ymin:0,ymax:1.5,xstep:.5,ystep:.5,xf:'0.0',yf:'0.0',legend:false,axisFont:21}));
 takeaway(s,'前后两段的相位误差能抵消多少，取决于它们是否相似');
}
// 22: explicit boundary, short and factual.
{
 const s=page('这次结果回答到了哪里','存活率和相干性分开说明',
 '左列是最新随机扩散存活标定实际输出。右列是仓库已有相干性扩展能够计算的量，但它没有与本次最终噪声条件一起完成新的联合验证。内部态模型的η=1.3e-4仍是简化假设。使用半经典路径，没有求完整运动波函数，也没有自旋相关力对运动的反馈。读者不能把两批结果合成一次已经完成的联合证明。',['noise','cli','phase','report']);
 table(s,[['本次存活标定','仓库中的相干性扩展'],['三维运动与随机扰动','沿轨迹累积差分光移'],['统计留下和丢失','统计相位分散'],['与四条实验曲线比较','比较回波与有限脉冲'],['尚有 400 μs 的明显差距','尚未与本次噪声联合验证']],[1,1],{h:364,size:28});
 takeaway(s,'同一套轨迹可以连接两类问题，当前证据仍需分别看');
}
// 23: the useful next physical measurements.
{
 const s=page('下一步需要哪些实验输入','先区分原因，再增加模型细节',
 '初态分布影响最初损失，AOD/SLM各自噪声谱影响随时间的能量扩散，真实位置和阱深控制关系影响操作激发。三类输入在原归档报告中均被列为优先独立测量。这里没有声称某一个是已证实的残差来源，也不建议直接给400微秒曲线添加独立损失倍率。',['report','params']);
 table(s,[['需要测量','要分清的问题'],['转移前的温度与能量分布','初态是否比假设更热或更偏'],['两条光路各自的噪声谱','随机加热到底来自哪里'],['位置、阱深和开关的实际波形','400 与 600 μs 为何响应不同']],[1.1,1.6],{h:340,size:29});
 takeaway(s,'更独立的物理输入，才能缩小剩余原因的范围');
}
// 24: final argument on one figure, not a word-heavy recap.
{
 const s=page('Monte Carlo 把单颗运动变成整体概率','当前共同模型已经改善比较，仍保留明确的未解释差距',
 '结束时回到第一页提出的问题。每条轨迹由初态、时变势和噪声决定，检查点将它转成留下或丢失，多次抽样给出存活曲线。最后不能只说拟合成功：三条明显改善，400微秒后期仍高估。右侧用四组真实实验末点与模拟末点对照，200取n30其余取n60。相干性则是对这些运动路径继续问内部态相位的问题，不在本次存活验收中。',['data','points','metrics','report']);
 text(s,'逐颗算运动',70,237,430,62,39,C.blue,true);
 text(s,'汇总成概率',70,342,430,62,39,C.blue,true);
 text(s,'用差距检验物理假设',70,451,490,78,32,C.text,true);
 const chosen=keys.map(k=>k===keys[0]?30:60);
 bars(s,['200 μs\nn=30','400 μs\nn=60','600 μs\nn=60','优化\nn=60'],[{name:'实验',y:keys.map((k,i)=>+points(k).find(r=>+r.n===chosen[i]).experiment),color:C.orange},{name:'模拟',y:keys.map((k,i)=>+curve(k).find(r=>+r.n===chosen[i]).survival),color:C.blue}],{x:582,y:215,w:641,h:379,ylabel:'末点存活比例',ymax:.8,percent:true,labels:false});
}
// Appendices: retain precision without interrupting the report.
// 25
{
 const s=page('附录：共同模型参数','四条曲线共用，同一组参数不能视为唯一测量结果',
 '最终温度19微开尔文，名义束腰1.11微米，AOD单边相对强度噪声谱5.75e-9/Hz。SLM短时间近似安静，声透镜在短途近似下关闭。束腰已经在探索范围边界。两批新种子各4096原子，机械步长0.0125微秒，噪声间隔0.05微秒。位点差异和对准抖动按原归档固定，不单独按曲线调节。',['params','report']);
 table(s,[['量','共同取值'],['初温 / 名义束腰','19 μK / 1.11 μm'],['静止阱 / 移动阱深','140 / 280 μK'],['AOD 强度噪声谱','5.75 × 10⁻⁹ /Hz'],['机械步长 / 噪声间隔','0.0125 / 0.05 μs'],['最终验证样本','4096 × 2 = 8192']],[1.1,1.5],{h:388,size:27});
 takeaway(s,'SLM 在短时间内近似安静，短途比较关闭声透镜');
}
// 26
{
 const s=page('附录：为什么使用三维势阱','同一束光，横向和光轴方向的约束不同',
 '曲线为束腰1.17微米、140微开尔文、波长1061纳米的说明性截面，与最终1.11微米参数分开。横向高斯衰减，沿光轴按瑞利长度变化。轴向较宽且曲率较低，所以束缚较弱。仅计算沿搬运方向的一维轨迹会漏掉横向和轴向能量交换或逃逸。空间运动采用经典粒子近似，不求完整运动波函数。',['engine','init']);
 plot(s,[ser('横向',D.axial.x,D.axial.radial,C.blue),ser('沿光轴',D.axial.x,D.axial.axial,C.orange)],{xlabel:'距焦点的位置（μm）',ylabel:'U / kB（μK）',xmin:-8,xmax:8,ymin:-150,ymax:0,xstep:2,ystep:50});
 takeaway(s,'初始位置与速度、受力和逃逸判断都在三个方向上计算');
}
// 27
{
 const s=page('附录：四个公式对应四个步骤','公式只补充物理含义',
 '第一行势能决定力，第二行牛顿方程决定轨迹，第三行目标阱能量给检查点束缚判据，第四行将保留轨迹计数变成存活比例。数值上用半步速度、一步位置、再半步速度反复推进，随机增量另按指定噪声间隔加入。每个时间步越小，越接近连续运动，但必须检查步长收敛。',['engine','noise']);
 table(s,[['物理步骤','表达式'],['势阱给出力','F = −∇U'],['力改变运动','m dv/dt = F，dr/dt = v'],['检查目标阱中的能量','E = ½m|v|² + U目标'],['统计存活','Sn = N留下 / N初始']],[1,1.5],{h:366,size:30});
 takeaway(s,'随机性来自初态、位点条件和有物理含义的噪声');
}
// 28
{
 const s=page('附录：样本误差与时间步误差','最终结果的精度来自两类不同的检查',
 '左图是p=0.5时二项独立抽样的近似95%半宽1.96 sqrt[p(1-p)/N]，不是实验误差也不包含模型误差。右图读取VALIDATION.md中的噪声间隔从0.05减为0.025微秒时最大存活率差，分别1.123、1.050、1.440、0.708个百分点。机械步长固定0.0125。噪声路径不做Brownian bridge耦合，所以是分布比较，不能要求单颗轨迹相同。机械步减半另有约0.68至2.12个百分点的差异，不能据此宣称亚百分点精度。',['valid']);
 label(s,'增加样本',127,180,475);label(s,'噪声间隔减半',761,180,455);
 bars(s,['256','4096','8192'],[{name:'半宽',y:[256,4096,8192].map(n=>1.96*Math.sqrt(.25/n)*100),color:C.blue}],{x:55,y:230,w:565,h:355,ylabel:'95% 区间半宽（百分点）',legend:false,ymax:7,labels:false});
 bars(s,['200','400','600','优化'],[{name:'最大差',y:[1.123,1.050,1.440,.708],color:C.orange}],{x:658,y:230,w:564,h:355,ylabel:'最大存活率差（百分点）',legend:false,ymax:2,ystep:.5,yf:'0.0',labels:false});
 takeaway(s,'增加样本不能消除时间步误差，更不能消除模型误差');
}
// 29
{
 const s=page('附录：长距离运输中的声透镜','光束形状变化会改变轴向束缚',
 '图为单束像散高斯光的轴向势能，两个正交方向的焦点从重合到分别偏移正负zR、正负2zR。它不是两束独立AOD光的势叠加。快速扫频可以造成柱面透镜效应，降低或改变束缚。最终短交接标定采用关闭透镜近似，长距离扩展属于另一协议，不能混入Fig6d短交接曲线比较。图参数为教学用280微开尔文、1.17微米、1055纳米。',['engine','report','paper']);
 plot(s,D.lensing.curves.map((v,i)=>ser(['焦点重合','焦点偏移 ±zR','焦点偏移 ±2zR'][i],D.lensing.x,v,cols[i])),{xlabel:'沿光轴的位置 z（μm）',ylabel:'U / kB（μK）',xmin:-12,xmax:12,ymin:-300,ymax:0,xstep:4,ystep:100});
 takeaway(s,'短交接与长距离运输，应分别确认对应的光束形状');
}
// 30
{
 const s=page('附录：来源与数据范围','讲解备注保留逐页出处',
 '所有定量实验对照直接读取归档CSV，没有重新运行全套标定。新增初态散点、势阱截面、噪声能量示意和相位示意在对应页面明确标为教学说明。历史单轨迹标明一维、600微秒和5微开尔文。数据来源冻结在4d5e200读取基准，原报告5bae8b0保留。源数据和构建脚本随报告提交。',['guide','report','data','points','metrics','paper','savard']);
 label(s,'思路参考',72,193,1080,C.blue,28);text(s,'仓库《蒙特卡洛模拟·大白话完整说明》',72,241,1110,58,29);
 label(s,'定量证据',72,328,1080,C.blue,28);text(s,'2026-09-16 共同标定归档与原始散点',72,376,1110,58,29);
 label(s,'物理依据',72,463,1080,C.blue,28);text(s,'Manetsch 等的光镊阵列论文；Savard 等的噪声加热理论',72,511,1110,88,27);
}

if(p.slides.items.length!==30)throw Error(`Expected 30, got ${p.slides.items.length}`);
await fs.writeFile(path.join(build,'presentation.json'),JSON.stringify(p.toProto()));
await fs.writeFile(path.join(here,'图解版讲解备注.md'),'# Monte Carlo 图解版讲解备注\n\n主报告24页，附录6页。\n\n'+notes.map(v=>`## ${v.n}. ${v.title}\n\n${v.caption}\n\n${v.note}\n\n来源：\n${v.sources.map(r=>`- ${r}`).join('\n')}\n`).join('\n'));
await fs.mkdir(path.join(here,'output'),{recursive:true});
const candidate=path.join(build,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
console.log(`Exported ${p.slides.items.length} slides.`);
const samples=(process.env.SAMPLE_SLIDES||'').split(',').filter(Boolean).map(Number);
for(const n of samples){const blob=await p.export({slide:p.slides.items[n-1],format:'png',scale:1});await fs.writeFile(path.join(build,`preview-${n}.png`),new Uint8Array(await blob.arrayBuffer()));}
if(process.env.DRAFT_ONLY!=='1'){
 const out=process.env.FINAL_PPTX||path.join(here,'output','Monte_Carlo_光镊原子_图解汇报版_修订.pptx');
 const tables=[...new Set(tableOwners)],charts=[...new Set(chartOwners)];
 const result=await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:out,pythonExecutable:python,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tables.flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:30,requiredNativeChartOwnerSlides:charts,requiredNativeTableOwnerSlides:tables,fontPolicy:{basis:'design',families:[FONT]},materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,
 receiptPath:path.join(build,`${path.basename(out)}.validation.json`)});
 console.log(JSON.stringify({path:out,sha256:result.finalSha256}));
}
