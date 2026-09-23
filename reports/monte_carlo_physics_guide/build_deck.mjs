import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';
const here=path.dirname(fileURLToPath(import.meta.url)), root=path.resolve(here,'../..');
const build=path.join(here,'.build');
const runtime=process.env.RUNTIME_NODE_MODULES || '/Users/imok/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.RUNTIME_NODE_MODULES=runtime;
const require=createRequire(path.join(runtime,'_resolve.cjs'));
const {Presentation,PresentationFile}=await import(pathToFileURL(require.resolve('@oai/artifact-tool')).href);
const skill=process.env.PRESENTATIONS_SKILL_DIR || '/Users/imok/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations';
const {applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const python=process.env.RUNTIME_PYTHON || '/Users/imok/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const D=JSON.parse(await fs.readFile(path.join(build,'data.json'),'utf8'));
const p=Presentation.create({slideSize:{width:1280,height:720}});
const C={text:'#172B3A',muted:'#526575',blue:'#17679B',orange:'#CB6C28',green:'#278071',red:'#AE3E3D',gray:'#8C98A3',grid:'#DEE5EA',bg:'#FFFFFF'};
const FONT='Hiragino Sans GB';
const notes=[],chartOwners=[],tableOwners=[];
const SOURCE={guide:'docs/monte_carlo_plain_guide.md', init:'simulation/level2_joint_transfer/src/continuous_transfer/initialization.py', protocol:'simulation/level2_joint_transfer/src/continuous_transfer/protocol.py', engine:'simulation/level2_joint_transfer/src/continuous_transfer/engine.py', noise:'simulation/level2_joint_transfer/src/continuous_transfer/noisy_survival.py', cli:'simulation/level2_joint_transfer/src/continuous_transfer/cli.py', report:'reports/movement_effects_v1_calibration_20260916/REPORT.md', physics:'reports/movement_effects_v1_calibration_20260916/PHYSICS.md', params:'reports/movement_effects_v1_calibration_20260916/PARAMETERS.md', valid:'reports/movement_effects_v1_calibration_20260916/VALIDATION.md', data:'reports/movement_effects_v1_calibration_20260916/survival_curves.csv', points:'reports/movement_effects_v1_calibration_20260916/point_residuals.csv', metrics:'reports/movement_effects_v1_calibration_20260916/metrics.csv', paper:'https://arxiv.org/html/2403.12021v4', savard:'https://doi.org/10.1103/PhysRevA.56.R1095'};
function tx(s,str,x,y,w,h,size=27,color=C.text,bold=false){
 const v=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 v.text=str;v.text.style={typeface:FONT,fontSize:size,color,bold,autoFit:'none',verticalAlignment:'top'};return v;
}
function page(title,sub,note,sources=['guide']){
 const s=p.slides.add();s.background.fill=C.bg;const n=p.slides.items.length;
 tx(s,title,60,38,1160,65,42,C.text,true);
 if(sub)tx(s,sub,62,113,1154,48,23,C.muted);
 tx(s,String(n).padStart(2,'0'),1170,675,60,25,16,C.gray);
 const refs=sources.map(k=>SOURCE[k]||k);
 s.speakerNotes.textFrame.setText(`讲解\n${note}\n\n来源与口径\n${refs.map(x=>x.startsWith('http')?x:`${x}\nhttps://github.com/Yifan2709/MonteCarlo-Atom-Tweezer/blob/4d5e200/${x}`).join('\n')}`);
 notes.push({n,title,sub,note,refs});return s;
}
function bottom(s,text,color=C.blue){tx(s,text,64,603,1150,65,26,color,true);}
function rows(s,items,{x=66,y=183,w=1145,gap=100,head=27,body=26}={}){
 items.forEach((it,i)=>{tx(s,it[0],x,y+i*gap,w,42,head,C.blue,true);tx(s,it[1],x,y+42+i*gap,w,70,body);});
}
function prose(title,sub,items,conclusion,note,refs=['guide']){
 const s=page(title,sub,note,refs);rows(s,items,{gap:items.length===4?98:130});if(conclusion)bottom(s,conclusion);return s;
}
function table(s,values,widths,{x=64,y=186,w=1152,h=365,size=25}={}){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnTracks:widths.map(value=>({mode:'fr',value})),values});
 t.borders.assign({style:'solid',fill:C.grid,width:.7});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){
  const cell=t.getCell(r,c);cell.fill=r===0?'#ECF3F7':'#FFFFFF';
  cell.text.style={typeface:FONT,fontSize:size,color:r===0?C.blue:C.text,bold:r===0,verticalAlignment:'middle',autoFit:'none'};
 }
 tableOwners.push(p.slides.items.indexOf(s)+1);return t;
}
function chart(s,xs,series,{x=55,y=184,w=800,h=391,xlabel='',ylabel='',xmin,xmax,ymin,ymax,xstep,ystep,format,percent=false,legend=true,kind='scatter',markers=false}={}){
 const colors=[C.blue,C.orange,C.green,C.gray,C.red];
 const rnd=v=>Number(Number(v).toPrecision(12));
 const xformat=format||((Math.max(...xs)-Math.min(...xs))<=16?'0.0':'0');
 const allY=series.flatMap(v=>v.y), yformat=(Math.max(...allY)-Math.min(...allY))<=5?'0.0':'0';
 const cc=s.charts.add(kind,{position:{left:x,top:y,width:w,height:h},
  categories:xs.map(v=>String(rnd(v))),
  series:series.map((v,i)=>({name:v.name,xValues:(v.x||xs).map(rnd),values:v.y.map(rnd),fill:v.color||colors[i],line:{fill:v.color||colors[i],width:v.pointsOnly?0:2.6},marker:{symbol:v.pointsOnly||markers?'circle':'none',size:v.pointsOnly?7:4},...v.extra})),
  scatterOptions:{style:'lineWithMarkers'},lineOptions:{smooth:false},barOptions:{direction:'column',grouping:'clustered',gapWidth:65},
  hasLegend:legend,legend:{position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:19,fill:C.text}},
  xAxis:{visible:true,title:{text:xlabel,textStyle:{typeface:FONT,fontSize:20,fill:C.muted}},min:xmin,max:xmax,majorUnit:xstep,numberFormatCode:xformat,position:'bottom',tickLabelPosition:'low',textStyle:{typeface:FONT,fontSize:18,fill:C.muted},line:{fill:C.gray,width:1},majorGridlines:null},
  yAxis:{visible:true,title:{text:ylabel,textStyle:{typeface:FONT,fontSize:20,fill:C.muted}},min:ymin,max:ymax,majorUnit:ystep,numberFormatCode:percent?'0%':yformat,position:'left',tickLabelPosition:'low',textStyle:{typeface:FONT,fontSize:18,fill:C.muted},line:{fill:C.gray,width:1},majorGridlines:{fill:C.grid,width:.7}},
  chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},plotAreaLine:{fill:'none',width:0}});
 applyPresentationChartFont(cc,{fontFamily:FONT});chartOwners.push(p.slides.items.indexOf(s)+1);return cc;
}
function side(s,items){rows(s,items.map(v=>[v[0],v[1].replace(/。$/,'')]),{x:870,y:192,w:346,gap:132,head:25,body:24});}
function graph(title,sub,xs,series,opt,items,conclusion,note,refs){
 const s=page(title,sub,note,refs);chart(s,xs,series,opt);side(s,items);bottom(s,conclusion);return s;
}
const labels={'manual_200us':'手工 200 μs','manual_400us':'手工 400 μs','manual_600us':'手工 600 μs','ml_400us':'优化 400 μs'};
const keys=Object.keys(labels);
const cv=k=>D.curves.filter(r=>r.group==='Joint model'&&r.waveform===k);
const pts=k=>D.points.filter(r=>r.group==='Joint model'&&r.waveform===k);
const mt=k=>D.metrics.find(r=>r.group==='Joint model'&&r.waveform===k);
// 01
{
 const s=page('光镊中的原子','', '本报告从原子所经历的物理过程出发。先讲光镊怎样约束原子，怎样随机准备一批初态，再讲拾取、等待和放回，最后解释存活统计和相干性。文档用于组织思路，具体公式、检查时点和结果以仓库实现与归档证据为准。',['guide','report']);
 tx(s,'Monte Carlo 模拟在做什么',64,210,1150,90,57,C.blue,true);
 tx(s,'从一颗原子的运动，到一批原子的存活概率',66,333,1120,75,32);
 tx(s,'物理原理与逐步实现',67,459,1050,50,29,C.muted);
 tx(s,'Yifan2709    2026 年 9 月',68,588,1080,50,24,C.muted);
}
// 02
prose('汇报的四个问题','按原子经历的顺序展开',[
 ['原子一开始是什么状态','光镊形状、温度和初始位置/速度。'],
 ['转移时，原子怎样运动','两口势阱怎样变化，力怎样改变原子的运动。'],
 ['随机性怎样影响结果','位点差异、初态差异和运动中的噪声。'],
 ['最后怎样判断与比较','存活、统计误差、相干性，以及与论文的差距。']],
 '贯穿全篇：每条轨迹都保留自己的历史',
 '可以把一次模拟看成一次虚拟实验。不是直接给原子一个成功概率，而是先计算轨迹，再判断成功。结果可信到什么程度，会在统计、步长检验和最后的模型边界中回答。');
// 03
prose('要回答的实验问题','同样的操作反复进行，最后还剩多少原子',[
 ['实验动作','用可移动的光镊，把原子从静止光镊中接走，再送回。'],
 ['主要观测','做完 n 次单程转移后，仍留在目标阱里的原子比例。'],
 ['模拟的任务','计算每颗原子的运动与能量，再把“留下”或“丢失”汇总。']],
 '更快、更多次的转移，是否仍能可靠地留住原子？',
 '这里对应论文 Fig. 6d 上图的短距离转移存活实验。n 数的是单程交接，拾取和放回各计一次。真实实验靠最后的成像统计，模拟用目标势阱中的束缚能近似这一观测。',['protocol','physics','paper']);
// 04
{
 const s=page('两口光镊，各自做什么','以铯-133 为例，把光镊理解为用激光形成的势阱',
 'SLM 是形成静止光镊阵列的器件，AOD 用来控制可移动光镊。这里不需要先学习器件工作原理，只需记住静止阱和移动阱这两个角色。表中的参数为讲解基准，最新共同模型把名义束腰改为1.11微米。',['guide','params']);
 table(s,[['物理量','静止光镊 SLM','移动光镊 AOD'],['作用','原子的起点和终点','接走、搬运、送回'],['名义阱深','140 μK','280 μK'],['激光波长','1061 nm','1055 nm'],['讲解用束腰半径','1.17 μm','1.17 μm'],['短距离交接','中心保持在原点','最多离开 2.4 μm']],[1.15,1.4,1.4],{h:365});
 bottom(s,'阱深用 μK 表示能量尺度，原子温度另行指定');
}
// 05
graph('势阱像一口有限深的碗','横轴是位置，纵轴是势能除以玻尔兹曼常数',D.well.x,[{name:'静止阱 140 μK',y:D.well.slm},{name:'移动阱 280 μK',y:D.well.aod}],{xlabel:'离阱中心的位置（μm）',ylabel:'U / kB（μK）',ymin:-300,ymax:20},[
 ['中心最低','原子离开中心，需要增加势能。'],['远处趋近零','阱的深度有限，足够高的总能量允许逃逸。'],['更深的阱','通常能容纳更大的运动能量。']],
 'U(x) = −D exp[−2(x − c)² / w²]',
 '曲线按仓库采用的高斯势公式计算，是物理示意，不是新一次完整模拟。D为正的阱深，c为中心，w为光强降到中心1/e²时的半径。远处U取零，所以阱内势能为负。不要把负能量理解为负温度。',['guide','engine']);
// 06
graph('原子受到的力来自势阱的坡度','静止阱的一维截面',D.well.x,[{name:'静止阱的力',y:D.well.force}],{xlabel:'离阱中心的位置（μm）',ylabel:'力（10⁻²¹ N）',legend:false},[
 ['左侧','力指向右边，把原子拉回中心。'],['右侧','力向左，把原子拉回中心。'],['阱中心','力为零，但原子可能仍有速度。']],
 'F = −∂U/∂x，随后用 F = ma 求运动变化',
 '零受力不等于零速度。原子经过阱底时速度通常最大，它会继续向另一侧运动，再被回复力拉回来。距离很远时，高斯势梯度又变小，已经跑远的原子不会被无限强的力拉回。',['engine','guide']);
// 07
prose('有限温度意味着原子一直在动','同一温度的原子，也不会处于同一个位置和速度',[
 ['位置有分布','有的在阱底附近，有的更靠近阱边。'],
 ['速度有分布','有的往左，有的往右，也有快慢差异。'],
 ['开始时的振动阶段不同','同样的光镊操作，会碰上不同的振动阶段。']],
 '只算一颗静止在阱底的原子，会漏掉热运动带来的差别',
 '热分布规定的是许多次制备的统计规律。对某一颗原子，位置与速度决定了当下的振动状态。即使没有任何后续噪声，不同初态也可以出现不同的最终能量与束缚结果。',['init','guide']);
// 08
prose('阱深、温度和原子能量','三个量相关，但含义不同',[
 ['阱深 D','从阱底到远处零势能的高度，本身是能量。'],
 ['初始温度 T','决定一批原子初始位置和速度的分布宽度。'],
 ['总能量 E','某颗原子当下的动能加势能，E = ½m|v|² + U。']],
 '三维谐振近似下，平均“高于阱底的能量”约为 3kBT',
 '三维谐振子有三个动能和三个势能的二次项，每项平均kBT/2，所以高于阱底的平均能量是3kBT。讲解基准15 μK对应45 kB·μK，但有限高斯阱、拒绝采样和非平衡运动都会使这一近似有偏离。末态平均总能量的变化也不能直接等同于温升。',['init','physics']);
// 09
graph('三维中，沿光轴的束缚更弱','同一束光的横向和轴向势能截面',D.axial.x,[{name:'横向 x',y:D.axial.radial},{name:'沿光轴 z',y:D.axial.axial}],{xlabel:'离焦点的距离（μm）',ylabel:'U / kB（μK）',ymin:-150,ymax:0},[
 ['横向较窄','原子受到较强的回复力。'],['轴向较宽','相同温度下，轴向分布更宽。'],['三维的必要性','原子可以沿 x、y 或 z 方向跑出阱。']],
 '只看搬运方向的一维模型，会漏掉其他方向的运动',
 '图使用140 μK、束腰1.17 μm、波长1061 nm。瑞利长度约4.05 μm。径向频率约25.46 kHz，轴向约5.20 kHz。这里横轴为了比较共用距离，蓝线远离中心后迅速接近零。',['init','engine']);
// 10
prose('这份模拟采用的物理近似','原子的空间运动与内部量子态分别描述',[
 ['空间运动','把原子当作有位置和速度的粒子，沿经典轨迹运动。'],
 ['内部量子态','需要讨论相干性时，再沿轨迹计算两能级状态的变化。'],
 ['适用范围','适合研究热运动、交接激发和损失趋势，未求完整运动波函数。']],
 '量子内部态可以保留，空间运动仍采用经典近似',
 '以SLM径向25.46 kHz为例，hf/kB约1.22 μK，小于讲解温度15 μK，因此经典热运动有一定依据。但这不是低温基态和隧穿问题的完整量子描述。本文不把它称作量子跳跃蒙特卡洛，也不引入多原子碰撞。',['engine','init']);
// 11
prose('Monte Carlo 怎样给出概率','按已知分布反复抽取样本，再逐颗计算',[
 ['准备很多次虚拟实验','每次得到一个不同的初态和位点条件。'],
 ['按相同物理规律向前演化','初态一旦给定，力和运动方程决定接下来的变化。'],
 ['把每次结局变成统计','留下记 1，丢失记 0，平均值就是模拟存活概率。']],
 '物理规律决定每条轨迹，随机抽样把它们组成一批实验',
 '随机性不是任意编造的成功或失败。分布需要来自热制备假设、装置差异或噪声模型。增加样本数可以降低抽样波动，但不能修正错误的势阱形状、错误波形或错误噪声假设。',['guide','noise']);
// 12
{
 const s=page('随机数分别代表什么','不同来源的随机性，有不同的更新时间',
 '四类随机量要区分。位点差异在一条轨迹内固定，初态只准备一次，段间对准抖动每轮每段变化，连续噪声才是在传播中不断产生的小速度扰动。旧的固定平均功率加热本身没有新的随机数。',['init','engine','noise','params']);
 table(s,[['随机量','物理含义','什么时候抽取'],['位点深度、束腰、偏差','每个光镊略有不同','每条轨迹开始前一次'],['位置和速度','热运动的初始状态','初态准备时一次'],['段间对准抖动','交接阶段的位置漂移','每轮、每段一次'],['运动中的噪声','光强和指向的快涨落','随机噪声模型中反复抽取']],[1.4,1.5,1.4],{h:355,size:25});
 bottom(s,'固定随机种子，便于重复同一组比较');
}
// 13
{
 const s=page('先固定一组可比较的物理条件','讲解基准与最终标定条件分开使用',
 '左列是文档采用的理解基准。右列是2026-09-16联合标定报告锁定的最终模型。四种操作共用右列参数，只允许操作时间和波形不同。19 μK与噪声幅值是模型条件下的估计，并不是独立测量。',['guide','params','report']);
 table(s,[['参数','讲解基准','最终共同模型'],['初始温度','15 μK','19 μK'],['SLM / AOD 名义阱深','140 / 280 μK','140 / 280 μK'],['名义束腰半径','1.17 μm','1.11 μm'],['加热描述','固定平均功率','随机力引起能量扩散'],['验证原子数','单批 4096','两批各 4096，合计 8192']],[1.5,1.45,1.8],{h:375,size:25});
 bottom(s,'先用基准数值理解步骤，再看共同模型的实际结果');
}
// 14
graph('先给每个位点一个实际阱深','示例：4096 个独立样本，SLM 深度相对标准差 11.4%',D.depthHist.x,[{name:'样本数',y:D.depthHist.y}],{kind:'bar',xlabel:'SLM 阱深（μK）',ylabel:'样本数',legend:false},[
 ['名义值 140 μK','各位点围绕它上下分布。'],['浅阱的影响','同一温度下，位置分布更宽，逃逸余量更小。'],['全程固定','同一颗原子始终属于同一个位点样本。']],
 '这些是独立模拟样本，不等于实物装置恰好有 4096 个位点',
 '图为本报告使用固定种子产生的教学抽样，非论文实测直方图。11.4%来自论文全阵列的深度散布，转移子阵列没有独立给出。AOD深度和束腰2%的散布、50 nm静态对准误差属于模型假设。',['params','init']);
// 15
{
 const s=page('初态的分布宽度来自温度和阱深','三维讲解基准：15 μK，140 μK SLM，束腰 1.17 μm',
 '标准差表示分布宽度，并不是每个原子都在这个位置。轴向频率低，所以z分布宽。速度的三个分量具有同样标准差。这里的数字是名义深度下的值，实际采样逐位点使用自己的SLM深度。',['init']);
 table(s,[['初始量','标准差','物理解释'],['x、y 位置','0.191 μm','横向束缚较紧'],['z 位置','0.938 μm','沿光轴束缚较弱'],['vx、vy、vz','3.06 cm/s','由温度和铯原子质量决定'],['径向 / 轴向振动频率','25.46 / 5.20 kHz','这里是 f，角频率 ω = 2πf']],[1.2,1.45,2.0],{h:338,size:25});
 bottom(s,'位置宽度 σx = √[kBT / (mωx²)]，速度宽度 σv = √(kBT / m)');
}
// 16
graph('一批初态在位置和速度上散开','从实际位点内准备的 4096 个样本中，展示前 180 个',D.sample.x,[{name:'一颗点代表一条初始轨迹',y:D.sample.v,pointsOnly:true}],{xlabel:'初始 x（μm）',ylabel:'初始 vx（cm/s）',legend:false,xmin:-.7,xmax:.7,ymin:-10,ymax:10},[
 ['同一位置','也可能具有不同的速度。'],['同一速度','也可能正处于振动的不同位置。'],['相同的转移','会把这批初态带向不同的结局。']],
 '三维时，每颗原子需要 3 个位置分量和 3 个速度分量',
 '本图直接调用仓库的实际位点束缚采样方法，再取前180个样本作展示。图中只画x与vx，其余四个分量仍存在。散点并非在模拟过程中彼此作用，原子之间的碰撞和相互作用没有加入。',['init']);
// 17
prose('初态必须已经被实际势阱束缚','先抽位点，再在这个位点里准备原子',[
 ['计算分布宽度','使用该位点的实际深度，而不是一律使用名义 140 μK。'],
 ['抽位置和速度','先用阱底近似的热分布提出一个初态。'],
 ['检查真实高斯势中的能量','只有 E < 0 才保留，否则在同一位点重新抽取。']],
 '重抽状态，保留位点，避免只剩下容易留住原子的深阱',
 '若先按名义阱准备，再把浅阱误差加进去，会混入初始就不束缚的原子。若直接删除不合格位点，又会偏向深阱。现有做法在同一位点重抽，保证每个位点一个已束缚初态。该方法是谐振分布提议加束缚条件，仍不等于有限深高斯阱的严格热平衡分布。',['init','report']);
// 18
prose('原子带着上一轮的状态继续运动','每一轮都保留位置、速度和已经累积的变化',[
 ['留下的原子','把本轮末态直接作为下一轮的起点。'],
 ['已经丢失的原子','后续各次统计一直记为丢失，不再补入新的原子。'],
 ['持续累积的影响','运动激发、噪声和位点差异会改变后续损失风险。']],
 '如果每轮重新抽热初态，就会抹掉重复转移的历史',
 '初始拒绝采样只发生在t=0，不能在后续每轮重新做。幸存者的分布经过筛选也会变化，因此第20次转移的原子群与第一次之前已经不同。这是存活曲线通常不能简单写成固定单次概率的原因之一。',['engine','noise']);
// 19
{
 const s=page('一轮短转移的时间安排','拾取一次、等待一次、放回一次，合计两次单程转移',
 '匹配Fig.6d的协议只做短程交接。一轮时长为2T+100 μs，三种手工操作的30轮总时长分别为15、27、39 ms。优化400 μs与手工400 μs总时长相同。n=60表示30轮。等待段纳入奇数检查点对应的拾取观察。',['protocol','params']);
 table(s,[['单程操作时间 T','拾取','等待','放回','30 轮总时长'],['200 μs','200 μs','100 μs','200 μs','15 ms'],['400 μs','400 μs','100 μs','400 μs','27 ms'],['600 μs','600 μs','100 μs','600 μs','39 ms']],[1.45,1,1,1,1.5],{h:302,size:25});
 bottom(s,'n = 60 是 60 次单程交接，也就是 30 个来回');
}
// 20
prose('拾取时，移动阱先变深，再离开','手工操作的前 48% 与后 52% 分别承担不同作用',[
 ['前 48%：在原位接住','AOD 中心与 SLM 重合，AOD 深度按时间的平方升高。'],
 ['后 52%：搬离原位','AOD 深度保持 280 μK，中心移动 2.4 μm。'],
 ['SLM 一直开着','原子在两束光形成的总势中运动，不预先指定属于哪一口阱。']],
 '原子能否跟随移动阱，要由它实际受到的力决定',
 '对400 μs操作，升深持续192 μs，移动持续208 μs。原子不会因为“拾取开始”这个名字就自动跳到AOD。两阱重合时总势加深，随后分开，原子能否随AOD离开取决于当时的位置、速度和总势形状。',['protocol','physics']);
// 21
{
 const s=page('手工 400 μs 的两条控制曲线','一条规定阱中心在哪里，一条规定阱有多深',
 '图为按目前仓库已纠正的论文手工波形计算的控制量。深度为二次升高，位置为单段3s²−2s³。它在移动段两端速度为零，但加速度不连续，不能把它与另一种端点加速度为零的分段S曲线混淆。',['protocol','physics']);
 chart(s,D.wave.t,[{name:'AOD 中心',y:D.wave.manual_pos}],{x:50,w:570,h:350,xlabel:'时间（μs）',ylabel:'中心位置（μm）',legend:false,ymin:0,ymax:2.6});
 chart(s,D.wave.t,[{name:'AOD 深度',y:D.wave.manual_depth,color:C.orange}],{x:660,w:570,h:350,xlabel:'时间（μs）',ylabel:'阱深（μK）',legend:false,ymin:0,ymax:300});
 tx(s,'升深用时 192 μs，移动用时 208 μs',66,560,1130,42,27);
 bottom(s,'移动段：c / 2.4 μm = 3s² − 2s³，其中 s 从 0 变到 1');
}
// 22
graph('拾取过程中，总势阱不断改变','一维截面示意，忽略位点偏差和噪声',D.snapshots.x,[{name:'起初，仅 SLM',y:D.snapshots.curves[0]},{name:'两阱重合',y:D.snapshots.curves[1]},{name:'分开 1.2 μm',y:D.snapshots.curves[2]},{name:'分开 2.4 μm',y:D.snapshots.curves[3]}],{xlabel:'位置 x（μm）',ylabel:'总势能 / kB（μK）',ymin:-450,ymax:0},[
 ['重合时','两束光相加，使中心势阱更深。'],['逐渐分开','势阱形状和最低点都在改变。'],['原子需要跟上','只移动光斑，不保证原子同步移动。']],
 '总势 U = USLM + UAOD，总力也由两束光相加得到',
 '四条线代表不同时刻的势能快照，并非同一时刻的四条轨迹。模型按两束光的势能直接相加，忽略光场干涉结构。这里以一维截面帮助理解，正式传播仍使用三维势。',['engine','protocol']);
// 23
prose('等待段检验能否独立留在 AOD','拾取结束后，关闭 SLM 并等待 100 μs',[
 ['关光的瞬间','原子位置和速度连续，势能因为少了一束光而改变。'],
 ['等待期间','原子只受 AOD 的束缚，继续振动并受到所选噪声影响。'],
 ['等待结束','用“动能 + AOD 势能”检查是否还处于束缚状态。']],
 '开关一束光会改变能量，不能凭空改变原子的位置或速度',
 '这是一个容易弄错的物理边界。关闭SLM并不是给原子额外设定一个速度，而是改变哈密顿量中的势能项。实现将开关前后的势能差计入切换功。在等待结束发现的损失，很多在关闭SLM后立刻就已经是正能。',['protocol','noise','physics']);
// 24
prose('放回时，光镊操作按时间倒放','先把 AOD 移回原处，再降低 AOD 深度',[
 ['SLM 重新打开','它再次参与总势，为最后接住原子提供束缚。'],
 ['AOD 反向操作','中心回到原点，深度最终降到零。'],
 ['原子继续正常向前运动','它保留当前速度和此前的振动，噪声也不会倒着消失。']],
 '倒放控制曲线，不会自动把原子带回原来的运动状态',
 '严格动力学时间反演还需要反转动量并反转完整动力学条件。本模拟只反转控制随时间的函数，不反转原子速度，也不逆转噪声。一个来回结束后可以有剩余激发，这些激发继续影响下一轮。',['protocol','engine']);
// 25
{
 const s=page('优化 400 μs 同时调整位置和深度','对照论文图中提取的控制曲线',
 'ML代表论文的优化波形。本报告用归档PDF控制点回放：13个内部点加起止点，共15个节点、14个时间间隔，再作普通三次插值。它不是原始仪器驱动记录，因此仍有数字化与校准的不确定性。本报告没有重新训练或优化波形。',['physics','reports/movement_effects_v1_calibration_20260916/ml_vector_nodes.json']);
 chart(s,D.wave.t,[{name:'手工',y:D.wave.manual_pos},{name:'论文优化',y:D.wave.ml_pos}],{x:50,w:570,h:350,xlabel:'时间（μs）',ylabel:'中心位置（μm）',ymin:0,ymax:2.6});
 chart(s,D.wave.t,[{name:'手工',y:D.wave.manual_depth},{name:'论文优化',y:D.wave.ml_depth}],{x:660,w:570,h:350,xlabel:'时间（μs）',ylabel:'阱深（μK）',ymin:0,ymax:300});
 tx(s,'更合适的操作时序，可以改变原子被激发的程度',66,560,1150,42,27);
 bottom(s,'比较操作效果时，其他物理参数和初始样本保持一致');
}
// 26
prose('每个极短时间步里发生什么','运动方程拆成小步，近似连续时间中的真实运动',[
 ['1  算当前位置受到的力','用当前时刻的两口势阱，得到三个方向的加速度。'],
 ['2  速度先更新半步，位置再走一步','原子带着已有速度前进，力逐渐改变它的运动。'],
 ['3  在新位置重新算力，再更新半步速度','这样能跟上势阱的移动和深度变化。'],
 ['4  加入所选噪声，再继续下一步','需要量子态时，还沿这段运动计算内部态变化。']],
 '力决定加速度，速度决定位置随时间的变化',
 '这对应velocity-Verlet算法，但口头报告不必强调算法名称。机械步长0.0125 μs小于阱振动周期，仍需通过进一步减半来检验误差。最新存活模型的噪声步长另为0.05 μs，机械与噪声时间间隔分别验证。',['engine','noise']);
// 27
graph('一颗原子怎样跟随移动阱','仓库的一维 600 μs 放回示例，初始温度 5 μK',D.trajectory.time_us,[{name:'AOD 中心',y:D.trajectory.aod_center_um,color:C.gray},{name:'原子位置',y:D.trajectory.x_um,color:C.blue}],{xlabel:'时间（μs）',ylabel:'位置（μm）',xmin:0,xmax:600},[
 ['总体跟随','原子从 2.4 μm 附近移向原点。'],['叠加振动','实际原子位置在阱中心附近来回摆动。'],['末态仍会运动','中心停下后，原子仍会振动。']],
 '光镊中心的轨迹与原子自己的轨迹，是两个不同的量',
 '这张图直接来自仓库Level1历史演示的thermal_example_trajectory.csv，并非最新三维校准结果。只做等间隔抽点以减小PPT数据体积。位置轨迹正确保留原始单位。该示例用于看懂动力学，不用于推断最新存活概率。',['simulation/level1_transfer_1d/outputs/level1_transfer_1d/demo/thermal_example_trajectory.csv']);
// 28
prose('能量改变，要分清来源','运动中的势阱会对原子做功',[
 ['移动中心、改变深度','即使完全没有噪声，时变势也可能给原子增加激发。'],
 ['光噪声和光子反冲','额外改变动量，使许多原子的能量分布发生变化。'],
 ['光镊静止且没有噪声时','总能量应近似守恒，可以检查计算是否准确。']],
 '能量变化 = 操作做功 + 噪声做功 + 数值误差',
 'dE/dt=∂U/∂t适用于无额外噪声的显式时变势。阱深变浅会把势能整体抬高，所以Ef−Ei不能一概称为加热量。描述末态激发应说明相对哪一个阱底，并把外部做功与噪声能量账目分开。',['physics','simulation/level2_joint_transfer/src/level2_joint_transfer/work_energy.py']);
// 29
prose('更慢的操作也有代价','机械激发与噪声累积，会随操作时间以不同方式变化',[
 ['操作太快','势阱来不及平稳地带动原子，结束后可能留下较大振动。'],
 ['操作更慢','机械激发可能减小，但受到噪声作用的时间更长。'],
 ['重复很多次','残余振动、噪声和逃逸筛选共同改变后续分布。']],
 '不存在只由“总时间越长越好”决定的统一答案',
 '在谐振近似的随动坐标q=x−c中，有q̈+ω²q=−c̈。光镊加速度相当于驱动，而激发取决于加速度时间形状与原子振动的配合。有限高斯势和双阱交接更复杂，所以需要实际传播轨迹。',['physics','protocol']);
// 30
{
 const s=page('三种加热来源，各自改变什么','这里的噪声幅值需要实验测量或明确的模型假设',
 '强度涨落改变势阱深度与回复力，指向涨落移动势阱位置，散射光子带来反冲。模型中的白噪声幅值并非装置已测得的谱。段间10 nm偏移属于另一种慢漂移近似，不应与连续指向噪声谱混为一谈。',['physics','noise','savard']);
 table(s,[['来源','直观图像','主要影响'],['光强涨落','碗忽深忽浅','改变回复力，带来能量扩散'],['光斑指向涨落','碗的位置轻微抖动','不断轻推原子'],['光子散射','原子吸收和发射光子','获得很小的反冲动量']],[1.2,1.4,1.8],{h:305,size:26});
 bottom(s,'噪声的物理类型与大小，共同决定模拟中的加热');
}
// 31
prose('最容易理解的加热近似','每一步给原子增加一个固定的平均能量',[
 ['规定平均功率 P','经过小段时间 Δt，增加动能 PΔt。'],
 ['把速度的长度放大一点','三维中沿原速度方向缩放，使增加的动能正好等于 PΔt。'],
 ['保留了什么','可以跟踪平均注入能量，便于理解和做基础诊断。']],
 '这个近似没有产生每一步的新随机扰动',
 '旧连续模型采用vnew=v·sqrt(1+2PΔt/(m|v|²))，速度为零时单独处理。这保证平均注能账目清楚，但它不能描述噪声使不同原子能量分布变宽。最新标定不再用一个Tref来设定整段固定功率。',['engine','physics']);
// 32
prose('更新的做法：随机力推动原子','相同的噪声强度，也会给不同原子不同的瞬时作用',[
 ['光强噪声改变力','随机速度变化与原子当下感受到的那束光的力有关。'],
 ['局部运动决定响应','位置、方向、阱深和波形变了，噪声的实际作用也会变。'],
 ['单颗原子有升有降','某一步可能增能，也可能减能，整批平均仍表现为加热。']],
 '最新存活标定采用随机传播，直接得到能量分布的变化',
 '对单边每Hz的强度噪声S，速度增量为(F/m)sqrt(S h/2) ξ，ξ为标准正态随机数。同一光路的标量涨落同时作用于xyz三个力分量，不是三方向各自独立的强度噪声。两束光按独立光路处理。噪声步长h不是拟合出来的物理带宽。',['noise','physics']);
// 33
{
 const s=page('平均能量一样，逃逸比例仍可不同','100 个原子的简化示意，能量都以“高于阱底”计',
 '这两组柱形数据是明确构造的教学例子，都有100个原子，平均激发能都是80 kB·μK。窄分布没有达到140，宽分布有20个达到或超过140。在深度140 μK的简化阱中，这意味着不同的逃逸比例。图不是模拟或实验结果。',['physics']);
 chart(s,D.toyEnergy.x,[{name:'较窄的分布',y:D.toyEnergy.narrow},{name:'较宽的分布',y:D.toyEnergy.wide}],{kind:'bar',x:60,y:183,w:1155,h:365,xlabel:'高于阱底的能量 / kB（μK）',ylabel:'原子数',ymin:0,ymax:65});
 tx(s,'两组平均都是 80 μK；达到 140 μK 的分别是 0 个与 20 个',65,561,1150,42,26);
 bottom(s,'存活率特别关心高能尾部，只有平均加热功率通常不够');
}
// 34
prose('噪声为什么会与振动频率有关','先用阱底附近的谐振运动理解',[
 ['光强噪声','在约两倍振动频率附近的起伏，更容易改变振动能量。'],
 ['指向噪声','在振动频率附近推动阱中心，更容易激发原子。'],
 ['真正传播时','使用完整高斯势的局部力，不把整个过程固定成一个谐振频率。']],
 '噪声谱告诉我们：不同频率附近的涨落有多强',
 '单边每Hz谱定义下，谐振极限参量平均增长率Γ=π² f² S_I(2f)，指向平均功率P=mω⁴Sx(f)/4。旧代码的相应系数分别少4倍与2倍，最新报告已纠正。公式的系数依赖谱定义，不能混用每Hz与每角频率、单边与双边谱。',['physics','savard']);
// 35
{
 const s=page('什么时候判断原子已经丢失','匹配 Fig. 6d 的短交接协议，每轮检查两次',
 '需要精确区分传播与观测：力在每个微小时间步更新，但短协议并不是每一步都按能量丢弃原子。等待末检查AOD，放回末检查SLM。E≥0是一种吸收边界近似，正能但还没跑远的原子在真实系统中可能再捕获。',['engine','noise','physics']);
 table(s,[['检查时点','用哪一口阱的势能','判定'],['等待 100 μs 结束','只用 AOD','E = K + UAOD < 0 为束缚'],['放回结束','只用 SLM','E = K + USLM < 0 为束缚'],['一旦 E ≥ 0','记录首次失去束缚','后续统计持续记为丢失']],[1.25,1.3,2],{h:299,size:25});
 bottom(s,'计算每一小步的运动，按明确检查点记录存活');
}
// 36
prose('发现损失的阶段，不一定是损失原因','等待末发现丢失，可能早在拾取结束时就已失去束缚',[
 ['原子先被激发','拾取时受到的机械激发，会带入等待段。'],
 ['SLM 关闭后显露出来','失去原先的附加束缚，原子在 AOD 中可能已经 E ≥ 0。'],
 ['最终标定中的记录','等待末损失中，约 98.6%–99.9% 在等待开始时已为正能。']],
 '分析原因，需要同时看检查前后的能量和光镊切换',
 '百分比来自最终两个独立验证种子的归档诊断：手工200/400/600和优化400分别约98.6/99.5/99.9/99.4%。它说明“pickup记录为零损失”不能证明拾取无激发。这是该模型的诊断结果，并不是直接测得的实验机制。',['report']);
// 37
prose('存活曲线就是每个检查点数人头','每颗初始原子的权重相同',[
 ['定义存活率','Sn = 第 n 次转移后仍在阱中的原子数 ÷ 初始原子数。'],
 ['一个数值例子','初始 4096 个，若还有 3072 个束缚，存活率就是 75%。'],
 ['整条轨迹一起保留','同一颗原子在多个时刻的存活记录彼此有关。']],
 'S0 = 1，采用丢失后不返回的规则时，Sn 不会上升',
 '例子中的3072只是为了说明计数，不是仓库某条曲线的测量值。实际统计使用alive数组逐检查点求和。绘制置信区间或重抽样时，需要保留一个原子的完整历史，不能把61个检查点当成相互独立的新实验。',['noise','cli']);
// 38
{
 const s=page('样本数决定抽样波动有多大','以下按独立原子样本、存活率约 50% 估计',
 '二项分布标准误差约为sqrt(S(1−S)/N)。表中95%半宽采用1.96倍标准误差作直观估计，正式曲线归档使用Wilson区间。真实实验可能有同一光路的共同噪声，独立样本计算不能替代那种相关误差。',['cli','valid']);
 table(s,[['原子样本数 N','标准误差约为','95% 区间半宽约为'],['256','3.13 个百分点','6.13 个百分点'],['4096','0.78 个百分点','1.53 个百分点'],['8192','0.55 个百分点','1.08 个百分点']],[1,1.5,1.7],{h:303,size:27});
 bottom(s,'样本数增加 4 倍，抽样误差大约减半');
}
// 39
prose('怎样公平比较四种操作','把差异集中在操作时间和控制曲线上',[
 ['共用物理条件','温度、名义阱深、位点差异与噪声参数保持相同。'],
 ['共用初始样本','相同的初始原子分别经历四种操作，减小无关抽样差异。'],
 ['选定参数后另取新样本','最终验证使用两个未参与选参的新种子，各 4096 个原子。']],
 '四条曲线共用一组参数，不能分别调温度来追每条曲线',
 '参数探索和独立验证分开。当前最终验证种子为26091621与26091622，合计8192个原子，不能把早期基准的23003当作最终随机模型的验证池。相同随机数有助于比较，但不同操作持续时间不同，并不等于每颗原子受到完全相同的物理历史。',['report','valid']);
// 40
{
 const s=page('最新共同模型采用的条件','用于随后四条存活曲线，属于模型条件下的标定',
 '最终模型假设短时间窗口SLM光路近似安静，AOD存在附加白噪声。这个光路区分受论文独立静态寿命约束，但噪声幅值仍没有独立装置测量。短途近似关闭声透镜。1.11 μm已经到所探索束腰范围边界，不能当作精确测量。',['params','report','physics']);
 table(s,[['共同条件','取值或做法'],['初始温度 / 名义束腰','19 μK / 1.11 μm'],['SLM / AOD 阱深','140 / 280 μK'],['AOD 相对强度噪声谱','5.75 × 10⁻⁹ /Hz，单边谱'],['SLM 技术噪声','在这段短时间内近似安静'],['机械 / 噪声时间间隔','0.0125 / 0.05 μs'],['传播内容','三维运动与随机噪声，专门统计存活']],[1.55,2.5],{h:385,size:25});
 bottom(s,'这组存活结果没有同时给出新的量子态保真度');
}
// 41
{
 const s=page('重复转移后的模拟存活曲线','最终共同模型，合计 8192 个独立原子样本',
 '本图来自归档最终共同模型。200 μs仍计算到n=60，但论文该组散点只到n=30，所以后半段是模型预测，不用它做区间外实验验收。其他三条与实验比较到n=60。线条用于连接离散检查点。',['data','report']);
 chart(s,cv(keys[0]).map(r=>+r.n),keys.map((k,i)=>({name:labels[k],y:cv(k).map(r=>+r.survival),color:[C.red,C.blue,C.orange,C.green][i]})),{x:61,w:1152,h:375,xlabel:'单程转移次数 n',ylabel:'存活比例 Sn',xmin:0,xmax:60,ymin:0,ymax:1,percent:true,xstep:10,ystep:.2});
 tx(s,'各条曲线的差异来自波形、作用时间和演化后的原子分布',66,562,1150,45,26);
 bottom(s,'200 μs 的 n > 30 部分仅为模拟预测，没有对应实验散点');
}
// 42--45
for(const [j,k] of keys.entries()){
 const r=cv(k), ref=pts(k), max=k==='manual_200us'?30:60, rr=r.filter(v=>+v.n<=max), m=mt(k);
 const comments=[
  ['快速操作','曲线误差从 19.11 降到 2.76 个百分点。','剩余差异','前期略低估、末段略高估存活。'],
  ['最明显的缺口','后期仍高估存活，平均偏差约 10 个百分点。','整体误差','误差 6.54 个百分点，\n未比原基线改善。'],
  ['较慢操作','曲线误差为 3.44 个百分点。','剩余差异','多数实验点处，模拟存活偏低。'],
  ['优化操作','曲线误差为 3.19 个百分点。','剩余差异','模拟仍低估一部分实验存活。']
 ][j];
 const s=page(`${labels[k]}：模拟与实验`, '圆点是论文实验散点，蓝线是最终共同模型',
  `数据直接取自归档survival_curves与point_residuals。RMSE为实验散点处的均方根差，以百分点计。${k==='manual_200us'?'该组有7个实验散点，只比较n=2到30。':'该组有10个实验散点，比较至n=60。'} 实验误差条的精确定义及相关性未明，本页用散点中心值作比较。当前原生图仅显示中心值，误差条上下限完整保留于源CSV，不能把这张图当成统计显著性检验。`,['data','points','metrics','report']);
 chart(s,rr.map(v=>+v.n),[{name:'共同模型',y:rr.map(v=>+v.survival)},{name:'论文散点',x:ref.map(v=>+v.n),y:ref.map(v=>+v.experiment),pointsOnly:true,color:C.orange}],{xlabel:'单程转移次数 n',ylabel:'存活比例 Sn',xmin:0,xmax:max,ymin:0,ymax:1,percent:true,xstep:10,ystep:.2});
 side(s,[[comments[0],comments[1]],[comments[2],comments[3]],['最大点偏差',`${(+m.max_abs_pp).toFixed(2)} 个百分点。`]]);
 bottom(s,'本页显示实验散点中心值；尚未全部进入 ±2 个百分点参照带');
}
// 46
{
 const s=page('共同模型改善了哪些曲线','与全部对应实验散点比较，误差单位为百分点',
 '联合RMSE采用每条曲线先计算均方误差，再对四条等权平均并开方，不是简单平均四个RMSE。最终联合值约4.26 pp，而基线约14.31 pp。200/600/优化400明显改善，400并没有整体改善。±2 pp为预先规定的比较参照，不是由论文误差条推导的严格统计验收线。',['metrics','report']);
 table(s,[['曲线','原基线 RMSE','最终 RMSE','最大散点偏差'],...keys.map(k=>[labels[k],(+D.metrics.find(r=>r.group==='Baseline'&&r.waveform===k).rmse_pp).toFixed(2),(+mt(k).rmse_pp).toFixed(2),(+mt(k).max_abs_pp).toFixed(2)])],[1.6,1.4,1.3,1.55],{h:335,size:26});
 tx(s,'RMSE：把各点误差平方后求平均，再开方，用一个数概括整条曲线的差距',66,550,1150,43,22,C.muted);
 bottom(s,'四条曲线仍未全部达到 ±2 个百分点的比较目标',C.red);
}
// 47
{
 const s=page('最明确的剩余差异在 400 与 600 μs','到第 60 次转移时，实验与模型对时长的反应不同',
 '这是一项结构性差异。模型给400和600几乎相同的末点，而实验600明显高于400。继续统一加大噪声会同时改变四条曲线，报告中的局部对照显示改善400会恶化600和优化400，因此不能给400单独加任意损失倍率来掩盖。',['points','metrics','report']);
 table(s,[['n = 60','论文实验','最终模型'],...['manual_400us','manual_600us'].map(k=>[labels[k],`${(100*+pts(k).find(r=>+r.n===60).experiment).toFixed(1)}%`,`${(100*+mt(k).S60).toFixed(1)}%`])],[1.6,1.3,1.3],{h:260,size:29});
 tx(s,'模型仍没有解释清楚：为何两种手工操作在实验中相差这么多',66,511,1150,66,28);
 bottom(s,'参数接近数据，不代表已经找到了唯一的真实损失机制');
}
// 48
prose('数值结果要经得起进一步细分时间','抽样误差之外，还有时间离散带来的误差',[
 ['机械时间步减半','四条曲线的最大存活差约 0.68–2.12 个百分点。'],
 ['600 μs 再做更细检查','进一步减半后，最大差约 0.98 个百分点。'],
 ['噪声时间间隔也减半','四条曲线的最大差约 0.71–1.44 个百分点。']],
 '已有检查支持总体趋势，但不支持声称小于 1 个百分点的精度',
 '机械步长0.0125到0.00625 μs，600组另到0.003125 μs。噪声间隔0.05到0.025 μs。报告还显示不少逐原子的逃逸标签发生变化，集合存活曲线接近不意味着每条轨迹都已收敛。参数锁定后做这些检查，不能挑最贴实验的步长。',['valid','report']);
// 49
prose('长距离运输是另一个物理问题','扩展协议在拾取之后加入远距离去程和回程',[
 ['短交接','每轮只搬离 2.4 μm，再送回，直接用于 Fig. 6d 对照。'],
 ['长距离扩展','另外去程 375 μm、远端停 54 μs、回程 375 μm。'],
 ['连续保留状态','长途之前已有的振动和内部态变化，也要带入长途过程。']],
 '扩展协议每轮约 754.8 μm，不能直接当作短交接实验来验收',
 '扩展的去回程各1450 μs，400 μs交接条件下每轮共3854 μs，30轮约115.62 ms，名义总路程22.644 mm。扩展模型有额外的长程段末检查，并非每一小步都执行逃逸吸收。最新Fig.6d随机噪声标定只覆盖matched短协议。',['protocol','noise']);
// 50
graph('快速扫描可能改变沿光轴的势阱','两横向方向的焦点错开时，轴向束缚会变形',D.lensing.x,[{name:'焦点重合',y:D.lensing.curves[0]},{name:'各偏移 1 个瑞利长度',y:D.lensing.curves[1]},{name:'各偏移 2 个瑞利长度',y:D.lensing.curves[2]}],{xlabel:'沿光轴位置 z（μm）',ylabel:'AOD 势能 / kB（μK）',ymin:-300,ymax:0},[
 ['声透镜效应','扫描时，声波调制可能产生附加聚焦作用。'],['焦点错开','束缚会变弱，也可能出现两个轴向最低点。'],['参数需标定','示意图不代表本装置已测得的焦点偏移。']],
 '它描述一束散光光镊的两横向焦点，不是两束独立 AOD 光',
 '曲线使用仓库a、b形式的高斯势，设s1=+δ、s2=−δ并取三档偏移作教学演示。实际偏移由两个扫描方向的速度与有效尺度vs决定。vs不是AOD中的声速。短程最终校准关闭了这一效应，所以不能把该图当作短交接剩余差异的既定解释。',['engine','physics']);
// 51
prose('原子留下来以后，还要看内部态','原子的位置状态与存储的量子信息分别受到影响',[
 ['两个内部态','可以把它们理解为量子比特的两个基础状态。'],
 ['光会改变两态的相对能量','原子走过不同位置，感受到的光强历史不同。'],
 ['相位逐渐不同','不同原子的内部“指针”可能转过不同角度。']],
 '原子没有丢失，也可能已经积累了不同的相位误差',
 '这里讨论铯钟态的差分光移。半经典近似用δ(t)=ηU(r(t),t)/ħ描述失谐，η=1.3×10⁻⁴在仓库中作为假设参数。重点是原子运动提供了光强历史，内部态沿这个历史演化。该模型忽略内部态对机械轨迹的反馈。',['engine','guide']);
// 52
graph('相位分散会降低整批原子的相干性','示意：相位服从宽度为 σφ 的高斯分布',D.coherence.x,[{name:'相干对比度',y:D.coherence.y}],{xlabel:'相位分布宽度 σφ（rad）',ylabel:'相干对比度 C',ymin:0,ymax:1,legend:false},[
 ['相位接近','平均后仍有清楚的共同方向。'],['相位分散','不同方向相互抵消，集体信号变弱。'],['共同相位','所有原子一起转同一角度，并不会降低这个 C。']],
 'C = |平均 exp(iφ)|，高斯示意中 C = exp(−σφ²/2)',
 '图为解析教学曲线，非仓库某一组实际相干结果。相干对比度取平均复相位的模。共同确定相位不降低C，却可能偏离指定目标门，造成保真度误差，所以C和平均门保真度不能互换。计算时还要说明是否只对幸存原子平均。',['cli','engine']);
// 53
graph('回波脉冲可以抵消缓慢相位误差','理想示意：一颗原子受到恒定失谐，半途施加一次翻转',D.echo.x,[{name:'不加回波',y:D.echo.free},{name:'理想回波',y:D.echo.echo}],{xlabel:'总时间的比例',ylabel:'归一化累积相位',ymin:0,ymax:1},[
 ['前半程','相位误差向一个方向积累。'],['翻转以后','对同样的慢失谐，后半程的相位贡献反号。'],['末尾抵消','恒定误差可以被理想回波消除。']],
 'XY4 使用四次、沿不同方向的翻转脉冲，减小慢变化误差',
 '图是理想单回波的toggling-frame相位，不是实际Bloch矢量角度。它帮助说明误差抵消。实际XY4还涉及X/Y轴交替以及有限脉冲期间的失谐，不是简单把真实运动或已经散射的光子倒回去。',['simulation/level2_joint_transfer/src/level5_finite_pulse_irb/dd_finite_pulses.py','engine']);
// 54
prose('真实微波脉冲需要有限时间','脉冲作用期间，原子仍在运动',[
 ['典型驱动频率','本模型采用 Ω / 2π = 24.611 kHz。'],
 ['一次 π 翻转','持续时间约为 20.3 μs，不能总当作瞬间完成。'],
 ['一起计算两种变化','运动改变光移，同时微波改变内部量子态。']],
 '有限脉冲效果需要完整内部态演化，不能只靠一个累积相位',
 'τπ=π/Ω=1/(2×24611 Hz)。统一旧连续模型用有限驱动计算SU(2)演化；无脉冲时可先累积z相位，但脉冲期间失谐与驱动不对易，需要共同传播。最新随机存活专用模型明确关闭DD和内部态，因此本节讲仓库的相干扩展能力，而非声称8192原子标定已同步给出相干结果。',['engine','protocol','noise']);
// 55
{
 const s=page('存活率与量子态质量，分别报告','首先明确分母和比较对象',
 '存活率以全部初始原子为分母。条件保真度只对幸存者构成的通道平均，目标门需要明确。把S×F作为损失记零下的平均得分时应说明定义，不能把它直接叫作条件量子通道保真度。论文IRB定义也不能由经典存活曲线替代。',['cli','simulation/level2_joint_transfer/src/level5_finite_pulse_irb/channel_tomography.py','paper']);
 table(s,[['量','回答什么','如何理解'],['存活率 S','原子还在不在目标阱里','以全部初始原子为分母'],['相干对比度 C','幸存原子的相位是否整齐','共同相位与相位分散不同'],['条件平均保真度 F','幸存原子的操作离目标多近','依赖目标操作和内部态演化'],['同时考虑损失的指标','留下且状态可用的总体表现','需要明确损失如何计分']],[1.4,1.75,1.9],{h:355,size:25});
 bottom(s,'存活曲线改善，不能直接推出量子比特保真度达到某个数值');
}
// 56
prose('模型目前能够解释到哪里','数值算得细，与物理模型足够完整，是两件事',[
 ['已经纳入','三维热初态、时变双阱、跨轮连续运动、位点差异和噪声。'],
 ['仍需外部测量','转移前温度、原子平面的噪声谱、真实波形和焦点分布。'],
 ['未完整描述','精确运动量子态、真实成像过程、全部光学像差与共同噪声。']],
 '增加模拟原子数，不能替代对实际装置的测量',
 '还要注意初态是近似制备分布，损失判据是检查点吸收边界，白噪声谱形尚未测定，位点之间的相关噪声没有作为本次MC区间的实验误差模型。报告中静态寿命检查已经排除一种过强的SLM白噪声解释，这说明独立实验约束对选模型很重要。',['report','physics','valid']);
// 57
prose('这套 Monte Carlo 的完整物理链条','每一步都对应一个明确的物理问题',[
 ['准备一批原子','先确定实际势阱，再抽取已束缚的热初态。'],
 ['让原子经历真实操作','不断计算力、更新运动，并加入有明确含义的噪声。'],
 ['保留历史并统计结局','连续做完多轮，在指定检查点判断，再汇总存活概率。']],
 '最终得到许多可能轨迹的统计，其可信度取决于物理输入',
 '结束主报告时可以强调：随机抽样只是把实验中的差异放进来，主要物理仍在势阱、受力、时间安排和噪声里。结果告诉我们哪种操作更容易激发和丢失，也告诉我们目前还有哪些实验行为没有被解释。',['guide','report']);
// 58
{
 const s=page('补充：常用公式与单位','符号只保留理解模型所必需的部分',
 '本页供答疑。μK形式的阱深要乘kB才能成为焦耳，公式里的D使用能量单位。公式σr和σz属于单束、轴向对称高斯阱底的谐振近似，s1、s2有偏移或两阱同时作用时不能直接套用作精确热分布。',['init','engine','physics']);
 table(s,[['物理量','公式或关系'],['势能与受力','U = −D exp[−2(x−c)²/w²]，F = −∇U'],['总能量','E = ½m|v|² + U'],['束缚判据','远处势能取零，E < 0 表示有能量束缚'],['初态速度宽度','σv = √(kBT/m)'],['初态位置宽度','σr = √(kBTw²/4D)，σz = √(kBTzR²/2D)'],['光束与频率','zR = πw²/λ，ω = 2πf']],[1.15,3.5],{h:395,size:24});
}
// 59
{
 const s=page('补充：讲解文档中需澄清的表述','依据当前物理实现与归档结果补充',
 '这些澄清不会改变文档作为思路参考的价值，但会影响物理结论。ML最新归档15个节点、14个间隔与旧描述14控制点不同。逐比特复现也需相同环境，固定种子本身不保证跨平台严格相同。最新随机存活校准与旧有限脉冲连续引擎是不同计算路径。',['guide','protocol','noise','physics','report']);
 table(s,[['容易混淆的说法','准确理解'],['每小步都检查存活','每小步算运动，短协议在等待末和放回末检查'],['末态能量减初态能量就是加热','还包含操作做功，应区分势能基准与激发能'],['随机加热与固定功率是一回事','两者对能量分布和高能尾部的描述不同'],['优化曲线有 14 个采样点','最新回放有 15 个节点，构成 14 个时间间隔'],['最新标定已同时算出保真度','该次最终随机模型专门验证存活']],[1.7,2.65],{h:382,size:23});
}
// 60
prose('补充：存活曲线怎样作经验拟合','拟合用于概括曲线形状，轨迹计算本身不依赖这个公式',[
 ['常见形式','Sn = p0 × pⁿ × [1 − exp(−1/(bn))]，用于 n > 0。'],
 ['参数怎样理解','p0 控制整体幅度，p 描述每次操作的指数因子，b 描述弯曲程度。'],
 ['解释的边界','它们共同描述损失曲线，p 不自动等于实际每一步的条件存活率。']],
 '不确定度重抽样时，要把同一原子的整段记录一起抽取',
 '这是仓库Level2C使用的经验拟合形式。n=0需要单独定义或用n趋于0的极限，不能直接把0代入1/(bn)。S0=1来自已束缚初态制备，p0是否自由应看具体拟合约定。本报告主要结果直接比较37个实验散点，不把密集拟合线当作独立测量，也没有重新拟合这些参数。',['simulation/level2_joint_transfer/src/level2c_pickup_survival/survival_fit.py','report']);
// 61
prose('补充：量子态计算的三个层次','公式用于答疑，物理含义与前面的运动过程相连',[
 ['沿轨迹累积光移相位','φ = ∫ ηU[r(t), t] / ℏ dt，光强历史决定相位历史。'],
 ['脉冲期间计算两态演化','H = ℏ[Δσz + Ω(cos θ σx + sin θ σy)] / 2。'],
 ['与指定目标操作比较','Favg = {平均 |Tr(Uideal†U)|² + 2} / 6，只对幸存原子平均。']],
 '相位、有限脉冲与目标操作，都要说明清楚才能解释保真度',
 '第一行是没有理想翻转记号的裸相位积分。理想瞬时DD可通过y(t)乘在积分中理解，有限驱动则使用第二行的哈密顿量。θ为微波驱动相位，Δ为失谐，σ为泡利矩阵。第三行适用于本模型的单比特随机幺正通道，平均在明确的幸存样本集合上完成。',['engine','cli']);
// 62
{
 const s=page('补充：一颗原子的第一小步','一维静止阱的计算例子，用来理解数值更新',
 '以下为单独构造的教学初态，非完整校准轨迹。Cs质量取2.2069469514e−25 kg，阱深140 kB·μK，束腰1.17 μm，初始x=0.2 μm、v=0.03 m/s，步长0.0125 μs。直接按高斯力和Verlet公式计算。三维正式传播是同时对三个方向做同样的更新。',['engine','init']);
 const kb=1.380649e-23,m=132.90545196*1.66053906660e-27,w=1.17e-6,d=140*kb*1e-6,dt=.0125e-6;
 const f=x=>-4*d*x/(w*w)*Math.exp(-2*x*x/(w*w));const x0=.2e-6,v0=.03,a0=f(x0)/m,vhalf=v0+a0*dt/2,x1=x0+vhalf*dt,v1=vhalf+f(x1)/m*dt/2;
 table(s,[['当前或更新的量','数值','这一小步的含义'],['初始位置 / 速度','0.2000 μm / 3.0000 cm/s','在阱中心右侧，并向右运动'],['初始加速度',`${a0.toFixed(0)} m/s²`,'力指向左侧，开始减速'],['先更新半步速度',`${(vhalf*100).toFixed(4)} cm/s`,'位置还没有更新'],['更新位置',`${(x1*1e6).toFixed(6)} μm`,'因为仍向右运动，位置略增'],['再更新半步速度',`${(v1*100).toFixed(4)} cm/s`,'使用新位置处的力']],[1.3,1.7,1.8],{h:375,size:24});
 bottom(s,'一次只走很短的一步，重复许多步才得到完整轨迹');
}
// 63
prose('补充：为什么需要这么多计算','400 μs 短交接为例，忽略少量额外脉冲边界',[
 ['一轮有多长','400 + 100 + 400 = 900 μs。'],
 ['一颗原子需要多少步','每步 0.0125 μs，一轮约 72000 步，30 轮约 216 万步。'],
 ['一批原子需要多少步','4096 颗原子全部跑满，约 88.5 亿个机械时间步。']],
 '丢失的原子会提前结束，但统计时仍保留它的丢失记录',
 '这里是没有提前丢失时的上限估算：4096×30×72000=8,847,360,000。每步通常两次使用力信息，具体实现可复用上一步的末力。最新存活专用模型不加入微波脉冲边界，旧连续量子模型会插入精确脉冲边沿。算得更多主要用于减小数值和抽样误差，无法让错误的物理假设自动变正确。',['protocol','engine','noise']);
// 64
{
 const s=page('资料来源与阅读顺序','详细出处和数据口径保留在各页备注中',
 '本报告以指定分支4d5e200的内容为读取起点。主线文档仅作为思路，主要物理细节核对初始化、协议、连续传播和随机存活实现；定量结果直接读取2026-09-16共同标定的CSV，没有把教学抽样当作新验收结果。原始论文与Savard论文用于核对物理概念，所有新增图表都是原生可编辑数值图表。',['guide','report','params','physics','valid','paper','savard']);
 rows(s,[['入门说明','《蒙特卡洛模拟·大白话完整说明》，2026-09-21。'],['本仓库的定量证据','2026-09-16 共同标定报告、参数表、逐点曲线和数值验证。'],['原始论文','Manetsch 等，A tweezer array with 6100 highly coherent atomic qubits。'],['噪声加热理论','Savard 等，Physical Review A 56, R1095（1997）。']],{gap:102,head:26,body:25});
 bottom(s,'实验事实、模型假设、教学示意和模拟结果均分别标注');
}

if(p.slides.items.length!==64)throw new Error(`Expected 64, got ${p.slides.items.length}`);
await fs.writeFile(path.join(build,'presentation.json'),JSON.stringify(p.toProto()));
await fs.writeFile(path.join(here,'讲解提纲与来源.md'),`# 光镊中的原子：Monte Carlo 物理讲解\n\n共64页。读取基准：docs/monte-carlo-plain-guide，4d5e200。\n\n`+notes.map(r=>`## ${r.n}. ${r.title}\n\n${r.sub}\n\n${r.note}\n\n来源：\n${r.refs.map(v=>`- ${v}`).join('\n')}\n`).join('\n'));
await fs.writeFile(path.join(build,'owners.json'),JSON.stringify({chartOwners:[...new Set(chartOwners)],tableOwners:[...new Set(tableOwners)]}));
const candidate=path.join(build,'candidate.pptx');
await(await PresentationFile.exportPptx(p)).save(candidate);
console.log(`Exported ${p.slides.items.length} slides; charts on ${new Set(chartOwners).size} slides.`);
const renderOnly=process.env.RENDER_ONLY==='1';
const sample=process.env.SAMPLE_SLIDES?process.env.SAMPLE_SLIDES.split(',').map(Number):[];
if(sample.length){
 for(const n of sample){const blob=await p.export({slide:p.slides.items[n-1],format:'png',scale:1});await fs.writeFile(path.join(build,`preview-${String(n).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));}
 console.log('Sample previews complete.');
}
if(!renderOnly){
 const out=process.env.FINAL_PPTX||path.join(here,'output','Monte_Carlo_光镊原子_物理原理与逐步实现_完整版.pptx');
 const result=await finalizePresentation({workspaceDir:root,candidatePath:candidate,finalPath:out,pythonExecutable:python,
  integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...[...new Set(tableOwners)].flatMap(n=>['--require-native-table-slide',String(n)])],
  explicitTotalSlideCount:64,requiredNativeChartOwnerSlides:[...new Set(chartOwners)],requiredNativeTableOwnerSlides:[...new Set(tableOwners)],
  fontPolicy:{basis:'design',families:[FONT]},materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,
  receiptPath:path.join(build,`${path.basename(out)}.validation.json`)});
 console.log(JSON.stringify({finalPath:result.finalPath,slides:64,charts:result.packageIntegrity.chart_count,sha256:result.finalSha256}));
}
