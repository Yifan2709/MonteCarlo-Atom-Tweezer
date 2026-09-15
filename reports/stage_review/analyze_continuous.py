"""Build the v10 report from completed continuous simulations, with no MC rerun."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
DATA = REPO/"simulation/level2_joint_transfer/outputs/continuous_level2c_345/run_20260914"
OUT = Path(__file__).parent/"output"
STEM = "AOD_SLM_continuous_level2c_345_2026-09-14_v10"
KEYS = ["manual_200us","manual_400us","manual_600us","ml_400us"]
LABELS = ["手工 200 μs","手工 400 μs","手工 600 μs","论文 ML 400 μs"]
BEST = dict(zip(KEYS,[5,15,25,25]))
rows = []
for file in sorted(DATA.glob("*.json")):
    row = json.loads(file.read_text(encoding="utf-8"))
    if "case_id" in row:
        row["result_file"] = str(file.relative_to(REPO))
        rows.append(row)


def pick(stage,t,arm,key,shots,dt):
    found = [r for r in rows if (r["stage"],r["temperature_uK"],r["arm"],r["waveform"],r["shots"],r["dt_us"]) == (stage,t,arm,key,shots,dt)]
    if len(found) != 1:
        raise ValueError(f"Missing/duplicate case {stage,t,arm,key,shots,dt}")
    return found[0]


def endpoint(row):
    cp = row["comparison_raw"]
    n = cp["n"][-1]
    return n, row["survival"][n], cp["reference"][-1]


def table(headers,data):
    safe = lambda x: str(x).replace("|", "\\|")
    return "\n".join(["| " + " | ".join(map(safe,headers)) + " |",
                       "| " + " | ".join(["---"]*len(headers)) + " |"] +
                      ["| " + " | ".join(map(safe,row)) + " |" for row in data])


scan = []
for t in [5,10,15,20,25,30,35,40]:
    group = [pick("scan",t,"matched",k,256,.05) for k in KEYS]
    scan.append({"temperature_uK":t,"pooled_rmse":float(np.sqrt(np.mean([r["comparison_raw"]["rmse"]**2 for r in group]))),
                 "max_deviation":max(r["comparison_raw"]["max_abs_deviation"] for r in group),
                 "all_curves_pass":all(r["comparison_raw"]["matches_reference"] for r in group)})
assert min(scan,key=lambda r:r["pooled_rmse"])["temperature_uK"] == 15
base = {str(t):{a:[pick("validation",t,a,k,512,.05) for k in KEYS]
                for a in ("baseline_1d","matched_lensing_off","matched","extended")} for t in (5,40)}
precision = {str(t):[pick("precision",t,"matched",k,4096,.0125) for k in KEYS] for t in (5,15,25)}
selected = [pick("precision",BEST[k],"matched",k,4096,.0125) for k in KEYS]
convergence = []
for t in (5,15,25):
    for key in KEYS:
        a,b = (pick("precision",t,"matched",key,4096,dt) for dt in (.025,.0125))
        na = np.load(REPO/Path(a["result_file"]).with_suffix(".npz"))
        nb = np.load(REPO/Path(b["result_file"]).with_suffix(".npz"))
        np.testing.assert_array_equal(na["initial_pos"],nb["initial_pos"])
        np.testing.assert_array_equal(na["initial_vel"],nb["initial_vel"])
        delta = np.asarray(b["survival"])-a["survival"]
        convergence.append({"temperature_uK":t,"waveform":key,"max_survival_change":float(max(abs(delta))),
                            "endpoint_change":float(delta[endpoint(b)[0]]),
                            "final_label_disagreement":float(np.mean(na["alive"][:,-1]!=nb["alive"][:,-1])),
                            "aggregate_within_0p02":bool(max(abs(delta))<=.02)})
coarse = []
for t in (5,40):
    for arm in ("baseline_1d","matched_lensing_off","matched","extended"):
        for key in KEYS:
            a,b = pick("validation",t,arm,key,128,.05),pick("convergence",t,arm,key,128,.025)
            coarse.append({"temperature_uK":t,"arm":arm,"waveform":key,
                           "max_survival_change":float(max(abs(np.asarray(a["survival"])-b["survival"])))})

confidence = []
for row in selected:
    cp = row["comparison_raw"]
    nr,sr = np.array(cp["n"]),np.array(cp["reference"])
    lo,hi = np.array(row["wilson95_low"])[nr],np.array(row["wilson95_high"])[nr]
    beyond = (lo>sr+.02)|(hi<sr-.02)
    worst = int(np.argmax(abs(np.array(cp["deviation"]))))
    confidence.append({"waveform":row["waveform"],"temperature_uK":row["temperature_uK"],
                       "worst_n":cp["n"][worst],"max_deviation":cp["max_abs_deviation"],
                       "points_with_pointwise95_interval_outside_band":int(beyond.sum()),
                       "matches_raw":cp["matches_reference"],
                       "matches_v8_anchors":row["comparison_v8_anchors"]["matches_reference"]})

manifest = []
for row in rows:
    for file in (REPO/row["result_file"],(REPO/row["result_file"]).with_suffix(".npz")):
        manifest.append({"path":str(file.relative_to(REPO)),"sha256":hashlib.sha256(file.read_bytes()).hexdigest()})
evidence = {"stem":STEM,"base":base,"temperature_scan":scan,"precision":precision,
            "selected":selected,"confidence":confidence,"convergence":convergence,
            "original_step_check":coarse,"result_case_count":len(rows),"result_manifest":manifest,
            "verdict":"No tested shared setting or per-waveform scan selection passes every Fig.6d point within ±0.02",
            "protocol_warning":"Extended 375 μm out/back protocol is not the Fig.6d experiment"}
evidence["source_snapshot"] = {"path":str((DATA/"reproducibility_sources.zip").relative_to(REPO)),
                               "sha256":hashlib.sha256((DATA/"reproducibility_sources.zip").read_bytes()).hexdigest()}
assert all(not r["matches_raw"] and not r["matches_v8_anchors"] for r in confidence)
assert all(r["points_with_pointwise95_interval_outside_band"] > 0 for r in confidence)
OUT.mkdir(exist_ok=True)
(OUT/f"{STEM}.evidence.json").write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding="utf-8")

# Standalone scientific plots for the written report. PPT uses native charts.
fig,axes = plt.subplots(2,2,figsize=(13,8),sharey=True)
for i,(key,ax) in enumerate(zip(KEYS,axes.flat)):
    row = selected[i]
    n,ref = np.array(row["comparison_raw"]["n"]),np.array(row["comparison_raw"]["reference"])
    ax.fill_between(n,np.maximum(ref-.02,0),np.minimum(ref+.02,1),color="black",alpha=.12,label="Paper ±0.02")
    ax.plot(n,ref,"k.",ms=3)
    for data,label,color,style in [(base["5"]["baseline_1d"][i],"1D T=5", "#667085","--"),
                                    (base["5"]["matched"][i],"3D T=5","#1570EF","-"),
                                    (row,f"Selected T={BEST[key]}, N=4096","#D92D20","-")]:
        ax.plot(n,np.asarray(data["survival"])[n],style,color=color,lw=1.8,label=label)
    ax.set_title(key.replace("_"," "))
    ax.set_xlabel("Short one-way transfer count n")
    ax.set_ylabel("Survival fraction")
    ax.set_ylim(0,1.03)
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT/f"{STEM}_curves.png",dpi=180)
plt.close(fig)

lines = ["# Level 2C 接入 Level 3/4/5 连续全流程重算", "", "2026-09-14 · v10", "",
"**结论：接入三维运动、AOD 柱面透镜、连续往返及有限脉冲后，曲线发生明显变化，但本次扫描与独立复核中，没有找到四条波形同时满足原 Level 2C 全曲线 ±0.02 标准的共同参数。即使分别使用每条波形在扫描池中最接近的参数，4096 个独立样本的整条曲线仍未通过。**", "",
"这是否定本次已测试场景的复现结果，不是证明任何参数或任何补充模型都不可能拟合。波形被冻结，未重新优化为其他形状。", "",
f"![四条存活曲线]({STEM}_curves.png)","",
"## 连续接入了什么", "",
"- Level 2C：原手工 200/400/600 μs 波形、论文 Fig. 6c 的原始 15 节点普通三次样条回放、SLM 开关、100 μs 等待、原加热公式和无序幅度。",
"- Level 3：三维高斯束缚热初态、三维速度 Verlet、crossed-AOD 柱面透镜，随速度移动的两个轴向焦点。",
"- Level 4 对应的连续流程：每个原子的 r、v 从上一段/上一轮末态接着推进。拾取、等待、放回或附加长运输均使用同一个原子。失阱为吸收事件。",
"- Level 5：沿同一条实时轨迹取差分光移，连续积累 SU(2) 传播子。有限 XY4 脉冲期间继续机械运动，脉冲起止精确切分积分步。输出无 DD、有限 DD 的条件对比度、相对理想脉冲序列的平均保真度及 PTM。",
"- 新连续入口是 `continuous_transfer`。原 Level 5 的 Clifford RB、47/195 站点和独立轨迹池仍是另外的实验协议。本次没有把其历史 IRB 拟合冒充连续模型结果。", "",
"## 协议、距离与计数", "",
table(["项目","matched：Fig. 6d 对照","extended：额外长运输"],[
["一轮","拾取 + 100 μs SLM 关闭等待 + 放回","拾取 + 等待 + 375 μm 去程 + 54 μs 保持 + 375 μm 回程 + 放回"],
["SLM","拾取/放回时开，等待时关","长运输及远端保持时也关"],
["400 μs 波形单轮时长","900 μs","3854 μs"],
["手工单轮名义路程","4.8 μm","754.8 μm"],
["n=60","30 轮，约 144 μm","30 轮，约 22.644 mm"],
["验收用途","可沿用 Level 2C 参考曲线","协议发生变化，只作扩展预测"]]),"",
"ML 数字化样条保留原始细小反向起伏：单程端点差约 2.400115 μm，积分实际路程约 2.430070 μm。文件中的 `distance_per_round_um` 是端点差对应的往返路程；以上 4.8/754.8 μm 是手工波形与名义协议口径。短拾取沿 x 轴，长运输沿 xy 对角线，总长 375 μm，每轴约 265.17 μm。", "",
"## 原 Level 2C 参数下的重算", "",
"以下使用独立复核池 512 shots、seed=23002、dt=0.05 μs。手工 200 μs 的论文最后有效点是 n=31，其余使用 n=60。", ""]
for t in (5,40):
    lines += [f"### T = Tref = {t} μK", "", table(["波形/末点","论文","1D 新基线","3D 透镜关闭","3D 透镜开启","另加长运输"],
        [[LABELS[i]+f" / n={endpoint(base[str(t)]['matched'][i])[0]}",f"{endpoint(base[str(t)]['matched'][i])[2]:.4f}"]+
         [f"{endpoint(base[str(t)][a][i])[1]:.4f}" for a in ("baseline_1d","matched_lensing_off","matched","extended")]
         for i in range(4)]),""]
lines += ["1D 新基线也重新抽取逐段抖动：随机流按原始 shot ID 分配，避免不同波形发生损失后改变后续随机样本对应关系。其数值不要求与 v8 的旧随机流逐点相同。保持了相同分布与种子约定。三维加热保持总功率，未将其乘以 3。", "",
"对准偏移和逐段抖动沿原 Level 2C 的 x 轴施加，未追加未经标定的 y/z 随机扰动。深度、束腰的静态无序则作用于整个三维势。", "",
"## 共同参数扫描与冻结后复核", "",
"扫描池固定 256 shots、seed=23001。候选 T 为 5、10、15、20、25、30、35、40 μK。按照原 Level 2C 默认行为，Tref 随 T 同时改变；这是一维的耦合参数场景扫描，不能把拟合所得数值称为独立测得的初温。四条曲线等权，目标是各曲线均方误差平均值的平方根。", "",
table(["T=Tref / μK","共同 RMSE","四曲线全部通过"],[[r["temperature_uK"],f"{r['pooled_rmse']:.5f}","否"] for r in scan]),"",
"共同最优扫描格点是 15 μK。下面给出 4096 shots、全新 seed=23003、dt=0.0125 μs 的冻结后复核；该样本池未参与参数选择。", "",
table(["波形","共同 T=15 的末点","论文末点","最大整曲线偏差","通过 ±0.02"],
      [[LABELS[i],f"{endpoint(r)[1]:.4f}",f"{endpoint(r)[2]:.4f}",f"{r['comparison_raw']['max_abs_deviation']:.4f}","否"] for i,r in enumerate(precision["15"])]),"",
"如果允许每条波形使用不同的扫描最优参数，结果如下。这仅用于诊断参数矛盾，不能当作同一个实验条件的复现。", "",
table(["波形","扫描选择 T","复核末点","论文末点","最大偏差 / 最差 n","点在 ±0.02 内的比例"],
      [[LABELS[i],r["temperature_uK"],f"{endpoint(r)[1]:.4f}",f"{endpoint(r)[2]:.4f}",
        f"{r['comparison_raw']['max_abs_deviation']:.4f} / {confidence[i]['worst_n']}",
        f"{r['comparison_raw']['fraction_within_band']:.1%}"] for i,r in enumerate(selected)]),"",
"原始数字化数据与 v8 固定锚点的单调处理结果均单独保留，未为了通过验收改变参考数据。每条精细复核曲线至少有一个点，其逐点 Wilson 95% 区间完全落在论文 ±0.02 带外。这里没有把逐点区间当作整曲线的同时置信带。", "",
"## 步长与 Monte Carlo 不确定性", "",
f"原 128-shot、0.05→0.025 μs 配对检查的最大曲线变化为 {max(r['max_survival_change'] for r in coarse):.4f}，超过 0.02。故仅靠原配置的小样本检查不足以宣称严格的 2% 数值精度。", "",
"新增精细检查对每个格点使用完全相同的 4096 个初态、静态无序和逐段抖动，比较 0.025 与 0.0125 μs。下表分别列出系综曲线差和最后一个 checkpoint 的逐原子存活标签差；混沌轨迹的标签不一致并不等于系综存活率同幅度变化。", "",
table(["T","波形","最大 |ΔS(n)|","末态标签不一致比例","曲线差≤0.02"],
      [[r["temperature_uK"],r["waveform"],f"{r['max_survival_change']:.4f}",f"{r['final_label_disagreement']:.2%}","是" if r["aggregate_within_0p02"] else "否"] for r in convergence]),"",
"4096 shots 在 S≈0.5 处的单点 95% 区间半宽约 0.0153。模型未通过论文验收的结论，依据整条曲线偏差及更细步长的重复结果，而不是只比较一位小数的末点。", "",
"## 有限脉冲与内部态", "",
"ηSLM=ηAOD=1.3×10⁻⁴，Ω/2π=24.611 kHz，bare π 脉冲约 20.316 μs。短协议每轮全程对称 XY4；长协议在每个长运输段各放一个 XY4。均未额外延长协议时间。Favg 相对于理想驱动序列计算，且只对当时存活原子取平均。", "",
table(["400 μs / 512 shots","S(60)","C 无 DD","C 有限 DD","Favg 有限 DD"],
      [[f"{a}, T={t}",f"{r['survival'][-1]:.4f}",f"{q['contrast_none']:.4f}",f"{q['contrast_dd']:.4f}",f"{q['favg_dd']:.4f}"]
       for t in (5,40) for a in ("matched","extended")
       for r in [base[str(t)][a][1]] for q in [r["quantum"]["checkpoints"][-1]]]),"",
"该模型没有自旋依赖机械力或自发 Raman 丢失，因此无 DD 与 DD 的存活率严格相同。DD 改变相位和相干性，不能通过纯幺正内部态演化补出实验的原子丢失。对比度高也不自动等于相对指定目标门的保真度高。", "",
"## 差距的可能来源与证据边界", "",
"1. **维度升级改变了可用热能和束缚判据。** 相同的每自由度温度在三维中具有更多运动能量；40 μK 场景发生明显过度损失。不同波形需要不同拟合参数，说明只调共同 T/Tref 不能解决形状差异。",
"2. **透镜效应对快速波形更敏感。** 200 μs 波形的 5 μK 末点从 3D 无透镜的约 0.467 降至约 0.270，而 600 μs 与论文 ML 基本不丢失，不能产生论文所需的统一衰减。v_s=0.53917 m/s 是旧 Level 3 severity=1 的相对场景，尚非器件独立标定。",
"3. **长运输增加演化时间和机械激发。** 它可显著降低部分存活率，但 600 μs 与 ML 的低温末点仍接近 1；而且 extended 与 Fig. 6d 的实验协议不同。",
"4. **原加热模型仍是常数平均功率近似。** 初温与 Tref 耦合、RIN/指向谱及系数约定仍待独立核对。本次刻意保留原公式以隔离维度影响，没有用额外损失项把曲线强行拉到论文上。",
"5. **波形与装置标定仍有缺项。** 手工位置曲线的具体解析形状、ML 数字化误差、实际 AWG/RF/光强响应和绝对 AOD 透镜参数可能影响结果。本次没有引入未标定的真空、成像、态制备或 Raman 损失。", "",
"## 复现与文件", "",
"```bash\npython -m pip install -e './simulation/level2_joint_transfer[dev,continuous]'\npython run_simulation.py run full --stage preflight\npython run_simulation.py run full --stage all --output-dir simulation/level2_joint_transfer/outputs/continuous_level2c_345/new_run\n```", "",
"中间温度扫描、同样本量粗细步长配对及冻结后的 4096-shot 精细确认命令见 `reports/stage_review/README.md`。`continuous_precision.py` 明确锁定扫描选择、seed 和两档步长。", "",
f"本次读取 {len(rows)} 个实际完成的计算格点，每点保存 JSON 指标和 NPZ 逐原子记录。NPZ 包含所有 61 个 checkpoint 的存活掩码、六维状态、相位及 SU(2) 系数，以及损失阶段、停止时间和实际累计加热。数值源文件与结果 SHA-256 位于 `{STEM}.evidence.json`。", "",
"新增回归覆盖三维力与 Level 3 一致、SU(2) 与 Level 5 一致、一维极限与旧 Level 2C 一致、第二轮从第一轮末态继续、损失吸收与加热停止、DD 无机械反馈、原始 ML 回放及有限脉冲时间边界。"]
(OUT/f"{STEM}.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
print(json.dumps({"cases":len(rows),"convergence_max":max(r["max_survival_change"] for r in convergence),
                  "precision_acceptance":[r["matches_raw"] for r in confidence]},ensure_ascii=False))
