"""Audit archived early losses; run bounded, paired diagnostic ablations.

Production inputs/results are not modified. All ablations use the archived
4096-shot initial states and disorder, except the explicitly scaled cold state.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
from dataclasses import replace

import numpy as np

from continuous_transfer.cli import case_inputs
from continuous_transfer.engine import run_continuous
from continuous_transfer.protocol import waveform_for, build_protocol
from level3_3d_transport_lensing.gaussian_3d import static_potential

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "simulation/level2_joint_transfer"
CASE = "precision_T25_matched_manual_600us_N4096_dt0.0125"
SOURCE = PROJECT / "outputs/continuous_level2c_345/run_20260914"
OUT = Path(__file__).resolve().parent / "output/manual600_early_audit"
K = 1.380649e-29


def summarize_result(result):
    survival = result["alive"].mean(axis=0)
    losses = np.where(result["lost_stage"] < 0, 999,
                      2 * result["lost_round"] + np.where(result["lost_stage"] == 1, 1, 2))
    return {
        "survival": survival.tolist(),
        "one_way_conditional_loss": (1 - survival[1:] / survival[:-1]).tolist(),
        "round_conditional_loss": (1 - survival[2::2] / survival[:-2:2]).tolist(),
        "losses_first_10": int((losses <= 10).sum()),
        "losses_by_checkpoint": {str(n): int((losses == n).sum()) for n in range(1, len(survival))},
    }


def render_report(report):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.unicode_minus": False,
                         "font.size": 11, "svg.fonttype": "none"})
    s = np.array(report["archived"]["survival"])
    paper = report["paper_reference"]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.05, 1]})
    ax = axes[0]
    ax.plot(paper["display_n"], paper["display_fit"], color="#245f9e", lw=2.2,
            label="论文原拟合线（PDF 矢量）")
    points = [p for p in paper["markers"] if p["n"] <= 16]
    ax.errorbar([p["n"] for p in points], [p["survival"] for p in points],
                yerr=[[p["survival"] - p["lower"] for p in points],
                      [p["upper"] - p["survival"] for p in points]],
                fmt="D", ms=4, capsize=3, color="#245f9e", label="论文散点及原误差条")
    ax.plot(np.arange(len(s)), s, color="#be4d36", lw=2.1, marker="o", ms=3,
            label="当前 MC：25 μK")
    ax.plot(np.arange(len(s)), report["constant_fraction_diagnostic"]["survival"],
            color="#848b94", lw=1.6, ls="--", label="同首末点的纯指数参照")
    ax.set(xlim=(0, 16), ylim=(.85, 1.008), xlabel="单程转移次数 n（2 次 = 1 往返）",
           ylabel="存活率", title="600 μs：初期损失超出纯指数参照")
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.grid(alpha=.18)
    names = ["baseline", "no_mean_heating", "no_lensing", "ideal_sites_and_alignment",
             "cold_initial_state_same_heating"]
    labels = ["原模型", "关闭平均\n加热", "关闭柱面\n透镜", "均匀阱深\n且理想对准", "压低初态能量\n加热不变"]
    losses = [1 - report["ablations"][name]["survival"][10] for name in names]
    ax = axes[1]
    ax.bar(np.arange(5), losses, color=["#be4d36", "#8898ae", "#8898ae", "#8898ae", "#31886a"], width=.63)
    for i, loss in enumerate(losses):
        ax.text(i, loss + .003, f"{loss:.2%}", ha="center", va="bottom", fontsize=11)
    ax.set_xticks(np.arange(5), labels, fontsize=10)
    ax.set(ylim=(0, .12), ylabel="前 5 个往返的累计丢失比例",
           title="前 5 个往返的机制对照")
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.grid(axis="y", alpha=.18)
    ax.set_axisbelow(True)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("初期快速下降主要关联初始高能尾部，关闭平均加热后仍然存在", fontsize=15, y=.98)
    fig.text(.05, .025,
             "诊断：4096 原子，最大步长 0.0125 μs，5 个往返；冷初态仅缩放初始位置与速度，不代表实测温度。\n"
             "虚线以 MC 的 S(0)、S(60) 定义，不是论文拟合；论文参考未平滑、未重拟合。", fontsize=9, color="#505b69")
    fig.subplots_adjust(left=.085, right=.985, top=.85, bottom=.20, wspace=.24)
    fig.savefig(OUT / "manual600_early_diagnosis.png", dpi=180)
    fig.savefig(OUT / "manual600_early_diagnosis.svg")
    plt.close(fig)
    ablation_rows = "\n".join(
        f"| {label.replace(chr(10), '')} | {report['ablations'][name]['survival'][10]:.4%} | "
        f"{report['ablations'][name]['losses_first_10']} / 4096 |"
        for name, label in zip(names, labels))
    group_rows = "\n".join(
        f"| 第 {b['after_round'] + 1}–{b['through_round']} 往返 | {b['geometric_conditional_loss']:.3%} |"
        for b in report["round_hazard_blocks"])
    text = f"""# 手工 600 μs 初期存活率诊断

2026-09-15。当前展示条件 T = Tref = 25 μK；Fig. 6 短转移协议，无长距离运输。

## 结论

初始高能尾部在转移过程中优先丢失，是当前模型初期快速下降的主要证据指向；初始采样与实际格点阱深不一致还引入了本不应计入初始存活池的原子。不能把初期快速下降主要归因于平均噪声加热或柱面透镜。

这不是对实验真实温度或真实损失机制的测量。本轮未修改生产模型、既有参考数据或 PPT。

## 指数下降的判断

纯指数 S(n)=S(0)p^n 在普通坐标中本来就前陡后缓；其特征是条件损失比例不随次数变化，而不是每次减少同样数量。两种单程操作的检查点不同，因此这里按完整往返比较。

| 往返区间 | 几何等效每往返条件损失 |
|---|---:|
{group_rows}

MC 的 S(2)={s[2]:.4%}、S(10)={s[10]:.4%}、S(60)={s[60]:.4%}。以相同首末点构造纯指数参照 S_exp(n)=S_MC(60)^(n/60)，其 S_exp(10)={report['constant_fraction_diagnostic']['survival'][10]:.4%}，比 MC 高约 2.17 个百分点。第 10–60 次单程的 log(S) 线性描述约为 S(n)=0.97451×0.992745^n，仅作形状诊断，不作为独立实验拟合或修正曲线。

论文 Methods 采用 S_n=p0·p^n·[1-exp(-1/(b·n))]（n>0），本身也不是纯指数；作者将指数分量可能归于实验不完善。当前读取的论文原拟合线在 n=10 约为 97.2901%，并非该位置的独立测量点。

## 配对短程对照

使用同一 4096 原子的归档初态、同一格点参数与逐程抖动序列，重跑前 5 个往返。基线逐原子存活状态与归档前缀完全一致。最大步长 0.0125 μs。

| 改动 | 第 10 次单程存活率 | 已丢失原子 |
|---|---:|---:|
{ablation_rows}

关闭平均加热：保留位点差异、逐程抖动和运动势阱。关闭透镜：其余输入不变。均匀格点诊断同时去掉阱深/束腰散布、静态对准和抖动，不能进一步拆分各项贡献。冷初态诊断仅将初始位置、速度各乘 sqrt(5/25)，保留原有 Tref=25 μK 的加热；它不是精确有限阱深平衡分布，也不是实验温度标定。

关闭平均加热和透镜的 S(10) 变化分别为 +0.317、+0.293 个百分点，配对标准误分别约 0.356、0.377 个百分点；这种小差异不足以确认对应改善，但足以看到关闭它们后快速下降仍然存在。机制之间存在相互作用，不能将存活率差直接解释为可相加的损失份额。

## 初态证据与初始化问题

前 10 次单程丢失的 394 个原子，其平均初始激发能 E/kB≈113.65 μK；同期幸存者约 56.58 μK。这里是相对各自 SLM 阱底的能量单位，不是温度。初始激发能小于 40 μK 的 1097 个原子在前 10 次没有损失；大于等于 100 μK 的 520 个原子中约 59.6% 丢失。

代码先按名义 140 μK SLM 阱采样并拒绝 E≥0 的提议，然后独立抽取实际格点深度；引擎却无条件把初始 alive 设为 True。重新按实际格点计算，50/4096=1.2207% 的原子初始 E≥0，其中 42 个在第一个往返被判丢失，50 个均在前 10 次单程内被判丢失。

仅在归档中排除这 50 个初始不束缚原子，并以剩余 4046 个作为分母：第一往返损失由 2.5391% 降为约 1.5324%；第 10 次单程存活率由 90.3809% 升为约 91.4978%，仍明显低于论文原拟合约 97.29%。这只是条件化重统计，不是重新生成正确初态后的模拟。

合理的初始化修正是先确定实际格点势阱，再按对应势阱准备与筛选初态，并明确 S(0) 的实验归一化。修正初始化并不能单独补齐与论文的约 6.91 个百分点初期偏差；高能分布、转移激发及后续加热仍须分别约束。降低温度可以移除初期损失，但现有扫描已显示其可能让后期损失不足，不能据此直接选温度来强制匹配。

## 追溯

- [逐原子统计与对照结果](diagnostics.json)
- [诊断脚本](../../diagnose_manual600_early.py)
- 原始归档：`{report['source']}`
- NPZ SHA-256：`{report['source_sha256']}`
- 初始化：`src/continuous_transfer/cli.py` 的 `case_inputs`；`src/level3_3d_transport_lensing/thermal_sampling_3d.py`。
- 存活初始化与检查点：`src/continuous_transfer/engine.py` 的 `integrate`。
- 论文：根目录 `2403.12021v4.md` Methods，Fig. 6 段落。

![初期下降诊断](manual600_early_diagnosis.png)
"""
    (OUT / "手工600us初期损失核查.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ablations", action="store_true")
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True, parents=True)
    if args.render_only:
        render_report(json.loads((OUT / "diagnostics.json").read_text(encoding="utf-8")))
        return
    metadata_path, archive_path = SOURCE / (CASE + ".json"), SOURCE / (CASE + ".npz")
    m = json.loads(metadata_path.read_text(encoding="utf-8"))
    z = np.load(archive_path)
    cfg, physics, pos, vel, draws, powers, sampling = case_inputs(
        m["settings"] | {"initial_sampling": "legacy_nominal"}, 25, 4096, m["seed"], "matched")
    assert np.array_equal(pos, z["initial_pos"])
    assert np.array_equal(vel, z["initial_vel"])
    w = physics["slm_waist_m"]
    zr = np.pi * w * w / 1061e-9
    depths = physics["slm_depth_j"] * draws.slm_depth_factor
    potential = static_potential(pos[:, 0], pos[:, 1], pos[:, 2], depths, w, zr)
    kinetic = .5 * physics["mass_kg"] * np.sum(vel ** 2, axis=1)
    energy = (kinetic + potential) / K
    excitation = energy + depths / K
    losses = np.where(z["lost_stage"] < 0, 999,
                      2 * z["lost_round"] + np.where(z["lost_stage"] == 1, 1, 2))
    groups = {}
    for name, mask in [("all", losses >= 0), ("lost_by_2", losses <= 2),
                       ("lost_by_10", losses <= 10), ("survive_10", losses > 10),
                       ("survive_60", losses > 60)]:
        groups[name] = {
            "shots": int(mask.sum()),
            "initial_excitation_mean_uK": float(excitation[mask].mean()),
            "initial_excitation_quantiles_10_50_90_uK": np.quantile(excitation[mask], [.1, .5, .9]).tolist(),
            "initial_binding_energy_mean_uK": float(energy[mask].mean()),
            "slm_depth_mean_uK": float(depths[mask].mean() / K),
            "initial_z_rms_um": float(np.sqrt(np.mean(pos[mask, 2] ** 2)) * 1e6),
        }
    bins = []
    for lo, hi in [(0, 40), (40, 60), (60, 80), (80, 100), (100, 1000)]:
        mask = (excitation >= lo) & (excitation < hi)
        bins.append({"excitation_interval_uK": [lo, hi], "shots": int(mask.sum()),
                     "loss_probability_by_10": float(np.mean(losses[mask] <= 10))})
    vectors = json.loads((PROJECT / "reference/fig6d_vector_v2/fig6d_vectors.json").read_text(encoding="utf-8"))
    paper = vectors["series"]["manual_600us"]
    report = {
        "source": str(archive_path.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest(),
        "archived": summarize_result(z), "groups": groups, "initial_energy_bins": bins,
        "initial_unbound_actual_site": int((energy >= 0).sum()),
        "initial_unbound_lost_by_10": int(((energy >= 0) & (losses <= 10)).sum()),
        "heating_power_uK_per_s": {k: v / K for k, v in powers.items()},
        "heat_per_full_round_uK": powers["total"] / K * .0013,
        "paper_reference": paper,
        "limitations": ["Initial excitation is relative to each atom's initial SLM well bottom.",
                        "Energy/loss correlations do not by themselves identify a causal heating channel.",
                        "Ablations cover only five full round trips, not the full 60-transfer experiment.",
                        "Cold diagnostic rescales archived positions/velocities; it is not an exact equilibrium ensemble."],
    }
    s = np.asarray(report["archived"]["survival"])
    n = np.arange(len(s))
    late_fit = np.polyfit(n[10:], np.log(s[10:]), 1)
    report["initial_bound_only_survival"] = z["alive"][energy < 0].mean(axis=0).tolist()
    report["initial_unbound_loss_counts"] = {
        str(i): int(((energy >= 0) & (losses == i)).sum()) for i in range(1, 11)}
    report["constant_fraction_diagnostic"] = {
        "definition": "S_exp(n)=S_MC(60)**(n/60), anchored to S0=1 and S60; not a paper fit",
        "survival": (s[-1] ** (n / 60)).tolist(),
        "late_log_fit_intercept": float(np.exp(late_fit[1])),
        "late_log_fit_one_way_factor": float(np.exp(late_fit[0])),
    }
    report["round_hazard_blocks"] = [
        {"after_round": a, "through_round": b,
         "geometric_conditional_loss": float(1 - (s[2 * b] / s[2 * a]) ** (1 / (b - a)))}
        for a, b in [(0, 1), (1, 5), (5, 10), (10, 20), (20, 30)]]
    print(json.dumps({k: v for k, v in report.items() if k != "paper_reference"}, ensure_ascii=False), flush=True)
    report_path = OUT / "diagnostics.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if not args.ablations:
        return
    import numba
    numba.set_num_threads(12)
    waveform = waveform_for("manual_600us", cfg, PROJECT / m["settings"]["ml_waveform"])
    protocol = build_protocol(waveform, cfg.roundtrip.wait_s, m["dt_us"] * 1e-6,
                              arm="matched", omega=2 * np.pi * m["rabi_frequency_hz"])
    base = dict(rounds=5, physics=physics, draws=draws,
                jitter_sigma_m=cfg.disorder.shot_jitter_alignment_sigma_um * 1e-6,
                power_w=powers["total"], vs_m_s=m["vs_m_s"], eta=m["eta"])
    ideal = replace(draws, slm_depth_factor=np.ones(4096), aod_depth_factor=np.ones(4096),
                    aod_waist_factor=np.ones(4096), alignment_offset_m=np.zeros(4096))
    cases = [
        ("baseline", pos, vel, {}),
        ("no_mean_heating", pos, vel, {"power_w": 0.0}),
        ("no_lensing", pos, vel, {"vs_m_s": 0.0}),
        ("ideal_sites_and_alignment", pos, vel, {"draws": ideal, "jitter_sigma_m": 0.0}),
        ("cold_initial_state_same_heating", pos * np.sqrt(5 / 25), vel * np.sqrt(5 / 25), {}),
    ]
    report["ablations"] = {}
    report["ablation_settings"] = {"shots": 4096, "seed": m["seed"], "dt_us": m["dt_us"],
                                     "rounds": 5, "heating_Tref_uK": 25,
                                     "same_original_atom_ids_and_jitter_streams": True}
    for name, p, v, changes in cases:
        start = time.perf_counter()
        result = run_continuous(p, v, protocol, **(base | changes))
        if name == "baseline":
            assert np.array_equal(result["alive"], z["alive"][:, :11]), "Archived prefix mismatch"
        summary = summarize_result(result)
        summary["runtime_s"] = time.perf_counter() - start
        changed = result["alive"][:, 10].astype(int) - z["alive"][:, 10].astype(int)
        summary["paired_S10_change"] = float(changed.mean())
        summary["paired_S10_change_standard_error"] = float(changed.std(ddof=1) / np.sqrt(len(changed)))
        summary["rescued_vs_baseline_at_10"] = int((changed == 1).sum())
        summary["newly_lost_vs_baseline_at_10"] = int((changed == -1).sum())
        report["ablations"][name] = summary
        np.savez_compressed(OUT / (name + ".npz"), **result)
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(name, json.dumps(summary, ensure_ascii=False), flush=True)
    render_report(report)


if __name__ == "__main__":
    main()
