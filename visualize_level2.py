# -*- coding: utf-8 -*-
"""Level 2 全流程可视化：从仓库输出数据生成完整过程图。

读取 outputs/level2_joint_transfer/demo/ 下的全部数据文件，按 Level 2 的
实际执行顺序（预检查 -> 候选探索 -> 精修选择冻结 -> 独立验证 -> 鲁棒性）
生成 8 张详细图，保存到 level2_figures/。

用法：
    python3 visualize_level2.py [数据目录] [输出目录]
默认数据目录为脚本所在仓库的
    simulation/level2_joint_transfer/outputs/level2_joint_transfer/demo
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei",
                        "Noto Sans CJK SC", "Arial Unicode MS", "sans-serif"],
    "axes.unicode_minus": False,
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "axes.grid": True,
    "grid.alpha": 0.3,
})

KB = 1.380649e-23                       # 玻尔兹曼常数 J/K
D_SLM_UK = 140.0                        # SLM 阱深 μK（与配置一致）
D_AOD_UK = 280.0
W_SLM_UM = 1.17                         # 两阱腰宽假设相同
W_AOD_UM = 1.17
D0_UM = 2.4                             # AOD 初始中心 μm
MASS_U = 132.90545196                   # Cs-133
WAVEFORMS = ["sequential_600us", "sequential_400us", "overlapped_optimized_400us"]
WF_LABEL = {
    "sequential_600us": "顺序基线 600 μs",
    "sequential_400us": "压缩顺序基线 400 μs",
    "overlapped_optimized_400us": "优化同步波形 400 μs",
}
WF_COLOR = {
    "sequential_600us": "#1f77b4",
    "sequential_400us": "#ff7f0e",
    "overlapped_optimized_400us": "#d62728",
}


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def gaussian_well(x_um, depth_uk, waist_um, center_um):
    """一维高斯光阱势（单位：μK），与仓库物理模型一致。"""
    return -depth_uk * np.exp(-2.0 * (x_um - center_um) ** 2 / waist_um ** 2)


# ---------------------------------------------------------------- 通用
def load_all(data_dir: Path) -> dict:
    data = {}
    data["metrics"] = json.loads((data_dir / "metrics.json").read_text())
    data["preflight"] = json.loads((data_dir / "preflight.json").read_text())
    data["frozen"] = yaml.safe_load((data_dir / "best_waveform.yaml").read_text())
    data["candidates"] = read_csv_rows(data_dir / "candidates.csv")
    data["history"] = read_csv_rows(data_dir / "optimization_history.csv")
    data["validation"] = read_csv_rows(data_dir / "validation_shots.csv")
    data["robustness"] = read_csv_rows(data_dir / "robustness.csv")
    with (data_dir / "waveforms.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = np.array([[to_float(v) for v in row] for row in reader])
    data["waveforms"] = {name: rows[:, header.index(name)] for name in header}
    return data


# ---------------------------------------------------------------- 图 1：流程总览
def fig1_pipeline(data: dict, out: Path):
    fig, ax = plt.subplots(figsize=(15, 8.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Level 2 全流程总览：同步移动-降深波形优化与鲁棒蒙特卡洛验证",
                 fontsize=15, fontweight="bold", pad=16)

    metrics = data["metrics"]
    counts = metrics.get("optimization") or {}
    val = metrics["validation"]["summaries"]
    frozen = data["frozen"]
    lhs_n = sum(1 for r in data["candidates"] if r["candidate_id"].startswith("lhs_"))
    later_n = counts.get("total", 0) - lhs_n

    boxes = [
        (0.3, 7.6, "① 预检查 preflight\n10 项数值自检\n(静态极限/功-能/步长收敛…)",
         f"全部通过\n耗时 {data['preflight'].get('elapsed_s', 0):.0f} s"),
        (2.75, 7.6, "② 阶段1 候选探索\noptimization_pool 128 shots\n拉丁超立方 + dt=0.10 μs",
         f"LHS 探索 {lhs_n} 候选\n后续入围/精修/样条 {max(later_n, 0)} 候选\n"
         f"（有效 {counts.get('ok', '?')} / 无效 {counts.get('invalid', '?')}"
         f" / 失败 {counts.get('failed', '?')}）"),
        (5.2, 7.6, "③ 阶段2 精修与选择\nselection_pool 512 shots\ndt=0.05 μs + PCHIP 样条",
         "唯一最终波形冻结\n→ best_waveform.yaml"),
        (7.65, 7.6, "④ 独立验证 validation_pool\n2000 shots × 3 波形\n同一初态数组(CRN)",
         "俘获率 + Wilson 95% CI\n配对 bootstrap 比较"),
        (7.65, 4.4, "⑤ 鲁棒性检查\n温度 3/5/7/10 μK\n对准 ±0.10 μm\n步长 0.10/0.05/0.025 μs",
         "冻结参数后不再调优"),
    ]
    for x, y, title, sub in boxes:
        ax.add_patch(FancyBboxPatch((x, y), 2.05, 1.9, boxstyle="round,pad=0.08",
                                    fc="#eaf3fb", ec="#1f77b4", lw=1.6))
        ax.text(x + 1.02, y + 1.52, title, ha="center", va="center", fontsize=10,
                fontweight="bold")
        ax.text(x + 1.02, y + 0.55, sub, ha="center", va="center", fontsize=9, color="#333")

    for x0, y0, x1, y1 in [(2.35, 8.55, 2.75, 8.55), (4.8, 8.55, 5.2, 8.55),
                           (7.25, 8.55, 7.65, 8.55), (8.67, 7.6, 8.67, 6.3)]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                     mutation_scale=22, lw=2, color="#555"))

    # 防过拟合三池隔离示意
    ax.add_patch(FancyBboxPatch((0.3, 4.4), 6.6, 2.4, boxstyle="round,pad=0.1",
                                fc="#fdf3e7", ec="#ff7f0e", lw=1.5))
    ax.text(3.6, 6.4, "防止蒙特卡洛过拟合：三池隔离 + 公共随机数(CRN)", ha="center",
            fontsize=11, fontweight="bold")
    for i, (pool, seed, shots, use) in enumerate([
            ("optimization_pool", 21001, 128, "只用于搜索候选"),
            ("selection_pool", 21002, 512, "只用于选择最终波形"),
            ("validation_pool", 21003, 2000, "只报告一次，禁止回调参数")]):
        x = 0.6 + i * 2.15
        color = "#8c8c8c" if i == 2 else "#2ca02c"
        ax.add_patch(FancyBboxPatch((x, 4.65), 1.95, 1.2, boxstyle="round,pad=0.06",
                                    fc="white", ec=color, lw=1.5))
        ax.text(x + 0.97, 5.55, pool, ha="center", fontsize=9, fontweight="bold", color=color)
        ax.text(x + 0.97, 5.1, f"seed={seed}\n{shots} shots\n{use}", ha="center",
                va="center", fontsize=8)

    # 冻结参数
    ax.add_patch(FancyBboxPatch((0.3, 1.4), 4.6, 2.3, boxstyle="round,pad=0.1",
                                fc="#eefaee", ec="#2ca02c", lw=1.6))
    ax.text(2.6, 3.35, "冻结的优化波形参数（5 维受约束族）", ha="center", fontsize=11,
            fontweight="bold")
    ax.text(2.6, 2.35,
            f"移动窗口 s ∈ [{frozen['move_start_fraction']:.3f}, "
            f"{frozen['move_end_fraction']:.3f}]\n"
            f"降深窗口 s ∈ [{frozen['ramp_start_fraction']:.3f}, "
            f"{frozen['ramp_end_fraction']:.3f}]\n"
            f"降深幂次 p = {frozen['ramp_power']:.3f}"
            f"（选中候选 {frozen['provenance']['candidate_id']}）",
            ha="center", va="center", fontsize=9.5)

    # 验证结论
    ax.add_patch(FancyBboxPatch((5.3, 1.4), 4.4, 2.3, boxstyle="round,pad=0.1",
                                fc="#fdeaea", ec="#d62728", lw=1.6))
    ax.text(7.5, 3.35, "独立验证结论（2000 shots）", ha="center", fontsize=11,
            fontweight="bold")
    lines = []
    for name in WAVEFORMS:
        s = val[name]
        lines.append(f"{WF_LABEL[name]}: {s['captured']}/{s['shots']} "
                     f"= {s['capture_fraction']:.4f}")
    ax.text(7.5, 2.25, "\n".join(lines) + "\n\n三波形均 100%（Wilson CI 下限 0.9981）\n"
            "→ 处于饱和区，未观察到显著提升（如实报告）",
            ha="center", va="center", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(out / "01_流程总览.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 2：波形
def fig2_waveforms(data: dict, out: Path):
    w = data["waveforms"]
    t_us = w["time_us"]
    frozen = data["frozen"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle("图2 三类主波形：AOD 中心与深度随时间演化（含导数）", fontsize=14,
                 fontweight="bold")

    ax = axes[0, 0]
    for name in WAVEFORMS:
        ax.plot(t_us[t_us <= 600], w[f"{name}_aod_center_um"][t_us <= 600],
                color=WF_COLOR[name], label=WF_LABEL[name], lw=2)
    ax.axvspan(frozen["move_start_fraction"] * 400, frozen["move_end_fraction"] * 400,
               color="#1f77b4", alpha=0.06, label="移动窗口")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("AOD 中心 (μm)")
    ax.set_title("AOD 中心位置 c(t)：从 2.4 μm 移动到 SLM 中心 0")
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    for name in WAVEFORMS:
        ax.plot(t_us[t_us <= 600], w[f"{name}_aod_depth_uK"][t_us <= 600],
                color=WF_COLOR[name], label=WF_LABEL[name], lw=2)
    ax.axvspan(frozen["ramp_start_fraction"] * 400, frozen["ramp_end_fraction"] * 400,
               color="#9467bd", alpha=0.08, label="降深窗口")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("AOD 深度 (μK)")
    ax.set_title("AOD 深度 D(t)：从 280 μK 降到 0")
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    m = t_us <= 400
    for name in WAVEFORMS:
        ax.plot(t_us[m], w[f"{name}_center_velocity_m_per_s"][m] * 1e3,
                color=WF_COLOR[name], label=WF_LABEL[name], lw=1.8)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("中心速度 (mm/s)")
    ax.set_title("中心一阶导数 ċ(t)（400 μs 段）：端点为零，无猛启停")
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    for name in WAVEFORMS:
        ax.plot(t_us[m], w[f"{name}_depth_rate_J_per_s"][m] / KB / 1e6,
                color=WF_COLOR[name], label=WF_LABEL[name], lw=1.8)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("深度变化率 (μK/s)")
    ax.legend(fontsize=9)

    ov = max(0.0, min(frozen["move_end_fraction"], frozen["ramp_end_fraction"])
             - max(frozen["move_start_fraction"], frozen["ramp_start_fraction"]))
    ov_t0 = max(frozen["move_start_fraction"], frozen["ramp_start_fraction"]) * 400
    ov_t1 = min(frozen["move_end_fraction"], frozen["ramp_end_fraction"]) * 400
    for ax in axes.flat:
        ax.axvspan(ov_t0, ov_t1, color="#d62728", alpha=0.07, zorder=0)
    axes[0, 0].text((ov_t0 + ov_t1) / 2, 1.5, "重叠区", ha="center", fontsize=8,
                    color="#d62728")
    axes[0, 0].text(120, 1.5, "← 仅移动", ha="center", fontsize=8, color="#666")
    axes[1, 1].set_title("深度一阶导数 Ḋ(t)：浅红竖带为移动/降深重叠区")
    fig.text(0.5, 0.005,
             f"优化波形关键特征：移动窗口几乎全程（s∈[{
                 frozen['move_start_fraction']:.3f}, 1]），降深在 s≥"
             f"{frozen['ramp_start_fraction']:.3f} 才开始 —— 移动与降深重叠约 "
             f"{ov * 100:.0f}% 的时间，这就是『同步波形』",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out / "02_波形对比.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 3：势能演化
def fig3_potentials(data: dict, out: Path):
    w = data["waveforms"]
    t_us = w["time_us"]
    name = "overlapped_optimized_400us"
    times = [0, 80, 160, 240, 300, 360, 400]
    x = np.linspace(-4, 5, 800)
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), sharey=True)
    fig.suptitle("图3 优化同步波形下总势场 U(x,t)=U_SLM+U_AOD 的演化（颜色越深势越深）",
                 fontsize=14, fontweight="bold")
    for ax, t in zip(axes.flat, times):
        idx = int(np.argmin(np.abs(t_us - t)))
        c = w[f"{name}_aod_center_um"][idx]
        d = w[f"{name}_aod_depth_uK"][idx]
        u_slm = gaussian_well(x, D_SLM_UK, W_SLM_UM, 0.0)
        u_aod = gaussian_well(x, d, W_AOD_UM, c) if d > 1e-9 else np.zeros_like(x)
        ax.plot(x, u_slm, "--", color="#7f7f7f", lw=1, label="SLM (140 μK)")
        if d > 1e-9:
            ax.plot(x, u_aod, ":", color="#d62728", lw=1.2, label=f"AOD ({d:.0f} μK)")
        ax.plot(x, u_slm + u_aod, color="#1f77b4", lw=2.2, label="总势")
        ax.axvline(c, color="#d62728", lw=0.8, alpha=0.5)
        ax.axvline(0, color="#7f7f7f", lw=0.8, alpha=0.5)
        ax.set_title(f"t = {t} μs   (AOD 中心 {c:+.2f} μm, 深度 {d:.0f} μK)",
                     fontsize=10)
        ax.set_xlabel("位置 (μm)")
    axes[0, 0].set_ylabel("势能 (μK)")
    axes[1, 0].set_ylabel("势能 (μK)")
    axes[0, 0].legend(fontsize=8)
    fig.text(0.5, 0.01, "原子初始呆在右侧红色 AOD 阱里；AOD 阱边移动边变浅，原子被逐步"
             "『倒』进左侧静止的灰色 SLM 阱", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out / "03_势能剖面演化.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 4：优化过程
def fig4_optimization(data: dict, out: Path):
    cand = data["candidates"]
    history = data["history"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9.5))
    fig.suptitle("图4 优化过程：候选探索 → 精修 → 样条 → 冻结（全部候选均保留）",
                 fontsize=14, fontweight="bold")

    ax = axes[0, 0]
    status_color = {"ok": "#1f77b4", "invalid": "#d62728", "failed": "k"}
    order = {"lhs_": 0, "sel_": 1, "ref_": 1, "spl_": 2}
    groups = {}
    for i, row in enumerate(cand):
        key = row["candidate_id"].split("_")[0] + "_"
        groups.setdefault(key, ([], []))
        groups[key][0].append(i)
        groups[key][1].append(to_float(row["objective_total"]))
    marker = {"lhs_": "o", "sel_": "s", "ref_": "s", "spl_": "^"}
    label = {"lhs_": "阶段1 LHS 探索", "sel_": "阶段2 入围复评", "ref_": "阶段2 局部精修",
             "spl_": "样条精修"}
    for key, (idxs, scores) in groups.items():
        ax.scatter(idxs, scores, s=26, marker=marker[key],
                   c=[status_color[cand[i]["status"]] for i in idxs], label=label[key])
    best = [r for r in cand if r["is_selected"] == "True"]
    if best:
        i_best = cand.index(best[0])
        ax.scatter([i_best], [to_float(best[0]["objective_total"])], marker="*",
                   s=350, color="#ffd700", edgecolor="k", zorder=5, label="最终选中")
    ax.set_xlabel("候选评估顺序")
    ax.set_ylabel("目标总分（俘获率为主 + 激发 + 光滑度惩罚）")
    ax.set_title("每个候选的得分（蓝=有效 红=违反约束 黑=运行失败）")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    scores = np.array([to_float(r["objective_total"]) for r in history
                       if r["status"] == "ok"])
    if len(scores):
        ax.plot(np.arange(1, len(scores) + 1), np.maximum.accumulate(scores),
                color="#2ca02c", lw=2, label="历史最优（单调不降）")
        ax.plot(np.arange(1, len(scores) + 1), scores, ".", color="#1f77b4",
                ms=3, alpha=0.5, label="单次评估")
    ax.set_xlabel("评估序号")
    ax.set_ylabel("目标总分")
    ax.set_title("最优值演化：先全局探索（LHS），后局部精修收敛")
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    ok = [r for r in cand if r["status"] == "ok" and r["candidate_id"].startswith("lhs_")]
    ms = [to_float(r["move_start_fraction"]) for r in ok]
    rs = [to_float(r["ramp_start_fraction"]) for r in ok]
    sc = [to_float(r["objective_total"]) for r in ok]
    p = ax.scatter(ms, rs, c=sc, cmap="viridis", s=40, edgecolor="k", lw=0.3)
    if best and best[0]["candidate_id"].startswith(("lhs", "sel", "ref")):
        ax.scatter([to_float(best[0]["move_start_fraction"])],
                   [to_float(best[0]["ramp_start_fraction"])], marker="*", s=350,
                   color="#ffd700", edgecolor="k", zorder=5)
    ax.set_xlabel("move_start_fraction（移动开始时刻 / T）")
    ax.set_ylabel("ramp_start_fraction（降深开始时刻 / T）")
    ax.set_title("阶段1 候选参数空间（颜色=得分）：越晚降深越好")
    plt.colorbar(p, ax=ax, label="目标总分")

    ax = axes[1, 1]
    ax.hist([to_float(r["ramp_power"]) for r in ok], bins=25,
            color="#9467bd", edgecolor="white")
    if best:
        ax.axvline(to_float(best[0]["ramp_power"]), color="#ffd700", lw=2.5,
                   label=f"选中 p={to_float(best[0]['ramp_power']):.2f}")
        ax.legend(fontsize=9)
    ax.set_xlabel("ramp_power（降深幂次）")
    ax.set_ylabel("候选数")
    ax.set_title("降深幂次分布：搜索覆盖 [0.5, 4]，最终选中接近线性")
    fig.tight_layout()
    fig.savefig(out / "04_优化过程.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 5：验证结果
def fig5_validation(data: dict, out: Path):
    shots = data["validation"]
    metrics = data["metrics"]
    fig = plt.figure(figsize=(16, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.28)
    fig.suptitle("图5 独立验证（validation_pool, 2000 shots, 三波形共用同一初态）",
                 fontsize=14, fontweight="bold")

    # (a) 俘获率 + Wilson CI：三波形均 1.0，改为展示更有信息量的 CI 下限柱
    ax = fig.add_subplot(gs[0, 0])
    for i, name in enumerate(WAVEFORMS):
        s = metrics["validation"]["summaries"][name]
        lo = s["wilson_low"]
        ax.bar(i, 1.0 - lo, bottom=lo, width=0.55, color=WF_COLOR[name], alpha=0.85)
        ax.annotate(f"{s['captured']}/{s['shots']}\nCI 下限 {lo:.4f}", (i, lo),
                    textcoords="offset points", xytext=(0, 8), fontsize=8,
                    ha="center", color="#333")
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.text(2.42, 1.0, " 1.0", va="center", fontsize=8, color="gray")
    ax.set_xticks(range(3))
    ax.set_xticklabels([WF_LABEL[n] for n in WAVEFORMS], rotation=12, fontsize=8.5)
    ax.set_ylim(0.9975, 1.002)
    ax.set_ylabel("俘获率（柱体 = Wilson 95% CI 范围）")
    ax.set_title("(a) 俘获率 1.0000，最坏 CI 下限 0.9981：三波形无显著差异")

    # (b) 配对计数
    ax = fig.add_subplot(gs[0, 1])
    comps = metrics["validation"]["comparisons"]
    labels, both, only_a, only_b = [], [], [], []
    for key, c in comps.items():
        labels.append(key.replace("_vs_", "\nvs\n").replace("sequential_", "seq_")
                      .replace("overlapped_optimized_400us", "optimized"))
        pc = c["paired_counts"]
        both.append(pc["both_captured"])
        only_a.append(pc["only_a_captured"])
        only_b.append(pc["only_b_captured"])
    x = np.arange(len(labels))
    ax.bar(x, both, 0.55, label="都成功", color="#2ca02c")
    ax.bar(x, only_a, 0.55, bottom=both, label="仅前者成功", color="#ff7f0e")
    ax.bar(x, only_b, 0.55, bottom=np.array(both) + np.array(only_a),
           label="仅后者成功", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("shots 数")
    ax.set_title("(b) 配对计数：2000 对初态全部『都成功』，差=0")
    ax.legend(fontsize=8)

    # (c) 俘获裕量分布
    ax = fig.add_subplot(gs[0, 2])
    for name in WAVEFORMS:
        m = [to_float(r["capture_margin"]) for r in shots if r["waveform"] == name]
        ax.hist(m, bins=40, alpha=0.45, color=WF_COLOR[name], label=WF_LABEL[name],
                density=True)
    ax.set_xlabel("归一化俘获裕量 M = -E_f / D_SLM（>0 即被俘获）")
    ax.set_ylabel("概率密度")
    ax.set_title("(c) 俘获裕量分布：三波形几乎重合")
    ax.legend(fontsize=8)

    # (d) 末态相空间
    ax = fig.add_subplot(gs[1, 0])
    for name, alpha in zip(WAVEFORMS, (0.25, 0.25, 0.6)):
        xs = np.array([to_float(r["final_x_um"]) for r in shots if r["waveform"] == name])
        vs = np.array([to_float(r["final_v_m_s"]) for r in shots if r["waveform"] == name])
        cap = np.array([r["captured"] == "True" for r in shots if r["waveform"] == name])
        ax.scatter(xs[cap], vs[cap] * 1e3, s=6, alpha=alpha, color=WF_COLOR[name])
        if not cap.all():
            ax.scatter(xs[~cap], vs[~cap] * 1e3, s=30, marker="x", color="k")
    mass_kg = MASS_U * 1.66053906660e-27
    d_j = D_SLM_UK * KB * 1e-6
    xline = np.linspace(-3, 3, 400)
    v_sep = np.sqrt(np.maximum(2 * (d_j + gaussian_well(xline, D_SLM_UK, W_SLM_UM, 0)
                                    * KB * 1e-6) / mass_kg, 0)) * 1e3
    ax.plot(xline, v_sep, "k--", lw=1.2, label="E_SLM=0 分界线")
    ax.plot(xline, -v_sep, "k--", lw=1.2)
    ax.set_xlabel("末态位置 (μm)")
    ax.set_ylabel("末态速度 (mm/s)")
    ax.set_title("(d) 末态相空间：全部点位于分界线内（全部俘获）")
    ax.legend(fontsize=8)

    # (e) 功-能核算
    ax = fig.add_subplot(gs[1, 1])
    names, res_med, res_max, hold = [], [], [], []
    for name in WAVEFORMS:
        rows = [r for r in shots if r["waveform"] == name]
        res = [to_float(r["max_abs_work_energy_residual_over_depth"]) for r in rows]
        names.append(name)
        res_med.append(np.nanmedian(res))
        res_max.append(np.nanmax(res))
        hold.append(np.nanmax([to_float(r["hold_peak_to_peak_energy_error_over_depth"])
                               for r in rows]))
    x = np.arange(3)
    ax.bar(x - 0.22, res_med, 0.2, label="功-能残差 中位数", color="#1f77b4")
    ax.bar(x, res_max, 0.2, label="功-能残差 最大值", color="#ff7f0e")
    ax.bar(x + 0.22, hold, 0.2, label="100μs 保持段能量误差 峰峰值", color="#2ca02c")
    ax.axhline(1e-3, color="r", ls="--", lw=1, label="验收阈值 1e-3")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([WF_LABEL[n] for n in names], rotation=12, fontsize=8.5)
    ax.set_ylabel("归一化误差（相对阱深）")
    ax.set_title("(e) 功-能平衡与保持段：远低于验收阈值")
    ax.legend(fontsize=7.5)

    # (f) 激发与外部功
    ax = fig.add_subplot(gs[1, 2])
    for name in WAVEFORMS:
        rows = [r for r in shots if r["waveform"] == name]
        ef = [to_float(r["final_excitation_over_depth"]) for r in rows]
        ax.hist(ef, bins=40, alpha=0.45, color=WF_COLOR[name], label=WF_LABEL[name],
                density=True)
    ax.set_xlabel("归一化末态激发 e_f（0=SLM 阱底, 1=逸出阈值）")
    ax.set_ylabel("概率密度")
    ax.set_title("(f) 末态激发分布：均值≈0.025，仅达阱深 2.5%")
    ax.legend(fontsize=8)
    fig.savefig(out / "05_独立验证结果.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 6：初态与轨迹
def fig6_initial_states(data: dict, out: Path):
    shots = data["validation"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    fig.suptitle("图6 蒙特卡洛初态抽样与末态（简谐提议 + 束缚拒绝，5 μK，2000 shots）",
                 fontsize=14, fontweight="bold")

    ax = axes[0]
    x0 = np.array([to_float(r["x0_um"]) for r in shots if r["waveform"] == WAVEFORMS[0]])
    v0 = np.array([to_float(r["v0_m_s"]) for r in shots if r["waveform"] == WAVEFORMS[0]])
    ax.scatter(x0, v0 * 1e3, s=4, alpha=0.4, color="#1f77b4")
    ax.axvline(D0_UM, color="#d62728", ls="--", lw=1, label="AOD 中心 2.4 μm")
    ax.set_xlabel("初态位置 x₀ (μm)")
    ax.set_ylabel("初态速度 v₀ (mm/s)")
    ax.set_title(f"初态相空间（n={len(x0)}，验收率见 metrics）")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.hist(x0, bins=50, color="#1f77b4", edgecolor="white")
    ax.axvline(D0_UM, color="#d62728", ls="--", lw=1.5)
    ax.set_xlabel("初态位置 x₀ (μm)")
    ax.set_ylabel("样本数")
    ax.set_title(f"初态位置分布：σ≈{np.std(x0):.3f} μm（热分布于 AOD 阱底附近）")

    ax = axes[2]
    for name in WAVEFORMS:
        xf = [to_float(r["final_x_um"]) for r in shots if r["waveform"] == name]
        ax.hist(xf, bins=50, alpha=0.45, color=WF_COLOR[name], label=WF_LABEL[name],
                density=True)
    ax.set_xlabel("末态位置 x_f (μm)")
    ax.set_ylabel("概率密度")
    ax.set_title("末态位置分布：集中在 SLM 阱（0 μm）内")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "06_初态与末态分布.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 7：鲁棒性
def fig7_robustness(data: dict, out: Path):
    rows = data["robustness"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    fig.suptitle("图7 鲁棒性敏感性（冻结参数后不再调优）：温度 / 末端对准 / 时间步长",
                 fontsize=14, fontweight="bold")

    for ax, scan, xlabel in [(axes[0], "temperature", "初态温度 (μK)"),
                             (axes[1], "alignment", "末端对准偏差 (μm)")]:
        for name in WAVEFORMS:
            pts = [(to_float(r["condition"]), to_float(r["capture_fraction"]),
                    to_float(r["wilson_low"]), to_float(r["wilson_high"]))
                   for r in rows if r["scan"] == scan and r["waveform"] == name]
            if not pts:
                continue
            pts.sort()
            c = [p[0] for p in pts]
            f = [p[1] for p in pts]
            lo = [p[2] for p in pts]
            hi = [p[3] for p in pts]
            ax.errorbar(c, f, yerr=[np.maximum(0, np.array(f) - np.array(lo)),
                                    np.maximum(0, np.array(hi) - np.array(f))],
                        fmt="-o", ms=6, capsize=4, color=WF_COLOR[name], label=WF_LABEL[name])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("俘获率")
        lo_all = min(to_float(r["wilson_low"]) for r in rows
                     if r["scan"] == scan and r.get("wilson_low"))
        ax.set_ylim(min(0.95, lo_all - 0.01), 1.005)
        ax.set_title({"temperature": "(a) 温度扫描：3–10 μK 每点 500 shots",
                      "alignment": "(b) 对准偏差：±0.10 μm 每点 500 shots"}[scan])
        ax.legend(fontsize=8)

    ax = axes[2]
    width = 0.25
    for i, name in enumerate(WAVEFORMS):
        pts = [(to_float(r["condition"]), to_float(r["final_energy_max_abs_diff_uK"]),
                to_float(r["label_flips_vs_reference"]))
               for r in rows if r["scan"] == "timestep" and r["waveform"] == name
               and r.get("final_energy_max_abs_diff_uK")]
        if not pts:
            continue
        pts.sort()
        dts = [p[0] for p in pts]
        diff = [p[1] for p in pts]
        flips = [p[2] for p in pts]
        bars = ax.bar(np.arange(len(dts)) + (i - 1) * width, diff, width,
                      color=WF_COLOR[name], label=WF_LABEL[name])
        for bar, flip in zip(bars, flips):
            if flip == flip and flip:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f"翻转{int(flip)}", ha="center", va="bottom", fontsize=8,
                        color="#d62728")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["dt=0.10 μs", "dt=0.05 μs", "dt=0.025 μs"])
    ax.set_ylabel("末态能量最大差 (μK, 相对最细步长)")
    ax.set_title("(c) 步长收敛：连续量随步长减小而收敛，俘获标签翻转≈0")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "07_鲁棒性.png")
    plt.close(fig)


# ---------------------------------------------------------------- 图 8：预检
def fig8_preflight(data: dict, out: Path):
    pf = data["preflight"]
    fig, ax = plt.subplots(figsize=(13, 7.5))
    ax.axis("off")
    ax.set_title("图8 Preflight 预检查：10 项数值自检，全部通过才允许开始优化",
                 fontsize=14, fontweight="bold")
    # 键名必须与 preflight.json 实际字段一致；waveform_constraints / ideal_runs
    # 无顶层 passed 字段，需从子项汇总。
    describe = {
        "existing_tests": ("existing_tests", "运行现有 Level 0/1 测试套件"),
        "static_limit": ("static_limit", "时变积分器在静态极限下与 Level 0 一致"),
        "level1_reproduction": ("level1_reproduction", "固定种子复现 Level 1 sequential_600us"),
        "waveform_constraints": ("waveform_constraints", "三类波形端点/范围/单调性/导数检查"),
        "ideal_runs": ("ideal_runs", "理想初态 (x₀=c_AOD, v₀=0) 下运行三波形"),
        "ideal_work_energy": ("ideal_work_energy", "理想轨迹功-能平衡残差 < 5e-4"),
        "ideal_timestep_convergence": ("ideal_timestep_convergence", "步长 0.05→0.025 μs 二阶收敛"),
        "mc_timestep_paired": ("mc_timestep_paired", "128 shots 配对步长检查（标签不一致率 ≤ 2%）"),
        "mc_work_energy_residual": ("mc_work_energy_residual", "MC 子集功-能残差中位数 < 1e-3"),
        "hold_segment": ("hold_segment", "100 μs 纯 SLM 保持段能量误差有界"),
    }
    rows = []
    for key, desc in describe.values():
        entry = pf.get(key, {})
        if isinstance(entry, bool):
            passed, detail = entry, ""
        elif key == "waveform_constraints":
            passed = bool(pf.get("waveform_constraints_all_valid",
                                 all(v.get("valid", False) for v in entry.values())))
            detail = "三波形均满足约束"
        elif key == "ideal_runs":
            passed = all(v.get("captured", False) for v in entry.values())
            detail = "三波形理想轨迹均被 SLM 俘获"
        else:
            passed = bool(entry.get("passed", False))
            detail = entry.get("summary_line", "")
        rows.append((desc, passed, detail))
    y = 0.92
    for desc, passed, detail in rows:
        color = "#2ca02c" if passed else "#d62728"
        ax.text(0.02, y, "✓" if passed else "✗", fontsize=16, color=color,
                fontweight="bold", va="center", ha="center")
        ax.text(0.07, y, desc, fontsize=11, va="center")
        if detail:
            ax.text(0.62, y, detail, fontsize=9, va="center", color="#666")
        y -= 0.085
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.text(0.5, 0.04,
             f"总耗时 {pf.get('elapsed_s', 0):.0f} s —— 任何一项失败都会终止流程，"
             f"禁止生成『看似完整』的最优结果", ha="center", fontsize=10)
    fig.savefig(out / "08_预检查清单.png")
    plt.close(fig)


def main():
    repo = Path(__file__).resolve().parent
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        repo / "simulation" / "level2_joint_transfer" / "outputs" /
        "level2_joint_transfer" / "demo")
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else repo / "level2_figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    data = load_all(data_dir)
    figures = [
        ("01_流程总览.png", fig1_pipeline),
        ("02_波形对比.png", fig2_waveforms),
        ("03_势能剖面演化.png", fig3_potentials),
        ("04_优化过程.png", fig4_optimization),
        ("05_独立验证结果.png", fig5_validation),
        ("06_初态与末态分布.png", fig6_initial_states),
        ("07_鲁棒性.png", fig7_robustness),
        ("08_预检查清单.png", fig8_preflight),
    ]
    for name, func in figures:
        func(data, out_dir)
        print(f"已生成 {out_dir / name}")
    print(f"\n共 {len(figures)} 张图，数据来自 {data_dir}")


if __name__ == "__main__":
    main()
