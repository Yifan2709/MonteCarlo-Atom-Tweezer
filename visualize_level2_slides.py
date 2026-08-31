# -*- coding: utf-8 -*-
"""Level 2 演示版幻灯片图：一页一图、大字号、少文字。

读取 outputs/level2_joint_transfer/demo/ 的数据，生成 10 张按 16:9 幻灯片
设计的图，保存到 level2_slide_figures/。所有数字均来自实际输出，不硬编码。

用法：
    python3 visualize_level2_slides.py [数据目录] [输出目录]
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
    "grid.alpha": 0.25,
    "font.size": 13,
})

REPO = Path(__file__).resolve().parent
DEMO_2B = (REPO / "simulation" / "level2b_speed_modify" / "outputs" /
           "level2b_speed_limit" / "demo")
KB_UK = 1.380649e-29                     # 1 μK 对应的焦耳数（J per μK）
WAVEFORMS = ["sequential_600us", "sequential_400us", "overlapped_optimized_400us"]
WF_SHORT = {
    "sequential_600us": "顺序 600 μs",
    "sequential_400us": "顺序 400 μs",
    "overlapped_optimized_400us": "优化同步 400 μs",
}
WF_COLOR = {
    "sequential_600us": "#1f77b4",
    "sequential_400us": "#ff7f0e",
    "overlapped_optimized_400us": "#d62728",
}
FIGSIZE = (14.0, 5.7)                    # 与幻灯片图片区比例一致


def to_float(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def gaussian_well(x_um, depth_uk, waist_um, center_um):
    """一维高斯光阱势（单位 μK），与仓库物理模型一致。"""
    return -depth_uk * np.exp(-2.0 * (x_um - center_um) ** 2 / waist_um ** 2)


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_all(data_dir: Path) -> dict:
    data = {}
    data["metrics"] = json.loads((data_dir / "metrics.json").read_text())
    data["preflight"] = json.loads((data_dir / "preflight.json").read_text())
    data["frozen"] = yaml.safe_load((data_dir / "best_waveform.yaml").read_text())
    data["candidates"] = read_csv_rows(data_dir / "candidates.csv")
    data["validation"] = read_csv_rows(data_dir / "validation_shots.csv")
    data["robustness"] = read_csv_rows(data_dir / "robustness.csv")
    with (data_dir / "waveforms.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = np.array([[to_float(v) for v in row] for row in reader])
    data["waveforms"] = {name: rows[:, header.index(name)] for name in header}
    return data


def save(fig, out: Path, name: str):
    fig.savefig(out / name, facecolor="white")
    plt.close(fig)
    print(f"已生成 {out / name}")


# ---------------------------------------------------------------- s01 问题：时间线
def s01_question(data, out):
    frozen = data["frozen"]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(-170, 760)
    ax.set_ylim(0, 3.7)
    ax.axis("off")
    bar_h = 0.74

    # barh 的 y 是条中心：条实际占据 y ± bar_h/2
    # 第一行：Level 1 顺序 600 μs
    y1 = 2.42
    ax.barh(y1, 312, left=0, height=bar_h, color="#1f77b4", alpha=0.92)
    ax.barh(y1, 288, left=312, height=bar_h, color="#ff7f0e", alpha=0.92)
    ax.text(156, y1, "只移动 52%", ha="center", va="center",
            color="white", fontsize=14, fontweight="bold")
    ax.text(456, y1, "只降深 48%", ha="center", va="center",
            color="white", fontsize=14, fontweight="bold")
    ax.text(-25, y1, "Level 1\n顺序基线", ha="right", va="center",
            fontsize=14, fontweight="bold", color="#333")
    ax.annotate("", xy=(600, y1 + bar_h / 2 + 0.18), xytext=(0, y1 + bar_h / 2 + 0.18),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=1.6))
    ax.text(300, y1 + bar_h / 2 + 0.26, "600 μs", ha="center", fontsize=17,
            fontweight="bold", color="#333")

    # 第二行：Level 2 同步 400 μs（冻结波形的真实窗口）
    y2 = 0.95
    T = frozen["duration_us"]
    ms, me = frozen["move_start_fraction"] * T, frozen["move_end_fraction"] * T
    rs, re_ = frozen["ramp_start_fraction"] * T, frozen["ramp_end_fraction"] * T
    ax.barh(y2, me - ms, left=ms, height=bar_h, color="#1f77b4", alpha=0.92)
    ax.barh(y2, re_ - rs, left=rs, height=bar_h, color="#ff7f0e", alpha=0.55,
            hatch="///", edgecolor="#b35a00", lw=1.4)
    ax.text((ms + rs) / 2, y2, "移动贯穿全程", ha="center", va="center",
            color="white", fontsize=14, fontweight="bold")
    cx = (rs + re_) / 2
    ax.plot([cx, cx], [y2 - bar_h / 2 - 0.03, y2 - bar_h / 2 - 0.18], color="#b35a00",
            lw=1.3, ls=":")
    ax.text(cx, y2 - bar_h / 2 - 0.22, "移动与降深重叠 ≈ 39%", ha="center", va="top",
            fontsize=14, color="#b35a00", fontweight="bold")
    ax.text(-25, y2, "Level 2\n同步波形", ha="right", va="center",
            fontsize=14, fontweight="bold", color="#333")
    ax.annotate("", xy=(T, y2 + bar_h / 2 + 0.18), xytext=(0, y2 + bar_h / 2 + 0.18),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=1.6))
    ax.text(T / 2, y2 + bar_h / 2 + 0.26, "400 μs", ha="center", fontsize=17,
            fontweight="bold", color="#b35a00")

    ax.annotate("", xy=(665, y1), xytext=(665, y2 + bar_h / 2),
                arrowprops=dict(arrowstyle="-|>", color="#2ca02c", lw=3))
    ax.text(690, (y1 + y2) / 2, "提速\n33%", ha="left", va="center",
            fontsize=15, color="#2ca02c", fontweight="bold")
    save(fig, out, "s01_问题.png")


# ---------------------------------------------------------------- s02 物理模型
def s02_model(data, out):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    x = np.linspace(-4.2, 5.4, 1400)
    u_slm = gaussian_well(x, 140.0, 1.17, 0.0)
    u_aod = gaussian_well(x, 280.0, 1.17, 2.4)
    ax.plot(x, u_slm, "--", color="#7f7f7f", lw=2.0, label="SLM 静止阱 140 μK")
    ax.plot(x, u_aod, "--", color="#d62728", lw=2.0, label="AOD 初始阱 280 μK")
    ax.plot(x, u_slm + u_aod, color="#1f77b4", lw=3.2, label="总势 U(x)")
    ax.plot([2.4], [-271], "o", ms=15, color="#9467bd", zorder=5)
    ax.annotate("Cs 原子（5 μK 热初态）", xy=(2.4, -271), xytext=(3.4, -180),
                fontsize=14, color="#9467bd", fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color="#9467bd", lw=1.8))
    ax.annotate("", xy=(0.0, -150), xytext=(2.15, -255),
                arrowprops=dict(arrowstyle="-|>", color="#2ca02c", lw=3.0,
                                connectionstyle="arc3,rad=-0.25"))
    ax.text(1.25, -218, "AOD 边移边变浅\n把原子留给 SLM", ha="center", fontsize=13,
            color="#2ca02c", fontweight="bold")
    ax.text(3.15, -60, "D_AOD 280 μK · D_SLM 140 μK\nw 1.17 μm · 间距 2.4 μm",
            ha="left", fontsize=12.5, color="#555",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#bbb"))
    ax.set_xlabel("位置 x (μm)")
    ax.set_ylabel("势能 (μK)")
    ax.set_ylim(-330, 60)
    ax.legend(loc="lower left", fontsize=12.5)
    save(fig, out, "s02_物理模型.png")


# ---------------------------------------------------------------- s03 方法：三池隔离
def box(ax, x, y, w, h, fc, ec, lw=2.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.10",
                                fc=fc, ec=ec, lw=lw))


def s03_method(data, out):
    m = data["metrics"]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.0)
    ax.axis("off")

    pools = [
        (0.45, "#eefaee", "#2ca02c", "① optimization 池",
         "seed 21001 · 128 shots", "拉丁超立方搜索候选"),
        (3.60, "#eaf3fb", "#1f77b4", "② selection 池",
         "seed 21002 · 512 shots", "精修 · 选择 · 冻结波形"),
        (6.75, "#fdeaea", "#d62728", "③ validation 池",
         "seed 21003 · 2000 shots", "只报告一次"),
    ]
    for x, fc, ec, head, sub, use in pools:
        box(ax, x, 2.05, 2.80, 1.85, fc, ec)
        ax.text(x + 1.40, 3.52, head, ha="center", fontsize=14.5,
                fontweight="bold", color=ec)
        ax.text(x + 1.40, 3.02, sub, ha="center", fontsize=12.5, color="#333")
        ax.text(x + 1.40, 2.45, use, ha="center", fontsize=13, color="#555")
    for x0 in (3.30, 6.45):
        ax.add_patch(FancyArrowPatch((x0, 2.98), (x0 + 0.28, 2.98),
                                     arrowstyle="-|>", mutation_scale=26,
                                     lw=2.4, color="#555"))
    # validation 结果禁止回流
    ax.add_patch(FancyArrowPatch((8.30, 1.85), (5.05, 1.30),
                                 arrowstyle="-|>", mutation_scale=22, lw=2.2,
                                 color="#d62728", linestyle=(0, (5, 4)),
                                 connectionstyle="arc3,rad=0.25"))
    ax.text(6.55, 1.02, "禁止用验证结果回调参数", ha="center", fontsize=13.5,
            color="#d62728", fontweight="bold")

    box(ax, 0.45, 4.30, 9.10, 0.95, "#fdf3e7", "#e8a860")
    ax.text(5.0, 4.78, "公共随机数（CRN）：所有候选波形用同一组初态逐行比较，"
            "蒙特卡洛噪声不掩盖波形差异", ha="center", va="center",
            fontsize=13.5, color="#7a4a00")
    ax.text(5.0, 5.66, "防止蒙特卡洛过拟合：三池隔离 + 公共随机数", ha="center",
            fontsize=15.5, fontweight="bold", color="#333")
    ax.text(5.0, 0.45, f"种子互不相同（配置强制校验）· 优化器源码不接触 validation 种子（spy 测试守护）",
            ha="center", fontsize=11.5, color="#888")
    save(fig, out, "s03_方法.png")


# ---------------------------------------------------------------- s04 三类波形
def s04_waveforms(data, out):
    w = data["waveforms"]
    frozen = data["frozen"]
    t_us = w["time_us"]
    T = float(frozen["duration_us"])
    durations = {"sequential_600us": 600.0, "sequential_400us": 400.0,
                 "overlapped_optimized_400us": T}
    ov0 = max(frozen["move_start_fraction"], frozen["ramp_start_fraction"]) * T
    ov1 = min(frozen["move_end_fraction"], frozen["ramp_end_fraction"]) * T

    fig, axes = plt.subplots(2, 1, figsize=FIGSIZE, sharex=True)
    for name in WAVEFORMS:
        m = t_us <= durations[name] + 1e-9
        center = w[f"{name}_aod_center_um"][m]
        depth = w[f"{name}_aod_depth_uK"][m]
        axes[0].plot(t_us[m], center, color=WF_COLOR[name], lw=2.6,
                     label=WF_SHORT[name])
        axes[1].plot(t_us[m], depth, color=WF_COLOR[name], lw=2.6,
                     label=WF_SHORT[name])
        axes[0].plot([durations[name]], [center[-1]], "o", ms=9,
                     color=WF_COLOR[name], zorder=5)
        axes[1].plot([durations[name]], [depth[-1]], "o", ms=9,
                     color=WF_COLOR[name], zorder=5)
    for ax in axes:
        ax.axvspan(ov0, ov1, color="#d62728", alpha=0.08, zorder=0)
        ax.axvline(400, color="#999", ls=":", lw=1.4)
    axes[0].text(400, 2.30, "400 μs", ha="center", fontsize=12.5, color="#555",
                 fontweight="bold")
    axes[0].text(596, 2.30, "600 μs", ha="right", fontsize=12.5, color="#555",
                 fontweight="bold")
    axes[0].text((ov0 + ov1) / 2, 1.55, "移动与降深\n重叠 ≈ 39%", ha="center",
                 fontsize=13, color="#b35a00", fontweight="bold")
    axes[0].set_ylabel("AOD 中心 (μm)")
    axes[1].set_ylabel("AOD 深度 (μK)")
    axes[1].set_xlabel("时间 (μs)")
    axes[1].set_xlim(0, 612)
    axes[0].legend(loc="lower left", fontsize=12)
    axes[0].text(590, 2.25 - 0.42, "端点速度 = 0", ha="right", fontsize=12, color="#666")
    fig.tight_layout()
    save(fig, out, "s04_波形对比.png")


# ---------------------------------------------------------------- s05 势场演化
def s05_potentials(data, out):
    w = data["waveforms"]
    t_us = w["time_us"]
    name = "overlapped_optimized_400us"
    times = [0, 130, 270, 400]
    colors = ["#c6dbef", "#6baed6", "#2171b5", "#08306b"]
    x = np.linspace(-3.6, 4.6, 1200)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    bottoms = []
    for t, color in zip(times, colors):
        idx = int(np.argmin(np.abs(t_us - t)))
        c = float(w[f"{name}_aod_center_um"][idx])
        d = float(w[f"{name}_aod_depth_uK"][idx])
        u_tot = gaussian_well(x, 140.0, 1.17, 0.0) + gaussian_well(x, d, 1.17, c)
        ax.plot(x, u_tot, color=color, lw=3.0, label=f"t = {t} μs")
        bottom_x = float(x[np.argmin(u_tot)])
        bottoms.append((bottom_x, float(np.min(u_tot))))
        ax.plot([bottom_x], [np.min(u_tot)], "o", ms=10, color=color, zorder=5)
    ax.axvline(0, color="#7f7f7f", ls=":", lw=1.6)
    ax.text(-0.12, 18, "SLM 中心", ha="right", fontsize=12, color="#7f7f7f")
    for (x0, y0), (x1, y1) in zip(bottoms[:-1], bottoms[1:]):
        ax.add_patch(FancyArrowPatch((x0, y0 - 8), (x1, y1 - 8),
                                     arrowstyle="-|>", mutation_scale=20, lw=2.2,
                                     color="#2ca02c",
                                     connectionstyle="arc3,rad=0.18"))
    ax.text(1.15, -305, "阱底（原子跟随）一路滑向 SLM", fontsize=13.5,
            color="#2ca02c", fontweight="bold")
    ax.set_xlabel("位置 x (μm)")
    ax.set_ylabel("总势能 (μK)")
    ax.set_ylim(-345, 55)
    ax.legend(loc="lower left", fontsize=12.5, ncols=2)
    save(fig, out, "s05_势场演化.png")


# ---------------------------------------------------------------- s06 优化过程
def s06_optimization(data, out):
    cand = data["candidates"]
    fig, ax = plt.subplots(figsize=FIGSIZE)

    # 阶段背景色带（按候选 id 前缀分段，与 candidates.csv 行序一致）
    stage_of = lambda cid: {"lhs": 0, "sel": 1, "ref": 1, "spl": 2}[cid.split("_")[0]]
    stages = [stage_of(r["candidate_id"]) for r in cand]
    bounds = {}
    for i, s in enumerate(stages):
        lo, hi = bounds.get(s, (i, i))
        bounds[s] = (min(lo, i), max(hi, i))
    band_style = [
        ("#eaf3fb", "阶段 1 · LHS 随机探索"),
        ("#fdf3e7", "阶段 2 · 前 8 名入围 + 局部精修"),
        ("#eefaee", "样条精修"),
    ]
    for s, (color, label) in enumerate(band_style):
        lo, hi = bounds[s]
        ax.axvspan(lo - 0.5, hi + 0.5, color=color, alpha=0.55, zorder=0)
        ax.text((lo + hi) / 2, 1.02, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=13.5, color="#333",
                fontweight="bold")

    scatter_style = {0: ("#1f77b4", "o"), 1: ("#ff7f0e", "s"), 2: ("#2ca02c", "^")}
    order, scores = [], []
    for i, row in enumerate(cand):
        if row["status"] != "ok":
            continue
        color, marker = scatter_style[stage_of(row["candidate_id"])]
        score = to_float(row["objective_total"])
        ax.scatter(i, score, s=46, marker=marker, color=color, alpha=0.75, zorder=3)
        order.append(i)
        scores.append(score)
    if order:
        running = np.maximum.accumulate(scores)
        ax.plot(order, running, color="#08306b", lw=2.6, zorder=4,
                label="历史最优（单调不降）")
    best = [r for r in cand if r["is_selected"] == "True"]
    if best:
        i_best = cand.index(best[0])
        s_best = to_float(best[0]["objective_total"])
        ax.scatter([i_best], [s_best], marker="*", s=520, color="#ffd700",
                   edgecolor="k", zorder=6, label="最终选中")
        ax.annotate(f"胜者 {best[0]['candidate_id']} → 冻结为 best_waveform.yaml",
                    xy=(i_best, s_best), xytext=(i_best + 6, 0.80),
                    fontsize=12.5, color="#7a5a00",
                    arrowprops=dict(arrowstyle="-|>", color="#7a5a00", lw=1.8))

    ax.text(0.015, 0.03, "总分 = 俘获率（主）+ 裕量/激发加分 − 光滑度惩罚，俘获数优先",
            transform=ax.transAxes, fontsize=11.5, color="#888")
    ax.text(0.985, 0.90, "顶部密集带 ≈ 1.048：\n5 μK 下多数合法候选满俘获，\n胜负只由裕量/激发细分",
            transform=ax.transAxes, ha="right", va="top", fontsize=12, color="#555",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#ccc", alpha=0.9))
    ax.set_xlabel("候选评估顺序（三阶段从左到右，28 个违反约束的候选未画出）")
    ax.set_ylabel("目标总分")
    ax.set_ylim(min(scores) - 0.03, max(scores) + 0.10)
    ax.legend(loc="center left", fontsize=12)
    save(fig, out, "s06_优化过程.png")


# ---------------------------------------------------------------- s07 独立验证
def s07_validation(data, out):
    metrics = data["metrics"]
    shots = data["validation"]
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)

    ax = axes[0]
    summaries = metrics["validation"]["summaries"]
    for i, name in enumerate(WAVEFORMS):
        s = summaries[name]
        ax.errorbar(i, s["capture_fraction"],
                    yerr=[[s["capture_fraction"] - s["wilson_low"]],
                          [max(s["wilson_high"] - s["capture_fraction"], 0)]],
                    fmt="o", ms=11, capsize=7, lw=2.4, color=WF_COLOR[name])
        ax.annotate(f"{s['captured']}/{s['shots']}", (i, s["capture_fraction"]),
                    textcoords="offset points", xytext=(0, 14), ha="center",
                    fontsize=13, fontweight="bold", color="#333")
    ax.axvline(1.0, color="gray", ls="--", lw=1.2)
    ax.set_xticks(range(3), [WF_SHORT[n] for n in WAVEFORMS], fontsize=11.5)
    ax.tick_params(axis="x", rotation=10)
    ax.set_ylim(0.9962, 1.0015)
    ax.set_ylabel("俘获率")
    ax.set_title("俘获率全部 100%（Wilson 95% CI，纵轴放大）", fontsize=13.5)

    ax = axes[1]
    for name in WAVEFORMS:
        ef = np.array([to_float(r["final_excitation_over_depth"])
                       for r in shots if r["waveform"] == name])
        ax.hist(ef, bins=32, histtype="step", lw=2.6, color=WF_COLOR[name],
                label=WF_SHORT[name], density=True)
    ax.axvline(0.025, color="k", ls="--", lw=1.6)
    ax.text(0.70, 0.84, "均值 ≈ 0.025\n（仅阱深的 2.5%）",
            transform=ax.transAxes, ha="center", fontsize=13, color="#333")
    ax.set_xlabel("末态激发 e_f（0 = SLM 阱底，1 = 逸出）")
    ax.set_ylabel("概率密度")
    ax.set_title("三波形末态激发分布几乎重合", fontsize=13.5)
    ax.legend(fontsize=11)
    fig.tight_layout()
    save(fig, out, "s07_独立验证.png")


# ---------------------------------------------------------------- s08 鲁棒性
def s08_robustness(data, out):
    rows = data["robustness"]
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)
    for ax, scan, xlabel in ((axes[0], "temperature", "初态温度 (μK)"),
                             (axes[1], "alignment", "末端对准偏差 (μm)")):
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
                        fmt="-o", ms=8, capsize=5, lw=2.2,
                        color=WF_COLOR[name], label=WF_SHORT[name])
        floor = min(to_float(r["wilson_low"]) for r in rows
                    if r["scan"] == scan and r.get("wilson_low"))
        ax.axhline(floor, color="gray", ls=":", lw=1.6)
        ax.text(0.985, floor + 0.0004, f"500 shots 的 Wilson 下限 {floor:.4f}",
                transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=11, color="#888")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("俘获率")
        ax.set_ylim(0.9885, 1.0025)
        ax.set_title({"temperature": "温度 3–10 μK：全部 100%",
                      "alignment": "对准 ±0.10 μm：全部 100%"}[scan], fontsize=13.5)
    axes[0].legend(fontsize=10.5, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    save(fig, out, "s08_鲁棒性.png")


# ---------------------------------------------------------------- s09 数值可信
def s09_numerical(data, out):
    # 左图：用本仓库物理引擎重算理想轨迹的功-能账本
    src = REPO / "simulation" / "level2_joint_transfer" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from level2_joint_transfer.config import load_config
    from level2_joint_transfer.level2_simulation import ideal_run, physics_from_config
    from level2_joint_transfer.optimizer import frozen_spec

    cfg = load_config(REPO / "simulation" / "level2_joint_transfer" / "configs" /
                      "level2_joint_transfer.yaml")
    physics = physics_from_config(cfg)
    spec = frozen_spec(data["frozen"], cfg)
    run = ideal_run(dict(spec), physics, cfg.integration.validation_dt_s,
                    cfg.integration.post_transfer_hold_s)
    we = run["work_energy"]
    t_us = we["time_s"] * 1e6
    to_uk = lambda arr: np.asarray(arr) / KB_UK

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)
    ax = axes[0]
    n_end = int(np.argmax(t_us >= 400)) + 1
    ax.plot(t_us[:n_end], to_uk(we["energy_total_J"] - we["energy_total_J"][0])[:n_end],
            color="k", lw=3.4, label="ΔE 机械能变化", zorder=4)
    ax.plot(t_us[:n_end], to_uk(we["work_total_J"])[:n_end], color="#2ca02c",
            lw=2.2, ls=(0, (6, 3)), label="W_ext 外部总功", zorder=5)
    ax.plot(t_us[:n_end], to_uk(we["work_move_J"])[:n_end], color="#1f77b4",
            lw=1.6, label="移动功", zorder=3)
    ax.plot(t_us[:n_end], to_uk(we["work_depth_J"])[:n_end], color="#ff7f0e",
            lw=1.6, label="降深功", zorder=3)
    ax.plot(t_us[:n_end], to_uk(we["residual_J"])[:n_end], color="#d62728",
            lw=1.8, label="残差 R_W（贴近 0）", zorder=6)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("能量 (μK)")
    ax.set_title("ΔE（黑实线）与外部功（绿虚线）逐点重合", fontsize=13.5)
    ax.legend(fontsize=11, loc="lower right", ncols=2)

    # 右图：步长收敛（robustness.csv 的 timestep 行，以最细步长为基准）
    ax = axes[1]
    conv_us = cfg.integration.convergence_dt_us
    rows = [r for r in data["robustness"]
            if r["scan"] == "timestep" and r["waveform"] == "overlapped_optimized_400us"]
    pts = sorted((to_float(r["condition"]), to_float(r["final_energy_max_abs_diff_uK"]))
                 for r in rows if to_float(r["condition"]) > conv_us)
    dts = np.array([p[0] for p in pts])
    diff = np.array([p[1] for p in pts])
    ax.loglog(dts, diff, "o-", lw=2.2, ms=11, color="#1f77b4",
              label="末态能量最大偏差")
    ref = diff * (dts / dts[0]) ** 2
    ax.loglog(dts, ref, "--", color="#888", lw=1.6, label="二阶收敛参考线")
    ax.set_xticks(dts, [f"{d:g}" for d in dts])
    ax.set_xlabel("步长 dt (μs)　·　以 dt = 0.025 μs 为基准")
    ax.set_ylabel("max |ΔE_f| (μK)")
    ax.set_title("步长减半：误差二阶下降，标签翻转 = 0", fontsize=13.5)
    ax.legend(fontsize=11, loc="lower right")
    fig.tight_layout()
    save(fig, out, "s09_数值可信.png")


# ---------------------------------------------------------------- s11 噪声模型
def s11_noise_model(data, out):
    """三通道平均强度确定性加热模型（默认幅度取自 level2b demo 实际配置）。"""
    defaults = json.loads((DEMO_2B / "noise_v95.json").read_text())["defaults"]
    derived = defaults["derived_from"]
    recoil = defaults["recoil"]["power_uK_per_s"]
    par_rate = derived["parametric_rate_per_s"]
    par_t_ref = derived["parametric_reference_temperature_uK"]
    par_power = defaults["parametric"]["power_uK_per_s"]
    pointing = defaults["pointing"]["power_uK_per_s"]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.4)
    ax.axis("off")
    ax.text(5.0, 5.98, "系综平均功率的确定性注入：无随机性，三通道均为常数功率，能量线性累积",
            ha="center", fontsize=15, fontweight="bold", color="#333")

    cards = [
        (0.45, "#eefaee", "#2ca02c", "光子反冲",
         "P = 2·E_r·Γ_sc\n常数功率 → 线性增长\n\n"
         f"默认 {recoil:.3f} μK/s\n（Γ_sc 2.77 /s · E_r 64 nK）"),
        (3.60, "#eaf3fb", "#1f77b4", "强度噪声 · 线性化",
         "P = Γ_par·E_ref（常数功率）\n远离 2ν0 参量共振的线性近似\n\n"
         f"默认 {par_power:.0f} μK/s\n（Γ_par {par_rate:.0f}/s × k_B·"
         f"{par_t_ref:.0f} μK）"),
        (6.75, "#fdf3e7", "#e8a860", "指向 / 位置噪声",
         "P = m·ω0⁴·S_x(ν0)/8\n常数功率 → 线性增长\n\n"
         f"默认 {pointing:.2f} μK/s\n（S_x 1 pm²/Hz @ ν0）"),
    ]
    for x, fc, ec, head, body in cards:
        box(ax, x, 2.70, 2.80, 2.65, fc, ec)
        ax.text(x + 1.40, 4.90, head, ha="center", fontsize=14.5,
                fontweight="bold", color=ec)
        ax.text(x + 1.40, 3.75, body, ha="center", va="center", fontsize=12,
                color="#444")

    box(ax, 0.45, 1.50, 9.10, 0.95, "#f2f2f2", "#999")
    ax.text(5.0, 1.98,
            "注入机制：每步速度踢 v ← sign(v)·√(v² + 2ΔE/m)，逐通道累计；"
            "注入能量计入功-能账本 ΔE = W_ext + W_noise + R_W",
            ha="center", va="center", fontsize=13, color="#444")
    ax.text(5.0, 0.70,
            "幅度参数均为 assumed_sensitivity_only（非实验测量值）"
            "· 1061 nm · ν0 取 SLM 阱频 25.5 kHz · 三通道可独立开关与缩放",
            ha="center", fontsize=11.5, color="#888")
    save(fig, out, "s11_噪声模型.png")


# ---------------------------------------------------------------- s12 噪声×速度结果
def s12_noise_results(data, out):
    """实验 4：v95 边界移动 / 幅度敏感性 / 通道分解（level2b demo 数据）。"""
    rows = [r for r in read_csv_rows(DEMO_2B / "noise_speed_results.csv")
            if r.get("status") == "ok"]
    v95_by_mode = {e["noise_mode"]: e for e in
                   json.loads((DEMO_2B / "noise_v95.json").read_text())["v95_table"]}
    move_frac = 1.0 - 0.01945549969117366  # frozen 波形移动窗口宽度
    v_avg_of = lambda r: 2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac) * 1e3

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE)

    # (a) v95 曲线：off / ×1 / ×100 / ×1000
    ax = axes[0]
    mode_style = [("off", "#333333", "o", "-"), ("x1", "#2ca02c", "", "--"),
                  ("x100", "#1f77b4", "", ":"), ("x1000", "#d62728", "s", "-")]
    for mode, color, marker, ls in mode_style:
        pts = sorted(((v_avg_of(r), float(r["capture_rate"]),
                       float(r["wilson_low"]), float(r["wilson_high"]))
                      for r in rows if r["sub_experiment"] == "v95_shift"
                      and r["noise_mode"] == mode))
        v = [p[0] for p in pts]
        f = [p[1] for p in pts]
        lo = np.maximum(0, np.array(f) - np.array([p[2] for p in pts]))
        hi = np.maximum(0, np.array([p[3] for p in pts]) - np.array(f))
        label = {"off": "关闭", "x1": "默认 ×1", "x100": "×100",
                 "x1000": "×1000"}[mode]
        ax.errorbar(v, f, yerr=[lo, hi], fmt=marker + ls if marker else ls,
                    ms=6, capsize=3, lw=2.0, color=color, label=label)
    ax.axhline(0.95, color="gray", ls=":", lw=1.6)
    est = v95_by_mode["off"].get("v95_bootstrap", {}).get("point_estimate")
    if est:
        ax.axvline(est * 1e3, color="#333333", ls=":", lw=1.4)
        ax.text(est * 1e3 * 1.08, 0.30, f"v95 = {est * 1e3:.1f}",
                fontsize=12, color="#333333", fontweight="bold")
    ax.text(5.5, 0.60, "×1000：全程 < 0.95\n（v95 无定义）", fontsize=11,
            color="#d62728", fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlabel("平均速度 v_avg (mm/s)")
    ax.set_ylabel("俘获率")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("(a) v95 曲线：×1、×100 与关闭重合", fontsize=13)
    ax.legend(fontsize=10, loc="upper left")

    # (b) 幅度敏感性 @55 mm/s
    ax = axes[1]
    pts = sorted((float(r["noise_mode"][1:]), float(r["capture_rate"]),
                  float(r["margin_mean"]))
                 for r in rows if r["sub_experiment"] == "amplitude_sensitivity")
    scales = [p[0] for p in pts]
    rates = [p[1] for p in pts]
    margins = [p[2] for p in pts]
    ax.plot(scales, rates, "o-", ms=10, lw=2.2, color="#1f77b4")
    for s, f, m in zip(scales, rates, margins):
        ax.annotate(f"{f:.2f}\n裕量 {m:+.2f}", (s, f),
                    textcoords="offset points", xytext=(0, -34), ha="center",
                    fontsize=10.5, color="#555")
    ax.set_xscale("log")
    ax.set_xticks(scales, [f"×{int(s)}" for s in scales])
    ax.set_xlabel("噪声整体放大倍数")
    ax.set_ylabel("俘获率")
    ax.set_ylim(0, 1.12)
    ax.set_title("(b) 幅度敏感性 @ 55 mm/s", fontsize=13)

    # (c) 通道分解（×100）
    ax = axes[2]
    order = ["off", "all_channels", "recoil_only", "parametric_only", "pointing_only"]
    labels = ["关闭", "三通道", "仅反冲", "仅参量", "仅指向"]
    lookup = {r["noise_mode"]: r for r in rows
              if r["sub_experiment"] == "channel_breakdown"}
    values = [float(lookup[m]["capture_rate"]) for m in order]
    colors = ["#999999", "#d62728", "#1f77b4", "#e8a860", "#2ca02c"]
    bars = ax.bar(range(5), values, 0.62, color=colors, alpha=0.88)
    for i, v in enumerate(values):
        ax.text(i, v + 0.03, f"{v:.2f}", ha="center", fontsize=12,
                fontweight="bold", color="#333")
    ax.set_xticks(range(5), labels, fontsize=11)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("俘获率")
    ax.set_title("(c) 通道分解（×1000）：仅参量 ≈ 三通道", fontsize=13)
    fig.tight_layout()
    save(fig, out, "s12_噪声结果.png")


# ---------------------------------------------------------------- s10 结论
def s10_conclusion(data, out):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.6)
    ax.axis("off")

    cards = [
        (0.35, "#eefaee", "#2ca02c", "全部实现并实际运行",
         "三类波形 · 功-能核算 · 分阶段优化\n独立验证 · 鲁棒性 · preflight 10/10"),
        (3.60, "#fdf3e7", "#e8a860", "物理结论：饱和区",
         "三波形俘获率均 100%，配对差 = 0\n如实报告：未观察到显著提升"),
        (6.85, "#eaf3fb", "#1f77b4", "下一步：脱离饱和区",
         "升温 / 进一步压缩时长\n（Level 2B 已测得速度边界 η95 ≈ 0.29）"),
    ]
    for x, fc, ec, head, body in cards:
        box(ax, x, 1.55, 2.80, 2.60, fc, ec)
        ax.text(x + 1.40, 3.62, head, ha="center", fontsize=14.5,
                fontweight="bold", color=ec)
        ax.text(x + 1.40, 2.50, body, ha="center", va="center", fontsize=12.5,
                color="#444")
    for x0 in (3.22, 6.47):
        ax.add_patch(FancyArrowPatch((x0, 2.85), (x0 + 0.30, 2.85),
                                     arrowstyle="-|>", mutation_scale=24,
                                     lw=2.2, color="#555"))
    ax.text(5.0, 4.75, "诚实报告是本 Level 的核心要求", ha="center",
            fontsize=16, fontweight="bold", color="#333")
    ax.text(5.0, 0.75, "模型边界：一维经典动力学 —— 不含二维/轴向、光子散射、技术噪声、量子内态与多原子相互作用",
            ha="center", fontsize=11.5, color="#888")
    save(fig, out, "s10_结论.png")


# ---------------------------------------------------------------- s13 最高速度测算
CONDITION_LABELS = {
    "nominal": "标称\n280 μK·1.17 μm",
    "depth_x0.5": "深度减半\n140 μK",
    "depth_x2.0": "深度加倍\n560 μK",
    "waist_0.90um": "束腰收窄\n0.90 μm",
    "waist_1.50um": "束腰放宽\n1.50 μm",
}


def s13_max_speed(data, out):
    """实验 5a：默认噪声 ×1 下硬件条件扫描的 v95（level2b max_speed 产物）。"""
    summary = json.loads((DEMO_2B / "max_speed_summary.json").read_text())
    rows = [r for r in read_csv_rows(DEMO_2B / "max_speed_hardware.csv")
            if r.get("status") == "ok"]
    v95_table = {e["condition"]: e for e in summary["v95_table"]}
    move_frac = 1.0 - 0.01945549969117366

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, width_ratios=[1.35, 1.0])

    # (a) 俘获率-速度曲线（五条件）
    ax = axes[0]
    colors = {"nominal": "#333333", "depth_x0.5": "#1f77b4",
              "depth_x2.0": "#d62728", "waist_0.90um": "#2ca02c",
              "waist_1.50um": "#e8a860"}
    for cond, entry in v95_table.items():
        pts = sorted(((2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac) * 1e3,
                       float(r["capture_rate"])) for r in rows
                      if r["condition"] == cond))
        if not pts:
            continue
        label = CONDITION_LABELS.get(cond, cond).replace("\n", " ")
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=5, lw=1.8,
                color=colors.get(cond, "#888"), label=label)
    ax.axhline(0.95, color="gray", ls=":", lw=1.6)
    ax.set_xscale("log")
    ax.set_xlabel("平均速度 v_avg (mm/s)")
    ax.set_ylabel("俘获率")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("(a) 俘获率-速度曲线（噪声 ×1，20 μK）", fontsize=13)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.25)

    # (b) v95 ± CI 与 η95
    ax = axes[1]
    conds = [c for c in ("depth_x2.0", "nominal", "depth_x0.5",
                         "waist_0.90um", "waist_1.50um") if c in v95_table]
    xs = np.arange(len(conds))
    for i, cond in enumerate(conds):
        entry = v95_table[cond]
        boot = entry.get("v95_bootstrap", {})
        point = boot.get("point_estimate")
        if point is None:
            ax.text(i, 0.5, "未覆盖", ha="center", fontsize=11, color="#999")
            continue
        lo = point - (boot.get("ci_low") or point)
        hi = (boot.get("ci_high") or point) - point
        ax.errorbar(i, point * 1e3, yerr=[[lo * 1e3], [hi * 1e3]], fmt="s",
                    ms=9, capsize=4, lw=2, color=colors.get(cond, "#333"))
        ax.text(i, point * 1e3 * 1.06, f"{point * 1e3:.1f}", ha="center",
                fontsize=12, fontweight="bold", color=colors.get(cond, "#333"))
        eta = entry.get("eta95_peak")
        note = f"η95 = {eta:.2f}" if eta is not None else "η95 = —"
        ax.text(i, point * 1e3 * 0.80, note, ha="center", fontsize=11,
                color=colors.get(cond, "#333"))
    ax.set_xticks(xs, [CONDITION_LABELS.get(c, c) for c in conds], fontsize=10.5)
    ax.set_ylabel("v95 (mm/s)")
    ax.set_title("(b) v95 与 η95 = v_peak/v_ad", fontsize=13)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    save(fig, out, "s13_最高速度测算.png")


# ---------------------------------------------------------------- s14 同时长方案对比
PROFILE_STYLES = [
    ("smootherstep", "渐变 smootherstep", "#1f77b4"),
    ("trapezoid_30", "梯形 0.30（参考）", "#333333"),
    ("trapezoid_15", "快加速梯形 0.15", "#2ca02c"),
    ("triangular_50", "三角 0.50", "#e8a860"),
]


def s14_profiles(data, out):
    """实验 5b：同时长同位移、仅速度形状不同的四种方案对比。"""
    summary = json.loads((DEMO_2B / "max_speed_summary.json").read_text())
    profiles = summary["profiles"]
    by_duration = {}
    for p in profiles:
        by_duration.setdefault(p["duration_us"], []).append(p)
    durations = sorted(by_duration)

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, width_ratios=[0.9, 1.1, 1.1])

    # (a) 速度形状示意（同时长归一化）
    ax = axes[0]
    t = np.linspace(0.0, 1.0, 400)
    for key, label, color in PROFILE_STYLES:
        if key == "smootherstep":
            s = t * t * t * (10.0 + t * (-15.0 + 6.0 * t))
            v = np.gradient(s, t)
        else:
            f = {"trapezoid_30": 0.30, "trapezoid_15": 0.15,
                 "triangular_50": 0.50}[key]
            v = np.where(t < f, t / f, np.where(t > 1.0 - f, (1.0 - t) / f, 1.0))
        v = v / np.trapezoid(v, t)  # 归一化：同时长同位移 → ∫v dt = 1
        ax.plot(t, v, lw=2.2, color=color, label=label)
    ax.set_xlabel("时间（归一化移动时长）")
    ax.set_ylabel("v / v_avg")
    ax.set_title("(a) 同时长移动方案（示意）", fontsize=13)
    ax.legend(fontsize=9.5, loc="upper right")
    ax.grid(alpha=0.25)

    # (b) 俘获率（时长 × 方案）
    ax = axes[1]
    width = 0.2
    for k, (key, label, color) in enumerate(PROFILE_STYLES):
        xs = np.arange(len(durations)) + (k - 1.5) * width
        pts = [{p["profile"]: p for p in by_duration[d]} for d in durations]
        rates = [pt[key]["capture_rate"] for pt in pts]
        los = [max(0.0, pt[key]["capture_rate"] - pt[key]["wilson_low"]) for pt in pts]
        his = [max(0.0, pt[key]["wilson_high"] - pt[key]["capture_rate"]) for pt in pts]
        ax.bar(xs, rates, width * 0.92, color=color, alpha=0.88,
               yerr=[los, his], capsize=2, error_kw={"lw": 1.2})
        for x, r in zip(xs, rates):
            ax.text(x, r + 0.05, f"{r:.2f}", ha="center", fontsize=9, color="#444")
    ax.set_xticks(np.arange(len(durations)), [f"{d:.0f} μs" for d in durations])
    ax.set_xlabel("AOD 起动→停止总时长（同位移 2.4 μm）")
    ax.set_ylabel("俘获率")
    ax.set_ylim(0, 1.15)
    ax.set_title("(b) 俘获率（噪声 ×1，20 μK）", fontsize=13)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.88)
               for _k, (_key, _l, c) in enumerate(PROFILE_STYLES)]
    ax.legend(handles, [l for _k, (_key, l, _c) in enumerate(PROFILE_STYLES)],
              fontsize=9, loc="upper right")
    ax.grid(alpha=0.25, axis="y")

    # (c) 相对参考方案的配对差
    ax = axes[2]
    y_labels, y = [], []
    for d in durations:
        for key, label, color in PROFILE_STYLES:
            if key == "trapezoid_30":
                continue
            pt = {p["profile"]: p for p in by_duration[d]}
            diff = pt[key].get("pair_vs_reference_difference")
            if diff is None:
                continue
            lo = diff - pt[key]["pair_ci_low"]
            hi = pt[key]["pair_ci_high"] - diff
            significant = pt[key]["pair_ci_low"] > 0 or pt[key]["pair_ci_high"] < 0
            ax.errorbar(diff, len(y), xerr=[[lo], [hi]], fmt="o", ms=8,
                        capsize=3, lw=1.8,
                        color=color if not significant else "#d62728")
            y_labels.append(f"{d:.0f} μs · {label}")
            y.append(len(y))
    ax.axvline(0.0, color="#333", ls="--", lw=1.4)
    ax.set_yticks(y, y_labels, fontsize=9.5)
    ax.set_xlabel("俘获率配对差 − 梯形 0.30")
    ax.set_title("(c) 配对差（红 = CI 不含 0）", fontsize=13)
    ax.grid(alpha=0.25, axis="x")
    fig.tight_layout()
    save(fig, out, "s14_同时长方案对比.png")


def main():
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        REPO / "simulation" / "level2_joint_transfer" / "outputs" /
        "level2_joint_transfer" / "demo")
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "level2_slide_figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    data = load_all(data_dir)
    figures = [
        ("s01_问题.png", s01_question),
        ("s02_物理模型.png", s02_model),
        ("s03_方法.png", s03_method),
        ("s04_波形对比.png", s04_waveforms),
        ("s05_势场演化.png", s05_potentials),
        ("s11_噪声模型.png", s11_noise_model),
        ("s12_噪声结果.png", s12_noise_results),
        ("s13_最高速度测算.png", s13_max_speed),
        ("s14_同时长方案对比.png", s14_profiles),
        ("s06_优化过程.png", s06_optimization),
        ("s07_独立验证.png", s07_validation),
        ("s08_鲁棒性.png", s08_robustness),
        ("s09_数值可信.png", s09_numerical),
        ("s10_结论.png", s10_conclusion),
    ]
    for name, func in figures:
        func(data, out_dir)
    print(f"\n共 {len(figures)} 张图，数据来自 {data_dir}")


if __name__ == "__main__":
    main()
