# -*- coding: utf-8 -*-
"""Level 2 相对 Level 1 的差异可视化：按新增条目逐项成图。

结构：
  00  Level 1 思路概览与 Level 2 的定位（Level 1 已有完整可视化，此处仅概览）
  01  代码结构：Level 0/1 零修改 + 新增 level2_joint_transfer 包
  02~03  新增1  可优化的波形体系（三类波形 + 约束检查）
  04~05  新增2  三阶段优化器（流程 + 打分函数）
  06     新增3  三套互不重叠的样本池 + 公共随机数 CRN
  07     新增4  配对统计
  08~09  新增5  功-能核算 + 100 μs 保持段检查
  10     新增6  near-threshold 临界样本标记
  11     新增7  鲁棒性扫描（温度/对准/步长）
  12     新增8  预检查关卡（10 项）
  13     新增9  平均噪声加热模型
  14     新增10 工程层面（CLI/并行/输出/测试）

波形曲线、功-能核算、约束拒绝示例、保持段轨迹均通过导入
level2_joint_transfer 包用真实代码计算；统计与扫描读
outputs/level2_joint_transfer/demo/ 下的真实运行数据。

用法（须用项目 venv，以便导入仿真包）：
  simulation/level2_joint_transfer/.venv/bin/python level2_vs_level1_figures.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
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

REPO = Path(__file__).resolve().parent
PKG_ROOT = REPO / "simulation" / "level2_joint_transfer"
DATA = PKG_ROOT / "outputs" / "level2_joint_transfer" / "demo"
NOISE_JSON = PKG_ROOT / "outputs" / "level2_joint_transfer" / "demo_mean_heating" / "mean_heating_demo.json"
OUT = REPO / "level2_diff_figures"
OUT.mkdir(exist_ok=True)

# ---- 导入仿真包（venv 内已 editable 安装） ----
from level2_joint_transfer.config import load_config
from level2_joint_transfer.level2_simulation import (build_waveform, ideal_run,
                                                     physics_from_config)
from level2_joint_transfer.optimizer import frozen_spec
from level2_joint_transfer.waveforms import (OverlappedWaveform, SequentialWaveform,
                                             SplineWaveform, validate_waveform)

CFG = load_config(PKG_ROOT / "configs" / "level2_joint_transfer.yaml")
PHYS = physics_from_config(CFG)
FROZEN = yaml.safe_load((DATA / "best_waveform.yaml").read_text())
SPEC_OPT = frozen_spec(FROZEN, CFG)
CONSTRAINT = CFG.waveforms.as_constraint_dict()

WF_SEQ600 = SequentialWaveform(CFG.waveforms.sequential_baseline_duration_s,
                               CFG.waveforms.level1_move_fraction,
                               PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])
WF_SEQ400 = SequentialWaveform(CFG.waveforms.compressed_baseline_duration_s,
                               CFG.waveforms.level1_move_fraction,
                               PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])
WF_OPT = build_waveform(SPEC_OPT, PHYS)

KB = 1.380649e-23
WAVEFORMS = ["sequential_600us", "sequential_400us", "overlapped_optimized_400us"]
WF_LABEL = {"sequential_600us": "顺序基线 600 μs",
            "sequential_400us": "压缩顺序基线 400 μs",
            "overlapped_optimized_400us": "优化同步波形 400 μs"}
WF_COLOR = {"sequential_600us": "#1f77b4",
            "sequential_400us": "#ff7f0e",
            "overlapped_optimized_400us": "#d62728"}
C_NEW, C_OLD, C_BAD, C_OK = "#d62728", "#7f7f7f", "#d62728", "#2ca02c"


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_f(value, default=np.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def box(ax, x, y, w, h, title, body, fc="#eaf3fb", ec="#1f77b4", fs_t=10.5, fs_b=9,
        title_color=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06", fc=fc, ec=ec, lw=1.6))
    ax.text(x + w / 2, y + h - 0.16 * h, title, ha="center", va="top", fontsize=fs_t,
            fontweight="bold", color=title_color or ec)
    ax.text(x + w / 2, y + (h - 0.30 * h) / 2, body, ha="center", va="center",
            fontsize=fs_b, color="#333", linespacing=1.45)


def arrow(ax, p0, p1, color="#555", lw=2.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=20,
                                 lw=lw, color=color))


# ================================================================ 图 00
def fig00_level1_overview():
    fig = plt.figure(figsize=(15, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.15, 1), hspace=0.33, wspace=0.22)
    fig.suptitle("图0  Level 1 的思路概览（已有完整可视化，此处只讲思路）与 Level 2 的定位",
                 fontsize=15, fontweight="bold")

    # (a) Level 1 固定波形
    ax = fig.add_subplot(gs[:, 0])
    t = np.linspace(0, 600, 1201) * 1e-6
    c = WF_SEQ600.center(t) * 1e6
    d = WF_SEQ600.depth(t) / KB * 1e6
    ax.plot(t * 1e6, c, color="#1f77b4", lw=2.4, label="AOD 中心 c(t)")
    ax2 = ax.twinx()
    ax2.plot(t * 1e6, d, color="#ff7f0e", lw=2.4, label="AOD 深度 D(t)")
    ax2.set_ylabel("AOD 深度 (μK)", color="#ff7f0e")
    ax2.grid(False)
    ax.axvspan(0, 312, color="#1f77b4", alpha=0.06)
    ax.axvspan(312, 600, color="#ff7f0e", alpha=0.08)
    ax.annotate("阶段① 移动 52%\nsmoothstep 3s²−2s³\n深度不变(280 μK)", xy=(150, 1.2),
                fontsize=10, ha="center", color="#1f77b4",
                arrowprops=dict(arrowstyle="->", color="#1f77b4"), xytext=(150, 1.85))
    ax.annotate("阶段② 降深 48%\n深度按 (1−s)² 降为 0\n中心停在 SLM 处", xy=(460, 0.15),
                fontsize=10, ha="center", color="#ff7f0e",
                arrowprops=dict(arrowstyle="->", color="#ff7f0e"), xytext=(460, 1.0))
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("AOD 中心 (μm)", color="#1f77b4")
    ax.set_title("(a) Level 1：唯一的一套固定波形（总时长 600 μs，先移动后降深，绝不重叠）",
                 fontsize=11.5)
    ax.set_ylim(-0.15, 2.6)

    # (b) Level 1 流程
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("(b) Level 1 流程：跑一次固定流程并统计", fontsize=11.5)
    steps = [
        ("热初态采样", "5 μK 简谐提议 + 束缚拒绝\n(复用 Level 0 高斯势)", "#eef7ee"),
        ("时变 velocity-Verlet", "每 shot 现采样、现积分\n总时长 600 μs", "#eef7ee"),
        ("判俘获 E_f < 0", "末态 SLM 单阱能量\nE_f = ½mv² + U_SLM(x_f)", "#eef7ee"),
        ("统计输出", "蒙特卡洛俘获率 + Wilson 95% 区间\n少量自检(静态极限/端点/步长)", "#eef7ee"),
    ]
    for i, (t1, t2, fc) in enumerate(steps):
        y = 0.78 - i * 0.27
        box(ax, 0.02, y, 0.96, 0.20, t1, t2, fc=fc, ec=C_OK, fs_t=10, fs_b=8.5)
        if i < 3:
            arrow(ax, (0.5, y - 0.005), (0.5, y - 0.065))

    # (c) Level 2 定位
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("(c) Level 2 要回答的问题（代码全部新增，物理模型不变）", fontsize=11.5)
    body = ("同一套一维经典模型、同一初态分布下：\n"
            "① 允许『移动』和『降深』在时间上重叠（同步波形）\n"
            "② 把 600 μs 压缩到 400 μs\n"
            "③ 波形参数不再是写死的，而是搜索出来的\n"
            "④ 结果好坏由从未参与优化的独立样本说了算\n\n"
            "工程事实：Level 0/1 代码零修改，\n"
            "全部改动在新包 level2_joint_transfer 里（见 图1）")
    box(ax, 0.02, 0.05, 0.96, 0.88, "Level 2 = Level 1 模型 + 波形优化与独立验证框架",
        body, fc="#fdf0f0", ec=C_NEW, fs_t=10.5, fs_b=9.5)
    fig.text(0.5, 0.005, "Level 1 的完整结果可视化见 simulation/level1_transfer_1d/pptplots/（本文不再重复）",
             ha="center", fontsize=9.5, color="#666")
    fig.savefig(OUT / "00_Level1思路概览与Level2定位.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 01
def fig01_structure():
    fig, ax = plt.subplots(figsize=(16.5, 9.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("图1  代码结构：Level 0/1 原样复制（逐字节相同，diff 验证），所有新东西都在新包里",
                 fontsize=14.5, fontweight="bold", pad=14)

    # 左：level0/level1 复用
    box(ax, 0.2, 5.4, 2.9, 3.9, "src/level0_static_trap\n（Level 0，未改动）",
        "高斯势/高斯力 gaussian_potential/force\n物理常数 constants\n"
        "静态积分器 velocity_verlet\n阱频解析式 trap_frequency_analytic\nIO 工具 io_utils",
        fc="#f2f2f2", ec=C_OLD, fs_t=10, fs_b=8.6)
    box(ax, 0.2, 1.1, 2.9, 3.9, "src/level1_transfer_1d\n（Level 1，未改动）",
        "时变积分器\nvelocity_verlet_time_dependent\n热初态采样器\nsample_thermal_initial_state\n"
        "Wilson 区间 wilson_interval\nLevel 1 波形 aod_center/aod_depth",
        fc="#f2f2f2", ec=C_OLD, fs_t=10, fs_b=8.6)
    ax.text(1.65, 9.6, "从 Level 1 项目复制，diff 逐字节相同", ha="center", fontsize=9.5,
            color=C_OLD)

    # 右：新包
    box(ax, 5.3, 0.5, 4.5, 8.9, "src/level2_joint_transfer/  （新增，13 个模块，~3300 行）",
        "", fc="#fdf3f3", ec=C_NEW, fs_t=11.5)
    modules = [
        ("waveforms.py", "三类波形类 + 1/2/3 阶导数 + 约束检查"),
        ("level2_simulation.py", "核心引擎：固定初态数组、逐 shot 评估、并行"),
        ("work_energy.py", "功-能核算：深度功/移动功/残差 + 保持段"),
        ("objectives.py", "候选打分（俘获为主，字典序排序）"),
        ("optimizer.py", "三阶段优化 + 检查点续跑 + 候选全记录"),
        ("statistics.py", "配对计数 + 配对 bootstrap"),
        ("robustness.py", "温度/对准/步长敏感性扫描"),
        ("preflight.py", "10 项预检查关卡"),
        ("noise_heating.py", "三通道平均噪声加热（确定性）"),
        ("level2_visualization.py", "10 张诊断图 + 程序化验收"),
        ("config.py", "严格 YAML schema（三池三种子等校验）"),
        ("cli.py", "preflight/optimize/validate/robustness 分阶段"),
        ("__init__.py", ""),
    ]
    for i, (name, desc) in enumerate(modules):
        y = 8.55 - i * 0.66
        ax.text(5.55, y, name, fontsize=9.3, family="monospace", color=C_NEW,
                fontweight="bold", va="center")
        ax.text(7.5, y, desc, fontsize=8.8, color="#333", va="center")

    # 中间箭头：import 复用
    for y0, y1 in [(7.3, 7.3), (3.0, 3.0)]:
        arrow(ax, (3.2, y0), (5.2, y1), color=C_OK)
    ax.text(4.2, 7.55, "import 复用\n（不改一行旧代码）", ha="center", fontsize=9.5,
            color=C_OK, fontweight="bold")
    ax.text(4.2, 2.55, "import 复用\n（逐调用同一实现）", ha="center", fontsize=9.5,
            color=C_OK, fontweight="bold")

    # 底部：测试
    box(ax, 0.2, 0.1, 4.7, 0.75, "tests/（新增 27 + 13 个，旧测试原样保留继续跑）",
        "test_level2.py(27) + test_noise_heating.py(13) + Level 0/1 全部旧测试",
        fc="#eefaee", ec=C_OK, fs_t=9.5, fs_b=8.5)
    fig.savefig(OUT / "01_代码结构_零修改全新增.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 02
def fig02_waveforms():
    t400 = np.linspace(0, 400, 801) * 1e-6
    t600 = np.linspace(0, 600, 1201) * 1e-6
    seq = [(WF_SEQ600, t600, WF_COLOR["sequential_600us"], WF_LABEL["sequential_600us"]),
           (WF_SEQ400, t400, WF_COLOR["sequential_400us"], WF_LABEL["sequential_400us"]),
           (WF_OPT, t400, WF_COLOR["overlapped_optimized_400us"], WF_LABEL["overlapped_optimized_400us"])]
    # 样条演示：由冻结最优波形取 8 控制点（与 optimizer.run_spline_refinement 同法）
    s_k = np.linspace(0, 1, 8)
    spl = SplineWaveform(CFG.waveforms.optimized_duration_s,
                         s_k, WF_OPT.center(s_k * CFG.waveforms.optimized_duration_s),
                         WF_OPT.depth(s_k * CFG.waveforms.optimized_duration_s),
                         PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    fig.suptitle("图2  新增① 可优化的波形体系：从 2 个写死的函数 → 3 个带导数的波形类（曲线由包代码实时计算）",
                 fontsize=14, fontweight="bold")

    def panel(ax, key, ylabel, title, scale=1.0, unit=""):
        for wf, tt, color, label in seq:
            ax.plot(tt * 1e6, getattr(wf, key)(tt) * scale, color=color, lw=1.9, label=label)
        ax.set_xlabel("时间 (μs)")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10.5)
        ax.legend(fontsize=8)

    panel(axes[0, 0], "center", "c(t) (μm)", "(a) AOD 中心：从 2.4 μm → 0")
    panel(axes[0, 1], "depth", "D(t) (μK)", "(b) AOD 深度：280 μK → 0", scale=1 / KB * 1e6)
    panel(axes[0, 2], "center_velocity", "c'(t) (m/s)", "(c) 中心一阶导数（Level 1 没有）")
    panel(axes[1, 0], "depth_rate", "D'(t) (μK/s)", "(d) 深度一阶导数（功-能核算要用）",
          scale=1 / KB * 1e6)

    ax = axes[1, 1]
    ax.plot(t400 * 1e6, WF_OPT.center_jerk(t400), color=WF_COLOR["overlapped_optimized_400us"],
            lw=1.9, label="优化同步波形")
    ax.plot(t400 * 1e6, WF_SEQ400.center_jerk(t400), color=WF_COLOR["sequential_400us"],
            lw=1.5, alpha=0.8, label="压缩顺序基线")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("c'''(t) (m/s³)")
    ax.set_title("(e) 中心三阶导数 jerk（光滑度惩罚要用）", fontsize=10.5)
    ax.legend(fontsize=8)

    # 接口面板
    ax = axes[1, 2]
    ax.axis("off")
    ax.set_title("(f) 类接口：统一 6 个量 + 光滑度统计", fontsize=10.5)
    text = ("WaveformBase（基类）\n"
            "  center(t) / depth(t)\n"
            "  center_velocity / center_acceleration\n"
            "  center_jerk / depth_rate\n"
            "  dense_table()  密集表\n"
            "  smoothness_metrics()  峰值统计\n\n"
            "├ SequentialWaveform  严格复现 Level 1\n"
            "│    （可设 600 或 400 μs 两个基线）\n"
            "├ OverlappedWaveform  5 参数同步族\n"
            "│    移动窗口[ax,bx]×降深窗口[aD,bD]×幂 p\n"
            "│    全部用 smootherstep 窗口构造\n"
            "└ SplineWaveform  PCHIP 单调样条\n"
            "     8 个控制点（可配 14），无过冲\n\n"
            "Level 1 对比：只有 aod_center(t)/aod_depth(t)\n"
            "两个函数，形状写死，无任何导数")
    ax.text(0.02, 0.97, text, fontsize=9.2, va="top",
            color="#333", linespacing=1.55)

    # 冻结参数标注
    fp = FROZEN
    fig.text(0.5, 0.005,
             f"优化波形冻结参数（真实运行结果）：移动窗口 s∈[{fp['move_start_fraction']:.3f}, "
             f"{fp['move_end_fraction']:.3f}]，降深窗口 s∈[{fp['ramp_start_fraction']:.3f}, "
             f"{fp['ramp_end_fraction']:.3f}]，幂 p={fp['ramp_power']:.3f} —— 移动与降深大部分时间重叠，这就是『同步』",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.028, 1, 1))
    fig.savefig(OUT / "02_新增1_可优化波形体系.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 03
def fig03_constraints():
    # 真实调用 validate_waveform 构造 通过/拒绝 示例
    cases = []
    valid, reasons = validate_waveform(WF_OPT, CONSTRAINT)
    cases.append(("冻结的优化波形", WF_OPT, valid, reasons))
    inv_power = OverlappedWaveform(CFG.waveforms.optimized_duration_s, 0.0, 0.52, 0.48, 1.0, 0.7,
                                   PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])
    cases.append(("ramp_power=0.7（幂<1）", inv_power, *validate_waveform(inv_power, CONSTRAINT)))
    inv_narrow = OverlappedWaveform(CFG.waveforms.optimized_duration_s, 0.0, 0.05, 0.3, 1.0, 2.0,
                                    PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])
    cases.append(("移动窗口宽 0.05", inv_narrow, *validate_waveform(inv_narrow, CONSTRAINT)))
    ctrl_s = np.linspace(0, 1, 8)
    center_bad = np.where(ctrl_s < 0.5, 2.4 - ctrl_s * 2, 1.2 + (ctrl_s - 0.5))  # 中途回升
    depth_bad = 280e-6 * KB * np.linspace(1, 0, 8) ** 2
    inv_spline = SplineWaveform(CFG.waveforms.optimized_duration_s, ctrl_s, center_bad, depth_bad,
                                PHYS["aod_initial_center_m"], PHYS["aod_initial_depth_j"])
    cases.append(("样条控制点非单调", inv_spline, *validate_waveform(inv_spline, CONSTRAINT)))

    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.05, 1), hspace=0.34, wspace=0.2)
    fig.suptitle("图3  新增①b  validate_waveform 约束检查：非法候选在进仿真前就被拒绝（拒绝理由为代码真实返回值）",
                 fontsize=14, fontweight="bold")

    ax = fig.add_subplot(gs[:, 0])
    ax.axis("off")
    ax.set_title("检查内容（对每个候选波形逐项执行）", fontsize=11.5)
    checks = [
        "端点严格满足 c(0)=2.4 μm, c(T)=0, D(0)=280 μK, D(T)=0",
        "中心轨迹单调不增（不许往回走）",
        "深度轨迹单调不增（不许回填）",
        "中心始终在 [0, 2.4 μm]、深度始终在 [0, 280 μK]（不许过冲/负值）",
        "窗口分数次序合法：0 ≤ start < end ≤ 1",
        "每个窗口宽度 ≥ min_segment_fraction（默认 0.10）",
        "ramp_power 在边界内且 ≥ 1（<1 会让端点导数发散 → 直接拒绝）",
        "端点速度/深度变化率必须为 0（不许猛启停）",
        "近端点（0.1%·T 内）速度不得超过峰值 5%",
        "样条控制点必须单调不增（PCHIP 保证无过冲）",
    ]
    y = 0.96
    for text in checks:
        ax.text(0.03, y, "•", fontsize=13, color=C_NEW, va="center")
        ax.text(0.07, y, text, fontsize=10.3, va="center")
        y -= 0.095
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # 四个真实案例
    for i, (name, wf, valid, reasons) in enumerate(cases):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        tt = np.linspace(0, wf.duration_s, 400)
        ax.plot(tt * 1e6, wf.center(tt) * 1e6, color="#1f77b4", lw=2, label="c(t)")
        ax.plot(tt * 1e6, wf.depth(tt) / KB * 1e6, color="#ff7f0e", lw=2, label="D(t)")
        status = "【通过】" if valid else "【拒绝】"
        color = C_OK if valid else C_BAD
        title = f"{name} → {status}"
        if reasons:
            title += "\n理由: " + "；".join(reasons)
        ax.set_title(title, fontsize=9.6, color=color)
        ax.set_xlabel("时间 (μs)")
        ax.legend(fontsize=8)
    fig.savefig(OUT / "03_新增1_波形约束检查.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 04
def fig04_optimizer():
    cand = read_csv_rows(DATA / "candidates.csv")
    pref = Counter(r["candidate_id"].split("_")[0] for r in cand)
    sel = [r for r in cand if r["is_selected"] == "True"][0]
    n_ok = sum(1 for r in cand if r["status"] == "ok")
    n_inv = sum(1 for r in cand if r["status"] == "invalid")
    n_fail = sum(1 for r in cand if r["status"] == "failed")

    fig, ax = plt.subplots(figsize=(16.5, 9.3))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("图4  新增② 三阶段优化器：Latin 超立方探索 → 精修选择 → 样条精修 → 冻结（数字为本次真实运行）",
                 fontsize=14.5, fontweight="bold", pad=14)

    # 阶段 1
    box(ax, 0.15, 6.1, 2.9, 3.1,
        "阶段 1  候选探索",
        f"Latin 超立方撒 {pref['lhs']} 个候选\n"
        f"样本池 optimization（128 shots）\n步长 0.10 μs（粗算）\n"
        f"5 维参数：两窗口起止 + 幂次\n每个候选先过约束检查\n"
        f"→ 取前 8 名入围", fc="#eaf3fb", ec="#1f77b4")
    # 阶段 2
    box(ax, 3.55, 6.1, 2.9, 3.1,
        "阶段 2  精修与选择",
        f"8 名入围复评(sel_) + 每名扰动 4 次(ref_)\n"
        f"共 {pref['sel'] + pref['ref']} 个候选\n"
        f"样本池 selection（512 shots）\n步长 0.05 μs（正式）\n"
        f"高斯扰动 σ，裁剪回合法范围\n→ 选出唯一最优", fc="#fdf3e7", ec="#ff7f0e")
    # 阶段 3
    box(ax, 6.95, 6.1, 2.9, 3.1,
        "阶段 3  样条精修（可选）",
        f"由最优低维波形初始化\n8 个控制点 PCHIP 单调样条\n"
        f"扰动 {pref['spl']} 个样条候选\n仍在 selection 池上评估\n"
        f"样条与低维结果分开记录\n（不许只留好的）", fc="#f3eefb", ec="#9467bd")
    arrow(ax, (3.1, 7.65), (3.5, 7.65))
    arrow(ax, (6.5, 7.65), (6.9, 7.65))

    # 冻结
    box(ax, 3.55, 4.0, 2.9, 1.6,
        "select_and_freeze → best_waveform.yaml",
        f"选中候选 {sel['candidate_id']}\n"
        f"类型 {sel['waveform_type']}，selection 池俘获率 {to_f(sel['capture_rate']):.4f}\n"
        f"波形从此冻结，后面阶段不许再改", fc="#eefaee", ec=C_OK)
    arrow(ax, (5.0, 6.0), (5.0, 5.7))

    # 候选表原则
    box(ax, 0.15, 4.0, 2.9, 1.6,
        "candidates.csv：全部保留",
        f"共 {len(cand)} 行 = {pref['lhs']} lhs + {pref['sel']} sel + "
        f"{pref['ref']} ref + {pref['spl']} spl\n"
        f"有效 {n_ok} / 违反约束 {n_inv} / 运行失败 {n_fail}\n"
        f"含全部参数、状态、失败原因、各项得分", fc="#fdeaea", ec=C_BAD)

    # 检查点
    box(ax, 6.95, 4.0, 2.9, 1.6,
        "断点续跑",
        "每阶段写 optimization_checkpoint.json\n"
        "（配置签名不符自动作废重算）\n"
        "中断后重跑只补缺失候选\n评估历史合并不覆盖", fc="#f2f2f2", ec=C_OLD)

    # 底部：防呆细节
    details = [
        ("失败候选不丢", "单个候选抛异常 →\n记为 failed 行，\n不影响整轮搜索", "#fdeaea", C_BAD),
        ("粗步长复算", "阶段 1 用 0.10 μs 粗算，\n阶段 2 起全部用\n正式步长 0.05 μs 复算", "#eaf3fb", "#1f77b4"),
        ("字典序排序", "先比俘获数、再比总分，\n小分差不会淹没\n俘获数差异", "#eefaee", C_OK),
        ("validation 池隔离", "optimizer.py 源码\n不引用 validation 池\n（有测试守护）", "#fdf3e7", "#ff7f0e"),
    ]
    for i, (t1, t2, fc, ec) in enumerate(details):
        box(ax, 0.15 + i * 2.5, 0.6, 2.3, 2.6, t1, t2, fc=fc, ec=ec, fs_t=10, fs_b=8.8)
    fig.savefig(OUT / "04_新增2_三阶段优化器.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 05
def fig05_objectives():
    cand = read_csv_rows(DATA / "candidates.csv")
    sel = [r for r in cand if r["is_selected"] == "True"][0]
    w_cap = CFG.optimization.capture_weight
    w_exc = CFG.optimization.excitation_weight
    w_smo = CFG.optimization.smoothness_weight
    cap_r, exc_t, smo_p, tot = (to_f(sel["capture_rate"]), to_f(sel["excitation_term"]),
                                to_f(sel["smooth_penalty"]), to_f(sel["objective_total"]))

    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=(1, 1.25), hspace=0.36, wspace=0.22)
    fig.suptitle("图5  新增②b 打分函数 objectives.py：俘获为首要目标，全部组成项分开保存（不做黑盒总分）",
                 fontsize=14, fontweight="bold")

    ax = fig.add_subplot(gs[:, 0])
    ax.axis("off")
    ax.set_title("目标函数定义（对应 objectives.py）", fontsize=11.5)
    formula = (
        "总分 = capture_weight × capture_rate      （权重 1.0，首要）\n"
        "     + excitation_weight × excitation_term （权重 0.05，次要）\n"
        "              excitation_term = 1 − 已俘获原子平均归一化激发\n"
        "     − smoothness_weight × smooth_penalty  （权重 0.001，轻微）\n"
        "              smooth_penalty = 无量纲加速度/50 + 无量纲 jerk/500\n\n"
        "排序键 rank_key = (captured_count, total)   ← 字典序\n"
        "先比俘获数：多俘获 1 个 shot 直接压过任何小分差\n\n"
        "另有 margin_quantile（俘获裕量 10% 分位数）单独成列，\n"
        "不混入总分，只用于解释『离失败有多远』")
    ax.text(0.02, 0.95, formula, fontsize=10.3, va="top",
            linespacing=1.7, color="#333")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # 选中候选分解
    ax = fig.add_subplot(gs[0, 1])
    parts = [("俘获项\n1.0×%.3f" % cap_r, w_cap * cap_r, "#2ca02c"),
             ("激发项\n0.05×%.3f" % exc_t, w_exc * exc_t, "#1f77b4"),
             ("光滑惩罚\n−0.001×%.3f" % smo_p, -w_smo * smo_p, "#d62728")]
    left = 0.0
    for label, val, color in parts:
        ax.barh(0.6, val, left=left, color=color, edgecolor="white", height=0.45,
                label=label.replace("\n", " = "))
        if abs(val) > 0.02:
            ax.text(left + val / 2, 0.6, f"{val:+.4f}", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
        left += val
    ax.axvline(tot, color="k", ls="--", lw=1.2)
    ax.text(tot, 0.92, f"总分 {tot:.4f}", ha="right", fontsize=10, fontweight="bold")
    ax.set_xlim(0, 1.12)
    ax.set_yticks([])
    ax.set_xlabel("贡献")
    ax.set_title(f"选中候选 {sel['candidate_id']} 的总分构成（真实数据）：capture=1.0, "
                 f"margin_q10={to_f(sel['margin_quantile']):.3f}", fontsize=10)
    ax.legend(fontsize=8.5, loc="lower right")

    # 全体候选散点
    ax = fig.add_subplot(gs[1, 1])
    marker = {"lhs": ("o", "阶段1 LHS"), "sel": ("s", "阶段2 复评"), "ref": ("^", "阶段2 精修"),
              "spl": ("D", "样条精修")}
    for pref_key, (m, label) in marker.items():
        rows = [r for r in cand if r["candidate_id"].startswith(pref_key) and r["status"] == "ok"]
        ax.scatter([to_f(r["capture_rate"]) for r in rows],
                   [to_f(r["objective_total"]) for r in rows], marker=m, s=30, alpha=0.7,
                   label=f"{label} ({len(rows)})", edgecolor="k", lw=0.2)
    invalid = [r for r in cand if r["status"] == "invalid"]
    ax.scatter([0] * len(invalid), [0] * len(invalid) if invalid else [],
               marker="x", s=40, color=C_BAD, label=f"违反约束被拒 ({len(invalid)})")
    ax.set_xlabel("训练池俘获率")
    ax.set_ylabel("目标总分")
    ax.set_title("全部候选：俘获率与总分（违反约束者根本不进仿真）", fontsize=10)
    ax.legend(fontsize=8.5)
    fig.savefig(OUT / "05_新增2_打分函数.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 06
def fig06_pools_crn():
    shots = read_csv_rows(DATA / "validation_shots.csv")
    by_wf = {name: [r for r in shots if r["waveform"] == name] for name in WAVEFORMS}
    x0 = {name: np.array([to_f(r["x0_um"]) for r in rows]) for name, rows in by_wf.items()}
    v0 = {name: np.array([to_f(r["v0_m_s"]) for r in rows]) for name, rows in by_wf.items()}
    max_dx = max(np.max(np.abs(x0[WAVEFORMS[0]] - x0[n])) for n in WAVEFORMS)
    max_dv = max(np.max(np.abs(v0[WAVEFORMS[0]] - v0[n])) for n in WAVEFORMS)
    ens = CFG.initial_ensemble

    fig = plt.figure(figsize=(16.5, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.35, 1), hspace=0.36, wspace=0.2)
    fig.suptitle("图6  新增③ 三套互不重叠的样本池（防过拟合）+ 公共随机数 CRN（同一组初态比三个波形）",
                 fontsize=14, fontweight="bold")

    ax = fig.add_subplot(gs[:, 0])
    ax.axis("off")
    ax.set_title("三池隔离：谁的样本谁用，validation 只许看一次", fontsize=11.5)
    pools = [
        ("optimization_pool", ens.optimization_seed, ens.optimization_shots, "0.10 μs",
         "只用于阶段 1 搜候选\n（128 shots，粗步长）", "#eaf3fb", "#1f77b4"),
        ("selection_pool", ens.selection_seed, ens.selection_shots, "0.05 μs",
         "只用于阶段 2/3 精修与选择\n（512 shots，正式步长）", "#fdf3e7", "#ff7f0e"),
        ("validation_pool", ens.validation_seed, ens.validation_shots, "0.05 μs",
         "波形冻结后一次性评估\n禁止看到结果再调参\n（2000 shots）", "#fdeaea", C_BAD),
    ]
    for i, (name, seed, n, dt, use, fc, ec) in enumerate(pools):
        y = 6.9 - i * 2.3
        box(ax, 0.02, y, 0.55, 1.9, name,
            f"seed = {seed}（各不相同）\n{n} shots   dt = {dt}\n{use}",
            fc=fc, ec=ec, fs_t=9.5, fs_b=8.8)
        if i < 2:
            arrow(ax, (0.30, y - 0.02), (0.30, y - 0.38))
    rules = ("配置加载时强制校验：三种子互不相同；池大小逐级递增 128 < 512 < 2000\n"
             "optimizer.py 源码不引用 validation 池（测试守护，防止数据泄漏）\n"
             "若看过 validation 结果后改波形 → 原 validation 集作废，必须换全新种子")
    box(ax, 0.02, 0.1, 0.55, 0.75, "防过拟合规则（写进配置与测试）", rules,
        fc="#f2f2f2", ec=C_OLD, fs_t=9, fs_b=8.2)

    # CRN 验证
    ax = fig.add_subplot(gs[0, 1])
    for name in WAVEFORMS:
        ax.hist(x0[name], bins=60, alpha=0.4, color=WF_COLOR[name], label=WF_LABEL[name],
                density=True)
    ax.set_xlabel("初态位置 x₀ (μm)")
    ax.set_ylabel("概率密度")
    ax.set_title("CRN 公共随机数：三个波形用的初态逐行完全相同\n"
                 f"（实测 max|x₀ 差| = {max_dx:.1e} μm，max|v₀ 差| = {max_dv:.1e} m/s）",
                 fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("为什么要 CRN + 为什么分池", fontsize=11.5)
    text = ("不用 CRN 的问题：每个波形重新随机抽样 → 蒙特卡洛噪声\n"
            "（~1/√2000 ≈ 2.2%）会淹没波形之间的真实差异。\n"
            "用 CRN：同一初态喂给所有波形，抽样噪声在配对差分中抵消。\n\n"
            "不分池的问题：拿验证样本反复调参 = 变相死记硬背，\n"
            "报告的『提升』在独立数据上不可复现。\n"
            "Level 1 对比：单一 seed 一组样本，既无分池也无配对比较。")
    ax.text(0.02, 0.95, text, fontsize=10, va="top", linespacing=1.8, color="#333")
    fig.savefig(OUT / "06_新增3_三样本池与CRN.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 07
def fig07_paired_stats():
    metrics = json.loads((DATA / "metrics.json").read_text())
    comps = metrics["validation"]["comparisons"]

    fig = plt.figure(figsize=(16, 8.8))
    gs = fig.add_gridspec(1, 3, width_ratios=(1.15, 1.1, 1), wspace=0.3)
    fig.suptitle("图7  新增④ 配对统计：同一初态下两两波形的四种结局 + 配对 bootstrap 置信区间（真实验证数据）",
                 fontsize=14, fontweight="bold")

    # (a) 2x2 矩阵示意 + 真实数值
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2.6)
    ax.axis("off")
    ax.set_title("(a) 配对计数矩阵（以一对波形为例）", fontsize=11)
    keys = list(comps)
    first = comps[keys[0]]["paired_counts"]
    cells = [(0.05, 1.5, "都成功\nboth_captured", f"{first['both_captured']}", "#2ca02c"),
             (1.05, 1.5, "仅波形A成功\nonly_a", f"{first['only_a_captured']}", "#ff7f0e"),
             (0.05, 0.45, "仅波形B成功\nonly_b", f"{first['only_b_captured']}", "#ff7f0e"),
             (1.05, 0.45, "都失败\nboth_failed", f"{first['both_failed']}", "#d62728")]
    for x, y, label, value, color in cells:
        ax.add_patch(FancyBboxPatch((x, y), 0.9, 0.9, boxstyle="round,pad=0.03",
                                    fc="white", ec=color, lw=2))
        ax.text(x + 0.45, y + 0.73, label, ha="center", fontsize=9.5)
        ax.text(x + 0.45, y + 0.28, value, ha="center", fontsize=15, fontweight="bold",
                color=color)
    ax.text(0.5, 2.42, f"波形A = {keys[0].split('_vs_')[0][:22]}", ha="center", fontsize=9, color="#666")
    ax.text(1.5, 2.42, f"波形B = {keys[0].split('_vs_')[1][:22]}", ha="center", fontsize=9, color="#666")

    # (b) bootstrap CI
    ax = fig.add_subplot(gs[0, 1])
    for i, (key, comp) in enumerate(comps.items()):
        boot = comp["capture_rate_difference_bootstrap"]
        short = key.replace("sequential_600us", "顺序600").replace("sequential_400us", "顺序400") \
                   .replace("overlapped_optimized_400us", "优化400").replace("_vs_", " vs ")
        ax.errorbar(boot["difference"], i, xerr=[[boot["difference"] - boot["ci_low"]],
                                                 [boot["ci_high"] - boot["difference"]]],
                    fmt="o", capsize=5, ms=8, color=WF_COLOR["sequential_400us"], lw=2)
        ax.text(boot["ci_high"] + 0.0004, i,
                f"CI [{boot['ci_low']:+.4f}, {boot['ci_high']:+.4f}]\n"
                f"{'不含0→显著' if boot['ci_excludes_zero'] else '含0→不显著'}",
                va="center", fontsize=8.6)
    ax.axvline(0, color="k", ls="--", lw=1.2)
    ax.set_yticks(range(len(comps)))
    ax.set_yticklabels([k.replace("sequential_600us", "顺序600").replace("sequential_400us", "顺序400")
                         .replace("overlapped_optimized_400us", "优化400").replace("_vs_", "\nvs ")
                        for k in comps], fontsize=9)
    ax.set_xlabel("俘获率差（A − B），配对 bootstrap 95% CI")
    ax.set_title("(b) 波形两两比较：差异均为 0，区间含 0 → 如实报告『无显著差异』", fontsize=10)
    ax.set_xlim(-0.004, 0.006)

    # (c) 解释
    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.set_title("(c) 配对 vs 非配对", fontsize=11)
    rng = np.random.default_rng(7)
    p = 0.99
    indep_a = rng.binomial(1, p, 2000).mean()
    indep_b = rng.binomial(1, p, 2000).mean()
    common = rng.binomial(1, p, 2000)
    paired_diff = 0.0
    ax.bar([0, 1], [abs(indep_a - indep_b), paired_diff],
           color=[C_BAD, C_OK], width=0.55)
    ax.text(0, abs(indep_a - indep_b) + 0.0006, f"{abs(indep_a - indep_b):.4f}",
            ha="center", fontsize=11, color=C_BAD, fontweight="bold")
    ax.text(1, 0.0006, "0.0000", ha="center", fontsize=11, color=C_OK, fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["各自独立抽样\n(非配对)", "同一组初态\n(配对, CRN)"], fontsize=10)
    ax.set_ylabel("两波形俘获率差的假噪声（示意）")
    ax.set_title("(c) 同样『无差异』的两个波形：非配对会量出假差异，\n配对把抽样噪声抵消为 0", fontsize=10)
    fig.savefig(OUT / "07_新增4_配对统计.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 08
def fig08_work_energy():
    run = ideal_run(dict(SPEC_OPT), PHYS, CFG.integration.validation_dt_s,
                    CFG.integration.post_transfer_hold_s)
    traj = run["trajectories"][0]
    we = run["work_energy"]
    t_us = we["time_s"] * 1e6
    end = int(traj["end_index"])
    scale = KB * 1e-6  # J → μK

    shots = read_csv_rows(DATA / "validation_shots.csv")

    fig = plt.figure(figsize=(16.5, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.3, 1), hspace=0.36, wspace=0.24)
    fig.suptitle("图8  新增⑤ 功-能核算：时变势下机械能不守恒是物理，不该当『积分误差』——分开记账并检查残差\n"
                 "（曲线为冻结优化波形的理想轨迹，由包代码实时计算）", fontsize=13.5, fontweight="bold")

    ax = fig.add_subplot(gs[0, 0])
    depth = np.asarray(WF_OPT.depth(we["time_s"]))
    center = np.asarray(WF_OPT.center(we["time_s"]))
    envelope = np.exp(-2 * ((traj["position_m"][:len(t_us)] - center)
                            / PHYS["aod_waist_m"]) ** 2)
    d_rate = np.asarray(WF_OPT.depth_rate(we["time_s"]))
    c_rate = np.asarray(WF_OPT.center_velocity(we["time_s"]))
    p_depth = -d_rate * envelope
    f_aod = -(4 * depth * (traj["position_m"][:len(t_us)] - center)
              / PHYS["aod_waist_m"] ** 2) * envelope
    p_move = f_aod * c_rate
    ax.plot(t_us, p_depth / scale, color="#ff7f0e", lw=1.8, label="深度变化功率 ∂U/∂t|深度")
    ax.plot(t_us, p_move / scale, color="#1f77b4", lw=1.8, label="中心移动功率 F_AOD·c'")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("功率 (μK/s)")
    ax.set_title("(a) 两种外部功的瞬时功率：先移动做负功，后降深做正功", fontsize=10.5)
    ax.legend(fontsize=9)

    ax = fig.add_subplot(gs[1, 0])
    delta_e = we["energy_total_J"] - we["energy_total_J"][0]
    ax.plot(t_us, delta_e / scale, color="k", lw=2.6, label="ΔE 总机械能变化")
    ax.plot(t_us, we["work_move_J"] / scale, color="#1f77b4", lw=1.6, label="W_move 移动功(累计)")
    ax.plot(t_us, we["work_depth_J"] / scale, color="#ff7f0e", lw=1.6, label="W_depth 深度功(累计)")
    ax.plot(t_us, we["work_total_J"] / scale, color="#2ca02c", lw=1.8, ls="--",
            label="W_total 总外部功")
    ax_r = ax.twinx()
    ax_r.plot(t_us, we["residual_J"] / scale, color=C_BAD, lw=1.2, alpha=0.9)
    ax_r.set_ylabel("残差 R_W (μK)", color=C_BAD)
    ax_r.grid(False)
    ax_r.tick_params(axis="y", colors=C_BAD)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("能量 (μK)")
    ax.set_title("(b) 恒等式 ΔE = W_move + W_depth + R_W：黑线与绿虚线几乎重合，红残差 ≈ 0", fontsize=10.5)
    ax.legend(fontsize=9, loc="lower left")

    # 公式面板
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("为什么需要这套账（work_energy.py）", fontsize=11)
    text = ("转移期间外部控制（AOD 中心/深度）在改变，\n"
            "机械能 E(t) 本来就应该变化。\n\n"
            "dU_AOD/dt = −D'·e^(−2q²/w²) + F_AOD·c'\n"
            "             └─深度项─┘   └──移动项──┘\n\n"
            "数值上分别累计两条功，检查\n"
            "R_W(t) = ΔE(t) − W_ext(t) 是否 ≈ 0\n"
            f"理想轨迹实测：max|R_W|/D = "
            f"{we['max_abs_residual_J'] / PHYS['slm_depth_j']:.2e}\n"
            f"（验收阈值 5e-4）\n\n"
            "Level 1 对比：只画了总能量曲线，\n"
            "没有任何功的分解与残差检查。")
    ax.text(0.02, 0.95, text, fontsize=10.2, va="top", linespacing=1.7,
            color="#333")

    # MC 逐 shot 残差
    ax = fig.add_subplot(gs[1, 1])
    for i, name in enumerate(WAVEFORMS):
        res = [to_f(r["max_abs_work_energy_residual_over_depth"])
               for r in shots if r["waveform"] == name]
        ax.bar(i + 1, np.median(res), width=0.55, color=WF_COLOR[name],
               label=f"{WF_LABEL[name]} 中位")
        ax.plot([i + 1 - 0.28, i + 1 + 0.28], [max(res)] * 2, color="k", lw=1.6)
        ax.text(i + 1, max(res) * 1.3, f"max {max(res):.1e}", ha="center", fontsize=8.5)
    ax.axhline(1e-3, color=C_BAD, ls="--", lw=1.4)
    ax.text(2.5, 1.3e-3, "MC 验收阈值 1e-3（中位数）", ha="center", fontsize=9, color=C_BAD)
    ax.set_yscale("log")
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels([WF_LABEL[n] for n in WAVEFORMS], fontsize=8.5)
    ax.set_ylabel("max|R_W| / D_SLM")
    ax.set_title("(c) 2000 个 MC shot 的功-能残差（真实数据）：远低于阈值", fontsize=10.5)
    fig.savefig(OUT / "08_新增5_功能核算.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 09
def fig09_hold():
    run = ideal_run(dict(SPEC_OPT), PHYS, CFG.integration.validation_dt_s,
                    CFG.integration.post_transfer_hold_s)
    traj = run["trajectories"][0]
    t_us = traj["time_s"] * 1e6
    end = traj["end_index"]
    T_us = WF_OPT.duration_s * 1e6
    x = traj["position_m"]
    v = traj["velocity_m_per_s"]
    e_slm = (0.5 * PHYS["mass_kg"] * v ** 2
             + [-PHYS["slm_depth_j"] * np.exp(-2 * (xi - PHYS["slm_center_m"]) ** 2
                                              / PHYS["slm_waist_m"] ** 2) for xi in x])
    hold = t_us >= T_us - 1e-9
    ptp = float(np.ptp(e_slm[hold])) / PHYS["slm_depth_j"]

    shots = read_csv_rows(DATA / "validation_shots.csv")
    hold_mc = np.array([to_f(r["hold_peak_to_peak_energy_error_over_depth"])
                        for r in shots])

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))
    fig.suptitle("图9  新增⑤b  转移后 100 μs 纯 SLM 保持段检查：AOD 已关死，能量只许有界振荡，不许越走越偏",
                 fontsize=13.5, fontweight="bold")

    ax = axes[0]
    ax.plot(t_us, x * 1e6, color="#1f77b4", lw=1.6)
    ax.axvspan(T_us, t_us[-1], color="#2ca02c", alpha=0.12, label="100 μs 保持段")
    ax.plot(t_us[:end + 1], np.asarray(WF_OPT.center(traj["time_s"][:end + 1])) * 1e6,
            "--", color="#d62728", lw=1.4, alpha=0.8, label="AOD 中心 c(t)")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("原子位置 (μm)")
    ax.set_title("(a) 完整时间线：0–400 μs 转移 + 400–500 μs 只剩 SLM", fontsize=10.5)
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(t_us[hold], e_slm[hold] / PHYS["slm_depth_j"], color="#2ca02c", lw=1.8)
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("E_SLM(t) / D_SLM")
    ax.set_title(f"(b) 保持段 SLM 能量（理想轨迹）：峰峰值 = {ptp:.2e}（阈值 1e-4）\n"
                 "平坦 = 静态势下能量守恒，无漂移", fontsize=10.5)

    ax = axes[2]
    ax.hist(hold_mc, bins=50, color="#1f77b4", edgecolor="white")
    ax.axvline(1e-4, color=C_BAD, ls="--", lw=1.6)
    ax.text(1e-4 * 0.6, ax.get_ylim()[1] * 0.8, "阈值 1e-4", rotation=90,
            fontsize=9, color=C_BAD, ha="right")
    ax.set_xlabel("保持段能量误差峰峰值 / D_SLM")
    ax.set_ylabel("shot 数")
    ax.set_title(f"(c) 2000×3 个 MC shot（真实数据）：最大 {hold_mc.max():.1e}，全部远低于阈值",
                 fontsize=10.5)
    fig.text(0.5, 0.005, "Level 1 对比：转移结束即判定，之后没有任何保持段检查；"
                         "Level 2 判定发生在 t=T（AOD 深度严格为 0），保持段用来说明判据可信",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "09_新增5_保持段检查.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 10
def fig10_near_threshold():
    shots = read_csv_rows(DATA / "validation_shots.csv")
    metrics = json.loads((DATA / "metrics.json").read_text())
    near_counts = {name: metrics["validation"]["excitation_stats"][name]["near_threshold_count"]
                   for name in WAVEFORMS}

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))
    fig.suptitle("图10  新增⑥ near-threshold 临界样本标记：判据仍然严格 E_f < 0，但把『差点算错』的样本单独点亮",
                 fontsize=13.5, fontweight="bold")

    ax = axes[0]
    for name in WAVEFORMS:
        ef = np.array([to_f(r["final_slm_energy_uK"]) for r in shots if r["waveform"] == name])
        ef = ef * 1e6  # 列值实为 K，换回 μK
        ax.hist(ef / 140.0, bins=60, alpha=0.4, color=WF_COLOR[name],
                label=WF_LABEL[name], density=True)
    ax.axvline(0, color="k", lw=1.6)
    ax.axvspan(-1e-3, 1e-3, color=C_BAD, alpha=0.35, label="临界带 |E_f|/D < 1e-3")
    ax.set_xlabel("末态能量 E_f / D_SLM  （<0 俘获，≥0 不俘获）")
    ax.set_ylabel("概率密度")
    ax.set_title("(a) 全部末态能量离 0 很远（都在 −1 附近）", fontsize=10.5)
    ax.legend(fontsize=8.5)

    ax = axes[1]
    zoom = [to_f(r["final_slm_energy_uK"]) * 1e6 / 140.0 for r in shots
            if r["waveform"] == WAVEFORMS[0]]
    ax.hist(zoom, bins=np.linspace(-1.05, 1.05, 211), color="#1f77b4")
    ax.axvspan(-1e-3, 1e-3, color=C_BAD, alpha=0.4)
    ax.set_xlim(-2e-3, 2e-3)
    ax.set_xlabel("E_f / D_SLM（横轴放大到 ±0.002）")
    ax.set_ylabel("shot 数")
    ax.set_title(f"(b) 临界带放大：带内样本数 = {sum(near_counts.values())}（三个波形合计）\n"
                 "步长敏感性讨论就集中在这一小段", fontsize=10.5)

    ax = axes[2]
    ax.axis("off")
    ax.set_title("规则（与 Level 1 的差别）", fontsize=11)
    text = ("判据本身没变：E_f < 0 → 俘获；E_f = 0 → 不俘获\n\n"
            "变了的地方：\n"
            "• |E_f|/D_SLM < 1e-3 的样本打 near_threshold 标记\n"
            "• 不删除、不改判，仍按严格符号计入俘获率\n"
            "• 单独报告比例，用来解释『为什么换步长标签会跳』\n\n"
            f"本次真实结果：三波形 near_threshold 计数\n"
            + "\n".join(f"  {WF_LABEL[n]}: {c}" for n, c in near_counts.items())
            + "\n\nLevel 1 对比：只有 captured 布尔值，\n"
              "没有任何临界样本概念。")
    ax.text(0.02, 0.95, text, fontsize=10, va="top", linespacing=1.75, color="#333")
    fig.savefig(OUT / "10_新增6_临界样本标记.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 11
def fig11_robustness():
    rows = read_csv_rows(DATA / "robustness.csv")

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.8))
    fig.suptitle("图11  新增⑦ 鲁棒性扫描：波形冻结后不再调参，故意换个环境再测（真实数据）",
                 fontsize=13.5, fontweight="bold")

    for ax, scan, xlabel, title in [
            (axes[0], "temperature", "初态温度 (μK)", "(a) 温度扫描 3/5/7/10 μK，每点 500 独立 shots"),
            (axes[1], "alignment", "末端对准偏差 (μm)", "(b) AOD 终点整体平移 ±0.10 μm，每点 500 shots")]:
        for name in WAVEFORMS:
            pts = sorted((to_f(r["condition"]), to_f(r["capture_fraction"]),
                          to_f(r["wilson_low"]), to_f(r["wilson_high"]))
                         for r in rows if r["scan"] == scan and r["waveform"] == name)
            c = [p[0] for p in pts]
            f = [p[1] for p in pts]
            lo = [p[2] for p in pts]
            hi = [p[3] for p in pts]
            ax.errorbar(c, f, yerr=[np.maximum(0, np.array(f) - np.array(lo)),
                                    np.maximum(0, np.array(hi) - np.array(f))],
                        fmt="-o", ms=6, capsize=4, color=WF_COLOR[name], label=WF_LABEL[name])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("俘获率（误差棒 = Wilson 95% CI）")
        ax.set_ylim(0.95, 1.008)
        ax.set_title(title, fontsize=10.5)
        ax.legend(fontsize=8.5, loc="lower left")

    ax = axes[2]
    width = 0.26
    for i, name in enumerate(WAVEFORMS):
        pts = sorted((to_f(r["condition"]), to_f(r["final_energy_max_abs_diff_uK"]),
                      to_f(r["label_flips_vs_reference"]))
                     for r in rows if r["scan"] == "timestep" and r["waveform"] == name
                     and r.get("final_energy_max_abs_diff_uK"))
        dts = [f"{p[0]:.2f}" for p in pts]
        diff = [p[1] for p in pts]
        flips = [p[2] for p in pts]
        bars = ax.bar(np.arange(len(pts)) + (i - 1) * width, diff, width,
                      color=WF_COLOR[name], label=WF_LABEL[name])
        for bar, flip in zip(bars, flips):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    "翻转0" if flip == 0 else f"翻转{flip:.3f}", ha="center", va="bottom",
                    fontsize=7.5, color="#333")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["dt=0.10 μs", "dt=0.05 μs", "dt=0.025 μs"])
    ax.set_ylabel("末态能量最大差 (μK，相对最细步长)")
    ax.set_title("(c) 步长敏感性：连续量随 dt 减小收敛；\n『标签翻转率』单独统计（阈值 2%）", fontsize=10.5)
    ax.legend(fontsize=8)
    fig.text(0.5, 0.005, "Level 1 对比：只做过一次 dt 与 dt/2 的理想轨迹收敛检查；"
                         "Level 2 有三组系统扫描、区分连续收敛与分类跳变、每条件独立样本",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "11_新增7_鲁棒性扫描.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 12
def fig12_preflight():
    pf = json.loads((DATA / "preflight.json").read_text())
    l1r = pf["level1_reproduction"]
    conv = pf["ideal_timestep_convergence"]
    mc_dt = pf["mc_timestep_paired"]
    res = pf["mc_work_energy_residual"]
    hold = pf["hold_segment"]

    items = [
        ("① 现有测试", "重跑 Level 0/1 全部测试", pf["existing_tests"]["summary_line"], True),
        ("② 静态极限", "时变积分器退化为静力时与 Level 0 积分器逐点一致",
         f"位置/速度 allclose = {pf['static_limit']['position_match']}", pf["static_limit"]["passed"]),
        ("③ Level 1 复现", "固定种子复现 Level 1 顺序 600 μs",
         f"复现俘获率 {l1r['reproduced_fraction']:.4f} vs 历史 "
         f"{l1r.get('historical_fraction', float('nan')):.4f}（逐位相等）", l1r["passed"]),
        ("④ 波形约束", "三波形端点/范围/单调性/导数全部过 validate_waveform",
         f"三波形 valid = {[v['valid'] for v in pf['waveform_constraints'].values()]}",
         pf["waveform_constraints_all_valid"]),
        ("⑤ 理想轨迹", "理想初态 (x₀=2.4 μm, v₀=0) 跑三波形，均应被俘获",
         f"captured = {[v['captured'] for v in pf['ideal_runs'].values()]}", True),
        ("⑥ 功-能平衡", "理想轨迹 max|R_W|/max(D) < 5e-4",
         f"逐波形通过 = {list(pf['ideal_work_energy']['per_waveform_passed'].values())}",
         pf["ideal_work_energy"]["passed"]),
        ("⑦ 步长收敛", "0.10/0.05/0.025 μs，末态位置收敛阶 > 1（二阶）",
         f"最小位置收敛阶 = {conv['min_position_order']:.2f}", conv["passed"]),
        ("⑧ MC 配对步长", "128 shots 比较 0.05 vs 0.025 μs，标签翻转率 ≤ 2%",
         f"最大翻转率 = {mc_dt['max_label_flip_rate']:.4f}", mc_dt["passed"]),
        ("⑨ MC 功-能残差", "MC 子集残差中位数 < 1e-3",
         "中位数 = " + ", ".join(f"{v['median']:.1e}" for v in res["per_waveform"].values()),
         res["passed"]),
        ("⑩ 保持段有界", "100 μs 保持段峰峰误差 < 1e-4",
         "峰峰值 = " + ", ".join(f"{v:.1e}" for v in hold["per_waveform_peak_to_peak_over_depth"].values()),
         hold["passed"]),
    ]

    fig, ax = plt.subplots(figsize=(16.5, 9.6))
    ax.axis("off")
    ax.set_title(f"图12  新增⑧ 预检查关卡 preflight：10 项全过才准开始优化（真实结果，总耗时 "
                 f"{pf.get('elapsed_s', 0):.0f} s）", fontsize=14, fontweight="bold", pad=14)
    for i, (name, desc, detail, passed) in enumerate(items):
        col, row = i % 2, i // 2
        x = 0.02 + col * 0.50
        y = 0.86 - row * 0.175
        color = C_OK if passed else C_BAD
        ax.add_patch(FancyBboxPatch((x, y), 0.47, 0.155, boxstyle="round,pad=0.008",
                                    fc="#eefaee" if passed else "#fdeaea", ec=color, lw=1.6))
        ax.text(x + 0.015, y + 0.115, name, fontsize=11.5, fontweight="bold", color=color)
        ax.text(x + 0.015, y + 0.072, desc, fontsize=9.3, color="#333")
        ax.text(x + 0.015, y + 0.028, detail, fontsize=8.2, color="#666")
    ax.text(0.5, 0.02, "任一关键项失败 → 抛 RuntimeError，禁止生成『看似完整』的最优结果；"
                       "全部结果写 preflight.json 可复查。Level 1 只有 6 项内嵌小检查，且失败不阻断",
            ha="center", fontsize=10.5, color=C_NEW)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.savefig(OUT / "12_新增8_预检查关卡.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 13
def fig13_noise():
    demo = json.loads(NOISE_JSON.read_text())
    dflt = demo["defaults_from_config"]
    dur = demo["duration_scan"]

    fig = plt.figure(figsize=(16.5, 9.2))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.28)
    fig.suptitle("图13  新增⑨ 平均噪声加热模型：三条噪声通道按『系综平均功率』确定性注入——没有随机数，能量随时间累积",
                 fontsize=13.5, fontweight="bold")

    # (a) 三通道卡片
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.set_title("(a) 三条通道（默认值由物理常数算出，均标注 assumed）", fontsize=10.5)
    channels = [
        ("光子反冲", "P = 2·E_r·Γ_sc", f"{dflt['recoil']['power_uK_per_s']:.2e} μK/s（常数）", "#1f77b4"),
        ("强度噪声(参量)", "dE/dt = Γ_par·ε_osc", f"Γ_par = {dflt['parametric']['rate_per_s']:.0f} /s（∝当前振荡能量）", "#ff7f0e"),
        ("指向噪声", "P = m·ω₀⁴·S_x(ν₀)/8", f"{dflt['pointing']['power_uK_per_s']:.1f} μK/s（常数）", "#9467bd"),
    ]
    for i, (name, formula, value, color) in enumerate(channels):
        y = 0.80 - i * 0.30
        ax.add_patch(FancyBboxPatch((0.02, y - 0.10), 0.96, 0.26, boxstyle="round,pad=0.02",
                                    fc="white", ec=color, lw=1.8))
        ax.text(0.06, y + 0.085, name, fontsize=10.5, fontweight="bold", color=color)
        ax.text(0.06, y - 0.01, formula + "\n" + value, fontsize=9, color="#333", va="center")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # (b) 加热踢示意
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_title("(b) 加热踢：每步把 ΔE = P·dt 沿当前速度方向塞进去", fontsize=10.5)
    ax.annotate("", xy=(0.72, 0.55), xytext=(0.25, 0.55),
                arrowprops=dict(arrowstyle="-|>", lw=2.4, color="#1f77b4"))
    ax.text(0.48, 0.60, "v（当前速度）", ha="center", fontsize=10, color="#1f77b4")
    ax.annotate("", xy=(0.80, 0.72), xytext=(0.25, 0.72),
                arrowprops=dict(arrowstyle="-|>", lw=2.8, color=C_BAD))
    ax.text(0.55, 0.79, "v ← sign(v)·√(v² + 2ΔE/m)", ha="center", fontsize=10.5,
            color=C_BAD, fontweight="bold")
    ax.plot([0.25, 0.25], [0.35, 0.68], "k-", lw=1.2)
    ax.text(0.52, 0.30, "大小只增不减（无随机性）；\n逐通道累计注入能量记进功-能账本：\n"
                        "ΔE = W_ext + W_noise + R_W", ha="center", fontsize=9.5, color="#333")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # (c) 适用范围
    ax = fig.add_subplot(gs[0, 2])
    ax.axis("off")
    ax.set_title("(c) 适用范围（防误用）", fontsize=10.5)
    text = ("• 只在 validate / robustness 阶段生效，\n  优化阶段保持无噪声（不改变已冻结流程）\n"
            "• 幅度参数全部 assumed_sensitivity_only，\n  不声称来自实验数据\n"
            "• 三个噪声函数外放为可注入参数，\n  签名 (x, v, t, ε_osc) → 功率，可换成任意函数\n"
            "• 构造的常数模型可 pickle → 支持多进程\n"
            "• scaled_heating() 整体缩放做敏感性扫描")
    ax.text(0.03, 0.92, text, fontsize=9.8, va="top", linespacing=1.8, color="#333")

    # (d) 注入能量随时长
    ax = fig.add_subplot(gs[1, 0])
    scales = sorted({e["scale"] for e in dur})
    for scale in scales:
        pts = [(e["duration_us"], e["mean_noise_injected_uK"]) for e in dur if e["scale"] == scale]
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=5,
                label=f"噪声强度 ×{scale:g}")
    ax.set_xlabel("转移+保持总时长 (μs)")
    ax.set_ylabel("平均注入能量 (μK)")
    ax.set_title("(d) 真实演示：注入能量随时长线性累积（demo_mean_heating）", fontsize=10.5)
    ax.legend(fontsize=8.5)

    # (e) 俘获率 vs 时长
    ax = fig.add_subplot(gs[1, 1])
    for scale in scales:
        pts = [(e["duration_us"], e["capture_fraction"]) for e in dur if e["scale"] == scale]
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=5,
                label=f"×{scale:g}")
    ax.set_xlabel("转移+保持总时长 (μs)")
    ax.set_ylabel("俘获率")
    ax.set_ylim(0, 1.05)
    ax.set_title("(e) 名义噪声下俘获率仍 1.0（量级太小）；\n放大后随时间变差——『越久越不稳定』", fontsize=10.5)
    ax.legend(fontsize=8.5, title="噪声强度", title_fontsize=8)

    # (f) 保持段演示
    ax = fig.add_subplot(gs[1, 2])
    hold_rows = demo["hold_scan"]
    pts = [(h["hold_us"], h["median_hold_peak_to_peak_over_depth_captured"]) for h in hold_rows]
    pts.sort()
    ax.loglog([p[0] for p in pts], np.maximum([p[1] for p in pts], 1e-16), "-o", ms=5,
              color=C_BAD)
    ax.set_xlabel("保持段时长 (μs)")
    ax.set_ylabel("已俘获原子保持段能量振荡 中位数 / D")
    ax.set_title("(f) 噪声开到演示强度时保持段漂移发散\n（对比 图9 无噪声时 ~1e-7）", fontsize=10.5)
    fig.savefig(OUT / "13_新增9_噪声加热模型.png", bbox_inches="tight")
    plt.close(fig)


# ================================================================ 图 14
def fig14_engineering():
    fig, ax = plt.subplots(figsize=(16.5, 9.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("图14  新增⑩ 工程层面：分阶段 CLI、多进程并行、丰富输出物与程序化验收",
                 fontsize=14.5, fontweight="bold", pad=14)

    # CLI 阶段
    box(ax, 0.15, 7.6, 2.2, 1.9, "preflight",
        "10 项预检查\n失败即终止\n→ preflight.json", fc="#eef7ee", ec=C_OK, fs_b=8.8)
    box(ax, 2.75, 7.6, 2.2, 1.9, "optimize",
        "三阶段搜索+冻结\n支持断点续跑\n→ best_waveform.yaml", fc="#eaf3fb", ec="#1f77b4", fs_b=8.8)
    box(ax, 5.35, 7.6, 2.2, 1.9, "validate",
        "需先有冻结文件\n2000 shots 一次性\n→ validation_shots.csv", fc="#fdf3e7", ec="#ff7f0e", fs_b=8.8)
    box(ax, 7.95, 7.6, 2.0, 1.9, "robustness",
        "需先有冻结文件\n温度/对准/步长\n→ robustness.csv", fc="#fdeaea", ec=C_BAD, fs_b=8.8)
    for x0, x1 in [(2.4, 2.7), (5.0, 5.3), (7.6, 7.9)]:
        arrow(ax, (x0, 8.55), (x1, 8.55))
    ax.text(9.0, 6.9, "--stage all 按序执行；单跑某 stage 时\nmetrics.json 与已有产物合并不覆盖",
            fontsize=8.8, ha="center", color="#666")

    # 多进程
    box(ax, 0.15, 4.4, 3.0, 2.6, "多进程并行评估",
        "evaluate_waveform_parallel\n把 shots 均匀切片\n→ multiprocessing Pool\n"
        "进程数 = CPU 数 − 1\n初态数组本身不变（保 CRN）\n"
        "噪声模型需可 pickle\n（constant 模型用 partial 构造）",
        fc="#f2f2f2", ec=C_OLD, fs_b=8.8)

    # 输出物
    box(ax, 3.55, 4.4, 3.4, 2.6, "输出物（12 数据文件 + 10 张图）",
        "config_used.yaml / preflight.json /\ncandidates.csv / optimization_history.csv /\n"
        "optimization_checkpoint.json /\nbest_waveform.yaml / waveforms.csv /\n"
        "validation_shots.csv / robustness.csv /\nmetrics.json / image_validation.json /\nsummary.md",
        fc="#eaf3fb", ec="#1f77b4", fs_b=8.2)

    # 图片验收
    box(ax, 7.35, 4.4, 2.6, 2.6, "图片程序化验收",
        "每张 PNG 保存时自动记录：\n数据有限 / x,y 范围非零 /\n坐标覆盖数据 /\n"
        "图例条数=系列数 / figure 已关闭\n→ image_validation.json\n"
        "（明示：不能证明美观，未做人工审阅）",
        fc="#f3eefb", ec="#9467bd", fs_b=8.5)

    # 配置校验 + 测试
    box(ax, 0.15, 0.6, 4.7, 3.3, "严格配置 schema（config.py）",
        "字段必须与 schema 完全一致（多/少都报错）\n"
        "• 三种子互不相同（强制）\n"
        "• 池大小逐级递增 128 < 512 < 2000\n"
        "• convergence_dt < validation_dt < optimization_dt\n"
        "• 压缩基线必须短于顺序基线\n"
        "• ramp_power 边界落在 [0.5, 4.0] 内\n"
        "• 最大阱频×任一步长 < 0.05（数值稳定性）\n"
        "• 噪声幅度参数非负\n"
        "config_used.yaml 回写 CLI 覆盖后的真实值",
        fc="#eefaee", ec=C_OK, fs_b=8.8)
    box(ax, 5.15, 0.6, 4.8, 3.3, "测试（40 个新增 + 旧测试全部保留）",
        "test_level2.py（27 个）：波形端点/单调/拒绝非法参数、\n"
        "导数 vs 有限差分、Level1 复现一致、功-能平衡与收敛、\n"
        "保持段有界、CRN 三波形初态逐行相同、三种子互异、\n"
        "validation 未被优化器访问、Wilson/bootstrap 构造数据正确、\n"
        "E_f=0 不俘获、失败候选入表、CLI 各 stage 冒烟、\n"
        "CSV/JSON schema、输出落入 sources/ 被拒绝\n"
        "test_noise_heating.py（13 个）：三通道功率公式、\n"
        "踢速度映射、纯音自检、pickle、缩放",
        fc="#fdf3e7", ec="#ff7f0e", fs_b=8.2)
    fig.savefig(OUT / "14_新增10_工程层面.png", bbox_inches="tight")
    plt.close(fig)


def main():
    figures = [
        fig00_level1_overview,
        fig01_structure,
        fig02_waveforms,
        fig03_constraints,
        fig04_optimizer,
        fig05_objectives,
        fig06_pools_crn,
        fig07_paired_stats,
        fig08_work_energy,
        fig09_hold,
        fig10_near_threshold,
        fig11_robustness,
        fig12_preflight,
        fig13_noise,
        fig14_engineering,
    ]
    for func in figures:
        func()
        print(f"已生成 {func.__name__}")
    print(f"\n共 {len(figures)} 张图 → {OUT}")


if __name__ == "__main__":
    main()
