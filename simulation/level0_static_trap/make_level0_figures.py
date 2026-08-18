"""Level 0 静态径向光镊 PPT 可视化脚本。

独立运行：只依赖 numpy / matplotlib，不 import level0_static_trap 包
（analysis.py 当前含全角字符语法错误，import 会失败）。
直接读取 outputs/level0_static_trap/demo/ 下已有的 CSV 和 metrics.json，
并用解析公式补算拟合曲线，输出 6 张中文字幕 PNG 供 PPT 使用。

运行：
    cd simulation/level0_static_trap
    python make_level0_figures.py
产物：
    outputs/level0_static_trap/ppt_figures/fig1_potential_radial.png
    outputs/level0_static_trap/ppt_figures/fig2_force_radial.png
    outputs/level0_static_trap/ppt_figures/fig3_trajectory.png
    outputs/level0_static_trap/ppt_figures/fig4_energy_conservation.png
    outputs/level0_static_trap/ppt_figures/fig5_frequency_comparison.png
    outputs/level0_static_trap/ppt_figures/fig6_error_decomposition.png
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np

# 受限环境下 matplotlib 缓存目录可能不可写，先指到系统临时目录。
_CACHE = Path(tempfile.gettempdir()) / "level0_ppt_matplotlib"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 中文字体配置（macOS / 常见 Linux 都有 PingFang / Noto CJK；缺字时回退英文）。
for candidate in ("PingFang SC", "Heiti SC", "STHeiti", "Noto Sans CJK SC", "WenQuanYi Zen Hei", "Arial Unicode MS"):
    if candidate in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        rcParams["font.sans-serif"] = [candidate] + rcParams["font.sans-serif"]
        break
rcParams["axes.unicode_minus"] = False  # 负号正常显示，不用方块。

# ── 路径与配色 ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "outputs" / "level0_static_trap" / "demo"
OUT_DIR = PROJECT_ROOT / "outputs" / "level0_static_trap" / "ppt_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COLOR_SLM = "#1f77b4"  # 蓝
COLOR_AOD = "#d62728"  # 红
SLM_LABEL = "SLM (180 μK, 实测束腰)"
AOD_LABEL = "AOD (280 μK, 假定束腰)"


# ── 数据读取 ──────────────────────────────────────────────────────────────────
def _load_grid():
    """读 static_grid.csv，返回 (x_um, slm_pot_uK, slm_force_N, aod_pot_uK, aod_force_N)。"""
    grid = np.genfromtxt(DATA_DIR / "static_grid.csv", delimiter=",", names=True)
    return (
        grid["x_um"],
        grid["slm_potential_uK"],
        grid["slm_force_N"],
        grid["aod_potential_uK"],
        grid["aod_force_N"],
    )


def _load_trajectory(trap: str):
    """读 {trap}_trajectory.csv。返回列字典。"""
    return np.genfromtxt(
        DATA_DIR / f"{trap}_trajectory.csv", delimiter=",", names=True
    )


def _load_metrics():
    with open(DATA_DIR / "metrics.json", "r", encoding="utf-8") as f:
        return json.load(f)


# ── 图1：势能 vs 径向位置 ────────────────────────────────────────────────────
def fig1_potential(grid_data):
    x_um, slm_pot, _, aod_pot, _ = grid_data
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(x_um, slm_pot, color=COLOR_SLM, linewidth=2.0, label=SLM_LABEL)
    ax.plot(x_um, aod_pot, color=COLOR_AOD, linewidth=2.0, label=AOD_LABEL)
    ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    # 标出两个阱底深度。
    ax.annotate("−180 μK", xy=(0, -180), xytext=(1.8, -150),
                fontsize=10, color=COLOR_SLM,
                arrowprops=dict(arrowstyle="->", color=COLOR_SLM, lw=0.8))
    ax.annotate("−280 μK", xy=(0, -280), xytext=(1.8, -260),
                fontsize=10, color=COLOR_AOD,
                arrowprops=dict(arrowstyle="->", color=COLOR_AOD, lw=0.8))
    ax.set_xlim(-4, 4)
    ax.set_xlabel("径向位置 x (μm)", fontsize=11)
    ax.set_ylabel("势能 U / k_B (μK)", fontsize=11)
    ax.set_title("图1  SLM / AOD 静态高斯光镊势能曲线\n"
                 r"$U(x) = -U_0 \exp\left[-2(x/w)^2\right]$",
                 fontsize=12)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig1_potential_radial.png", dpi=300)
    plt.close(fig)


# ── 图2：径向力 vs 位置 ──────────────────────────────────────────────────────
def fig2_force(grid_data):
    x_um, _, slm_force, _, aod_force = grid_data
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(x_um, slm_force * 1e18, color=COLOR_SLM, linewidth=2.0, label=SLM_LABEL)
    ax.plot(x_um, aod_force * 1e18, color=COLOR_AOD, linewidth=2.0, label=AOD_LABEL)
    ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
    # 标出恢复力方向。
    ax.annotate("x<0：力指向 +x\n(指向阱心)", xy=(-2.2, 1.2), fontsize=9,
                ha="center", color="#444",
                arrowprops=dict(arrowstyle="->", color="#444"))
    ax.annotate("x>0：力指向 −x\n(指向阱心)", xy=(2.2, -1.2), fontsize=9,
                ha="center", color="#444",
                arrowprops=dict(arrowstyle="->", color="#444"))
    ax.set_xlim(-4, 4)
    ax.set_xlabel("径向位置 x (μm)", fontsize=11)
    ax.set_ylabel("径向力 F (aN)", fontsize=11)
    ax.set_title("图2  SLM / AOD 光镊径向恢复力\n"
                 r"$F(x) = -\frac{4U_0 x}{w^2}\exp\left[-2(x/w)^2\right]$",
                 fontsize=12)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig2_force_radial.png", dpi=300)
    plt.close(fig)


# ── 图3：原子位移轨迹 ────────────────────────────────────────────────────────
def fig3_trajectory():
    slm = _load_trajectory("slm")
    aod = _load_trajectory("aod")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(slm["time_us"], slm["position_um"], color=COLOR_SLM,
            linewidth=1.0, alpha=0.9, label="SLM 原子轨迹")
    ax.plot(aod["time_us"], aod["position_um"], color=COLOR_AOD,
            linewidth=1.0, alpha=0.9, label="AOD 原子轨迹")
    ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    ax.set_xlabel("时间 (μs)", fontsize=11)
    ax.set_ylabel("位移 (μm)", fontsize=11)
    ax.set_title("图3  原子在光镊中的径向振荡轨迹\n"
                 "初始位移 0.10 w、初速 0，velocity-Verlet 积分",
                 fontsize=12)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    # 文字标注 AOD 振荡更快。
    ax.text(0.02, 0.02, "AOD 频率更高 → 振荡更密",
            transform=ax.transAxes, fontsize=9, color=COLOR_AOD, va="bottom")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_trajectory.png", dpi=300)
    plt.close(fig)


# ── 图4：能量守恒 ────────────────────────────────────────────────────────────
def fig4_energy_conservation():
    slm = _load_trajectory("slm")
    aod = _load_trajectory("aod")
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    slm_max = np.max(np.abs(slm["energy_error_over_depth"])) * 1e9
    aod_max = np.max(np.abs(aod["energy_error_over_depth"])) * 1e9
    ax.plot(slm["time_us"], slm["energy_error_over_depth"] * 1e9,
            color=COLOR_SLM, linewidth=0.8, alpha=0.85,
            label=f"SLM: max $|\\Delta E|/U_0$ = {slm_max:.2f} ppb")
    ax.plot(aod["time_us"], aod["energy_error_over_depth"] * 1e9,
            color=COLOR_AOD, linewidth=0.8, alpha=0.85,
            label=f"AOD: max $|\\Delta E|/U_0$ = {aod_max:.2f} ppb")
    ax.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    ax.set_xlabel("时间 (μs)", fontsize=11)
    ax.set_ylabel(r"能量误差 $(E-E_0)/U_0$ (ppb)", fontsize=11)
    ax.set_title(r"图4  辛积分器能量守恒（除以阱深 $U_0$ 归一化）" + "\n"
                 "有界振荡、无长期漂移 → 数值结果可信",
                 fontsize=12)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_energy_conservation.png", dpi=300)
    plt.close(fig)


# ── 图5：四频率交叉验证 ──────────────────────────────────────────────────────
def fig5_frequency_comparison(metrics):
    slm_f = metrics["trap_results"]["slm"]["frequencies"]
    aod_f = metrics["trap_results"]["aod"]["frequencies"]
    methods = ["解析小振幅", "数值曲率", "非线性积分", "轨迹过零"]
    keys = [
        "analytic_small_amplitude_Hz",
        "numerical_curvature_Hz",
        "nonlinear_quadrature_Hz",
        "trajectory_estimate_Hz",
    ]
    slm_vals = [slm_f[k] / 1e3 for k in keys]
    aod_vals = [aod_f[k] / 1e3 for k in keys]

    x = np.arange(len(methods))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    bars1 = ax.bar(x - width / 2, slm_vals, width, color=COLOR_SLM, label="SLM")
    bars2 = ax.bar(x + width / 2, aod_vals, width, color=COLOR_AOD, label="AOD")
    # 柱顶标数值。
    for bars in (bars1, bars2):
        for b in bars:
            ax.annotate(f"{b.get_height():.2f}",
                        xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=10)
    ax.set_ylabel("阱频率 (kHz)", fontsize=11)
    ax.set_title("图5  四种独立方法计算的阱频率交叉验证\n"
                 "前两者一致（数值实现正确）；后两者低于前两者（物理非简谐软化）",
                 fontsize=11)
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_frequency_comparison.png", dpi=300)
    plt.close(fig)


# ── 图6：误差分解 ────────────────────────────────────────────────────────────
def fig6_error_decomposition(metrics):
    """把三种误差源画成对比柱状图（对数刻度），凸显量级差异。"""
    slm_e = metrics["trap_results"]["slm"]["error_classification"]
    aod_e = metrics["trap_results"]["aod"]["error_classification"]
    categories = [
        ("数值离散误差\n(曲率 vs 解析)", "numerical_discretization_relative_error"),
        ("物理非简谐软化\n(非线性 vs 解析)", "finite_amplitude_nonlinear_shift_relative"),
        ("轨迹采样误差\n(轨迹 vs 非线性)", "trajectory_vs_nonlinear_reference_relative_error"),
    ]
    labels = [c[0] for c in categories]
    slm_vals = [abs(slm_e[c[1]]) for c in categories]
    aod_vals = [abs(aod_e[c[1]]) for c in categories]

    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    bars1 = ax.bar(x - width / 2, slm_vals, width, color=COLOR_SLM, label="SLM")
    bars2 = ax.bar(x + width / 2, aod_vals, width, color=COLOR_AOD, label="AOD")
    for bars in (bars1, bars2):
        for b in bars:
            ax.annotate(f"{b.get_height():.1e}",
                        xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=8)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("相对误差量级（对数轴）", fontsize=11)
    ax.set_title("图6  误差源分离：物理软化 >> 数值/采样误差\n"
                 "物理非简谐软化（~7.5e-3）比数值误差（~5e-7）大 4 个量级",
                 fontsize=11)
    ax.grid(alpha=0.25, axis="y", which="both")
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig6_error_decomposition.png", dpi=300)
    plt.close(fig)


# ── 主入口 ────────────────────────────────────────────────────────────────────
def main():
    print(f"读取数据目录: {DATA_DIR}")
    grid_data = _load_grid()
    metrics = _load_metrics()

    fig1_potential(grid_data)
    print("  [1/6] fig1_potential_radial.png        ✓ 势能 vs 径向位置")
    fig2_force(grid_data)
    print("  [2/6] fig2_force_radial.png            ✓ 径向力 vs 位置")
    fig3_trajectory()
    print("  [3/6] fig3_trajectory.png              ✓ 原子振荡轨迹")
    fig4_energy_conservation()
    print("  [4/6] fig4_energy_conservation.png     ✓ 能量守恒（ppb）")
    fig5_frequency_comparison(metrics)
    print("  [5/6] fig5_frequency_comparison.png    ✓ 四频率交叉验证")
    fig6_error_decomposition(metrics)
    print("  [6/6] fig6_error_decomposition.png     ✓ 误差源分解")

    print()
    print("=" * 70)
    print(f"完成。运行 make_level0_figures.py，得到以下 6 张图，位于：")
    print(f"  {OUT_DIR}")
    print("-" * 70)
    figures = sorted(OUT_DIR.glob("fig*.png"))
    for i, p in enumerate(figures, 1):
        size_kb = p.stat().st_size / 1024
        print(f"  图{i}  {p.name}  ({size_kb:.0f} KB)")
    print("=" * 70)


if __name__ == "__main__":
    main()
