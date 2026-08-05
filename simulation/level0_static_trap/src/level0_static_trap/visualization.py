"""绘图模块：势/力/轨迹/能量/频率比较图。

约定：英文坐标与标题（避免中文字体缺失）；非交互后端（在 ``__init__`` 设置 Agg）；
>=200 dpi；不调用 ``plt.show()``；保存后关闭 figure。每个绘图函数返回一个元数据
字典，供程序化图像检查使用（系列数、图例标签数、数据范围）。
"""

from __future__ import annotations

import math
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

DPI = 200


def _metadata(fig, series_count: int, legend_count: int, x_range: tuple, y_range: tuple) -> dict:
    return {
        "series_count": int(series_count),
        "legend_label_count": int(legend_count),
        "x_data_range": (float(x_range[0]), float(x_range[1])),
        "y_data_range": (float(y_range[0]), float(y_range[1])),
        "x_data_nonzero": bool(abs(x_range[1] - x_range[0]) > 0.0),
        "y_data_nonzero": bool(abs(y_range[1] - y_range[0]) > 0.0),
        "figure_size_inches": (float(fig.get_size_inches()[0]), float(fig.get_size_inches()[1])),
    }


def plot_potential_comparison(
    path: str,
    x_um: np.ndarray,
    traps: dict[str, dict],
) -> dict:
    """绘制势能比较图。``traps`` 形如 ``{name: {"uK": array, "color": str}}``。

    返回图像元数据用于程序化检查。
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    series = 0
    legends = 0
    ymin = math.inf
    ymax = -math.inf
    for name, info in traps.items():
        ax.plot(x_um, info["uK"], color=info.get("color"), lw=1.6, label=f"{name} potential")
        series += 1
        legends += 1
        ymin = min(ymin, float(np.min(info["uK"])))
        ymax = max(ymax, float(np.max(info["uK"])))
    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel("Radial position x (um)")
    ax.set_ylabel("Potential U (uK)")
    ax.set_title("Static radial Gaussian optical dipole potential")
    if legends:
        ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    meta = _metadata(fig, series, legends, (float(x_um[0]), float(x_um[-1])), (ymin, ymax))
    plt.close(fig)
    return meta


def plot_force_comparison(
    path: str,
    x_um: np.ndarray,
    traps: dict[str, dict],
) -> dict:
    """绘制力比较图。``traps`` 形如 ``{name: {"force_N": array, "color": str}}``。"""
    fig, ax = plt.subplots(figsize=(8, 5))
    series = 0
    legends = 0
    ymin = math.inf
    ymax = -math.inf
    for name, info in traps.items():
        f = info["force_N"]
        ax.plot(x_um, f, color=info.get("color"), lw=1.6, label=f"{name} force")
        series += 1
        legends += 1
        ymin = min(ymin, float(np.min(f)))
        ymax = max(ymax, float(np.max(f)))
    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel("Radial position x (um)")
    ax.set_ylabel("Force F (N)")
    ax.set_title("Static radial force F = -dU/dx")
    if legends:
        ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    meta = _metadata(fig, series, legends, (float(x_um[0]), float(x_um[-1])), (ymin, ymax))
    plt.close(fig)
    return meta


def plot_single_trajectory(
    path: str,
    times_us: np.ndarray,
    displacement_um: np.ndarray,
    trap_name: str,
) -> dict:
    """绘制单势阱经典轨迹（位移随时间）。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(times_us, displacement_um, color="C0", lw=1.0, label=f"{trap_name} x(t) - x_c")
    ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Displacement from trap center (um)")
    ax.set_title(f"{trap_name} classical trajectory (velocity-Verlet)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    meta = _metadata(
        fig,
        1,
        1,
        (float(times_us[0]), float(times_us[-1])),
        (float(np.min(displacement_um)), float(np.max(displacement_um))),
    )
    plt.close(fig)
    return meta


def plot_energy_evolution(
    path: str,
    times_us: np.ndarray,
    dE_over_U0: np.ndarray,
    trap_name: str,
    energy_tol: float | None = None,
) -> dict:
    """绘制归一化能量误差 ``(E(t)-E0)/U0`` 随时间的演化。

    刻意只画能量误差，避免动能/势能/负总能在同一尺度上掩盖误差。
    """
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(times_us, dE_over_U0, color="C3", lw=1.0, label="(E(t)-E0)/U0")
    ax.axhline(0.0, color="0.4", lw=0.8, ls="--")
    if energy_tol is not None:
        ax.axhline(energy_tol, color="0.6", lw=0.8, ls=":", label=f"+tol ({energy_tol:g})")
        ax.axhline(-energy_tol, color="0.6", lw=0.8, ls=":", label=f"-tol ({energy_tol:g})")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Normalized energy error (E-E0)/U0")
    ax.set_title(f"{trap_name} energy error evolution")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    legend_count = 1 + (2 if energy_tol is not None else 0)
    meta = _metadata(
        fig,
        1 + (2 if energy_tol is not None else 0),
        legend_count,
        (float(times_us[0]), float(times_us[-1])),
        (float(np.min(dE_over_U0)), float(np.max(dE_over_U0))),
    )
    plt.close(fig)
    return meta


def plot_frequency_comparison(
    path: str,
    trap_metrics: dict[str, dict],
) -> dict:
    """绘制四个频率（analytic / numerical curvature / nonlinear quadrature /
    trajectory estimated）的分组柱状比较图。"""
    names = list(trap_metrics.keys())
    categories = [
        "analytic_small_amp_Hz",
        "numerical_curvature_Hz",
        "nonlinear_quadrature_Hz",
        "trajectory_estimated_Hz",
    ]
    cat_labels = ["Analytic small-amp", "Numerical curvature", "Nonlinear quadrature", "Trajectory estimated"]
    fig, ax = plt.subplots(figsize=(9, 5))
    n_groups = len(names)
    n_cat = len(categories)
    width = 0.8 / n_cat
    x = np.arange(n_groups)
    colors = ["C0", "C1", "C2", "C3"]
    series = 0
    legends = 0
    ymin = math.inf
    ymax = -math.inf
    for j, cat in enumerate(categories):
        vals = []
        for name in names:
            f = trap_metrics[name]["frequencies"].get(cat, float("nan"))
            vals.append(float(f) if math.isfinite(f) else 0.0)
            if math.isfinite(f):
                ymin = min(ymin, f)
                ymax = max(ymax, f)
        ax.bar(x + (j - (n_cat - 1) / 2.0) * width, vals, width, label=cat_labels[j], color=colors[j])
        series += 1
        legends += 1
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("Trap frequency comparison")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    if not math.isfinite(ymin):
        ymin, ymax = 0.0, 1.0
    meta = _metadata(fig, series, legends, (float(x[0]) - 0.5, float(x[-1]) + 0.5), (ymin, ymax))
    plt.close(fig)
    return meta
