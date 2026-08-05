"""使用非交互后端生成可程序化验收的英文标注 PNG。"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from typing import Iterable

# 某些受限运行环境的默认用户缓存目录不可写；缓存明确落到系统临时目录。
_MPL_CACHE = Path(tempfile.gettempdir()) / "level0_static_trap_matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _finish_plot(fig, ax, path: Path, series_count: int, arrays: Iterable[np.ndarray]) -> dict:
    finite = all(np.all(np.isfinite(np.asarray(array))) for array in arrays)
    legend = ax.get_legend()
    legend_count = len(legend.get_texts()) if legend is not None else 0
    x_limits = tuple(float(v) for v in ax.get_xlim())
    y_limits = tuple(float(v) for v in ax.get_ylim())
    lines = ax.get_lines()
    x_values = np.concatenate([np.asarray(line.get_xdata(), dtype=float) for line in lines])
    y_values = np.concatenate([np.asarray(line.get_ydata(), dtype=float) for line in lines])
    data_x_range = (float(np.min(x_values)), float(np.max(x_values)))
    data_y_range = (float(np.min(y_values)), float(np.max(y_values)))
    axes_cover_data = (
        x_limits[0] <= data_x_range[0] <= data_x_range[1] <= x_limits[1]
        and y_limits[0] <= data_y_range[0] <= data_y_range[1] <= y_limits[1]
    )
    fig.tight_layout()
    figure_number = fig.number
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return {
        "path": str(path),
        "data_finite": bool(finite),
        "data_x_range_nonzero": data_x_range[1] > data_x_range[0],
        "data_y_range_nonzero": data_y_range[1] > data_y_range[0],
        "axes_cover_data": bool(axes_cover_data),
        "series_count": series_count,
        "legend_label_count": legend_count,
        "legend_matches_series": legend_count == series_count,
        "figure_closed": not plt.fignum_exists(figure_number),
    }


def plot_potential_comparison(grid_m, traps: dict, path: Path) -> dict:
    """绘制所选静态势阱的势能曲线。"""

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    arrays = [np.asarray(grid_m)]
    for key, result in traps.items():
        values = result["arrays"]["potential_grid_J"] / 1.380649e-29
        arrays.append(values)
        source_label = "assumed waist" if result["config"].waist_is_assumed else "measured waist"
        ax.plot(np.asarray(grid_m) * 1e6, values, label=f"{key.upper()} ({source_label})")
    ax.set(xlabel="Radial position x (um)", ylabel="Potential U/kB (uK)", title="Static Gaussian trap potential")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    return _finish_plot(fig, ax, path, len(traps), arrays)


def plot_force_comparison(grid_m, traps: dict, path: Path) -> dict:
    """绘制所选静态势阱的径向力。"""

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    arrays = [np.asarray(grid_m)]
    for key, result in traps.items():
        values = result["arrays"]["force_grid_N"] * 1e18
        arrays.append(values)
        ax.plot(np.asarray(grid_m) * 1e6, values, label=key.upper())
    ax.set(xlabel="Radial position x (um)", ylabel="Force (aN)", title="Static Gaussian trap force")
    ax.grid(alpha=0.25)
    ax.legend()
    return _finish_plot(fig, ax, path, len(traps), arrays)


def plot_single_trajectory(key: str, result: dict, path: Path) -> dict:
    """绘制单个势阱的径向位移轨迹。"""

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    time_us = result["arrays"]["time_s"] * 1e6
    displacement_um = (result["arrays"]["position_m"] - result["config"].center_m) * 1e6
    ax.plot(time_us, displacement_um, linewidth=1.0, label=f"{key.upper()} trajectory")
    ax.set(xlabel="Time (us)", ylabel="Displacement from center (um)", title=f"{key.upper()} static-trap trajectory")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    return _finish_plot(fig, ax, path, 1, [time_us, displacement_um])


def plot_energy_evolution(key: str, result: dict, path: Path) -> dict:
    """绘制以正势阱深度归一化的总能量误差。"""

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    time_us = result["arrays"]["time_s"] * 1e6
    error = result["arrays"]["energy_error_over_depth"]
    error_ppb = error * 1.0e9
    ax.plot(time_us, error_ppb, linewidth=0.7, alpha=0.85, label="(E - E0) / U0 x 1e9")
    ax.axhline(0.0, color="black", linewidth=0.7, alpha=0.45)
    ax.set(xlabel="Time (us)", ylabel="Energy error / U0 (ppb)", title=f"{key.upper()} bounded energy error")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    zoom = ax.inset_axes([0.67, 0.64, 0.30, 0.29])
    zoom_end = min(float(time_us[-1]), 100.0)
    mask = time_us <= zoom_end
    zoom.plot(time_us[mask], error_ppb[mask], linewidth=0.8)
    zoom.axhline(0.0, color="black", linewidth=0.5, alpha=0.4)
    zoom.set_title("First 100 us", fontsize=8)
    zoom.tick_params(labelsize=7)
    zoom.grid(alpha=0.2)
    return _finish_plot(fig, ax, path, 1, [time_us, error_ppb])


def plot_frequency_comparison(traps: dict, path: Path) -> dict:
    """比较解析、数值曲率、非线性积分和轨迹频率。"""

    fig, ax = plt.subplots(figsize=(7.4, 4.7))
    labels = ["Analytic", "Curvature", "Nonlinear", "Trajectory"]
    x = np.arange(len(labels), dtype=float)
    arrays: list[np.ndarray] = [x]
    for key, result in traps.items():
        frequencies = result["metrics"]["frequencies"]
        values = np.array([
            frequencies["analytic_small_amplitude_Hz"], frequencies["numerical_curvature_Hz"],
            frequencies["nonlinear_quadrature_Hz"], frequencies["trajectory_estimate_Hz"],
        ]) / 1e3
        arrays.append(values)
        ax.plot(x, values, marker="o", linewidth=1.5, label=key.upper())
    ax.set_xticks(x, labels)
    ax.set(xlabel="Frequency estimate", ylabel="Frequency (kHz)", title="Trap-frequency consistency and finite-amplitude shift")
    ax.grid(alpha=0.25)
    ax.legend()
    return _finish_plot(fig, ax, path, len(traps), arrays)
