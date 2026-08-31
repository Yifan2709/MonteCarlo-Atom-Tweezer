"""Level 2 诊断图：非交互后端，保存后关闭，附带程序化验收记录。"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

_MPL_CACHE = Path(tempfile.gettempdir()) / "level2_joint_transfer_matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_AOD = "#1f77b4"
C_SLM = "#d62728"
C_OPT = "#2ca02c"
C_BASE = "#7f7f7f"


def _finish(fig, ax, path, series_count, arrays):
    """保存 PNG 并生成与 level0.io_utils.validate_pngs 兼容的验收记录。

    图例与数据范围统计遍历整个 figure 的所有子图，适配多面板图。
    """
    finite = all(np.all(np.isfinite(np.asarray(array))) for array in arrays)
    legend_count = sum(len(axis.get_legend().get_texts()) for axis in fig.axes
                       if axis.get_legend() is not None)
    cover, x_nonzero, y_nonzero = True, True, True
    for axis in fig.axes:
        # 只统计 data 坐标下的曲线；axhline/axvline 参考线使用混合坐标，跳过
        lines = [line for line in axis.get_lines()
                 if line.get_transform() is axis.transData]
        if not lines:
            continue
        x_values = np.concatenate([np.asarray(line.get_xdata(), dtype=float) for line in lines])
        y_values = np.concatenate([np.asarray(line.get_ydata(), dtype=float) for line in lines])
        if not (np.all(np.isfinite(x_values)) and np.all(np.isfinite(y_values))):
            continue
        x_limits = tuple(float(v) for v in axis.get_xlim())
        y_limits = tuple(float(v) for v in axis.get_ylim())
        cover = cover and (x_limits[0] <= float(np.min(x_values)) and float(np.max(x_values)) <= x_limits[1]
                           and y_limits[0] <= float(np.min(y_values)) and float(np.max(y_values)) <= y_limits[1])
        x_nonzero = x_nonzero and (float(np.max(x_values)) > float(np.min(x_values)))
        y_nonzero = y_nonzero and (float(np.max(y_values)) > float(np.min(y_values)))
    fig.tight_layout()
    number = fig.number
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return {
        "path": str(path), "data_finite": bool(finite),
        "data_x_range_nonzero": bool(x_nonzero), "data_y_range_nonzero": bool(y_nonzero),
        "axes_cover_data": bool(cover), "series_count": series_count,
        "legend_label_count": legend_count, "legend_matches_series": legend_count == series_count,
        "figure_closed": not plt.fignum_exists(number),
    }


def plot_waveform_comparison(waveform_tables: dict, path: Path):
    """三类主波形的中心与深度随时间变化，标出重叠区。"""
    fig, axs = plt.subplots(2, 1, figsize=(9.5, 7.0), sharex=True)
    arrays = []
    series = 0
    for name, table in waveform_tables.items():
        t_us = table["time_s"] * 1e6
        axs[0].plot(t_us, table["aod_center_m"] * 1e6, label=name)
        axs[1].plot(t_us, table["aod_depth_J"] / 1.380649e-29, label=name)
        arrays += [t_us, table["aod_center_m"], table["aod_depth_J"]]
        series += 1
    axs[0].set_ylabel("AOD center (um)")
    axs[1].set_ylabel("AOD depth (uK)")
    axs[1].set_xlabel("Time (us)")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend()
    axs[0].set_title("Waveform comparison (move/ramp overlap highlighted in optimizer)")
    return _finish(fig, axs[1], path, 2 * series, arrays)


def plot_potential_evolution(spec, physics, path, num_snapshots=6):
    """优化波形下若干关键时刻的总势剖面。"""
    from .level2_simulation import build_waveform
    from level0_static_trap.potentials import gaussian_potential
    waveform = build_waveform(spec, physics)
    times = np.linspace(0.0, waveform.duration_s, num_snapshots)
    x = np.linspace(-3e-6, 3.2e-6, 1201)
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    arrays = [x * 1e6]
    series = 0
    for t in times:
        center = float(waveform.center(t))
        depth = float(waveform.depth(t))
        total = (gaussian_potential(x, physics["slm_depth_j"], physics["slm_waist_m"], physics["slm_center_m"])
                 + gaussian_potential(x, depth, physics["aod_waist_m"], center))
        ax.plot(x * 1e6, total / 1.380649e-29, label=f"t = {t*1e6:.0f} us")
        arrays.append(total)
        series += 1
    ax.axvline(0.0, color=C_SLM, ls=":", lw=1.2, label="SLM center")
    ax.set(xlabel="x (um)", ylabel="Total potential / kB (uK)", title="Potential evolution (optimized waveform)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    return _finish(fig, ax, path, series + 1, arrays)


def plot_representative_trajectories(trajectories: dict, path: Path):
    """成功、失败与阈值附近代表轨迹，叠加 AOD 中心。"""
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    arrays = []
    series = 0
    for label, trajectory in trajectories.items():
        time_us = trajectory["time_s"] * 1e6
        end = trajectory["end_index"]
        ax.plot(time_us, trajectory["position_m"] * 1e6, label=f"atom ({label})")
        ax.plot(time_us, trajectory["aod_center_m"] * 1e6, ls="--", alpha=0.55)
        arrays += [time_us, trajectory["position_m"] * 1e6]
        series += 2
    ax.set(xlabel="Time (us)", ylabel="Position (um)", title="Representative trajectories vs AOD center")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    return _finish(fig, ax, path, series // 2, arrays)


def plot_final_phase_space(records: list, physics, path):
    """末态 (x,v) 散点按 captured 着色，并画 E_SLM=0 分界线。"""
    from level0_static_trap.potentials import gaussian_potential
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    captured = [r for r in records if r["captured"]]
    missed = [r for r in records if not r["captured"]]
    arrays = []
    series = 0
    for group, color, label in ((captured, C_OPT, "captured"), (missed, C_SLM, "not captured")):
        if group:
            xs = np.array([r["final_x_um"] for r in group])
            vs = np.array([r["final_v_m_s"] for r in group])
            ax.scatter(xs, vs, s=14, alpha=0.6, color=color, label=label)
            arrays += [xs, vs]
            series += 1
    v_grid = np.linspace(-0.12, 0.12, 400)
    x_grid = np.linspace(-3.0, 3.0, 400)
    boundary = []
    for v in v_grid:
        kinetic = 0.5 * physics["mass_kg"] * v**2
        target = -kinetic
        depths = gaussian_potential(x_grid * 1e-6, physics["slm_depth_j"], physics["slm_waist_m"], 0.0)
        index = int(np.argmin(np.abs(depths - target)))
        boundary.append(x_grid[index])
    ax.plot(boundary, v_grid, "k--", lw=1.2, label="E_SLM = 0")
    ax.set(xlabel="Final x (um)", ylabel="Final v (m/s)", title="Final phase space")
    ax.grid(alpha=0.3)
    ax.legend()
    return _finish(fig, ax, path, series + 1, arrays + [v_grid, np.array(boundary)])


def plot_capture_rate_comparison(validation_summaries: dict, training_summaries: dict, path):
    """validation 俘获率与 Wilson 区间；训练结果以浅色区分。"""
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    names = list(validation_summaries)
    x = np.arange(len(names))
    arrays = [x]
    series = 0
    for offset, summaries, alpha, label in ((-0.16, training_summaries, 0.45, "training/selection"),
                                            (0.16, validation_summaries, 1.0, "independent validation")):
        positions, rates, lows, highs = [], [], [], []
        for position, name in enumerate(names):
            if name not in summaries:
                continue
            entry = summaries[name]
            if entry.get("capture_fraction") is None or entry.get("wilson_low") is None:
                continue
            positions.append(position + offset)
            rates.append(entry["capture_fraction"])
            lows.append(max(entry["capture_fraction"] - entry["wilson_low"], 0.0))
            highs.append(max(entry["wilson_high"] - entry["capture_fraction"], 0.0))
        if not positions:
            continue
        ax.bar(positions, rates, width=0.3, alpha=alpha, label=label,
               yerr=[lows, highs], capsize=4,
               color=[C_BASE, C_AOD, C_OPT][:len(names)])
        arrays += [np.array(positions), np.array(rates)]
        series += 1
    ax.set_xticks(x, names, fontsize=9)
    ax.set(ylabel="Capture fraction", ylim=(0, 1.1), title="Capture rate: training vs independent validation")
    ax.grid(alpha=0.3, axis="y")
    ax.legend()
    return _finish(fig, ax, path, series, arrays)


def plot_excitation_and_margin(records_by_waveform: dict, path):
    """已俘获原子的末态激发、俘获裕量与激发能变化分布。"""
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6))
    arrays = []
    series = 0
    legend_expected = 0
    for panel_index, (ax, key, label) in enumerate((
            (axs[0], "final_excitation_over_depth", "Final excitation / depth"),
            (axs[1], "capture_margin", "Capture margin"),
            (axs[2], "excitation_change_uK", "Excitation change (uK)"))):
        for name, records in records_by_waveform.items():
            values = [r[key] for r in records if r["captured"]]
            if values:
                ax.hist(values, bins=24, alpha=0.55, label=name)
                arrays.append(np.array(values))
                series += 1
                if panel_index == 0:
                    legend_expected += 1
        ax.set(xlabel=label, ylabel="Shots")
        ax.grid(alpha=0.2)
    axs[0].legend(fontsize=9)
    fig.suptitle("Captured-shot excitation and margin distributions", fontsize=13)
    return _finish(fig, axs[2], path, legend_expected, arrays)


def plot_optimization_history(history_rows: list, path):
    """候选评估顺序、最佳值演化与无效候选。"""
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    scores = np.array([r.get("objective_total") if r.get("objective_total") is not None else np.nan
                       for r in history_rows], dtype=float)
    index = np.arange(len(history_rows))
    valid = np.isfinite(scores)
    # marker-only Line2D（而非 scatter）：候选点进入 _finish 的数据范围检查，
    # 否则 running-best 为水平线时 data_y_range_nonzero 会误报失败
    ax.plot(index[valid], scores[valid], "o", ms=5, alpha=0.6, color=C_AOD,
            label="valid candidate")
    ax.plot(index[~valid], np.zeros(int((~valid).sum())), "o", ms=5, alpha=0.6, color=C_SLM,
            label="invalid/failed")
    if valid.any():
        # 从首个有效候选起画累计最优，避免前导无效候选产生 -inf 污染 data_finite 检查
        first_valid = int(np.argmax(valid))
        tail_valid = valid[first_valid:]
        tail_scores = scores[first_valid:]
        running = np.maximum.accumulate(np.where(tail_valid, tail_scores, -np.inf))
        ax.plot(index[first_valid:], running, color=C_OPT, lw=1.8, label="running best")
        arrays = [index, running]
    else:
        arrays = [index]
    ax.set(xlabel="Evaluation order", ylabel="Objective total", title="Optimization history")
    ax.grid(alpha=0.3)
    ax.legend()
    return _finish(fig, ax, path, 3, arrays + [index[valid], scores[valid]])


def plot_work_energy_balance(work_energy_results: dict, path):
    """多条轨迹的 DeltaE、移动功、深度功、总功与残差随时间。"""
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    scale = 1.380649e-29
    colors = ["#4c72b0", "#dd8452"]
    arrays = []
    series = 0
    for index, (label, we) in enumerate(work_energy_results.items()):
        t_us = we["time_s"] * 1e6
        delta_e = we["energy_total_J"] - we["energy_total_J"][0]
        series_items = [
            (f"ΔE ({label})", delta_e, "black", 1.8),
            (f"move work ({label})", we["work_move_J"], C_AOD, 1.2),
            (f"depth work ({label})", we["work_depth_J"], "#ff7f0e", 1.2),
            (f"total work ({label})", we["work_total_J"], C_OPT, 1.2),
            (f"residual R_W ({label})", we["residual_J"], C_SLM, 1.2),
        ]
        arrays.append(t_us)
        for name, values, color, lw in series_items:
            ax.plot(t_us, values / scale, color=color, lw=lw, ls="-" if index == 0 else "--",
                    label=name, alpha=0.55 + 0.45 * (index == 0))
            arrays.append(values)
            series += 1
    ax.set(xlabel="Time (us)", ylabel="Energy / kB (uK)", title="Work-energy balance")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return _finish(fig, ax, path, series, arrays)


def plot_timestep_convergence(timestep_rows: list, path):
    """连续量末能差与标签翻转率随步长变化。"""
    by_waveform = {}
    for row in timestep_rows:
        by_waveform.setdefault(row["waveform"], []).append(row)
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.8))
    arrays = []
    series = 0
    for name, rows in by_waveform.items():
        rows.sort(key=lambda r: r["condition"])
        dts = np.array([r["condition"] for r in rows])
        energy_diff = np.array([r["final_energy_max_abs_diff_uK"] for r in rows])
        flips = np.array([r["label_flips_vs_reference"] for r in rows])
        axs[0].semilogy(dts, np.clip(energy_diff, 1e-12, None), "o-", label=name)
        axs[1].plot(dts, flips, "s-", label=name)
        arrays += [dts, energy_diff, flips]
        series += 1
    all_dts = np.array(sorted({r["condition"] for r in timestep_rows}))
    axs[1].plot([all_dts.min(), all_dts.max()], [0.02, 0.02], color="gray", ls=":", lw=1.2)  # 翻转容差参考线
    axs[0].set(xlabel="dt (us)", ylabel="Max |dE_final| (uK)", title="Continuous convergence")
    axs[1].set(xlabel="dt (us)", ylabel="Label flip rate", title="Classification jumps")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    return _finish(fig, axs[1], path, 2 * series, arrays)


def plot_robustness_heatmaps(robustness_rows: list, path):
    """温度与末端对准偏差下的俘获率两组敏感性图。"""
    fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.8))
    arrays = []
    series = 0
    for ax, scan, xlabel in ((axs[0], "temperature", "Temperature (uK)"),
                             (axs[1], "alignment", "Final alignment offset (um)")):
        rows = [r for r in robustness_rows if r["scan"] == scan]
        waveforms = sorted({r["waveform"] for r in rows})
        conditions = sorted({r["condition"] for r in rows})
        lookup = {(r["waveform"], r["condition"]): r["capture_fraction"] for r in rows}
        for index, name in enumerate(waveforms):
            rates = [lookup.get((name, condition), np.nan) for condition in conditions]
            ax.plot(conditions, rates, "o-", label=name)
            arrays.append(np.array(rates))
            series += 1
        wilson_floor = min(float(r["wilson_low"]) for r in rows if str(r["wilson_low"]) not in ("", "None"))
        ax.plot([min(conditions), max(conditions)], [wilson_floor, wilson_floor],
                color="gray", ls=":", lw=1.2)  # 该样本量下的 Wilson 下界参考线
        ax.set(xlabel=xlabel, ylabel="Capture fraction", ylim=(0, 1.08),
               title=f"{scan} sensitivity")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    return _finish(fig, axs[1], path, series, arrays + [np.array(conditions)])
