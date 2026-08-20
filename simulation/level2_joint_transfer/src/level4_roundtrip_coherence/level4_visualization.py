"""Level 4 诊断图：非交互后端、保存后关闭、附程序化验收记录。

在 prompt §15 要求的 13 张图基础上扩展：初态采样、上游冻结来源、preflight
数值诊断、basin 归属、checkpoint 能量演化、分段相位分解、配对比较、长尾
生存、usable coherent fraction、每轮失败事件与噪声模型均有对应图；每一步
计算都可在图中核查。
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

_MPL_CACHE = Path(tempfile.gettempdir()) / "level4_matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Hiragino Sans GB", "Heiti SC",
                        "Microsoft YaHei", "Noto Sans CJK SC", "Arial Unicode MS",
                        "sans-serif"],
    "axes.unicode_minus": False,
})

from level1_transfer_1d.analysis import wilson_interval
from level3_3d_transport_lensing.gaussian_3d import KB
from .dd_sequences import toggling_function

C_A = "#1f77b4"
C_B = "#ff7f0e"
C_C = "#2ca02c"
C_D = "#d62728"
C_E = "#9467bd"
DD_COLOR = {"none": C_A, "spin_echo": C_C, "xy4_per_long_move": C_D}
BASIN_COLORS = {"bound_to_slm": C_A, "bound_to_aod": C_D,
                "shared_or_ambiguous": C_E, "unbound": "0.45"}
BASIN_LABELS = {"bound_to_slm": "SLM 束缚", "bound_to_aod": "AOD 束缚",
                "shared_or_ambiguous": "shared/ambiguous", "unbound": "失束缚"}
STAGES = ["pickup", "split", "outbound_transport", "remote_hold",
          "inbound_transport", "merge", "dropoff", "final_hold"]

# 输出图片清单（文件名 → 用途说明），供 summary.md 与 metrics 引用。
FIGURE_CATALOG = {
    "protocol_timeline.png": "Protocol A/B 分段时间线、checkpoint 与 DD toggling/脉冲",
    "trap_path_and_depths.png": "Protocol A/B 深度、AOD 中心/速度与 lensing 焦点偏移",
    "initial_ensemble.png": "初态六维采样分布与接受率/温度/种子元数据",
    "upstream_provenance.png": "上游冻结输入（路径、SHA-256、时间）一览",
    "preflight_checks.png": "20 项 preflight 状态与关键数值诊断",
    "representative_roundtrips.png": "成功/各失败阶段/recaptured 代表轨迹（位置+能量+功）",
    "checkpoint_survival.png": "四个条件的逐 checkpoint 条件/累计成功概率",
    "basin_classification.png": "各 checkpoint basin 归属、重叠区 E_SLM/E_AOD 散点",
    "checkpoint_energy_evolution.png": "checkpoint 总能量分位数带随协议阶段演变",
    "failure_stage_breakdown.png": "单轮 first-failure 阶段堆叠（按条件）",
    "segment_phase_breakdown.png": "各 DD 模式的分段相位贡献分解",
    "paired_comparisons.png": "配对四格计数与 bootstrap 差异区间",
    "survival_vs_round.png": "10 轮 absorbing/非吸收生存与拟合",
    "long_tail_survival.png": "128-shot 长尾（至 40 轮）absorbing 生存曲线",
    "failure_events_by_round.png": "重复往返每轮新增失败事件的阶段堆叠",
    "heating_vs_round.png": "留阱 shots 的 SLM 激发分位数随轮次",
    "phase_and_contrast.png": "相位分布与各 DD 模式对比度随轮次",
    "usable_coherent_fraction.png": "usable_coherent_fraction = P(success)×C_cond 随轮次",
    "dd_toggling_and_phase.png": "代表轨迹的 y(t)、δω 分解与三种 DD 累计相位",
    "energy_ledger.png": "分段功分解、ΔE–W_ext 对照与残差",
    "noise_model_diagnostic.png": "assumed 强度噪声示例轨迹与统计",
    "ablation_comparison.png": "消融的成功率、激发与 C(none)/C(xy4)",
    "timestep_convergence.png": "状态/功/相位/成功率随步长收敛",
    "sensitivity_panels.png": "温度、深度、对准、lensing、η、噪声与脉冲时刻敏感性",
}


def _finish(fig, path, series_count, arrays, y_nonzero_any_axis=False,
            x_nonzero_any_axis=False):
    """保存 PNG 并生成与 io_utils.validate_pngs 兼容的验收记录。

    y/x_nonzero_any_axis：允许“至少一个子图数据范围非零”（单点、单
    checkpoint 或全同数据的合法场景），默认要求每个含线子图范围非零。
    """
    finite = all(np.all(np.isfinite(np.asarray(a, dtype=float)))
                 for a in arrays if a is not None and np.size(a))
    legend_count = sum(len(ax.get_legend().get_texts()) for ax in fig.axes
                       if ax.get_legend() is not None)
    cover, x_nonzero, y_nonzero = True, True, True
    any_x_nonzero, any_nonzero, saw_lines = False, False, False
    for ax in fig.axes:
        lines = [ln for ln in ax.get_lines() if ln.get_transform() is ax.transData]
        if not lines:
            continue
        xs = np.concatenate([np.asarray(ln.get_xdata(), dtype=float)
                             for ln in lines])
        ys = np.concatenate([np.asarray(ln.get_ydata(), dtype=float)
                             for ln in lines])
        if xs.size == 0 or ys.size == 0:
            continue
        saw_lines = True
        xl, yl = tuple(map(float, ax.get_xlim())), tuple(map(float, ax.get_ylim()))
        cover = cover and xl[0] <= xs.min() and xs.max() <= xl[1] \
            and yl[0] <= ys.min() and ys.max() <= yl[1]
        ax_x_nonzero = xs.max() > xs.min()
        ax_nonzero = ys.max() > ys.min()
        x_nonzero = x_nonzero and ax_x_nonzero
        y_nonzero = y_nonzero and ax_nonzero
        any_x_nonzero = any_x_nonzero or ax_x_nonzero
        any_nonzero = any_nonzero or ax_nonzero
    if y_nonzero_any_axis:
        y_nonzero = any_nonzero or not saw_lines
    if x_nonzero_any_axis:
        x_nonzero = any_x_nonzero or not saw_lines
    fig.tight_layout()
    number = fig.number
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return {"path": str(path), "data_finite": bool(finite),
            "data_x_range_nonzero": bool(x_nonzero),
            "data_y_range_nonzero": bool(y_nonzero), "axes_cover_data": bool(cover),
            "series_count": series_count, "legend_label_count": legend_count,
            "legend_matches_series": legend_count == series_count,
            "figure_closed": not plt.fignum_exists(number)}


# --------------------------------------------------------------- 协议与控制
def plot_protocol_timeline(protocols, checkpoints, path, dd_pulse_sets=None):
    """图1：分段时间线、checkpoint、DD toggling 与脉冲时刻。"""
    n = len(protocols)
    fig, axes = plt.subplots(n, 1, figsize=(13, 3.4 * n), squeeze=False)
    arrays = []
    series = 0
    dd_pulse_sets = dd_pulse_sets or {}
    for row, (name, timeline) in enumerate(protocols.items()):
        ax = axes[row][0]
        for i, seg in enumerate(timeline.segments):
            color = plt.cm.viridis(i / max(1, len(timeline.segments) - 1))
            ax.axvspan(timeline.segment_starts[i] * 1e6,
                       (timeline.segment_starts[i] + seg.duration_s) * 1e6,
                       color=color, alpha=0.35)
            ax.text((timeline.segment_starts[i] + seg.duration_s / 2) * 1e6,
                    0.0, f"{seg.name}\n{seg.duration_s * 1e6:.0f}μs",
                    rotation=90, ha="center", va="center", fontsize=6.5)
        for ck in checkpoints.get(name, []):
            ax.axvline(ck["time_s"] * 1e6, color="k", ls="--", lw=0.8)
        pulses = dd_pulse_sets.get(name, {})
        xy4 = np.asarray(pulses.get("xy4_per_long_move", np.empty(0)))
        if xy4.size:
            grid_t = np.arange(0.0, timeline.total_duration_s + 1e-12, 0.05e-6)
            grid_t = np.clip(grid_t, 0.0, timeline.total_duration_s)
            y = toggling_function(grid_t, xy4)
            ax.plot(grid_t * 1e6, y, color=C_D, lw=0.9, alpha=0.85,
                    label="XY4 toggling y(t)")
            ax.plot(xy4 * 1e6, np.full(xy4.size, 1.18), "^", color=C_D, ms=5,
                    ls="none", label="XY4 π 脉冲")
            arrays += [grid_t, y]
            if row == 0:
                series += 2
        echo = np.asarray(pulses.get("spin_echo", np.empty(0)))
        if echo.size:
            ax.plot(echo * 1e6, np.full(echo.size, -1.18), "v", color=C_C, ms=6,
                    ls="none", label="spin-echo π 脉冲")
            arrays.append(echo)
            if row == 0:
                series += 1
        ax.set_ylim(-1.45, 1.45)
        ax.set_yticks([-1.0, 0.0, 1.0])
        ax.set_title(f"{name}（总时长 {timeline.total_duration_s*1e6:.0f} μs；"
                     "虚线=checkpoint，色带=分段）")
        ax.set_xlabel("时间 (μs)")
        ax.set_xlim(0, timeline.total_duration_s * 1e6)
        if row == 0 and series:
            ax.legend(fontsize=7, loc="upper left")
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_trap_path_and_depths(grid, path, grid_b=None):
    """图2：深度、AOD 中心/速度与 lensing 偏移（可选 Protocol A/B 对比）。"""
    grids = [("Protocol A", grid)] + ([("Protocol B", grid_b)] if grid_b else [])
    n_rows = len(grids)
    fig, axes = plt.subplots(n_rows, 4, figsize=(15, 3.6 * n_rows), squeeze=False)
    arrays = []
    series = 0
    for r, (label, grid) in enumerate(grids):
        t_us = grid["times"] * 1e6
        axes[r][0].plot(t_us, grid["d_slm"] / KB * 1e6, color=C_A,
                        label="SLM 深度")
        axes[r][0].plot(t_us, grid["d_aod"] / KB * 1e6, color=C_D,
                        label="AOD 深度")
        axes[r][0].set_ylabel("深度 (μK)")
        axes[r][1].plot(t_us, grid["c1"] * 1e6, color=C_A, label="c1")
        axes[r][1].plot(t_us, grid["c2"] * 1e6, color=C_C, ls="--", label="c2")
        axes[r][1].set_ylabel("AOD 中心 (μm)")
        axes[r][2].plot(t_us, grid["c1d"], color=C_A, label="ċ1")
        axes[r][2].plot(t_us, grid["c2d"], color=C_C, ls="--", label="ċ2")
        axes[r][2].set_ylabel("AOD 速度 (m/s)")
        axes[r][3].plot(t_us, grid["zs1"] * 1e6, color=C_E, label="zs1")
        axes[r][3].plot(t_us, grid["zs2"] * 1e6, color=C_B, ls="--", label="zs2")
        axes[r][3].set_ylabel("lensing 焦点偏移 (μm)")
        arrays += [grid["d_slm"], grid["d_aod"], grid["c1"], grid["c1d"],
                   grid["zs1"]]
        series += 8
        for c in range(4):
            axes[r][c].set_xlim(t_us[0], t_us[-1])
            axes[r][c].legend(fontsize=7)
            if r == n_rows - 1:
                axes[r][c].set_xlabel("时间 (μs)")
            else:
                axes[r][c].set_xticklabels([])
        axes[r][0].set_title(f"{label}：SLM/AOD 深度", fontsize=10)
        axes[r][1].set_title(f"{label}：AOD 中心", fontsize=10)
        axes[r][2].set_title(f"{label}：AOD 速度", fontsize=10)
        axes[r][3].set_title(f"{label}：lensing 焦点偏移", fontsize=10)
    return _finish(fig, path, series, arrays)


def plot_initial_ensemble(states, path):
    """新增：初态六维采样分布 + 接受率/温度/种子元数据。"""
    pos = np.asarray(states["pos_m"])
    vel = np.asarray(states["vel_m_per_s"])
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6),
                             layout="constrained")
    arrays = []
    series = 0
    specs = [("x", pos[:, 0] * 1e6, "μm"), ("y", pos[:, 1] * 1e6, "μm"),
             ("z", pos[:, 2] * 1e6, "μm"), ("vx", vel[:, 0], "m/s"),
             ("vy", vel[:, 1], "m/s"), ("vz", vel[:, 2], "m/s")]
    for i, (name, data, unit) in enumerate(specs):
        ax = axes[i // 3][i % 3]
        ax.hist(data, bins=40, color=C_A, alpha=0.85)
        mean = float(np.mean(data))
        std = float(np.std(data))
        ax.axvline(mean, color=C_D, ls="--", lw=1.3,
                   label=f"均值 {mean:+.3g}，std {std:.3g}")
        ax.set_title(f"{name} 初态分布（{states['shots']} shots）")
        ax.set_xlabel(f"{name} ({unit})")
        ax.set_ylabel("shot 数")
        ax.legend(fontsize=7)
        arrays.append(data)
        series += 1
    acc = states.get("acceptance_rate")
    acc_str = f"{acc:.3f}" if isinstance(acc, (int, float)) else "n/a"
    fig.suptitle(
        f"初态采样（validation pool 复现）：seed={states.get('seed')}，"
        f"T={states.get('temperature_uK', states.get('temperature'))} μK，"
        f"阱深={states.get('trap_depth_uK', 'n/a')} μK，接受率={acc_str}\n"
        f"{states.get('sampling_note', '')}",
        fontsize=8)
    return _finish(fig, path, series, arrays)


def plot_upstream_provenance(manifest, path):
    """新增：上游冻结输入清单（角色、文件、SHA-256 前缀、修改时间）。"""
    files = manifest.get("files", {})
    roles = list(files)
    fig, ax = plt.subplots(figsize=(13, max(2.2, 0.62 * len(roles) + 1.6)))
    ax.axis("off")
    cell_text = [[role,
                  Path(files[role]["absolute_path"]).name,
                  files[role].get("sha256", "")[:12] + "…",
                  str(files[role].get("size_bytes", "")),
                  (files[role].get("mtime_utc", "") or "")[:19],
                  files[role].get("note", "")]
                 for role in roles]
    table = ax.table(cellText=cell_text,
                     colLabels=["角色", "文件", "SHA-256 前缀", "字节数",
                                "修改时间 (UTC)", "说明"],
                     colWidths=[0.16, 0.24, 0.12, 0.08, 0.20, 0.24],
                     loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.5)
    recon = manifest.get("reconstructed_from_config")
    ax.set_title(f"上游冻结输入 manifest（生成于 "
                 f"{(manifest.get('generated_utc') or '')[:19]} UTC；"
                 f"reconstructed_from_config={recon}）", fontsize=10)
    return _finish(fig, path, 0, [np.arange(4.0)])


def plot_preflight_checks(checks, path):
    """新增：20 项检查状态/耗时 + 关键数值诊断（回归率、时间反向、相位验证）。"""
    names = [n for n in checks.get("critical_checks", []) if n in checks]
    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.0])
    series = 0
    arrays = []
    # (a) 检查耗时与状态
    ax_a = fig.add_subplot(gs[0, 0])
    elapsed = [float(checks[n].get("elapsed_s", 0.0)) for n in names]
    passed = [bool(checks[n].get("passed", False)) for n in names]
    ax_a.barh(np.arange(len(names)), elapsed,
              color=[C_C if p else C_D for p in passed], label="通过")
    ax_a.barh(np.arange(len(names)),
              [0.0 if p else max(e, 0.05) for p, e in zip(passed, elapsed)],
              color=C_D, label="失败")
    ax_a.set_yticks(np.arange(len(names)))
    ax_a.set_yticklabels([f"{'PASS' if p else 'FAIL'} {n}"
                          for n, p in zip(names, passed)], fontsize=7)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("耗时 (s)")
    ax_a.set_title(f"(a) {len(names)} 项检查（总 "
                   f"{checks.get('elapsed_s_total', 0):.0f} s，"
                   f"{'全部通过' if checks.get('all_passed') else '存在失败'}）",
                   fontsize=10)
    ax_a.legend(fontsize=7)
    arrays.append(np.asarray(elapsed, dtype=float))
    series += 2
    # (b) 控制时间反向误差 vs 步长
    ax_b = fig.add_subplot(gs[0, 1])
    rev = checks.get("control_time_reversal", {}).get("errors", {})
    if rev:
        dts = sorted(float(k) for k in rev)
        pos_err = [float(rev[k]["pos_err_m"]) * 1e6 for k in sorted(rev, key=float)]
        vel_err = [float(rev[k]["vel_err_m_s"]) for k in sorted(rev, key=float)]
        ax_b.loglog(dts, pos_err, "o-", color=C_A, label="往返位置误差 (μm)")
        ax_b.loglog(dts, vel_err, "s--", color=C_C, label="往返速度误差 (m/s)")
        arrays += [np.asarray(pos_err), np.asarray(vel_err)]
        series += 2
    ax_b.set_xlabel("dt (μs)")
    ax_b.set_ylabel("反转回初态的误差")
    ax_b.set_title("(b) 控制时间反向 + 速度反转（Verlet 可逆性）", fontsize=10)
    ax_b.legend(fontsize=7)
    # (c) L2/L3 回归率与 contrast 构造案例
    ax_c = fig.add_subplot(gs[1, 0])
    items, vals = [], []
    l2 = checks.get("level2_regression", {})
    l3 = checks.get("level3_regression", {})
    cc = checks.get("ensemble_contrast_cases", {})
    if l2:
        items.append("L2 三维切片\n俘获率")
        vals.append(float(l2.get("capture_rate", np.nan)))
    if l3:
        items.append("L3 长运输\n留阱率")
        vals.append(float(l3.get("retention_fraction", np.nan)))
    for key, label in (("identical_contrast", "contrast 全同相"),
                       ("uniform_contrast", "contrast 均匀相位"),
                       ("two_point_contrast", "contrast 两点相位")):
        if key in cc and cc[key] is not None:
            items.append(label)
            vals.append(float(cc[key]))
    if items:
        x = np.arange(len(items))
        ax_c.bar(x, vals, color=C_A, label="检查值")
        ax_c.axhline(0.99, color=C_D, ls="--", lw=1.0, label="回归阈值 0.99")
        if cc.get("two_point_expected") is not None:
            ax_c.scatter([len(items) - 1], [float(cc["two_point_expected"])],
                         marker="*", s=120, color="k", zorder=3,
                         label="两点案例解析值")
            series += 1
        ax_c.set_xticks(x)
        ax_c.set_xticklabels(items, fontsize=7)
        arrays.append(np.nan_to_num(np.asarray(vals, dtype=float)))
        series += 2
    ax_c.set_ylim(-0.02, 1.12)
    ax_c.set_ylabel("概率 / 对比度")
    ax_c.set_title("(c) Level 2/3 回归率与 ensemble contrast 构造案例", fontsize=10)
    ax_c.legend(fontsize=7)
    # (d) 相位与功-能数值误差（log 刻度条形）
    ax_d = fig.add_subplot(gs[1, 1])
    d_items, d_vals = [], []
    ap = checks.get("analytic_phase", {})
    if ap.get("max_relative_error") is not None:
        d_items.append("解析相位\n相对误差")
        d_vals.append(float(ap["max_relative_error"]))
    pg = checks.get("pulse_grid_handling", {})
    if pg.get("relative_error") is not None:
        d_items.append("脉冲网格拆分\n相对误差")
        d_vals.append(float(pg["relative_error"]))
    dd = checks.get("dd_cancellation", {}).get("per_mode", {})
    for mode, rec in dd.items():
        none_phi = abs(float(rec.get("none_phi", np.nan)))
        dd_phi = abs(float(rec.get("dd_phase", 0.0)))
        if none_phi > 0:
            d_items.append(f"DD 消除残差比\n({mode})")
            d_vals.append(max(dd_phi / none_phi, 1e-16))
    we = checks.get("work_energy_residual", {})
    if we.get("relative_to_depth") is not None:
        d_items.append("功-能残差\n/阱深")
        d_vals.append(abs(float(we["relative_to_depth"])))
    sh = checks.get("static_hold_bounded", {})
    if sh.get("max_w_ext_uK") is not None:
        d_items.append("静态 hold\nw_ext (μK)")
        d_vals.append(max(float(sh["max_w_ext_uK"]), 1e-12))
    if d_items:
        x = np.arange(len(d_items))
        ax_d.bar(x, d_vals, color=C_E, label="数值")
        ax_d.set_xticks(x)
        ax_d.set_xticklabels(d_items, fontsize=7)
        arrays.append(np.asarray(d_vals, dtype=float))
        series += 1
    ax_d.set_yscale("log")
    ax_d.set_ylabel("相对误差 / 比值 / μK")
    ax_d.set_title("(d) 相位解析验证、DD 消除与功-能残差（对数刻度）", fontsize=10)
    ax_d.legend(fontsize=7)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


# ------------------------------------------------------- 单轮 validation 图
def plot_checkpoint_survival(checkpoint_stats, path):
    """图4：逐 checkpoint 条件/累计成功概率（多条件小倍数，含 Wilson 带）。"""
    if isinstance(checkpoint_stats, dict):
        conds = list(checkpoint_stats)
    else:
        conds = ["validation"]
        checkpoint_stats = {conds[0]: checkpoint_stats}
    n = len(conds)
    n_cols = 2 if n > 1 else 1
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.2 * n_cols, 4.4 * n_rows),
                             squeeze=False)
    arrays = []
    series = 0
    for i, cond in enumerate(conds):
        ax = axes[i // n_cols][i % n_cols]
        rows = checkpoint_stats[cond]
        names = [row["checkpoint"] for row in rows]
        x = np.arange(len(names))
        cum = np.array([row["cumulative"] for row in rows], dtype=float)
        cond_p = np.array([row["conditional"] for row in rows], dtype=float)
        low = np.array([row.get("cumulative_wilson_low", row["cumulative"] - 0.05)
                        for row in rows], dtype=float)
        high = np.array([row.get("cumulative_wilson_high",
                                 row["cumulative"] + 0.05)
                         for row in rows], dtype=float)
        ax.fill_between(x, low, high, alpha=0.2, color=C_A, label="Wilson 95% CI")
        ax.plot(x, cum, "o-", color=C_A, label="累计成功")
        ax.plot(x, low, "-", color=C_A, alpha=0.45, lw=0.9, label="Wilson 下界")
        ax.plot(x, cond_p, "s--", color=C_D, label="条件成功（给定前一 ckpt）")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=7)
        ax.set_ylim(-0.02, 1.05)
        ax.set_ylabel("概率")
        ax.set_title(cond, fontsize=9)
        ax.legend(fontsize=7)
        arrays += [cum, cond_p, low, high]
        series += 4
    for j in range(n, n_rows * n_cols):
        fig.delaxes(axes[j // n_cols][j % n_cols])
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True,
                   x_nonzero_any_axis=True)


def plot_failure_stage_breakdown(stage_counts, path):
    """图5：first-failure stage 堆叠图（按协议）。"""
    protocols = list(stage_counts)
    fig, ax = plt.subplots(figsize=(11, 5))
    bottoms = np.zeros(len(protocols))
    x = np.arange(len(protocols))
    series = 0
    arrays = []
    for stage in STAGES:
        vals = np.array([stage_counts[p].get(stage, 0) for p in protocols],
                        dtype=float)
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, bottom=bottoms, label=stage)
        bottoms += vals
        arrays.append(vals)
        series += 1
    none_vals = np.array([max(0.0, stage_counts[p].get("none", 0))
                          for p in protocols], dtype=float)
    ax.bar(x, none_vals, bottom=bottoms, label="none（无失败）", color="lightgray")
    arrays.append(none_vals)
    series += 1
    ax.set_xticks(x)
    ax.set_xticklabels(protocols, rotation=15)
    ax.set_ylabel("shot 数")
    ax.legend(fontsize=8)
    return _finish(fig, path, series, arrays)


def plot_basin_classification(records, path):
    """新增：checkpoint basin 归属堆叠、重叠 checkpoint 的 E_SLM/E_AOD 散点。"""
    from .trap_assignment import BASINS
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.15, 1.0, 0.75])
    names = [rec["checkpoint"]["name"] for rec in records]
    x = np.arange(len(records))
    series = 0
    arrays = []
    # (a) 各 checkpoint 的 basin 计数堆叠
    ax_a = fig.add_subplot(gs[0, :])
    bottoms = np.zeros(len(records))
    for basin in BASINS:
        vals = np.array([np.sum(np.asarray(rec["classification"]["basin"])
                                == basin) for rec in records], dtype=float)
        ax_a.bar(x, vals, bottom=bottoms, color=BASIN_COLORS[basin],
                 label=BASIN_LABELS[basin])
        bottoms += vals
        arrays.append(vals)
        series += 1
    # 期望归属：在对应堆叠柱顶部放星标
    star_x, star_y = [], []
    cum = np.zeros(len(records))
    for k, basin in enumerate(BASINS):
        vals = np.array([np.sum(np.asarray(rec["classification"]["basin"])
                                == basin) for rec in records], dtype=float)
        for i, rec in enumerate(records):
            if rec["checkpoint"]["expected_basin"] == basin:
                star_x.append(i)
                star_y.append(cum[i] + vals[i] + 0.02 * bottoms.max())
        cum += vals
    ax_a.plot(star_x, star_y, "*", color="k", ms=11, ls="none",
              label="期望归属")
    series += 1
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax_a.set_ylabel("shot 数")
    ax_a.set_title("(a) 各 checkpoint 的 basin 归属计数（Protocol A nominal，"
                   "2000 shots）", fontsize=10)
    ax_a.legend(fontsize=8)
    # (b) 两个关键重叠/转移 checkpoint 的 E_AOD vs E_SLM 散点
    scatter_idx = []
    for want in ("split_end", "merge_end", "dropoff_end", "pickup_end"):
        for i, rec in enumerate(records):
            if rec["checkpoint"]["name"] == want and i not in scatter_idx:
                scatter_idx.append(i)
    scatter_idx = scatter_idx[:2] or [0, len(records) - 1]
    for k, i in enumerate(scatter_idx):
        ax_s = fig.add_subplot(gs[1, k])
        rec = records[i]
        cl = rec["classification"]
        e_slm = np.asarray(cl["e_slm"]) / KB * 1e6
        e_aod = np.asarray(cl["e_aod"]) / KB * 1e6
        basins = np.asarray(rec["classification"]["basin"])
        sub = np.arange(len(e_slm))[:: max(1, len(e_slm) // 700)]
        for basin in BASINS:
            mask = basins[sub] == basin
            if not mask.any():
                continue
            ax_s.plot(e_slm[sub][mask], e_aod[sub][mask], "o", ms=2.5,
                      alpha=0.6, color=BASIN_COLORS[basin])
        ax_s.axhline(0, color="0.5", lw=0.8)
        ax_s.axvline(0, color="0.5", lw=0.8)
        arrays += [e_slm, e_aod]
        ax_s.set_xlabel("E_SLM (μK)")
        ax_s.set_ylabel("E_AOD (μK)")
        ax_s.set_title(f"(b{k + 1}) {rec['checkpoint']['name']} "
                       f"（阱间距 {cl['separation_m'] * 1e6:.2f} μm，"
                       f"颜色=basin，与 (a) 图例一致）", fontsize=9)
    # (c) ambiguous / near-threshold 比例
    ax_c = fig.add_subplot(gs[2, :])
    amb = np.array([float(np.mean(rec["classification"]["ambiguous"]))
                    for rec in records])
    near = np.array([float(np.mean(
        rec["classification"]["near_slm_threshold"]
        | rec["classification"]["near_aod_threshold"])) for rec in records])
    ax_c.plot(x, amb, "o-", color=C_E, label="ambiguous 比例")
    ax_c.plot(x, near, "s--", color=C_B, label="near-threshold 比例")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax_c.set_ylabel("比例")
    ax_c.set_title("(c) 无法稳定归类与接近阈值的比例（不得静默删除）", fontsize=10)
    ax_c.legend(fontsize=8)
    arrays += [amb, near]
    series += 2
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_checkpoint_energy_evolution(records, path):
    """新增：checkpoint 总能量分位数带随协议阶段演变。"""
    fig, ax = plt.subplots(figsize=(13, 5.5))
    names = [rec["checkpoint"]["name"] for rec in records]
    expected = [rec["checkpoint"]["expected_basin"] for rec in records]
    x = np.arange(len(records))
    e = np.stack([np.asarray(rec["e_total"]) / KB * 1e6 for rec in records])
    q05, q25, q50, q75, q95 = [np.quantile(e, q, axis=1)
                               for q in (0.05, 0.25, 0.50, 0.75, 0.95)]
    ax.fill_between(x, q05, q95, alpha=0.18, color=C_A, label="5–95% 分位带")
    ax.fill_between(x, q25, q75, alpha=0.35, color=C_A, label="25–75% 分位带")
    ax.plot(x, q50, "o-", color=C_D, lw=1.6, label="中位数")
    ax.axhline(0, color="k", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n({BASIN_LABELS.get(e_, e_)})"
                        for n, e_ in zip(names, expected)],
                       rotation=25, ha="right", fontsize=7.5)
    ax.set_xlabel("checkpoint（括号内为期望归属）")
    ax.set_ylabel("总机械能 E_tot (μK)")
    ax.set_title("checkpoint 总能量分布随协议阶段演变（Protocol A nominal；"
                 "E=0 为连续统阈值）", fontsize=10)
    ax.legend(fontsize=8)
    return _finish(fig, path, 3, [q05, q50, q95], y_nonzero_any_axis=True)


def plot_segment_phase_breakdown(seg_phases, seg_names, success_mask, path):
    """新增：各 DD 模式按分段的相位贡献（success 条件下 mean|φ| ± std）。"""
    modes = list(seg_phases)
    fig, axes = plt.subplots(1, len(modes), figsize=(6.0 * len(modes), 4.6),
                             squeeze=False)
    arrays = []
    series = 0
    mask = None if success_mask is None else np.asarray(success_mask, dtype=bool)
    for k, mode in enumerate(modes):
        ax = axes[0][k]
        mat = np.asarray(seg_phases[mode], dtype=float)
        if mask is not None and mask.size == mat.shape[1]:
            mat = mat[:, mask]
        mean_abs = np.mean(np.abs(mat), axis=1)
        std = np.std(mat, axis=1)
        x = np.arange(len(seg_names))
        ax.bar(x, mean_abs, color=DD_COLOR.get(mode, C_A), alpha=0.85,
               label="mean |φ_seg|")
        ax.errorbar(x, mean_abs, yerr=std, fmt="none", ecolor="k", capsize=2,
                    lw=1.0, label="std(φ_seg)")
        ax.set_xticks(x)
        ax.set_xticklabels(seg_names, rotation=35, ha="right", fontsize=7)
        ax.set_xlabel("分段")
        ax.set_ylabel("分段相位 (rad)")
        ax.set_title(f"DD = {mode}", fontsize=10)
        ax.legend(fontsize=8)
        arrays += [mean_abs, std]
        series += 2
    fig.suptitle("分段相位贡献分解（条件于 roundtrip success；η_SLM=η_AOD=nominal）",
                 fontsize=10)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_paired_comparisons(comparisons, path):
    """新增：同初态配对差异（bootstrap 95% CI）与四格配对计数。"""
    items = list(comparisons)
    n = len(items)
    if not n:
        return None
    fig = plt.figure(figsize=(12, 3.2 * n + 2.6))
    gs = GridSpec(n + 1, 1, figure=fig,
                  height_ratios=[2.2] + [1.0] * n)
    arrays = []
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(n)
    diffs = np.array([float(comparisons[k]["bootstrap"]["difference"])
                      for k in items])
    lows = np.array([float(comparisons[k]["bootstrap"]["ci_low"])
                     for k in items])
    highs = np.array([float(comparisons[k]["bootstrap"]["ci_high"])
                      for k in items])
    yerr = np.vstack([np.clip(diffs - lows, 0.0, None),
                      np.clip(highs - diffs, 0.0, None)])
    y_all = np.concatenate([diffs, lows, highs])
    degenerate = not (y_all.max() > y_all.min())
    if not degenerate:
        ax.bar(x, diffs, yerr=yerr, capsize=4, color=C_A, alpha=0.85,
               label="配对差（bootstrap 95% CI）")
        ax.plot([x[0] - 0.5, x[-1] + 0.5], [0.0, 0.0], color="0.4", ls="--",
                lw=1.0)
    else:
        # 全部配对完全平局（差与 CI 相同）：以文字说明代替零高度柱
        ax.text(0.5, 0.5,
                f"全部配对差与 CI 相同（Δ={diffs[0]:+.4f}，"
                f"CI [{lows[0]:+.4f}, {highs[0]:+.4f}]）",
                transform=ax.transAxes, ha="center", va="center", fontsize=11)
        ax.set_ylim(-1.0, 1.0)
    labels = []
    for i, name in enumerate(items):
        comp = comparisons[name]
        a, b = comp["pair"]
        sig = comp["bootstrap"].get("ci_excludes_zero")
        labels.append(f"{name}\nP({a})−P({b})"
                      f"{'*' if sig else ''}")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("配对成功概率差")
    ax.set_title("同初态配对比较（* = 95% CI 不含 0，统计显著）", fontsize=10)
    if not degenerate:
        ax.legend(fontsize=8)
    arrays.append(diffs)
    for i, name in enumerate(items):
        axt = fig.add_subplot(gs[i + 1, 0])
        axt.axis("off")
        counts = comparisons[name]["paired_counts"]
        a, b = comparisons[name]["pair"]
        cells = [[str(counts["both_captured"]), str(counts["only_a_captured"])],
                 [str(counts["only_b_captured"]), str(counts["both_failed"])]]
        table = axt.table(cellText=cells,
                          rowLabels=[a, b],
                          colLabels=[f"{a} 成功", f"{a} 失败"],
                          cellLoc="center", loc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(7.5)
        table.scale(1.0, 1.4)
        axt.set_title(f"{name}：四格配对计数（{counts['shots']} shots）",
                      fontsize=8)
    # 单组比较（n=1）时误差帽 x 恒定属合法退化，放宽 x 范围检查
    return _finish(fig, path, 0 if degenerate else 1, arrays,
                   x_nonzero_any_axis=True)


# ------------------------------------------------------------ 代表轨迹图
def plot_representative_roundtrips(trajs, path):
    """图3：成功/各失败阶段/recaptured 代表轨迹（位置 + 能量/累积功）。"""
    n = max(1, len(trajs))
    has_energy = any("energy_uK" in t for t in trajs.values())
    fig, axes = plt.subplots(n, 2 if has_energy else 1,
                             figsize=(13 if has_energy else 12, 2.6 * n),
                             squeeze=False)
    arrays = []
    series = 0
    for i, (label, traj) in enumerate(trajs.items()):
        t = traj["time_s"] * 1e6
        pos = traj["position_m"]
        ax = axes[i][0]
        ax.plot(t, pos[:, 0, 0] * 1e6, color=C_A, label="x1")
        ax.plot(t, pos[:, 0, 1] * 1e6, color=C_C, label="x2")
        ax.plot(t, pos[:, 0, 2] * 1e6, color=C_D, label="z")
        ax.set_title(label, fontsize=9)
        ax.set_ylabel("位置 (μm)")
        ax.legend(fontsize=7)
        arrays += [pos[:, 0, 0], pos[:, 0, 1], pos[:, 0, 2]]
        series += 3
        if has_energy:
            axe = axes[i][1]
            if "energy_uK" in traj:
                energy = np.asarray(traj["energy_uK"], dtype=float)
                work = np.asarray(traj.get("work_uK",
                                           np.full(len(energy), np.nan)),
                                  dtype=float)
                axe.plot(t, energy, color=C_E, label="总机械能")
                axe.plot(t, work, color=C_B, ls="--", label="累积外部功")
                axe.axhline(0, color="0.5", lw=0.8)
                arrays += [energy, np.nan_to_num(work)]
                series += 2
            axe.set_ylabel("能量 / 功 (μK)")
            axe.legend(fontsize=7)
    axes[-1][0].set_xlabel("时间 (μs)")
    if has_energy:
        axes[-1][1].set_xlabel("时间 (μs)")
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_dd_toggling(toggling_data, path):
    """图9：三种 DD 模式的 y(t)、δω 的 SLM/AOD 分解与累计 φ 对比。"""
    t = toggling_data["time_s"] * 1e6
    modes = toggling_data.get("modes", {})
    fig, axes = plt.subplots(3, 1, figsize=(12, 9.5), sharex=True)
    arrays = []
    series = 0
    ax0 = axes[0]
    for mode, rec in modes.items():
        y = np.asarray(rec["y"], dtype=float)
        color = DD_COLOR.get(mode, None)
        ax0.step(t, y, where="post", color=color, lw=1.0, alpha=0.9,
                 label=f"y(t)（{mode}）")
        for pulse in np.asarray(rec.get("pulses_s", []), dtype=float):
            ax0.axvline(pulse * 1e6, color=color, ls=":", lw=0.7, alpha=0.6)
        arrays.append(y)
        series += 1
    ax0.set_ylabel("toggling y(t)")
    ax0.set_ylim(-1.35, 1.35)
    ax0.legend(fontsize=8, loc="right")
    ax1 = axes[1]
    dw = np.asarray(toggling_data["delta_omega"], dtype=float)
    dw_slm = np.asarray(toggling_data["delta_omega_slm"], dtype=float)
    dw_aod = np.asarray(toggling_data["delta_omega_aod"], dtype=float)
    ax1.plot(t, dw * 1e3, color="k", lw=1.2, label="δω 总计")
    ax1.plot(t, dw_slm * 1e3, color=C_A, lw=0.9, alpha=0.8, label="η·U_SLM/ħ")
    ax1.plot(t, dw_aod * 1e3, color=C_D, lw=0.9, alpha=0.8, label="η·U_AOD/ħ")
    ax1.set_ylabel("δω (mrad/s)")
    ax1.legend(fontsize=8)
    arrays += [dw, dw_slm, dw_aod]
    series += 3
    ax2 = axes[2]
    for mode, rec in modes.items():
        phi = np.asarray(rec["cumulative_phi"], dtype=float)
        ax2.plot(t, phi, color=DD_COLOR.get(mode, None), lw=1.2,
                 label=f"累计 φ（{mode}）")
        arrays.append(phi)
        series += 1
    ax2.set_ylabel("φ (rad)")
    ax2.set_xlabel("时间 (μs)")
    ax2.legend(fontsize=8)
    return _finish(fig, path, series, arrays)


# ------------------------------------------------------------ 多轮统计图
def plot_survival_vs_round(round_data, fits, path):
    """图6：absorbing 与 non-absorbing survival、Wilson 区间与可选拟合。"""
    fig, ax = plt.subplots(figsize=(10, 5.5))
    series = 0
    arrays = []
    for label, data in round_data.items():
        rounds = np.array(data["rounds"], dtype=float)
        absorb = np.array(data["absorbing"], dtype=float)
        nonabs = np.array(data["non_absorbing"], dtype=float)
        low = np.array(data["absorbing_low"], dtype=float)
        high = np.array(data["absorbing_high"], dtype=float)
        color = DD_COLOR.get(label.split("|")[-1], C_A) if "|" in label else C_A
        ax.fill_between(rounds, low, high, alpha=0.15, color=color)
        ax.plot(rounds, absorb, "o-", color=color,
                label=f"{label} absorbing")
        ax.plot(rounds, low, "-", color=color, alpha=0.4, lw=0.8)
        ax.plot(rounds, nonabs, "s--", color=color, alpha=0.55,
                label=f"{label} 非吸收")
        arrays += [absorb, nonabs, low, high]
        series += 2
        fit = fits.get(label)
        if fit and fit.get("fits"):
            best_name = min(fit["fits"], key=lambda k: fit["fits"][k]["aic"])
            params = fit["fits"][best_name]["params"]
            grid_n = np.linspace(rounds.min(), rounds.max(), 200)
            if best_name == "exponential":
                curve = params[0] ** grid_n
            elif best_name == "clipped_boltzmann":
                curve = params[1] + (1 - params[1]) * np.exp(-grid_n / params[0])
            else:
                curve = (params[2] + (1 - params[2]) * (params[0] ** grid_n)
                         * np.exp(-grid_n / params[1]))
            ax.plot(grid_n, curve, ":", color=color, alpha=0.8,
                    label=f"{label} {best_name}")
            arrays.append(curve)
            series += 1
    ax.set_xlabel("轮次")
    ax.set_ylabel("生存概率")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize=7)
    return _finish(fig, path, series, arrays)


def plot_long_tail_survival(long_tail, path):
    """新增：128-shot 长尾（至 40 轮）absorbing survival 与 Wilson 区间。"""
    rounds = np.asarray(long_tail["rounds"], dtype=float)
    surv = np.asarray(long_tail["survival"], dtype=float)
    n = int(long_tail.get("shots", np.nan))
    counts = np.round(surv * n).astype(int) if np.isfinite(n) else None
    fig, ax = plt.subplots(figsize=(10, 5))
    arrays = [surv]
    series = 1
    label = (f"absorbing survival（{long_tail.get('shots')} shots，"
             f"dt={long_tail.get('dt_us', 'n/a')} μs）")
    if counts is not None:
        low = np.empty(len(rounds))
        high = np.empty(len(rounds))
        for i, k in enumerate(counts):
            lo, hi = wilson_interval(int(k), int(n))
            low[i], high[i] = lo, hi
        ax.fill_between(rounds, low, high, alpha=0.2, color=C_A,
                        label="Wilson 95% CI")
        ax.plot(rounds, low, "-", color=C_A, alpha=0.45, lw=0.9,
                label="Wilson 下界")
        arrays += [low, high]
        series += 2
    ax.plot(rounds, surv, "o-", color=C_A, ms=4, label=label)
    ax.set_xlabel("轮次")
    ax.set_ylabel("absorbing 生存概率")
    ax.set_ylim(-0.02, 1.05)
    ax.set_title(f"长尾重复往返（最多 {long_tail.get('max_rounds')} 轮，"
                 "Protocol A；仅在检查点轮次评估）", fontsize=10)
    ax.legend(fontsize=8)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_failure_events_by_round(failure_rows, protocols, path):
    """新增：重复往返中每轮新增 first-failure 事件的阶段堆叠。"""
    fig, axes = plt.subplots(len(protocols), 1,
                             figsize=(11, 4.2 * len(protocols)), squeeze=False)
    arrays = []
    series = 0
    for r, pname in enumerate(protocols):
        ax = axes[r][0]
        rows = [f for f in failure_rows if f["protocol"] == pname]
        max_round = max((int(f["round"]) for f in rows), default=0)
        bottoms = np.zeros(max_round)
        ax_series = 0
        for stage in STAGES:
            vals = np.array([sum(1 for f in rows
                                 if int(f["round"]) == rd
                                 and f["stage"] == stage)
                             for rd in range(1, max_round + 1)], dtype=float)
            if vals.sum() == 0:
                continue
            ax.bar(np.arange(1, max_round + 1), vals, bottom=bottoms,
                   label=stage)
            bottoms += vals
            arrays.append(vals)
            series += 1
            ax_series += 1
        ax.set_xlabel("轮次")
        ax.set_ylabel("本轮新增失败 shot 数")
        ax.set_title(f"{pname}：每轮新增 first-failure 事件（吸收语义，"
                     "每 shot 至多记一次）", fontsize=10)
        if ax_series:
            ax.legend(fontsize=7, ncol=4)
    return _finish(fig, path, series, arrays)


def plot_heating_vs_round(heating_data, path):
    """图7：仅成功/留阱 shots 的激发分布随轮次。"""
    fig, ax = plt.subplots(figsize=(10, 5.5))
    series = 0
    arrays = []
    for label, data in heating_data.items():
        rounds = np.array(data["rounds"], dtype=float)
        med = np.array(data["median"], dtype=float)
        q90 = np.array(data["q90"], dtype=float)
        q95 = np.array(data["q95"], dtype=float)
        ax.plot(rounds, med, "o-", label=f"{label} 中位")
        ax.fill_between(rounds, med - np.array(data["q10"], dtype=float),
                        q95, alpha=0.2)
        ax.plot(rounds, q90, "--", alpha=0.6, label=f"{label} 90%")
        arrays += [med, q90]
        series += 2
    ax.set_xlabel("轮次")
    ax.set_ylabel("SLM 激发 (μK)")
    ax.legend(fontsize=7)
    return _finish(fig, path, series, arrays)


def plot_phase_and_contrast(phase_data, contrast_data, path):
    """图8：相位分布与 none/echo/XY4 contrast 随轮次。"""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    series = 0
    arrays = []
    ax = axes[0]
    for label, phases in phase_data.items():
        phi = np.asarray(phases, dtype=float)
        ax.hist(phi, bins=40, alpha=0.5, color=DD_COLOR.get(label, None),
                label=f"{label}（C="
                      f"{abs(np.exp(1j * phi).mean()):.3f}）")
        arrays.append(phi)
        series += 1
    ax.set_xlabel("相位 φ (rad)")
    ax.set_ylabel("shot 数")
    ax.legend(fontsize=8)
    ax = axes[1]
    for label, data in contrast_data.items():
        rounds = np.array(data["rounds"], dtype=float)
        contrast = np.array(data["contrast"], dtype=float)
        ax.plot(rounds, contrast, "o-",
                color=DD_COLOR.get(label, None), label=label)
        arrays.append(contrast)
        series += 1
    ax.set_xlabel("轮次")
    ax.set_ylabel("条件相干对比度 C_cond")
    ax.set_ylim(0, 1.02)
    if contrast_data:
        ax.legend(fontsize=8)
    return _finish(fig, path, series, arrays,
                   y_nonzero_any_axis=bool(contrast_data))


def plot_usable_coherent_fraction(contrast_curve, survival, path):
    """新增：usable_coherent_fraction = P(roundtrip success)×C_cond 随轮次。"""
    protocols = sorted({k.split("|")[0] for k in contrast_curve})
    fig, axes = plt.subplots(1, len(protocols),
                             figsize=(7.0 * len(protocols), 4.8), squeeze=False)
    arrays = []
    series = 0
    for i, pname in enumerate(protocols):
        ax = axes[0][i]
        surv = survival.get(f"{pname}|none")
        if not surv:
            continue
        rounds = np.asarray(surv["rounds"], dtype=float)
        absorb = np.asarray(surv["absorbing"], dtype=float)
        for mode in ("none", "spin_echo", "xy4_per_long_move"):
            data = contrast_curve.get(f"{pname}|{mode}")
            if not data:
                continue
            contrast = np.asarray(data["contrast"], dtype=float)
            usable = contrast * absorb
            ax.plot(rounds, usable, "o-",
                    color=DD_COLOR.get(mode, None), label=mode)
            arrays.append(usable)
            series += 1
        ax.set_xlabel("轮次")
        ax.set_ylabel("usable coherent fraction")
        ax.set_ylim(0, 1.02)
        ax.set_title(pname, fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("usable_coherent_fraction = P(roundtrip success) × C_cond"
                 "（工程代理，非量子过程保真度）", fontsize=10)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


# ------------------------------------------------------------- 账本与消融
def plot_energy_ledger(ledger_summary, path):
    """图10：分段功分解、ΔE–W_ext 对照与残差（max 与 mean）。"""
    segs = [row["segment"] for row in ledger_summary]
    x = np.arange(len(segs))
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    w = 0.22
    series = 0
    arrays = []
    axes[0].bar(x - 1.5 * w, [r["mean_w_move_uK"] for r in ledger_summary],
                width=w, label="移动功", color=C_A)
    axes[0].bar(x - 0.5 * w, [r["mean_w_shift_uK"] for r in ledger_summary],
                width=w, label="lensing 偏移功", color=C_C)
    axes[0].bar(x + 0.5 * w, [r["mean_w_aod_depth_uK"] for r in ledger_summary],
                width=w, label="AOD 深度功", color=C_D)
    axes[0].bar(x + 1.5 * w, [r["mean_w_slm_depth_uK"] for r in ledger_summary],
                width=w, label="SLM 深度功", color=C_E)
    axes[0].plot(x, [r["mean_delta_e_uK"] for r in ledger_summary], "k_",
                 ms=12, label="ΔE")
    axes[0].set_ylabel("均值 (μK)")
    axes[0].set_title("(a) 各分段外部功四通道分解与 ΔE（均值）", fontsize=10)
    axes[0].legend(fontsize=8)
    arrays += [[r["mean_w_move_uK"] for r in ledger_summary],
               [r["mean_delta_e_uK"] for r in ledger_summary]]
    series += 5
    axes[1].bar(x - 0.18, [r["mean_delta_e_uK"] for r in ledger_summary],
                width=0.36, label="ΔE（能量变化）", color=C_A)
    axes[1].bar(x + 0.18, [r["mean_w_ext_uK"] for r in ledger_summary],
                width=0.36, label="W_ext（外部功）", color=C_C)
    axes[1].set_ylabel("均值 (μK)")
    axes[1].set_title("(b) 功-能平衡：ΔE 与 W_ext 应逐段一致", fontsize=10)
    axes[1].legend(fontsize=8)
    arrays += [[r["mean_delta_e_uK"] for r in ledger_summary],
               [r["mean_w_ext_uK"] for r in ledger_summary]]
    series += 2
    axes[2].bar(x - 0.18, [r["max_abs_residual_uK"] for r in ledger_summary],
                width=0.36, label="max|R_W|", color="0.55")
    if all("mean_abs_residual_uK" in r for r in ledger_summary):
        axes[2].bar(x + 0.18, [r["mean_abs_residual_uK"] for r in ledger_summary],
                    width=0.36, label="mean|R_W|", color=C_E)
        arrays.append([r["mean_abs_residual_uK"] for r in ledger_summary])
        series += 2
    else:
        arrays.append([r["max_abs_residual_uK"] for r in ledger_summary])
        series += 1
    axes[2].set_ylabel("|残差| (μK)")
    axes[2].set_yscale("log")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(segs, rotation=25, ha="right", fontsize=8)
    axes[2].set_title("(c) 分段残差 R_W = ΔE − W_ext（对数刻度；静态 hold 段"
                      "应只有数值误差）", fontsize=10)
    axes[2].legend(fontsize=8)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_ablation_comparison(rows, path):
    """图11：消融对比（成功率、激发、C(none)、C(xy4) 四联图）。"""
    names = [row["condition"] for row in rows]

    def col(key):
        return np.array([row.get(key) if row.get(key) is not None else np.nan
                         for row in rows], dtype=float)

    vals = col("success_fraction")
    lows = np.clip(vals - col("wilson_low"), 0.0, None)
    highs = np.clip(col("wilson_high") - vals, 0.0, None)
    lows = np.nan_to_num(lows)
    highs = np.nan_to_num(highs)
    exc = col("median_excitation_uK")
    c_none = col("contrast_none")
    c_xy4 = col("contrast_xy4")
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.5))
    x = np.arange(len(names))
    series = 0
    arrays = []
    axes[0][0].bar(x, np.nan_to_num(vals), yerr=[lows, highs], capsize=3,
                   color=C_A, label="roundtrip success（Wilson 95% CI）")
    axes[0][0].set_ylabel("roundtrip success 概率")
    axes[0][0].set_ylim(0, 1.05)
    axes[0][0].legend(fontsize=8)
    axes[0][0].set_title("(a) 协议成功率", fontsize=10)
    arrays.append(np.nan_to_num(vals))
    series += 1
    axes[0][1].bar(x, np.nan_to_num(exc), color=C_C,
                   label="末态激发中位数")
    axes[0][1].set_ylabel("SLM 激发中位数 (μK)")
    axes[0][1].legend(fontsize=8)
    axes[0][1].set_title("(b) 留阱激发（加热）", fontsize=10)
    arrays.append(np.nan_to_num(exc))
    series += 1
    axes[1][0].bar(x, np.nan_to_num(c_none), color=C_A,
                   label="C(none)")
    axes[1][0].set_ylim(0, 1.05)
    axes[1][0].set_ylabel("C_cond")
    axes[1][0].legend(fontsize=8)
    axes[1][0].set_title("(c) 无 DD 条件相干对比度", fontsize=10)
    arrays.append(np.nan_to_num(c_none))
    series += 1
    axes[1][1].bar(x, np.nan_to_num(c_xy4), color=C_D,
                   label="C(XY4)")
    axes[1][1].set_ylim(0, 1.05)
    axes[1][1].set_ylabel("C_cond")
    axes[1][1].legend(fontsize=8)
    axes[1][1].set_title("(d) XY4-per-long-move 条件相干对比度", fontsize=10)
    arrays.append(np.nan_to_num(c_xy4))
    series += 1
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xticks(x)
            ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_timestep_convergence(rows, path):
    """图12：状态、功、相位与 success 标签随步长收敛。"""
    all_rows = list(rows)
    diff_rows = [r for r in all_rows if r.get("max_pos_diff_um") is not None]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    series = 0
    arrays = []
    if diff_rows:
        dts = np.array([row["dt_us"] for row in diff_rows], dtype=float)
        axes[0].plot(dts, [r["max_pos_diff_um"] for r in diff_rows], "o-",
                     color=C_A, label="max|Δpos| (μm)")
        axes[0].plot(dts, [r["max_vel_diff"] for r in diff_rows], "s--",
                     color=C_C, label="max|Δvel| (m/s)")
        axes[0].set_yscale("log")
        axes[0].legend(fontsize=8)
        arrays += [dts]
        series += 2
    axes[0].set_xlabel("dt (μs)")
    axes[0].set_title("(a) 末态状态差（相对 dt=0.025）", fontsize=10)
    if diff_rows:
        dts = np.array([row["dt_us"] for row in diff_rows], dtype=float)
        axes[1].plot(dts, [max(r["max_work_diff_uK"], 1e-9) for r in diff_rows],
                     "o-", color=C_D, label="max|ΔW| (μK)")
        axes[1].set_yscale("log")
        axes[1].legend(fontsize=8)
        arrays.append(dts)
        series += 1
    axes[1].set_xlabel("dt (μs)")
    axes[1].set_title("(b) 外部功差", fontsize=10)
    if diff_rows:
        dts = np.array([row["dt_us"] for row in diff_rows], dtype=float)
        axes[2].plot(dts, [max(r["max_phase_diff_rad"], 1e-12) for r in diff_rows],
                     "o-", color=C_E, label="max|Δφ| (rad)")
        axes[2].set_yscale("log")
        axes[2].legend(fontsize=8)
        arrays.append(dts)
        series += 1
    axes[2].set_xlabel("dt (μs)")
    axes[2].set_title("(c) 相位差", fontsize=10)
    dts_all = np.array([r["dt_us"] for r in all_rows], dtype=float)
    ok_all = np.array([r.get("success_rate") if r.get("success_rate")
                       is not None else np.nan for r in all_rows], dtype=float)
    axes[3].plot(dts_all, ok_all, "o-", color=C_A, label="success 率")
    for row in diff_rows:
        if row.get("label_flips") is not None:
            axes[3].annotate(f"{row['label_flips']} 翻转",
                             (row["dt_us"], np.nan_to_num(
                                 row.get("success_rate") or 0.0) + 0.012),
                             fontsize=7, ha="center")
    axes[3].set_ylim(min(0.0, np.nanmin(ok_all) - 0.02), 1.05)
    axes[3].set_xlabel("dt (μs)")
    axes[3].set_ylabel("success 率")
    axes[3].legend(fontsize=8)
    axes[3].set_title("(d) success 率与标签翻转数", fontsize=10)
    arrays.append(np.nan_to_num(ok_all))
    series += 1
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_noise_diagnostic(traces, path):
    """新增：assumed 强度噪声示例（OU 轨迹、quasistatic 分布、SLM/AOD 相关）。"""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    arrays = []
    series = 0
    ou = traces.get("ou")
    if ou:
        t_us = np.asarray(ou["times_s"]) * 1e6
        eps_s = np.asarray(ou["eps_slm"])
        eps_a = np.asarray(ou["eps_aod"])
        rms = float(ou["rms"])
        k = min(6, eps_s.shape[1] if eps_s.ndim > 1 else 1)
        for j in range(k):
            axes[0].plot(t_us, eps_s[:, j], color=C_A, lw=0.8, alpha=0.7,
                         label="ε_SLM(t)" if j == 0 else None)
            axes[0].plot(t_us, eps_a[:, j], color=C_D, lw=0.8, alpha=0.7,
                         ls="--", label="ε_AOD(t)" if j == 0 else None)
        axes[0].axhline(rms, color="0.3", lw=0.9, ls=":")
        axes[0].axhline(-rms, color="0.3", lw=0.9, ls=":")
        axes[0].set_xlabel("时间 (μs)")
        axes[0].set_ylabel("分数深度偏移 ε")
        axes[0].set_title(f"(a) OU 轨迹示例（rms={rms:g}，τ="
                          f"{ou['tau_s'] * 1e6:.0f} μs，前 {k} shots）",
                          fontsize=9)
        axes[0].legend(fontsize=8)
        arrays += [eps_s[:, 0], eps_a[:, 0]]
        series += 2
    quasi = traces.get("quasi")
    if quasi:
        qs = np.asarray(quasi["eps_slm"])
        qa = np.asarray(quasi["eps_aod"])
        axes[1].hist(qs, bins=40, alpha=0.6, color=C_A,
                     label=f"ε_SLM（std={np.std(qs):.3f}）")
        axes[1].hist(qa, bins=40, alpha=0.6, color=C_D,
                     label=f"ε_AOD（std={np.std(qa):.3f}）")
        axes[1].set_xlabel("quasistatic 分数偏移 ε")
        axes[1].set_ylabel("shot 数")
        axes[1].set_title(f"(b) quasistatic 每 shot 偏移分布"
                          f"（rms={quasi['rms']:g}）", fontsize=9)
        axes[1].legend(fontsize=8)
        arrays += [qs, qa]
        series += 2
        axes[2].plot(qs, qa, "o", ms=2.5, alpha=0.5, color=C_E,
                     label=f"quasistatic（ρ={np.corrcoef(qs, qa)[0, 1]:+.2f}）")
        arrays += [qs, qa]
        series += 1
    if ou:
        f_s = np.asarray(ou["eps_slm"])[:, -1]
        f_a = np.asarray(ou["eps_aod"])[:, -1]
        axes[2].plot(f_s, f_a, "s", ms=2.5, alpha=0.5, color=C_B,
                     label=f"OU 末态（ρ={np.corrcoef(f_s, f_a)[0, 1]:+.2f}）")
        arrays += [f_s, f_a]
        series += 1
    axes[2].set_xlabel("ε_SLM")
    axes[2].set_ylabel("ε_AOD")
    axes[2].set_title("(c) SLM/AOD 噪声相关性（假设 cross_correlation）",
                      fontsize=9)
    axes[2].legend(fontsize=8)
    fig.suptitle("强度噪声模型示例（参数 assumed，仅敏感性模式启用；"
                 "种子与初态采样隔离）", fontsize=10)
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


def plot_sensitivity_panels(rows, path):
    """图13：温度、深度、对准、lensing、η、噪声与 pulse timing 敏感性。"""
    factors = sorted({row["factor"] for row in rows})
    n = len(factors)
    n_cols = 3
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.8 * n_cols, 3.8 * n_rows),
                             squeeze=False)
    arrays = []
    series = 0
    for i, factor in enumerate(factors):
        ax = axes[i // n_cols][i % n_cols]
        sub = sorted([r for r in rows if r["factor"] == factor],
                     key=lambda r: r["value"] if not isinstance(r["value"], str)
                     else str(r["value"]))
        xlabels = [str(r["value"]) for r in sub]
        x = np.arange(len(sub))
        vals = [r["success_fraction"] if r["success_fraction"] is not None
                else np.nan for r in sub]
        has_prob = any(r["success_fraction"] is not None for r in sub)
        if has_prob:
            errs = np.clip(
                [[(r["success_fraction"] - r["wilson_low"])
                  if r["success_fraction"] is not None else 0.0 for r in sub],
                 [(r["wilson_high"] - r["success_fraction"])
                  if r["success_fraction"] is not None else 0.0 for r in sub]],
                0.0, None)
            ax.bar(x, np.nan_to_num(vals), yerr=errs, capsize=3, color=C_A,
                   alpha=0.85)
        else:
            ax.bar(x, np.zeros(len(sub)), color="none",
                   edgecolor="0.7", hatch="//")
        c_none = [r.get("contrast_none") for r in sub]
        c_xy4 = [r.get("contrast_xy4") for r in sub]
        if any(v is not None for v in c_none) or \
                any(v is not None for v in c_xy4):
            ax2 = ax.twinx()
            if any(v is not None for v in c_none):
                ax2.plot(x, [v if v is not None else np.nan for v in c_none],
                         "o--", color=C_D, label="C(none)")
                arrays.append(np.nan_to_num(np.array(
                    [v if v is not None else 0.0 for v in c_none],
                    dtype=float)))
                series += 1
            if any(v is not None for v in c_xy4):
                ax2.plot(x, [v if v is not None else np.nan for v in c_xy4],
                         "s-", color=DD_COLOR["xy4_per_long_move"],
                         label="C(XY4)")
                arrays.append(np.nan_to_num(np.array(
                    [v if v is not None else 0.0 for v in c_xy4],
                    dtype=float)))
                series += 1
            ax2.set_ylabel("C_cond", fontsize=8)
            ax2.set_ylim(0, 1.02)
            ax2.legend(fontsize=7, loc="lower left")
        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, fontsize=7, rotation=20)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("success 概率", fontsize=8)
        ax.set_title(factor, fontsize=10)
        arrays.append(np.nan_to_num(np.array(vals, dtype=float)))
    for j in range(n, n_rows * n_cols):
        fig.delaxes(axes[j // n_cols][j % n_cols])
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)
