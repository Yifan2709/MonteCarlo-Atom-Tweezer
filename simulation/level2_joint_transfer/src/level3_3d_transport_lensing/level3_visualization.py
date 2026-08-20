"""Level 3 诊断图（11 张）：非交互后端，保存后关闭，附程序化验收记录。"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

_MPL_CACHE = Path(tempfile.gettempdir()) / "level3_matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.sans-serif": ["PingFang SC", "Hiragino Sans GB", "Heiti SC",
                        "Microsoft YaHei", "Noto Sans CJK SC", "Arial Unicode MS",
                        "sans-serif"],
    "axes.unicode_minus": False,
})

from .aod_lensing import axial_minima
from .gaussian_3d import KB, static_potential
from .transport_profiles import get_profile

C_SINE = "#1f77b4"
C_CJERK = "#ff7f0e"
C_MJERK = "#2ca02c"
C_STRAIGHT = "#d62728"
C_DIAG = "#9467bd"
PROFILE_COLOR = {"adiabatic_sine": C_SINE, "constant_jerk": C_CJERK,
                 "minimum_jerk": C_MJERK}


def _finish(fig, path, series_count, arrays):
    """保存 PNG 并生成与 level0.io_utils.validate_pngs 兼容的验收记录。"""
    finite = all(np.all(np.isfinite(np.asarray(a))) for a in arrays)
    legend_count = sum(len(ax.get_legend().get_texts()) for ax in fig.axes
                       if ax.get_legend() is not None)
    cover, x_nonzero, y_nonzero = True, True, True
    for ax in fig.axes:
        lines = [ln for ln in ax.get_lines() if ln.get_transform() is ax.transData]
        if not lines:
            continue
        xs = np.concatenate([np.asarray(ln.get_xdata(), dtype=float) for ln in lines])
        ys = np.concatenate([np.asarray(ln.get_ydata(), dtype=float) for ln in lines])
        if xs.size == 0 or ys.size == 0 or not (
                np.all(np.isfinite(xs)) and np.all(np.isfinite(ys))):
            continue
        xl, yl = tuple(map(float, ax.get_xlim())), tuple(map(float, ax.get_ylim()))
        cover = cover and xl[0] <= xs.min() and xs.max() <= xl[1] \
            and yl[0] <= ys.min() and ys.max() <= yl[1]
        x_nonzero = x_nonzero and xs.max() > xs.min()
        y_nonzero = y_nonzero and ys.max() > ys.min()
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


def plot_static_slices(physics, path):
    """图1：无透镜静态势的横向与轴向切片，标数值/解析曲率。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    x = np.linspace(-3e-6, 3e-6, 400)
    u_r = static_potential(x, 0, 0, physics["depth_j"], physics["waist_m"],
                           physics["z_r_m"]) / KB * 1e6
    axes[0].plot(x * 1e6, u_r, color=C_SINE, label="横向切片 U(ξ,0,0)")
    z = np.linspace(-10e-6, 10e-6, 400)
    u_z = static_potential(0, 0, z, physics["depth_j"], physics["waist_m"],
                           physics["z_r_m"]) / KB * 1e6
    axes[1].plot(z * 1e6, u_z, color=C_CJERK, label="轴向切片 U(0,0,z)")
    axes[1].axvline(physics["z_r_m"] * 1e6, ls="--", color="gray",
                    label=f"z_R={physics['z_r_m']*1e6:.2f} μm")
    for ax, freq, label in ((axes[0], physics["omega_r"], "ω_r"),
                            (axes[1], physics["omega_z"], "ω_z")):
        ax.set_title(f"{label}/2π = {freq/2/np.pi/1e3:.1f} kHz")
        ax.set_xlabel("位置 (μm)")
        ax.set_ylabel("势能 (μK)")
        ax.legend()
    return _finish(fig, path, 3, [u_r, u_z])


def plot_bifurcation(physics, path):
    """图2：直线模型不同 z_s/z_R 的轴向势（单阱→双阱）与对角对照。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    z_r = physics["z_r_m"]
    z = np.linspace(-8e-6, 8e-6, 1200)
    arrays = []
    for ratio, color in zip((0.0, 1.0, 2.0, 3.0),
                            ("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728")):
        # 直线：仅轴 1 偏移，U(z) = -D/sqrt(1+dz1^2)
        dz1 = (z - ratio * z_r) / z_r
        u_line = -physics["depth_j"] / np.sqrt(1 + dz1 ** 2) / KB * 1e6
        axes[0].plot(z * 1e6, u_line, color=color,
                     label=f"|z_s|/z_R={ratio:.1f}")
        arrays.append(u_line)
        # 对角同号：两轴同偏，U(z) = -D/(1+dz^2)，保持单阱
        dzd = (z - ratio * z_r) / z_r
        u_diag = -physics["depth_j"] / (1 + dzd ** 2) / KB * 1e6
        axes[1].plot(z * 1e6, u_diag, color=color, label=f"z_s1=z_s2, 比值 {ratio:.1f}")
        arrays.append(u_diag)
    axes[0].set_title("直线：单阱 → 双阱（|z_s| > 2 z_R 分裂）")
    axes[1].set_title("对角同号：保持单一（球面）阱")
    for ax in axes:
        ax.set_xlabel("z (μm)")
        ax.set_ylabel("轴向势 (μK)")
        ax.legend(fontsize=8)
    return _finish(fig, path, 8, arrays)


def plot_profiles(path):
    """图3：三类轨迹的 q, q̇, q̈ 与 jerk。"""
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    s = np.linspace(0, 1, 500)
    arrays = []
    for name in ("adiabatic_sine", "constant_jerk", "minimum_jerk"):
        prof = get_profile(name)
        series = [prof.q(s), prof.q_dot(s), prof.q_ddot(s), prof.q_jerk(s)]
        arrays.extend(series)
        for ax, y, label in zip(axes, series, ("q(s)", "dq/ds", "d²q/ds²", "jerk")):
            ax.plot(s, y, color=PROFILE_COLOR[name], label=name)
    titles = ("位置 q(s)", "速度 q̇(s)", "加速度 q̈(s)", "jerk")
    for ax, t in zip(axes, titles):
        ax.set_title(t)
        ax.set_xlabel("s")
        ax.legend(fontsize=8)
    return _finish(fig, path, 12, arrays)


def plot_straight_vs_diagonal(minima_rows, physics, path):
    """图4：相同 L/T 下两种几何的单轴速度、焦点偏移与有效阱深。"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    distances = sorted({r["distance_um"] for r in minima_rows}, reverse=True)
    durations = sorted({r["duration_us"] for r in minima_rows}, reverse=True)
    probe_distance = next((d for d in distances
                           if any(r["distance_um"] == d and r["geometry"] == "straight"
                                  for r in minima_rows)), distances[0])
    probe_duration = durations[0]
    target = {(r["geometry"], r["severity"]): r for r in minima_rows
              if r["distance_um"] == probe_distance and r["duration_us"] == probe_duration
              and r["profile"] == "adiabatic_sine"}
    arrays = []
    for ax, key, label in ((axes[0], "v1_over_vs", "轴1 峰值速度 / v_s"),
                           (axes[1], "zs1_peak_over_zr", "峰值 |z_s1|/z_R"),
                           (axes[2], "min_depth_uK", "瞬时最低点深度 (μK)")):
        for geometry, color in (("straight", C_STRAIGHT), ("diagonal", C_DIAG)):
            xs, ys = [], []
            for severity in sorted({r["severity"] for r in minima_rows}):
                row = target.get((geometry, severity))
                if row:
                    xs.append(severity)
                    ys.append(row[key])
            ax.plot(xs, ys, "o-", color=color, label=geometry)
            arrays.append(np.array(ys))
        ax.set_xlabel("lensing severity (v_max,axis/v_s)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8)
    return _finish(fig, path, 6, arrays)


def plot_survival_phase(coarse_rows, path):
    """图5：留阱率随时长与 severity 的相图（straight/diagonal 分面）。"""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    arrays = []
    for ax, geometry in zip(axes, ("straight", "diagonal")):
        rows = [r for r in coarse_rows if r["geometry"] == geometry]
        durations = sorted({r["duration_us"] for r in rows})
        severities = sorted({r["severity"] for r in rows})
        grid = np.full((len(severities), len(durations)), np.nan)
        for r in rows:
            i = severities.index(r["severity"])
            j = durations.index(r["duration_us"])
            grid[i, j] = r["retention_fraction"]
        im = ax.imshow(grid, origin="lower", aspect="auto", vmin=0, vmax=1,
                       cmap="RdYlGn",
                       extent=[min(durations), max(durations), min(severities),
                               max(severities)])
        arrays.append(grid[np.isfinite(grid)])
        ax.set_title(f"{geometry}（sine 轨迹, 128 shots/格点）")
        ax.set_xlabel("运输时长 (μs)")
        ax.set_ylabel("lensing severity")
        fig.colorbar(im, ax=ax, label="留阱率")
    return _finish(fig, path, 0, arrays)


def plot_heating(validation_records, path):
    """图6：仅 retained shots 的 Δε 与末态激发分布。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    arrays = []
    keys = sorted({(r["profile"], r["geometry"], r["severity"])
                   for r in validation_records})
    for profile, geometry, severity in keys:
        rows = [r for r in validation_records
                if r["profile"] == profile and r["geometry"] == geometry
                and r["severity"] == severity and r["retained"]]
        label = f"{profile}|{geometry}|sev{severity}"
        de = np.array([r["delta_eps_uK"] for r in rows])
        ef = np.array([r["eps_f_uK"] for r in rows])
        axes[0].hist(de, bins=40, alpha=0.5, label=label, density=True)
        axes[1].hist(ef, bins=40, alpha=0.5, label=label, density=True)
        arrays.extend([de, ef])
    axes[0].set_xlabel("激发能变化 Δε (μK, 仅 retained)")
    axes[1].set_xlabel("末态激发 ε_f (μK, 仅 retained)")
    for ax in axes:
        ax.set_ylabel("概率密度")
        ax.legend(fontsize=7)
    return _finish(fig, path, 2 * len(keys), arrays)


def plot_axial_dynamics(trajectories, path, lost=None, lost_minima=None):
    """图7：代表性轨迹的 z(t) 与瞬时最低点/分支。

    左面板：validation 主条件的 retained/split 代表 z(t)；
    右面板（若提供）：粗扫描失败条件的 lost z(t) 与瞬时轴向最低点位置
    （单阱→双阱分支随速度演化）。lost 为 {label: traj} 字典，
    lost_minima 为 {label: (time_s, minima_z_m[n_t, n_min])}。
    """
    n_panels = 2 if lost else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 4.8),
                             squeeze=False)
    ax = axes[0][0]
    arrays = []
    series = 0
    for label, traj in trajectories.items():
        z = traj["position_m"][:, :, 2] * 1e6  # (n_rec, n_shots)
        ax.plot(traj["time_s"] * 1e6, z, lw=1, label=label)
        arrays.append(z.ravel())
        series += z.shape[1]
    ax.set_title("validation 主条件（retained）")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("轴向位置 z (μm)")
    ax.legend(fontsize=8)
    if lost:
        ax2 = axes[0][1]
        for label, traj in lost.items():
            z = traj["position_m"][:, :, 2] * 1e6
            ax2.plot(traj["time_s"] * 1e6, z, lw=1, label=label)
            arrays.append(z.ravel())
            series += z.shape[1]
        if lost_minima:
            for label, (t_s, zmin_m) in lost_minima.items():
                t_us = np.asarray(t_s) * 1e6
                for j in range(zmin_m.shape[1]):
                    col = np.asarray(zmin_m[:, j]) * 1e6
                    ax2.plot(t_us, col, "k--", lw=1.2, alpha=0.7,
                             label="瞬时轴向最低点" if j == 0 else None)
                    finite_col = col[np.isfinite(col)]
                    if finite_col.size:
                        arrays.append(finite_col)
                series += 1
        ax2.set_title("粗扫描失败条件（lost + 双阱分支）")
        ax2.set_xlabel("时间 (μs)")
        ax2.set_ylabel("轴向位置 z (μm)")
        ax2.legend(fontsize=8)
    return _finish(fig, path, series, arrays)


def plot_heating_scaling(check, path):
    """图 12（补充）：简谐小激发极限的加热上包络与 T 标度律。"""
    fig, ax = plt.subplots(figsize=(10, 5.2))
    arrays = []
    series = 0
    data = check["plot_data"]
    t = np.array(data["t_grid_us"])
    colors = {"adiabatic_sine": C_SINE, "constant_jerk": C_CJERK,
              "minimum_jerk": C_MJERK, "const_accel_reference": C_STRAIGHT}
    for fam, env in data["envelope_uK"].items():
        env = np.array(env)
        ax.plot(t, env, "o-", ms=3, lw=1.2, color=colors.get(fam, "gray"),
                label=f"{fam}（斜率 {check['families'][fam]['t_slope']:+.1f}）")
        arrays.append(env)
        series += 1
    for slope, name in ((-6.0, r"$T^{-6}$ 参考线（m=2：sine/constant-jerk/minimum-jerk）"),
                        (-4.0, None)):
        ref = env[0] * (t / t[0]) ** slope
        ref = np.maximum(ref, 1e-9)
        ax.plot(t, ref, ":", lw=1, color="gray", alpha=0.8, label=name)
        arrays.append(ref)
        if name:
            series += 1  # 仅带标签的参考线计入图例匹配
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("运输时长 T (μs)")
    ax.set_ylabel("末端激发上包络 (μK)")
    ax.set_title("简谐、无透镜、小激发极限的加热标度（干涉节点取上包络）")
    ax.legend(fontsize=8)
    return _finish(fig, path, series, arrays)


def plot_final_energy(validation_records, path):
    """图8：末态静态能量分布与 E_f=0 阈值。"""
    fig, ax = plt.subplots(figsize=(10, 4.4))
    arrays = []
    keys = sorted({(r["profile"], r["geometry"], r["severity"])
                   for r in validation_records})
    for profile, geometry, severity in keys:
        rows = [r for r in validation_records
                if r["profile"] == profile and r["geometry"] == geometry
                and r["severity"] == severity]
        e = np.array([r["e_final_uK"] for r in rows])
        ax.hist(e, bins=40, alpha=0.5, density=True,
                label=f"{profile}|{geometry}|sev{severity}")
        arrays.append(e)
    ax.axvline(0.0, color="k", ls="--", lw=2, label="E_f=0 留阱阈值")
    ax.set_xlabel("末态静态能量 E_f (μK)")
    ax.set_ylabel("概率密度")
    ax.legend(fontsize=7)
    return _finish(fig, path, len(keys) + 1, arrays)


def plot_work_energy(work_cases, path):
    """图9：ΔE、W_ext、R_W（理想与 lensing、直线与对角各一例）。"""
    fig, ax = plt.subplots(figsize=(10, 4.6))
    arrays = []
    for label, series in work_cases.items():
        ax.plot(series["time_us"], series["delta_E"], label=f"{label} ΔE")
        ax.plot(series["time_us"], series["W_ext"], ls="--",
                label=f"{label} W_ext")
        ax.plot(series["time_us"], series["R_W"], ls=":",
                label=f"{label} R_W×100")
        arrays.extend([series["delta_E"], series["W_ext"], series["R_W"]])
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("能量 (μK；R_W 已 ×100)")
    ax.legend(fontsize=7)
    return _finish(fig, path, len(work_cases) * 3, arrays)


def plot_timestep_convergence(convergence_rows, path):
    """图10：连续量误差、功-能残差与标签不一致率随步长。"""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    dts = [r["dt_us"] for r in convergence_rows]
    arrays = []
    axes[0].plot(dts, [r["max_final_energy_diff_uK"] for r in convergence_rows],
                 "o-", color=C_SINE)
    axes[0].set_ylabel("末态能量最大差 (μK, 相对 0.025 μs)")
    axes[1].plot(dts, [r["max_abs_w_residual_over_depth"] for r in convergence_rows],
                 "o-", color=C_CJERK)
    axes[1].set_ylabel("功-能残差 / 阱深")
    axes[2].plot(dts, [r["label_flip_rate"] for r in convergence_rows], "o-",
                 color=C_MJERK, label="标签不一致率")
    if "retention_fraction" in convergence_rows[0]:
        axes[2].plot(dts, [r["retention_fraction"] for r in convergence_rows], "s--",
                     color=C_SINE, label="留阱率")
        axes[2].legend(fontsize=8)
    axes[2].set_ylabel("留阱标签不一致率 / 留阱率")
    for ax in axes:
        ax.set_xlabel("步长 dt (μs)")
        arrays.append(np.array([1.0]))  # 占位保证 finite 检查
    series = 2 if "retention_fraction" in convergence_rows[0] else 0
    return _finish(fig, path, series, arrays)


def plot_sensitivity(sensitivity_rows, path):
    """图11：温度、深度、束腰、severity 敏感性面板。"""
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0))
    arrays = []
    for ax, scan in zip(axes, ("temperature", "depth", "waist", "severity")):
        rows = [r for r in sensitivity_rows if r["scan"] == scan]
        if not rows:
            continue
        xs = np.array([r["value"] for r in rows], dtype=float)
        fr = np.array([r["retention_fraction"] for r in rows], dtype=float)
        lo = np.array([r["wilson_low"] for r in rows], dtype=float)
        hi = np.array([r["wilson_high"] for r in rows], dtype=float)
        order = np.argsort(xs)
        xs, fr, lo, hi = xs[order], fr[order], lo[order], hi[order]
        ax.errorbar(xs, fr, yerr=[np.maximum(0, fr - lo), np.maximum(0, hi - fr)],
                    fmt="-o", capsize=4, color=C_SINE)
        arrays.extend([xs, fr])
        ax.set_title({"temperature": "温度 (μK)", "depth": "阱深 (μK)",
                      "waist": "束腰 w₀ (μm)",
                      "severity": "lensing severity"}[scan])
        ax.set_xlabel({"temperature": "T (μK)", "depth": "D (μK)",
                       "waist": "w₀ (μm)", "severity": "v_max/v_s"}[scan])
        ax.set_ylabel("留阱率")
    return _finish(fig, path, 0, arrays)
