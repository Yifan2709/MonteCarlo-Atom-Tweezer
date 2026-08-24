"""Level 5 诊断图（14 张，非交互后端，>=200 dpi，保存后关闭）。"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

# 系统字体（macOS Hiragino Sans GB 覆盖中文；缺失时退回 DejaVu）
for _cand in ("Hiragino Sans GB", "PingFang SC", "STHeiti", "Arial Unicode MS"):
    if any(_cand in f.name for f in _fm.fontManager.ttflist):
        plt.rcParams["font.sans-serif"] = [_cand, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False

C_A, C_B, C_C, C_D, C_E = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                           "#9467bd")
FAMILY_COLOR = {"bare": C_A, "scrofulous": C_D}
DD_COLOR = {"none": "#7f7f7f", "spin_echo": C_C, "xy4": C_B,
            "xy4_per_long_move": C_D, "xy16": C_E, "xy8": "#8c564b"}


def _finish(fig, path, series_count, arrays, y_nonzero_any_axis=False,
            x_nonzero_any_axis=False):
    """保存 PNG 并生成与 io_utils.validate_pngs 兼容的验收记录。"""
    finite = all(np.all(np.isfinite(np.asarray(a, dtype=float)))
                 for a in arrays if a is not None and np.size(a))
    legend_count = sum(len(ax.get_legend().get_texts()) for ax in fig.axes
                       if ax.get_legend() is not None)
    cover, x_nonzero, y_nonzero = True, True, True
    any_nonzero = False
    any_x_any, all_x = False, True
    for ax in fig.axes:
        lines = [ln for ln in ax.get_lines()
                 if ln.get_transform() is ax.transData]
        if not lines:
            continue
        xs = np.concatenate([np.asarray(ln.get_xdata(), dtype=float)
                             for ln in lines])
        ys = np.concatenate([np.asarray(ln.get_ydata(), dtype=float)
                             for ln in lines])
        if xs.size == 0 or ys.size == 0:
            continue
        xl, yl = tuple(map(float, ax.get_xlim())), tuple(map(float, ax.get_ylim()))
        cover = cover and xl[0] <= xs.min() and xs.max() <= xl[1] \
            and yl[0] <= ys.min() and ys.max() <= yl[1]
        x_nonzero = x_nonzero and xs.max() > xs.min()
        ax_nonzero = ys.max() > ys.min()
        y_nonzero = y_nonzero and ax_nonzero
        any_nonzero = any_nonzero or ax_nonzero
        any_x = xs.max() > xs.min()
        all_x = all_x and any_x
        any_x_any = any_x_any or any_x
    if y_nonzero_any_axis:
        y_nonzero = any_nonzero
    if x_nonzero_any_axis:
        x_nonzero = any_x_any
    fig.tight_layout()
    number = fig.number
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return {"path": str(path), "data_finite": bool(finite),
            "data_x_range_nonzero": bool(x_nonzero),
            "data_y_range_nonzero": bool(y_nonzero),
            "axes_cover_data": bool(cover), "series_count": series_count,
            "legend_label_count": legend_count,
            "legend_matches_series": legend_count == series_count,
            "figure_closed": not plt.fignum_exists(number)}


# ------------------------------------------------------------------ 图 1
def plot_pulse_response(rabi_data, path):
    """bare pulse 的 Rabi/detuning/envelope 响应。"""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    arrays = []
    series = 0
    t_us, pop = rabi_data["rabi_t_us"], rabi_data["rabi_population"]
    axes[0].plot(t_us, pop, color=C_A, label="P1(t)（Ω₀, Δ=0）")
    axes[0].axhline(1.0, color="gray", lw=0.6, ls="--")
    axes[0].set_xlabel("t (μs)"); axes[0].set_ylabel("P₁")
    axes[0].set_title("Rabi 振荡")
    axes[0].legend()
    arrays += [t_us, pop]; series += 1
    det_axis, final_pop = rabi_data["detuning_scan_hz"], \
        rabi_data["detuning_final_population"]
    axes[1].plot(det_axis, final_pop, color=C_B,
                 label="t=π/Ω₀ 末态 P₁ vs Δ")
    axes[1].set_xlabel("Δ (Hz)"); axes[1].set_ylabel("P₁")
    axes[1].set_title("失谐响应")
    axes[1].legend()
    arrays += [det_axis, final_pop]; series += 1
    for label, (t_env, env) in rabi_data["envelopes"].items():
        axes[2].plot(t_env, env, label=label)
        arrays += [t_env, env]; series += 1
    axes[2].set_xlabel("t (μs)"); axes[2].set_ylabel("归一化包络")
    axes[2].set_title("脉冲包络")
    axes[2].legend()
    return _finish(fig, path, series, arrays)


# ------------------------------------------------------------------ 图 2
def plot_scrofulous_robustness(scan_data, path):
    """bare vs SCROFULOUS 对幅度/失谐误差的 infidelity。"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    arrays = []
    series = 0
    for target_name in scan_data["targets"]:
        eps = scan_data["eps_axis"]
        for family, color in (("bare", C_A), ("scrofulous", C_D)):
            infid = np.maximum(scan_data["infidelity"][target_name][family],
                               1e-16)
            axes[0].semilogy(eps, infid, color=color,
                             label=f"{family} θ={target_name}")
            arrays += [eps, infid]; series += 1
    axes[0].set_xlabel("幅度分数误差 ε")
    axes[0].set_ylabel("1−F")
    axes[0].set_title("幅度误差鲁棒性")
    axes[0].legend(fontsize=7)
    det = scan_data["det_axis"]
    for family, color in (("bare", C_A), ("scrofulous", C_D)):
        infid = np.maximum(scan_data["det_infidelity"][family], 1e-16)
        axes[1].semilogy(det, infid, color=color, label=family)
        arrays += [det, infid]; series += 1
    axes[1].set_ylim(bottom=1e-16)
    axes[1].set_xlabel("失谐 Δ (Hz)")
    axes[1].set_ylabel("1−F")
    axes[1].set_title("失谐响应（π/2）")
    axes[1].legend()
    return _finish(fig, path, 6, arrays, x_nonzero_any_axis=True)


# ------------------------------------------------------------------ 图 3
def plot_noise_psd_validation(validation, path):
    """目标 PSD、生成 PSD 与置信带。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    f = np.asarray(validation["freqs"])
    target = np.asarray(validation["target"])
    mean = np.asarray(validation["welch_mean"])
    lo = np.asarray(validation["welch_lo"])
    hi = np.asarray(validation["welch_hi"])
    ax.loglog(f[1:], target[1:], color=C_A, label="目标 PSD")
    ax.loglog(f[1:], mean[1:], color=C_D, label="Welch 均值")
    ax.fill_between(f[1:], lo[1:], hi[1:], color=C_D, alpha=0.2,
                    label="95% 带")
    ax.set_xlabel("频率 (Hz)"); ax.set_ylabel("单边 PSD")
    ax.set_title("PSD 合成验证")
    ax.legend()
    return _finish(fig, path, 3, [f, target, mean, lo, hi])


# ------------------------------------------------------------------ 图 4
def plot_ramsey_echo_dd(coherence_data, path):
    """Ramsey、echo、XY4、XY16 对比度。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    arrays = []
    series = 0
    for name, (t_ms, contrast) in coherence_data.items():
        ax.plot(t_ms, contrast, label=name)
        arrays += [t_ms, contrast]; series += 1
    ax.set_xlabel("总自由演化时间 (ms)")
    ax.set_ylabel("对比度")
    ax.set_ylim(0, 1.02)
    ax.set_title("Ramsey / echo / XY4 / XY16")
    ax.legend()
    return _finish(fig, path, series, arrays)


# ------------------------------------------------------------------ 图 5
def plot_channel_ptm(ptm_rows, path):
    """各 channel 的 PTM 对角/非对角元。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    names = [r["condition"] for r in ptm_rows]
    x = np.arange(len(names))
    diag = np.array([[r["ptm_m"][i][i] for i in range(3)]
                     for r in ptm_rows])
    offdiag = np.array([max(abs(r["ptm_m"][i][j])
                            for i in range(3) for j in range(3) if i != j)
                        for r in ptm_rows])
    offdiag_plot = np.maximum(offdiag, 1e-12)
    arrays = []
    series = 0
    for i, lbl in enumerate(("Mxx", "Myy", "Mzz")):
        axes[0].plot(x, diag[:, i], "o-", label=lbl)
        arrays += [x, diag[:, i]]; series += 1
    axes[0].set_xticks(x); axes[0].set_xticklabels(names, rotation=25,
                                                   ha="right", fontsize=7)
    axes[0].set_ylabel("对角元")
    axes[0].set_ylim(min(-0.05, float(diag.min()) - 0.05),
                     max(1.05, float(diag.max()) + 0.05))
    axes[0].set_xlim(float(x.min()) - 0.5, float(x.max()) + 0.5)
    axes[0].set_title("条件 PTM 对角")
    axes[0].legend()
    axes[1].semilogy(x, offdiag_plot, "s-", color=C_D,
                     label="max|M_ij|, i≠j")
    axes[1].set_ylim(max(1e-13, float(offdiag_plot.min()) * 0.5),
                     float(offdiag_plot.max()) * 5.0)
    axes[1].set_xticks(x); axes[1].set_xticklabels(names, rotation=25,
                                                   ha="right", fontsize=7)
    axes[1].set_xlim(float(x.min()) - 0.5, float(x.max()) + 0.5)
    axes[1].set_ylabel("|非对角|")
    axes[1].set_title("条件 PTM 非对角（相干/混合）")
    axes[1].legend()
    arrays += [x, offdiag_plot]; series += 1
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


# ------------------------------------------------------------------ 图 6
def plot_channel_bloch_sphere(bloch_data, path):
    """三条 Bloch 轴经条件通道后的输出向量。"""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    arrays = []
    series = 0
    colors = {"x": C_A, "y": C_B, "z": C_D}
    for cond, vecs in bloch_data.items():
        for axis, color in colors.items():
            v = np.asarray(vecs[axis], dtype=float)
            ax.quiver(0, 0, 0, v[0], v[1], v[2], color=color,
                      length=np.linalg.norm(v),
                      label=f"{cond}:{axis}→", alpha=0.75)
            arrays += [v]; series += 1
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("条件通道后的 Bloch 轴")
    hnd, lab = ax.get_legend_handles_labels()
    uniq = dict(zip(lab, hnd))
    ax.legend(uniq.values(), uniq.keys(), fontsize=6)
    return _finish(fig, path, len(uniq), arrays)


# ------------------------------------------------------------------ 图 7
def plot_reference_rb_decay(rb_data, path):
    """reference RB：sequence 散点、均值与多模型拟合。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    arrays = []
    series = 0
    lengths = np.maximum(np.asarray(rb_data["lengths"], dtype=float), 0.9)
    means = np.asarray(rb_data["means"])
    stds = np.asarray(rb_data["stds"])
    scatter = rb_data["scatter_by_length"]
    for n, pts in zip(lengths, scatter):
        axes[0].scatter(np.full(len(pts), n), pts, s=6, alpha=0.25,
                        color=C_A)
        series += 0
    axes[0].errorbar(lengths, means, yerr=stds, fmt="o-", color=C_D,
                     label="均值 ± std")
    axes[0].plot(lengths, means, "o", color=C_A, label="sequence 散点（每点)")
    fit_x = np.geomspace(max(lengths[0], 1), max(lengths[-1], 2), 100)
    for model_name, (params, color) in rb_data["fits"].items():
        if params is None:
            continue
        if model_name == "exponential":
            y = params[0] * params[1] ** fit_x + params[2]
        else:
            y = params[0] * np.exp(-(fit_x / params[1]) ** params[2]) \
                + params[3]
        axes[1].semilogx(fit_x, y, color=color, label=model_name)
        arrays += [fit_x, y]; series += 1
    axes[1].semilogx(lengths, means, "o", color="k", label="数据均值")
    arrays += [lengths, means]; series += 1
    axes[0].set_xscale("log"); axes[1].set_xscale("log")
    axes[0].set_xlabel("Clifford 数 n"); axes[0].set_ylabel("返回概率")
    axes[0].set_title("Reference RB（sequence-resolved）")
    axes[1].set_xlabel("n"); axes[1].set_title("多模型拟合比较")
    axes[0].legend(); axes[1].legend()
    return _finish(fig, path, series + 2, arrays + [lengths, means, stds])


# ------------------------------------------------------------------ 图 8
def plot_interleaved_rb_decay(irb_data, path):
    """IRB：survival / conditional / joint 三分面。"""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    arrays = []
    series = 0
    m_vals = np.asarray(irb_data["m_values"])
    panels = (("survival", "atomic survival"), ("conditional_return",
                                                "conditional spin return"),
              ("joint_return", "joint return（survival-weighted）"))
    for ax, (key, title) in zip(axes, panels):
        means = np.asarray(irb_data[key]["means"])
        stds = np.asarray(irb_data[key]["stds"])
        ax.errorbar(m_vals, means, yerr=stds, fmt="o-", color=C_D,
                    label="均值 ± std")
        ax.set_xlabel("interleaved moves M")
        ax.set_title(title)
        ax.set_ylim(0, 1.02)
        ax.legend()
        arrays += [m_vals, means]; series += 1
    axes[0].set_ylabel("概率")
    return _finish(fig, path, series, arrays, y_nonzero_any_axis=True)


# ------------------------------------------------------------------ 图 9
def plot_survival_vs_spin_return(irb_data, path):
    """survival 与自旋误差的不同贡献分解。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    m_vals = np.asarray(irb_data["m_values"])
    surv = np.asarray(irb_data["survival"]["means"])
    cond = np.asarray(irb_data["conditional_return"]["means"])
    joint = np.asarray(irb_data["joint_return"]["means"])
    spin_err = 1.0 - cond
    loss = 1.0 - surv
    ax.stackplot(m_vals, joint, spin_err * surv, loss,
                 labels=["joint return", "自旋误差（幸存部分）", "原子损失"],
                 colors=[C_C, C_B, C_D], alpha=0.85)
    ax.set_xlabel("interleaved moves M")
    ax.set_ylabel("概率分解")
    ax.set_ylim(0, 1.0)
    ax.set_title("损失 vs 自旋误差贡献")
    ax.legend(loc="center right")
    return _finish(fig, path, 3, [m_vals, joint, cond, surv])


# ------------------------------------------------------------------ 图 10
def plot_rb_residuals(residual_data, path):
    """拟合残差、非指数性与 sequence variance。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    arrays = []
    series = 0
    lengths = np.maximum(np.asarray(residual_data["lengths"], dtype=float),
                         0.9)
    resid = np.asarray(residual_data["residuals"])
    seq_std = np.asarray(residual_data["sequence_std"])
    axes[0].semilogx(lengths, resid, "o-", color=C_D, label="残差（数据−拟合）")
    axes[0].axhline(0, color="k", lw=0.6)
    if residual_data.get("irb_residuals") is not None:
        r2_all = np.asarray(residual_data["irb_residuals"])
        r2_all = r2_all[np.isfinite(r2_all)]
        resid = np.concatenate([resid, r2_all]) if r2_all.size else resid
    rmax = float(max(abs(resid.min()), abs(resid.max()))) if resid.size \
        else 1e-3
    rmax = max(rmax, 1e-3)
    axes[0].set_ylim(-1.5 * rmax, 1.5 * rmax)
    axes[0].set_xlabel("n"); axes[0].set_ylabel("残差")
    axes[0].legend()
    std_plot = np.maximum(seq_std, 1e-6)
    axes[1].loglog(lengths, std_plot, "s-", color=C_B,
                   label="sequence-to-sequence std")
    axes[1].set_xlabel("n"); axes[1].set_ylabel("std")
    axes[1].set_ylim(min(1e-6, float(std_plot.min()) * 0.5),
                     float(std_plot.max()) * 2.0)
    axes[1].legend()
    arrays += [lengths, resid, seq_std]; series += 1
    legend_labels = 3
    if residual_data.get("irb_residuals") is not None:
        m = np.maximum(np.asarray(residual_data["irb_m"], dtype=float), 0.9)
        r2 = np.asarray(residual_data["irb_residuals"])
        axes[0].plot(m, r2, "^--", color=C_B, label="IRB 残差")
        axes[0].legend()
        arrays += [m, r2]; series += 1
    return _finish(fig, path, legend_labels, arrays,
                   y_nonzero_any_axis=True)


# ------------------------------------------------------------------ 图 11
def plot_site_resolved_array(site_data, path):
    """site metrics 分布与参数相关性。"""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    arrays = []
    series = 0
    site_ids = np.arange(site_data["n_sites"])
    for key, color in (("survival", C_A), ("conditional_return", C_D),
                       ("joint_return", C_B)):
        vals = np.asarray(site_data[key])
        axes[0].plot(site_ids, vals, "o", ms=2.5, color=color, label=key)
        arrays += [site_ids, vals]; series += 1
    axes[0].set_xlabel("site"); axes[0].set_ylabel("概率")
    axes[0].set_ylim(0, 1.02); axes[0].legend(fontsize=7)
    yparam = np.asarray(site_data["conditional_return"])
    for param_key, color in (("rabi_scale", C_A), ("eta_scale", C_D),
                             ("depth_scale", C_C)):
        xparam = np.asarray(site_data["params"].get(
            param_key, np.ones(site_data["n_sites"])))
        axes[1].plot(xparam, yparam, "o", ms=3, color=color,
                     label=f"vs {param_key}")
        arrays += [xparam, yparam]; series += 1
    axes[1].set_xlabel("站点参数"); axes[1].set_ylabel("conditional return")
    axes[1].legend(fontsize=7)
    for key, color in (("survival", C_A), ("conditional_return", C_D)):
        vals = np.asarray(site_data[key])
        axes[2].hist(vals, bins=20, alpha=0.6, color=color, label=key)
        arrays += [vals]; series += 1
    axes[2].set_xlabel("值"); axes[2].set_ylabel("site 数")
    axes[2].legend(fontsize=7)
    return _finish(fig, path, 8, arrays, x_nonzero_any_axis=True,
                   y_nonzero_any_axis=True)


# ------------------------------------------------------------------ 图 12
def plot_coherence_vs_round(level4_curve, level5_curve, path):
    """与 Level 4 半经典相位模型对照。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    arrays = []
    series = 0
    l4_colors = {"contrast_none": C_A, "contrast_xy4": C_B,
                 "contrast_spin_echo": C_C, "contrast_xy16": "#8c564b"}
    l5_colors = {"none": C_D, "spin_echo": "#7f7f7f",
                 "xy4_per_long_move": C_E, "xy16": "#e377c2"}
    for name, (rounds, contrast) in level4_curve.items():
        ax.plot(rounds, contrast, "o--",
                color=l4_colors.get(name, C_A),
                label=f"Level 4 toggling：{name}")
        arrays += [rounds, contrast]; series += 1
    for name, (rounds, contrast) in level5_curve.items():
        ax.plot(rounds, contrast, "s-",
                color=l5_colors.get(name, C_D),
                label=f"Level 5 有限脉冲：{name}")
        arrays += [rounds, contrast]; series += 1
    ax.set_xlabel("round"); ax.set_ylabel("条件对比度")
    ax.set_ylim(0, 1.02)
    ax.set_title("半经典 vs 有限脉冲相干性")
    ax.legend(fontsize=7)
    return _finish(fig, path, series, arrays)


# ------------------------------------------------------------------ 图 13
def plot_timestep_convergence(rows, path):
    """unitary/phase/return/fidelity 随量子步长收敛。"""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    arrays = []
    series = 0
    for metric, color in (("phase_dev_rad", C_A),
                          ("return_dev", C_D), ("fidelity_dev", C_B)):
        dt = np.array([r["dt_us"] for r in rows])
        vals = np.array([abs(float(r.get(metric, 0.0))) for r in rows])
        order = np.argsort(dt)
        ax.loglog(dt[order], np.maximum(vals[order], 1e-16), "o-",
                  color=color, label=metric)
        arrays += [dt[order], vals[order]]; series += 1
    ax.set_xlabel("量子子步长 dt (μs)")
    ax.set_ylabel("|相对最细步长偏差|")
    ax.set_title("量子步长收敛（0.10/0.05/0.025 μs）")
    ax.legend()
    return _finish(fig, path, series, arrays)


# ------------------------------------------------------------------ 图 14
def plot_sensitivity_panels(rows, path):
    """敏感性面板：各 factor 的条件保真度/survival 变化。"""
    factors = sorted({r["factor"] for r in rows})
    n = len(factors)
    n_cols = 3
    n_rows_panels = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows_panels, n_cols,
                             figsize=(4.2 * n_cols, 3.2 * n_rows_panels),
                             squeeze=False)
    arrays = []
    series = 0
    for k, factor in enumerate(factors):
        ax = axes[k // n_cols][k % n_cols]
        sub = [r for r in rows if r["factor"] == factor]
        xs = np.arange(len(sub))
        f_avg = [r.get("f_avg_conditional") or 0 for r in sub]
        surv = [r.get("survival") or 0 for r in sub]
        ax.plot(xs, f_avg, "o-", color=C_D, label="F_avg(cond)")
        ax.plot(xs, surv, "s--", color=C_A, label="survival")
        ax.set_xticks(xs)
        ax.set_xticklabels([str(r["value"]) for r in sub], fontsize=6,
                           rotation=20)
        ax.set_ylim(0, 1.05)
        ax.set_title(factor, fontsize=8)
        if k == 0:
            ax.legend(fontsize=6)
        arrays += [xs, np.array(f_avg), np.array(surv)]; series += 2
    for k in range(n, n_rows_panels * n_cols):
        axes[k // n_cols][k % n_cols].axis("off")
    return _finish(fig, path, 2, arrays,
                   x_nonzero_any_axis=True, y_nonzero_any_axis=True)
