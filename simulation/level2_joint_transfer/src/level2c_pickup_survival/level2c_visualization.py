"""Level 2C 诊断图：非交互后端，复用 Level 2 的程序化验收记录。"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# macOS CJK 字体回退，避免中文标签显示为方框
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS", "Hiragino Sans GB",
                                   "STHeiti", "PingFang SC", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

from level2_joint_transfer.level2_visualization import _finish

K_B_UK = 1.380649e-29

C_ML = "#1f77b4"
C_600 = "#606060"
C_400 = "#a0a0a0"
C_200 = "#d0d0d0"
C_SIM = "#d62728"


def plot_pickup_waveforms(waveforms: dict, path: Path):
    """各时长 pick-up 波形的中心与深度（实线）及时间反演 drop-off（虚线）。"""
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    arrays = []
    for label, wf in waveforms.items():
        table = wf.dense_table(1501)
        t_us = table["time_s"] * 1e6
        axes[0].plot(t_us, table["aod_center_m"] * 1e6, label=label)
        axes[1].plot(t_us, table["aod_depth_J"] / K_B_UK, label=label)
        arrays += [t_us, table["aod_center_m"], table["aod_depth_J"]]
    axes[0].set_ylabel("AOD 中心 (µm)")
    axes[1].set_ylabel("AOD 深度 (µK)")
    axes[1].set_xlabel("时间 (µs)")
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    axes[0].set_title("pick-up 手工轨迹（二次升深 + constant-jerk 移动）")
    return _finish(fig, axes, path, 2 * len(waveforms), arrays)


def plot_roundtrip_program(legs, wait_us, path: Path):
    """一个往返的完整势场程序：AOD 中心/深度与 SLM 开关。"""
    fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    n_pu, n_wait, n_do = legs.steps_pickup, legs.steps_wait, legs.steps_dropoff
    t = np.concatenate([np.arange(n_pu + 1),
                        n_pu + np.arange(n_wait + 1),
                        n_pu + n_wait + np.arange(n_do + 1)])
    center = np.concatenate([legs.center_pickup, legs.center_wait, legs.center_dropoff])
    depth = np.concatenate([legs.depth_pickup, legs.depth_wait, legs.depth_dropoff])
    factor = legs.slm_factor
    axes[0].plot(t, center * 1e6)
    axes[1].plot(t, depth / K_B_UK)
    axes[2].step(t[:factor.size], factor, where="post")
    axes[0].set_ylabel("AOD 中心 (µm)")
    axes[1].set_ylabel("AOD 深度 (µK)")
    axes[2].set_ylabel("SLM 开关")
    axes[2].set_xlabel(f"步（往返 = pick-up + {wait_us:.0f} µs 等待 + drop-off）")
    return _finish(fig, axes, path, 0, [t, center, depth, factor])


def plot_survival_vs_transfers(curves_sim: dict, reference: dict, tolerance_band: float,
                               path: Path, title_suffix=""):
    """核心对比图：模拟 survival(n)（实线带标记）叠加论文数字化曲线（灰阶）。

    curves_sim: {label: (n, survival)}；reference: {series: (n, survival)}。
    """
    fig, ax = plt.subplots(figsize=(8, 5.5))
    arrays = []
    ref_style = {"ml_400us": ("#6baed6", "论文 0.4 ms (ML)"),
                 "manual_600us": ("#737373", "论文 0.6 ms 手工"),
                 "manual_400us": ("#bdbdbd", "论文 0.4 ms 手工"),
                 "manual_200us": ("#e2e2e2", "论文 0.2 ms 手工")}
    for key, (n, s) in reference.items():
        color, label = ref_style.get(key, ("#999999", key))
        ax.plot(n, s, "o--", color=color, ms=3, lw=1.2, label=label)
        arrays += [n, s]
    for label, (n, s) in curves_sim.items():
        ax.plot(n, s, "-s", ms=3, lw=1.6, label=label)
        arrays += [np.asarray(n), np.asarray(s)]
    ax.set_xlabel("单程转移次数 n")
    ax.set_ylabel("存活率")
    ax.set_ylim(0.0, 1.05)
    ax.legend(fontsize=8)
    ax.set_title(f"atom survival vs 单程转移次数：模拟 vs 论文 Fig. 6d 顶图{title_suffix}")
    return _finish(fig, ax, path, len(reference) + len(curves_sim), arrays)


def plot_survival_zoom(curves: dict, path: Path, reference: dict | None = None,
                       ylim=(0.3, 1.05)):
    """保守模型 survival(n) 与论文参考曲线同图（诚实展示数值差距）。

    论文曲线以浅灰虚线作背景参照；模拟曲线预计钉在 1 附近。
    reference 提供时一并绘制（同时保证图内 y 数据范围非零）。
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    arrays = []
    series = 0
    if reference:
        ref_style = {"ml_400us": ("#6baed6", "论文 0.4 ms (ML)"),
                     "manual_600us": ("#737373", "论文 0.6 ms 手工"),
                     "manual_400us": ("#bdbdbd", "论文 0.4 ms 手工"),
                     "manual_200us": ("#e2e2e2", "论文 0.2 ms 手工")}
        for key, (n, s) in reference.items():
            color, label = ref_style.get(key, ("#999999", key))
            ax.plot(n, s, "o--", color=color, ms=3, lw=1.0, alpha=0.7, label=label)
            arrays += [n, s]
            series += 1
    for label, (n, s) in curves.items():
        ax.plot(n, s, "-o", ms=3, label=label)
        arrays += [np.asarray(n), np.asarray(s)]
        series += 1
    ax.set_xlabel("单程转移次数 n")
    ax.set_ylabel("存活率")
    ax.set_ylim(*ylim)
    ax.legend(fontsize=8)
    ax.set_title("保守模型 survival(n)（无损耗机制 ≈1；灰阶虚线为论文 Fig. 6d 参考）")
    return _finish(fig, ax, path, series, arrays)


def plot_ml_optimization(rows: list, path: Path):
    """ML 候选评估历史：存活率与状态随评估序号。"""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    idx = np.arange(len(rows))
    survival = np.array([r.get("survival_at_objective", np.nan) for r in rows], dtype=float)
    status = [r.get("status") for r in rows]
    ok = np.array([s == "ok" for s in status])
    ax.plot(idx[ok], survival[ok], "o", ms=4, label="有效候选")
    invalid = np.array([s == "invalid" for s in status])
    ax.plot(idx[invalid], np.zeros(int(invalid.sum())), "x", ms=4, label="无效（过冲/约束）")
    failed = np.array([s == "failed" for s in status])
    ax.plot(idx[failed], np.full(int(failed.sum()), -0.02), "v", ms=4, label="运行失败")
    ax.axhline(1.0, color="k", lw=0.5)
    ax.set_xlabel("评估序号")
    ax.set_ylabel(f"目标存活率")
    ax.legend(fontsize=8)
    ax.set_title("ML 轨迹优化历史（LHS + 局部精修）")
    finite_arrays = [idx[ok], survival[ok], idx[invalid], idx[failed],
                     np.zeros(int(invalid.sum())), np.full(int(failed.sum()), -0.02)]
    return _finish(fig, ax, path, 3, finite_arrays)


def plot_ml_trajectory(manual_wf, ml_wf, path: Path):
    """ML 冻结轨迹 vs 手工 400 µs：中心/深度及其控制点。"""
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    t_m = manual_wf.dense_table(1501)["time_s"] * 1e6
    axes[0].plot(t_m, manual_wf.center(t_m * 1e-6) * 1e6, label="手工 400 µs")
    axes[1].plot(t_m, manual_wf.depth(t_m * 1e-6) / K_B_UK, label="手工 400 µs")
    t_ml = ml_wf.dense_table(1501)["time_s"] * 1e6
    s_ctrl = ml_wf.control_s * ml_wf.duration_s * 1e6
    axes[0].plot(t_ml, ml_wf.center(t_ml * 1e-6) * 1e6, label="ML 400 µs（冻结）")
    axes[0].plot(s_ctrl, ml_wf.center_values_m * 1e6, "ko", ms=4)
    axes[1].plot(t_ml, ml_wf.depth(t_ml * 1e-6) / K_B_UK, label="ML 400 µs（冻结）")
    axes[1].plot(s_ctrl, ml_wf.depth_values_j / K_B_UK, "ko", ms=4)
    axes[0].set_ylabel("AOD 中心 (µm)")
    axes[1].set_ylabel("AOD 深度 (µK)")
    axes[1].set_xlabel("时间 (µs)")
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    return _finish(fig, axes, path, 4,
                   [t_m, t_ml, s_ctrl, ml_wf.center_values_m, ml_wf.depth_values_j])


def plot_heating_accumulation(variants: dict, path: Path):
    """加热/差异系综变体的能量累积：按波形分 2×2 面板，每面板 4 个损耗变体。

    锯齿为奇/偶判定点系统差异（奇数次在 AOD 阱判定、偶数次在 SLM 阱判定，
    两阱深度不同归一化不同）——真实物理，非噪声。
    """
    keys = sorted(variants)
    waveforms = sorted({k.split("|")[1] for k in keys})
    variant_names = sorted({k.split("|")[0] for k in keys})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    axes = axes.ravel()
    arrays = []
    series_total = 0
    for ax, wname in zip(axes, waveforms):
        for vname in variant_names:
            key = f"{vname}|{wname}"
            if key not in variants:
                continue
            energies = variants[key]
            max_len = max((len(e) for e in energies), default=0)
            if max_len == 0:
                continue
            mean_curve, ns = [], []
            for k in range(max_len):
                vals = [e[k] for e in energies if len(e) > k]
                mean_curve.append(float(np.mean(vals)))
                ns.append(k + 1)
            ax.plot(ns, mean_curve, "-", lw=1.2, label=vname)
            arrays += [np.asarray(ns), np.asarray(mean_curve)]
            series_total += 1
        ax.set_yscale("log")
        ax.set_title(wname, fontsize=10)
        ax.legend(fontsize=7)
    for ax in axes[2:]:
        ax.set_xlabel("单程转移次数 n")
    for ax in (axes[0], axes[2]):
        ax.set_ylabel("判定时刻平均归一化激发 ε/D")
    fig.suptitle("重复转移中的能量累积（存活 shot 系综平均；锯齿=奇偶判定阱不同）")
    return _finish(fig, axes, path, series_total, arrays)


def plot_loss_attribution(result: dict, path: Path):
    """丢失归因：逐转移丢失计数按 leg（pick-up 段 / drop-off 段）。"""
    lost_leg = result["lost_leg"]
    lost_transfer = result["lost_transfer"]
    legs = sorted(set(x for x in lost_leg if x))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    arrays = []
    for leg in legs:
        ns = [t for l, t in zip(lost_leg, lost_transfer) if l == leg and t > 0]
        if ns:
            ax.hist(ns, bins=range(1, max(ns) + 2), alpha=0.6, label=f"{leg} 段丢失")
            arrays.append(np.asarray(ns, dtype=float))
    ax.set_xlabel("丢失发生的单程转移序号")
    ax.set_ylabel("丢失 shot 数")
    ax.legend(fontsize=8)
    ax.set_title("丢失归因（按段）")
    return _finish(fig, ax, path, len(legs), arrays if arrays else [np.array([0.0])])
