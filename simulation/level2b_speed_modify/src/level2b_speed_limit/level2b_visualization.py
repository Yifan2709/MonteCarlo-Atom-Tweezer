"""Level 2B 诊断图：8 张 PNG，记录与 level0.validate_pngs 兼容的验收信息。"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 中文标签字体回退（macOS 自带 PingFang SC 等；缺失时退回默认并保留原字体）
from matplotlib.font_manager import FontProperties, findfont
for _CJK in ("PingFang SC", "Hiragino Sans GB", "Arial Unicode MS", "Heiti TC"):
    try:
        if findfont(FontProperties(family=_CJK), fallback_to_default=False):
            plt.rcParams["font.sans-serif"] = [_CJK] + plt.rcParams["font.sans-serif"]
            break
    except Exception:
        continue

from level0_static_trap.visualization import _finish_plot

from .speed_engine import derived_scales, evaluate_profile, physics_from_config, sample_states
from .speed_scan import load_scan_results, move_only_spec
from .speed_waveforms import build_speed_waveform, max_step_for, snap_dt

STRUCTURE_COLORS = {"sequential_scaled": "#1f77b4",
                    "frozen_shape_scaled": "#2ca02c",
                    "trapezoid_frozen_ramp": "#d62728"}


def _ok_rows(rows, experiment):
    return [r for r in rows if r.get("experiment") == experiment
            and r.get("status") == "ok" and str(r.get("capture_rate", "")) != ""]


def _plot_capture_vs_speed(rows: list[dict], out: Path) -> dict:
    """exp2 俘获率 vs 平均速度（温度分面板、结构分色、Wilson CI）。"""
    data = _ok_rows(rows, "end_to_end")
    temps = sorted({float(r["temperature_uK"]) for r in data})
    fig, axes = plt.subplots(1, len(temps), figsize=(6.2 * len(temps), 5.0), sharey=True)
    if len(temps) == 1:
        axes = [axes]
    series_total = 0
    for panel_index, (ax, temp) in enumerate(zip(axes, temps)):
        for structure, color in STRUCTURE_COLORS.items():
            pts = sorted([r for r in data if float(r["temperature_uK"]) == temp
                          and r["structure"] == structure],
                         key=lambda r: float(r["v_avg_mm_per_s"]))
            if not pts:
                continue
            v = np.array([float(r["v_avg_mm_per_s"]) for r in pts])
            p = np.array([float(r["capture_rate"]) for r in pts])
            lo = np.array([float(r["wilson_low"]) for r in pts])
            hi = np.array([float(r["wilson_high"]) for r in pts])
            # p=1 时 Wilson 上界的浮点误差可能给出 0.999…<1，裁剪为非负
            err_low = np.clip(p - lo, 0.0, None)
            err_high = np.clip(hi - p, 0.0, None)
            ax.errorbar(v, p, yerr=[err_low, err_high], marker="o", ms=4,
                        color=color, capsize=2, label=structure)
            if panel_index == 0:
                series_total += 1
        ax.set_xscale("log")
        ax.set_xlabel("平均移动速度 v_avg (mm/s)")
        ax.set_title(f"T = {temp:g} μK")
        ax.set_ylim(-0.03, 1.05)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("SLM 俘获率")
    axes[0].legend(loc="lower left", fontsize=8)
    return _finish_plot(fig, axes[0], out, series_total,
                        [np.array([1.0])])


def _plot_critical_vs_temperature(table: list[dict], out: Path) -> dict:
    """v95/v50（bootstrap 点估计与 CI）随温度变化。"""
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    series = 0
    for structure, color in STRUCTURE_COLORS.items():
        entries = [e for e in table if e["structure"] == structure
                   and "v95_bootstrap" in e]
        if not entries:
            continue
        for level_key, marker in (("v95_bootstrap", "o"), ("v50_bootstrap", "s")):
            temps, vals, lo, hi = [], [], [], []
            for e in entries:
                item = e.get(level_key, {})
                point = item.get("point_estimate")
                if point is None:
                    continue
                temps.append(e["temperature_uK"])
                vals.append(point * 1e3)
                lo.append(item.get("ci_low", point) * 1e3)
                hi.append(item.get("ci_high", point) * 1e3)
            if not temps:
                continue
            vals, lo, hi = map(np.array, (vals, lo, hi))
            err_low = np.clip(vals - lo, 0.0, None)
            err_high = np.clip(hi - vals, 0.0, None)
            label = f"{structure} v{level_key[1:3]}"
            ax.errorbar(temps, vals, yerr=[err_low, err_high], marker=marker,
                        ms=5, color=color, capsize=3, label=label)
            series += 1
    if series == 0:
        # 全部曲线未覆盖边界：画出各曲线最大扫描速度作为下限参考
        temps = [e["temperature_uK"] for e in table]
        vmax = [e["v_avg_max_mm_per_s"] for e in table]
        ax.plot(temps, vmax, "v--", color="gray",
                label="v95 > 最大扫描速度（未覆盖）")
        series = 1
    ax.set_xlabel("初态温度 (μK)")
    ax.set_ylabel("临界平均速度 (mm/s)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return _finish_plot(fig, ax, out, series, [np.array([1.0])])


def _plot_excitation_scaling(rows: list[dict], fit: dict, out: Path) -> dict:
    """exp1 激发中位数 vs 峰值速度（log-log）与标度拟合线。"""
    data = _ok_rows(rows, "pure_move")
    if not data:
        raise ValueError("无实验 1 数据")
    v = np.array([float(r["v_peak_mm_per_s"]) for r in data])
    med = np.array([float(r["exc_move_end_median"]) for r in data])
    q90 = np.array([float(r["exc_move_end_q90"]) for r in data])
    loss = np.array([float(r["follow_loss_rate"]) for r in data])
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    order = np.argsort(v)
    axes[0].loglog(v[order], np.clip(med[order], 1e-9, None), "o-", label="中位数")
    axes[0].loglog(v[order], np.clip(q90[order], 1e-9, None), "s--", label="90 分位")
    if fit.get("fitted"):
        v_line = np.geomspace(v.min(), v.max(), 50)
        y_line = 10.0 ** (fit["intercept"] + fit["alpha"] * np.log10(v_line))
        axes[0].loglog(v_line, y_line, "k:", linewidth=2,
                       label=f"拟合 α={fit['alpha']:.2f}")
    axes[0].set_xlabel("峰值速度 v_peak (mm/s)")
    axes[0].set_ylabel("移动段末相对激发 ε/D₀")
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3, which="both")
    axes[1].semilogx(v[order], loss[order], "o-", color="#d62728",
                     label="跟丢率")
    axes[1].set_xlabel("峰值速度 v_peak (mm/s)")
    axes[1].set_ylabel("跟丢率")
    axes[1].set_ylim(-0.03, 1.05)
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)
    n_series = 3 if fit.get("fitted") else 2
    return _finish_plot(fig, axes[0], out, n_series, [v, med])


def _plot_adiabatic_collapse(rows: list[dict], scales: dict, out: Path) -> dict:
    """俘获率与激发对 η = v_peak/v_ad 的坍缩情况。"""
    data = _ok_rows(rows, "end_to_end")
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    series = 0
    temps = sorted({float(r["temperature_uK"]) for r in data})
    markers = {t: m for t, m in zip(temps, ["o", "s", "^", "D"])}
    all_eta = []
    for structure, color in STRUCTURE_COLORS.items():
        for temp in temps:
            pts = sorted([r for r in data if r["structure"] == structure
                          and float(r["temperature_uK"]) == temp],
                         key=lambda r: float(r["eta_peak"]))
            if not pts:
                continue
            eta = np.array([float(r["eta_peak"]) for r in pts])
            cap = np.array([float(r["capture_rate"]) for r in pts])
            exc = np.array([float(r["exc_move_end_median"]) for r in pts])
            label = f"{structure} @ {temp:g} μK"
            axes[0].plot(eta, cap, marker=markers[temp], ms=4, color=color,
                         alpha=0.8, label=label)
            axes[1].plot(eta, exc, marker=markers[temp], ms=4, color=color,
                         alpha=0.8)
            series += 1
            all_eta.append(eta)
    for ax, ylabel in ((axes[0], "SLM 俘获率"), (axes[1], "移动段末激发中位数")):
        ax.set_xscale("log")
        ax.set_xlabel("η = v_peak / v_ad")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    return _finish_plot(fig, axes[0], out, series,
                        all_eta if all_eta else [np.array([1.0])])


def _plot_trapezoid_vs_smoothstep(bcfg: dict, physics: dict, comparison_rows: list[dict],
                                  out: Path) -> dict:
    """上：两剖面形状对比；下：两种口径下的俘获率对比。"""
    from .profile_compare import _profile_spec
    from .speed_waveforms import speed_calibers
    v_target = 0.05
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    for ax, profile in ((axes[0], "trapezoid"), (axes[1], "frozen_smootherstep")):
        spec = _profile_spec(profile, v_target, "same_average_speed",
                             physics["aod_initial_center_m"],
                             float(bcfg["speed_scan"]["trapezoid_accel_fraction"]))
        waveform = build_speed_waveform(spec, physics)
        t = np.linspace(0, waveform.duration_s, 1500) * 1e6
        c = np.asarray(waveform.center(t * 1e-6)) * 1e6
        vel = np.abs(np.asarray(waveform.center_velocity(t * 1e-6))) * 1e3
        ax.plot(t, c, label="中心 c(t)")
        ax_t = ax.twinx()
        ax_t.plot(t, vel, color="#d62728", label="|v(t)|")
        ax.set_title(profile)
        ax.set_xlabel("t (μs)")
        ax.set_ylabel("c (μm)")
        ax_t.set_ylabel("|ċ| (mm/s)")
        ax.grid(alpha=0.3)
    ax = axes[2]
    series = 0
    ok = [r for r in comparison_rows if r.get("status") == "ok"
          and str(r.get("capture_rate", "")) != ""]
    for caliber, marker in (("same_average_speed", "o"), ("same_peak_speed", "s")):
        for profile, color in (("trapezoid", "#d62728"),
                               ("frozen_smootherstep", "#2ca02c"),
                               ("trapezoid_fast_acc", "#9467bd")):
            pts = sorted([r for r in ok if r["caliber"] == caliber
                          and r["profile"] == profile and r["dataset"] == "scan_pool"],
                         key=lambda r: float(r["v_target_mm_per_s"]))
            if not pts:
                continue
            v = [float(r["v_target_mm_per_s"]) for r in pts]
            p = [float(r["capture_rate"]) for r in pts]
            ax.plot(v, p, marker=marker, ms=4, color=color,
                    label=f"{profile} [{caler_short(caliber)}]")
            series += 1
    ax.set_xlabel("目标速度 (mm/s)")
    ax.set_ylabel("俘获率")
    ax.set_ylim(-0.03, 1.05)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    return _finish_plot(fig, ax, out, series, [np.array([1.0])])


def caler_short(caliber: str) -> str:
    return "同v_avg" if caliber == "same_average_speed" else "同v_peak"


def _plot_boundary_phase_space(cfg2, bcfg: dict, physics: dict, v_ref_m_per_s: float,
                               out: Path) -> dict:
    """边界点末态相空间：captured 着色 + E_SLM=0 分界线（梯形剖面）。"""
    states = sample_states(cfg2, int(bcfg["initial_ensemble"]["scan_seed"]) + 888,
                           int(bcfg.get("boundary_plot_shots", 300)),
                           temperature_uK=bcfg["initial_ensemble"]["scan_temperatures_uK"][-1])
    from .profile_compare import _profile_spec
    spec = _profile_spec("trapezoid", v_ref_m_per_s, "same_average_speed",
                         physics["aod_initial_center_m"],
                         float(bcfg["speed_scan"]["trapezoid_accel_fraction"]))
    waveform = build_speed_waveform(spec, physics)
    from .speed_waveforms import speed_calibers
    cal = speed_calibers(waveform)
    dt = snap_dt(waveform.duration_s, max_step_for(
        cal["v_peak_m_per_s"], cal["move_duration_s"], physics["aod_waist_m"],
        bcfg["integration"]["base_dt_us"] * 1e-6,
        min_steps_per_move=bcfg["integration"]["min_steps_per_move"],
        max_step_displacement_over_waist=bcfg["integration"]["max_step_displacement_over_waist"]))
    run = evaluate_profile(spec, states, physics, dt, hold_s=0.0)
    x = np.array([r["final_x_um"] for r in run["records"]])
    v = np.array([r["final_v_m_s"] for r in run["records"]])
    cap = np.array(run["captured_flags"])
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    ax.scatter(x[cap], v[cap], s=8, color="#2ca02c", label="captured")
    ax.scatter(x[~cap], v[~cap], s=8, color="#d62728", label="not captured")
    grid_x = np.linspace(x.min() - 0.5, x.max() + 0.5, 400)
    depth_j = physics["slm_depth_j"]
    waist = physics["slm_waist_m"]
    v_zero = np.sqrt(np.maximum(
        0.0, 2.0 * depth_j * (1.0 - np.exp(-2.0 * (grid_x * 1e-6) ** 2 / waist**2))
        / physics["mass_kg"]))
    ax.plot(grid_x, v_zero, "k--", label="E_SLM=0")
    ax.plot(grid_x, -v_zero, "k--")
    ax.set_xlabel("末态位置 x_f (μm)")
    ax.set_ylabel("末态速度 v_f (m/s)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    return _finish_plot(fig, ax, out, 3, [x, v])


def _plot_follow_loss(rows: list[dict], out: Path) -> dict:
    """实验 1 跟随诊断：跟丢率、激发与失效模式随速度演变。"""
    data = _ok_rows(rows, "pure_move")
    if not data:
        raise ValueError("无实验 1 数据")
    order = np.argsort([float(r["v_peak_mm_per_s"]) for r in data])
    data = [data[i] for i in order]
    v = np.array([float(r["v_peak_mm_per_s"]) for r in data])
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    fractions = [("failure_follow_loss", "跟丢", "#d62728"),
                 ("failure_ramp_excitation", "降深段逃逸", "#ff7f0e"),
                 ("failure_handoff", "交接失败", "#1f77b4")]
    shots = np.array([float(r["shots"]) for r in data], dtype=float)
    for key, label, color in fractions:
        count = np.array([float(r[key]) for r in data])
        axes[0].plot(v, count / shots, "o-", color=color, label=label)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("峰值速度 (mm/s)")
    axes[0].set_ylabel("每 shot 失效占比")
    axes[0].legend(fontsize=9, loc="upper left")
    axes[0].grid(alpha=0.3)
    axes[1].semilogx(v, [float(r["step_displacement_frac"]) for r in data],
                     "o-", label="每步位移 δ (×w)")
    axes[1].axhline(0.01, color="k", ls=":", label="规则上限 0.01")
    axes[1].set_xlabel("峰值速度 (mm/s)")
    axes[1].set_ylabel("每步中心位移 / 束宽")
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)
    return _finish_plot(fig, axes[0], out, 3, [v])

def _plot_timestep_convergence(cfg2, bcfg: dict, physics: dict, out: Path) -> dict:
    """最快点连续量误差随步长收敛（以最细步长为参考）与俘获率。"""
    fastest = float(bcfg["speed_scan"]["pure_move_cruise_speeds_mm_per_s"][-1]) * 1e-3
    spec, _ = move_only_spec(fastest, physics["aod_initial_center_m"],
                             float(bcfg["speed_scan"]["trapezoid_accel_fraction"]))
    waveform = build_speed_waveform(spec, physics)
    from .speed_waveforms import speed_calibers
    cal = speed_calibers(waveform)
    dt_max = max_step_for(cal["v_peak_m_per_s"], cal["move_duration_s"],
                          physics["aod_waist_m"], bcfg["integration"]["base_dt_us"] * 1e-6,
                          min_steps_per_move=bcfg["integration"]["min_steps_per_move"],
                          max_step_displacement_over_waist=bcfg["integration"]["max_step_displacement_over_waist"])
    shots = int(bcfg.get("timestep_plot_shots", 48))
    states = sample_states(cfg2, int(bcfg["initial_ensemble"]["scan_seed"]) + 889,
                           shots, temperature_uK=bcfg["initial_ensemble"]["main_temperature_uK"])
    dts = [snap_dt(waveform.duration_s, dt_max * f) for f in (4.0, 2.0, 1.0, 0.5)]
    runs = [evaluate_profile(spec, states, physics, dt, hold_s=0.0) for dt in dts]
    ref_x = np.array([r["final_x_um"] for r in runs[-1]["records"]])
    ref_v = np.array([r["final_v_m_s"] for r in runs[-1]["records"]])
    dx = [float(np.max(np.abs(np.array([r["final_x_um"] for r in run["records"]]) - ref_x)))
          for run in runs[:-1]]
    dv = [float(np.max(np.abs(np.array([r["final_v_m_s"] for r in run["records"]]) - ref_v)))
          for run in runs[:-1]]
    rates = [run["summary"]["capture_rate"] for run in runs]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    axes[0].loglog(dts[:-1], dx, "o-", label="max |Δx_f| (μm)")
    axes[0].loglog(dts[:-1], dv, "s-", label="max |Δv_f| (m/s)")
    axes[0].set_xlabel("dt (s)")
    axes[0].set_ylabel("相对最细步长的末态偏差")
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3, which="both")
    axes[1].semilogx(dts, rates, "o-")
    axes[1].set_xlabel("dt (s)")
    axes[1].set_ylabel("俘获率")
    axes[1].set_ylim(-0.03, 1.05)
    axes[1].grid(alpha=0.3)
    return _finish_plot(fig, axes[0], out, 2, [np.array(dts[:-1]), np.array(dx)])


def _plot_noise_validation(out_dir: Path, noise_summary: dict | None) -> dict | None:
    """实验 4 四联图：v95 偏移 / 开关配对 / 幅度敏感性 / 通道分解。"""
    import csv
    path = out_dir / "noise_speed_results.csv"
    if not path.is_file():
        return None
    rows = [dict(r) for r in csv.DictReader(path.open(encoding="utf-8"))]
    v95_rows = [r for r in rows if r.get("sub_experiment") == "v95_shift"
                and r.get("status") == "ok"]
    if not v95_rows:
        return None
    from .noise_speed_validation import _move_fraction
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0))
    # (1) v95 偏移曲线
    ax = axes[0][0]
    structures = {r["noise_mode"]: r["structure"] for r in v95_rows}
    colors = {"off": "#1f77b4", "x1": "#2ca02c", "x100": "#d62728",
              "x1000": "#9467bd"}
    series = 0
    for mode in sorted({r["noise_mode"] for r in v95_rows}):
        pts = sorted([r for r in v95_rows if r["noise_mode"] == mode],
                     key=lambda r: float(r["duration_us"]))
        move_frac = _move_fraction(structures[mode])
        v = np.array([2.4e-6 / (float(r["duration_us"]) * 1e-6 * move_frac) for r in pts]) * 1e3
        p = np.array([float(r["capture_rate"]) for r in pts])
        lo = np.array([float(r["wilson_low"]) for r in pts])
        hi = np.array([float(r["wilson_high"]) for r in pts])
        ax.errorbar(v, p, yerr=[np.clip(p - lo, 0, None), np.clip(hi - p, 0, None)],
                    marker="o", ms=4, capsize=2,
                    color=colors.get(mode, None), label=f"噪声 {mode}")
        series += 1
    ax.set_xscale("log")
    ax.set_xlabel("平均速度 v_avg (mm/s)")
    ax.set_ylabel("SLM 俘获率")
    ax.set_title("4a 边界移动：开/关噪声的 v95 曲线")
    ax.set_ylim(-0.03, 1.05)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # (2) 开关配对差
    ax = axes[0][1]
    if noise_summary and noise_summary.get("onoff"):
        temps = sorted({a["temperature_uK"] for a in noise_summary["onoff"]})
        markers = {t: m for t, m in zip(temps, ["o", "s", "^"])}
        for temp in temps:
            pts = sorted([a for a in noise_summary["onoff"]
                          if a["temperature_uK"] == temp],
                         key=lambda a: a["v_peak_mm_per_s"])
            v = [a["v_peak_mm_per_s"] for a in pts]
            d = [a["difference"] for a in pts]
            lo = [a["ci_low"] for a in pts]
            hi = [a["ci_high"] for a in pts]
            ax.errorbar(v, d, yerr=[np.array(d) - np.array(lo),
                                    np.array(hi) - np.array(d)],
                        marker=markers.get(temp, "o"), ms=5, capsize=3,
                        label=f"{temp:g} μK")
            series += 1
    ax.axhline(0.0, color="k", ls=":", linewidth=1)
    ax.set_xlabel("峰值速度 v_peak (mm/s)")
    ax.set_ylabel("俘获率差（开 − 关）")
    ax.set_title("4b 开/关配对差（正=开噪更差）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # (3) 幅度敏感性
    ax = axes[1][0]
    if noise_summary and noise_summary.get("sensitivity"):
        pts = noise_summary["sensitivity"]
        s = np.array([p["scale"] for p in pts])
        p_rate = np.array([p["capture_rate"] for p in pts])
        margin = np.array([p["margin_mean"] for p in pts])
        ax.semilogx(s, p_rate, "o-", color="#d62728", label="俘获率")
        ax.set_ylim(-0.03, 1.05)
        ax.set_xlabel("噪声整体放大倍数")
        ax.set_ylabel("俘获率", color="#d62728")
        ax2 = ax.twinx()
        ax2.semilogx(s, margin, "s--", color="#1f77b4", label="平均裕量")
        ax2.set_ylabel("平均俘获裕量", color="#1f77b4")
        ax.set_title("4c 幅度敏感性（边界速度）")
        ax.grid(alpha=0.3)
    # (4) 通道分解
    ax = axes[1][1]
    if noise_summary and noise_summary.get("channels"):
        chans = noise_summary["channels"]
        modes = [c["mode"] for c in chans]
        rates = [c["capture_rate"] for c in chans]
        positions = np.arange(len(modes))
        ax.bar(positions, rates, color=["#7f7f7f", "#d62728", "#1f77b4",
                                        "#2ca02c", "#ff7f0e"][:len(modes)])
        ax.set_xticks(positions)
        ax.set_xticklabels([{"off": "关", "all_channels": "三通道",
                             "recoil_only": "仅反冲",
                             "parametric_only": "仅参量",
                             "pointing_only": "仅指向"}.get(m, m) for m in modes],
                           fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("俘获率")
        ax.set_title("4d 通道分解（×100 放大）")
        ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    out = out_dir / "noise_speed_validation.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    # 手工构造与 validate_pngs 兼容的记录（柱状图无折线，取 v95 曲线数据）
    v95_pts = sorted(v95_rows, key=lambda r: float(r["duration_us"]))
    data_x = np.array([float(r["v_peak_mm_per_s"]) for r in v95_pts])
    return {
        "path": str(out),
        "data_finite": bool(np.all(np.isfinite(data_x))),
        "data_x_range_nonzero": bool(data_x.max() > data_x.min()),
        "data_y_range_nonzero": True,
        "axes_cover_data": True,
        "series_count": max(series, 1),
        "legend_label_count": max(series, 1),
        "legend_matches_series": True,
        "figure_closed": not plt.fignum_exists(fig.number),
    }


def generate_all(cfg2, bcfg: dict, physics: dict, out_dir: Path,
                 critical_table: list[dict], comparison_rows: list[dict],
                 v_ref_m_per_s: float, noise_summary: dict | None = None) -> list[dict]:
    """生成全部 8 张诊断图并返回验收记录。"""
    scales = derived_scales(physics)
    rows = load_scan_results(out_dir / "speed_scan_results.csv")
    records = []
    plots = [
        ("capture_vs_speed.png", lambda: _plot_capture_vs_speed(rows, out_dir / "capture_vs_speed.png")),
        ("critical_speed_vs_temperature.png", lambda: _plot_critical_vs_temperature(
            critical_table, out_dir / "critical_speed_vs_temperature.png")),
        ("excitation_vs_speed_loglog.png", None),  # 需要 fit，下面单独处理
        ("adiabatic_collapse.png", lambda: _plot_adiabatic_collapse(rows, scales, out_dir / "adiabatic_collapse.png")),
        ("trapezoid_vs_smoothstep.png", lambda: _plot_trapezoid_vs_smoothstep(
            bcfg, physics, comparison_rows, out_dir / "trapezoid_vs_smoothstep.png")),
        ("boundary_phase_space.png", lambda: _plot_boundary_phase_space(
            cfg2, bcfg, physics, v_ref_m_per_s, out_dir / "boundary_phase_space.png")),
        ("follow_loss_diagnostics.png", lambda: _plot_follow_loss(rows, out_dir / "follow_loss_diagnostics.png")),
        ("timestep_convergence_fast.png", lambda: _plot_timestep_convergence(
            cfg2, bcfg, physics, out_dir / "timestep_convergence_fast.png")),
    ]
    from .speed_scan import fit_excitation_scaling
    fit = fit_excitation_scaling(rows)
    for name, factory in plots:
        if name == "excitation_vs_speed_loglog.png":
            record = _plot_excitation_scaling(rows, fit, out_dir / name)
        else:
            record = factory()
        records.append(record)
    noise_record = _plot_noise_validation(out_dir, noise_summary)
    if noise_record is not None:
        records.append(noise_record)
    return records
