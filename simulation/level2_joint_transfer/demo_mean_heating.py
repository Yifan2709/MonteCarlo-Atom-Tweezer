"""平均强度确定性噪声加热演示：机制验证 + 随时间不稳定的定量结果。

运行：.venv/bin/python demo_mean_heating.py
输出：outputs/level2_joint_transfer/demo_mean_heating/
  - mean_heating_demo.png（机制面板 + 俘获率随转移时长衰减 + 保持段漂移）
  - mean_heating_demo.json（全部数值结果）

三通道（默认常数，均为 assumed）：
  recoil    P = 2·E_r·Γ_sc ≈ 0.36 μK/s（1061 nm、280 μK 两能级近似）
  parametric P = Γ_par·E_ref，Γ_par = (π²/4)·ν₀²·S_RIN(2ν₀) ≈ 16 /s，
    E_ref = k_B·T_ref（缺省取初始系综温度）——远离 2ν₀ 参量共振的线性化常数功率
  pointing  P = m·ω₀⁴·S_x(ν₀)/8 ≈ 1.3 μK/s（S_x = 1 pm²/Hz @25.5 kHz）
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from level0_static_trap.constants import microkelvin_to_joule
from level0_static_trap.potentials import gaussian_force, gaussian_potential
from level2_joint_transfer.config import load_config
from level2_joint_transfer.level2_simulation import (evaluate_waveform_on_states,
                                                     evaluate_waveform_parallel,
                                                     physics_from_config,
                                                     sample_initial_states)
from level2_joint_transfer.noise_heating import (build_mean_heating, heating_from_config,
                                                 scaled_heating,
                                                 velocity_verlet_time_dependent_heating)

PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "level2_joint_transfer" / "demo_mean_heating"
DT = 0.05e-6
# 记录 *_uK 字段已统一为真实 μK（J / 1.380649e-29）


def mechanism_panel(ax, physics):
    """静态 SLM 阱中三通道的解析 vs 数值增长曲线。"""
    mass, depth, waist = physics["mass_kg"], physics["slm_depth_j"], physics["slm_waist_m"]
    total = 500e-6
    times = np.arange(round(total / DT) + 1) * DT
    force = lambda x, t: gaussian_force(x, depth, waist, 0.0)

    def e_osc_fn(x, v, t):
        return max(0.0, 0.5 * mass * v * v + float(gaussian_potential(x, depth, waist, 0.0)) + depth)

    energy = lambda x, v: (0.5 * mass * v ** 2 + gaussian_potential(x, depth, waist, 0.0)) / microkelvin_to_joule(1.0)

    # 常数功率（反冲/指向类）：P = 100 μK/s，线性
    power = microkelvin_to_joule(100.0)
    _t, x, v, _a, _w = velocity_verlet_time_dependent_heating(
        force, mass, 0.2e-6, 0.0, times, build_mean_heating(recoil_power_w=power), e_osc_fn)
    numeric = energy(x, v) - energy(x, v)[0]
    ax.plot(times * 1e6, numeric, label="常数功率 100 μK/s（数值）", color="tab:blue")
    ax.plot(times * 1e6, 100.0 * times, "--", color="tab:blue", alpha=0.6,
            label="解析 P·t")

    # 参量通道（线性化常数功率）：P = Γ_par·E_ref = 1000/s × k_B·10 μK
    rate, t_ref_uK = 1000.0, 10.0
    p_par = microkelvin_to_joule(rate * t_ref_uK)  # 10 mK/s
    _t, x, v, _a, _w = velocity_verlet_time_dependent_heating(
        force, mass, 0.2e-6, 0.0, times, build_mean_heating(parametric_power_w=p_par), e_osc_fn)
    numeric_par = energy(x, v) - energy(x, v)[0]
    ax.plot(times * 1e6, numeric_par, label="参量 P=Γ·E_ref=10 mK/s（数值）", color="tab:red")
    ax.plot(times * 1e6, rate * t_ref_uK * times, "--", color="tab:red", alpha=0.6,
            label="解析 P·t")
    ax.set_xlabel("时间 (μs)")
    ax.set_ylabel("能量 (μK)")
    ax.set_title("机制：三通道均为常数功率（线性；参量为 Γ·E_ref 线性化）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    return {"constant_power_max_rel_error": float(abs(numeric[-1] / (100.0 * times[-1]) - 1.0)),
            "parametric_linear_max_rel_error": float(abs(
                numeric_par[-1] / (rate * t_ref_uK * times[-1]) - 1.0))}


def duration_scan(cfg, physics, base_heating):
    """顺序波形时长 × 噪声放大倍数：俘获率随时间衰减（系统越发不稳定）。"""
    states = sample_initial_states(cfg, seed=31337, shots=256)
    durations_us = [400.0, 600.0, 800.0, 1200.0, 1600.0]
    scales = [1, 100, 1000]
    table = []
    for scale in scales:
        heating = scaled_heating(base_heating, float(scale))
        for duration_us in durations_us:
            spec = {"type": "sequential", "duration_s": duration_us * 1e-6,
                    "move_fraction": 0.52, "name": f"seq_{duration_us:.0f}us_x{scale}"}
            evaluation = evaluate_waveform_parallel(
                dict(spec), states, physics, DT, 100e-6, heating=heating)
            records = evaluation["records"]
            table.append({
                "scale": scale, "duration_us": duration_us,
                "capture_fraction": float(np.mean([r["captured"] for r in records])),
                "mean_margin_captured": float(np.mean([r["capture_margin"] for r in records
                                                       if r["captured"]])) if any(r["captured"] for r in records) else None,
                "mean_noise_injected_uK": float(np.mean(
                    [r["noise_total_energy_uK"] for r in records])),
            })
    return table


def hold_scan(cfg, physics, base_heating):
    """固定 400 μs 转移，延长保持段：注入能量与 SLM 能量漂移随时长增长。"""
    states = sample_initial_states(cfg, seed=31338, shots=128)
    heating = scaled_heating(base_heating, 1000.0)
    rows = []
    for hold_us in (100.0, 500.0, 1000.0, 2000.0):
        spec = {"type": "sequential", "duration_s": 400e-6, "move_fraction": 0.52}
        evaluation = evaluate_waveform_on_states(dict(spec), states, physics, DT, hold_us * 1e-6,
                                                 heating=heating)
        records = evaluation["records"]
        captured = [r for r in records if r["captured"]]
        rows.append({
            "hold_us": hold_us,
            "captured_fraction": len(captured) / len(records),
            # 仅统计被俘获原子，避免逃逸原子失控能量污染均值
            "median_hold_peak_to_peak_over_depth_captured":
                float(np.median([r["hold_peak_to_peak_energy_error_over_depth"] for r in captured]))
                if captured else None,
            "mean_noise_total_at_hold_end_uK":
                float(np.mean([r["noise_total_energy_uK"] for r in records])),
        })
    return rows


def main():
    cfg = load_config(PROJECT_ROOT / "configs" / "level2_joint_transfer.yaml")
    physics = physics_from_config(cfg)
    base = heating_from_config(cfg) or build_mean_heating(
        recoil_power_w=2 * 64e-9 * microkelvin_to_joule(1.0) * 2.8,
        parametric_power_w=16.0 * microkelvin_to_joule(5.0),
        pointing_power_w=1.3 * microkelvin_to_joule(1.0))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))

    errors = mechanism_panel(axes[0], physics)

    duration_table = duration_scan(cfg, physics, base)
    for scale in (1, 100, 1000):
        pts = [row for row in duration_table if row["scale"] == scale]
        axes[1].plot([p["duration_us"] for p in pts], [p["capture_fraction"] for p in pts],
                     marker="o", label=f"噪声 ×{scale}")
    axes[1].set_xlabel("转移总时长 (μs)")
    axes[1].set_ylabel("SLM 俘获率（256 shots）")
    axes[1].set_title("俘获率随转移时长衰减（常数噪声持续注入）")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    axes[1].set_ylim(-0.03, 1.03)

    hold_rows = hold_scan(cfg, physics, base)
    axes[2].plot([r["hold_us"] for r in hold_rows],
                 [r["median_hold_peak_to_peak_over_depth_captured"] for r in hold_rows], marker="s",
                 color="tab:purple", label="保持段 SLM 能量峰峰/D_SLM（已俘获中位数）")
    axes[2].set_xlabel("保持段时长 (μs)")
    axes[2].set_ylabel("峰峰能量误差 / D_SLM")
    axes[2].set_yscale("log")
    axes[2].set_title("保持段能量漂移增长（噪声 ×1000）")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)

    fig.suptitle("Level 2 平均强度确定性噪声加热：反冲/参量/指向（常数、无随机性）", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mean_heating_demo.png", dpi=200)
    plt.close(fig)

    result = {
        "defaults_from_config": base.describe()["channels"],
        "mechanism_errors": errors,
        "duration_scan": duration_table,
        "hold_scan": hold_rows,
        "notes": {
            "units": "duration/hold 表中 mean_noise_injected_uK 即记录字段原值（真实 μK）",
            "provenance": "assumed_sensitivity_only，默认参数推导见 configs/level2_joint_transfer.yaml 注释",
        },
    }
    (OUTPUT_DIR / "mean_heating_demo.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("mechanism_errors", "duration_scan", "hold_scan")},
                     indent=2, ensure_ascii=False))
    print(f"输出目录：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
