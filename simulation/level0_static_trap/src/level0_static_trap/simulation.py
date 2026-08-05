"""单势阱模拟编排与 ``--trap`` 模式输出控制。

每个被选中的势阱都遵循同一流程：建网格 -> 计算势/力 -> 求平衡 -> 解析/数值
曲率 -> 四个频率 -> 时间步与束缚检查 -> velocity-Verlet 轨迹 -> 能量与频率
-> 结构化结果。``--trap`` 严格决定生成哪些输出，绝不伪造未选中的势阱。
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Callable

import numpy as np

from . import io_utils, visualization
from .analysis import (
    estimate_frequency_from_trajectory,
    find_equilibrium_position,
    gaussian_potential,
    summarize_static_trap,
)
from .config import Config, TrapConfig
from .constants import j_to_uk
from .integrators import velocity_verlet
from .potentials import gaussian_force


@dataclass
class TrapResult:
    """单势阱的运行结果（含原始数组与已汇总指标）。"""

    name: str
    trap: TrapConfig
    mass_kg: float
    grid_m: np.ndarray
    potential_J: np.ndarray
    force_N: np.ndarray
    times_s: np.ndarray
    positions_m: np.ndarray
    velocities_m_per_s: np.ndarray
    accelerations_m_per_s2: np.ndarray
    kinetic_energy_J: np.ndarray
    potential_energy_J: np.ndarray
    total_energy_J: np.ndarray
    energy_error_J: np.ndarray
    energy_error_over_depth: np.ndarray
    abs_energy_error_over_depth: np.ndarray
    metrics: dict


def build_grid_for_trap(trap: TrapConfig, cfg: Config) -> np.ndarray:
    """为单势阱构建网格：``center + [min,max]*waist``。"""
    lo = trap.center_m + cfg.min_offset_waist * trap.waist_m
    hi = trap.center_m + cfg.max_offset_waist * trap.waist_m
    return np.linspace(lo, hi, cfg.num_points)


def build_shared_grid(traps: list[TrapConfig], cfg: Config) -> np.ndarray:
    """构建覆盖每个势阱 ``center+[min,max]*waist`` 范围的共享网格。

    取所有势阱范围的并集 [min_lo, max_hi]，点数沿用配置。
    """
    los = [t.center_m + cfg.min_offset_waist * t.waist_m for t in traps]
    his = [t.center_m + cfg.max_offset_waist * t.waist_m for t in traps]
    return np.linspace(min(los), max(his), cfg.num_points)


def _time_array(cfg: Config) -> np.ndarray:
    """等间隔时间网格。"""
    # 用 np.arange 可能因浮点累积导致末端少一点；显式构造保证等间隔。
    n_steps = int(round(cfg.duration_s / cfg.dt_s))
    if n_steps < 1:
        raise ValueError(f"duration/dt 过小：n_steps={n_steps}")
    times = np.linspace(0.0, n_steps * cfg.dt_s, n_steps + 1)
    return times


def run_single_trap(trap: TrapConfig, cfg: Config) -> TrapResult:
    """对单个势阱运行完整流程并返回 :class:`TrapResult`。"""
    mass = cfg.mass_kg

    # 1) 网格与势/力
    grid = build_grid_for_trap(trap, cfg)
    U = gaussian_potential(grid, trap.depth_J, trap.waist_m, trap.center_m)
    F = gaussian_force(grid, trap.depth_J, trap.waist_m, trap.center_m)

    # 2) 平衡位置（与汇总内重复一次，但这里也用于诊断）
    _ = find_equilibrium_position(
        lambda x: gaussian_potential(x, trap.depth_J, trap.waist_m, trap.center_m),
        lambda x: gaussian_force(x, trap.depth_J, trap.waist_m, trap.center_m),
        trap.center_m,
        trap.waist_m,
    )

    # 3) 初始条件：x(0) = x_c + q0*w
    x0 = trap.center_m + cfg.initial_displacement_waist * trap.waist_m
    v0 = cfg.initial_velocity_m_per_s

    # 4) 束缚检查：E0 = U(x0) + 0.5 m v0^2 必须 < 0
    e0_pot = float(gaussian_potential(x0, trap.depth_J, trap.waist_m, trap.center_m))
    e0_kin = 0.5 * mass * v0 * v0
    e0 = e0_pot + e0_kin
    if e0 >= 0.0:
        # 给出清晰中文错误并停止该势阱的轨迹频率分析
        raise BoundTrajectoryError(
            f"[{trap.name}] 初始能量 E0={e0:g} J (={j_to_uk(e0):g} μK) >= 0，"
            f"轨迹非束缚（势能 {j_to_uk(e0_pot):g} μK + 动能 {j_to_uk(e0_kin):g} μK）。"
            f"请减小 initial_displacement_waist 或 initial_velocity_m_per_s。"
        )

    # 5) 时间步 & 周期数检查（已在 config 校验过 omega*dt；这里再验周期数）
    from .potentials import trap_frequency_analytic

    f0 = trap_frequency_analytic(trap.depth_J, trap.waist_m, mass)
    times = _time_array(cfg)
    est_cycles = f0 * (times[-1] - times[0])
    if est_cycles < cfg.minimum_complete_cycles:
        raise ValueError(
            f"[{trap.name}] 预估完整周期数 {est_cycles:g} < 最小要求 "
            f"{cfg.minimum_complete_cycles}（f0={f0:g} Hz, 时长 {times[-1]:g} s）。"
        )

    # 6) 积分
    F_fn: Callable = lambda x: gaussian_force(x, trap.depth_J, trap.waist_m, trap.center_m)
    _, x_traj, v_traj, a_traj = velocity_verlet(F_fn, mass, x0, v0, times)

    # 7) 能量
    ke = 0.5 * mass * v_traj * v_traj
    pe = gaussian_potential(x_traj, trap.depth_J, trap.waist_m, trap.center_m)
    total = ke + pe
    e0_actual = float(total[0])
    dE = total - e0_actual
    dE_over = dE / trap.depth_J

    # 8) 汇总指标
    metrics = summarize_static_trap(
        name=trap.name,
        depth_J=trap.depth_J,
        waist_m=trap.waist_m,
        center_m=trap.center_m,
        mass_kg=mass,
        times=times,
        positions=x_traj,
        velocities=v_traj,
        accelerations=a_traj,
        step_waist=cfg.curvature_step_waist,
        frequency_estimator=cfg.frequency_estimator,
        minimum_complete_cycles=cfg.minimum_complete_cycles,
        curv_tol=cfg.curvature_relative_tolerance,
        harm_tol=cfg.harmonic_frequency_relative_tolerance,
        traj_tol=cfg.trajectory_frequency_relative_tolerance,
        energy_tol=cfg.max_energy_error_over_depth,
        sources={
            "depth_source": trap.depth_source,
            "waist_source": trap.waist_source,
            "isotope": cfg.isotope,
        },
        waist_assumed=trap.waist_assumed,
    )
    # 补回 mass_u（汇总内未填）
    metrics["inputs"]["mass_u"] = cfg.mass_u

    return TrapResult(
        name=trap.name,
        trap=trap,
        mass_kg=mass,
        grid_m=grid,
        potential_J=U,
        force_N=F,
        times_s=times,
        positions_m=x_traj,
        velocities_m_per_s=v_traj,
        accelerations_m_per_s2=a_traj,
        kinetic_energy_J=ke,
        potential_energy_J=pe,
        total_energy_J=total,
        energy_error_J=dE,
        energy_error_over_depth=dE_over,
        abs_energy_error_over_depth=np.abs(dE_over),
        metrics=metrics,
    )


class BoundTrajectoryError(ValueError):
    """初始条件非束缚。"""


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------

TRAJECTORY_COLUMNS = [
    "time_s",
    "time_us",
    "position_m",
    "displacement_from_center_m",
    "position_um",
    "velocity_m_per_s",
    "acceleration_m_per_s2",
    "kinetic_energy_J",
    "potential_energy_J",
    "total_energy_J",
    "energy_error_J",
    "energy_error_over_depth",
    "abs_energy_error_over_depth",
]


def _trajectory_columns(res: TrapResult) -> dict[str, np.ndarray]:
    return {
        "time_s": res.times_s,
        "time_us": res.times_s * 1.0e6,
        "position_m": res.positions_m,
        "displacement_from_center_m": res.positions_m - res.trap.center_m,
        "position_um": res.positions_m * 1.0e6,
        "velocity_m_per_s": res.velocities_m_per_s,
        "acceleration_m_per_s2": res.accelerations_m_per_s2,
        "kinetic_energy_J": res.kinetic_energy_J,
        "potential_energy_J": res.potential_energy_J,
        "total_energy_J": res.total_energy_J,
        "energy_error_J": res.energy_error_J,
        "energy_error_over_depth": res.energy_error_over_depth,
        "abs_energy_error_over_depth": res.abs_energy_error_over_depth,
    }


def write_outputs(
    cfg: Config,
    trap_choice: str,
    results: dict[str, TrapResult],
) -> dict:
    """根据 ``trap_choice`` 写出全部输出，并返回图像检查结果字典。

    ``results`` 的键是已实际运行的势阱名（小写：slm / aod）。
    """
    out_dir = io_utils.ensure_output_dir(cfg.output_dir)
    image_results: dict[str, dict] = {}

    # ---- config_used.yaml ----
    from .config import config_to_used_yaml

    used = config_to_used_yaml(cfg)
    used["runtime"]["trap_choice"] = trap_choice
    io_utils.write_yaml(os.path.join(out_dir, "config_used.yaml"), used)

    # ---- static_grid.csv ----
    if trap_choice == "both":
        shared = build_shared_grid([cfg.slm, cfg.aod], cfg)
        cols: dict[str, np.ndarray] = {
            "x_m": shared,
            "x_um": shared * 1.0e6,
        }
        if "slm" in results:
            slm = results["slm"]
            cols["slm_potential_J"] = gaussian_potential(shared, slm.trap.depth_J, slm.trap.waist_m, slm.trap.center_m)
            cols["slm_potential_uK"] = j_to_uk(cols["slm_potential_J"])
            cols["slm_force_N"] = gaussian_force(shared, slm.trap.depth_J, slm.trap.waist_m, slm.trap.center_m)
        if "aod" in results:
            aod = results["aod"]
            cols["aod_potential_J"] = gaussian_potential(shared, aod.trap.depth_J, aod.trap.waist_m, aod.trap.center_m)
            cols["aod_potential_uK"] = j_to_uk(cols["aod_potential_J"])
            cols["aod_force_N"] = gaussian_force(shared, aod.trap.depth_J, aod.trap.waist_m, aod.trap.center_m)
        if cfg.save_csv:
            io_utils.write_csv(os.path.join(out_dir, "static_grid.csv"), cols)
    else:
        # 单势阱：只保留公共坐标列与该势阱列
        key = "slm" if trap_choice == "slm" else "aod"
        if key not in results:
            raise ValueError(f"trap_choice={trap_choice} 但缺少运行结果 {key}")
        res = results[key]
        grid = res.grid_m
        cols = {
            "x_m": grid,
            "x_um": grid * 1.0e6,
            f"{key}_potential_J": res.potential_J,
            f"{key}_potential_uK": j_to_uk(res.potential_J),
            f"{key}_force_N": res.force_N,
        }
        if cfg.save_csv:
            io_utils.write_csv(os.path.join(out_dir, "static_grid.csv"), cols)

    # ---- 轨迹 CSV ----
    if cfg.save_csv:
        for key, res in results.items():
            io_utils.write_csv(
                os.path.join(out_dir, f"{key}_trajectory.csv"),
                _trajectory_columns(res),
            )

    # ---- 图像 ----
    if cfg.save_png:
        # 势能比较（both 模式：共享网格）
        if trap_choice == "both":
            shared = build_shared_grid([cfg.slm, cfg.aod], cfg)
            traps_meta = {}
            if "slm" in results:
                traps_meta["SLM"] = {
                    "uK": j_to_uk(gaussian_potential(shared, results["slm"].trap.depth_J, results["slm"].trap.waist_m, results["slm"].trap.center_m)),
                    "color": "C0",
                }
            if "aod" in results:
                traps_meta["AOD"] = {
                    "uK": j_to_uk(gaussian_potential(shared, results["aod"].trap.depth_J, results["aod"].trap.waist_m, results["aod"].trap.center_m)),
                    "color": "C1",
                }
            p = os.path.join(out_dir, "potential_comparison.png")
            meta = visualization.plot_potential_comparison(p, shared * 1.0e6, traps_meta)
            image_results["potential_comparison.png"] = io_utils.validate_image(p, meta)

            traps_force = {}
            if "slm" in results:
                traps_force["SLM"] = {
                    "force_N": gaussian_force(shared, results["slm"].trap.depth_J, results["slm"].trap.waist_m, results["slm"].trap.center_m),
                    "color": "C0",
                }
            if "aod" in results:
                traps_force["AOD"] = {
                    "force_N": gaussian_force(shared, results["aod"].trap.depth_J, results["aod"].trap.waist_m, results["aod"].trap.center_m),
                    "color": "C1",
                }
            pf = os.path.join(out_dir, "force_comparison.png")
            meta = visualization.plot_force_comparison(pf, shared * 1.0e6, traps_force)
            image_results["force_comparison.png"] = io_utils.validate_image(pf, meta)
        else:
            key = "slm" if trap_choice == "slm" else "aod"
            res = results[key]
            label = "SLM" if key == "slm" else "AOD"
            p = os.path.join(out_dir, "potential_comparison.png")
            meta = visualization.plot_potential_comparison(
                p, res.grid_m * 1.0e6, {label: {"uK": j_to_uk(res.potential_J), "color": "C0"}}
            )
            image_results["potential_comparison.png"] = io_utils.validate_image(p, meta)
            pf = os.path.join(out_dir, "force_comparison.png")
            meta = visualization.plot_force_comparison(
                pf, res.grid_m * 1.0e6, {label: {"force_N": res.force_N, "color": "C0"}}
            )
            image_results["force_comparison.png"] = io_utils.validate_image(pf, meta)

        # 单势阱轨迹与能量图
        for key, res in results.items():
            label = "SLM" if key == "slm" else "AOD"
            pt = os.path.join(out_dir, f"{key}_trajectory.png")
            meta = visualization.plot_single_trajectory(
                pt,
                res.times_s * 1.0e6,
                (res.positions_m - res.trap.center_m) * 1.0e6,
                label,
            )
            image_results[f"{key}_trajectory.png"] = io_utils.validate_image(pt, meta)

            pe = os.path.join(out_dir, f"{key}_energy.png")
            meta = visualization.plot_energy_evolution(
                pe,
                res.times_s * 1.0e6,
                res.energy_error_over_depth,
                label,
                energy_tol=cfg.max_energy_error_over_depth,
            )
            image_results[f"{key}_energy.png"] = io_utils.validate_image(pe, meta)

        # 频率比较图（基于实际运行结果）
        if results:
            freq_metrics = {("SLM" if k == "slm" else "AOD"): v.metrics for k, v in results.items()}
            pfc = os.path.join(out_dir, "frequency_comparison.png")
            meta = visualization.plot_frequency_comparison(pfc, freq_metrics)
            image_results["frequency_comparison.png"] = io_utils.validate_image(pfc, meta)

    # ---- image_validation.json ----
    io_utils.write_image_validation(os.path.join(out_dir, "image_validation.json"), {
        "images": image_results,
        "all_pass": bool(image_results) and all(img["checks"]["overall_pass"] for img in image_results.values()),
    })

    # ---- metrics.json ----
    if cfg.save_json:
        io_utils.write_json(
            os.path.join(out_dir, "metrics.json"),
            {
                "trap_choice": trap_choice,
                "traps": {k: v.metrics for k, v in results.items()},
            },
        )

    return image_results


def build_summary_md(
    cfg: Config,
    trap_choice: str,
    results: dict[str, TrapResult],
    image_results: dict[str, dict],
) -> str:
    """生成简洁中文 Markdown 汇总报告。"""
    lines: list[str] = []
    lines.append("# Level 0 静态势阱自洽检查汇总报告\n")
    lines.append("## 一、范围与边界\n")
    lines.append("- `x` 为高斯光束焦平面径向坐标，**不**是光轴方向；不模拟轴向 5.64 kHz 运动。")
    lines.append("- 仅静态、一维、经典、单原子；不含 SLM-AOD 交接、深度时变、运输、二维、噪声、")
    lines.append("  量子内态、加热、损失或多原子。")
    lines.append("- SLM 与 AOD 调用同一高斯势实现，仅由配置参数与来源区分。")
    lines.append("- 数值曲率 vs 解析曲率为**代码一致性检查**，非独立实验验证。")
    lines.append("- 默认 SLM(180 μK) 与 AOD(280 μK) 为两个独立静态势阱的 Level 0 比较，")
    lines.append("  不代表论文 Fig.6 任意时刻的实际重叠势（该处 SLM 约 140 μK）。\n")

    lines.append("## 二、参数\n")
    lines.append(f"- 原子：{cfg.isotope}，质量 {cfg.mass_u:g} u = {cfg.mass_kg:g} kg。")
    lines.append(f"- 网格：{cfg.num_points} 点，偏移 [{cfg.min_offset_waist:g}, {cfg.max_offset_waist:g}] w。")
    lines.append(f"- 轨迹：初始位移 {cfg.initial_displacement_waist:g} w，初速 {cfg.initial_velocity_m_per_s:g} m/s，"
                 f"时长 {cfg.duration_s*1e6:g} μs，dt {cfg.dt_s*1e6:g} μs。\n")

    lines.append("## 三、主要结果\n")
    header = "| 势阱 | 深度(μK) | 束腰(μm,来源) | f_解析(kHz) | f_数值曲率(kHz) | f_非线性(kHz) | f_轨迹(kHz) | 离散化误差 | 非简谐频移 | max|ΔE|/U0 |"
    sep = "|---|---|---|---|---|---|---|---|---|---|"
    lines.append(header)
    lines.append(sep)
    for key, res in results.items():
        m = res.metrics
        f = m["frequencies"]
        e = m["error_classification"]
        en = m["energy"]
        waist_src = f"{res.trap.waist_um:g}" + ("(assumed)" if res.trap.waist_assumed else "")
        lines.append(
            f"| {res.name} | {res.trap.depth_uK:g} | {waist_src} | "
            f"{f['analytic_small_amp_Hz']/1e3:.4f} | {f['numerical_curvature_Hz']/1e3:.4f} | "
            f"{f['nonlinear_quadrature_Hz']/1e3:.4f} | {f['trajectory_estimated_Hz']/1e3:.4f} | "
            f"{e['discretization_error_relative']:.2e} | {e['anharmonic_shift_relative']:.3e} | "
            f"{en['max_abs_energy_error_over_depth']:.2e} |"
        )
    lines.append("")

    lines.append("## 四、误差分类说明\n")
    lines.append("- **数值离散化误差**：数值曲率频率 vs 解析小振幅频率（目标 < 1e-4）。")
    lines.append("- **有限振幅非简谐频移**：小振幅频率 vs 非线性积分频率（约 1–3%，物理软化，非积分误差）。")
    lines.append("- **轨迹积分误差**：轨迹频率 vs 同能量非线性积分频率（严格数值容差）。\n")

    lines.append("## 五、验收状态\n")
    all_ok = True
    for key, res in results.items():
        v = res.metrics["validation"]
        lines.append(f"- **{res.name}**：曲率容差={v['curvature_relative_error_within_tolerance']}，"
                     f"谐振频率容差={v['harmonic_frequency_within_tolerance']}，"
                     f"轨迹频率容差(3%)={v['trajectory_frequency_within_tolerance']}，"
                     f"能量误差容差={v['max_energy_error_within_tolerance']}，"
                     f"束缚={v['is_bound']}，周期数充足={v['cycle_count_sufficient']}。")
        all_ok = all_ok and all(v.values())
    lines.append(f"\n全部验收项通过：{all_ok}\n")

    lines.append("## 六、图像检查\n")
    n_img = len(image_results)
    n_pass = sum(1 for r in image_results.values() if r["checks"].get("overall_pass"))
    lines.append(f"- 程序化图像检查：{n_pass}/{n_img} 张通过。")
    lines.append("- 已完成程序化图像检查。视觉审阅状态见 `image_validation.json` 中的 "
                 "`visual_review_performed` 字段（默认 false，除非显式记录）。\n")

    lines.append("## 七、尚未覆盖的物理\n")
    lines.append("- 论文轴向 5.64 kHz 运动、SLM-AOD 交接、时变深度、长距离运输、")
    lines.append("  AOD 柱面透镜效应、二维运动、噪声、量子内态演化、加热、损失、多原子。\n")

    return "\n".join(lines)
