"""SLM/AOD 共用的一维静态势阱模拟流水线。"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .analysis import (
    calculate_relative_error,
    estimate_frequency_from_trajectory,
    estimate_nonlinear_frequency_by_quadrature,
    find_equilibrium_position,
    numerical_curvature,
    numerical_trap_frequency,
    summarize_static_trap,
)
from .config import SimulationConfig, TrapConfig
from .constants import BOLTZMANN_CONSTANT
from .integrators import velocity_verlet
from .potentials import gaussian_curvature_analytic, gaussian_force, gaussian_potential, trap_frequency_analytic


SCOPE_LIMITATIONS = [
    "x 是高斯光束焦平面内的一维径向坐标，不是光轴方向。",
    "不模拟论文报告的轴向 5.64 kHz 运动。",
    "不模拟 SLM-AOD 交接、时变深度、长距离移动、AOD 柱面透镜效应或二维运动。",
    "不包含噪声、量子内态、加热、损失或多原子效应。",
    "SLM 与 AOD 调用同一高斯势函数，仅由配置参数及来源区分。",
    "解析/数值曲率比较是代码一致性检查，不是独立实验验证。",
]


def _selected_traps(config: SimulationConfig, selection: str) -> dict[str, TrapConfig]:
    if selection == "slm":
        return {"slm": config.slm}
    if selection == "aod":
        return {"aod": config.aod}
    if selection == "both":
        return {"slm": config.slm, "aod": config.aod}
    raise ValueError("trap 必须是 slm、aod 或 both")


def _shared_grid(config: SimulationConfig, traps: dict[str, TrapConfig]) -> np.ndarray:
    lower = min(t.center_m + config.grid.min_offset_waist * t.waist_m for t in traps.values())
    upper = max(t.center_m + config.grid.max_offset_waist * t.waist_m for t in traps.values())
    return np.linspace(lower, upper, config.grid.num_points)


def _run_trap(config: SimulationConfig, trap: TrapConfig, grid: np.ndarray) -> dict[str, Any]:
    depth, waist, center, mass = trap.depth_j, trap.waist_m, trap.center_m, config.atom.mass_kg
    potential = gaussian_potential(grid, depth, waist, center)
    force = gaussian_force(grid, depth, waist, center)
    equilibrium = find_equilibrium_position(depth, waist, center)
    analytic_curvature = gaussian_curvature_analytic(depth, waist)
    curvature = numerical_curvature(
        lambda x: gaussian_potential(x, depth, waist, center), center, config.grid.curvature_step_waist * waist
    )
    analytic_omega, analytic_frequency = trap_frequency_analytic(depth, waist, mass)
    numerical_omega, numerical_frequency = numerical_trap_frequency(curvature["curvature_h_over_2_J_per_m2"], mass)

    trajectory_cfg = config.trajectory
    initial_position = center + trajectory_cfg.initial_displacement_waist * waist
    initial_velocity = trajectory_cfg.initial_velocity_m_per_s
    initial_energy = 0.5 * mass * initial_velocity**2 + float(gaussian_potential(initial_position, depth, waist, center))
    if initial_energy >= 0.0:
        raise ValueError(f"{trap.name} 初始状态不束缚：E0={initial_energy:.6e} J，必须满足 E0<0")
    nonlinear = estimate_nonlinear_frequency_by_quadrature(depth, waist, center, mass, initial_energy)
    n_steps = int(round(trajectory_cfg.duration_s / trajectory_cfg.dt_s))
    times = np.arange(n_steps + 1, dtype=float) * trajectory_cfg.dt_s
    if not np.all(np.diff(times) > 0.0) or not np.allclose(np.diff(times), trajectory_cfg.dt_s):
        raise ValueError("构造的时间数组不是严格递增等间隔数组")
    omega_dt = analytic_omega * trajectory_cfg.dt_s
    if omega_dt >= config.validation.max_omega_dt:
        raise ValueError(
            f"{trap.name} 时间步长检查失败：ωΔt={omega_dt:.6g}，阈值={config.validation.max_omega_dt:.6g}"
        )
    estimated_cycles = trajectory_cfg.duration_s * analytic_frequency
    if estimated_cycles < trajectory_cfg.minimum_complete_cycles:
        raise ValueError(
            f"{trap.name} 预计完整周期数 {estimated_cycles:.3f} 小于要求 {trajectory_cfg.minimum_complete_cycles}"
        )

    time, position, velocity, acceleration = velocity_verlet(
        lambda x: gaussian_force(x, depth, waist, center), mass, initial_position, initial_velocity, times
    )
    kinetic = 0.5 * mass * velocity**2
    potential_trajectory = gaussian_potential(position, depth, waist, center)
    total = kinetic + potential_trajectory
    energy_error = total - initial_energy
    energy_error_over_depth = energy_error / depth
    frequency = estimate_frequency_from_trajectory(time, position, center)

    numeric_frequency_error = calculate_relative_error(numerical_frequency, analytic_frequency)
    finite_amplitude_shift = (nonlinear["frequency_Hz"] - analytic_frequency) / analytic_frequency
    trajectory_to_nonlinear = calculate_relative_error(frequency["frequency_Hz"], nonlinear["frequency_Hz"])
    trajectory_to_harmonic = calculate_relative_error(frequency["frequency_Hz"], analytic_frequency)
    max_energy_error = float(np.max(np.abs(energy_error_over_depth)))
    rms_energy_error = float(np.sqrt(np.mean(energy_error_over_depth**2)))
    peak_to_peak_energy_error = float(np.ptp(energy_error_over_depth))
    force_scale = depth / waist

    validation = {
        "equilibrium_center_passed": abs(equilibrium["position_m"] - center) / waist < 1.0e-8,
        "force_residual_passed": equilibrium["force_residual_N"] / force_scale < 1.0e-10,
        "curvature_convergence_passed": curvature["relative_convergence_difference"]
        < config.validation.curvature_relative_tolerance,
        "harmonic_frequency_consistency_passed": numeric_frequency_error
        < config.validation.harmonic_frequency_relative_tolerance,
        "trajectory_finite_amplitude_band_passed": trajectory_to_harmonic
        < config.validation.trajectory_frequency_relative_tolerance,
        "trajectory_nonlinear_reference_passed": trajectory_to_nonlinear
        < min(3.0e-3, config.validation.trajectory_frequency_relative_tolerance / 5.0),
        "energy_error_passed": max_energy_error < config.validation.max_energy_error_over_depth,
        "time_step_passed": omega_dt < config.validation.max_omega_dt,
        "cycle_count_passed": estimated_cycles >= trajectory_cfg.minimum_complete_cycles,
        "initial_state_bound_passed": initial_energy < 0.0,
    }
    validation["all_passed"] = all(validation.values())

    metrics = summarize_static_trap(
        input_parameters={
            "name": trap.name,
            "depth": {"value": trap.depth_uK, "unit": "uK", "si_value_J": depth, "source": trap.depth_source},
            "waist": {
                "value": trap.waist_um, "unit": "um", "si_value_m": waist, "source": trap.waist_source,
                "assumed": trap.waist_is_assumed,
            },
            "center": {"value": trap.center_um, "unit": "um", "si_value_m": center},
            "mass": {"value": config.atom.mass_u, "unit": "u", "si_value_kg": mass, "isotope": config.atom.isotope},
        },
        equilibrium={
            **equilibrium,
            "center_offset_m": equilibrium["position_m"] - center,
            "force_residual_over_depth_per_waist": equilibrium["force_residual_N"] / force_scale,
        },
        curvature={
            "analytic_J_per_m2": analytic_curvature,
            **curvature,
            "h_over_2_relative_error_vs_analytic": calculate_relative_error(
                curvature["curvature_h_over_2_J_per_m2"], analytic_curvature
            ),
        },
        frequencies={
            "analytic_small_amplitude_Hz": analytic_frequency,
            "analytic_omega_rad_per_s": analytic_omega,
            "numerical_curvature_Hz": numerical_frequency,
            "numerical_omega_rad_per_s": numerical_omega,
            "nonlinear_quadrature_Hz": nonlinear["frequency_Hz"],
            "trajectory_estimate_Hz": frequency["frequency_Hz"],
            "trajectory_period_s": frequency["period_s"],
            "fft_peak_Hz": frequency["fft_peak_frequency_Hz"],
            "fft_raw_resolution_Hz": frequency["fft_raw_resolution_Hz"],
            "zero_crossing_complete_cycle_intervals": frequency["complete_cycle_intervals"],
            "estimated_duration_cycles": estimated_cycles,
        },
        error_classification={
            "numerical_discretization_relative_error": numeric_frequency_error,
            "finite_amplitude_nonlinear_shift_relative": finite_amplitude_shift,
            "trajectory_vs_nonlinear_reference_relative_error": trajectory_to_nonlinear,
            "trajectory_vs_harmonic_relative_difference": trajectory_to_harmonic,
            "note": "前两项分别是数值离散化误差与物理非简谐频移，不得合并解释。",
        },
        trajectory={
            "initial_position_m": initial_position,
            "initial_displacement_waist": trajectory_cfg.initial_displacement_waist,
            "initial_velocity_m_per_s": initial_velocity,
            "initial_total_energy_J": initial_energy,
            "bound": initial_energy < 0.0,
            "nonlinear_turning_amplitude_m": nonlinear["turning_amplitude_m"],
            "max_abs_displacement_m": float(np.max(np.abs(position - center))),
            "max_abs_velocity_m_per_s": float(np.max(np.abs(velocity))),
            "omega_dt": omega_dt,
        },
        energy_error={
            "normalization": "trap depth U0 (positive)",
            "max_abs_energy_error_over_depth": max_energy_error,
            "rms_energy_error_over_depth": rms_energy_error,
            "peak_to_peak_energy_error_over_depth": peak_to_peak_energy_error,
        },
        validation=validation,
    )
    arrays = {
        "time_s": time,
        "position_m": position,
        "velocity_m_per_s": velocity,
        "acceleration_m_per_s2": acceleration,
        "kinetic_energy_J": kinetic,
        "potential_energy_J": potential_trajectory,
        "total_energy_J": total,
        "energy_error_J": energy_error,
        "energy_error_over_depth": energy_error_over_depth,
        "potential_grid_J": potential,
        "force_grid_N": force,
    }
    return {"config": trap, "metrics": metrics, "arrays": arrays}


def run_simulation(config: SimulationConfig, selection: str = "both") -> dict[str, Any]:
    """对所选势阱执行完全相同的静态分析和经典轨迹流程。"""

    traps = _selected_traps(config, selection)
    grid = _shared_grid(config, traps)
    results = {key: _run_trap(config, trap, grid) for key, trap in traps.items()}
    return {
        "selection": selection,
        "coordinate": config.model.coordinate,
        "scope_and_limitations": SCOPE_LIMITATIONS,
        "grid_m": grid,
        "traps": results,
    }


def build_summary(simulation: dict[str, Any], image_validation: dict[str, Any]) -> str:
    """生成包含边界、数值结果和验证状态的简洁中文 Markdown 报告。"""

    lines = [
        "# Level 0 静态径向光镊汇总",
        "",
        "## 模型边界",
        "",
    ]
    lines.extend(f"- {item}" for item in simulation["scope_and_limitations"])
    lines.extend(["", "SLM 180 μK 与 AOD 280 μK 是两个独立静态势阱的比较，不代表实际重叠转移势。AOD 的 1.17 μm 束腰是与 SLM 后孔径匹配的假设值。", "", "## 主要结果", ""])
    lines.append("| Trap | Analytic (kHz) | Curvature (kHz) | Nonlinear (kHz) | Trajectory (kHz) | Numerical error | Nonlinear shift | Max abs(dE)/U0 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for key, result in simulation["traps"].items():
        metric = result["metrics"]
        freq = metric["frequencies"]
        errors = metric["error_classification"]
        energy = metric["energy_error"]
        lines.append(
            f"| {key.upper()} | {freq['analytic_small_amplitude_Hz']/1e3:.6f} | "
            f"{freq['numerical_curvature_Hz']/1e3:.6f} | {freq['nonlinear_quadrature_Hz']/1e3:.6f} | "
            f"{freq['trajectory_estimate_Hz']/1e3:.6f} | {errors['numerical_discretization_relative_error']:.3e} | "
            f"{errors['finite_amplitude_nonlinear_shift_relative']:.3%} | {energy['max_abs_energy_error_over_depth']:.3e} |"
        )
    lines.extend([
        "",
        "数值曲率误差属于数值离散化一致性检查；有限振幅频移是高斯势非简谐软化的物理效应，二者没有合并成单一‘相对误差’。能量误差以正的势阱深度 U0 归一化，报告的是有界误差幅度，不称为长期漂移。",
        "",
        "## 验证",
        "",
    ])
    for key, result in simulation["traps"].items():
        status = "通过" if result["metrics"]["validation"]["all_passed"] else "失败"
        lines.append(f"- {key.upper()} 数值验收：{status}。")
    image_status = "通过" if image_validation["all_passed"] else "失败"
    lines.append(f"- PNG 程序化图像检查：{image_status}（{len(image_validation['images'])} 张）。")
    if image_validation.get("visual_review_performed"):
        findings = "；".join(image_validation.get("visual_review_findings", [])) or "未发现明显布局问题"
        lines.append(f"- 已执行模型视觉审阅：{findings}。")
    else:
        lines.append("- 已完成程序化图像检查，未进行人工或模型视觉审阅。")
    lines.extend([
        "",
        "## 未包含的后续物理",
        "",
        "本 Level 0 不包含 SLM-AOD 交接、时变势深、运输、AOD 柱面透镜效应、二维/轴向运动、噪声、量子内态、加热、损失和多原子效应。若扩展到 F(x,t)，必须重新审视积分器接口、能量解释和辛结构。",
        "",
    ])
    return "\n".join(lines)
