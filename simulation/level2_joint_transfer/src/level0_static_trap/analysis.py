"""平衡位置、曲率、轨迹频率和有限振幅基准分析。"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .potentials import gaussian_force, gaussian_potential


def calculate_relative_error(value: float, reference: float) -> float:
    """计算相对误差；参考量为零时显式拒绝病态定义。"""

    if reference == 0.0:
        raise ValueError("参考量为零，无法定义相对误差")
    return abs(float(value) - float(reference)) / abs(float(reference))


def find_equilibrium_position(depth: float, waist: float, center: float) -> dict[str, float]:
    """在无量纲束腰坐标内有界搜索势阱平衡位置。"""

    # 用 F/(U0/w) 消除焦耳和牛顿的微小绝对量级，避免优化器默认绝对容差主导结果。
    normalized_force = lambda q: float(
        gaussian_force(center + q * waist, depth, waist, center) / (depth / waist)
    )
    root_q = brentq(normalized_force, -2.0, 2.0, xtol=1.0e-15, rtol=1.0e-15)
    position = center + root_q * waist
    return {
        "position_m": position,
        "potential_J": float(gaussian_potential(position, depth, waist, center)),
        "force_residual_N": abs(float(gaussian_force(position, depth, waist, center))),
    }


def numerical_curvature(
    potential_function: Callable[[float], float], center: float, step: float
) -> dict[str, float]:
    """用 h 和 h/2 中心差分计算曲率并报告收敛差异。"""

    if step <= 0.0 or not math.isfinite(step):
        raise ValueError("曲率差分步长必须是有限正数")

    def second_derivative(h: float) -> float:
        return (
            float(potential_function(center + h))
            - 2.0 * float(potential_function(center))
            + float(potential_function(center - h))
        ) / h**2

    curvature_h = second_derivative(step)
    curvature_h_over_2 = second_derivative(step / 2.0)
    if curvature_h_over_2 <= 0.0:
        raise ValueError(f"数值曲率必须为正，实际为 {curvature_h_over_2:.6e} J/m²")
    return {
        "step_h_m": step,
        "step_h_over_2_m": step / 2.0,
        "curvature_h_J_per_m2": curvature_h,
        "curvature_h_over_2_J_per_m2": curvature_h_over_2,
        "relative_convergence_difference": calculate_relative_error(curvature_h, curvature_h_over_2),
    }


def numerical_trap_frequency(curvature: float, mass: float) -> tuple[float, float]:
    """由正的数值曲率返回角频率和普通频率。"""

    if curvature <= 0.0 or mass <= 0.0:
        raise ValueError("曲率和质量必须为正")
    omega = math.sqrt(curvature / mass)
    return omega, omega / (2.0 * math.pi)


def estimate_frequency_from_trajectory(times, positions, center: float) -> dict[str, float | int]:
    """用同方向线性插值过零点回归估计轨迹频率，并给出 FFT 诊断。"""

    time = np.asarray(times, dtype=float)
    displacement = np.asarray(positions, dtype=float) - center
    if time.ndim != 1 or displacement.shape != time.shape or time.size < 4:
        raise ValueError("轨迹时间和位置必须是一维、等长且至少含 4 点")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(displacement)):
        raise ValueError("轨迹含 NaN 或 Infinity")
    upward = np.flatnonzero((displacement[:-1] <= 0.0) & (displacement[1:] > 0.0))
    crossings: list[float] = []
    for index in upward:
        y0, y1 = displacement[index], displacement[index + 1]
        fraction = -y0 / (y1 - y0)
        crossings.append(float(time[index] + fraction * (time[index + 1] - time[index])))
    if len(crossings) < 2:
        raise ValueError("同方向过零点不足，无法估计完整周期")
    crossing_array = np.asarray(crossings)
    cycle_index = np.arange(crossing_array.size, dtype=float)
    period, intercept = np.polyfit(cycle_index, crossing_array, 1)
    del intercept
    if period <= 0.0:
        raise ValueError("过零点回归得到非正周期")
    dt = float(time[1] - time[0])
    duration = float(time[-1] - time[0])
    spectrum = np.abs(np.fft.rfft(displacement - np.mean(displacement)))
    frequencies = np.fft.rfftfreq(time.size, dt)
    spectrum[0] = 0.0
    peak_index = int(np.argmax(spectrum))
    return {
        "frequency_Hz": 1.0 / float(period),
        "period_s": float(period),
        "crossing_count": len(crossings),
        "complete_cycle_intervals": len(crossings) - 1,
        "fft_peak_frequency_Hz": float(frequencies[peak_index]),
        "fft_raw_resolution_Hz": 1.0 / duration,
    }


def estimate_nonlinear_frequency_by_quadrature(
    depth: float, waist: float, center: float, mass: float, initial_energy: float
) -> dict[str, float]:
    """以同一初始能量的一维积分计算有限振幅非线性周期和频率。"""

    if initial_energy >= 0.0:
        raise ValueError(f"初始能量 {initial_energy:.6e} J 非负，轨迹不束缚")
    minimum = -depth
    if initial_energy < minimum * (1.0 + 1.0e-12):
        raise ValueError("初始能量低于势阱最低点，物理上不可能")
    if math.isclose(initial_energy, minimum, rel_tol=1.0e-14, abs_tol=0.0):
        omega = math.sqrt(4.0 * depth / (mass * waist**2))
        return {"turning_amplitude_m": 0.0, "period_s": 2.0 * math.pi / omega, "frequency_Hz": omega / (2.0 * math.pi)}

    amplitude = waist * math.sqrt(-0.5 * math.log(-initial_energy / depth))
    turning_slope = (4.0 * depth * amplitude / waist**2) * math.exp(-2.0 * (amplitude / waist) ** 2)
    endpoint_limit = math.sqrt(2.0 * amplitude / turning_slope)

    def transformed_integrand(theta: float) -> float:
        cosine = math.cos(theta)
        if abs(cosine) < 1.0e-10:
            return endpoint_limit
        position = center + amplitude * math.sin(theta)
        available = initial_energy - float(gaussian_potential(position, depth, waist, center))
        if available <= 0.0:
            return endpoint_limit
        return amplitude * abs(cosine) / math.sqrt(available)

    integral, _ = quad(transformed_integrand, -math.pi / 2.0, math.pi / 2.0, epsrel=1.0e-11, limit=200)
    period = math.sqrt(2.0 * mass) * integral
    return {"turning_amplitude_m": amplitude, "period_s": period, "frequency_Hz": 1.0 / period}


def summarize_static_trap(**metrics) -> dict:
    """将单势阱分析量整理为可序列化的结构化摘要。"""

    return dict(metrics)
