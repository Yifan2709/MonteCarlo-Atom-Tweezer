import math

import numpy as np
import pytest

from level0_static_trap.analysis import (
    calculate_relative_error,
    estimate_frequency_from_trajectory,
    estimate_nonlinear_frequency_by_quadrature,
    find_equilibrium_position,
    numerical_curvature,
    numerical_trap_frequency,
)
from level0_static_trap.constants import atomic_mass_to_kg, microkelvin_to_joule
from level0_static_trap.integrators import velocity_verlet
from level0_static_trap.potentials import (
    gaussian_curvature_analytic,
    gaussian_force,
    gaussian_potential,
    trap_frequency_analytic,
)


DEPTH = microkelvin_to_joule(180.0)
WAIST = 1.17e-6
MASS = atomic_mass_to_kg(132.90545196)
CENTER = 0.43e-6


def _trajectory(amplitude_waist: float, duration=1.0e-3, dt=0.1e-6):
    times = np.arange(round(duration / dt) + 1) * dt
    return velocity_verlet(
        lambda x: gaussian_force(x, DEPTH, WAIST, CENTER), MASS,
        CENTER + amplitude_waist * WAIST, 0.0, times,
    )


def test_nonzero_center_equilibrium_and_force_residual():
    result = find_equilibrium_position(DEPTH, WAIST, CENTER)
    assert abs(result["position_m"] - CENTER) / WAIST < 1e-8
    assert result["force_residual_N"] / (DEPTH / WAIST) < 1e-10


def test_two_step_curvature_converges_to_analytic_value():
    result = numerical_curvature(
        lambda x: gaussian_potential(x, DEPTH, WAIST, CENTER), CENTER, 0.002 * WAIST
    )
    analytic = gaussian_curvature_analytic(DEPTH, WAIST)
    error_h = calculate_relative_error(result["curvature_h_J_per_m2"], analytic)
    error_h2 = calculate_relative_error(result["curvature_h_over_2_J_per_m2"], analytic)
    assert error_h2 < error_h
    assert error_h2 < 3e-6  # 严格的曲率离散化容差
    assert result["relative_convergence_difference"] < 1e-5


def test_numerical_curvature_frequency_matches_analytic_frequency():
    result = numerical_curvature(
        lambda x: gaussian_potential(x, DEPTH, WAIST, CENTER), CENTER, 0.002 * WAIST
    )
    _, numerical = numerical_trap_frequency(result["curvature_h_over_2_J_per_m2"], MASS)
    _, analytic = trap_frequency_analytic(DEPTH, WAIST, MASS)
    assert calculate_relative_error(numerical, analytic) < 2e-6  # 严格数值一致性容差


def test_small_amplitude_trajectory_matches_harmonic_frequency():
    time, position, _, _ = _trajectory(0.01)
    estimated = estimate_frequency_from_trajectory(time, position, CENTER)["frequency_Hz"]
    analytic = trap_frequency_analytic(DEPTH, WAIST, MASS)[1]
    assert calculate_relative_error(estimated, analytic) < 2e-4  # 含有限采样和 0.01w 软化


def test_finite_amplitude_shift_is_physical_and_trajectory_matches_quadrature():
    time, position, velocity, _ = _trajectory(0.10)
    energy = 0.5 * MASS * velocity[0] ** 2 + gaussian_potential(position[0], DEPTH, WAIST, CENTER)
    nonlinear = estimate_nonlinear_frequency_by_quadrature(DEPTH, WAIST, CENTER, MASS, float(energy))
    trajectory = estimate_frequency_from_trajectory(time, position, CENTER)["frequency_Hz"]
    harmonic = trap_frequency_analytic(DEPTH, WAIST, MASS)[1]
    shift = (nonlinear["frequency_Hz"] - harmonic) / harmonic
    assert -0.02 < shift < -0.001  # 非简谐软化，不是积分器错误
    assert calculate_relative_error(trajectory, nonlinear["frequency_Hz"]) < 3e-4


def test_fft_reports_only_coarse_raw_resolution():
    time, position, _, _ = _trajectory(0.01, duration=0.5e-3)
    estimate = estimate_frequency_from_trajectory(time, position, CENTER)
    assert estimate["fft_raw_resolution_Hz"] == pytest.approx(2000.0)
    assert estimate["complete_cycle_intervals"] >= 10


def test_zero_reference_relative_error_is_explicitly_rejected():
    with pytest.raises(ValueError, match="参考量为零"):
        calculate_relative_error(1.0, 0.0)

