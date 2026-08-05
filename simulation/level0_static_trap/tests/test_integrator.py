import math

import numpy as np
import pytest

from level0_static_trap.constants import atomic_mass_to_kg, microkelvin_to_joule
from level0_static_trap.integrators import velocity_verlet
from level0_static_trap.potentials import gaussian_force, gaussian_potential, trap_frequency_analytic


def test_velocity_verlet_phase_error_has_second_order_convergence():
    omega = 2.0 * math.pi * 1000.0
    mass = 2.0
    duration = 4.3 / 1000.0

    def endpoint_error(steps_per_period: int) -> float:
        dt = 1.0 / (1000.0 * steps_per_period)
        times = np.arange(round(duration / dt) + 1) * dt
        _, position, _, _ = velocity_verlet(lambda x: -mass * omega**2 * x, mass, 1.0, 0.0, times)
        exact = math.cos(omega * times[-1])
        return abs(position[-1] - exact)

    coarse = endpoint_error(50)
    fine = endpoint_error(100)
    assert 3.5 < coarse / fine < 4.5


def test_gaussian_energy_error_scales_second_order_when_dt_halves():
    depth = microkelvin_to_joule(180.0)
    waist = 1.17e-6
    mass = atomic_mass_to_kg(132.90545196)
    duration = 200e-6

    def peak_to_peak(dt: float) -> float:
        times = np.arange(round(duration / dt) + 1) * dt
        _, x, v, _ = velocity_verlet(
            lambda position: gaussian_force(position, depth, waist, 0.0), mass, 0.1 * waist, 0.0, times
        )
        energy = 0.5 * mass * v**2 + gaussian_potential(x, depth, waist, 0.0)
        return float(np.ptp((energy - energy[0]) / depth))

    coarse = peak_to_peak(0.2e-6)
    fine = peak_to_peak(0.1e-6)
    assert coarse < 1e-3 and fine < 1e-3
    assert 3.5 < coarse / fine < 4.5


def test_integrator_rejects_nonuniform_or_nonincreasing_times():
    with pytest.raises(ValueError, match="等间隔"):
        velocity_verlet(lambda x: -x, 1.0, 1.0, 0.0, [0.0, 0.1, 0.21])
    with pytest.raises(ValueError, match="严格递增"):
        velocity_verlet(lambda x: -x, 1.0, 1.0, 0.0, [0.0, 0.1, 0.1])

