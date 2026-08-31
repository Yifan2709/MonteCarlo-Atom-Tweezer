"""Bound thermal initial states in the exact static AOD Gaussian."""
from dataclasses import dataclass
import math
from level0_static_trap.constants import BOLTZMANN_CONSTANT, MICROKELVIN_TO_KELVIN
from level0_static_trap.potentials import gaussian_potential
@dataclass(frozen=True)
class ThermalInitialState:
    x0_m: float
    v0_m_s: float
    aod_total_energy_J: float
    aod_excitation_energy_J: float
    aod_energy_over_depth: float
def thermal_scales(config):
    omega = math.sqrt(4*config.aod.depth_j/(config.atom.mass_kg*config.aod.waist_m**2)); kbt = BOLTZMANN_CONSTANT*config.thermal.temperature_uK*MICROKELVIN_TO_KELVIN
    return math.sqrt(kbt/(config.atom.mass_kg*omega**2)), math.sqrt(kbt/config.atom.mass_kg)
def sample_thermal_initial_state(config, rng, max_attempts=100000, return_attempts=False):
    sigma_x, sigma_v = thermal_scales(config)
    for attempts in range(1, max_attempts+1):
        x = config.aod.initial_center_m + rng.normal(0, sigma_x); v = rng.normal(0, sigma_v)
        energy = 0.5*config.atom.mass_kg*v**2 + float(gaussian_potential(x, config.aod.depth_j, config.aod.waist_m, config.aod.initial_center_m))
        if energy < 0:
            excitation = energy + config.aod.depth_j
            state = ThermalInitialState(x, v, energy, excitation, energy/config.aod.depth_j)
            return (state, attempts) if return_attempts else state
    raise RuntimeError("could not sample bound state")
