"""Composition of the verified Level 0 Gaussian trap functions."""
from level0_static_trap.potentials import gaussian_force, gaussian_potential
from .waveforms import aod_center, aod_depth
def slm_potential(x, config): return gaussian_potential(x, config.slm.depth_j, config.slm.waist_m, config.slm.center_m)
def aod_potential(x, t, config): return gaussian_potential(x, aod_depth(t, config), config.aod.waist_m, aod_center(t, config))
def total_potential(x, t, config): return slm_potential(x, config) + aod_potential(x, t, config)
def total_force(x, t, config): return gaussian_force(x, config.slm.depth_j, config.slm.waist_m, config.slm.center_m) + gaussian_force(x, aod_depth(t, config), config.aod.waist_m, aod_center(t, config))
