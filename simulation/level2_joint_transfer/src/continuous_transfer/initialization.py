"""Prepare an initially bound ensemble in the actual, already drawn SLM sites.

Harmonic proposals are conditioned on binding in the full Gaussian potential.
This remains an approximate preparation model, not a finite-depth canonical
ensemble or an independent calibration of experimental temperature.
"""
from __future__ import annotations

import numpy as np

KB = 1.380649e-23


def slm_energies(pos, vel, physics, depth_factors):
    depth = physics["slm_depth_j"] * np.asarray(depth_factors)
    waist = physics["slm_waist_m"]
    zr = np.pi * waist ** 2 / 1061e-9
    x = pos[:, 0] - physics["slm_center_m"]
    radial = x ** 2 + pos[:, 1] ** 2
    factor = 1 + (pos[:, 2] / zr) ** 2
    potential = -depth / factor * np.exp(-2 * radial / (waist ** 2 * factor))
    return .5 * physics["mass_kg"] * np.sum(vel ** 2, axis=1) + potential


def sample_site_bound_harmonic(seed, temperature_uK, physics, depth_factors,
                               *, dimensions=3, max_attempts=200):
    """Keep one atom per original site; retry rejected states at that SAME site."""
    factors = np.asarray(depth_factors, dtype=float)
    if factors.ndim != 1 or not len(factors) or np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        raise ValueError("SLM site depths must be a nonempty positive finite vector")
    if not np.isfinite(temperature_uK) or temperature_uK < 0 or dimensions not in (1, 3):
        raise ValueError("Expected nonnegative finite temperature and dimensions 1 or 3")
    if max_attempts < 1:
        raise ValueError("max_attempts must be positive")
    depth = physics["slm_depth_j"] * factors
    waist = physics["slm_waist_m"]
    zr = np.pi * waist ** 2 / 1061e-9
    thermal = KB * temperature_uK * 1e-6
    sigma_r = np.sqrt(thermal * waist ** 2 / (4 * depth))
    sigma_z = np.sqrt(thermal * zr ** 2 / (2 * depth))
    sigma_v = np.sqrt(thermal / physics["mass_kg"])
    rng = np.random.default_rng(seed)
    pos = np.zeros((len(factors), 3))
    vel = np.zeros_like(pos)
    attempts = np.zeros(len(factors), dtype=int)
    pending = np.arange(len(factors))
    for _ in range(max_attempts):
        proposals = rng.standard_normal((len(pending), 2 * dimensions))
        pos[pending, 0] = physics["slm_center_m"] + sigma_r[pending] * proposals[:, 0]
        if dimensions == 3:
            pos[pending, 1] = sigma_r[pending] * proposals[:, 1]
            pos[pending, 2] = sigma_z[pending] * proposals[:, 2]
        vel[pending, :dimensions] = sigma_v * proposals[:, dimensions:]
        attempts[pending] += 1
        energy = slm_energies(pos[pending], vel[pending], physics, factors[pending])
        pending = pending[(energy >= 0) | ~np.isfinite(energy)]
        if not len(pending):
            break
    if len(pending):
        raise RuntimeError(f"Could not prepare bound atoms at {len(pending)} SLM sites")
    return {"pos_m": pos, "vel_m_per_s": vel,
            "initial_energy_j": slm_energies(pos, vel, physics, factors),
            "attempts_per_site": attempts,
            "sampler": f"site_bound_harmonic_{dimensions}d_v1"}
