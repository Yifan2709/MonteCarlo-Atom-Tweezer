"""Site disorder draws and bound initial-state preparation (F3).

Generalized from the research initialization: site parameters are drawn
first, then every atom is prepared *bound in its actual site* (harmonic
proposal conditioned on the drawn depth, rejection in the full Gaussian
potential, retry at the same site, site identity never discarded).
"""
from __future__ import annotations

import numpy as np

from ..errors import ToolkitError, E_BAD_VALUE
from ..units import K_B
from .field import BeamContext, beam_field


def draw_site_disorder(prep, shots: int, seed: int) -> dict:
    """Draw per-site parameters. Disabled disorder -> exact unit factors."""
    if not prep.site_disorder_enabled:
        return {"slm_depth_factor": np.ones(shots),
                "aod_depth_factor": np.ones(shots),
                "aod_waist_factor": np.ones(shots),
                "alignment_offset_m": np.zeros(shots),
                "enabled": False}
    rng = np.random.default_rng(seed)
    return {"slm_depth_factor": 1.0 + prep.slm_depth_rel_sigma * rng.standard_normal(shots),
            "aod_depth_factor": 1.0 + prep.aod_depth_rel_sigma * rng.standard_normal(shots),
            "aod_waist_factor": 1.0 + prep.aod_waist_rel_sigma * rng.standard_normal(shots),
            "alignment_offset_m": prep.alignment_sigma_m * rng.standard_normal(shots),
            "enabled": True}


def slm_site_energies(pos, vel, ctx: BeamContext, slm_depth_factor):
    """Total energy in each atom's actual SLM site (J)."""
    depth = ctx.slm_depth_j * np.asarray(slm_depth_factor)
    u = beam_field(pos[:, 0] - ctx.slm_center_m[0], pos[:, 1] - ctx.slm_center_m[1],
                   pos[:, 2], depth, ctx.slm_waist_m, ctx.slm_zr_m)[0]
    return 0.5 * ctx.mass_kg * np.sum(np.asarray(vel) ** 2, axis=1) + u


def sample_site_bound(ctx: BeamContext, prep, shots: int, seed: int,
                      slm_depth_factor, *, max_attempts: int = 200) -> dict:
    """3D site-bound harmonic preparation with in-place rejection.

    Proposals use the per-site actual depth; rejection checks E < 0 in the
    full Gaussian potential; rejected atoms are redrawn AT THE SAME site so
    the site census never changes. This remains an approximate preparation
    model, not a finite-depth canonical ensemble.
    """
    factors = np.asarray(slm_depth_factor, dtype=float)
    if factors.ndim != 1 or len(factors) != shots or np.any(factors <= 0):
        raise ToolkitError(E_BAD_VALUE, "slm depth factors must be positive, one per shot",
                           "preparation")
    depth = ctx.slm_depth_j * factors
    thermal = K_B * prep.temperature_uK * 1e-6
    sigma_r = np.sqrt(thermal * ctx.slm_waist_m ** 2 / (4.0 * depth))
    sigma_z = np.sqrt(thermal * ctx.slm_zr_m ** 2 / (2.0 * depth))
    sigma_v = np.sqrt(thermal / ctx.mass_kg)
    rng = np.random.default_rng(seed)
    pos = np.zeros((shots, 3))
    vel = np.zeros((shots, 3))
    attempts = np.zeros(shots, dtype=int)
    pending = np.arange(shots)
    for _ in range(max_attempts):
        proposals = rng.standard_normal((len(pending), 6))
        pos[pending, 0] = ctx.slm_center_m[0] + sigma_r[pending] * proposals[:, 0]
        pos[pending, 1] = ctx.slm_center_m[1] + sigma_r[pending] * proposals[:, 1]
        pos[pending, 2] = sigma_z[pending] * proposals[:, 2]
        vel[pending, :] = sigma_v * proposals[:, 3:]
        attempts[pending] += 1
        energy = slm_site_energies(pos[pending], vel[pending], ctx, factors[pending])
        pending = pending[(energy >= 0) | ~np.isfinite(energy)]
        if not len(pending):
            break
    if len(pending):
        raise ToolkitError(E_BAD_VALUE,
                           f"could not prepare bound atoms at {len(pending)} sites "
                           f"(temperature too high for the shallowest site)",
                           "preparation.temperature_uK")
    return {"pos_m": pos, "vel_m_per_s": vel,
            "initial_energy_j": slm_site_energies(pos, vel, ctx, factors),
            "attempts_per_site": attempts,
            "sampler": "site_bound_harmonic_3d_v1"}
