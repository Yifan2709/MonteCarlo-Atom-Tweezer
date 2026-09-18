"""Gaussian-beam fields, vectorized over atoms (F3 physics core).

Math is the closed-form 3D Gaussian tweezer with optional axial focus
shifts s1, s2 per scanning axis (acousto lensing: s = 2*z_R*v_cmd/v_s),
extracted from the research engine and verified against it analytically in
tests. Wavelengths, waists and depths enter as parameters — there are no
hardcoded 1055/1061 nm constants in this module.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def beam_field(x, y, z, depth, waist, zr, s1=0.0, s2=0.0):
    """Potential energy (J) and force components (N) of one Gaussian beam.

    u = -D/sqrt(a b) exp(-2 (x^2/a + y^2/b)/w^2),  a = 1+((z-s1)/zr)^2 etc.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    a = 1.0 + ((z - s1) / zr) ** 2
    b = 1.0 + ((z - s2) / zr) ** 2
    w2 = waist * waist
    u = -depth / np.sqrt(a * b) * np.exp(-2.0 * (x * x / a + y * y / b) / w2)
    fx = u * 4.0 * x / (w2 * a)
    fy = u * 4.0 * y / (w2 * b)
    fz = -u * ((z - s1) / (zr * zr) * (-1.0 / a + 4.0 * x * x / (w2 * a * a))
               + (z - s2) / (zr * zr) * (-1.0 / b + 4.0 * y * y / (w2 * b * b)))
    return u, fx, fy, fz


def beam_x_gradient(x, y, z, depth, waist, zr, s1=0.0, s2=0.0):
    """Field plus d(Fx,Fy,Fz)/dx for one beam (pointing-noise coupling).

    A beam-position displacement delta_x along x changes the force by
    (dFx/dx, dFy/dx, dFz/dx) * delta_x; these are the g-components used by
    the diffusion kicks. Extracted verbatim from the research kernel.
    """
    u, fx, fy, fz = beam_field(x, y, z, depth, waist, zr, s1, s2)
    a = 1.0 + ((np.asarray(z) - s1) / zr) ** 2
    w2 = waist * waist
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    ell = -4.0 * x / (w2 * a)
    gx = fx * ell + 4.0 * u / (w2 * a)
    gy = fy * ell
    gz = fz * ell - 8.0 * u * x * (z - s1) / (w2 * zr * zr * a * a)
    return u, fx, fy, fz, gx, gy, gz


@dataclass(frozen=True)
class BeamContext:
    """Per-run beam parameters (SI), including per-site disorder factors."""
    mass_kg: float
    slm_waist_m: float
    slm_zr_m: float
    slm_depth_j: float            # nominal (before site factor)
    slm_center_m: tuple[float, float]
    aod_waist_m: float            # nominal (before site factor)
    aod_zr_m: float
    vs_m_s: float                 # 0 disables lensing shifts

    @staticmethod
    def from_device(device, aod_waist_m=None) -> "BeamContext":
        slm, aod = device.slm, device.aod
        zr_s = np.pi * slm.waist_m ** 2 / slm.wavelength_m
        waist_a = aod_waist_m if aod_waist_m is not None else aod.waist_m
        zr_a = np.pi * waist_a ** 2 / aod.wavelength_m
        vs = device.lensing_vs_m_s if device.lensing_mode == "velocity_shift" else 0.0
        return BeamContext(device.mass_kg, slm.waist_m, zr_s, slm.nominal_depth_j,
                           (slm.center_m[0], slm.center_m[1]), waist_a, zr_a, vs)


def total_field(x, y, z, ctrl, align_x, ctx: BeamContext, site):
    """Sum of SLM and AOD beams at one control row.

    ``ctrl`` = (t, cx, cy, vx, vy, aod_depth_j, slm_factor).
    ``site`` = (slm_depth_factor, aod_depth_factor, aod_waist_factor) with
    array shape (shots,). ``align_x`` is the per-shot AOD x misalignment (m).
    """
    _, cx, cy, vx, vy, depth_a, factor_s = ctrl
    s_factor = site[0]
    a_factor = site[1]
    waist_a = ctx.aod_waist_m * site[2]
    zr_a = ctx.aod_zr_m * site[2] ** 2          # zr scales with waist^2
    s1 = 2.0 * zr_a * vx / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
    s2 = 2.0 * zr_a * vy / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
    ua, ax, ay, az = beam_field(x - cx - align_x, y - cy, z, depth_a * a_factor,
                                waist_a, zr_a, s1, s2)
    us, bx, by, bz = beam_field(x - ctx.slm_center_m[0], y - ctx.slm_center_m[1], z,
                                ctx.slm_depth_j * s_factor * factor_s,
                                ctx.slm_waist_m, ctx.slm_zr_m)
    return us + ua, bx + ax, by + ay, bz + az


def total_potential(x, y, z, ctrl, align_x, ctx, site):
    return total_field(x, y, z, ctrl, align_x, ctx, site)[0]


def beam_energy(x, y, z, vx, vy, vz, mass, depth, waist, zr, center=(0.0, 0.0)):
    """Kinetic + single-beam potential (binding judge)."""
    u = beam_field(x - center[0], y - center[1], z, depth, waist, zr)[0]
    return 0.5 * mass * (vx * vx + vy * vy + vz * vz) + u
