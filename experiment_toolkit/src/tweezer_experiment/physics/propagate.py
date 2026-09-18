"""Velocity-Verlet propagation with optional force diffusion (F3 core).

Semantics extracted from the research kernel and kept identical:

- Force is evaluated at the left and right control of each step (segments
  own duplicated boundary nodes, so an SLM switch changes the potential
  without teleporting the atom).
- Binding judges run at checkpoints (segment ends carrying ``observe``):
  energy in the named single beam; E >= 0 or non-finite fails. With
  ``absorb=True`` the atom stays dead afterwards (absorbing semantics);
  with ``absorb=False`` only the first failure is recorded and propagation
  continues (recapture-allowed diagnosis).
- Diffusion kicks (model=diffusion) apply every ``noise_dt`` of accumulated
  time and at each segment end: one scalar RIN increment per beam shared by
  that beam's three force components, one pointing increment per beam
  coupling through dF/dx, and isotropic recoil with rate proportional to
  the local total depth. local_mean applies the expected energy of the same
  channels deterministically (diagnostic); legacy_constant_power is the
  historical fixed-power heating kept only as an explicitly named legacy
  mode.

Vectorization note: dead atoms continue to be integrated (cheap, masked
out of all records); their random draws keep the per-atom streams
independent, so surviving-atom statistics are unaffected.
"""
from __future__ import annotations

import numpy as np

from ..errors import ToolkitError, E_NOT_DIVISIBLE
from .field import BeamContext, beam_field, beam_x_gradient, beam_energy
from .initialization import slm_site_energies


def propagate(compiled, pos0, vel0, ctx: BeamContext, site, jitter,
              noise, noise_dt_s: float, seed: int, *, absorb: bool = True,
              save_states: str = "checkpoints") -> dict:
    controls = compiled.controls
    bounds = compiled.segment_bounds
    checkpoints = compiled.checkpoints
    dt = compiled.dt_s
    if noise.model in ("diffusion", "local_mean"):
        stride = int(round(noise_dt_s / dt))
        if stride < 1 or not np.isclose(stride * dt, noise_dt_s, rtol=0, atol=1e-15):
            raise ToolkitError(E_NOT_DIVISIBLE,
                               f"noise_dt {noise_dt_s * 1e6:g} us must be an integer multiple "
                               f"of dt {dt * 1e6:g} us", "numerics.noise_dt_us")
    else:
        stride = None

    shots = len(pos0)
    x, y, z = np.asarray(pos0, dtype=float).T.copy()
    vx, vy, vz = np.asarray(vel0, dtype=float).T.copy()
    mass = ctx.mass_kg
    slm_factor_site = np.asarray(site[0])
    aod_factor = np.asarray(site[1])
    align_site = np.asarray(site[3])
    alive_now = np.ones(shots, dtype=bool)
    n_ck = len(checkpoints)
    alive = np.zeros((shots, n_ck + 1), dtype=bool)
    alive[:, 0] = True
    ck_energy = np.full((shots, n_ck + 1), np.nan)
    keep_states = save_states in ("checkpoints", "all")
    ck_states = np.full((shots, n_ck + 1, 6), np.nan) if keep_states else None
    n_inst = len(bounds)
    stage_energy = np.full((shots, n_inst), np.nan)
    switch_work = np.full((shots, n_inst), np.nan)
    noise_ledger = np.zeros((shots, 2))          # realized, expected
    lost = np.full((shots, 2), -1, dtype=np.int64)
    rng = np.random.default_rng(seed)
    mode = {"none": 0, "legacy_constant_power": 0, "diffusion": 1, "local_mean": 2}[noise.model]
    power = noise.legacy_power_w or 0.0
    rin, sx = noise.aod_rin_psd_per_hz, noise.aod_pointing_psd_m2_per_hz
    slm_rin, slm_sx = noise.slm_rin_psd_per_hz, noise.slm_pointing_psd_m2_per_hz
    recoil_c = noise.recoil_per_joule

    def fields(ctrl, align_x):
        """Total field with per-atom disorder and alignment."""
        s1 = 2.0 * ctx.aod_zr_m * site[2] ** 2 * ctrl[3] / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
        s2 = 2.0 * ctx.aod_zr_m * site[2] ** 2 * ctrl[4] / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
        waist_a = ctx.aod_waist_m * site[2]
        zr_a = ctx.aod_zr_m * site[2] ** 2
        ua, ax, ay, az = beam_field(x - ctrl[1] - align_x, y - ctrl[2], z,
                                    ctrl[5] * aod_factor, waist_a, zr_a, s1, s2)
        us, bx, by, bz = beam_field(x - ctx.slm_center_m[0], y - ctx.slm_center_m[1], z,
                                    ctx.slm_depth_j * slm_factor_site * ctrl[6],
                                    ctx.slm_waist_m, ctx.slm_zr_m)
        return us + ua, bx + ax, by + ay, bz + az

    prev_u = beam_field(x - ctx.slm_center_m[0], y - ctx.slm_center_m[1], z,
                        ctx.slm_depth_j * slm_factor_site, ctx.slm_waist_m, ctx.slm_zr_m)[0]
    ck_index = {(label, idx): k + 1 for k, (label, idx, _trap) in enumerate(checkpoints)}
    # initial record (checkpoint index 0): the prepared state itself
    ck_energy[:, 0] = 0.5 * mass * (vx * vx + vy * vy + vz * vz) + prev_u
    if keep_states:
        ck_states[:, 0] = np.column_stack([x, y, z, vx, vy, vz])

    for inst, (i0, i1) in enumerate(bounds):
        align = align_site + jitter[inst]
        u0, fx, fy, fz = fields(controls[i0], align)
        switch_work[:, inst] = u0 - prev_u
        accumulated = 0.0
        for i in range(i0, i1):
            step = controls[i + 1, 0] - controls[i, 0]
            if step <= 0:
                _, fx, fy, fz = fields(controls[i + 1], align)
                continue
            vx += 0.5 * step * fx / mass
            vy += 0.5 * step * fy / mass
            vz += 0.5 * step * fz / mass
            x += step * vx
            y += step * vy
            z += step * vz
            u, fx, fy, fz = fields(controls[i + 1], align)
            vx += 0.5 * step * fx / mass
            vy += 0.5 * step * fy / mass
            vz += 0.5 * step * fz / mass
            accumulated += step
            if mode == 0 and power > 0.0:
                k0 = 0.5 * mass * (vx * vx + vy * vy + vz * vz)
                scale = np.sqrt(1.0 + power * step / np.maximum(k0, 1e-300))
                vx, vy, vz = vx * scale, vy * scale, vz * scale
                noise_ledger[:, 0] += power * step
                noise_ledger[:, 1] += power * step
            elif mode > 0 and ((i - i0 + 1) % stride == 0 or i == i1 - 1):
                h = accumulated
                accumulated = 0.0
                ctrl = controls[i + 1]
                s1 = 2.0 * ctx.aod_zr_m * site[2] ** 2 * ctrl[3] / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
                s2 = 2.0 * ctx.aod_zr_m * site[2] ** 2 * ctrl[4] / ctx.vs_m_s if ctx.vs_m_s > 0 else 0.0
                waist_a = ctx.aod_waist_m * site[2]
                zr_a = ctx.aod_zr_m * site[2] ** 2
                ua, ax, ay, az, agx, agy, agz = beam_x_gradient(
                    x - ctrl[1] - align, y - ctrl[2], z, ctrl[5] * aod_factor, waist_a, zr_a, s1, s2)
                us, bx, by, bz, bgx, bgy, bgz = beam_x_gradient(
                    x - ctx.slm_center_m[0], y - ctx.slm_center_m[1], z,
                    ctx.slm_depth_j * slm_factor_site * ctrl[6], ctx.slm_waist_m, ctx.slm_zr_m)
                rec = np.maximum(0.0, -(ua + us)) * recoil_c
                expected = h * (rin * (ax * ax + ay * ay + az * az)
                                + slm_rin * (bx * bx + by * by + bz * bz)
                                + sx * (agx * agx + agy * agy + agz * agz)
                                + slm_sx * (bgx * bgx + bgy * bgy + bgz * bgz)) / (4.0 * mass) \
                    + rec * h
                k0 = 0.5 * mass * (vx * vx + vy * vy + vz * vz)
                if mode == 1:
                    ra = np.sqrt(rin * h / 2.0) * rng.standard_normal(shots) / mass
                    rb = np.sqrt(slm_rin * h / 2.0) * rng.standard_normal(shots) / mass
                    pa = np.sqrt(sx * h / 2.0) * rng.standard_normal(shots) / mass
                    pb = np.sqrt(slm_sx * h / 2.0) * rng.standard_normal(shots) / mass
                    rr = np.sqrt(2.0 * rec * h / (3.0 * mass))
                    r3 = rng.standard_normal((shots, 3))
                    vx += ax * ra + bx * rb - agx * pa - bgx * pb + rr * r3[:, 0]
                    vy += ay * ra + by * rb - agy * pa - bgy * pb + rr * r3[:, 1]
                    vz += az * ra + bz * rb - agz * pa - bgz * pb + rr * r3[:, 2]
                else:  # local_mean: deterministic expected-energy kick
                    scale = np.sqrt(1.0 + expected / np.maximum(k0, 1e-300))
                    vx, vy, vz = vx * scale, vy * scale, vz * scale
                noise_ledger[:, 0] += 0.5 * mass * (vx * vx + vy * vy + vz * vz) - k0
                noise_ledger[:, 1] += expected
        stage_energy[:, inst] = 0.5 * mass * (vx * vx + vy * vy + vz * vz) + u
        prev_u = u
        # checkpoint judge for this instance, if any
        for label, idx, trap in checkpoints:
            if idx != inst:
                continue
            ctrl_end = controls[i1]
            if trap == "slm":
                depth = ctx.slm_depth_j * slm_factor_site * ctrl_end[6]
                energy = beam_energy(x, y, z, vx, vy, vz, mass, depth,
                                     ctx.slm_waist_m, ctx.slm_zr_m, ctx.slm_center_m[:2])
            else:
                waist_a = ctx.aod_waist_m * site[2]
                zr_a = ctx.aod_zr_m * site[2] ** 2
                energy = beam_energy(x, y, z, vx, vy, vz, mass,
                                     ctrl_end[5] * aod_factor, waist_a, zr_a,
                                     (ctrl_end[1] + align, ctrl_end[2]))
            k = ck_index[(label, idx)]
            first = alive_now & (lost[:, 0] < 0)
            ck_energy[first, k] = energy[first]
            if keep_states:
                ck_states[first, k] = np.column_stack([x, y, z, vx, vy, vz])[first]
            failed = first & (~np.isfinite(energy) | (energy >= 0))
            if np.any(failed):
                lost[failed, 0] = inst
                lost[failed, 1] = i1
                if absorb:
                    alive_now &= ~failed
            passed = alive_now & (lost[:, 0] < 0)
            alive[passed, k] = True

    return {"alive": alive, "checkpoint_energy_j": ck_energy,
            "checkpoint_states": ck_states, "stage_energy_j": stage_energy,
            "switch_work_j": switch_work, "noise_ledger_j": noise_ledger,
            "lost_instance": lost[:, 0], "lost_control_index": lost[:, 1],
            "final_states": np.column_stack([x, y, z, vx, vy, vz]),
            "initial_energy_j": slm_site_energies(np.asarray(pos0), np.asarray(vel0),
                                                  ctx, slm_factor_site)}
