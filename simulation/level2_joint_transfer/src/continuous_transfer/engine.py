"""Continuous 6D velocity-Verlet and finite-duration SU(2) propagation.

Numba is an optional acceleration, not a different physics model. Its scalar
field kernel is regression checked against Level 3's vector implementation.
Each original shot owns its complete history: no resampling between moves.
"""
import math
import numpy as np
try:
    from numba import njit, prange
    ACCELERATED = True
except ImportError:
    ACCELERATED = False
    prange = range
    def njit(*args, **kwargs):
        return lambda f: f


@njit(cache=True)
def field(x, y, z, d, w, zr, s1, s2):
    a = 1 + ((z - s1) / zr) ** 2
    b = 1 + ((z - s2) / zr) ** 2
    w2 = w * w
    u = -d / math.sqrt(a * b) * math.exp(-2 * (x*x/a + y*y/b) / w2)
    fx = u * 4 * x / (w2 * a)
    fy = u * 4 * y / (w2 * b)
    fz = -u * ((z-s1)/(zr*zr) * (-1/a + 4*x*x/(w2*a*a)) +
               (z-s2)/(zr*zr) * (-1/b + 4*y*y/(w2*b*b)))
    return u, fx, fy, fz


@njit(cache=True)
def total_field(x, y, z, control, alignment, ds, da_factor, ws, wa, zrs, zra, vs):
    s1 = 2 * zra * control[3] / vs if vs > 0 else 0.0
    s2 = 2 * zra * control[4] / vs if vs > 0 else 0.0
    ua, ax, ay, az = field(x-control[1]-alignment, y-control[2], z,
                          control[5]*da_factor, wa, zra, s1, s2)
    us, sx, sy, sz = field(x, y, z, ds*control[6], ws, zrs, 0.0, 0.0)
    return us+ua, sx+ax, sy+ay, sz+az


@njit(cache=True)
def spin_step(a, b, delta, omega, phi, dt):
    """First row (a,b) of SU(2) [[a,b],[-b*,a*]], left multiplication."""
    q = math.sqrt(delta*delta + omega*omega)
    c = math.cos(q*dt/2)
    s = math.sin(q*dt/2)/q if q > 1e-30 else dt/2
    aa = c - 1j*s*delta
    bb = -1j*s*omega * complex(math.cos(phi), -math.sin(phi))
    return aa*a - bb*b.conjugate(), aa*b + bb*a.conjugate()


@njit(cache=True, parallel=True)
def integrate(pos, vel, controls, segments, rounds, mass, ds_nom, ws, wa_nom,
              slm_wavelength, aod_wavelength, disorder, jitter, power, vs, eta):
    shots = len(pos)
    count = 2*rounds + 1
    alive = np.zeros((shots, count), dtype=np.bool_)
    phase = np.full((shots, count), np.nan)
    spin = np.full((shots, count, 2), complex(np.nan, np.nan), dtype=np.complex128)
    states = np.full((shots, count, 6), np.nan)
    energies = np.full((shots, count), np.nan)
    elapsed = np.zeros(shots)
    injected = np.zeros(shots)
    lost_stage = np.full(shots, -1, dtype=np.int64)
    lost_round = np.full(shots, -1, dtype=np.int64)
    final = np.zeros((shots, 6))
    for shot in prange(shots):
        x,y,z = pos[shot]
        vx,vy,vz = vel[shot]
        ds = ds_nom * disorder[shot, 0]
        wa = wa_nom * disorder[shot, 2]
        zrs, zra = math.pi*ws*ws/slm_wavelength, math.pi*wa*wa/aod_wavelength
        a, b, pending, phase0 = 1+0j, 0+0j, 0.0, 0.0
        alive[shot,0], phase[shot,0] = True, 0.0
        spin[shot,0,0], spin[shot,0,1] = a,b
        states[shot,0,:3], states[shot,0,3:] = pos[shot], vel[shot]
        dead = False
        for r in range(rounds):
            for leg in range(len(segments)):
                start, end, stream, judge = segments[leg]
                alignment = disorder[shot,3] + jitter[r,stream,shot]
                u,fx,fy,fz = total_field(x,y,z,controls[start],alignment,ds,disorder[shot,1],ws,wa,zrs,zra,vs)
                for i in range(start, end):
                    dt = controls[i+1,0] - controls[i,0]
                    vx += .5*dt*fx/mass
                    vy += .5*dt*fy/mass
                    vz += .5*dt*fz/mass
                    x += dt*vx
                    y += dt*vy
                    z += dt*vz
                    un,fx,fy,fz = total_field(x,y,z,controls[i+1],alignment,ds,disorder[shot,1],ws,wa,zrs,zra,vs)
                    vx += .5*dt*fx/mass
                    vy += .5*dt*fy/mass
                    vz += .5*dt*fz/mass
                    if power > 0:
                        speed2 = vx*vx + vy*vy + vz*vz
                        if speed2 > 1e-300:
                            scale = math.sqrt(1 + 2*power*dt/(mass*speed2))
                            vx,vy,vz = vx*scale,vy*scale,vz*scale
                        else:
                            vx = math.sqrt(2*power*dt/mass)
                        injected[shot] += power*dt
                    delta = eta * .5*(u+un) / 1.054571817e-34
                    phase0 += delta*dt
                    om = controls[i,7]
                    if om != 0:
                        if pending != 0:
                            rotation = complex(math.cos(pending/2), -math.sin(pending/2))
                            a,b = rotation*a, rotation*b
                            pending = 0.0
                        a,b = spin_step(a,b,delta,om,controls[i,8],dt)
                    else:
                        pending += delta*dt
                    u = un
                    elapsed[shot] += dt
                if judge:
                    if judge == 2:
                        bound = field(x,y,z,ds,ws,zrs,0.0,0.0)[0]
                    else:
                        bound = field(x-controls[end,1]-alignment,y-controls[end,2],z,
                                      controls[end,5]*disorder[shot,1],wa,zra,0.0,0.0)[0]
                    energy = .5*mass*(vx*vx+vy*vy+vz*vz) + bound
                    cp = 2*r + (1 if judge == 1 else 2)
                    energies[shot,cp] = energy
                    if not math.isfinite(energy) or energy >= 0:
                        dead = True
                        lost_stage[shot], lost_round[shot] = leg,r
                        break
                    if judge != 3:
                        rotation = complex(math.cos(pending/2), -math.sin(pending/2))
                        a,b = rotation*a, rotation*b
                        pending = 0.0
                        alive[shot,cp], phase[shot,cp] = True,phase0
                        spin[shot,cp,0], spin[shot,cp,1] = a,b
                        states[shot,cp] = np.array([x,y,z,vx,vy,vz])
            if dead:
                break
        final[shot] = np.array([x,y,z,vx,vy,vz])
    return alive, phase, spin, states, energies, elapsed, injected, lost_stage, lost_round, final


def run_continuous(pos, vel, protocol, *, rounds, physics, draws, jitter_sigma_m,
                   power_w, vs_m_s, eta=1.3e-4):
    if rounds < 1 or len(pos) == 0 or np.shape(pos) != np.shape(vel) or np.shape(pos)[1] != 3:
        raise ValueError("Expected nonempty (shots,3) states and positive rounds")
    disorder = np.column_stack((draws.slm_depth_factor, draws.aod_depth_factor,
                                draws.aod_waist_factor, draws.alignment_offset_m))
    if np.any(disorder[:,:3] <= 0) or power_w < 0 or jitter_sigma_m < 0:
        raise ValueError("Nonphysical depths, waists or noise")
    # Six independent streams indexed by ORIGINAL shot ID, never compacted pool.
    # The same PU/wait/DO samples apply in both protocol arms and all waveforms.
    jitter = np.stack([np.random.default_rng(np.random.SeedSequence([draws.jitter_seed,k]))
                       .normal(0,jitter_sigma_m,(rounds,len(pos))) for k in range(6)], axis=1)
    values = integrate(np.asarray(pos), np.asarray(vel), protocol.controls, protocol.segments,
                       rounds, physics["mass_kg"], physics["slm_depth_j"], physics["slm_waist_m"],
                       physics["aod_waist_m"],1061e-9,1055e-9,disorder,jitter,
                       power_w,vs_m_s,eta)
    keys = ("alive","phase_rad","spin_ab","states","energy_j","elapsed_s","injected_j",
            "lost_stage","lost_round","final_states")
    return dict(zip(keys,values))
