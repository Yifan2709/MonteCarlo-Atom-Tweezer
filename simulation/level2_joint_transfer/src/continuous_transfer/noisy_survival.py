"""Classical full-Gaussian force diffusion, one-sided PSD per Hz.

C(t)=S delta(t)/2, so a noisy beam contributes dv=F/m sqrt(S/2)dW.
The same scalar intensity increment acts on all three axes of that beam;
the 1061/1055 nm beams are assumed independent. This is a *model* of flat
weak noise, not a measured Cs noise spectrum. Pointing is along the transport
x axis, with dv=-dF/dx /m sqrt(Sx/2)dW, independently for each beam.
No harmonic energy estimator, reference temperature, or transport heating
term enters these kicks. Noise step and mechanical step converge separately.
"""
import math
import numpy as np
from .engine import njit, prange, field, total_field


@njit(cache=True)
def beam_noise(x,y,z,d,w,zr,s1,s2):
    u,fx,fy,fz = field(x,y,z,d,w,zr,s1,s2)
    a = 1+((z-s1)/zr)**2
    # dF/dx, evaluated at the actual atom position (including mixed x-z term).
    l = -4*x/(w*w*a)
    gx = fx*l + 4*u/(w*w*a)
    gy = fy*l
    gz = fz*l - 8*u*x*(z-s1)/(w*w*zr*zr*a*a)
    return u,fx,fy,fz,gx,gy,gz


@njit(cache=True, parallel=True)
def propagate(pos,vel,controls,segments,rounds,mass,ds_nom,ws,wa_nom,disorder,
              jitter,vs,mode,power,rin,sx,noise_stride,seed,recoil_per_joule,
              slm_rin,slm_sx,absorb_at_checkpoints):
    shots = len(pos)
    alive = np.zeros((shots,2*rounds+1),np.bool_)
    states = np.full((shots,2*rounds+1,6),np.nan)
    energy = np.full((shots,2*rounds+1),np.nan)
    stage_energy = np.full((shots,rounds,len(segments)),np.nan)
    switch_work = np.full((shots,rounds,len(segments)),np.nan)
    injected = np.zeros((shots,3)) # realized noise work, expected noise work, elapsed
    lost = np.full((shots,2),-1,np.int64)
    for k in prange(shots):
        np.random.seed((seed+104729*k) % 4294967291)
        x,y,z = pos[k]
        vx,vy,vz = vel[k]
        ds,da,wa = ds_nom*disorder[k,0],disorder[k,1],wa_nom*disorder[k,2]
        zrs,zra = math.pi*ws*ws/1061e-9,math.pi*wa*wa/1055e-9
        alive[k,0] = True
        states[k,0,:3],states[k,0,3:] = pos[k],vel[k]
        prev_u = field(x,y,z,ds,ws,zrs,0.,0.)[0]
        energy[k,0] = .5*mass*(vx*vx+vy*vy+vz*vz)+prev_u
        dead = False
        for r in range(rounds):
            for leg in range(len(segments)):
                start,end,stream,judge = segments[leg]
                align = disorder[k,3]+jitter[r,stream,k]
                u,fx,fy,fz = total_field(x,y,z,controls[start],align,ds,da,ws,wa,zrs,zra,vs)
                switch_work[k,r,leg] = u-prev_u
                accumulated_dt = 0.
                for i in range(start,end):
                    dt = controls[i+1,0]-controls[i,0]
                    vx += .5*dt*fx/mass
                    vy += .5*dt*fy/mass
                    vz += .5*dt*fz/mass
                    x += dt*vx
                    y += dt*vy
                    z += dt*vz
                    c = controls[i+1]
                    u,fx,fy,fz = total_field(x,y,z,c,align,ds,da,ws,wa,zrs,zra,vs)
                    vx += .5*dt*fx/mass
                    vy += .5*dt*fy/mass
                    vz += .5*dt*fz/mass
                    accumulated_dt += dt
                    if mode == 0 and power > 0:
                        k0 = .5*mass*(vx*vx+vy*vy+vz*vz)
                        if k0 > 1e-300:
                            fac = math.sqrt(1+power*dt/k0)
                            vx,vy,vz = vx*fac,vy*fac,vz*fac
                        else:
                            vx = math.sqrt(2*power*dt/mass)
                        injected[k,0] += power*dt
                        injected[k,1] += power*dt
                    elif mode > 0 and ((i-start+1)%noise_stride == 0 or i == end-1):
                        h = accumulated_dt
                        accumulated_dt = 0.
                        f1,f2 = 0.,0.
                        if vs > 0:
                            f1,f2 = 2*zra*c[3]/vs,2*zra*c[4]/vs
                        ua,ax,ay,az,agx,agy,agz = beam_noise(x-c[1]-align,y-c[2],z,c[5]*da,wa,zra,f1,f2)
                        us,bx,by,bz,bgx,bgy,bgz = beam_noise(x,y,z,ds*c[6],ws,zrs,0.,0.)
                        # Isotropic recoil diffusion: sum_i m Var(dv_i)/2 = P dt.
                        rec = max(0.,-ua-us)*recoil_per_joule
                        expected = h*(rin*(ax*ax+ay*ay+az*az)+slm_rin*(bx*bx+by*by+bz*bz)
                                      + sx*(agx*agx+agy*agy+agz*agz)
                                      + slm_sx*(bgx*bgx+bgy*bgy+bgz*bgz))/(4*mass)+rec*h
                        k0 = .5*mass*(vx*vx+vy*vy+vz*vz)
                        if mode == 1:
                            ra = math.sqrt(rin*h/2)*np.random.randn()/mass
                            rb = math.sqrt(slm_rin*h/2)*np.random.randn()/mass
                            pa = math.sqrt(sx*h/2)*np.random.randn()/mass
                            pb = math.sqrt(slm_sx*h/2)*np.random.randn()/mass
                            rr = math.sqrt(2*rec*h/(3*mass))
                            vx += ax*ra+bx*rb-agx*pa-bgx*pb+rr*np.random.randn()
                            vy += ay*ra+by*rb-agy*pa-bgy*pb+rr*np.random.randn()
                            vz += az*ra+bz*rb-agz*pa-bgz*pb+rr*np.random.randn()
                        else:
                            # Diagnostic only: same local mean drift, suppress diffusion.
                            if k0 > 1e-300:
                                fac = math.sqrt(1+expected/k0)
                                vx,vy,vz = vx*fac,vy*fac,vz*fac
                            else:
                                vx = math.sqrt(2*expected/mass)
                        injected[k,0] += .5*mass*(vx*vx+vy*vy+vz*vz)-k0
                        injected[k,1] += expected
                    injected[k,2] += dt
                prev_u = u
                stage_energy[k,r,leg] = .5*mass*(vx*vx+vy*vy+vz*vz)+u
                if judge:
                    if judge == 2:
                        bound = field(x,y,z,ds,ws,zrs,0.,0.)[0]
                    else:
                        bound = field(x-controls[end,1]-align,y-controls[end,2],z,
                                      controls[end,5]*da,wa,zra,0.,0.)[0]
                    cp = 2*r+(1 if judge == 1 else 2)
                    e = .5*mass*(vx*vx+vy*vy+vz*vz)+bound
                    energy[k,cp] = e
                    states[k,cp] = np.array([x,y,z,vx,vy,vz])
                    if not math.isfinite(e) or e >= 0:
                        if lost[k,0] < 0:
                            lost[k,0],lost[k,1] = r,leg
                        if absorb_at_checkpoints or not math.isfinite(e):
                            dead = True
                            break
                        continue
                    alive[k,cp] = True
            if dead:
                break
    return alive,states,energy,stage_energy,switch_work,injected,lost


def run_noisy_survival(pos,vel,protocol,*,rounds,physics,draws,jitter_sigma_m,
                       vs_m_s,mode="diffusion",power_w=0.,rin=0.,pointing=0.,
                       noise_step_s=0.05e-6,seed=26091611,recoil_per_joule=0.,
                       slm_rin=None,slm_pointing=None,absorb_at_checkpoints=True):
    if mode not in ("constant","diffusion","local_mean"):
        raise ValueError(mode)
    if rounds < 1 or len(pos)==0 or np.shape(pos)!=np.shape(vel) or np.shape(pos)[1:]!=(3,):
        raise ValueError("Expected nonempty (shots,3) states and positive rounds")
    if protocol.arm != "matched":
        raise ValueError("This survival-only kernel supports the Fig6d matched protocol")
    slm_rin = rin if slm_rin is None else slm_rin
    slm_pointing = pointing if slm_pointing is None else slm_pointing
    amplitudes=np.array([power_w,rin,pointing,jitter_sigma_m,recoil_per_joule,slm_rin,slm_pointing])
    if not np.all(np.isfinite(amplitudes)) or np.any(amplitudes<0):
        raise ValueError("Noise amplitudes must be nonnegative")
    if protocol.dd_mode != "none":
        raise ValueError("Survival-only protocol requires dd=False (no spin force)")
    steps = np.diff(protocol.controls[:,0])
    dt = np.median(steps[steps>0])
    stride = int(round(noise_step_s/dt))
    if stride < 1 or not np.isclose(stride*dt,noise_step_s,rtol=0,atol=1e-14):
        raise ValueError("Noise step must be an integer multiple of mechanical step")
    disorder=np.column_stack((draws.slm_depth_factor,draws.aod_depth_factor,
                              draws.aod_waist_factor,draws.alignment_offset_m))
    jitter=np.stack([np.random.default_rng(np.random.SeedSequence([draws.jitter_seed,k]))
                     .normal(0,jitter_sigma_m,(rounds,len(pos))) for k in range(6)],axis=1)
    values=propagate(np.asarray(pos),np.asarray(vel),protocol.controls,protocol.segments,rounds,
        physics['mass_kg'],physics['slm_depth_j'],physics['slm_waist_m'],physics['aod_waist_m'],
        disorder,jitter,vs_m_s,{"constant":0,"diffusion":1,"local_mean":2}[mode],power_w,
        rin,pointing,stride,seed,recoil_per_joule,slm_rin,slm_pointing,absorb_at_checkpoints)
    return dict(zip(('alive','states','energy_j','stage_energy_j','switch_work_j',
                     'noise_ledger','lost_round_stage'),values))
