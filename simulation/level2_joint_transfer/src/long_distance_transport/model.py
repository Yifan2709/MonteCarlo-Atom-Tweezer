"""SI dynamics. No empirical losses; observations only after complete rests.

Production reuses audited Gaussian/astigmatic field and its noise derivative.
Each atom carries its original site, RNG stream and 6D state across all moves.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from numba import njit, prange
from scipy.stats import beta
from continuous_transfer.engine import field
from continuous_transfer.noisy_survival import beam_noise
from paper_integration.exact_trajectories import polynomial

KB = 1.380649e-23
U = 1.66053906660e-27


@dataclass(frozen=True)
class Config:
    species: str = "Cs"
    mass: float = 132.90545196*U
    wavelength: float = 1055e-9
    slm_wavelength: float = 1061e-9
    waist: float = 1.11e-6
    depth_uK: float = 280.
    slm_depth_uK: float = 140.
    temperature_uK: float = 19.
    tau_us: float = 0.
    geometry: str = "x"
    trajectory: str = "rb_cubic"
    distance_um: float = 510.
    duration_us: float = 1000.
    protocol: str = "A"
    repeats: int = 1
    ramp_us: float = 200.
    hold_us: float = 100.
    disorder: float = 0.
    alignment_nm: float = 0.
    rin: float = 0.
    pointing: float = 0.
    dt_us: float = .1
    noise_us: float = .2
    shots: int = 512
    seed: int = 26091701


def interval(k, n):
    return (float(beta.ppf(.025,k,n-k+1)) if k else 0.,
            float(beta.ppf(.975,k+1,n-k)) if k<n else 1.)


def metrics(c):
    p = polynomial(c.trajectory)
    extrema = []
    for order in (1,2,3):
        d = p.deriv(order)
        roots = d.deriv().roots()
        loc = [0.,1.] + [r.real for r in roots if abs(r.imag)<1e-10 and 0<r.real<1]
        extrema.append(float(np.max(np.abs(d(loc)))))
    L,T = c.distance_um*1e-6,c.duration_us*1e-6
    axis = (1.,0.) if c.geometry=="x" else (2**-.5,2**-.5)
    return dict(average_m_s=L/T, peak_m_s=extrema[0]*L/T,
        peak_x_m_s=axis[0]*extrema[0]*L/T,peak_y_m_s=axis[1]*extrema[0]*L/T,
        peak_accel_m_s2=extrema[1]*L/T**2,interior_peak_jerk_m_s3=extrema[2]*L/T**3,
        acceleration_join_jump_m_s2=abs(float(p.deriv(2)(0)))*L/T**2,
        complete_oneway_us=T*1e6+c.hold_us+(2*c.ramp_us if c.protocol=="B" else 0),
        omega_r_rad_s=np.sqrt(4*c.depth_uK*1e-6*KB/(c.mass*c.waist**2)),
        maximum_static_force_acceleration=2*c.depth_uK*1e-6*KB*np.exp(-.5)/(c.mass*c.waist),
        tau_over_T=c.tau_us/c.duration_us)


def controls(c):
    """One outward and one return leg, reused without resetting physical state.

    Columns t,cx,cy,vx,vy,zsx,zsy,Da,Ds,sx,sy (SI). At stage boundaries
    duplicate times carry explicit left/right potentials; never integrate dt=0.
    """
    p = polynomial(c.trajectory)
    direction = np.array([1.,0.]) if c.geometry=="x" else np.array([1.,1.])/np.sqrt(2)
    signs = np.array([1.,-1.]) if c.geometry=="opposite" else np.ones(2)
    end = c.distance_um*1e-6*direction
    zr = np.pi*c.waist**2/c.wavelength
    da,ds = c.depth_uK*1e-6*KB,c.slm_depth_uK*1e-6*KB
    arr, segments, labels, offset, clock = [],[],[],0,0.
    for leg in range(2):
        start,stop = (np.zeros(2),end) if leg==0 else (end,np.zeros(2))
        kinds = ["pickup","move","dropoff","hold"] if c.protocol=="B" else ["move","hold"]
        for kind in kinds:
            length = (c.duration_us if kind=="move" else c.hold_us if kind=="hold" else c.ramp_us)*1e-6
            n = max(1,int(round(length/(c.dt_us*1e-6))))
            t = np.linspace(0,length,n+1)
            s = t/length
            a = np.zeros((n+1,11)); a[:,0] = clock+t
            if kind=="move":
                a[:,1:3] = start+p(s)[:,None]*(stop-start)
                a[:,3:5] = p.deriv()(s)[:,None]*(stop-start)/length
                a[:,5:7] = c.tau_us*1e-6/c.waist*zr*a[:,3:5]*signs
                a[:,7] = da
            else:
                a[:,1:3] = start if kind=="pickup" else stop
                a[:,9:11] = a[:,1:3]
                if kind=="pickup":
                    a[:,7],a[:,8] = da*s,ds*(1-s)
                elif kind=="dropoff":
                    a[:,7],a[:,8] = da*(1-s),ds*s
                else:
                    a[:,8 if c.protocol=="B" else 7] = ds if c.protocol=="B" else da
            arr.append(a)
            segments.append((offset,offset+n,int(kind=="hold"),leg))
            labels.append(kind)
            offset += n+1; clock += length
    return np.vstack(arr),np.asarray(segments,dtype=np.int64),labels


def initial(c):
    """Draw site disorder BEFORE thermal proposals, reject in same actual site."""
    rng = np.random.default_rng(np.random.SeedSequence([c.seed,0]))
    factors = 1+c.disorder*rng.normal(size=(c.shots,2))
    if np.any(factors<=0):
        raise ValueError("nonpositive site")
    offsets = c.alignment_nm*1e-9*rng.normal(size=c.shots)
    d = (c.slm_depth_uK if c.protocol=="B" else c.depth_uK)*1e-6*KB
    dep = d*factors[:,1 if c.protocol=="B" else 0]
    zr = np.pi*c.waist**2/(c.slm_wavelength if c.protocol=="B" else c.wavelength)
    thermal = c.temperature_uK*1e-6*KB
    scales = np.column_stack((np.sqrt(thermal*c.waist**2/(4*dep)),
                              np.sqrt(thermal*c.waist**2/(4*dep)),
                              np.sqrt(thermal*zr**2/(2*dep))))
    pos,vel = np.zeros((c.shots,3)),np.zeros((c.shots,3))
    pending=np.arange(c.shots); attempts=np.zeros(c.shots,int)
    for _ in range(1000):
        r = rng.normal(size=(len(pending),6))
        pos[pending] = r[:,:3]*scales[pending]
        vel[pending] = r[:,3:]*np.sqrt(thermal/c.mass)
        attempts[pending]+=1
        q=1+(pos[pending,2]/zr)**2
        e=.5*c.mass*np.sum(vel[pending]**2,axis=1)-dep[pending]/q*np.exp(-2*np.sum(pos[pending,:2]**2,axis=1)/(c.waist**2*q))
        pending=pending[e>=0]
        if not len(pending):break
    if len(pending):raise RuntimeError("bound preparation failed")
    if c.protocol=="A":pos[:,0]+=offsets
    return pos,vel,factors,offsets,attempts


@njit(cache=True)
def forces(x,y,z,a,fa,fs,offset,w,zra,zrs):
    ua,fx,fy,fz = field(x-a[1]-offset,y-a[2],z,a[7]*fa,w,zra,a[5],a[6])
    us,gx,gy,gz = field(x-a[9],y-a[10],z,a[8]*fs,w,zrs,0.,0.)
    return ua+us,fx+gx,fy+gy,fz+gz


@njit(cache=True,parallel=True)
def integrate(pos,vel,factors,offsets,ctrl,segs,m,w,zra,zrs,repeats,rin,sx,stride,seed):
    N=len(pos); per=len(segs)//2
    states=np.full((N,repeats+1,6),np.nan)
    energies=np.full((N,repeats+1),np.nan)
    stage_energy=np.full((N,repeats,per),np.nan)
    alive=np.zeros((N,repeats+1),np.bool_)
    ledger=np.zeros((N,2))
    for k in prange(N):
        np.random.seed((seed+104729*k)%4294967291)
        x,y,z=pos[k]; vx,vy,vz=vel[k]
        fa,fs=factors[k]; off=offsets[k]
        states[k,0,:3]=pos[k];states[k,0,3:]=vel[k];alive[k,0]=True
        pot,fx,fy,fz=forces(x,y,z,ctrl[0],fa,fs,off,w,zra,zrs)
        energies[k,0]=.5*m*(vx*vx+vy*vy+vz*vz)+pot
        for r in range(repeats):
            for st in range(per):
                first,last,judge,leg=segs[(r%2)*per+st]
                pot,fx,fy,fz=forces(x,y,z,ctrl[first],fa,fs,off,w,zra,zrs)
                h=0.
                for i in range(first,last):
                    dt=ctrl[i+1,0]-ctrl[i,0]
                    vx+=.5*dt*fx/m;vy+=.5*dt*fy/m;vz+=.5*dt*fz/m
                    x+=dt*vx;y+=dt*vy;z+=dt*vz
                    a=ctrl[i+1]
                    pot,fx,fy,fz=forces(x,y,z,a,fa,fs,off,w,zra,zrs)
                    vx+=.5*dt*fx/m;vy+=.5*dt*fy/m;vz+=.5*dt*fz/m
                    h+=dt
                    if (rin>0 or sx>0) and ((i-first+1)%stride==0 or i==last-1):
                        _,ax,ay,az,gx,gy,gz=beam_noise(x-a[1]-off,y-a[2],z,a[7]*fa,w,zra,a[5],a[6])
                        dr=math.sqrt(rin*h/2)/m*np.random.randn()
                        dx=math.sqrt(sx*h/2)/m*np.random.randn()
                        kinetic=.5*m*(vx*vx+vy*vy+vz*vz)
                        vx+=ax*dr-gx*dx;vy+=ay*dr-gy*dx;vz+=az*dr-gz*dx
                        ledger[k,0]+=.5*m*(vx*vx+vy*vy+vz*vz)-kinetic
                        ledger[k,1]+=h*(rin*(ax*ax+ay*ay+az*az)+sx*(gx*gx+gy*gy+gz*gz))/(4*m)
                        h=0.
                # Diagnostic relative energy. Never absorb using this during motion.
                a=ctrl[last]
                stage_energy[k,r,st]=.5*m*((vx-a[3])**2+(vy-a[4])**2+vz*vz)+pot
            states[k,r+1]=np.array([x,y,z,vx,vy,vz])
            e=.5*m*(vx*vx+vy*vy+vz*vz)+pot
            energies[k,r+1]=e
            alive[k,r+1]=math.isfinite(e) and e<0
            if not alive[k,r+1]:break
    return states,energies,stage_energy,alive,ledger


def run(c):
    if c.protocol not in ('A','B') or c.repeats<1:raise ValueError(c)
    ctrl,segs,labels=controls(c)
    pos,vel,factors,offsets,attempts=initial(c)
    stride=max(1,int(round(c.noise_us/c.dt_us)))
    if c.rin>0 and not np.isclose(stride*c.dt_us,c.noise_us):raise ValueError("noise step not divisible")
    vals=integrate(pos,vel,factors,offsets,ctrl,segs,c.mass,c.waist,
        np.pi*c.waist**2/c.wavelength,np.pi*c.waist**2/c.slm_wavelength,
        c.repeats,c.rin,c.pointing,stride,c.seed)
    raw=dict(zip(('states','energy_j','stage_energy_j','alive','noise_work_j'),vals))
    raw.update(initial_pos=pos,initial_vel=vel,site_factors=factors,offsets=offsets,attempts=attempts)
    rows=[]
    for r in range(1,c.repeats+1):
        k=int(raw['alive'][:,r].sum());lo,hi=interval(k,c.shots)
        en=raw['energy_j'][:,r]/KB*1e6
        bound=en[raw['alive'][:,r]]
        rows.append(dict(repeat=r,alive=k,survival=k/c.shots,lower=lo,upper=hi,
            conditional_loss=1-k/max(1,int(raw['alive'][:,r-1].sum())),
            median_energy_uK=float(np.median(bound)) if len(bound) else None,
            energy_p95_uK=float(np.quantile(bound,.95)) if len(bound) else None))
    return dict(config=asdict(c),kinematics=metrics(c),observations=rows,stage_labels=labels[:len(labels)//2]),raw
