"""Independent physics references, all numerical residuals retained."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from numba import njit,prange
from long_distance_transport.model import Config,controls,initial,integrate,run,KB
from continuous_transfer.engine import field
from continuous_transfer.noisy_survival import beam_noise
from paper_integration.exact_trajectories import polynomial

ROOT=Path(__file__).resolve().parent

def potential(q,d,w,zr,sh):
    # Independently written optical intensity, no production force invocation.
    wx=w*np.sqrt(1+((q[2]-sh[0])/zr)**2)
    wy=w*np.sqrt(1+((q[2]-sh[1])/zr)**2)
    return -d*w*w/(wx*wy)*np.exp(-2*(q[0]/wx)**2-2*(q[1]/wy)**2)

def force_reference(q,d,w,zr,sh):
    h=w*2e-5
    return np.array([-(potential(q+np.eye(3)[j]*h,d,w,zr,sh)-potential(q-np.eye(3)[j]*h,d,w,zr,sh))/(2*h) for j in range(3)])

@njit(parallel=True,cache=True)
def oscillator_noise(n,dt,tend,omega,rin,sx,seed):
    e=np.zeros(n)
    # units mass=1, initial E=1; random phases, no initial energy variance.
    for k in prange(n):
        np.random.seed(seed+k)
        ph=2*np.pi*np.random.rand()
        x=np.sqrt(2)/omega*np.cos(ph);v=np.sqrt(2)*np.sin(ph)
        for j in range(int(round(tend/dt))):
            v-=.5*dt*omega**2*x;x+=dt*v;v-=.5*dt*omega**2*x
            v+=-omega**2*x*np.sqrt(rin*dt/2)*np.random.randn()+omega**2*np.sqrt(sx*dt/2)*np.random.randn()
        e[k]=.5*(v*v+omega**2*x*x)
    return e

def main():
    out={};c=Config(shots=32,temperature_uK=5,distance_um=0,duration_us=300,hold_us=100)
    d=c.depth_uK*1e-6*KB;w=c.waist;zr=np.pi*w*w/c.wavelength
    rng=np.random.default_rng(1971)
    errors=[];noise_errors=[]
    for sh in [(0.,0.),(zr,0.),(zr,zr),(zr,-zr)]:
        for _ in range(40):
            q=rng.normal(size=3)*[w,w,zr]*.4
            actual=np.array(field(*q,d,w,zr,*sh)[1:])
            ref=force_reference(q,d,w,zr,sh)
            errors.append(np.linalg.norm(ref-actual)/(d/w))
            h=w*1e-5
            fd=(np.array(field(*(q+[h,0,0]),d,w,zr,*sh)[1:])-np.array(field(*(q-[h,0,0]),d,w,zr,*sh)[1:]))/(2*h)
            noise_errors.append(np.linalg.norm(fd-np.array(beam_noise(*q,d,w,zr,*sh)[4:]))/(d/w**2))
    out['force_derivative_max_scaled']=max(errors)
    out['pointing_derivative_max_scaled']=max(noise_errors)
    out['static']=[];out['galilean']=[]
    for dt in [.1,.05,.025]:
        cc=Config(**{**c.__dict__,'dt_us':dt})
        result,raw=run(cc)
        out['static'].append(dict(dt_us=dt,max_energy_drift_over_depth=float(np.max(abs(raw['energy_j'][:,-1]-raw['energy_j'][:,0]))/d)))
        a,segs,_=controls(cc);p,v,fac,off,_=initial(cc)
        # Entire protocol uniformly translated, including nominal hold; bootstrap v+V.
        V=.6;a[:,1]=V*a[:,0];a[:,3]=V
        values=integrate(p,v+[V,0,0],fac,off,a,segs,cc.mass,w,zr,zr,1,0.,0.,1,cc.seed)
        state=values[0][:,-1]
        relative=state.copy();relative[:,0]-=V*a[segs[1,1],0];relative[:,3]-=V
        out['galilean'].append(dict(dt_us=dt,max_state_difference=float(np.max(abs(relative-raw['states'][:,-1])))))
    out['harmonic']=[]
    mass=c.mass;omega=2*np.pi*40000
    for tr in ['rb_cubic','zero_jerk','min_jerk']:
        poly=polynomial(tr);T=217e-6;L=55e-6
        f=lambda t:L/T**2*poly.deriv(2)(t/T)
        spec=quad(lambda t:f(t)*np.cos(omega*t),0,T,epsabs=1e-10)[0]+1j*quad(lambda t:f(t)*np.sin(omega*t),0,T,epsabs=1e-10)[0]
        ref=.5*mass*abs(spec)**2
        sol=solve_ivp(lambda t,y:[y[1],-omega**2*y[0]-f(t)],(0,T),[0.,0.],rtol=2e-11,atol=1e-14,method='DOP853')
        q,v=sol.y[:,-1];e=.5*mass*(v*v+omega**2*q*q)
        out['harmonic'].append(dict(trajectory=tr,energy_fourier_J=ref,energy_DOP853_J=e,relative_error=abs(e/ref-1)))
    # Independent full nonlinear reference, finite-difference gradient + DOP853.
    out['gaussian_reference']=[]
    for tau in [0.,2.]:
        cc=Config(shots=4,distance_um=270,duration_us=600,tau_us=tau,dt_us=.05,temperature_uK=19)
        result,raw=run(cc);T=cc.duration_us*1e-6;L=cc.distance_um*1e-6;p=polynomial(cc.trajectory)
        def rhs(t,y):
            s=min(t/T,1.);cx=L*p(s);vv=L/T*p.deriv()(s) if t<T else 0.
            q=y[:3]-[cx,0,0]
            f=force_reference(q,d,w,zr,(tau*1e-6/w*zr*vv,0.))
            return np.r_[y[3:],f/cc.mass]
        errs=[]
        for k in range(4):
            sol=solve_ivp(rhs,(0,T+cc.hold_us*1e-6),np.r_[raw['initial_pos'][k],raw['initial_vel'][k]],method='DOP853',rtol=2e-9,atol=1e-12,max_step=1e-6)
            delta=sol.y[:,-1]-raw['states'][k,-1]
            errs.append(float(np.linalg.norm(delta[:3])/w+np.linalg.norm(delta[3:])/np.sqrt(d/cc.mass)))
        out['gaussian_reference'].append(dict(tau_us=tau,max_scaled_state_error=max(errs)))
    out['noise']=[];noise_samples={}
    for mode in ['rin','pointing']:
        om=2*np.pi*30e3;t=200e-6;rin=5e-9 if mode=='rin' else 0.;sx=2e-18 if mode=='pointing' else 0.
        gamma=om**2*rin/4;power=om**4*sx/4
        mean=np.exp(gamma*t) if rin else 1+power*t
        var=np.exp(3*gamma*t)-np.exp(2*gamma*t) if rin else 2*power*t+(power*t)**2
        for dt in [.1,.05,.025]:
            vals=oscillator_noise(16384,dt*1e-6,t,om,rin,sx,971)
            noise_samples[f'{mode}_dt_{dt}']=vals
            fourth=float(np.mean((vals-vals.mean())**4))
            variance_se=np.sqrt((fourth-float(vals.var())**2)/len(vals))
            out['noise'].append(dict(mode=mode,dt_us=dt,n=len(vals),mean=float(vals.mean()),variance=float(vals.var()),reference_mean=mean,reference_variance=var,mean_z=(float(vals.mean())-mean)/np.sqrt(var/len(vals)),variance_SE=variance_se,variance_z=(float(vals.var())-var)/variance_se))
    np.savez_compressed(ROOT/'noise_reference_energies.npz',**noise_samples)
    # Exact published waveform derivatives form an external test, not copied coefficients.
    import pandas as pd
    tab=pd.read_csv(ROOT/'zhang_source/dataset/extended_data_fig2c.csv')
    out['published_trajectory_max_error']={}
    for tr,tag in [('zero_jerk','zero_jerk'),('min_jerk','min_jerk')]:
        for n,name in [(2,'acceleration'),(3,'jerk')]:
            col=f'{name}_{tag}_L_T{n}'
            out['published_trajectory_max_error'][col]=float(np.max(abs(polynomial(tr).deriv(n)(tab.t_over_T)-tab[col])))
    # One 2-leg call versus chained 1-leg calls through identical controls.
    cc=Config(shots=32,repeats=2,duration_us=1000,tau_us=2)
    a,segs,_=controls(cc);p,v,fac,off,_=initial(cc)
    args=(cc.mass,w,zr,np.pi*w*w/cc.slm_wavelength)
    both=integrate(p,v,fac,off,a,segs,*args,2,0.,0.,1,0)
    one=integrate(p,v,fac,off,a,segs,*args,1,0.,0.,1,0)
    # Swap segment halves to start with the inbound control.
    mid=len(segs)//2;rev=np.vstack([segs[mid:],segs[:mid]])
    two=integrate(one[0][:,-1,:3],one[0][:,-1,3:],fac,off,a,rev,*args,1,0.,0.,1,0)
    good=both[3][:,-1]
    out['continuity_max_state_difference']=float(np.max(abs(both[0][good,-1]-two[0][good,-1])))
    out['passed']=(max(errors)<1e-7 and max(noise_errors)<1e-7 and
        max(x['relative_error'] for x in out['harmonic'])<1e-6 and
        max(x['max_energy_drift_over_depth'] for x in out['static'])<1e-4 and
        max(x['max_state_difference'] for x in out['galilean'])<1e-6 and
        max(x['max_scaled_state_error'] for x in out['gaussian_reference'])<.02 and
        max(abs(x['mean_z']) for x in out['noise'])<5 and
        max(abs(x['variance_z']) for x in out['noise'])<5 and
        max(out['published_trajectory_max_error'].values())<1e-10 and
        out['continuity_max_state_difference']<1e-12)
    (ROOT/'benchmarks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2),flush=True)
    if not out['passed']:raise RuntimeError('Physics benchmark failed')

if __name__=='__main__':main()
