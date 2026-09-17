"""Phase-mixed harmonic benchmark; not fitted to survival data."""
from run_study import *
from numba import njit,prange

@njit(cache=True,parallel=True)
def harmonic(shots,steps,dt,omega,rin,pointing,seed):
    out=np.empty(shots)
    for k in prange(shots):
        np.random.seed((seed+104729*k)%4294967291)
        phase=2*np.pi*np.random.random()
        x=np.sqrt(2.)/omega*np.cos(phase);v=np.sqrt(2.)*np.sin(phase)
        for i in range(steps):
            v-=.5*dt*omega**2*x;x+=dt*v;v-=.5*dt*omega**2*x
            v+=-omega**2*x*np.sqrt(rin*dt/2)*np.random.randn()
            v+=omega**2*np.sqrt(pointing*dt/2)*np.random.randn()
        out[k]=.5*(v*v+omega**2*x*x)
    return out

def main():
    numba.set_num_threads(8)
    rows=[];w=2*np.pi*1000;t=.02;rin=1e-6;gamma=w*w*rin/4
    for kind in ('parametric','pointing'):
        for dt in (1e-5,5e-6,2.5e-6):
            # Choose Pt=0.2; harmonic mass and E0 are both one in this benchmark.
            sx=4*.2/t/w**4
            e=harmonic(30000,round(t/dt),dt,w,rin if kind=='parametric' else 0,
                       sx if kind=='pointing' else 0,26091631)
            mean=np.exp(gamma*t) if kind=='parametric' else 1.2
            var=np.exp(3*gamma*t)-np.exp(2*gamma*t) if kind=='parametric' else .44
            row=dict(kind=kind,dt_s=dt,mean=float(e.mean()),expected_mean=mean,
                     variance=float(e.var()),expected_variance=var,
                     mean_se=float(e.std()/np.sqrt(len(e))),shots=len(e),seed=26091631)
            row['mean_z']=(row['mean']-mean)/row['mean_se']
            rows.append(row);print(row,flush=True)
            assert abs(row['mean_z'])<5
            assert abs(row['variance']/var-1)<.08
    write(HERE/'noise_benchmarks.json',rows)

if __name__=='__main__':main()
