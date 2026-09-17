"""Independent hold prediction: Fig3a no-cooling lifetime is ~132 s."""
from run_study import *
from continuous_transfer.protocol import ContinuousProtocol
from level2c_pickup_survival.survival_engine import DisorderDraws

def main():
    numba.set_num_threads(8)
    n=512;dt=.1e-6;chunk=.001
    _,physics,*_=case_inputs(SETTINGS,15,n,26091631,'matched')
    physics=dict(physics,slm_depth_j=180*1.380649e-29)
    draws=DisorderDraws.disabled(n)
    initial=sample_site_bound_harmonic(26091631,15,physics,draws.slm_depth_factor)
    c=np.zeros((round(chunk/dt)+1,9));c[:,0]=np.linspace(0,chunk,len(c));c[:,6]=1
    p=ContinuousProtocol(c,np.array([[0,len(c)-1,0,2]]),['static'],[],chunk,0,'matched','none',(chunk/2,chunk))
    result=run_noisy_survival(initial['pos_m'],initial['vel_m_per_s'],p,rounds=200,
        physics=physics,draws=draws,jitter_sigma_m=0,vs_m_s=0,rin=1e-8,pointing=0,
        seed=26091631,noise_step_s=.1e-6)
    s=result['alive'][:,::2].mean(axis=0)
    row=dict(seed=26091631,shots=n,temperature_uK=15,depth_uK=180,rin=1e-8,
             time_s=np.arange(201)*chunk,survival=s,observed_1e_lifetime_s=132,
             caveat='Independent paper operating condition, not a fitted likelihood; excludes recoil and other heating, so added heating cannot rescue this white spectrum.')
    write(HERE/'static_hold.json',row)
    print('S(0.2s)=',s[-1],flush=True)

if __name__=='__main__':main()
