"""Exercise the actual production diffusion kernel in its harmonic limit."""
from run_study import *
from continuous_transfer.protocol import ContinuousProtocol
from level2c_pickup_survival.survival_engine import DisorderDraws

def main():
    numba.set_num_threads(8)
    n=16000;t=.002;dt=.05e-6
    mass=132.90545196*1.6605390666e-27;depth=140*1.380649e-29;w=1.17e-6;zr=np.pi*w*w/1061e-9
    omega=np.array([np.sqrt(4*depth/(mass*w*w))]*2+[np.sqrt(2*depth/(mass*zr*zr))])
    e0=.02*1.380649e-29
    phase=np.random.default_rng(26091631).uniform(0,2*np.pi,(n,3))
    pos=np.sqrt(2*e0/mass)/omega*np.cos(phase);vel=np.sqrt(2*e0/mass)*np.sin(phase)
    controls=np.zeros((round(t/dt)+1,9));controls[:,0]=np.linspace(0,t,len(controls));controls[:,6]=1
    protocol=ContinuousProtocol(controls,np.array([[0,len(controls)-1,0,2]]),['static'],[],t,0,'matched','none',(t/2,t))
    result=run_noisy_survival(pos,vel,protocol,rounds=1,physics=dict(mass_kg=mass,slm_depth_j=depth,slm_waist_m=w,aod_waist_m=w),
        draws=DisorderDraws.disabled(n),jitter_sigma_m=0,vs_m_s=0,rin=1e-8,pointing=0,seed=26091631)
    state=result['states'][:,-1]
    e=.5*mass*(state[:,3:]**2+omega**2*state[:,:3]**2)/e0
    gamma=omega**2*1e-8/4
    rows=[]
    for j in range(3):
        expected=np.exp(gamma[j]*t);variance=np.exp(3*gamma[j]*t)-np.exp(2*gamma[j]*t)
        se=e[:,j].std()/np.sqrt(n)
        row=dict(axis=j,frequency_hz=omega[j]/(2*np.pi),mean=e[:,j].mean(),expected_mean=expected,
                 mean_z=(e[:,j].mean()-expected)/se,variance=e[:,j].var(),expected_variance=variance)
        assert abs(row['mean_z'])<5 and abs(row['variance']/variance-1)<.1
        rows.append(row)
    ledger=result['noise_ledger']
    rows.append(dict(mean_realized_work=float(ledger[:,0].mean()),mean_expected_work=float(ledger[:,1].mean()),
                     work_difference_se=float((ledger[:,0]-ledger[:,1]).std()/np.sqrt(n))))
    write(HERE/'gaussian_benchmark.json',rows);print(json.dumps(clean(rows)),flush=True)

if __name__=='__main__':main()
