"""Full production noise increments vs quadratic moment prediction and work ledger."""
import json
from pathlib import Path
import numpy as np
from long_distance_transport.model import Config,run,KB,metrics
ROOT=Path(__file__).resolve().parent
rows=[]
for dt in [.1,.05,.025]:
    c=Config(shots=16384,distance_um=0,duration_us=400,hold_us=100,
             temperature_uK=.05,rin=5.75e-9,dt_us=dt,noise_us=.1,seed=26091761)
    result,r=run(c);depth=c.depth_uK*KB*1e-6
    initial=r['energy_j'][:,0]+depth;end=r['energy_j'][:,-1]+depth
    wr=metrics(c)['omega_r_rad_s'];zr=np.pi*c.waist*c.waist/c.wavelength
    wz=np.sqrt(2*depth/(c.mass*zr*zr));t=500e-6
    expected=KB*c.temperature_uK*1e-6*(2*np.exp(wr*wr*c.rin*t/4)+np.exp(wz*wz*c.rin*t/4))
    se=end.std(ddof=1)/np.sqrt(c.shots)
    residual=(end-initial)-r['noise_work_j'][:,0]
    difference=r['noise_work_j'][:,0]-r['noise_work_j'][:,1]
    rows.append(dict(dt_us=dt,n=c.shots,mean_energy_J=float(end.mean()),harmonic_reference_J=float(expected),
        mean_z=float((end.mean()-expected)/se),
        noise_work_mean_z=float(difference.mean()/(difference.std(ddof=1)/np.sqrt(c.shots))),
        ledger_rms_over_initial_mean=float(np.sqrt(np.mean(residual**2))/initial.mean()),
        survival=result['observations'][-1]['survival']))
(ROOT/'production_noise_benchmark.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(rows,indent=2))
assert all(abs(r['mean_z'])<5 and abs(r['noise_work_mean_z'])<5 for r in rows)
