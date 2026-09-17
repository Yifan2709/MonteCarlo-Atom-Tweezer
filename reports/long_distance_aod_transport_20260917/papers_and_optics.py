"""Extract vector markers, compare external data, test finite acoustic propagation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pdfplumber
from scipy.integrate import quad
from scipy.special import erf
from long_distance_transport.model import KB,U
from paper_integration.exact_trajectories import polynomial
ROOT=Path(__file__).resolve().parent

def main():
    page=pdfplumber.open('tmp/pdfs/paper_relation_20260915/bluvstein_2022.pdf').pages[16]
    markers=[]
    for c in page.curves:
        if c.get('non_stroking_color')==(0.,.625,.492,0.) and c['fill'] and 221<c['x0']<317:
            # Exclude legend point (x~224.4,y~140.3); retain all 24 experimental markers.
            x=(c['x0']+c['x1'])/2;y=(c['top']+c['bottom'])/2
            if 137<y<143 and x<230:continue
            split=.2+(x-221.402837155)*.8/(309.335726375-221.402837155)
            surv=(156.080486475-y)/(156.080486475-64.00503095)
            markers.append(dict(x_pdf_pt=x,y_pdf_pt=y,separation_speed_m_s=split,
                atom_average_speed_m_s=split/2,duration_us=110/split,survival=surv,
                extraction_error_survival=.05/(156.080486475-64.00503095),
                experimental_standard_error=np.nan))
    tab=pd.DataFrame(markers).sort_values('separation_speed_m_s')
    tab.to_csv(ROOT/'rb_ed2a_digitized.csv',index=False)
    (ROOT/'rb_digitization.json').write_text(json.dumps(dict(pdf='tmp/pdfs/paper_relation_20260915/bluvstein_2022.pdf',page_1based=17,
        method='PDF vector filled CMYK markers; panel-a tick affine calibration',
        points=len(tab),coordinate_uncertainty_pt=.05,experimental_error='not recoverable from rendered figure; not equal to digitization uncertainty',
        axes={'x_ticks_pt':[221.402837155,309.335726375],'x_values':[.2,1.],
              'y_ticks_pt':[156.080486475,64.00503095],'y_values':[0.,1.]},
        background='paper subtracts 0.007 background loss; reproduced markers are already corrected'),indent=2))
    # Harmonic exact spectrum vs asymptote at the same frequency distribution.
    mass=86.909180531*U;hbar=1.054571817e-34;D=55e-6
    rows=[]
    for T in np.geomspace(15e-6,600e-6,80):
        omegas=np.linspace(.85,1.15,401)*2*np.pi*40000
        k=omegas*T
        # Integral of cubic acceleration: 6D/T * integral_0^1(1-2s)e^(iks)ds.
        inte=(1+np.exp(1j*k))*1j/k+2*(1-np.exp(1j*k))/k**2
        exact=mass/(2*hbar*omegas)*abs(6*D/T*inte)**2
        asym=36*mass*D*D/(hbar*omegas**3*T**4)
        rows.append(dict(duration_us=T*1e6,omegaT_min=k.min(),exact_mean_N=float(exact.mean()),asym_mean_N=float(asym.mean()),
                         relative_difference=float(exact.mean()/asym.mean()-1)))
    pd.DataFrame(rows).to_csv(ROOT/'rb_harmonic_averaging.csv',index=False)
    # Pupil propagation: acoustic phase phi(u,t) = 2/(w*tau) int_0^(tau*u)[x(t)-x(t-r)]dr.
    # Quadratic coefficient tau*v/w gives Eq.S1. Higher orders are NOT free loss terms.
    # Compare complex pupil fields using full phase vs instantaneous quadratic lens.
    nodes,weights=np.polynomial.legendre.leggauss(40)
    u=np.linspace(-3,3,1201);weight=np.exp(-2*u*u);norm=np.trapezoid(weight,u)
    result=[]
    for tr in ['rb_cubic','zero_jerk','min_jerk']:
        p=polynomial(tr)
        for T_us in [100.,400.,1000.,1600.,3000.]:
            T=T_us*1e-6;L=510e-6;w=1.11e-6;tau=5e-6
            overlaps=[];maxr=[]
            for t in np.linspace(-3*tau,T+3*tau,101):
                r=tau*u[:,None]*(nodes+1)[None,:]/2
                integral=tau*u/2*np.sum(weights[None,:]*(L*p(np.clip(t/T,0,1))-L*p(np.clip((t-r)/T,0,1))),axis=1)
                full=2*integral/(w*tau)
                v=L/T*p.deriv()(t/T) if 0<t<T else 0.
                remainder=full-tau*v/w*u*u
                overlaps.append(abs(np.trapezoid(weight*np.exp(1j*remainder),u)/norm)**2)
                maxr.append(np.sqrt(np.trapezoid(weight*remainder**2,u)/norm))
            result.append(dict(trajectory=tr,duration_us=T_us,tau_us=5,distance_um=510,
                minimum_pupil_mode_overlap=min(overlaps),max_rms_phase_rad=max(maxr),
                warning='mode overlap is an optical approximation diagnostic, NOT atom survival'))
    pd.DataFrame(result).to_csv(ROOT/'acoustic_propagation_check.csv',index=False)
    # External data metadata, measurement waveform range and independent scatter exposure.
    a=pd.read_csv(ROOT/'zhang_source/dataset/extended_data_fig3a.csv')
    sc=pd.read_csv(ROOT/'zhang_source/dataset/extended_data_fig3c.csv')
    (ROOT/'paper_data_metadata.json').write_text(json.dumps(dict(zenodo='https://doi.org/10.5281/zenodo.19491381',
        md5='4701e05c2bb8abac833a586d8f09e428',license='CC-BY-4.0',
        measured_distance_um=float(a.position_um.max()),full_duration_ms=float(a.time_ms.max()),
        scattering_final_row=sc.iloc[-1].to_dict(),
        missing=['ED2 initial temperature','ED2 beam waist','ED2 depth specified independently of ED3','illuminated AOD radius and sound velocity','x/y optical geometry'],
        rb_source_search='publisher returned interstitial; no numerical source attached locally; digitized vector markers instead'),indent=2))

if __name__=='__main__':main()
