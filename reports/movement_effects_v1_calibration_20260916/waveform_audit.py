"""Audit control waveforms against PDF vectors; never use survival values."""
from run_study import *
import pdfplumber
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    cfg,*_=case_inputs(SETTINGS,15,1,26091611,'matched')
    raw=waveform_for('ml_400us',cfg,ROOT/SETTINGS['ml_waveform'])
    with pdfplumber.open(WORK/'2403.12021v4.pdf') as pdf:
        page=pdf.pages[6]
        # Tick calibration: 0..400 us, 0..0.25 mK at six left ticks;
        # 0.5..2.0 um at four right ticks, all read from page 7 vector objects.
        x0,x1=267.295,370.111
        depth_y=np.array([227.766,215.359,202.951,190.543,178.131,165.723])
        depth_coef=np.polyfit(depth_y,np.arange(6)*.05,1)
        pos_coef=np.polyfit([212.775,197.78,182.784,167.793],[.5,1.,1.5,2.],1)
        t=np.linspace(0,400,15)
        extracted={}
        for name,ids,coef,ends in [('pos_um',range(16,29),pos_coef,[0,2.4]),
                                  ('depth_mK',range(2,15),depth_coef,[0,.280])]:
            centers=np.array([[(page.curves[i]['x0']+page.curves[i]['x1'])/2,
                (page.curves[i]['top']+page.curves[i]['bottom'])/2] for i in ids])
            assert np.max(abs((centers[:,0]-x0)/(x1-x0)*400-t[1:-1]))<.02
            extracted[name]=np.r_[ends[0],np.polyval(coef,centers[:,1]),ends[1]]
        extracted['grid_us']=t
        extracted['provenance']='13 internal PDF black dots plus 0/400 us endpoints; 14 intervals; endpoint values from Methods. Not AWG commands.'
        write(HERE/'ml_vector_nodes.json',extracted)
        summary={};fig,axs=plt.subplots(2,2,figsize=(10,7))
        for j,(name,idx,coef,scale,fn) in enumerate([
            ('pos_um',15,pos_coef,1e6,raw.center),('depth_mK',1,depth_coef,1/1.380649e-26,raw.depth)]):
            path=np.array(page.curves[idx]['pts'])
            tt=(path[:,0]-x0)/(x1-x0)*400;yy=np.polyval(coef,path[:,1])
            spline=CubicSpline(t,extracted[name])
            discrepancy=fn(tt*1e-6)*scale-yy
            summary[name]=dict(raw_path_rmse=float(np.sqrt(np.mean(discrepancy**2))),
                raw_path_max=float(np.max(abs(discrepancy))),vector_nodes_path_rmse=float(np.sqrt(np.mean((spline(tt)-yy)**2))))
            ax=axs[j,0];ax.plot(tt,yy,'k.',label='Published vector path');dense=np.linspace(0,400,1201)
            ax.plot(dense,fn(dense*1e-6)*scale,label='Frozen raster replay')
            ax.plot(dense,spline(dense),'--',label='PDF marker cubic (audit)')
            ax.set(xlabel='Time (us)',ylabel=name);ax.legend(fontsize=7)
        for j,duration in enumerate([200,400,600]):
            w=waveform_for(f'manual_{duration}us',cfg,ROOT/SETTINGS['ml_waveform'])
            summary[f'manual_{duration}us']=dict(vmax=float(1.5*2.4e-6/(duration*.52e-6)),
                jerk=float(-12*2.4e-6/(duration*.52e-6)**3),lens_max_shift_over_2zr=float(1.5*2.4e-6/(duration*.52e-6)/SETTINGS['vs_m_s']))
        dense=np.linspace(0,400e-6,1201)
        for profile in ['legacy_piecewise','paper_cubic']:
            w=waveform_for('manual_400us',cfg,ROOT/SETTINGS['ml_waveform'],manual_profile=profile)
            axs[0,1].plot(dense*1e6,w.center(dense)*1e6,label=profile)
            axs[1,1].plot(dense*1e6,w.center_velocity(dense),label=profile)
        axs[0,1].set(xlabel='Time (us)',ylabel='Manual position (um)');axs[1,1].set(xlabel='Time (us)',ylabel='Manual speed (m/s)')
        for ax in axs[:,1]:ax.legend(fontsize=8)
        summary['raw_ml_endpoints']=dict(position_um=[float(raw.center(t))*1e6 for t in [0,400e-6]],
            velocity_m_s=[float(raw.center_velocity(t)) for t in [0,400e-6]],
            depth_uK=[float(raw.depth(t))/1.380649e-29 for t in [0,400e-6]])
        summary['pdf_sha256']=hashlib.sha256((WORK/'2403.12021v4.pdf').read_bytes()).hexdigest()
        write(HERE/'waveform_audit.json',summary)
        fig.tight_layout();fig.savefig(HERE/'waveform_audit.png',dpi=160)

if __name__=='__main__':main()
