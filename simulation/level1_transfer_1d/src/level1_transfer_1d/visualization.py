"""Compact diagnostics for the Level 1 run."""
from pathlib import Path
import os, tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"level1_transfer_matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from level0_static_trap.constants import BOLTZMANN_CONSTANT
from .dynamics import total_potential
from .waveforms import aod_center, aod_depth

def _save(fig,path): fig.tight_layout(); fig.savefig(path,dpi=220); plt.close(fig)
def plot_waveform(t,c,path):
    fig,axs=plt.subplots(3,1,figsize=(7,7),sharex=True); u=t*1e6
    axs[0].plot(u,aod_center(t,c)*1e6); axs[1].plot(u,aod_depth(t,c)/BOLTZMANN_CONSTANT*1e6); axs[2].plot(u,np.full_like(t,c.slm.depth_uK))
    for ax,y in zip(axs,("AOD center (um)","AOD depth (uK)","SLM depth (uK)")): ax.set_ylabel(y); ax.grid(alpha=.25)
    axs[-1].set_xlabel("Time (us)"); _save(fig,path)
def plot_potential_snapshot(c,t,path,title):
    x=np.linspace(-3,5,2001)*1e-6; y=total_potential(x,t,c)/BOLTZMANN_CONSTANT*1e6
    fig,ax=plt.subplots(figsize=(7,4.2)); ax.plot(x*1e6,y,label="total potential"); ax.set(xlabel="x (um)",ylabel="U/kB (uK)",title=title); ax.grid(alpha=.25); ax.legend(); _save(fig,path)
def plot_position(record,c,path,title):
    t=record["time_s"]*1e6; fig,ax=plt.subplots(figsize=(7,4.2)); ax.plot(t,record["x_m"]*1e6,label="atom"); ax.plot(t,record["aod_center_m"]*1e6,label="AOD center"); ax.axhline(c.slm.center_um,color="k",ls="--",label="SLM center"); ax.set(xlabel="Time (us)",ylabel="Position (um)",title=title); ax.grid(alpha=.25); ax.legend(); _save(fig,path)
def plot_energy(record,path,title):
    t=record["time_s"]*1e6; fig,ax=plt.subplots(figsize=(7,4.2)); ax.plot(t,record["total_energy_J"]/BOLTZMANN_CONSTANT*1e6,label="mechanical energy"); ax.set(xlabel="Time (us)",ylabel="Energy/kB (uK)",title=title); ax.grid(alpha=.25); ax.legend(); _save(fig,path)
def plot_phase(mc,path,initial):
    x=[r["x0_um" if initial else "final_x_um"] for r in mc]; v=[r["v0_m_s" if initial else "final_v_m_s"] for r in mc]
    fig,ax=plt.subplots(figsize=(6,4.5)); ax.scatter(x,v,s=22); ax.set(xlabel="x (um)",ylabel="v (m/s)",title=("Initial" if initial else "Final")+" phase space"); ax.grid(alpha=.25); _save(fig,path)
def plot_distribution(values,path,xlabel,threshold=None):
    fig,ax=plt.subplots(figsize=(6.5,4.2)); ax.hist(values,bins=12,edgecolor="white");
    if threshold is not None: ax.axvline(threshold,color="red",ls="--",label="escape threshold"); ax.legend()
    ax.set(xlabel=xlabel,ylabel="Shots"); ax.grid(alpha=.2); _save(fig,path)
