"""Escape-decision explainer figure: two real 140 uK atoms, captured vs escaped."""
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from level0_static_trap.constants import BOLTZMANN_CONSTANT
from level0_static_trap.potentials import gaussian_potential
from level1_transfer_1d.config import load_config
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level1_transfer_1d.dynamics import total_force, total_potential
from level1_transfer_1d.thermal import sample_thermal_initial_state
from level1_transfer_1d.waveforms import aod_center, aod_depth

HERE = Path(__file__).resolve().parent
CONFIG = HERE.parents[0] / "configs" / "transfer_1d.yaml"
OUT = HERE / "figures" / "escape_decision_example.png"
J_TO_UK = 1e6 / BOLTZMANN_CONSTANT

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "font.size": 13,
    "axes.titlesize": 15, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.labelweight": "bold",
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
    "axes.linewidth": 1.4, "lines.linewidth": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
})
C_AOD = "#1f6fb4"; C_SLM = "#d62728"; C_ATOM = "#2ca02c"; C_GREY = "#7f7f7f"


def propagate(config, x0, v0):
    times = np.arange(round(config.transfer.duration_s/config.trajectory.dt_s)+1)*config.trajectory.dt_s
    t, x, v, _ = velocity_verlet_time_dependent(
        lambda xx, tt: total_force(xx, tt, config),
        config.atom.mass_kg, x0, v0, times)
    return t, x, v


def find_pair(config):
    """Reproduce the exact two cases from the walkthrough."""
    rng = np.random.default_rng(20260811)
    m = config.atom.mass_kg
    cap = esc = None
    for _ in range(2000):
        s = sample_thermal_initial_state(config, rng)
        _, x, v = propagate(config, s.x0_m, s.v0_m_s)
        xf, vf = float(x[-1]), float(v[-1])
        U = float(gaussian_potential(xf, config.slm.depth_j, config.slm.waist_m, config.slm.center_m))
        E = 0.5*m*vf**2 + U
        if E < 0 and cap is None:
            cap = (s, t_, x_, v_) if False else (s, x, v)
        elif E >= 0 and esc is None:
            esc = (s, x, v)
        if cap and esc:
            break
    return cap, esc


def main():
    cfg = replace(load_config(CONFIG),
                  thermal=replace(load_config(CONFIG).thermal, temperature_uK=140.0))
    (s_cap, x_cap, v_cap), (s_esc, x_esc, v_esc) = find_pair(cfg)
    times = np.arange(round(cfg.transfer.duration_s/cfg.trajectory.dt_s)+1)*cfg.trajectory.dt_s
    t_us = times * 1e6

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.22,
                          height_ratios=[1.1, 1.0])

    # ===== Row 1: trajectories =====
    for col, (state, x, v, title, color) in enumerate([
        (s_cap, x_cap, v_cap, "CASE A - CAPTURED", C_ATOM),
        (s_esc, x_esc, v_esc, "CASE B - ESCAPED", C_SLM),
    ]):
        ax = fig.add_subplot(gs[0, col])
        ax.plot(t_us, x*1e6, color=C_ATOM, lw=2.5, label="Atom position")
        ax.plot(t_us, aod_center(times, cfg)*1e6, color=C_AOD, lw=2,
                label="AOD center (moving)")
        ax.axhline(0, color=C_SLM, ls="--", lw=2, label="SLM center (x=0)")
        ax.axvline(312, color=C_GREY, ls=":", lw=1.8)
        ax.text(312, ax.get_ylim()[1]*0.9, " move ends ", color=C_GREY,
                fontsize=10, ha="center")
        # mark initial state
        ax.scatter([0], [state.x0_m*1e6], s=90, color=C_ATOM,
                   edgecolor="black", zorder=5)
        ax.annotate(f"x0={state.x0_m*1e6:+.3f} um\nv0={state.v0_m_s:+.4f} m/s",
                    xy=(0, state.x0_m*1e6), xytext=(60, state.x0_m*1e6),
                    fontsize=10, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="black"))
        ax.set(xlabel="Time (us)", ylabel="Position x (um)", title=title,
               xlim=(0, 600))
        ax.legend(loc="best", fontsize=10)

    # ===== Row 2 left: final SLM potential + energy levels =====
    ax3 = fig.add_subplot(gs[1, 0])
    xx = np.linspace(-3, 3, 1201) * 1e-6
    U_slm = gaussian_potential(xx, cfg.slm.depth_j, cfg.slm.waist_m, 0) * J_TO_UK

    ax3.plot(xx*1e6, U_slm, color=C_SLM, lw=3, label="SLM potential U_SLM(x)")
    ax3.fill_between(xx*1e6, U_slm, 0, color=C_SLM, alpha=0.12)
    ax3.axhline(0, color="black", lw=2, label="Escape threshold (E = 0)")
    ax3.axhline(-140, color=C_GREY, ls=":", lw=1.8, label="Well bottom (-140 uK)")

    # Case A: captured, energy level below 0
    xf_c = float(x_cap[-1]); vf_c = float(v_cap[-1])
    U_c = float(gaussian_potential(xf_c, cfg.slm.depth_j, cfg.slm.waist_m, 0))*J_TO_UK
    KE_c = 0.5*cfg.atom.mass_kg*vf_c**2 * J_TO_UK
    E_c = KE_c + U_c
    ax3.axhline(E_c, color=C_ATOM, lw=2.5, ls="--")
    ax3.scatter([xf_c*1e6], [E_c], s=130, color=C_ATOM, edgecolor="black", zorder=5)
    ax3.text(2.0, E_c+3, f"E_A = {E_c:+.1f} uK\n(< 0, bound)",
             color=C_ATOM, fontsize=11, fontweight="bold")
    ax3.set(xlabel="x (um)", ylabel="Energy (uK)",
            title=f"Final state in SLM well  (CASE A captured)",
            xlim=(-3, 3), ylim=(-160, 20))
    ax3.legend(loc="upper right", fontsize=9)

    # ===== Row 2 right: case B on same well =====
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(xx*1e6, U_slm, color=C_SLM, lw=3, label="SLM potential U_SLM(x)")
    ax4.fill_between(xx*1e6, U_slm, 0, color=C_SLM, alpha=0.12)
    ax4.axhline(0, color="black", lw=2, label="Escape threshold (E = 0)")
    ax4.axhline(-140, color=C_GREY, ls=":", lw=1.8, label="Well bottom (-140 uK)")

    xf_e = float(x_esc[-1]); vf_e = float(v_esc[-1])
    # arrow: atom flew far away, so show it off-frame with annotation
    U_e = float(gaussian_potential(xf_e, cfg.slm.depth_j, cfg.slm.waist_m, 0))*J_TO_UK
    KE_e = 0.5*cfg.atom.mass_kg*vf_e**2 * J_TO_UK
    E_e = KE_e + U_e
    ax4.axhline(E_e, color=C_SLM, lw=2.5, ls="--")
    # the atom is at x=-5um (off-frame); draw an arrow pointing left
    ax4.annotate(f"atom flew to x={xf_e*1e6:+.2f} um\n(way outside the well)",
                 xy=(-3, E_e), xytext=(1.0, E_e+5),
                 fontsize=10, fontweight="bold", color=C_SLM,
                 arrowprops=dict(arrowstyle="->", color=C_SLM, lw=2))
    ax4.text(-2.9, E_e+3, f"E_B = {E_e:+.2f} uK\n(> 0, unbound)",
             color=C_SLM, fontsize=11, fontweight="bold")
    ax4.set(xlabel="x (um)", ylabel="Energy (uK)",
            title=f"Final state in SLM well  (CASE B escaped)",
            xlim=(-3, 3), ylim=(-160, 20))
    ax4.legend(loc="upper right", fontsize=9)

    fig.suptitle("How a 140 uK atom's capture is decided:  E_final = KE + U_SLM(x_f),  "
                 "capture  <=>  E < 0",
                 fontsize=15, fontweight="bold", y=0.997)
    fig.savefig(OUT, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
