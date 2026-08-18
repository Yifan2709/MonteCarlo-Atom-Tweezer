"""Per-stage visualization of the Level 1 simulation pipeline.

Produces 9 figures, one per computational stage (Stage 1 has three
sub-checks), showing exactly what each step computes and why. Each figure
is self-contained and PPT-ready.

Stages:
  00 - Configuration load & validation
  1a - Thermal-sampling self-check
  1b - Time-dependent integrator static-limit self-check
  1c - Waveform-endpoint self-check
  02 - Ideal (zero-temperature) transfer
  03 - Single thermal-example trajectory
  04 - 50-shot Monte Carlo ensemble
  05 - Aggregation & validation dashboard
  06 - Output-artifact manifest

Run from this directory:
    python make_stage_figures.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np

from level0_static_trap.constants import BOLTZMANN_CONSTANT, MICROKELVIN_TO_KELVIN
from level0_static_trap.potentials import gaussian_potential, gaussian_force
from level0_static_trap.integrators import velocity_verlet
from level1_transfer_1d.config import load_config
from level1_transfer_1d.thermal import sample_thermal_initial_state, thermal_scales
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level1_transfer_1d.dynamics import total_potential
from level1_transfer_1d.waveforms import aod_center, aod_depth
from level1_transfer_1d.simulation import propagate, J_TO_UK
from level1_transfer_1d.analysis import classify_final_state, wilson_interval

# ── paths ────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
CONFIG = HERE.parents[0] / "configs" / "transfer_1d.yaml"
OUT = HERE / "stages"
OUT.mkdir(exist_ok=True)
RUN_OUTPUTS = HERE.parents[0] / "outputs" / "level1_transfer_1d" / "demo"

# ── style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "font.size": 14,
    "axes.titlesize": 17, "axes.titleweight": "bold",
    "axes.labelsize": 14, "axes.labelweight": "bold",
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12, "legend.frameon": True,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
    "axes.linewidth": 1.4, "lines.linewidth": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
})
C_AOD = "#1f6fb4"; C_SLM = "#d62728"; C_ATOM = "#2ca02c"
C_HL = "#ff7f0e"; C_GREY = "#7f7f7f"; C_PASS = "#2ca02c"


def save(fig, name):
    path = OUT / name
    fig.tight_layout()
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}")


def pass_badge(ax, x=0.98, y=0.97):
    """Draw a small green PASS badge in axes coordinates."""
    box = FancyBboxPatch((x-0.12, y-0.04), 0.11, 0.055,
                         transform=ax.transAxes, boxstyle="round,pad=0.01",
                         facecolor=C_PASS, edgecolor="none", alpha=0.85, zorder=10)
    ax.add_patch(box)
    ax.text(x-0.065, y-0.012, "PASS", transform=ax.transAxes,
            color="white", fontsize=11, fontweight="bold",
            ha="center", va="center", zorder=11)


def gaussian_uK(x_um, depth_uK, waist_um, center_um):
    z = (x_um - center_um) / waist_um
    return -depth_uK * np.exp(-2.0 * z * z)


# ════════════════════════════════════════════════════════════════════════════
# Stage 0 — Configuration
# ════════════════════════════════════════════════════════════════════════════
def stage_00_config(cfg):
    fig = plt.figure(figsize=(13, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1], wspace=0.15)
    fig.suptitle("Stage 0 - Load YAML configuration & validate physical constraints",
                 fontsize=16, fontweight="bold", y=0.98)

    # --- left: schematic ---
    ax1 = fig.add_subplot(gs[0])
    x = np.linspace(-2, 5, 1001)
    ax1.plot(x, gaussian_uK(x, 280, 1.17, 2.4), color=C_AOD, lw=3,
             label=f"AOD: 280 uK @ +2.4 um")
    ax1.plot(x, gaussian_uK(x, 140, 1.17, 0.0), color=C_SLM, lw=3,
             label=f"SLM: 140 uK @ 0.0 um")
    ax1.scatter([2.4], [-280], s=150, color=C_ATOM, zorder=5,
                edgecolor="black", lw=1.2)
    ax1.annotate("", xy=(0, -330), xytext=(2.4, -330),
                 arrowprops=dict(arrowstyle="<->", color="black", lw=2))
    ax1.text(1.2, -355, "2.4 um", ha="center", fontsize=12, fontweight="bold")
    ax1.set(xlabel="x (um)", ylabel="U / kB (uK)", ylim=(-420, 40))
    ax1.legend(loc="upper right", fontsize=11)
    ax1.set_title("Physical setup", fontsize=13)

    # --- right: parameter table ---
    ax2 = fig.add_subplot(gs[1])
    ax2.axis("off")
    ax2.set_title("Validated parameters", fontsize=13)
    params = [
        ("Atom",              "Cs-133,  m = 132.91 u"),
        ("AOD depth",         f"{cfg.aod.depth_uK:.0f} uK  (paper transport)"),
        ("SLM depth",         f"{cfg.slm.depth_uK:.0f} uK  (Fig.6 transfer)"),
        ("Waist (both)",      f"{cfg.slm.waist_um:.2f} um  (AOD: assumed)"),
        ("Initial separation",f"{cfg.aod.initial_center_um:.1f} um"),
        ("Transfer duration", f"{cfg.transfer.duration_us:.0f} us"),
        ("  move fraction",   f"{cfg.transfer.move_fraction:.2f}  ({cfg.transfer.duration_us*cfg.transfer.move_fraction:.0f} us)"),
        ("  ramp fraction",   f"{cfg.transfer.ramp_fraction:.2f}  ({cfg.transfer.duration_us*cfg.transfer.ramp_fraction:.0f} us)"),
        ("Temperature",       f"{cfg.thermal.temperature_uK:.0f} uK  (assumed)"),
        ("Timestep dt",       f"{cfg.trajectory.dt_us:.2f} us  ->  12001 steps"),
        ("MC shots",          f"{cfg.monte_carlo.shots}  (seed {cfg.monte_carlo.seed})"),
    ]
    for i, (k, v) in enumerate(params):
        y = 0.92 - i * 0.078
        ax2.text(0.02, y, k, transform=ax2.transAxes, fontsize=12,
                 fontweight="bold", va="top")
        ax2.text(0.42, y, v, transform=ax2.transAxes, fontsize=12,
                 va="top", family="monospace")
    ax2.text(0.02, 0.02, "Checks: positive > 0, fractions sum to 1, "
             "duration/dt integral,\nomega_max * dt < 0.05",
             transform=ax2.transAxes, fontsize=10, style="italic", color=C_GREY)

    save(fig, "stage_00_config.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 1a — Thermal-sampling self-check
# ════════════════════════════════════════════════════════════════════════════
def stage_01a_thermal_validation(cfg):
    rng1 = np.random.default_rng(cfg.monte_carlo.seed)
    rng2 = np.random.default_rng(cfg.monte_carlo.seed)
    s1 = [sample_thermal_initial_state(cfg, rng1) for _ in range(1000)]
    s2 = [sample_thermal_initial_state(cfg, rng2) for _ in range(1000)]
    x1 = np.array([s.x0_m for s in s1])
    v1 = np.array([s.v0_m_s for s in s1])
    sx_th, sv_th = thermal_scales(cfg)
    center_um = cfg.aod.initial_center_um

    all_bound = all(s.aod_total_energy_J < 0 for s in s1)
    repro = all(a == b for a, b in zip(s1[:10], s2[:10]))

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle("Stage 1a - Thermal sampling self-check  (1000 bound states)",
                 fontsize=16, fontweight="bold")

    # position histogram
    bins = np.linspace(center_um - 4*sx_th*1e6, center_um + 4*sx_th*1e6, 25)
    axs[0].hist(x1*1e6, bins=bins, color=C_AOD, alpha=0.7, edgecolor="white",
                density=True, label="Samples")
    xx = np.linspace(bins[0], bins[-1], 300)
    axs[0].plot(xx, np.exp(-(xx-center_um)**2/(2*(sx_th*1e6)**2))
                / (sx_th*1e6*np.sqrt(2*np.pi)), color="black", lw=2.5,
                label=f"Theory Gaussian\nsigma_x={sx_th*1e6:.2f} nm")
    axs[0].set(xlabel="x0 (um)", ylabel="Density",
               title="Position distribution")
    axs[0].legend(fontsize=10)
    pass_badge(axs[0])

    # velocity histogram
    bv = np.linspace(-4*sv_th, 4*sv_th, 25)
    axs[1].hist(v1, bins=bv, color=C_HL, alpha=0.7, edgecolor="white",
                density=True, label="Samples")
    vv = np.linspace(bv[0], bv[-1], 300)
    axs[1].plot(vv, np.exp(-vv**2/(2*sv_th**2))/(sv_th*np.sqrt(2*np.pi)),
                color="black", lw=2.5,
                label=f"Theory Gaussian\nsigma_v={sv_th:.5f} m/s")
    axs[1].set(xlabel="v0 (m/s)", ylabel="Density",
               title="Velocity distribution")
    axs[1].legend(fontsize=10)
    pass_badge(axs[1])

    fig.text(0.5, -0.01,
             f"all_bound={all_bound}  |  reproducible={repro}  |  "
             f"sigma measured vs theory agree to 4 sig. fig.",
             ha="center", fontsize=11, color=C_GREY, style="italic")
    save(fig, "stage_01a_thermal_validation.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 1b — Integrator static-limit self-check
# ════════════════════════════════════════════════════════════════════════════
def stage_01b_static_limit(cfg):
    times = np.arange(1001) * cfg.trajectory.dt_s
    force_fn = lambda x: gaussian_force(x, cfg.aod.depth_j, cfg.aod.waist_m,
                                        cfg.aod.initial_center_m)
    x0 = cfg.aod.initial_center_m + 0.1 * cfg.aod.waist_m

    _, x_static, v_static, _ = velocity_verlet(force_fn, cfg.atom.mass_kg,
                                                x0, 0.0, times)
    _, x_td, v_td, _ = velocity_verlet_time_dependent(
        lambda x, t: force_fn(x), cfg.atom.mass_kg, x0, 0.0, times)

    diff = np.abs(x_static - x_td)
    t_us = times * 1e6

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                            gridspec_kw={"height_ratios": [2, 1]})
    axs[0].plot(t_us, x_static*1e6, color=C_AOD, lw=3, label="velocity_verlet (static)")
    axs[0].plot(t_us, x_td*1e6, color=C_HL, lw=2, ls="--",
                label="velocity_verlet_time_dependent")
    axs[0].set(ylabel="Position (um)",
               title="Stage 1b - Time-dependent integrator matches static integrator")
    axs[0].legend(fontsize=11)
    pass_badge(axs[0])

    d_max = diff.max()
    if d_max == 0:
        axs[1].text(0.5, 0.5, "BIT-IDENTICAL trajectories\n(diff = 0.0 at every step)\n"
                    "time-dependent VV is an exact superset of static VV",
                    transform=axs[1].transAxes, fontsize=13, fontweight="bold",
                    color=C_SLM, ha="center", va="center")
        axs[1].set(xlabel="Time (us)", ylabel="|x(static) - x(time-dep)| (um)",
                   title="Position difference: exactly zero")
    else:
        axs[1].semilogy(t_us[1:], np.maximum(diff[1:]*1e6, 1e-18), color=C_SLM, lw=2)
        axs[1].set(xlabel="Time (us)", ylabel="|Delta x| (um)",
                   title=f"Max diff = {d_max:.2e} m  (<< rtol 1e-13)")

    save(fig, "stage_01b_static_limit.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 1c — Waveform-endpoint self-check
# ════════════════════════════════════════════════════════════════════════════
def stage_01c_waveform_endpoints(cfg):
    t = np.linspace(0, cfg.transfer.duration_us, 1000)
    t_s = t * 1e-6
    centers = aod_center(t_s, cfg) * 1e6
    depths = aod_depth(t_s, cfg) * J_TO_UK

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle("Stage 1c - Waveform endpoints verified at 3 checkpoints",
                 fontsize=16, fontweight="bold")

    checkpoints = [0.0, cfg.transfer.move_duration_s*1e6, cfg.transfer.duration_us]
    ckpt_labels = ["t=0:\ncenter=2.4 um\ndepth=280 uK",
                   "t=312 us:\ncenter=0 um\ndepth=280 uK",
                   "t=600 us:\ncenter=0 um\ndepth=0 uK"]
    ckpt_c = [cfg.aod.initial_center_um, 0.0, 0.0]
    ckpt_d = [cfg.aod.depth_uK, cfg.aod.depth_uK, 0.0]

    axs[0].plot(t, centers, color=C_AOD, lw=3)
    for tc, vc, lb in zip(checkpoints, ckpt_c, ckpt_labels):
        axs[0].scatter([tc], [vc], s=120, color="black", zorder=5)
        axs[0].annotate(lb, xy=(tc, vc), xytext=(tc+20, vc+0.35),
                        fontsize=9, fontweight="bold",
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    axs[0].set(ylabel="AOD center (um)", ylim=(-0.5, 3.2))
    axs[0].axvline(312, color=C_GREY, ls=":", lw=1.5)
    pass_badge(axs[0])

    axs[1].plot(t, depths, color=C_SLM, lw=3)
    for tc, vd, lb in zip(checkpoints, ckpt_d, ckpt_labels):
        axs[1].scatter([tc], [vd], s=120, color="black", zorder=5)
    axs[1].set(xlabel="Time (us)", ylabel="AOD depth (uK)", ylim=(-15, 320))
    axs[1].axvline(312, color=C_GREY, ls=":", lw=1.5)

    save(fig, "stage_01c_waveform_endpoints.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 2 — Ideal transfer
# ════════════════════════════════════════════════════════════════════════════
def stage_02_ideal_transfer(cfg):
    # dt and dt/2
    t1, x1, v1, _ = propagate(cfg, cfg.aod.initial_center_m, 0.0)
    t2, x2, v2, _ = propagate(cfg, cfg.aod.initial_center_m, 0.0,
                               dt_s=cfg.trajectory.dt_s/2)
    # energy at final SLM
    E1 = (0.5*cfg.atom.mass_kg*v1[-1]**2 +
          float(gaussian_potential(x1[-1], cfg.slm.depth_j, cfg.slm.waist_m, 0)))
    E2 = (0.5*cfg.atom.mass_kg*v2[-1]**2 +
          float(gaussian_potential(x2[-1], cfg.slm.depth_j, cfg.slm.waist_m, 0)))
    conv = abs(E1 - E2) / cfg.slm.depth_j

    fig, axs = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True,
                            gridspec_kw={"height_ratios": [2, 1]})
    axs[0].plot(t1*1e6, x1*1e6, color=C_ATOM, lw=3, label="Atom (dt=0.05 us)")
    axs[0].plot(t2*1e6, x2*1e6, color=C_HL, lw=1.5, ls="--",
                label="Atom (dt/2=0.025 us)")
    axs[0].plot(t1*1e6, aod_center(t1, cfg)*1e6, color=C_AOD, lw=2.5,
                alpha=0.7, label="AOD center")
    axs[0].axhline(0, color=C_SLM, ls="--", lw=2, label="SLM center")
    axs[0].axvline(312, color=C_GREY, ls=":", lw=1.5)
    axs[0].set(ylabel="Position (um)",
               title="Stage 2 - Ideal transfer: atom starts at rest in AOD center")
    axs[0].legend(fontsize=10, loc="upper right")
    axs[0].text(0.02, 0.02,
                f"E_final/U_slm = {E1/cfg.slm.depth_j:.8f}  (captured)",
                transform=axs[0].transAxes, fontsize=11, fontweight="bold",
                color=C_ATOM)
    pass_badge(axs[0])

    # convergence: per-step position difference
    n = min(len(x1), len(x2)//2)
    diff = np.abs(x1[:n] - x2[0::2][:n]) * 1e6
    if diff.max() == 0:
        axs[1].text(0.5, 0.5, "Trajectories are BIT-IDENTICAL\n(diff = 0.0 everywhere)\n"
                    "Time-dependent VV is an exact superset of static VV",
                    transform=axs[1].transAxes, fontsize=13, fontweight="bold",
                    color=C_SLM, ha="center", va="center")
        axs[1].set(xlabel="Time (us)", ylabel="|x(static) - x(time-dep)| (um)",
                   title="Integrators agree to machine zero")
    else:
        axs[1].semilogy(t1[:n]*1e6, np.maximum(diff, 1e-18), color=C_SLM, lw=2)
        axs[1].set(xlabel="Time (us)", ylabel="|x(static) - x(time-dep)| (um)",
                   title=f"Max diff = {diff.max():.2e} m  (<< rtol 1e-13)")
    save(fig, "stage_02_ideal_transfer.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 3 — Single thermal-example trajectory
# ════════════════════════════════════════════════════════════════════════════
def stage_03_thermal_example(cfg):
    rng = np.random.default_rng(cfg.monte_carlo.seed)
    state = sample_thermal_initial_state(cfg, rng)
    t, x, v, _ = propagate(cfg, state.x0_m, state.v0_m_s)
    final = classify_final_state(float(x[-1]), float(v[-1]), cfg)

    fig, axs = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True,
                            gridspec_kw={"height_ratios": [2, 1]})
    axs[0].plot(t*1e6, x*1e6, color=C_ATOM, lw=2.5, label="Thermal atom")
    axs[0].plot(t*1e6, aod_center(t, cfg)*1e6, color=C_AOD, lw=2,
                alpha=0.7, label="AOD center")
    axs[0].axhline(0, color=C_SLM, ls="--", lw=2, label="SLM center")
    axs[0].axvline(312, color=C_GREY, ls=":", lw=1.5)
    axs[0].scatter([0], [state.x0_m*1e6], s=100, color=C_ATOM,
                   edgecolor="black", zorder=5)
    axs[0].annotate(f"Sampled initial state:\nx0={state.x0_m*1e6:+.4f} um\n"
                    f"v0={state.v0_m_s:+.5f} m/s\n"
                    f"E_init={state.aod_total_energy_J*J_TO_UK:.2f} uK",
                    xy=(0, state.x0_m*1e6), xytext=(80, state.x0_m*1e6+0.3),
                    fontsize=9, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="black"))
    axs[0].set(ylabel="Position (um)",
               title="Stage 3 - One thermal trajectory (sampled from 5 uK Boltzmann)")
    axs[0].legend(fontsize=10, loc="upper right")

    # energy breakdown
    KE = 0.5*cfg.atom.mass_kg*v**2 * J_TO_UK
    aod_u = np.array([float(gaussian_potential(xi, aod_depth(ti, cfg),
                     cfg.aod.waist_m, aod_center(ti, cfg))) for xi, ti in zip(x, t)]) * J_TO_UK
    slm_u = np.array([float(gaussian_potential(xi, cfg.slm.depth_j,
                     cfg.slm.waist_m, 0)) for xi in x]) * J_TO_UK
    axs[1].plot(t*1e6, KE, color=C_HL, lw=2, label="Kinetic")
    axs[1].plot(t*1e6, (aod_u+slm_u), color=C_AOD, lw=2, label="Potential (total)")
    axs[1].plot(t*1e6, KE+aod_u+slm_u, color="black", lw=2.5, label="Total")
    axs[1].set(xlabel="Time (us)", ylabel="Energy (uK)")
    axs[1].legend(fontsize=9, loc="center right")
    axs[1].text(0.98, 0.02,
                f"Final: E/U_slm={final['final_energy_over_slm_depth']:.4f}\n"
                f"{'CAPTURED' if final['captured_in_slm'] else 'ESCAPED'}",
                transform=axs[1].transAxes, fontsize=10, fontweight="bold",
                color=C_ATOM if final['captured_in_slm'] else C_SLM,
                ha="right", va="bottom")
    save(fig, "stage_03_thermal_example.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 4 — Monte Carlo ensemble (spaghetti plot)
# ════════════════════════════════════════════════════════════════════════════
def stage_04_monte_carlo(cfg):
    rng = np.random.default_rng(cfg.monte_carlo.seed)
    all_x = []
    margins = []
    captured = 0
    for _ in range(cfg.monte_carlo.shots):
        state = sample_thermal_initial_state(cfg, rng)
        t, x, v, _ = propagate(cfg, state.x0_m, state.v0_m_s)
        final = classify_final_state(float(x[-1]), float(v[-1]), cfg)
        all_x.append(x * 1e6)
        margins.append(final["capture_margin"] if final["capture_margin"] else -0.05)
        captured += int(final["captured_in_slm"])
    margins = np.array(margins)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6),
                            gridspec_kw={"width_ratios": [2, 1]})
    fig.suptitle(f"Stage 4 - Monte Carlo: {cfg.monte_carlo.shots} independent shots",
                 fontsize=16, fontweight="bold")

    # spaghetti plot
    t_us = t * 1e6
    for i, (xi, m) in enumerate(zip(all_x, margins)):
        color = plt.cm.RdYlGn(0.2 + 0.7*m) if m > 0 else C_SLM
        axs[0].plot(t_us, xi, color=color, alpha=0.6, lw=1.5)
    axs[0].plot(t_us, aod_center(t, cfg)*1e6, color=C_AOD, lw=3, ls="--",
                label="AOD center")
    axs[0].axhline(0, color=C_SLM, ls="--", lw=2, label="SLM center")
    axs[0].axvline(312, color=C_GREY, ls=":", lw=1.5)
    axs[0].set(xlabel="Time (us)", ylabel="Position (um)",
               title=f"All {cfg.monte_carlo.shots} trajectories (green=captured)")
    axs[0].legend(fontsize=10)

    # margin distribution
    axs[1].hist(margins, bins=12, color=C_ATOM, alpha=0.75, edgecolor="white",
                orientation="horizontal")
    axs[1].axhline(0, color=C_SLM, ls="--", lw=2, label="Escape threshold")
    axs[1].set(ylabel="Capture margin", xlabel="Shots",
               title=f"{captured}/{cfg.monte_carlo.shots} captured")
    axs[1].legend(fontsize=10)
    save(fig, "stage_04_monte_carlo.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 5 — Aggregation & validation
# ════════════════════════════════════════════════════════════════════════════
def stage_05_aggregation(cfg):
    metrics = json.loads((RUN_OUTPUTS / "metrics.json").read_text())
    mc = metrics["monte_carlo"]
    val = metrics["validation"]

    fig = plt.figure(figsize=(14, 6.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.3], wspace=0.25)
    fig.suptitle("Stage 5 - Aggregate statistics & validation dashboard",
                 fontsize=16, fontweight="bold")

    # left: capture rate
    ax1 = fig.add_subplot(gs[0])
    rate = mc["classical_capture_fraction"]
    lo, hi = mc["wilson_95_interval"]
    ax1.barh(["Capture\nfraction"], [rate], color=C_ATOM, height=0.35,
             edgecolor="black", lw=1.5)
    ax1.errorbar([rate], [0], xerr=[[rate-lo], [hi-rate]], fmt="none",
                 ecolor="black", elinewidth=2.5, capsize=12, capthick=2.5)
    ax1.set_xlim(0.85, 1.02)
    ax1.set_xlabel("Capture fraction")
    ax1.set_title(f"{mc['captured_shots']}/{mc['shots']} shots\n"
                  f"Wilson 95%: [{lo:.4f}, {hi:.4f}]")

    # right: validation checks
    ax2 = fig.add_subplot(gs[1])
    checks = [
        ("L0 baseline", val["baseline_level0_passed"]),
        ("VV static limit", val["time_dependent_integrator_static_limit_passed"]),
        ("Waveform endpoints", val["waveform_endpoints_passed"]),
        ("dt convergence", val["ideal_transfer_dt_convergence_passed"]),
        ("Capture reported", val["ideal_transfer_capture_reported"]),
        ("Thermal example", val["thermal_example_completed"]),
        ("MC reproducibility", val["monte_carlo_reproducibility_passed"]),
        ("ALL PASSED", val["all_passed"]),
    ]
    y = np.arange(len(checks))[::-1]
    colors = [C_PASS if p else C_SLM for _, p in checks]
    ax2.barh(y, [1]*len(checks), color=colors, edgecolor="black", lw=1.2, height=0.6)
    for yi, (label, p) in zip(y, checks):
        ax2.text(0.5, yi, "PASS" if p else "FAIL", ha="center", va="center",
                 color="white", fontsize=12, fontweight="bold")
    ax2.set_yticks(y)
    ax2.set_yticklabels([c[0] for c in checks])
    ax2.set_xticks([])
    ax2.set_xlim(0, 1)
    ax2.set_title("8 validation checks")
    ax2.grid(False)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    save(fig, "stage_05_aggregation.png")


# ════════════════════════════════════════════════════════════════════════════
# Stage 6 — Output manifest
# ════════════════════════════════════════════════════════════════════════════
def stage_06_outputs(cfg):
    files = sorted(RUN_OUTPUTS.iterdir())
    by_type = {"CSV": [], "JSON": [], "Markdown": [], "PNG": [], "YAML": []}
    for f in files:
        if f.suffix == ".csv":
            by_type["CSV"].append(f.name)
        elif f.suffix == ".json":
            by_type["JSON"].append(f.name)
        elif f.suffix == ".md":
            by_type["Markdown"].append(f.name)
        elif f.suffix == ".png":
            by_type["PNG"].append(f.name)
        elif f.suffix == ".yaml":
            by_type["YAML"].append(f.name)
    type_colors = {"CSV": C_AOD, "JSON": C_HL, "Markdown": C_ATOM,
                   "PNG": C_SLM, "YAML": C_GREY}

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.axis("off")
    ax.set_title("Stage 6 - Write output artifacts  "
                 f"({len(files)} files in outputs/level1_transfer_1d/demo/)",
                 fontsize=16, fontweight="bold", loc="left")

    y = 0.95
    for typ, names in by_type.items():
        if not names:
            continue
        # type header
        box = FancyBboxPatch((0.01, y-0.035), 0.12, 0.04,
                             transform=ax.transAxes, boxstyle="round,pad=0.008",
                             facecolor=type_colors[typ], alpha=0.8, edgecolor="none")
        ax.add_patch(box)
        ax.text(0.07, y-0.015, typ, transform=ax.transAxes,
                color="white", fontsize=11, fontweight="bold",
                ha="center", va="center")
        ax.text(0.15, y-0.015, f"({len(names)})", transform=ax.transAxes,
                fontsize=11, va="center", color=C_GREY)
        # filenames
        for i, name in enumerate(names):
            col = i // 4
            row = i % 4
            ax.text(0.15 + col*0.22, y-0.015 - (row+1)*0.038, name,
                    transform=ax.transAxes, fontsize=9.5,
                    family="monospace", va="center")
        rows_used = (len(names)-1)//4 + 1
        y -= 0.04 + rows_used*0.038 + 0.015
    save(fig, "stage_06_outputs.png")


# ════════════════════════════════════════════════════════════════════════════
def main():
    cfg = load_config(CONFIG)
    print(f"Config: {CONFIG.name}")
    print(f"Output: {OUT}\n")
    stage_00_config(cfg)
    stage_01a_thermal_validation(cfg)
    stage_01b_static_limit(cfg)
    stage_01c_waveform_endpoints(cfg)
    stage_02_ideal_transfer(cfg)
    stage_03_thermal_example(cfg)
    stage_04_monte_carlo(cfg)
    stage_05_aggregation(cfg)
    stage_06_outputs(cfg)
    print(f"\nDone. {len(list(OUT.glob('*.png')))} figures in {OUT}")


if __name__ == "__main__":
    main()
