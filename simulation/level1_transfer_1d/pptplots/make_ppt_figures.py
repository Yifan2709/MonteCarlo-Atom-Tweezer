"""Level 1 AOD->SLM transfer: PPT-ready visualization suite.

Reads the artifacts produced by `level1_transfer_1d.cli` from
outputs/level1_transfer_1d/demo/ and writes 12 self-contained, presentation-
ready PNGs into pptplots/figures/. Each figure tells one part of the story:

    01 setup          - the physical scenario (two traps, atom, separation)
    02 waveform       - the engineered AOD move + ramp schedule
    03 potential_evo  - three snapshots of the live total potential
    04 ideal_pos      - ideal atom riding the moving AOD into the SLM
    05 ideal_energy   - ideal mechanical energy (time-dependent, NOT conserved)
    06 thermal_dist   - how finite-T initial states are sampled
    07 thermal_pos    - one example thermal trajectory
    08 phase_space    - initial vs final MC phase space (cooling + focusing)
    09 mc_capture     - per-shot capture margin
    10 mc_excitation  - per-shot motional excitation change
    11 mc_summary     - capture rate + Wilson CI
    12 validation     - the 8-check validation dashboard

Style: white background, large fonts, thick lines, semantic colors
(AOD=blue, SLM=red, atom=green, escape threshold=dashed red).

Run from this directory:
    python make_ppt_figures.py
"""
from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

# ----------------------------------------------------------------------------
# Style: presentation-ready
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 15,
    "axes.titlesize": 18,
    "axes.titleweight": "bold",
    "axes.labelsize": 15,
    "axes.labelweight": "bold",
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 13,
    "legend.frameon": True,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "axes.linewidth": 1.5,
    "lines.linewidth": 3.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Semantic colors (consistent across every figure)
C_AOD = "#1f6fb4"      # blue  - AOD trap
C_SLM = "#d62728"      # red   - SLM trap
C_ATOM = "#2ca02c"     # green - atom
C_THRESHOLD = "#d62728"
C_NEUTRAL = "#7f7f7f"
C_HIGHLIGHT = "#ff7f0e"

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "figures"
OUT_DIR.mkdir(exist_ok=True)
DATA = HERE.parent / "outputs" / "level1_transfer_1d" / "demo"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def load_csv(name: str) -> list[dict]:
    with open(DATA / name, newline="") as f:
        return list(csv.DictReader(f))


def load_ideal() -> dict[str, np.ndarray]:
    rows = load_csv("ideal_transfer_trajectory.csv")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def load_thermal() -> dict[str, np.ndarray]:
    rows = load_csv("thermal_example_trajectory.csv")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def load_mc() -> list[dict]:
    rows = load_csv("mc_shots.csv")
    out = []
    for r in rows:
        rr = dict(r)
        for k in ("shot_id",):
            rr[k] = int(r[k])
        for k, v in r.items():
            if k == "captured_in_slm":
                rr[k] = v == "True"
            elif k != "shot_id":
                rr[k] = float(v)
        out.append(rr)
    return out


def save(fig, name: str):
    path = OUT_DIR / name
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.name}")


def gaussian_uK(x_um, depth_uK, waist_um, center_um):
    """1D Gaussian optical dipole potential in uK (sign kept: negative well)."""
    z = (x_um - center_um) / waist_um
    return -depth_uK * np.exp(-2.0 * z * z)


# ----------------------------------------------------------------------------
# 01 - Setup / scenario
# ----------------------------------------------------------------------------
def fig_01_setup():
    x = np.linspace(-3, 6, 1201)
    aod = gaussian_uK(x, 280.0, 1.17, 2.4)
    slm = gaussian_uK(x, 140.0, 1.17, 0.0)
    total = aod + slm

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(x, aod, color=C_AOD, label="AOD trap (280 uK @ x=2.4 um)")
    ax.plot(x, slm, color=C_SLM, label="SLM trap (140 uK @ x=0)")
    ax.plot(x, total, color=C_NEUTRAL, ls="--", lw=2, alpha=0.7,
            label="Sum (initial overlap is small)")

    # atom marker
    ax.scatter([2.4], [-280], s=220, color=C_ATOM, zorder=5,
               edgecolor="black", linewidth=1.5, label="Atom (starts in AOD)")
    # separation arrow
    ax.annotate("", xy=(0.0, -340), xytext=(2.4, -340),
                arrowprops=dict(arrowstyle="<->", color="black", lw=2))
    ax.text(1.2, -360, "2.4 um\nseparation", ha="center", va="top",
            fontsize=13, fontweight="bold")

    ax.set(xlabel="Radial position x (um)", ylabel="U / kB (uK)",
           title="Step 0 - Two independent traps before transfer",
           xlim=(-3, 6), ylim=(-430, 40))
    ax.axhline(0, color="black", lw=0.8)
    ax.legend(loc="upper right")
    save(fig, "01_setup_two_traps.png")


# ----------------------------------------------------------------------------
# 02 - Waveform schedule
# ----------------------------------------------------------------------------
def fig_02_waveform():
    wf = load_csv("transfer_waveform.csv")
    t = np.array([float(r["time_us"]) for r in wf])
    aod_c = np.array([float(r["aod_center_um"]) for r in wf])
    aod_d = np.array([float(r["aod_depth_uK"]) for r in wf])
    slm_d = np.array([float(r["slm_depth_uK"]) for r in wf])

    T_move = 312.0
    T_total = 600.0

    fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True)

    # phase shading
    for ax in axs:
        ax.axvspan(0, T_move, color=C_AOD, alpha=0.06)
        ax.axvspan(T_move, T_total, color=C_SLM, alpha=0.06)
    axs[0].text(T_move / 2, 2.55, "MOVE\n(312 us)", ha="center",
                fontweight="bold", color=C_AOD, fontsize=13)
    axs[0].text((T_move + T_total) / 2, 2.55, "RAMP\n(288 us)", ha="center",
                fontweight="bold", color=C_SLM, fontsize=13)

    axs[0].plot(t, aod_c, color=C_AOD)
    axs[0].set(ylabel="AOD center (um)", ylim=(-0.1, 2.8))
    axs[0].axhline(0, color=C_SLM, ls=":", lw=2, label="SLM center (target)")
    axs[0].legend(loc="upper right")

    axs[1].plot(t, aod_d, color=C_AOD)
    axs[1].set(ylabel="AOD depth (uK)", ylim=(-15, 310))

    axs[2].plot(t, slm_d, color=C_SLM)
    axs[2].set(ylabel="SLM depth (uK)", ylim=(130, 150))
    axs[2].set_xlabel("Time (us)")

    fig.suptitle("Step 1 - Engineered transfer waveform (smoothstep move, square ramp-down)",
                 fontsize=17, fontweight="bold", y=0.995)
    save(fig, "02_waveform_schedule.png")


# ----------------------------------------------------------------------------
# 03 - Potential evolution (3 snapshots)
# ----------------------------------------------------------------------------
def fig_03_potential_evolution():
    x = np.linspace(-3, 5, 1201)
    ideal = load_ideal()

    def total_at(t_us):
        aod_c = np.interp(t_us, ideal["time_us"], ideal["aod_center_um"])
        aod_d = np.interp(t_us, ideal["time_us"], ideal["aod_depth_uK"])
        slm_d = np.interp(t_us, ideal["time_us"], ideal["slm_depth_uK"])
        return (gaussian_uK(x, aod_d, 1.17, aod_c)
                + gaussian_uK(x, slm_d, 1.17, 0.0))

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    snaps = [(0.0, "t = 0 us\nAOD carrying atom"), (312.0, "t = 312 us\nTraps merged"),
             (600.0, "t = 600 us\nSLM-only")]
    for ax, (t_us, title) in zip(axs, snaps):
        u = total_at(t_us)
        ax.plot(x, u, color=C_NEUTRAL, lw=3)
        ax.fill_between(x, u, 0, color=C_NEUTRAL, alpha=0.15)
        ax.set(xlabel="x (um)", title=title, xlim=(-3, 5))
        ax.axhline(0, color="black", lw=0.8)
    axs[0].set_ylabel("U / kB (uK)", fontsize=15)
    fig.suptitle("Step 2 - Live potential during the transfer",
                 fontsize=17, fontweight="bold")
    save(fig, "03_potential_evolution.png")


# ----------------------------------------------------------------------------
# 04 - Ideal atom position
# ----------------------------------------------------------------------------
def fig_04_ideal_position():
    d = load_ideal()
    t = d["time_us"]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(t, d["aod_center_um"], color=C_AOD, label="AOD center (moving)")
    ax.plot(t, d["x_um"], color=C_ATOM, label="Atom position")
    ax.axhline(0.0, color=C_SLM, ls="--", lw=2.5, label="SLM center (x=0)")
    ax.axvline(312, color=C_NEUTRAL, ls=":", lw=2)
    ax.text(312, 2.45, " move ends ", color=C_NEUTRAL, fontsize=12,
            fontweight="bold", ha="center")

    ax.set(xlabel="Time (us)", ylabel="Position (um)",
           title="Step 3 - Ideal atom rides the moving AOD into the SLM",
           xlim=(0, 600), ylim=(-0.2, 2.7))
    ax.legend(loc="upper right")
    save(fig, "04_ideal_position.png")


# ----------------------------------------------------------------------------
# 05 - Ideal mechanical energy (time-dependent -> not conserved)
# ----------------------------------------------------------------------------
def fig_05_ideal_energy():
    d = load_ideal()
    t = d["time_us"]

    fig, axs = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axs[0].plot(t, d["kinetic_uK"], color=C_HIGHLIGHT, label="Kinetic")
    axs[0].plot(t, d["total_potential_uK"], color=C_AOD, label="Potential (AOD+SLM)")
    axs[0].plot(t, d["total_energy_uK"], color="black", lw=3.5, label="Total")
    axs[0].set(ylabel="Energy / kB (uK)", ylim=(-300, 30))
    axs[0].legend(loc="upper right")
    axs[0].set_title("Ideal atom: energy components (time-dependent system)")

    # zoom on total energy
    axs[1].plot(t, d["total_energy_uK"], color="black", lw=3.5)
    axs[1].set(xlabel="Time (us)", ylabel="Total energy (uK)",
               title="Zoom: total energy is NOT conserved - work is done on the atom")
    axs[1].axvline(312, color=C_NEUTRAL, ls=":", lw=2)
    save(fig, "05_ideal_energy.png")


# ----------------------------------------------------------------------------
# 06 - Thermal initial-state distribution
# ----------------------------------------------------------------------------
def fig_06_thermal_distribution():
    mc = load_mc()
    x0 = np.array([r["x0_um"] for r in mc])
    v0 = np.array([r["v0_m_s"] for r in mc])

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    axs[0].hist(x0, bins=10, color=C_AOD, alpha=0.75, edgecolor="white")
    axs[0].axvline(2.4, color="black", ls="--", lw=2, label="AOD center")
    axs[0].set(xlabel="Initial position x0 (um)", ylabel="Shots",
               title="Initial position spread (sigma_x ~ 78 nm @ 5 uK)")
    axs[0].legend()

    axs[1].hist(v0, bins=10, color=C_HIGHLIGHT, alpha=0.75, edgecolor="white")
    axs[1].axvline(0, color="black", ls="--", lw=2)
    axs[1].set(xlabel="Initial velocity v0 (m/s)", ylabel="Shots",
               title="Initial velocity spread (sigma_v ~ 0.018 m/s)")
    fig.suptitle("Step 4 - Sampling bound thermal initial states (5 uK, 50 shots)",
                 fontsize=16, fontweight="bold")
    save(fig, "06_thermal_distribution.png")


# ----------------------------------------------------------------------------
# 07 - Example thermal trajectory
# ----------------------------------------------------------------------------
def fig_07_thermal_position():
    d = load_thermal()
    t = d["time_us"]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(t, d["aod_center_um"], color=C_AOD, label="AOD center")
    ax.plot(t, d["x_um"], color=C_ATOM, label="Thermal atom (example)")
    ax.axhline(0.0, color=C_SLM, ls="--", lw=2.5, label="SLM center")
    ax.axvline(312, color=C_NEUTRAL, ls=":", lw=2)
    ax.set(xlabel="Time (us)", ylabel="Position (um)",
           title="Step 5 - Example thermal trajectory (still captured)",
           xlim=(0, 600))
    ax.legend(loc="upper right")
    save(fig, "07_thermal_position.png")


# ----------------------------------------------------------------------------
# 08 - Phase space: initial vs final
# ----------------------------------------------------------------------------
def fig_08_phase_space():
    mc = load_mc()
    x0 = np.array([r["x0_um"] for r in mc])
    v0 = np.array([r["v0_m_s"] for r in mc])
    xf = np.array([r["final_x_um"] for r in mc])
    vf = np.array([r["final_v_m_s"] for r in mc])

    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    axs[0].scatter(x0, v0, s=70, color=C_AOD, edgecolor="black", linewidth=0.8)
    axs[0].axhline(0, color="black", lw=0.8)
    axs[0].axvline(2.4, color=C_AOD, ls="--", lw=2)
    axs[0].set(xlabel="x (um)", ylabel="v (m/s)",
               title="INITIAL: spread around AOD @ 2.4 um")

    axs[1].scatter(xf, vf, s=70, color=C_SLM, edgecolor="black", linewidth=0.8)
    axs[1].axhline(0, color="black", lw=0.8)
    axs[1].axvline(0, color=C_SLM, ls="--", lw=2)
    axs[1].set(xlabel="x (um)", ylabel="v (m/s)",
               title="FINAL: clustered around SLM @ 0 (focused + cooled)")

    # keep same y-scale to make focusing visible
    ypad = 0.02
    ymin = min(v0.min(), vf.min()) - ypad
    ymax = max(v0.max(), vf.max()) + ypad
    for ax in axs:
        ax.set_ylim(ymin, ymax)
    fig.suptitle("Step 6 - Phase space collapses: 50 MC shots before vs after",
                 fontsize=16, fontweight="bold")
    save(fig, "08_phase_space_collapse.png")


# ----------------------------------------------------------------------------
# 09 - Per-shot capture margin
# ----------------------------------------------------------------------------
def fig_09_capture_margin():
    mc = load_mc()
    margins = np.array([r["capture_margin"] for r in mc])
    order = np.argsort(margins)
    m_sorted = margins[order]

    fig, ax = plt.subplots(figsize=(11, 6))
    colors = [C_SLM if m > 0 else C_THRESHOLD for m in m_sorted]
    ax.bar(range(len(m_sorted)), m_sorted, color=colors, edgecolor="white")
    ax.axhline(0, color=C_THRESHOLD, ls="--", lw=2.5,
               label="Escape threshold (margin=0)")
    ax.axhline(1.0, color=C_NEUTRAL, ls=":", lw=2,
               label="Ideal margin (E = -U_slm)")
    ax.set(xlabel="Shot (sorted by margin)", ylabel="Capture margin (-E_final / U_slm)",
           title=f"Step 7 - Per-shot capture margin (min = {m_sorted.min():.3f}, all captured)")
    ax.legend(loc="lower right")
    save(fig, "09_capture_margin_per_shot.png")


# ----------------------------------------------------------------------------
# 10 - Per-shot excitation energy change
# ----------------------------------------------------------------------------
def fig_10_excitation_change():
    mc = load_mc()
    dE = np.array([r["excitation_energy_change_uK"] for r in mc])

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.hist(dE, bins=12, color=C_HIGHLIGHT, alpha=0.8, edgecolor="white")
    ax.axvline(0, color="black", ls="--", lw=2.5,
               label="No heating / cooling")
    ax.axvline(dE.mean(), color=C_SLM, lw=3,
               label=f"Mean = {dE.mean():.2f} uK (net COOLING)")
    ax.set(xlabel="Excitation-energy change (uK)", ylabel="Shots",
           title="Step 8 - Motional excitation get colder on average")
    ax.legend(loc="upper left")
    save(fig, "10_excitation_change_hist.png")


# ----------------------------------------------------------------------------
# 11 - Capture-rate summary with Wilson CI
# ----------------------------------------------------------------------------
def fig_11_capture_summary():
    m = json.loads((DATA / "metrics.json").read_text())
    mc = m["monte_carlo"]
    captured = mc["captured_shots"]
    total = mc["shots"]
    rate = mc["classical_capture_fraction"]
    lo, hi = mc["wilson_95_interval"]

    fig, ax = plt.subplots(figsize=(11, 6))
    # bar
    ax.barh(["Classical\ncapture fraction"], [rate], color=C_ATOM,
            height=0.4, edgecolor="black", linewidth=1.5)
    # CI error bar
    err_lo = rate - lo
    err_hi = hi - rate
    ax.errorbar([rate], [0], xerr=[[err_lo], [err_hi]], fmt="none",
                ecolor="black", elinewidth=2.5, capsize=14, capthick=2.5,
                label=f"Wilson 95% CI: [{lo:.4f}, {hi:.4f}]")
    ax.axvline(1.0, color=C_NEUTRAL, ls=":", lw=2)

    ax.set_xlim(0.85, 1.02)
    ax.set_xlabel("Capture fraction")
    ax.set_title(f"Step 9 - Result: {captured}/{total} captured "
                 f"(rate = {rate:.3f})", fontsize=18)
    ax.legend(loc="lower right", fontsize=14)

    # caveat as subtitle
    fig.text(0.5, -0.02,
             "Classical mechanical capture - NOT experimental survival / "
             "transfer fidelity / quantum fidelity.",
             ha="center", fontsize=11, style="italic", color=C_NEUTRAL)
    save(fig, "11_capture_rate_summary.png")


# ----------------------------------------------------------------------------
# 12 - Validation dashboard
# ----------------------------------------------------------------------------
def fig_12_validation():
    m = json.loads((DATA / "metrics.json").read_text())
    v = m["validation"]
    checks = [
        ("Level 0 baseline", v["baseline_level0_passed"]),
        ("Time-dependent VV\nstatic-limit agrees", v["time_dependent_integrator_static_limit_passed"]),
        ("Waveform\nendpoints", v["waveform_endpoints_passed"]),
        ("Ideal transfer\ndt convergence", v["ideal_transfer_dt_convergence_passed"]),
        ("Ideal capture\nreported", v["ideal_transfer_capture_reported"]),
        ("Thermal example\ncompleted", v["thermal_example_completed"]),
        ("Monte Carlo\nreproducibility", v["monte_carlo_reproducibility_passed"]),
        ("ALL CHECKS", v["all_passed"]),
    ]
    labels = [c[0] for c in checks]
    passed = [c[1] for c in checks]
    colors = [C_ATOM if p else C_THRESHOLD for p in passed]

    fig, ax = plt.subplots(figsize=(13, 6.5))
    y = np.arange(len(labels))[::-1]
    ax.barh(y, [1] * len(labels), color=colors, edgecolor="black",
            linewidth=1.5, height=0.6)
    for yi, p in zip(y, passed):
        ax.text(0.5, yi, "PASS" if p else "FAIL", ha="center", va="center",
                color="white", fontsize=15, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("Step 10 - Validation dashboard (all green)", fontsize=18)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    save(fig, "12_validation_dashboard.png")


# ----------------------------------------------------------------------------
def main():
    print(f"Reading from: {DATA}")
    print(f"Writing to:   {OUT_DIR}\n")
    fig_01_setup()
    fig_02_waveform()
    fig_03_potential_evolution()
    fig_04_ideal_position()
    fig_05_ideal_energy()
    fig_06_thermal_distribution()
    fig_07_thermal_position()
    fig_08_phase_space()
    fig_09_capture_margin()
    fig_10_excitation_change()
    fig_11_capture_summary()
    fig_12_validation()
    print(f"\nDone. {len(list(OUT_DIR.glob('*.png')))} figures in {OUT_DIR}")


if __name__ == "__main__":
    main()
