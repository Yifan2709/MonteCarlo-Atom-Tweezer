"""Shared framework for Level 1 parameter scans.

Extends the pattern of run_temperature_scan.py: every scan sweeps one (or
two) parameters, runs a Monte Carlo at each point with a single shared
multiprocessing Pool, and writes raw per-shot CSV, an aggregate summary CSV,
a meta JSON, and a standard multi-panel detail figure.

Scan modes:
  - thermal          : finite-temperature sampling (default, like the demo)
  - fixed_excitation : deterministic initial excitation energy, position
                       sampled uniformly between the turning points of the
                       static AOD well, |v| from energy conservation, random
                       sign. NOTE: uniform-in-x on the energy shell, not a
                       microcanonical (uniform-in-phase) distribution.

Waveform variants (for waveform-shape / ramp-exponent scans) replace the
hardcoded smoothstep / quadratic drop-off in level1_transfer_1d.waveforms.
The baseline "smoothstep" + exponent 2.0 reproduces the library waveform
exactly (verified in test / smoke checks).

All scans keep the Level 1 honesty boundary: results are classical mechanics
for this configured model, not experimental survival or transfer fidelity.
"""
from __future__ import annotations

import csv
import json
import math
import os
import time
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from level0_static_trap.constants import BOLTZMANN_CONSTANT
from level0_static_trap.potentials import gaussian_force, gaussian_potential
from level1_transfer_1d.analysis import classify_final_state, wilson_interval
from level1_transfer_1d.config import load_config
from level1_transfer_1d.integrator import velocity_verlet_time_dependent
from level1_transfer_1d.simulation import time_grid
from level1_transfer_1d.thermal import sample_thermal_initial_state

J_TO_UK = 1e6 / BOLTZMANN_CONSTANT

HERE = Path(__file__).resolve().parent
CONFIG = HERE.parents[0] / "configs" / "transfer_1d.yaml"
NPROC = max(1, (os.cpu_count() or 2) - 1)
BASE_SEED = 20260811

# Presentation style (matches make_ppt_figures.py / run_temperature_scan.py)
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "font.size": 14,
    "axes.titlesize": 16, "axes.titleweight": "bold",
    "axes.labelsize": 14, "axes.labelweight": "bold",
    "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 11, "legend.frameon": True,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
    "axes.linewidth": 1.5, "lines.linewidth": 2.8,
    "axes.spines.top": False, "axes.spines.right": False,
})
C_AOD = "#1f6fb4"; C_SLM = "#d62728"; C_ATOM = "#2ca02c"
C_HL = "#ff7f0e"; C_NEUTRAL = "#7f7f7f"; C_ESC = "#9467bd"

# ----------------------------------------------------------------------------
# Waveform variants
# ----------------------------------------------------------------------------
MOVE_SHAPES = {
    "linear":       lambda s: s,
    "smoothstep":   lambda s: 3 * s**2 - 2 * s**3,          # baseline
    "smootherstep": lambda s: 6 * s**5 - 15 * s**4 + 10 * s**3,
    "sine":         lambda s: 0.5 * (1.0 - math.cos(math.pi * s)),
    "ease_in_quad": lambda s: s**2,
}


def move_h(name: str, s: float) -> float:
    return MOVE_SHAPES[name](s)


def make_waveform(config, move_shape="smoothstep", ramp_exponent=2.0):
    """Return (center_fn(t), depth_fn(t)) for a configurable waveform.

    Baseline (smoothstep, exponent 2.0) is identical to
    level1_transfer_1d.waveforms.aod_center / aod_depth.
    """
    T = config.transfer.duration_s
    tm = config.transfer.move_duration_s
    x0 = config.aod.initial_center_m
    U0 = config.aod.depth_j

    def center(t):
        t = min(max(float(t), 0.0), T)
        s = min(t / tm, 1.0)
        return x0 * (1.0 - move_h(move_shape, s))

    def depth(t):
        t = min(max(float(t), 0.0), T)
        s = min(max((t - tm) / (T - tm), 0.0), 1.0)
        return U0 * (1.0 - s) ** ramp_exponent

    return center, depth


def make_force(config, move_shape="smoothstep", ramp_exponent=2.0):
    center_fn, depth_fn = make_waveform(config, move_shape, ramp_exponent)
    sd, sw, sc = config.slm.depth_j, config.slm.waist_m, config.slm.center_m
    aw = config.aod.waist_m

    def force(x, t):
        return (gaussian_force(x, sd, sw, sc)
                + gaussian_force(x, depth_fn(t), aw, center_fn(t)))
    return force


# ----------------------------------------------------------------------------
# Initial-state samplers
# ----------------------------------------------------------------------------
def turning_half_width(config, excitation_uK):
    """Half-width of the classically allowed region in the static AOD well
    for a given excitation energy above the well bottom."""
    U0 = config.aod.depth_j
    w = config.aod.waist_m
    exc_J = excitation_uK / J_TO_UK
    E = -U0 + exc_J
    if E >= 0.0:
        raise ValueError("excitation meets or exceeds the AOD well depth")
    return w * math.sqrt(-math.log(-E / U0) / 2.0)


def sample_fixed_excitation(config, rng, excitation_uK):
    """One initial state at a deterministic excitation energy.

    x is uniform between the turning points of the static AOD well,
    |v| follows from energy conservation, sign random.
    Returns (x0, v0, E_total_J, excitation_J, half_width).
    """
    U0 = config.aod.depth_j
    w = config.aod.waist_m
    xc = config.aod.initial_center_m
    exc_J = excitation_uK / J_TO_UK
    E = -U0 + exc_J
    dx = turning_half_width(config, excitation_uK)
    x = xc + rng.uniform(-dx, dx)
    U = float(gaussian_potential(x, U0, w, xc))
    KE = max(E - U, 0.0)
    v = math.sqrt(2.0 * KE / config.atom.mass_kg) * float(rng.choice([-1.0, 1.0]))
    return x, v, E, exc_J, dx


# ----------------------------------------------------------------------------
# Core shot runner (single process) and Pool worker
# ----------------------------------------------------------------------------
def run_shots(config, shots, seed, mode="thermal", fixed_excitation_uK=None,
              move_shape="smoothstep", ramp_exponent=2.0, dt_s=None):
    rng = np.random.default_rng(seed)
    force = make_force(config, move_shape, ramp_exponent)
    times = time_grid(config, dt_s)
    mass = config.atom.mass_kg
    records = []
    for i in range(shots):
        if mode == "thermal":
            st = sample_thermal_initial_state(config, rng)
            x0, v0 = st.x0_m, st.v0_m_s
            E0, exc0 = st.aod_total_energy_J, st.aod_excitation_energy_J
            x_turn = None
        elif mode == "fixed_excitation":
            x0, v0, E0, exc0, dx = sample_fixed_excitation(
                config, rng, fixed_excitation_uK)
            x_turn = (x0 - config.aod.initial_center_m) / dx  # in [-1, 1]
        else:
            raise ValueError(f"unknown mode {mode}")
        _, x, v, _ = velocity_verlet_time_dependent(force, mass, x0, v0, times)
        fin = classify_final_state(float(x[-1]), float(v[-1]), config)
        records.append({
            "shot_id": i,
            "x0_um": x0 * 1e6, "v0_m_s": v0,
            "initial_aod_energy_uK": E0 * J_TO_UK,
            "initial_excitation_energy_uK": exc0 * J_TO_UK,
            "initial_phase_coord": x_turn,
            "final_x_um": float(x[-1]) * 1e6, "final_v_m_s": float(v[-1]),
            "final_slm_energy_uK": fin["final_slm_energy_J"] * J_TO_UK,
            "final_excitation_energy_uK": fin["final_excitation_energy_J"] * J_TO_UK,
            "excitation_energy_change_uK":
                (fin["final_excitation_energy_J"] - exc0) * J_TO_UK,
            "final_energy_over_slm_depth": fin["final_energy_over_slm_depth"],
            "captured_in_slm": fin["captured_in_slm"],
            "capture_margin": fin["capture_margin"],
        })
    return records


def ideal_shot(config, move_shape="smoothstep", ramp_exponent=2.0, dt_s=None):
    """Deterministic atom starting at rest at the AOD centre."""
    force = make_force(config, move_shape, ramp_exponent)
    times = time_grid(config, dt_s)
    _, x, v, _ = velocity_verlet_time_dependent(
        force, config.atom.mass_kg, config.aod.initial_center_m, 0.0, times)
    fin = classify_final_state(float(x[-1]), float(v[-1]), config)
    return {
        "captured_in_slm": fin["captured_in_slm"],
        "final_energy_over_slm_depth": fin["final_energy_over_slm_depth"],
        "final_excitation_energy_uK": fin["final_excitation_energy_J"] * J_TO_UK,
        "excitation_energy_change_uK": fin["final_excitation_energy_J"] * J_TO_UK,
    }


def worker(job):
    """Top-level Pool worker (macOS spawn-safe). Runs one chunk of shots."""
    cfg = load_config(CONFIG)
    for sec, kv in job.get("overrides", {}).items():
        cfg = replace(cfg, **{sec: replace(getattr(cfg, sec), **kv)})
    dt_s = job["dt_us"] * 1e-6 if job.get("dt_us") else None
    recs = run_shots(
        cfg, job["shots"], job["seed"],
        mode=job.get("mode", "thermal"),
        fixed_excitation_uK=job.get("fixed_excitation_uK"),
        move_shape=job.get("move_shape", "smoothstep"),
        ramp_exponent=job.get("ramp_exponent", 2.0),
        dt_s=dt_s)
    for r in recs:
        r["label"] = job["label"]
    return recs


def ideal_worker(job):
    cfg = load_config(CONFIG)
    for sec, kv in job.get("overrides", {}).items():
        cfg = replace(cfg, **{sec: replace(getattr(cfg, sec), **kv)})
    dt_s = job["dt_us"] * 1e-6 if job.get("dt_us") else None
    out = ideal_shot(cfg, move_shape=job.get("move_shape", "smoothstep"),
                     ramp_exponent=job.get("ramp_exponent", 2.0), dt_s=dt_s)
    out["label"] = job["label"]
    return out


# ----------------------------------------------------------------------------
# Scan driver
# ----------------------------------------------------------------------------
def run_scan(scan_name, jobs, ideal_jobs=None, x_field=None, out_dir=None):
    """Run all jobs in ONE shared Pool; group records by label.

    jobs: list of worker job dicts (already chunked, with unique seeds).
    Returns (records_by_label, ideal_by_label, meta_dict).
    Writes raw CSV / ideal CSV / meta JSON incrementally.
    """
    out_dir = Path(out_dir or (HERE / scan_name))
    out_dir.mkdir(exist_ok=True)
    t0 = time.time()
    print(f"[{scan_name}] {len(jobs)} chunk-jobs on {NPROC} procs -> {out_dir.name}/")
    with Pool(NPROC) as pool:
        chunked = pool.map(worker, jobs)
        ideals = pool.map(ideal_worker, ideal_jobs) if ideal_jobs else []
    records_by_label = {}
    for chunk in chunked:
        if not chunk:
            continue
        records_by_label.setdefault(chunk[0]["label"], []).extend(chunk)
    ideal_by_label = {r["label"]: r for r in ideals}
    meta = {
        "scan": scan_name, "x_field": x_field,
        "n_points": len(records_by_label),
        "total_shots": sum(len(v) for v in records_by_label.values()),
        "elapsed_s": time.time() - t0, "n_procs": NPROC, "base_seed": BASE_SEED,
        "slm_depth_uK": load_config(CONFIG).slm.depth_uK,
        "aod_depth_uK": load_config(CONFIG).aod.depth_uK,
        "temperature_uK": load_config(CONFIG).thermal.temperature_uK,
        "transfer_duration_us": load_config(CONFIG).transfer.duration_us,
    }
    # raw per-shot CSV
    all_recs = [r for v in records_by_label.values() for r in v]
    if all_recs:
        keys = list(all_recs[0].keys())
        with open(out_dir / "mc_shots.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in all_recs:
                w.writerow({k: ("" if r.get(k) is None else r[k]) for k in keys})
    if ideal_by_label:
        keys = list(next(iter(ideal_by_label.values())).keys())
        with open(out_dir / "ideal.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for lbl in sorted(ideal_by_label):
                w.writerow(ideal_by_label[lbl])
    (out_dir / "scan_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{scan_name}] done: {meta['total_shots']} shots in "
          f"{meta['elapsed_s']:.1f}s")
    return records_by_label, ideal_by_label, meta


def chunked_jobs(points, shots_per_point, job_factory, chunks=NPROC):
    """Build chunk-level jobs. job_factory(point, label) -> partial job dict
    (everything except shots/seed/label)."""
    jobs = []
    per_chunk = max(1, shots_per_point // chunks)
    for pi, point in enumerate(points):
        label = job_factory(point, None)["label"]
        for c in range(chunks):
            job = job_factory(point, label)
            job.update({"shots": per_chunk,
                        "seed": BASE_SEED + pi * 1000 + c, "label": label})
            jobs.append(job)
    return jobs


# ----------------------------------------------------------------------------
# Aggregation
# ----------------------------------------------------------------------------
def aggregate(records, x_value, x_field):
    n = len(records)
    cap = sum(r["captured_in_slm"] for r in records)
    lo, hi = wilson_interval(cap, n)
    margins = [r["capture_margin"] for r in records
               if r["capture_margin"] is not None]
    dE = [r["excitation_energy_change_uK"] for r in records]
    fin_exc = [r["final_excitation_energy_uK"] for r in records]
    return {
        x_field: x_value, "shots": n, "captured": cap,
        "capture_fraction": cap / n, "wilson_lo": lo, "wilson_hi": hi,
        "mean_capture_margin": float(np.mean(margins)) if margins else None,
        "median_capture_margin": float(np.median(margins)) if margins else None,
        "mean_excitation_change_uK": float(np.mean(dE)),
        "median_excitation_change_uK": float(np.median(dE)),
        "mean_final_excitation_uK": float(np.mean(fin_exc)),
    }


def write_summary(out_dir, summary_rows):
    out_dir = Path(out_dir)
    with open(out_dir / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        for row in summary_rows:
            w.writerow({k: ("" if row.get(k) is None else row[k])
                        for k in row})
    print(f"  summary -> {out_dir.name}/summary.csv")


# ----------------------------------------------------------------------------
# Standard detail figure (2 rows x 3 cols)
# ----------------------------------------------------------------------------
def detail_figure(out_dir, x, summary, records_by_label, ideal_by_label,
                  x_label, title, x_log=False, vlines=None, filename=None,
                  note_lines=None):
    """Six-panel scan figure:
    1 capture fraction + Wilson CI   2 excitation change (MC + ideal)
    3 capture-margin box plot        4 final-excitation histograms (4 points)
    5 initial vs final excitation scatter (escapes highlighted)
    6 text panel with key numbers
    """
    out_dir = Path(out_dir)
    x = np.asarray(x, dtype=float)
    rate = np.array([s["capture_fraction"] for s in summary])
    lo = np.array([s["wilson_lo"] for s in summary])
    hi = np.array([s["wilson_hi"] for s in summary])
    err_lo = np.clip(rate - lo, 0.0, None)
    err_hi = np.clip(hi - rate, 0.0, None)

    fig, axes = plt.subplots(2, 3, figsize=(19, 10.5))
    fig.suptitle(title, fontsize=19, fontweight="bold")

    # 1) capture fraction
    ax = axes[0, 0]
    ax.errorbar(x, rate, yerr=[err_lo, err_hi], fmt="-o", color=C_ATOM,
                ecolor=C_NEUTRAL, elinewidth=2, capsize=5, capthick=2,
                markersize=8, markeredgecolor="black", markeredgewidth=1,
                label="Capture fraction (Wilson 95% CI)")
    for xv, c, lab in (vlines or []):
        ax.axvline(xv, color=c, ls="--", lw=2, label=lab)
    ax.set(xlabel=x_label, ylabel="Classical capture fraction",
           ylim=(-0.03, 1.09), title="Capture fraction")
    if x_log:
        ax.set_xscale("log")
    ax.legend(loc="lower left", fontsize=10)

    # 2) excitation change
    ax = axes[0, 1]
    mean_dE = [s["mean_excitation_change_uK"] for s in summary]
    med_dE = [s["median_excitation_change_uK"] for s in summary]
    ax.plot(x, mean_dE, "o-", color=C_HL, markersize=7,
            markeredgecolor="black", markeredgewidth=0.8, label="MC mean")
    ax.plot(x, med_dE, "s--", color=C_AOD, markersize=6,
            markeredgecolor="black", markeredgewidth=0.8, label="MC median")
    if ideal_by_label:
        labels = [s.get("__label__") for s in summary]
        if all(l in ideal_by_label for l in labels):
            idE = [ideal_by_label[l]["excitation_energy_change_uK"] for l in labels]
            ax.plot(x, idE, "d-.", color=C_SLM, markersize=6,
                    label="Ideal atom (rest at AOD centre)")
    ax.axhline(0, color="black", ls="--", lw=1.8)
    ax.set(xlabel=x_label, ylabel="Excitation-energy change (uK)",
           title="Motional excitation (negative = net cooling)")
    if x_log:
        ax.set_xscale("log")
    ax.legend(loc="best", fontsize=10)

    # 3) capture-margin box plot
    ax = axes[0, 2]
    labels = [s.get("__label__") for s in summary]
    data = []
    for lbl in labels:
        vals = [(r["capture_margin"] if r["capture_margin"] is not None else -0.05)
                for r in records_by_label[lbl]]
        data.append(vals)
    bp = ax.boxplot(data, positions=range(len(x)), widths=0.6,
                    patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(C_ATOM); patch.set_alpha(0.5)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(2)
    ax.axhline(0, color=C_SLM, ls="--", lw=2, label="Escape threshold")
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels([f"{v:g}" for v in x], rotation=45, ha="right", fontsize=10)
    ax.set(xlabel=x_label, ylabel="Capture margin (-E_final/U_slm)",
           title="Capture-margin distribution", ylim=(-0.15, 1.1))
    ax.legend(loc="lower left", fontsize=10)

    # 4) final-excitation histograms at up to 4 points
    ax = axes[1, 0]
    idxs = np.unique(np.linspace(0, len(x) - 1, min(4, len(x))).astype(int))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(idxs)))
    for i, c in zip(idxs, colors):
        vals = [r["final_excitation_energy_uK"] for r in records_by_label[labels[i]]]
        ax.hist(vals, bins=40, histtype="step", lw=2.2, color=c,
                label=f"{x_label.split('(')[0].strip()} = {x[i]:g}")
    ax.set(xlabel="Final SLM excitation energy (uK)", ylabel="Shots",
           title="Final-energy histograms")
    ax.legend(fontsize=10)

    # 5) initial vs final excitation scatter (subsample, escapes highlighted)
    ax = axes[1, 1]
    all_r = [r for lbl in labels for r in records_by_label[lbl]]
    if len(all_r) > 4000:
        sel = np.random.default_rng(0).choice(len(all_r), 4000, replace=False)
        all_r = [all_r[i] for i in sel]
    ei = [r["initial_excitation_energy_uK"] for r in all_r]
    ef = [r["final_excitation_energy_uK"] for r in all_r]
    capd = np.array([r["captured_in_slm"] for r in all_r])
    ei = np.array(ei); ef = np.array(ef)
    ax.scatter(ei[capd], ef[capd], s=6, alpha=0.25, color=C_ATOM,
               label="captured", rasterized=True)
    ax.scatter(ei[~capd], ef[~capd], s=18, alpha=0.8, color=C_ESC,
               marker="x", label="escaped", rasterized=True)
    ax.axhline(load_config(CONFIG).slm.depth_uK, color=C_SLM, ls="--", lw=2,
               label="SLM depth (escape above)")
    ax.set(xlabel="Initial AOD excitation (uK)",
           ylabel="Final SLM excitation (uK)",
           title="Per-shot energy transfer")
    ax.legend(fontsize=10, loc="upper left")

    # 6) text panel
    ax = axes[1, 2]
    ax.axis("off")
    lines = list(note_lines or [])
    ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes,
            fontsize=11.5, va="top", family="monospace")
    ax.set_title("Key numbers")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    name = filename or "scan_detail.png"
    fig.savefig(out_dir / name, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig -> {out_dir.name}/{name}")
