"""Temperature scan: how capture fraction varies with atom temperature.

Sweeps atom temperature across the physically interesting range and runs a
Monte Carlo (1000 shots per temperature) at each point to map the capture
fraction. Parallelized across CPU cores. Writes raw per-shot data, an
aggregate summary, and four PPT-ready figures.

Physics intuition for the sweep range:
  - Baseline operating point: 5 uK  -> near-100% capture (this is what the
    Level 1 demo reports).
  - SLM depth is 140 uK. Once the atom's thermal energy is comparable to the
    SLM well depth, a meaningful fraction of initial states are no longer
    capturable -> expect the curve to roll off somewhere around the SLM
    depth. The scan is densified there to resolve the onset.

Run from this directory:
    python run_temperature_scan.py
"""
from __future__ import annotations

import csv
import json
import os
import time
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from level1_transfer_1d.analysis import wilson_interval
from level1_transfer_1d.config import load_config
from level1_transfer_1d.simulation import monte_carlo

# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
CONFIG = HERE.parents[0] / "configs" / "transfer_1d.yaml"
OUT_DIR = HERE / "scan_temperature"
OUT_DIR.mkdir(exist_ok=True)

NPROC = max(1, (os.cpu_count() or 2) - 1)  # leave one core for the OS
SHOTS_PER_TEMP = 1000
CHUNKS = NPROC  # split each temperature's shots into NPROC chunks
BASE_SEED = 20260811

# Presentation style (matches make_ppt_figures.py)
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "font.size": 15,
    "axes.titlesize": 18, "axes.titleweight": "bold",
    "axes.labelsize": 15, "axes.labelweight": "bold",
    "xtick.labelsize": 13, "ytick.labelsize": 13,
    "legend.fontsize": 13, "legend.frameon": True,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--",
    "axes.linewidth": 1.5, "lines.linewidth": 3.0,
    "axes.spines.top": False, "axes.spines.right": False,
})
C_AOD = "#1f6fb4"; C_SLM = "#d62728"; C_ATOM = "#2ca02c"
C_HL = "#ff7f0e"; C_NEUTRAL = "#7f7f7f"


# ----------------------------------------------------------------------------
# Temperature grid: densified in the transition region (40-160 uK)
# ----------------------------------------------------------------------------
TEMPS = [5, 10, 20, 30, 40, 50, 60, 70, 80, 95, 110, 125, 140, 160, 200, 250]


# ----------------------------------------------------------------------------
# Parallel worker: run one chunk of MC at a fixed temperature.
# Must be top-level for macOS spawn-based pickling.
# ----------------------------------------------------------------------------
def run_chunk(args):
    temp_uK, shots, seed = args
    cfg = replace(
        load_config(CONFIG),
        thermal=replace(load_config(CONFIG).thermal, temperature_uK=temp_uK),
    )
    return monte_carlo(cfg, shots=shots, seed=seed)


def run_one_temperature(temp_uK: float, idx: int):
    """Run SHOTS_PER_TEMP shots at one temperature, in parallel chunks."""
    per_chunk = SHOTS_PER_TEMP // CHUNKS
    jobs = [(temp_uK, per_chunk, BASE_SEED + idx * 1000 + c) for c in range(CHUNKS)]
    with Pool(NPROC) as pool:
        chunks = pool.map(run_chunk, jobs)
    # renumber shot ids globally
    records = []
    sid = 0
    for chunk in chunks:
        for r in chunk:
            r["shot_id"] = sid
            r["temperature_uK"] = temp_uK
            records.append(r)
            sid += 1
    return records


def aggregate(records):
    n = len(records)
    cap = sum(r["captured_in_slm"] for r in records)
    rate = cap / n
    lo, hi = wilson_interval(cap, n)
    margins = [r["capture_margin"] for r in records if r["captured_in_slm"]]
    dE = [r["excitation_energy_change_uK"] for r in records]
    return {
        "temperature_uK": records[0]["temperature_uK"],
        "shots": n,
        "captured": cap,
        "capture_fraction": rate,
        "wilson_lo": lo,
        "wilson_hi": hi,
        "mean_capture_margin": float(np.mean(margins)) if margins else None,
        "median_capture_margin": float(np.median(margins)) if margins else None,
        "mean_excitation_change_uK": float(np.mean(dE)),
        "median_excitation_change_uK": float(np.median(dE)),
    }


# ----------------------------------------------------------------------------
def main():
    print(f"Temperature scan: {len(TEMPS)} points x {SHOTS_PER_TEMP} shots "
          f"= {len(TEMPS)*SHOTS_PER_TEMP} total")
    print(f"Parallel: {NPROC} procs, {SHOTS_PER_TEMP//CHUNKS} shots/chunk")
    print(f"Output: {OUT_DIR}\n")
    print(f"{'T(uK)':>7} {'cap/n':>10} {'rate':>7} {'Wilson95% CI':>20}  elapsed")
    print("-" * 64)

    all_records = []
    summary = []
    t0 = time.time()
    for i, t in enumerate(TEMPS):
        recs = run_one_temperature(t, i)
        all_records.extend(recs)
        s = aggregate(recs)
        summary.append(s)
        print(f"{t:>7.1f} {s['captured']:>4}/{s['shots']:<4d} {s['capture_fraction']:>7.3f}"
              f"  [{s['wilson_lo']:.3f}, {s['wilson_hi']:.3f}]   {time.time()-t0:5.1f}s")

    elapsed = time.time() - t0
    total_shots = len(all_records)
    print("-" * 64)
    print(f"Done: {total_shots} shots in {elapsed:.1f}s "
          f"({1000*elapsed/total_shots:.1f} ms/shot effective)")

    # ---- save raw per-shot data ----
    raw_path = OUT_DIR / "mc_shots_by_temperature.csv"
    with open(raw_path, "w", newline="") as f:
        keys = ["temperature_uK"] + [k for k in all_records[0] if k != "temperature_uK"]
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in all_records:
            w.writerow({k: ("" if r.get(k) is None else r[k]) for k in keys})
    print(f"raw  -> {raw_path.name}  ({total_shots} rows)")

    # ---- save summary ----
    sum_path = OUT_DIR / "summary_by_temperature.csv"
    with open(sum_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print(f"summ -> {sum_path.name}")

    meta = {
        "shots_per_temperature": SHOTS_PER_TEMP,
        "n_temperatures": len(TEMPS),
        "total_shots": total_shots,
        "elapsed_s": elapsed,
        "ms_per_shot_effective": 1000 * elapsed / total_shots,
        "n_procs": NPROC,
        "base_seed": BASE_SEED,
        "slm_depth_uK": load_config(CONFIG).slm.depth_uK,
        "aod_depth_uK": load_config(CONFIG).aod.depth_uK,
    }
    (OUT_DIR / "scan_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"meta -> scan_meta.json")

    # ---- figures ----
    make_figures(summary, all_records, meta)
    print(f"\nAll figures written to {OUT_DIR}")


# ----------------------------------------------------------------------------
def make_figures(summary, all_records, meta):
    T = np.array([s["temperature_uK"] for s in summary])
    rate = np.array([s["capture_fraction"] for s in summary])
    lo = np.array([s["wilson_lo"] for s in summary])
    hi = np.array([s["wilson_hi"] for s in summary])
    slm_d = meta["slm_depth_uK"]
    aod_d = meta["aod_depth_uK"]

    def save(fig, name):
        path = OUT_DIR / name
        fig.tight_layout()
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"  fig -> {name}")

    # clip tiny negative widths caused by float roundoff at rate=1.0
    err_lo = np.clip(rate - lo, 0.0, None)
    err_hi = np.clip(hi - rate, 0.0, None)

    # 1) capture fraction vs temperature (the headline figure)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.errorbar(T, rate, yerr=[err_lo, err_hi], fmt="-o", color=C_ATOM,
                ecolor=C_NEUTRAL, elinewidth=2, capsize=6, capthick=2,
                markersize=9, markeredgecolor="black", markeredgewidth=1,
                label="Capture fraction (1000 shots/point, Wilson 95% CI)")
    ax.axvline(slm_d, color=C_SLM, ls="--", lw=2.5,
               label=f"SLM depth = {slm_d:.0f} uK")
    ax.axvline(aod_d, color=C_AOD, ls="--", lw=2.5,
               label=f"AOD depth = {aod_d:.0f} uK")
    ax.axvspan(5, 20, color=C_ATOM, alpha=0.08)
    ax.text(12, 0.06, "baseline\noperating zone", ha="center", fontsize=11,
            color=C_ATOM, fontweight="bold")
    ax.set(xlabel="Atom temperature (uK)", ylabel="Classical capture fraction",
           title="Capture fraction vs atom temperature",
           ylim=(0, 1.08), xlim=(0, max(T) + 10))
    ax.legend(loc="lower left")
    save(fig, "capture_fraction_vs_temperature.png")

    # 2) Wilson CI band (uncertainty view, log-x to stretch the low-T regime)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.fill_between(T, lo, hi, color=C_ATOM, alpha=0.25, label="Wilson 95% CI")
    ax.plot(T, rate, "o-", color=C_ATOM, markersize=9,
            markeredgecolor="black", markeredgewidth=1, label="Point estimate")
    ax.axvline(slm_d, color=C_SLM, ls="--", lw=2.5,
               label=f"SLM depth = {slm_d:.0f} uK")
    ax.set_xscale("log")
    ax.set(xlabel="Atom temperature (uK, log scale)",
           ylabel="Capture fraction",
           title="Statistical uncertainty across the temperature scan",
           ylim=(0, 1.05))
    ax.legend(loc="lower left")
    save(fig, "capture_fraction_uncertainty.png")

    # 3) excitation-energy change vs temperature
    mean_dE = np.array([s["mean_excitation_change_uK"] for s in summary])
    med_dE = np.array([s["median_excitation_change_uK"] for s in summary])
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(T, mean_dE, "o-", color=C_HL, markersize=9,
            markeredgecolor="black", markeredgewidth=1, label="Mean")
    ax.plot(T, med_dE, "s--", color=C_AOD, markersize=8,
            markeredgecolor="black", markeredgewidth=1, label="Median")
    ax.axhline(0, color="black", ls="--", lw=2)
    ax.axvline(slm_d, color=C_SLM, ls="--", lw=2.5,
               label=f"SLM depth = {slm_d:.0f} uK")
    ax.set(xlabel="Atom temperature (uK)", ylabel="Excitation change (uK)",
           title="Motional excitation vs temperature\n(negative = net cooling)",
           xlim=(0, max(T) + 10))
    ax.legend(loc="best")
    save(fig, "excitation_change_vs_temperature.png")

    # 4) capture-margin distribution shift (box plot per temperature)
    by_T = {}
    for r in all_records:
        by_T.setdefault(r["temperature_uK"], []).append(
            r["capture_margin"] if r["capture_margin"] is not None else -0.05
        )
    data = [by_T[t] for t in T]
    fig, ax = plt.subplots(figsize=(12, 6.5))
    bp = ax.boxplot(data, positions=range(len(T)), widths=0.6, patch_artist=True,
                    showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(C_ATOM); patch.set_alpha(0.5)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(2)
    ax.axhline(0, color=C_SLM, ls="--", lw=2, label="Escape threshold (margin=0)")
    ax.axhline(1, color=C_NEUTRAL, ls=":", lw=2, label="Ideal (E=-U_slm)")
    ax.set_xticks(range(len(T)))
    ax.set_xticklabels([f"{t:g}" for t in T])
    ax.set(xlabel="Atom temperature (uK)", ylabel="Capture margin (-E_final/U_slm)",
           title="Capture-margin distribution shifts as temperature rises",
           ylim=(-0.15, 1.1))
    ax.legend(loc="lower left")
    save(fig, "capture_margin_box_by_temperature.png")


if __name__ == "__main__":
    main()
