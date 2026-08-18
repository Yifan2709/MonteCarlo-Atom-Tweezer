"""Scan 14: Monte Carlo statistics (shots x seeds).

Runs the baseline configuration with N in {50..800} shots, repeated with 4
independent seeds. Shows how the capture-fraction estimate and its Wilson
interval converge with N, and the seed-to-seed scatter. Determines the
shots-per-point budget used across the other scans.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

NSHOTS = [50, 100, 200, 400, 800]
SEEDS = [11, 22, 33, 44]
SCAN = "scan_mc_stats"


def main():
    jobs = []
    for pi, n in enumerate(NSHOTS):
        for si, seed in enumerate(SEEDS):
            jobs.append({
                "overrides": {}, "mode": "thermal",
                "shots": n, "seed": sl.BASE_SEED + seed,
                "label": f"N={n},seed={seed}",
            })
    out = sl.HERE / SCAN
    recs, _, meta = sl.run_scan(SCAN, jobs, None, "shots", out)

    rows = []
    for n in NSHOTS:
        for seed in SEEDS:
            lbl = f"N={n},seed={seed}"
            row = sl.aggregate(recs[lbl], n, "shots")
            row["__label__"] = lbl
            row["seed"] = seed
            row["wilson_halfwidth"] = (row["wilson_hi"] - row["wilson_lo"]) / 2
            rows.append(row)
    sl.write_summary(out, rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    ax = axes[0]
    for seed in SEEDS:
        r = [row for row in rows if row["seed"] == seed]
        ax.errorbar([row["shots"] for row in r],
                    [row["capture_fraction"] for row in r],
                    yerr=[row["wilson_halfwidth"] for row in r],
                    fmt="o-", capsize=5, markersize=7, label=f"seed {seed}")
    ax.set_xscale("log")
    ax.set(xlabel="Shots per Monte Carlo", ylabel="Capture fraction",
           ylim=(0.9, 1.05), title="Estimate convergence with N")
    ax.legend(fontsize=11)
    ax = axes[1]
    mean_hw = [np.mean([row["wilson_halfwidth"] for row in rows
                        if row["shots"] == n]) for n in NSHOTS]
    ax.plot(NSHOTS, mean_hw, "o-", color=sl.C_AOD, markersize=9,
            markeredgecolor="black", markeredgewidth=1, label="empirical")
    nn = np.linspace(min(NSHOTS), max(NSHOTS), 100)
    ax.plot(nn, mean_hw[0] * np.sqrt(NSHOTS[0] / nn), "k--", lw=2,
            label="1/sqrt(N) scaling")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set(xlabel="Shots per Monte Carlo", ylabel="Wilson 95% CI half-width",
           title="Statistical uncertainty vs N")
    ax.legend(fontsize=11)
    fig.suptitle("Monte Carlo statistics: how many shots are enough?",
                 fontsize=17, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out / "mc_stats.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_mc_stats/mc_stats.png")


if __name__ == "__main__":
    main()
