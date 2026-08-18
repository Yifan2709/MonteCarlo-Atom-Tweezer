"""Scan 3: AOD centre-move waveform shape.

Compares five move profiles at fixed duration (600 us) and fractions
(0.52/0.48): linear, smoothstep (baseline 3s^2-2s^3), smootherstep
(6s^5-15s^4+10s^3), sine, ease-in-quadratic. Physics: profiles with lower
peak velocity/acceleration should excite less; the scan quantifies how much
the profile shape matters for capture and heating.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

SHAPES = ["linear", "smoothstep", "smootherstep", "sine", "ease_in_quad"]
SHOTS = 200
SCAN = "scan_waveform_shape"
X_FIELD = "move_shape"


def job_factory(shape, label):
    return {
        "overrides": {},
        "mode": "thermal",
        "move_shape": shape,
        "label": label if label is not None else shape,
    }


def main():
    jobs = sl.chunked_jobs(SHAPES, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(s, s), shots=0, seed=sl.BASE_SEED)
                  for s in SHAPES]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for i, s in enumerate(SHAPES):
        row = sl.aggregate(recs[s], i, X_FIELD)   # x = index for plotting
        row["__label__"] = s
        row["shape"] = s
        row["ideal_excitation_change_uK"] = ideals[s]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    x = np.arange(len(SHAPES))
    notes = [
        f"shots/shape : {SHOTS}",
        "baseline    : smoothstep (3s^2-2s^3)",
        "",
        "All profiles move the AOD centre",
        "2.4 um in 312 us; they differ in",
        "peak velocity / acceleration.",
    ]
    sl.detail_figure(out, x, summary, recs, ideals,
                     "Move profile (index)",
                     "Scan: AOD move waveform shape  (5 uK, 600 us)",
                     filename="scan_waveform_shape_detail.png",
                     note_lines=notes)

    # extra figure: the profiles themselves + capture comparison with labels
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    sgrid = np.linspace(0, 1, 400)
    for s in SHAPES:
        ax.plot(sgrid, [sl.move_h(s, v) for v in sgrid], lw=2.5, label=s)
    ax.set(xlabel="Normalized time s = t / t_move",
           ylabel="Moved fraction h(s)", title="Move profiles h(s)")
    ax.legend(fontsize=11)
    ax = axes[1]
    rate = [r["capture_fraction"] for r in summary]
    elo = [max(r["capture_fraction"] - r["wilson_lo"], 0) for r in summary]
    ehi = [max(r["wilson_hi"] - r["capture_fraction"], 0) for r in summary]
    ax.errorbar(x, rate, yerr=[elo, ehi], fmt="o", color=sl.C_ATOM,
                ecolor=sl.C_NEUTRAL, elinewidth=2, capsize=6, markersize=10,
                markeredgecolor="black", markeredgewidth=1)
    ax.set_xticks(x); ax.set_xticklabels(SHAPES, rotation=20, ha="right")
    ax.set(ylabel="Capture fraction", ylim=(0.9, 1.02),
           title="Capture fraction by profile")
    fig.suptitle("AOD move waveform: profiles and capture", fontsize=17,
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out / "waveform_shapes_comparison.png", dpi=170,
                bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_waveform_shape/waveform_shapes_comparison.png")


if __name__ == "__main__":
    main()
