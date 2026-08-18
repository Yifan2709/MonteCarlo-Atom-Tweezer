"""Scan 4: AOD depth-ramp exponent.

The AOD depth ramps down as U(s) = U0 (1-s)^p during the ramp phase.
Baseline p = 2. Sweeps p in [0.5 .. 6]. Physics: the exponent sets how long
the atom spends in a shallow AOD well while the SLM well already dominates;
very slow initial ramps (large p) vs abrupt late drops (small p) excite
differently.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

EXPONENTS = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
SHOTS = 200
SCAN = "scan_ramp_exponent"
X_FIELD = "ramp_exponent"


def job_factory(p, label):
    return {
        "overrides": {},
        "mode": "thermal",
        "ramp_exponent": float(p),
        "label": label if label is not None else f"p={p:g}",
    }


def main():
    jobs = sl.chunked_jobs(EXPONENTS, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(p, f"p={p:g}"), shots=0,
                       seed=sl.BASE_SEED) for p in EXPONENTS]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for p in EXPONENTS:
        lbl = f"p={p:g}"
        row = sl.aggregate(recs[lbl], p, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS}",
        "baseline    : p = 2 (quadratic)",
        "",
        "U_aod(s) = U0 (1-s)^p during",
        "the 288 us ramp phase.",
        "p<1 drops fast early; p>1 keeps",
        "the AOD deep, then falls steeply.",
    ]
    sl.detail_figure(out, EXPONENTS, summary, recs, ideals,
                     "Ramp exponent p",
                     "Scan: AOD depth-ramp exponent  (5 uK, 600 us)",
                     vlines=[(2.0, sl.C_NEUTRAL, "baseline p=2")],
                     filename="scan_ramp_exponent_detail.png",
                     note_lines=notes)

    # extra: ramp curves
    fig, ax = plt.subplots(figsize=(8, 5.5))
    sgrid = np.linspace(0, 1, 400)
    for p in EXPONENTS:
        ax.plot(sgrid, (1 - sgrid) ** p, lw=2.5, label=f"p={p:g}")
    ax.set(xlabel="Normalized ramp time s", ylabel="U_aod / U0",
           title="Depth-ramp laws compared")
    ax.legend(fontsize=11, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "ramp_laws.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_ramp_exponent/ramp_laws.png")


if __name__ == "__main__":
    main()
