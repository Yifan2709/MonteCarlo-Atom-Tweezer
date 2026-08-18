"""Scan 13: integrator timestep convergence.

Sweeps dt around the baseline 0.05 us. Two views:
  - deterministic: ideal-atom final SLM excitation vs the finest-dt reference
  - statistical: capture fraction (200 shots, identical seeds per dt)
The Level 1 validator requires omega_max*dt < 0.05; replace() bypasses the
loader check intentionally to probe where results actually degrade.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

DTS_US = [0.0125, 0.025, 0.05, 0.10, 0.15]
SHOTS = 200
SCAN = "scan_dt"
X_FIELD = "dt_us"


def job_factory(dt, label):
    return {
        "overrides": {"trajectory": {"dt_us": float(dt)}},
        "mode": "thermal",
        "label": label if label is not None else f"dt={dt:g}us",
    }


def main():
    jobs = sl.chunked_jobs(DTS_US, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(d, f"dt={d:g}us"), shots=0,
                       seed=sl.BASE_SEED) for d in DTS_US]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    ref_exc = ideals["dt=0.0125us"]["final_excitation_energy_uK"]
    summary = []
    for d in DTS_US:
        lbl = f"dt={d:g}us"
        row = sl.aggregate(recs[lbl], d, X_FIELD)
        row["__label__"] = lbl
        idl = ideals[lbl]
        row["ideal_excitation_change_uK"] = idl["excitation_energy_change_uK"]
        row["ideal_final_exc_rel_dev_vs_ref"] = (
            abs(idl["final_excitation_energy_uK"] - ref_exc) / abs(ref_exc)
            if ref_exc else None)
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS} (identical seeds)",
        "baseline    : dt = 0.05 us",
        "validator   : omega_max*dt < 0.05",
        "reference   : dt = 0.0125 us",
        "",
        "Same seeds at every dt -> capture",
        "differences reflect the integrator,",
        "not sampling noise.",
    ]
    sl.detail_figure(out, DTS_US, summary, recs, ideals,
                     "Timestep dt (us)",
                     "Scan: integrator timestep convergence  (5 uK, 600 us)",
                     x_log=True,
                     vlines=[(0.05, sl.C_NEUTRAL, "baseline 0.05 us")],
                     filename="scan_dt_detail.png",
                     note_lines=notes)

    # extra: ideal-atom relative deviation vs dt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    dev = [max(s["ideal_final_exc_rel_dev_vs_ref"] or 1e-16, 1e-16)
           for s in summary]
    ax.plot(DTS_US, dev, "o-", color=sl.C_SLM, markersize=9,
            markeredgecolor="black", markeredgewidth=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set(xlabel="Timestep dt (us)", ylabel="|E(dt)-E(ref)| / |E(ref)|",
           title="Ideal-atom final excitation vs finest-dt reference")
    fig.tight_layout()
    fig.savefig(out / "dt_convergence_ideal.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_dt/dt_convergence_ideal.png")


if __name__ == "__main__":
    main()
