"""Scan 5: initial AOD-SLM separation.

Sweeps the initial AOD centre position (final SLM position stays 0).
Physics: a longer haul means higher peak transport speed and stronger
non-adiabatic excitation; at very small separations the two wells overlap
strongly during the whole sequence. Maps capture vs haul distance.
"""
from __future__ import annotations

import scan_lib as sl

SEPARATIONS_UM = [0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.8, 4.5, 5.5]
SHOTS = 200
SCAN = "scan_separation"
X_FIELD = "initial_separation_um"


def job_factory(sep, label):
    return {
        "overrides": {"aod": {"initial_center_um": float(sep)}},
        "mode": "thermal",
        "label": label if label is not None else f"d={sep:g}um",
    }


def main():
    jobs = sl.chunked_jobs(SEPARATIONS_UM, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(s, f"d={s:g}um"), shots=0,
                       seed=sl.BASE_SEED) for s in SEPARATIONS_UM]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for s in SEPARATIONS_UM:
        lbl = f"d={s:g}um"
        row = sl.aggregate(recs[lbl], s, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS}",
        "baseline    : 2.4 um",
        f"waist       : 1.17 um (both wells)",
        "",
        "Larger separation = faster peak",
        "transport speed in the same 312 us",
        "move window = more excitation.",
        "Small separation = wells overlap",
        "early, gentler effective hand-off.",
    ]
    sl.detail_figure(out, SEPARATIONS_UM, summary, recs, ideals,
                     "Initial AOD-SLM separation (um)",
                     "Scan: initial trap separation  (5 uK, 600 us)",
                     vlines=[(2.4, sl.C_NEUTRAL, "baseline 2.4 um"),
                             (1.17, sl.C_HL, "waist 1.17 um")],
                     filename="scan_separation_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
