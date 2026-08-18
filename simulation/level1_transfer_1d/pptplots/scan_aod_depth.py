"""Scan 7: AOD (transport) well depth.

Sweeps the initial AOD depth. Physics: a deeper AOD well (a) confines the
5 uK thermal ensemble more tightly (smaller position spread, same velocity
spread) and (b) stiffens the trap during the move, changing how the moving
well drags the atom. Also changes the excitation reference (well bottom).
"""
from __future__ import annotations

import scan_lib as sl

DEPTHS_UK = [140, 200, 280, 380, 500, 650, 850]
SHOTS = 200
SCAN = "scan_aod_depth"
X_FIELD = "aod_depth_uK"


def job_factory(depth, label):
    return {
        "overrides": {"aod": {"depth_uK": float(depth)}},
        "mode": "thermal",
        "label": label if label is not None else f"aod={depth:g}uK",
    }


def main():
    jobs = sl.chunked_jobs(DEPTHS_UK, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(d, f"aod={d:g}uK"), shots=0,
                       seed=sl.BASE_SEED) for d in DEPTHS_UK]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for d in DEPTHS_UK:
        lbl = f"aod={d:g}uK"
        row = sl.aggregate(recs[lbl], d, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS}",
        "baseline    : 280 uK",
        "SLM depth   : 140 uK fixed",
        "",
        "Deeper AOD = tighter initial",
        "confinement + stiffer transport;",
        "also lowers the excitation",
        "reference (well bottom).",
    ]
    sl.detail_figure(out, DEPTHS_UK, summary, recs, ideals,
                     "AOD depth (uK)",
                     "Scan: AOD transport-well depth  (5 uK, 600 us)",
                     vlines=[(280, sl.C_NEUTRAL, "baseline 280 uK")],
                     filename="scan_aod_depth_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
