"""Scan 6: SLM (receiving) well depth.

Sweeps the static SLM depth. Physics: the SLM well is the final container -
a shallow well can only hold atoms whose excitation stays below its depth,
so the capture fraction should collapse once the SLM depth approaches the
atom's post-transfer excitation. Directly maps to the laser-power budget.
"""
from __future__ import annotations

import scan_lib as sl

DEPTHS_UK = [30, 50, 70, 90, 110, 140, 170, 210, 260, 320]
SHOTS = 200
SCAN = "scan_slm_depth"
X_FIELD = "slm_depth_uK"


def job_factory(depth, label):
    return {
        "overrides": {"slm": {"depth_uK": float(depth)}},
        "mode": "thermal",
        "label": label if label is not None else f"slm={depth:g}uK",
    }


def main():
    jobs = sl.chunked_jobs(DEPTHS_UK, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(d, f"slm={d:g}uK"), shots=0,
                       seed=sl.BASE_SEED) for d in DEPTHS_UK]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for d in DEPTHS_UK:
        lbl = f"slm={d:g}uK"
        row = sl.aggregate(recs[lbl], d, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    half = min(summary, key=lambda s: abs(s["capture_fraction"] - 0.5))
    notes = [
        f"shots/point : {SHOTS}",
        "baseline    : 140 uK",
        f"temperature : 5 uK",
        f"~50% capture: SLM depth ~ {half[X_FIELD]:g} uK",
        "",
        "Shallower receiving well = less",
        "room for post-transfer excitation.",
        "Steep rolloff expected well below",
        "the 5 uK atom temperature times",
        "the excitation growth factor.",
    ]
    sl.detail_figure(out, DEPTHS_UK, summary, recs, ideals,
                     "SLM depth (uK)",
                     "Scan: SLM receiving-well depth  (5 uK, 600 us)",
                     vlines=[(140, sl.C_NEUTRAL, "baseline 140 uK")],
                     filename="scan_slm_depth_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
