"""Scan 9: SLM-AOD alignment offset.

Sweeps the static SLM centre position; the AOD still moves to x = 0.
Physics: models experimental pointing / calibration error - the AOD hands
the atom off at x=0 while the receiving SLM well sits at x = offset. Maps
the capture tolerance to sub-waist misalignment.
"""
from __future__ import annotations

import scan_lib as sl

OFFSETS_UM = [-0.80, -0.50, -0.30, -0.15, 0.0, 0.15, 0.30, 0.50, 0.80]
SHOTS = 200
SCAN = "scan_alignment"
X_FIELD = "slm_center_offset_um"


def job_factory(off, label):
    return {
        "overrides": {"slm": {"center_um": float(off)}},
        "mode": "thermal",
        "label": label if label is not None else f"off={off:g}um",
    }


def main():
    jobs = sl.chunked_jobs(OFFSETS_UM, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(o, f"off={o:g}um"), shots=0,
                       seed=sl.BASE_SEED) for o in OFFSETS_UM]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for o in OFFSETS_UM:
        lbl = f"off={o:g}um"
        row = sl.aggregate(recs[lbl], o, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS}",
        "AOD endpoint: x = 0",
        "baseline    : perfectly aligned",
        f"waist       : 1.17 um",
        "",
        "Models pointing/calibration error:",
        "the AOD hands the atom off at x=0",
        "while the SLM well sits at offset.",
        "Tolerance to sub-waist offsets is",
        "an experimental error budget item.",
    ]
    sl.detail_figure(out, OFFSETS_UM, summary, recs, ideals,
                     "SLM centre offset (um)",
                     "Scan: SLM-AOD alignment offset  (5 uK, 600 us)",
                     vlines=[(0.0, sl.C_NEUTRAL, "aligned"),
                             (1.17, sl.C_HL, "1 waist"), (-1.17, sl.C_HL, "")],
                     filename="scan_alignment_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
