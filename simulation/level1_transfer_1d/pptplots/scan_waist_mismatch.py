"""Scan 8: AOD-SLM waist mismatch.

Sweeps the AOD waist while the SLM waist stays 1.17 um. The Level 1 model
assumes matched waists (AOD waist is an `assumed` parameter), so this scan
is a sensitivity study of that assumption: how much does mode mismatch
between the transport and receiving beams matter?
"""
from __future__ import annotations

import scan_lib as sl

WAISTS_UM = [0.85, 0.95, 1.05, 1.17, 1.30, 1.45, 1.65]
SHOTS = 200
SCAN = "scan_waist_mismatch"
X_FIELD = "aod_waist_um"


def job_factory(waist, label):
    return {
        "overrides": {"aod": {"waist_um": float(waist)}},
        "mode": "thermal",
        "label": label if label is not None else f"w={waist:g}um",
    }


def main():
    jobs = sl.chunked_jobs(WAISTS_UM, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(w, f"w={w:g}um"), shots=0,
                       seed=sl.BASE_SEED) for w in WAISTS_UM]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for w in WAISTS_UM:
        lbl = f"w={w:g}um"
        row = sl.aggregate(recs[lbl], w, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    notes = [
        f"shots/point : {SHOTS}",
        "SLM waist   : 1.17 um fixed",
        "baseline    : matched (1.17/1.17)",
        "",
        "Sensitivity study of the assumed",
        "AOD waist: curvature scales as",
        "U0/w^2, so smaller AOD waist =",
        "stiffer transport well.",
    ]
    sl.detail_figure(out, WAISTS_UM, summary, recs, ideals,
                     "AOD waist (um)",
                     "Scan: AOD-SLM waist mismatch  (SLM waist 1.17 um, 5 uK)",
                     vlines=[(1.17, sl.C_NEUTRAL, "matched 1.17 um")],
                     filename="scan_waist_mismatch_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
