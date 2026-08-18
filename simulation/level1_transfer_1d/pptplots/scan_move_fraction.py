"""Scan 2: move/ramp time allocation (move_fraction).

Sweeps the fraction of the total 600 us spent on the AOD centre move; the
remainder is the depth ramp. Physics: the move phase and the ramp phase
contribute differently to motional excitation - this answers "how should the
waveform budget its time?".
"""
from __future__ import annotations

import scan_lib as sl

FRACTIONS = [0.15, 0.25, 0.35, 0.45, 0.52, 0.60, 0.70, 0.80, 0.90]
SHOTS = 200
SCAN = "scan_move_fraction"
X_FIELD = "move_fraction"


def job_factory(frac, label):
    return {
        "overrides": {"transfer": {"move_fraction": float(frac),
                                   "ramp_fraction": float(1.0 - frac)}},
        "mode": "thermal",
        "label": label if label is not None else f"mf={frac:g}",
    }


def main():
    jobs = sl.chunked_jobs(FRACTIONS, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(f, f"mf={f:g}"), shots=0,
                       seed=sl.BASE_SEED) for f in FRACTIONS]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for f in FRACTIONS:
        lbl = f"mf={f:g}"
        row = sl.aggregate(recs[lbl], f, X_FIELD)
        row["__label__"] = lbl
        row["ideal_excitation_change_uK"] = ideals[lbl]["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    best = max(summary, key=lambda s: s["capture_fraction"])
    notes = [
        f"shots/point : {SHOTS}",
        f"duration    : 600 us fixed",
        f"baseline    : move_fraction = 0.52",
        f"best point  : {best[X_FIELD]:g} -> {best['capture_fraction']:.3f}",
        "",
        "move_fraction = share of the 600 us",
        "spent moving the AOD centre;",
        "the rest is the depth ramp.",
        "",
        "Classical result for this model -",
        "not experimental survival/fidelity.",
    ]
    sl.detail_figure(out, FRACTIONS, summary, recs, ideals,
                     "Move fraction (of 600 us)",
                     "Scan: move/ramp time allocation  (5 uK, 600 us total)",
                     vlines=[(0.52, sl.C_NEUTRAL, "baseline 0.52")],
                     filename="scan_move_fraction_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
