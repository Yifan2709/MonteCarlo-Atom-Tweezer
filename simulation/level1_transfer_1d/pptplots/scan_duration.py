"""Scan 1: transfer duration.

Sweeps the total AOD->SLM transfer duration at fixed move/ramp fractions
(0.52/0.48) and 5 uK temperature. Physics: shorter transfers are less
adiabatic (more motional excitation, eventual escape); longer transfers are
harmless in this noiseless classical model but costly in a real experiment
(spontaneous scattering, vacuum lifetime). The scan locates the onset of
non-adiabatic capture loss.

Run from pptplots/:
    PYTHONPATH=../src:. python scan_duration.py
"""
from __future__ import annotations

import numpy as np

import scan_lib as sl

DURATIONS_US = [100, 200, 300, 400, 500, 600, 800, 1000, 1400, 2000]
SHOTS = 200
SCAN = "scan_duration"
X_FIELD = "transfer_duration_us"


def job_factory(duration, label):
    return {
        "overrides": {"transfer": {"duration_us": float(duration)}},
        "mode": "thermal",
        "label": label if label is not None else f"T={duration:g}us",
    }


def main():
    jobs = sl.chunked_jobs(DURATIONS_US, SHOTS, job_factory)
    ideal_jobs = [dict(job_factory(d, f"T={d:g}us"), shots=0,
                       seed=sl.BASE_SEED) for d in DURATIONS_US]
    out = sl.HERE / SCAN
    recs, ideals, meta = sl.run_scan(SCAN, jobs, ideal_jobs, X_FIELD, out)

    summary = []
    for d in DURATIONS_US:
        lbl = f"T={d:g}us"
        row = sl.aggregate(recs[lbl], d, X_FIELD)
        row["__label__"] = lbl
        idl = ideals[lbl]
        row["ideal_captured"] = idl["captured_in_slm"]
        row["ideal_excitation_change_uK"] = idl["excitation_energy_change_uK"]
        summary.append(row)
    sl.write_summary(out, summary)

    worst = min(summary, key=lambda s: s["capture_fraction"])
    notes = [
        f"shots/point      : {SHOTS}",
        f"temperature      : {meta['temperature_uK']:g} uK",
        f"baseline duration: 600 us (capture = "
        f"{next(s['capture_fraction'] for s in summary if s[X_FIELD]==600):.3f})",
        f"worst point      : {worst[X_FIELD]:g} us -> {worst['capture_fraction']:.3f}",
        "",
        "Shorter transfer = less adiabatic:",
        "excitation grows and captures",
        "start to fail below ~300-400 us.",
        "Long transfer is benign in this",
        "noiseless classical model.",
        "",
        "Classical result for this model -",
        "not experimental survival/fidelity.",
    ]
    sl.detail_figure(out, DURATIONS_US, summary, recs, ideals,
                     "Transfer duration (us)",
                     "Scan: transfer duration  (5 uK, AOD 280 uK -> SLM 140 uK)",
                     x_log=True,
                     vlines=[(600, sl.C_NEUTRAL, "baseline 600 us")],
                     filename="scan_duration_detail.png",
                     note_lines=notes)


if __name__ == "__main__":
    main()
