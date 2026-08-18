"""Scan 11: fixed-excitation-energy boundary (classical separatrix).

Instead of thermal sampling, every shot at a scan point starts with the
SAME deterministic excitation energy above the AOD well bottom (position
uniform between the turning points, |v| from energy conservation, random
sign - uniform-in-x on the energy shell, NOT microcanonical). Maps the
capture fraction vs initial excitation and exposes the energy threshold at
which the transfer can no longer hold the atom.

Complements the temperature scan: temperature is a distribution; this scan
is the underlying per-energy response curve.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

EXCITATIONS_UK = [2, 5, 10, 20, 35, 50, 70, 90, 110, 130, 150, 180, 220]
SHOTS = 200
SCAN = "scan_energy_boundary"
X_FIELD = "initial_excitation_uK"


def job_factory(exc, label):
    return {
        "overrides": {},
        "mode": "fixed_excitation",
        "fixed_excitation_uK": float(exc),
        "label": label if label is not None else f"E={exc:g}uK",
    }


def main():
    jobs = sl.chunked_jobs(EXCITATIONS_UK, SHOTS, job_factory)
    out = sl.HERE / SCAN
    recs, _, meta = sl.run_scan(SCAN, jobs, None, X_FIELD, out)

    summary = []
    for e in EXCITATIONS_UK:
        lbl = f"E={e:g}uK"
        row = sl.aggregate(recs[lbl], e, X_FIELD)
        row["__label__"] = lbl
        summary.append(row)
    sl.write_summary(out, summary)

    half = min(summary, key=lambda s: abs(s["capture_fraction"] - 0.5))
    notes = [
        f"shots/point : {SHOTS}",
        "sampling    : fixed energy shell,",
        "  x uniform between turning pts,",
        "  |v| from energy conservation",
        "  (uniform-in-x, not microcanonical)",
        f"~50% capture: E_exc ~ {half[X_FIELD]:g} uK",
        f"SLM depth   : {meta['slm_depth_uK']:g} uK",
        "",
        "Deterministic-energy version of the",
        "temperature scan: shows the energy",
        "threshold where the transfer fails.",
    ]
    sl.detail_figure(out, EXCITATIONS_UK, summary, recs, None,
                     "Initial excitation energy (uK)",
                     "Scan: fixed-excitation capture boundary  (600 us)",
                     vlines=[(140, sl.C_SLM, "SLM depth 140 uK")],
                     filename="scan_energy_boundary_detail.png",
                     note_lines=notes)

    # extra: per-shot final excitation vs initial energy (ridge view)
    fig, ax = plt.subplots(figsize=(10, 6))
    for e in EXCITATIONS_UK:
        vals = [r["final_excitation_energy_uK"] for r in recs[f"E={e:g}uK"]]
        jit = np.random.default_rng(1).normal(0, 0.012 * max(EXCITATIONS_UK), len(vals))
        capd = np.array([r["captured_in_slm"] for r in recs[f"E={e:g}uK"]])
        vals = np.array(vals)
        ax.scatter(np.full(len(vals), e) + jit, vals,
                   c=[sl.C_ATOM if c else sl.C_ESC for c in capd],
                   s=7, alpha=0.45, rasterized=True)
    ax.axhline(140, color=sl.C_SLM, ls="--", lw=2, label="SLM depth 140 uK")
    ax.plot([0, max(EXCITATIONS_UK)], [0, max(EXCITATIONS_UK)], "k:", lw=1.5,
            label="no excitation change")
    ax.set(xlabel="Initial AOD excitation (uK)",
           ylabel="Final SLM excitation (uK)",
           title="Per-shot final excitation vs initial excitation\n"
                 "(green = captured, purple = escaped)")
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "energy_boundary_scatter.png", dpi=170,
                bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_energy_boundary/energy_boundary_scatter.png")


if __name__ == "__main__":
    main()
