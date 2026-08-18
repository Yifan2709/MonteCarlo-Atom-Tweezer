"""Scan 12: initial-orbit phase sensitivity.

Fixes the initial excitation energy (50 uK) and samples the initial phase
coordinate of the orbit in the static AOD well (normalized position
(x0-xc)/x_turn in [-1, 1], |v| from energy conservation, both velocity
branches kept and analyzed separately). Physics: if the capture outcome
depends on where in its orbit the atom sits when the AOD starts moving, the
experiment could in principle trigger the transfer synchronized to the atom
phase; if not, phase-averaging is justified.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

EXC_UK = 50.0
SHOTS = 1200
SCAN = "scan_phase"


def main():
    def job_factory(_, label):
        return {
            "overrides": {},
            "mode": "fixed_excitation",
            "fixed_excitation_uK": EXC_UK,
            "label": label if label is not None else f"E={EXC_UK:g}uK",
        }
    jobs = sl.chunked_jobs([EXC_UK], SHOTS, job_factory)
    out = sl.HERE / SCAN
    recs, _, meta = sl.run_scan(SCAN, jobs, None, "excitation_uK", out)
    recs = recs[f"E={EXC_UK:g}uK"]

    phase = np.array([r["initial_phase_coord"] for r in recs])
    vsign = np.array([np.sign(r["v0_m_s"]) for r in recs])
    capd = np.array([r["captured_in_slm"] for r in recs])
    dE = np.array([r["excitation_energy_change_uK"] for r in recs])

    bins = np.linspace(-1, 1, 21)
    centres = 0.5 * (bins[:-1] + bins[1:])

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(f"Initial-phase sensitivity at fixed excitation {EXC_UK:g} uK "
                 f"({SHOTS} shots)", fontsize=17, fontweight="bold")

    for ax, sign, name in [(axes[0, 0], 1, "v0 > 0"), (axes[0, 1], -1, "v0 < 0")]:
        m = vsign == sign
        frac, _ = np.histogram(phase[m & capd], bins=bins)
        tot, _ = np.histogram(phase[m], bins=bins)
        rate = np.divide(frac, tot, out=np.full_like(frac, np.nan, dtype=float),
                         where=tot > 0)
        ax.bar(centres, rate, width=0.045, color=sl.C_ATOM, alpha=0.75,
               edgecolor="black")
        ax.axhline(np.mean(capd[m]), color=sl.C_NEUTRAL, ls="--", lw=2,
                   label=f"mean {np.mean(capd[m]):.3f}")
        ax.set(xlabel="Initial phase coordinate (x0-xc)/x_turn",
               ylabel="Capture fraction", ylim=(0, 1.1),
               title=f"Capture vs initial phase, {name}")
        ax.legend(fontsize=11)

    ax = axes[1, 0]
    ax.scatter(phase[capd], dE[capd], s=9, alpha=0.35, color=sl.C_ATOM,
               label="captured", rasterized=True)
    ax.scatter(phase[~capd], dE[~capd], s=22, alpha=0.85, color=sl.C_ESC,
               marker="x", label="escaped", rasterized=True)
    ax.axhline(0, color="black", ls="--", lw=1.5)
    ax.set(xlabel="Initial phase coordinate (x0-xc)/x_turn",
           ylabel="Excitation-energy change (uK)",
           title="Excitation change vs initial phase")
    ax.legend(fontsize=11)

    ax = axes[1, 1]
    ax.hist(phase[vsign > 0], bins=bins, histtype="step", lw=2.2,
            color=sl.C_AOD, label="v0 > 0")
    ax.hist(phase[vsign < 0], bins=bins, histtype="step", lw=2.2,
            color=sl.C_SLM, label="v0 < 0")
    ax.set(xlabel="Initial phase coordinate (x0-xc)/x_turn", ylabel="Shots",
           title="Phase-coordinate sampling coverage")
    ax.legend(fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out / "phase_sensitivity.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_phase/phase_sensitivity.png")

    # summary text
    n_pos = int(np.sum(vsign > 0)); n_neg = int(np.sum(vsign < 0))
    lines = [
        f"fixed excitation       : {EXC_UK:g} uK",
        f"shots                  : {SHOTS} (v>0: {n_pos}, v<0: {n_neg})",
        f"capture fraction v0>0  : {np.mean(capd[vsign > 0]):.4f}",
        f"capture fraction v0<0  : {np.mean(capd[vsign < 0]):.4f}",
        f"escaped shots          : {int(np.sum(~capd))}",
    ]
    (out / "phase_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
