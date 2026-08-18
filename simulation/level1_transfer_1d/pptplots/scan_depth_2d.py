"""Scan 10: 2D depth phase diagram (SLM x AOD).

Coarse 2D map of capture fraction over (SLM depth, AOD depth) - the two
laser-power knobs. Produces a heatmap with the baseline operating point and
the depth-matched diagonal marked.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

SLM_DEPTHS = [50, 80, 110, 140, 180, 230]
AOD_DEPTHS = [140, 210, 280, 380, 500, 650]
SHOTS = 120
SCAN = "scan_depth_2d"
X_FIELD = "slm_depth_uK"


def main():
    points = [(s, a) for a in AOD_DEPTHS for s in SLM_DEPTHS]

    def job_factory(pt, label):
        s, a = pt
        return {
            "overrides": {"slm": {"depth_uK": float(s)},
                          "aod": {"depth_uK": float(a)}},
            "mode": "thermal",
            "label": label if label is not None else f"slm={s:g},aod={a:g}",
        }

    jobs = sl.chunked_jobs(points, SHOTS, job_factory)
    out = sl.HERE / SCAN
    recs, _, meta = sl.run_scan(SCAN, jobs, None, X_FIELD, out)

    rate = np.zeros((len(AOD_DEPTHS), len(SLM_DEPTHS)))
    dE = np.zeros_like(rate)
    lo = np.zeros_like(rate); hi = np.zeros_like(rate)
    summary = []
    for i, a in enumerate(AOD_DEPTHS):
        for j, s in enumerate(SLM_DEPTHS):
            lbl = f"slm={s:g},aod={a:g}"
            row = sl.aggregate(recs[lbl], s, X_FIELD)
            row["__label__"] = lbl
            row["aod_depth_uK"] = a
            summary.append(row)
            rate[i, j] = row["capture_fraction"]
            dE[i, j] = row["mean_excitation_change_uK"]
            lo[i, j] = row["wilson_lo"]; hi[i, j] = row["wilson_hi"]
    sl.write_summary(out, summary)

    # padded extent so edge-cell annotations are not clipped
    dx = (max(SLM_DEPTHS) - min(SLM_DEPTHS)) / (len(SLM_DEPTHS) - 1) / 2
    dy = (max(AOD_DEPTHS) - min(AOD_DEPTHS)) / (len(AOD_DEPTHS) - 1) / 2
    EXT = [min(SLM_DEPTHS) - dx, max(SLM_DEPTHS) + dx,
           min(AOD_DEPTHS) - dy, max(AOD_DEPTHS) + dy]

    # heatmap: capture fraction
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ax = axes[0]
    im = ax.imshow(rate, origin="lower", aspect="auto", cmap="RdYlGn",
                   vmin=0.0, vmax=1.0, extent=EXT)
    ax.plot([50, 230], [50, 230], "w--", lw=2, label="SLM = AOD")
    ax.plot(140, 280, marker="*", ms=22, color="black", mec="white",
            mew=1.5, ls="none", label="baseline (140, 280)")
    for i, a in enumerate(AOD_DEPTHS):
        for j, s in enumerate(SLM_DEPTHS):
            ax.text(s, a, f"{rate[i,j]:.2f}", ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="black" if 0.3 < rate[i, j] < 0.9 else "white")
    ax.set(xlabel="SLM depth (uK)", ylabel="AOD depth (uK)",
           title="Capture fraction (5 uK, 600 us)")
    ax.legend(loc="lower right", fontsize=10)
    fig.colorbar(im, ax=ax, label="capture fraction")

    ax = axes[1]
    ciw = (hi - lo) / 2
    im2 = ax.imshow(ciw, origin="lower", aspect="auto", cmap="magma",
                    extent=EXT)
    ax.set(xlabel="SLM depth (uK)", ylabel="AOD depth (uK)",
           title=f"Wilson 95% CI half-width ({SHOTS} shots/point)")
    fig.colorbar(im2, ax=ax, label="CI half-width")
    fig.suptitle("Depth phase diagram: capture fraction over the two power knobs",
                 fontsize=17, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out / "depth_2d_heatmap.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_depth_2d/depth_2d_heatmap.png")

    # cuts: capture vs SLM depth at fixed AOD, and vs AOD at fixed SLM=140
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    for i, a in enumerate(AOD_DEPTHS):
        ax.plot(SLM_DEPTHS, rate[i], "o-", label=f"AOD {a} uK")
    ax.set(xlabel="SLM depth (uK)", ylabel="Capture fraction",
           ylim=(-0.03, 1.05), title="Cuts at fixed AOD depth")
    ax.legend(fontsize=9, ncol=2)
    ax = axes[1]
    for j, s in enumerate(SLM_DEPTHS):
        ax.plot(AOD_DEPTHS, rate[:, j], "s-", label=f"SLM {s} uK")
    ax.set(xlabel="AOD depth (uK)", ylabel="Capture fraction",
           ylim=(-0.03, 1.05), title="Cuts at fixed SLM depth")
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "depth_2d_cuts.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_depth_2d/depth_2d_cuts.png")

    # heating heatmap: mean excitation change (the informative observable -
    # capture is 1.0 everywhere at 5 uK)
    fig, ax = plt.subplots(figsize=(8.5, 6))
    lim = float(np.max(np.abs(dE)))
    im = ax.imshow(dE, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-lim, vmax=lim, extent=EXT)
    ax.plot([50, 230], [50, 230], "k--", lw=2, label="SLM = AOD")
    ax.plot(140, 280, marker="*", ms=22, color="black", mec="white",
            mew=1.5, ls="none", label="baseline (140, 280)")
    for i, a in enumerate(AOD_DEPTHS):
        for j, s in enumerate(SLM_DEPTHS):
            ax.text(s, a, f"{dE[i,j]:+.1f}", ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="white" if abs(dE[i, j]) > 0.55 * lim else "black")
    ax.set(xlabel="SLM depth (uK)", ylabel="AOD depth (uK)")
    ax.legend(loc="upper left", fontsize=10)
    fig.colorbar(im, ax=ax, label="mean dE_exc (uK)")
    fig.suptitle("Depth phase diagram: mean excitation-energy change (uK)\n"
                 "red = net heating, blue = net cooling",
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out / "depth_2d_heating.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  fig -> scan_depth_2d/depth_2d_heating.png")


if __name__ == "__main__":
    main()
