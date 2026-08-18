"""Summary dashboard across all parameter scans.

Reads every scan_* directory's summary.csv and composes:
  1. a 3x4 headline dashboard (capture-fraction curve per scan)
  2. a sensitivity ranking bar chart (capture-fraction range across each scan)
  3. SCAN_RESULTS.md - tables + key findings per scan

Run after all scan scripts have completed:
    PYTHONPATH=../src:. python make_scan_dashboard.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import scan_lib as sl

HERE = Path(__file__).resolve().parent

# (scan_dir, x_field, x_label, title, x_log, baseline_x)
PANELS = [
    ("scan_duration",       "transfer_duration_us", "Duration (us)",      "1. Transfer duration",        True,  600),
    ("scan_move_fraction",  "move_fraction",        "Move fraction",      "2. Move/ramp allocation",     False, 0.52),
    ("scan_waveform_shape", None,                   None,                 "3. Move waveform shape",      False, None),
    ("scan_ramp_exponent",  "ramp_exponent",        "Ramp exponent p",    "4. Depth-ramp exponent",      False, 2.0),
    ("scan_separation",     "initial_separation_um","Separation (um)",    "5. Initial separation",       False, 2.4),
    ("scan_slm_depth",      "slm_depth_uK",         "SLM depth (uK)",     "6. SLM depth",                False, 140),
    ("scan_aod_depth",      "aod_depth_uK",         "AOD depth (uK)",     "7. AOD depth",                False, 280),
    ("scan_waist_mismatch", "aod_waist_um",         "AOD waist (um)",     "8. Waist mismatch",           False, 1.17),
    ("scan_alignment",      "slm_center_offset_um", "SLM offset (um)",    "9. Alignment offset",         False, 0.0),
    ("scan_energy_boundary","initial_excitation_uK","Init. excitation (uK)","10. Energy boundary",       False, None),
    ("scan_dt",             "dt_us",                "dt (us)",            "11. Timestep",                True,  0.05),
    ("scan_temperature",    "temperature_uK",       "Temperature (uK)",   "12. Temperature (prior work)",True,  5),
]


def read_summary(scan_dir):
    path = HERE / scan_dir / "summary.csv"
    if not path.is_file():
        alt = HERE / scan_dir / "summary_by_temperature.csv"
        path = alt if alt.is_file() else path
    if not path.is_file():
        return None
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    available = []
    for scan_dir, *_ in PANELS:
        rows = read_summary(scan_dir)
        if rows:
            available.append((scan_dir, rows))
    if not available:
        raise SystemExit("no completed scans found")

    # ---- dashboard ----
    n = len(PANELS)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, 5.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (scan_dir, x_field, x_label, title, x_log, base_x) in zip(axes, PANELS):
        rows = read_summary(scan_dir)
        if not rows:
            ax.text(0.5, 0.5, "not run", ha="center", va="center",
                    transform=ax.transAxes, fontsize=14, color=sl.C_NEUTRAL)
            ax.set_title(title)
            continue
        if scan_dir == "scan_waveform_shape":
            labels = [r.get("shape") or r.get("__label__") for r in rows]
            x = np.arange(len(rows))
            rate = [float(r["capture_fraction"]) for r in rows]
            lo = [float(r["wilson_lo"]) for r in rows]
            hi = [float(r["wilson_hi"]) for r in rows]
            ax.errorbar(x, rate, yerr=[np.clip(np.array(rate)-np.array(lo),0,None),
                                       np.clip(np.array(hi)-np.array(rate),0,None)],
                        fmt="o", color=sl.C_ATOM, ecolor=sl.C_NEUTRAL,
                        capsize=4, markersize=7, markeredgecolor="black")
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=9)
        else:
            x = np.array([float(r[x_field]) for r in rows])
            rate = np.array([float(r["capture_fraction"]) for r in rows])
            lo = np.array([float(r["wilson_lo"]) for r in rows])
            hi = np.array([float(r["wilson_hi"]) for r in rows])
            order = np.argsort(x)
            x, rate, lo, hi = x[order], rate[order], lo[order], hi[order]
            ax.errorbar(x, rate, yerr=[np.clip(rate-lo,0,None),
                                       np.clip(hi-rate,0,None)],
                        fmt="-o", color=sl.C_ATOM, ecolor=sl.C_NEUTRAL,
                        capsize=3, markersize=5, markeredgecolor="black",
                        markeredgewidth=0.6, elinewidth=1.2)
            if base_x is not None:
                ax.axvline(base_x, color=sl.C_NEUTRAL, ls="--", lw=1.8)
            if x_log:
                ax.set_xscale("log")
            ax.set_xlabel(x_label)
        ax.set(ylim=(-0.05, 1.1), title=title, ylabel="capture fraction")
    for ax in axes[len(PANELS):]:
        ax.axis("off")
    fig.suptitle("Level 1 parameter scans - capture fraction overview "
                 "(classical 1D model, not experimental fidelity)",
                 fontsize=20, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(HERE / "scan_dashboard.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("fig -> scan_dashboard.png")

    # ---- heating dashboard: mean excitation change (the informative
    # observable - capture fraction is flat ~1.0 at 5 uK in most scans) ----
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, 5.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (scan_dir, x_field, x_label, title, x_log, base_x) in zip(axes, PANELS):
        rows = read_summary(scan_dir)
        if not rows or scan_dir == "scan_waveform_shape":
            if not rows:
                ax.text(0.5, 0.5, "not run", ha="center", va="center",
                        transform=ax.transAxes, fontsize=14, color=sl.C_NEUTRAL)
            if scan_dir == "scan_waveform_shape" and rows:
                labels = [r.get("shape") or r.get("__label__") for r in rows]
                x = np.arange(len(rows))
                dEm = [float(r["mean_excitation_change_uK"]) for r in rows]
                dEi = [float(r.get("ideal_excitation_change_uK") or "nan")
                       for r in rows]
                ax.plot(x, dEm, "o-", color=sl.C_HL, label="MC mean")
                ax.plot(x, dEi, "d--", color=sl.C_SLM, label="ideal atom")
                ax.axhline(0, color="black", ls="--", lw=1.5)
                ax.set_xticks(x)
                ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=9)
                ax.legend(fontsize=9)
                ax.set_title(title)
            continue
        key = "temperature_uK" if scan_dir == "scan_temperature" else x_field
        x = np.array([float(r[key]) for r in rows])
        dEm = np.array([float(r["mean_excitation_change_uK"]) for r in rows])
        order = np.argsort(x)
        x, dEm = x[order], dEm[order]
        ax.plot(x, dEm, "o-", color=sl.C_HL, markersize=5,
                markeredgecolor="black", markeredgewidth=0.6,
                label="MC mean dE")
        if any(r.get("ideal_excitation_change_uK") for r in rows):
            dEi = np.array([float(r["ideal_excitation_change_uK"])
                            for r in rows])[order]
            ax.plot(x, dEi, "d--", color=sl.C_SLM, markersize=5,
                    label="ideal atom")
        ax.axhline(0, color="black", ls="--", lw=1.5)
        if base_x is not None:
            ax.axvline(base_x, color=sl.C_NEUTRAL, ls="--", lw=1.5)
        if x_log:
            ax.set_xscale("log")
        ax.set(xlabel=x_label, ylabel="mean dE_exc (uK)", title=title)
        ax.legend(fontsize=9)
    for ax in axes[len(PANELS):]:
        ax.axis("off")
    fig.suptitle("Level 1 parameter scans - motional excitation change "
                 "(positive = net heating; the sensitive observable at 5 uK)",
                 fontsize=20, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(HERE / "scan_heating_dashboard.png", dpi=160,
                bbox_inches="tight")
    plt.close(fig)
    print("fig -> scan_heating_dashboard.png")

    # ---- sensitivity ranking by heating range (capture is ~flat at 5 uK) ----
    ranking = []
    for scan_dir, x_field, _, title, _, _ in PANELS:
        rows = read_summary(scan_dir)
        if not rows:
            continue
        rates = [float(r["capture_fraction"]) for r in rows]
        dEs = [float(r["mean_excitation_change_uK"]) for r in rows
               if r.get("mean_excitation_change_uK") not in (None, "")]
        ranking.append({
            "scan": title,
            "scan_dir": scan_dir,
            "rate_min": min(rates), "rate_max": max(rates),
            "rate_range": max(rates) - min(rates),
            "dE_range": (max(dEs) - min(dEs)) if dEs else 0.0,
            "max_mean_heating_uK": max(dEs) if dEs else None,
        })
    ranking.sort(key=lambda r: -r["dE_range"])

    fig, ax = plt.subplots(figsize=(10, 6.5))
    names = [r["scan"] for r in ranking]
    vals = [r["dE_range"] for r in ranking]
    colors = [sl.C_SLM if v > 2.0 else sl.C_HL if v > 0.5 else sl.C_ATOM
              for v in vals]
    ax.barh(range(len(names)), vals, color=colors, alpha=0.85, edgecolor="black")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=11)
    ax.invert_yaxis()
    ax.set(xlabel="Mean excitation-change range across scan (uK)",
           title="Which parameters move the atom's energy most?\n"
                 "(red > 2 uK, orange > 0.5 uK; capture fraction is ~flat "
                 "at 5 uK except energy boundary)")
    for i, v in enumerate(vals):
        ax.text(v + 0.005, i, f"{v:.2f}", va="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "scan_sensitivity_ranking.png", dpi=170,
                bbox_inches="tight")
    plt.close(fig)
    print("fig -> scan_sensitivity_ranking.png")

    # ---- markdown report ----
    lines = ["# Level 1 parameter scan results", "",
             "Classical 1D AOD->SLM transfer model (Cs-133, 5 uK baseline "
             "temperature, 600 us baseline duration). All capture fractions "
             "are classical mechanics results for this configured model - "
             "not experimental survival, transfer fidelity, or quantum "
             "fidelity.", ""]
    for scan_dir, x_field, x_label, title, _, _ in PANELS:
        rows = read_summary(scan_dir)
        if not rows:
            continue
        meta_path = HERE / scan_dir / "scan_meta.json"
        meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        lines.append(f"## {title}")
        lines.append("")
        if x_field:
            hdr = f"| {x_field} | shots | captured | capture fraction | Wilson 95% CI | mean dE_exc (uK) |"
            lines.append(hdr)
            lines.append("|---|---|---|---|---|---|")
            key = "temperature_uK" if scan_dir == "scan_temperature" else x_field
            for r in sorted(rows, key=lambda r: float(r[key])):
                lines.append(
                    f"| {float(r[key]):g} | {r['shots']} | {r['captured']} "
                    f"| {float(r['capture_fraction']):.4f} "
                    f"| [{float(r['wilson_lo']):.3f}, {float(r['wilson_hi']):.3f}] "
                    f"| {float(r['mean_excitation_change_uK']):.2f} |")
        else:
            lines.append("| shape | shots | captured | capture fraction | Wilson 95% CI | mean dE_exc (uK) |")
            lines.append("|---|---|---|---|---|---|")
            for r in rows:
                nm = r.get("shape") or r.get("__label__")
                lines.append(
                    f"| {nm} | {r['shots']} | {r['captured']} "
                    f"| {float(r['capture_fraction']):.4f} "
                    f"| [{float(r['wilson_lo']):.3f}, {float(r['wilson_hi']):.3f}] "
                    f"| {float(r['mean_excitation_change_uK']):.2f} |")
        if meta:
            lines.append("")
            lines.append(f"_elapsed {meta.get('elapsed_s', 0):.0f}s, "
                         f"{meta.get('total_shots', '?')} shots_")
        lines.append("")
    (HERE / "SCAN_RESULTS.md").write_text("\n".join(lines))
    print("report -> SCAN_RESULTS.md")


if __name__ == "__main__":
    main()
