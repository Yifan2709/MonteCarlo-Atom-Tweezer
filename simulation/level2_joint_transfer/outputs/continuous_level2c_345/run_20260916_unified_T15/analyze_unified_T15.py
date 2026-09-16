"""Analysis for the unified T0=Tref=15 uK Fig. 6d rerun (run_20260916_unified_T15).

Verifies (a) zero initially unbound atoms, (b) common random numbers shared by
all four waveforms, (c) identical heating powers, (d) paired timestep stability,
then emits the comparison figure, per-transfer CSVs, S(n) table, loss rates,
reference metrics (published fit line and experimental markers kept separate),
and the archived pre-fix contrast. No parameter is fitted anywhere here.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # simulation/level2_joint_transfer
RUN14 = ROOT / "outputs/continuous_level2c_345/run_20260914"
REF = ROOT / "reference/fig6d_vector_v2"
WAVEFORMS = ["manual_200us", "manual_400us", "manual_600us", "ml_400us"]
LABEL = {"manual_200us": "manual 200 us", "manual_400us": "manual 400 us",
         "manual_600us": "manual 600 us", "ml_400us": "paper ML 400 us"}
FORMAL = {w: HERE / f"validation_T15_matched_{w}_N4096_dt0.0125_initSite_Tref15.json" for w in WAVEFORMS}
DTCHECK = {w: HERE / f"convergence_T15_matched_{w}_N4096_dt0.025_initSite_Tref15.json" for w in WAVEFORMS}
ARCHIVE_T = {"manual_200us": 5, "manual_400us": 15, "manual_600us": 25, "ml_400us": 25}


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_npz(path):
    return np.load(str(path.with_suffix(".npz")))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load_reference():
    rows = list(csv.DictReader((REF / "published_fit_reference.csv").open()))
    fit = {}
    for w in WAVEFORMS:
        pts = [(int(r["n_one_way_transfers"]), float(r[f"survival_{w}"]))
               for r in rows if r[f"survival_{w}"]]
        fit[w] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
    rows = list(csv.DictReader((REF / "experimental_marker_centers.csv").open()))
    markers = {}
    for w in WAVEFORMS:
        pts = [(int(r["n"]), float(r["survival"]),
                float(r["survival"]) - float(r["error_lower"]),
                float(r["error_upper"]) - float(r["survival"]))
               for r in rows if r["waveform"] == w]
        markers[w] = {k: np.array([p[i] for p in pts]) for i, k in enumerate(("n", "s", "elo", "ehi"))}
    return fit, markers


def metrics(n_sim, s_sim, n_ref, s_ref, band=0.02):
    common = n_ref[n_ref <= n_sim.max()]
    sim_at = np.interp(common, n_sim, s_sim)
    dev = sim_at - s_ref[n_ref <= n_sim.max()]
    return {"n_points": int(len(common)), "n_min": int(common.min()), "n_max": int(common.max()),
            "rmse": float(np.sqrt(np.mean(dev ** 2))), "max_abs_dev": float(np.max(np.abs(dev))),
            "argmax_n": int(common[int(np.argmax(np.abs(dev)))]),
            "fraction_within_band": float(np.mean(np.abs(dev) <= band))}


def per_transfer_loss(n, s, lo, hi):
    """Mean per-transfer conditional loss rate inside [lo, hi] (geometric mean)."""
    k = hi - lo
    return float(1.0 - (s[hi] / s[lo]) ** (1.0 / k)) if s[lo] > 0 and s[hi] > 0 else float("nan")


def main():
    report = {}

    # --- 1. integrity checks -------------------------------------------------
    formal = {w: load_json(p) for w, p in FORMAL.items()}
    dtc = {w: load_json(p) for w, p in DTCHECK.items()}
    unbound = {w: formal[w]["initial_unbound_actual_site"] for w in WAVEFORMS}
    assert all(v == 0 for v in unbound.values()), f"unbound atoms: {unbound}"
    npz = {w: load_npz(p) for w, p in FORMAL.items()}
    base = WAVEFORMS[0]
    crn = {w: bool(np.array_equal(npz[base]["initial_pos"], npz[w]["initial_pos"])
                   and np.array_equal(npz[base]["initial_vel"], npz[w]["initial_vel"])
                   and np.array_equal(npz[base]["initial_slm_depth_factor"],
                                      npz[w]["initial_slm_depth_factor"]))
           for w in WAVEFORMS}
    assert all(crn.values()), f"CRN broken: {crn}"
    energies_ok = bool(np.max([npz[w]["initial_energy_j"].max() for w in WAVEFORMS]) < 0)
    assert energies_ok, "some atoms are not bound in their actual site"
    powers = {w: formal[w]["heating_power_uK_per_s"]["total"] for w in WAVEFORMS}
    assert len(set(powers.values())) == 1, f"heating powers differ: {powers}"
    # survival recomputed from the per-atom NPZ record must match the metadata
    for w in WAVEFORMS:
        assert np.array_equal(npz[w]["alive"].sum(axis=0),
                              np.asarray(formal[w]["survival_counts"])), w
    report["integrity"] = {"initial_unbound_actual_site": unbound,
                           "common_random_numbers_across_waveforms": crn,
                           "all_initially_bound": energies_ok,
                           "heating_total_power_uK_per_s": powers[base],
                           "survival_recomputed_from_npz_matches": True}

    # --- 2. reference comparison (fit line and markers separate) -------------
    fit, markers = load_reference()
    comp = {}
    for w in WAVEFORMS:
        n = np.asarray(formal[w]["n"])
        s = np.asarray(formal[w]["survival"])
        m = markers[w]
        comp[w] = {"vs_published_fit": metrics(n, s, *fit[w]),
                   "vs_experimental_markers": metrics(n, s, m["n"], m["s"])}
    report["reference_comparison"] = comp

    # --- 3. S(n), loss rates ---------------------------------------------------
    sn = {}
    losses = {}
    for w in WAVEFORMS:
        n = np.asarray(formal[w]["n"]); s = np.asarray(formal[w]["survival"])
        sn[w] = {k: float(np.interp(k, n, s)) for k in (2, 10, 20, 40, 60)}
        sn[w]["paper_fit"] = {k: (float(np.interp(k, *fit[w])) if fit[w][0].min() <= k <= fit[w][0].max() else None)
                              for k in (2, 10, 20, 40, 60)}
        losses[w] = {"first_transfer": float(1.0 - s[1]),
                     "early_per_transfer_n1_10": per_transfer_loss(n, s, 0, 10),
                     "late_per_transfer_n31_60": per_transfer_loss(n, s, 30, 60)}
    report["sn"] = sn
    report["loss_rates"] = losses

    # --- 4. paired timestep check ---------------------------------------------
    dtpairs = {}
    for w in WAVEFORMS:
        a = npz[w]["alive"]
        b = load_npz(DTCHECK[w])["alive"]
        n = np.asarray(formal[w]["n"])
        ds = np.abs(np.asarray(formal[w]["survival"]) - np.asarray(dtc[w]["survival"]))
        dtpairs[w] = {"max_abs_dS": float(ds.max()), "argmax_n": int(n[int(np.argmax(ds))]),
                      "final_label_disagreement": float(np.mean(a[:, -1] != b[:, -1])),
                      "curve_difference_within_0.02": bool(ds.max() <= 0.02)}
    report["paired_dt_0.0125_vs_0.025"] = dtpairs

    # --- 5. archived contrasts (recomputed against the CURRENT fit reference) --
    contrast = {}
    for w in WAVEFORMS:
        row = {"new_unified_S60": float(formal[w]["survival"][-1]),
               "new_unified_metrics_vs_fit": comp[w]["vs_published_fit"]}
        for tag, path in (("archived_T15_prefix_init", RUN14 / f"precision_T15_matched_{w}_N4096_dt0.0125.json"),
                          ("archived_percurve_scanT", RUN14 / f"precision_T{ARCHIVE_T[w]}_matched_{w}_N4096_dt0.0125.json")):
            if path.is_file():
                old = load_json(path)
                n_old = np.asarray(old["n"]); s_old = np.asarray(old["survival"])
                row[tag] = {"temperature_uK": old["temperature_uK"], "Tref_uK": old["parametric_Tref_uK"],
                            "S60": float(s_old[-1]),
                            "metrics_vs_current_fit": metrics(n_old, s_old, *fit[w])}
        contrast[w] = row
    report["archive_contrast"] = contrast

    # --- 6. figure: four subplots ---------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2), sharex=True, sharey=True)
    for ax, w in zip(axes.flat, WAVEFORMS):
        n = np.asarray(formal[w]["n"]); s = np.asarray(formal[w]["survival"])
        lo = np.asarray(formal[w]["wilson95_low"]); hi = np.asarray(formal[w]["wilson95_high"])
        ax.fill_between(n, lo, hi, color="#1f77b4", alpha=0.18, label="MC Wilson 95%")
        ax.plot(n, s, color="#1f77b4", lw=1.8, label="MC unified T0=Tref=15 uK")
        fn, fs = fit[w]
        ax.plot(fn, fs, color="#7f7f7f", lw=1.3, ls="--", label="paper published fit")
        m = markers[w]
        ax.errorbar(m["n"], m["s"], yerr=(m["elo"], m["ehi"]), fmt="D", ms=5,
                    color="k", capsize=3, lw=1, label="paper data points")
        c = comp[w]["vs_published_fit"]
        ax.set_title(f"{LABEL[w]}   RMSE={c['rmse']:.3f}  max|dev|={c['max_abs_dev']:.3f} @ n={c['argmax_n']}",
                     fontsize=11)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, 60); ax.set_ylim(0, 1.02)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    for ax in axes[1]: ax.set_xlabel("one-way transfers n")
    for ax in axes[:, 0]: ax.set_ylabel("survival S(n)")
    fig.suptitle("Unified Fig. 6d rerun: common T0=15 uK, Tref=15 uK, site-bound init, "
                 "N=4096, dt=0.0125 us, seed=23003, matched arm", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(HERE / "fig6d_unified_T15_comparison.png", dpi=300)

    # --- 7. CSV outputs ---------------------------------------------------------
    with (HERE / "survival_curves_unified_T15.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["waveform", "n", "survivors", "survival", "wilson95_low", "wilson95_high"])
        for w in WAVEFORMS:
            d = formal[w]
            for i, n in enumerate(d["n"]):
                wr.writerow([w, n, d["survival_counts"][i], f"{d['survival'][i]:.6f}",
                             f"{d['wilson95_low'][i]:.6f}", f"{d['wilson95_high'][i]:.6f}"])
    with (HERE / "metrics_vs_reference.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["waveform", "comparison_target", "n_min", "n_max", "n_points",
                     "rmse", "max_abs_dev", "argmax_n", "fraction_within_0.02"])
        for w in WAVEFORMS:
            for target, key in (("published_fit_line", "vs_published_fit"),
                                ("experimental_markers", "vs_experimental_markers")):
                c = comp[w][key]
                wr.writerow([w, target, c["n_min"], c["n_max"], c["n_points"],
                             f"{c['rmse']:.5f}", f"{c['max_abs_dev']:.5f}", c["argmax_n"],
                             f"{c['fraction_within_band']:.3f}"])
    with (HERE / "sn_and_loss_rates.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["waveform", "S2", "S10", "S20", "S40", "S60",
                     "fit_S2", "fit_S10", "fit_S20", "fit_S40", "fit_S60",
                     "first_transfer_loss", "early_per_transfer_loss_n1_10", "late_per_transfer_loss_n31_60"])
        for w in WAVEFORMS:
            s, l, pf = sn[w], losses[w], sn[w]["paper_fit"]
            wr.writerow([w] + [f"{s[k]:.4f}" for k in (2, 10, 20, 40, 60)]
                        + ["" if pf[k] is None else f"{pf[k]:.4f}" for k in (2, 10, 20, 40, 60)]
                        + [f"{l['first_transfer']:.5f}", f"{l['early_per_transfer_n1_10']:.5f}",
                           f"{l['late_per_transfer_n31_60']:.5f}"])
    with (HERE / "paired_dt_check.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["waveform", "max_abs_dS_0.0125_vs_0.025", "argmax_n",
                     "final_label_disagreement", "curve_difference_within_0.02"])
        for w in WAVEFORMS:
            d = dtpairs[w]
            wr.writerow([w, f"{d['max_abs_dS']:.5f}", d["argmax_n"],
                         f"{d['final_label_disagreement']:.4f}", d["curve_difference_within_0.02"]])
    with (HERE / "archive_contrast.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["waveform", "case", "temperature_uK", "Tref_uK", "S60",
                     "rmse_vs_current_fit", "max_abs_dev_vs_current_fit"])
        for w in WAVEFORMS:
            for tag in ("archived_T15_prefix_init", "archived_percurve_scanT"):
                if tag in contrast[w]:
                    c = contrast[w][tag]; m = c["metrics_vs_current_fit"]
                    wr.writerow([w, tag, c["temperature_uK"], c["Tref_uK"], f"{c['S60']:.4f}",
                                 f"{m['rmse']:.5f}", f"{m['max_abs_dev']:.5f}"])
            c = contrast[w]
            m = c["new_unified_metrics_vs_fit"]
            wr.writerow([w, "new_unified_T15", 15, 15, f"{c['new_unified_S60']:.4f}",
                         f"{m['rmse']:.5f}", f"{m['max_abs_dev']:.5f}"])

    # --- 8. source hashes --------------------------------------------------------
    sources = [ROOT / "configs/continuous_fig6d_unified_T15.yaml",
               ROOT / "configs/level2c_pickup_survival.yaml",
               ROOT / "outputs/level2c_pickup_survival/scan_stage1/fig6c_ml_trajectory_digitized.json",
               ROOT / "reference/fig6d_vector_v2/published_fit_reference.csv",
               ROOT / "reference/fig6d_vector_v2/experimental_marker_centers.csv",
               ROOT / "src/continuous_transfer/engine.py",
               ROOT / "src/continuous_transfer/protocol.py",
               ROOT / "src/continuous_transfer/initialization.py",
               ROOT / "src/continuous_transfer/cli.py"]
    report["source_sha256_16"] = {str(p.relative_to(ROOT)): sha256(p) for p in sources}
    (HERE / "analysis_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")

    print(json.dumps({"integrity": report["integrity"],
                      "paired_dt": {w: dtpairs[w]["max_abs_dS"] for w in WAVEFORMS},
                      "rmse_fit": {w: round(comp[w]["vs_published_fit"]["rmse"], 4) for w in WAVEFORMS},
                      "S60": {w: round(sn[w][60], 4) for w in WAVEFORMS}},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
