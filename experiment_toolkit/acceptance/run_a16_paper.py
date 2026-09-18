"""A16: optional paper regression against the archived research study.

Uses the research repository ONLY as an explicit user-supplied resource:
the ML vector nodes (Fig. 6c) and the archived Joint-model survival curves
(2026-09-16 calibration study, 8192 shots). The toolkit re-runs the four
matched short-transfer curves with the SAME locked physics (diffusion,
T=19 uK, waist 1.11 um, AOD RIN 5.75e-9 /Hz, pointing 1e-24 m^2/Hz, SLM
noise zero, lensing off, jitter 10 nm, alignment 50 nm, paper cubic manual
48%, ML replay) and compares against the archived curves.

Criteria fixed BEFORE running:
  per curve: >= 90% of common n satisfy |dS| <= 3.5 * sigma_mc(S_mean)
             AND max |dS| <= 0.06
where sigma_mc combines both pools' binomial errors. This compares the
MIGRATED IMPLEMENTATION to the SAME MODEL's archived results — passing the
paper's own data is explicitly NOT a criterion here.

Notes on legitimate statistical differences (documented, not fitted away):
independent per-atom noise draws (vectorized streams) vs the archived
per-shot legacy RNG; per-occurrence jitter draws. Both are distributionally
equivalent; curves must still agree within Monte Carlo error.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOOLKIT = HERE.parent
sys.path.insert(0, str(TOOLKIT / "src"))

import tweezer_experiment as te  # noqa: E402
from tweezer_experiment.adapters import paper as paper_adapter  # noqa: E402
from tweezer_experiment import units as U  # noqa: E402

RESEARCH = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    TOOLKIT.parent / "simulation/level2_joint_transfer")
CALIB_DIR = TOOLKIT.parent / "reports/movement_effects_v1_calibration_20260916"
ML_NODES = CALIB_DIR / "ml_vector_nodes.json"
ARCHIVED_CURVES = CALIB_DIR / "survival_curves.csv"

SHOTS = 2048
DT_US = 0.0125
NOISE_DT_US = 0.05
SEED = 20260917
Z95 = 1.959963984540054

WAVEFORMS = {"manual_200us": 200.0, "manual_400us": 400.0, "manual_600us": 600.0}


def archived_joint_curves():
    rows = list(csv.DictReader(open(ARCHIVED_CURVES, encoding="utf-8-sig")))
    out = {}
    for row in rows:
        if row["group"] != "Joint model":
            continue
        out.setdefault(row["waveform"], {"n": [], "s": [], "shots": int(row["shots"])})
        out[row["waveform"]]["n"].append(int(row["n"]))
        out[row["waveform"]]["s"].append(float(row["survival"]))
    return {w: {"n": np.asarray(d["n"]), "s": np.asarray(d["s"]), "shots": d["shots"]}
            for w, d in out.items()}


def bundle_for(model: str):
    raw = {
        "device": {
            "species": "Cs133",
            "slm": {"wavelength_nm": 1061, "waist_um": 1.11, "depth_uK": 140},
            "aod": {"wavelength_nm": 1055, "waist_um": 1.11, "depth_uK": 280},
            "lensing": {"mode": "off"}},
        "preparation": {
            "temperature_uK": 19.0,
            "site_disorder": True,
            "slm_depth_rel_sigma": 0.114,
            "aod_depth_rel_sigma": 0.02,
            "aod_waist_rel_sigma": 0.02,
            "alignment_sigma_um": 0.05,
            "shot_jitter_alignment_sigma_nm": 10.0},
        "noise": {
            "model": "diffusion",
            "aod_rin_psd_per_hz": 5.75e-9,
            "aod_pointing_psd_m2_per_hz": 1.0e-24,
            "slm_rin_psd_per_hz": 0.0,
            "slm_pointing_psd_m2_per_hz": 0.0,
            "recoil_per_joule": 0.0},
        "numerics": {"dt_us": DT_US, "noise_dt_us": NOISE_DT_US,
                     "shots": SHOTS, "seed": SEED},
    }
    return raw


def protocol_for(name: str, ml_csv: Path | None):
    if name == "ml_400us":
        return paper_adapter.matched_roundtrip_spec(
            pickup_duration_us=400.0, ml_csv=str(ml_csv))
    return paper_adapter.matched_roundtrip_spec(
        pickup_duration_us=WAVEFORMS[name], family="paper_manual")


def main():
    if not ML_NODES.is_file() or not ARCHIVED_CURVES.is_file():
        print(json.dumps({"status": "NOT_RUN",
                          "reason": "paper resources not provided"}))
        return 0
    archived = archived_joint_curves()
    ml_csv = HERE / "evidence/a16_ml_replay.csv"
    nodes = paper_adapter.load_ml_nodes(ML_NODES)
    paper_adapter.write_ml_control_csv(nodes, ml_csv, aod_depth_uK=280.0)
    paper_adapter.reverse_control_csv(ml_csv, Path(str(ml_csv) + ".reversed.csv"))

    report = {"shots": SHOTS, "dt_us": DT_US, "noise_dt_us": NOISE_DT_US, "seed": SEED,
              "criteria": ">=90% of common n within 3.5*sigma_mc AND max|dS|<=0.06 per curve",
              "curves": {}}
    all_pass = True
    from dataclasses import replace as _replace
    placeholder = {"segments": [{"id": "x", "type": "static", "duration_us": 1,
                                 "aod_center_um": [0, 0], "aod_depth_uK": 0,
                                 "slm_factor": 1.0, "observe": "slm"}],
                   "sequence": ["x"]}
    for name in ("manual_200us", "manual_400us", "manual_600us", "ml_400us"):
        raw = bundle_for(name)
        raw["protocol"] = placeholder          # parse requirement; replaced below
        protocol = protocol_for(name, ml_csv if name == "ml_400us" else None)
        bundle = _replace(te.ExperimentBundle.parse(raw), protocol=protocol)
        result = te.simulate(bundle, protocol)
        arch = archived[name]
        n_arch, s_arch = arch["n"], arch["s"]
        common = n_arch[n_arch <= 60]
        sim_at = np.interp(common, np.arange(len(result.survival)), result.survival)
        ref_at = np.interp(common, n_arch, s_arch)
        d = sim_at - ref_at
        p_bar = np.clip((ref_at + sim_at) / 2, 1e-6, 1 - 1e-6)
        sigma = np.sqrt(p_bar * (1 - p_bar) / SHOTS + p_bar * (1 - p_bar) / arch["shots"])
        frac_within = float(np.mean(np.abs(d) <= 3.5 * np.maximum(sigma, 1e-4)))
        max_d = float(np.max(np.abs(d)))
        passed = frac_within >= 0.90 and max_d <= 0.06
        all_pass = all_pass and passed
        report["curves"][name] = {
            "n_compared": int(len(common)), "frac_within_3.5sigma": round(frac_within, 4),
            "max_abs_dS": round(max_d, 4),
            "argmax_n": int(common[int(np.argmax(np.abs(d)))]) if len(common) else None,
            "final_S_sim": float(result.survival[-1]), "final_S_archived": float(s_arch[-1]),
            "passed": passed}
    report["status"] = "PASS" if all_pass else "FAIL"
    out = HERE / "evidence/a16_paper_regression.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
