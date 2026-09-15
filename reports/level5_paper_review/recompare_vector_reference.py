"""Reassess archived Monte Carlo outputs against independent PDF vector data."""
from pathlib import Path
import copy
import hashlib
import json

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROJECT = ROOT / "simulation/level2_joint_transfer"
VECTOR_FILE = PROJECT / "reference/fig6d_vector_v2/fig6d_vectors.json"
OLD_FILE = HERE.parent / "stage_review/output/AOD_SLM_continuous_level2c_345_2026-09-14_v10.evidence.json"
DATA = PROJECT / "outputs/continuous_level2c_345/run_20260914"
OUT = HERE / "output"
KEYS = ["manual_200us", "manual_400us", "manual_600us", "ml_400us"]


def metrics(row, n, reference, basis):
    n, reference = np.asarray(n, dtype=int), np.asarray(reference)
    simulated = np.asarray(row["survival"])[n]
    delta = simulated - reference
    return {"basis": basis, "n": n.tolist(), "reference": reference.tolist(),
            "simulated": simulated.tolist(), "deviation": delta.tolist(),
            "max_abs_deviation": float(max(abs(delta))), "worst_n": int(n[np.argmax(abs(delta))]),
            "rmse": float(np.sqrt(np.mean(delta**2))), "matches_reference": bool(np.all(abs(delta) <= .02)),
            "fraction_within_band": float(np.mean(abs(delta) <= .02))}


def main():
    vectors = json.loads(VECTOR_FILE.read_text("utf-8"))
    old = json.loads(OLD_FILE.read_text("utf-8"))
    evidence = copy.deepcopy(old)
    def revise(row):
        s = vectors["series"][row["waveform"]]
        row["comparison_paper_fit"] = metrics(row, s["comparison_n"], s["comparison_fit"], "published PDF fit path; prior n range retained")
        row["comparison_paper_markers"] = metrics(row, [m["n"] for m in s["markers"]], [m["survival"] for m in s["markers"]], "published marker centers at measured n only")
        return row
    for temp in evidence["base"].values():
        for rows in temp.values():
            for row in rows:
                revise(row)
    for rows in evidence["precision"].values():
        for row in rows:
            revise(row)
    for row in evidence["selected"]:
        revise(row)
    cases = []
    for file in sorted(DATA.glob("*.json")):
        row = json.loads(file.read_text("utf-8"))
        if "case_id" not in row:
            continue
        row["result_file"] = str(file.relative_to(ROOT))
        cases.append(revise(row))
    scan = []
    for temp in [5, 10, 15, 20, 25, 30, 35, 40]:
        rows = [r for r in cases if r["stage"] == "scan" and r["temperature_uK"] == temp
                and r["arm"] == "matched" and r["shots"] == 256 and r["dt_us"] == .05]
        assert len(rows) == 4
        scan.append({"temperature_uK": temp,
                     "pooled_rmse": float(np.sqrt(np.mean([r["comparison_paper_fit"]["rmse"]**2 for r in rows]))),
                     "all_curves_pass": all(r["comparison_paper_fit"]["matches_reference"] for r in rows)})
    evidence["temperature_scan_original"] = evidence["temperature_scan"]
    evidence["temperature_scan"] = scan
    evidence["reference_vectors"] = vectors
    evidence["reference_revision"] = {"file": str(VECTOR_FILE.relative_to(ROOT)),
        "sha256": hashlib.sha256(VECTOR_FILE.read_bytes()).hexdigest(),
        "old_evidence": str(OLD_FILE.relative_to(ROOT)),
        "old_sha256": hashlib.sha256(OLD_FILE.read_bytes()).hexdigest(),
        "simulation_rerun": False, "selected_parameters": "Previously frozen T = 5/15/25/25 retained; no reselection for this presentation",
        "acceptance": "Same ±0.02 comparison tolerance, now around the actual published fit path. This tolerance is not the paper's error bar or a confidence interval."}
    summary = []
    for row in evidence["selected"]:
        cp = row["comparison_paper_fit"]
        n, ref = np.array(cp["n"]), np.array(cp["reference"])
        lo, hi = np.array(row["wilson95_low"])[n], np.array(row["wilson95_high"])[n]
        summary.append({"waveform": row["waveform"], "temperature_uK": row["temperature_uK"],
                        "old_max_deviation": row["comparison_raw"]["max_abs_deviation"],
                        "fit_max_deviation": cp["max_abs_deviation"],
                        "marker_max_deviation": row["comparison_paper_markers"]["max_abs_deviation"],
                        "fit_endpoint": cp["reference"][-1], "sim_endpoint": cp["simulated"][-1],
                        "worst_n": cp["worst_n"], "fit_pass": cp["matches_reference"],
                        "points_with_MC_95_interval_outside_fit_tolerance": int(((lo>ref+.02)|(hi<ref-.02)).sum())})
    changes = {}
    for row in evidence["selected"]:
        delta = np.array(row["comparison_paper_fit"]["reference"]) - row["comparison_raw"]["reference"]
        changes[row["waveform"]] = {"max_old_to_fit_change": float(max(abs(delta))),
                                  "rms_old_to_fit_change": float(np.sqrt(np.mean(delta**2)))}
    evidence["recomparison_summary"] = summary
    evidence["reference_change"] = changes
    evidence["verdict"] = "Corrected published fit references do not pass the same ±0.02 band for any of the four frozen high-precision simulations."
    # Legacy confidence was derived from flawed raster points, so preserve under
    # an explicit historical label instead of continuing to publish it as current.
    evidence["confidence_original"] = evidence.pop("confidence")
    evidence["confidence"] = summary
    OUT.mkdir(exist_ok=True)
    (OUT / "Level5_vector_reference_v2.evidence.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "Level5_vector_reference_v2.recomparison.json").write_text(json.dumps({"cases": [{"case_id": r["case_id"], "result_file": r["result_file"], "comparison_paper_fit": r["comparison_paper_fit"], "comparison_paper_markers": r["comparison_paper_markers"]} for r in cases]}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"summary": summary, "shared_best": min(scan,key=lambda r:r["pooled_rmse"]), "changes": changes, "reassessed_cases": len(cases)}, ensure_ascii=False))
    return evidence


if __name__ == "__main__":
    main()
