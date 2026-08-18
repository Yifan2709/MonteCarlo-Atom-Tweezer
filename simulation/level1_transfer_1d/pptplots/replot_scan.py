"""Re-render scan figures from the saved CSV (no re-run needed).

Usage:  python replot_scan.py
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

from run_temperature_scan import make_figures, OUT_DIR


def load_summary():
    with open(OUT_DIR / "summary_by_temperature.csv", newline="") as f:
        return list(csv.DictReader(f))


def load_raw():
    with open(OUT_DIR / "mc_shots_by_temperature.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        rr = {}
        for k, v in r.items():
            if v == "" or v is None:
                rr[k] = None
            elif k in ("shot_id",):
                rr[k] = int(v)
            elif k == "captured_in_slm":
                rr[k] = v == "True"
            elif k == "temperature_uK":
                rr[k] = float(v)
            else:
                rr[k] = float(v)
        out.append(rr)
    return out


def main():
    summary_raw = load_summary()
    # cast numeric fields to float (CSV stores strings)
    summary = []
    for s in summary_raw:
        ss = {}
        for k, v in s.items():
            try:
                ss[k] = float(v)
            except (ValueError, TypeError):
                ss[k] = v
        ss["shots"] = int(ss["shots"])
        ss["captured"] = int(ss["captured"])
        summary.append(ss)

    all_records = load_raw()
    meta = json.loads((OUT_DIR / "scan_meta.json").read_text())
    print(f"Loaded {len(all_records)} shots, {len(summary)} temps from {OUT_DIR}")
    make_figures(summary, all_records, meta)
    print("Done.")


if __name__ == "__main__":
    main()
