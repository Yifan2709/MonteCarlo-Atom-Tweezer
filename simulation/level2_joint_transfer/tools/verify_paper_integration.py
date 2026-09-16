"""验收核对工具：合同 ↔ 运行目录 ↔ 状态清单 三方核对。

用法：
  python tools/verify_paper_integration.py --contract configs/paper_integration/contract \
      --run-dir outputs/paper_integration/run_all [--tier core]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

CHECKS = {
    "B0": ["preflight.json"],
    "M1": ["preflight.json"],
    "R1": ["scan_rb2022.json"],
    "R2": ["scan_rb2022.json"],
    "R3": ["validation_rb2022.json"],
    "Y1": ["scan_yb2026.json"],
    "Y2": ["validation_yb2026.json"],
    "Y3": ["validation_yb2026.json"],
    "Y4": ["validation_yb2026.json"],
    "Y5": ["validation_yb2026.json"],
    "Y6": ["validation_yb2026.json"],
    "I1": ["integration.json"],
    "V1": ["requirements_status.csv", "report.md"],
}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--contract", type=Path, required=True)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--run-dir-yb", type=Path, default=None,
                   help="Y 侧运行目录（缺省用 <run-dir>_yb）")
    p.add_argument("--tier", default="core", choices=["core", "full_paper"])
    a = p.parse_args(argv)

    run_yb = a.run_dir_yb or a.run_dir.with_name(a.run_dir.name + "_yb")
    contract = json.loads((a.contract / "acceptance_contract.json")
                          .read_text("utf-8"))
    coverage = list(csv.DictReader(
        (a.contract / "paper_coverage.csv").open(encoding="utf-8")))
    problems: list[str] = []
    rows = []
    for req in coverage:
        rid = req["requirement_id"]
        files = CHECKS.get(rid, [])
        search = [a.run_dir, run_yb] if rid.startswith("Y") else [a.run_dir]
        present = [any((d / f).exists() for d in search) for f in files]
        evidence_ok = all(present) if files else False
        rows.append({
            "requirement_id": rid,
            "scope": req["scope"],
            "implementation": "implemented" if evidence_ok else "missing",
            "numeric": "ran" if evidence_ok else "not_run",
            "paper_match": "see_report" if evidence_ok else "missing",
            "independent_review": "SELF_CHECKED",
            "evidence": ";".join(str(a.run_dir / f) for f in files),
        })
        if not evidence_ok:
            problems.append(f"{rid}: 缺证据文件 {files}")
    manifest_files = list(a.run_dir.glob("run_manifest.json")) \
        + list(run_yb.glob("run_manifest.json"))
    if not manifest_files:
        problems.append("缺 run_manifest.json")
    for m in manifest_files:
        try:
            man = json.loads(m.read_text("utf-8"))
            if man.get("contract_sha256") != contract.get("contract_sha256"):
                problems.append(
                    f"{m}: 合同哈希不一致 {man.get('contract_sha256')}")
        except json.JSONDecodeError as e:
            problems.append(f"{m}: manifest 解析失败 {e}")

    out_csv = a.run_dir / "verification_matrix.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"核对完成：{len(rows)} 项，问题 {len(problems)} 个")
    for pr in problems:
        print(" -", pr)
    print("矩阵 ->", out_csv)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
