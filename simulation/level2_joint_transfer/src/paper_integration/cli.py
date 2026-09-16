"""paper_integration CLI：统一入口（注册于 simulation_workflow LEVELS）。

用法：
  python -m paper_integration.cli list
  python -m paper_integration.cli contract
  python -m paper_integration.cli run --stage preflight [--output-dir D]
  python -m paper_integration.cli run paper-rb2022 --stage scan --shots N --seed S
  python -m paper_integration.cli run paper-yb2026 --stage scan|validation ...
  python -m paper_integration.cli run --stage integration
"""
from __future__ import annotations

import argparse
from pathlib import Path

from . import contract_data, runs


def main(argv=None):
    p = argparse.ArgumentParser(prog="paper_integration")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    sub.add_parser("contract")
    r = sub.add_parser("run")
    r.add_argument("paper", nargs="?", choices=["paper-rb2022", "paper-yb2026"],
                   default=None)
    r.add_argument("--stage", required=True,
                   choices=["preflight", "scan", "validation", "integration",
                            "all"])
    r.add_argument("--shots", type=int, default=1024)
    r.add_argument("--trials", type=int, default=200000)
    r.add_argument("--seed", type=int, default=20260915)
    r.add_argument("--output-dir", type=Path, default=None)
    args = p.parse_args(argv)

    if args.cmd == "list":
        print("paper-rb2022: scan/validation | paper-yb2026: scan/validation "
              "| 无 paper: preflight/integration")
        return 0
    if args.cmd == "contract":
        print(contract_data.write_contract())
        return 0

    base = args.output_dir or (runs.OUTPUT_ROOT / f"run_{args.paper or 'base'}")
    base.mkdir(parents=True, exist_ok=True)
    if args.stage == "preflight":
        out = runs.stage_preflight(base)
        print("preflight ok:", out["preflight_ok"])
    elif args.stage == "scan":
        if args.paper == "paper-rb2022":
            runs.stage_scan_rb(base, args.shots, args.seed)
        elif args.paper == "paper-yb2026":
            runs.stage_scan_yb(base, args.shots, args.seed)
        else:
            raise SystemExit("scan 需要 paper-rb2022 或 paper-yb2026")
        print("scan done ->", base)
    elif args.stage == "validation":
        which = {"paper-rb2022": "rb2022",
                 "paper-yb2026": "yb2026"}.get(args.paper, "all")
        runs.stage_validation(base, args.trials, args.seed, which)
        print("validation done ->", base)
    elif args.stage == "integration":
        runs.stage_integration(base)
        print("integration done ->", base)
    elif args.stage == "all":
        runs.stage_preflight(base)
        runs.stage_scan_rb(base, args.shots, args.seed)
        runs.stage_scan_yb(base, args.shots, args.seed + 1)
        runs.stage_validation(base, args.trials, args.seed + 2, "rb2022")
        runs.stage_validation(base, args.trials, args.seed + 3, "yb2026")
        runs.stage_integration(base)
        print("all done ->", base)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
