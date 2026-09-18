"""tweezer-experiment command line (same application layer as the API).

Subcommands: validate | config | simulate | compare | optimize | export |
report | selftest. All commands accept plain experiment YAML files and do
not require any research-repository artifacts.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from . import __version__
from .errors import ToolkitError
from .provenance import write_json, json_safe
from .schemas import ExperimentBundle, validate_experiment, effective_config
from .simulation import simulate, save_result
from .search import search_controls
from .export import export_controls


def _load(path: str) -> ExperimentBundle:
    return ExperimentBundle.load(path)


def cmd_validate(args) -> int:
    try:
        bundle = _load(args.config)
    except ToolkitError as exc:
        report = type("R", (), {"ok": False, "errors": [f"{exc}"], "warnings": [],
                                "capability": {}, "effective_config": {}})()
        bundle = None
    else:
        report = validate_experiment(bundle)
    if args.json:
        print(json.dumps(json_safe(report.__dict__), indent=2, ensure_ascii=False))
    else:
        print(f"valid: {report.ok}")
        for err in report.errors:
            print(f"  ERROR {err}")
        for warn in report.warnings:
            print(f"  warn  {warn}")
        for key, value in report.capability.items():
            print(f"  {key}: {value}")
    return 0 if report.ok else 1


def cmd_config(args) -> int:
    bundle = _load(args.config)
    report = validate_experiment(bundle)
    if not report.ok:
        print(json.dumps(json_safe(report.__dict__), indent=2, ensure_ascii=False))
        return 1
    print(json.dumps(report.effective_config, indent=2, ensure_ascii=False))
    return 0


def cmd_simulate(args) -> int:
    bundle = _load(args.config)
    started = time.perf_counter()
    result = simulate(bundle, shots=args.shots, seed=args.seed)
    result.run.elapsed_s = time.perf_counter() - started
    out = Path(args.out) if args.out else Path(f"tweezer_runs/{result.run.run_id}")
    path = save_result(out, result)
    summary = {"run_id": result.run.run_id, "output": str(path),
               "elapsed_s": result.run.elapsed_s,
               "checkpoints": [{"label": label, "trap": trap,
                                "survival": float(s), "count": int(c)}
                               for (label, trap, _t), s, c in
                               zip(result.checkpoints, result.survival, result.survival_counts)],
               "initial_unbound": result.initial_unbound,
               "loss_histogram": result.loss_histogram,
               "warnings": result.warnings}
    print(json.dumps(json_safe(summary), indent=2, ensure_ascii=False))
    return 0


def cmd_compare(args) -> int:
    a, b = _load(args.config_a), _load(args.config_b)
    ra, rb = simulate(a, shots=args.shots), simulate(b, shots=args.shots)
    n = min(len(ra.survival), len(rb.survival))
    diff = {
        "config_a": args.config_a, "config_b": args.config_b,
        "common_checkpoints": int(n),
        "max_abs_dS": float(max((abs(ra.survival[i] - rb.survival[i]) for i in range(n)), default=0.0)),
        "final_capture": [float(ra.survival[-1]), float(rb.survival[-1])],
        "seeds": [ra.run.seeds["master"], rb.run.seeds["master"]]}
    print(json.dumps(json_safe(diff), indent=2, ensure_ascii=False))
    return 0


def cmd_optimize(args) -> int:
    bundle = _load(args.config)
    started = time.perf_counter()
    result = search_controls(bundle)
    result.run.elapsed_s = time.perf_counter() - started
    out = Path(args.out) if args.out else Path(f"tweezer_runs/{result.run.run_id}_search")
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "search_result.json", result.to_dict())
    rec = result.recommended
    print(json.dumps(json_safe({
        "candidates": len(result.candidates),
        "ok": sum(c.status == "ok" for c in result.candidates),
        "infeasible": sum(c.status == "infeasible" for c in result.candidates),
        "failed": sum(c.status == "failed" for c in result.candidates),
        "recommended": None if rec is None else rec.to_dict(),
        "notes": result.notes, "output": str(out / "search_result.json")}),
        indent=2, ensure_ascii=False))
    return 0


def cmd_export(args) -> int:
    bundle = _load(args.config)
    manifest = export_controls(bundle, args.out, verify=not args.no_verify)
    print(json.dumps(json_safe({
        "physical_table": manifest["physical_table"]["path"],
        "switch_events": len(manifest["physical_table"]["switch_events"]),
        "overshoot_reports": len(manifest["physical_table"]["audit"].get("overshoot", [])),
        "device_instructions": (manifest["device_instructions"] or {}).get("path")
        if manifest["device_instructions"] else None,
        "limits_ok": (manifest["device_instructions"] or {}).get("limits_ok", None)
        if manifest["device_instructions"] else None,
        "readback": manifest["readback"]}), indent=2, ensure_ascii=False))
    return 0


def cmd_report(args) -> int:
    path = Path(args.dir)
    result_file = path / "result.json"
    if not result_file.is_file():
        print(f"no result.json under {path}", file=sys.stderr)
        return 1
    data = json.loads(result_file.read_text(encoding="utf-8"))
    print(json.dumps(data, indent=2, ensure_ascii=False))
    return 0


def cmd_selftest(args) -> int:
    from .units import (temperature_uK_to_joule, joule_to_temperature_uK,
                        hz_to_rad_per_s, rad_per_s_to_hz, um_to_m, us_to_s)
    checks = [
        ("uK->J->uK roundtrip", abs(joule_to_temperature_uK(temperature_uK_to_joule(140.0)) - 140.0) < 1e-9),
        ("J anchor 1 uK = k_B*1e-6", abs(temperature_uK_to_joule(1.0) - 1.380649e-29) < 1e-40),
        ("Hz<->rad/s roundtrip", abs(rad_per_s_to_hz(hz_to_rad_per_s(24611.0)) - 24611.0) < 1e-6),
        ("um->m", abs(um_to_m(2.4) - 2.4e-6) < 1e-15),
        ("us->s", abs(us_to_s(400.0) - 4e-4) < 1e-15),
    ]
    ok = all(passed for _, passed in checks)
    for name, passed in checks:
        print(f"  {'PASS' if passed else 'FAIL'} {name}")
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="tweezer-experiment",
        description="Offline experiment planning for tweezer transport (phase 1: mechanics).")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate", help="validate an experiment YAML and classify capabilities")
    p.add_argument("config")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_validate)

    p = sub.add_parser("config", help="print the effective configuration with sources")
    p.add_argument("config")
    p.set_defaults(fn=cmd_config)

    p = sub.add_parser("simulate", help="run a forward prediction")
    p.add_argument("config")
    p.add_argument("--out")
    p.add_argument("--shots", type=int)
    p.add_argument("--seed", type=int)
    p.set_defaults(fn=cmd_simulate)

    p = sub.add_parser("compare", help="simulate two configs and diff outcomes")
    p.add_argument("config_a")
    p.add_argument("config_b")
    p.add_argument("--shots", type=int)
    p.set_defaults(fn=cmd_compare)

    p = sub.add_parser("optimize", help="run the constrained candidate search")
    p.add_argument("config")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_optimize)

    p = sub.add_parser("export", help="export control tables and verify read-back")
    p.add_argument("config")
    p.add_argument("out")
    p.add_argument("--no-verify", action="store_true")
    p.set_defaults(fn=cmd_export)

    p = sub.add_parser("report", help="print a stored result")
    p.add_argument("dir")
    p.set_defaults(fn=cmd_report)

    p = sub.add_parser("selftest", help="fast built-in unit/anchor checks")
    p.set_defaults(fn=cmd_selftest)

    args = parser.parse_args(argv)
    try:
        return args.fn(args)
    except ToolkitError as exc:
        print(f"error [{exc.code}] {exc.field}: {exc.message}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
