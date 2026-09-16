"""统一调用已有 CLI。list/check 不启动 Monte Carlo，也不导入数值库。"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LEVELS = {
    "full": {
        "title": "Level 2C→3/4/5 连续全流程",
        "module": "continuous_transfer.cli",
        "config": "continuous_level2c_345.yaml",
        "output": "continuous_level2c_345",
        "stages": ["preflight", "scan", "validation", "convergence", "all"],
        "unit": "60 短转移计数=30 连续往返；matched/extended 两种协议分别验收",
        "depends_on": ["Level 2C 原始波形与加热", "Level 3 三维透镜场", "有限脉冲 SU(2)"],
    },
    "2c": {
        "title": "一维短距离 pick-up/drop-off 存活率",
        "module": "level2c_pickup_survival.level2c_cli",
        "config": "level2c_pickup_survival.yaml",
        "output": "level2c_pickup_survival",
        "stages": ["preflight", "pickup", "roundtrip", "ml", "survival_scan", "all"],
        "unit": "单程 2.4 μm；60 单程 = 30 短往返；无长运输",
        "depends_on": [],
    },
    "3": {
        "title": "三维 AOD 长距离运输与柱面透镜",
        "module": "level3_3d_transport_lensing.level3_cli",
        "config": "level3_3d_transport.yaml",
        "output": "level3_3d_transport_lensing",
        "stages": ["preflight", "potential_scan", "coarse_mc", "validate", "sensitivity", "compensation", "all"],
        "unit": "单次 270/510/610 μm；无拾取、放回或重复往返",
        "depends_on": [],
    },
    "4": {
        "title": "完整往返与半经典相干性",
        "module": "level4_roundtrip_coherence.level4_cli",
        "config": "level4_roundtrip.yaml",
        "output": "level4_end_to_end_roundtrip",
        "stages": ["preflight", "single_round", "repeated_rounds", "ablations", "sensitivity", "all"],
        "unit": "每轮 2.4 + 375 + 375 + 2.4 = 754.8 μm；默认 10/40 轮",
        "depends_on": ["Level 2 冻结波形（Protocol B）", "Level 3 冻结透镜场景"],
    },
    "5": {
        "title": "有限脉冲、PTM 与多站点 RB/IRB",
        "module": "level5_finite_pulse_irb.level5_cli",
        "config": "level5_finite_pulse_irb.yaml",
        "output": "level5_finite_pulse_irb",
        "stages": ["preflight", "channel", "reference_rb", "interleaved_rb", "array_extension", "sensitivity", "all"],
        "unit": "默认 M 为完整组合往返数；每次独立抽取冻结轨迹",
        "depends_on": ["Level 4 协议与轨迹池"],
    },
,
    "paper-rb2022": {
        "title": "R 论文（Bluvstein 2022）Rb 纠缠运输与线路接入",
        "module": "paper_integration.cli",
        "config": "paper_integration/contract/acceptance_contract.json",
        "output": "paper_integration/run_paper-rb2022",
        "stages": ["preflight", "scan", "validation", "integration", "all"],
        "unit": "运输扫描按单次移动；线路 MC 按线路 trial",
        "depends_on": ["合同 configs/paper_integration/contract"],
    },
    "paper-yb2026": {
        "title": "Y 论文（Zhang 2026）Yb 擦除转换与逻辑比特接入",
        "module": "paper_integration.cli",
        "config": "paper_integration/contract/acceptance_contract.json",
        "output": "paper_integration/run_paper-yb2026",
        "stages": ["preflight", "scan", "validation", "integration", "all"],
        "unit": "运输扫描按单程；逻辑按线路 trial；三口径分母见合同",
        "depends_on": ["合同 configs/paper_integration/contract"],
    },
}


def inventory(root=PROJECT_ROOT):
    """配置、代码和历史指标分别核查；文件存在不等于物理验收通过。"""
    rows = []
    for level, spec in LEVELS.items():
        files = {
            "source": root / "src" / (spec["module"].replace(".", "/") + ".py"),
            "config": root / "configs" / spec["config"],
            "historical_metrics": root / "outputs" / spec["output"] / "demo/metrics.json",
        }
        if level == "full":
            files["historical_metrics"] = (root / "outputs/continuous_level2c_345/run_20260914/"
                                           "validation_T5_matched_manual_400us_N512_dt0.05.json")
        details = {}
        for kind, file in files.items():
            details[kind] = {"path": str(file), "exists": file.is_file()}
            if file.is_file():
                details[kind]["sha256"] = hashlib.sha256(file.read_bytes()).hexdigest()
        rows.append({"level": level, **spec, "files": details,
                     "runnable_files_present": all(details[k]["exists"] for k in ("source", "config")),
                     "result_status": ("continuous_results_rerun_20260914" if level == "full"
                                       and files["historical_metrics"].is_file() else
                                       "historical_irb_and_array_require_rerun" if level == "5"
                                       else "historical_not_rerun_in_report_update")})
    return rows


def command_for(level, stage, config=None, output_dir=None, root=PROJECT_ROOT):
    spec = LEVELS[level]
    if stage not in spec["stages"]:
        raise ValueError(f"Level {level} 不支持 {stage}；可用阶段：{', '.join(spec['stages'])}")
    cfg = Path(config).resolve() if config else root / "configs" / spec["config"]
    if not cfg.is_file():
        raise FileNotFoundError(cfg)
    out = (Path(output_dir).resolve() if output_dir else root / "outputs" /
           spec["output"] / ("run_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")))
    return [sys.executable, "-m", spec["module"], "--config", str(cfg),
            "--stage", stage, "--output-dir", str(out)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("list", help="层级、操作单位和依赖关系")
    check = sub.add_parser("check", help="核查本地源代码、配置、历史指标及 SHA-256")
    check.add_argument("--json", action="store_true")
    run = sub.add_parser("run", help="运行单个层级，不自动拼接不同协议")
    run.add_argument("level", choices=LEVELS)
    run.add_argument("--stage", default="preflight")
    run.add_argument("--config")
    run.add_argument("--output-dir")
    run.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.action in ("list", "check"):
        rows = inventory()
        if args.action == "check" and args.json:
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        else:
            for row in rows:
                status = "代码/配置齐全" if row["runnable_files_present"] else "缺少文件"
                print(f"Level {row['level']}: {row['title']} [{status}]\n  {row['unit']}")
            print("full 入口接入 Level 2C 连续轨迹；原有独立模块仍保留各自协议与历史结果。")
        return 0 if all(r["runnable_files_present"] for r in rows) else 1
    try:
        command = command_for(args.level, args.stage, args.config, args.output_dir)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    print(json.dumps({"cwd": str(PROJECT_ROOT), "argv": command}, ensure_ascii=False), flush=True)
    if args.dry_run:
        return 0
    # 模块从本项目 src 查找；调用解释器与当前环境一致。
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = str(PROJECT_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(command, cwd=PROJECT_ROOT, env=env, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
