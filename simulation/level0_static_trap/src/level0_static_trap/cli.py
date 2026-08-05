"""命令行入口。

用法：

    python -m level0_static_trap.cli \
        --config configs/default.yaml \
        --trap both \
        --output-dir outputs/level0_static_trap/demo

退出码 0 表示成功；非 0 表示配置或运行错误。
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

from . import __init__  # noqa: F401  触发 Agg 后端设置
from .config import ConfigError, load_config
from .simulation import (
    BoundTrajectoryError,
    build_summary_md,
    run_single_trap,
    write_outputs,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="level0-static-trap",
        description="Level 0：一维径向静态高斯光偶极阱自洽检查。",
    )
    p.add_argument("--config", required=True, help="YAML 配置文件路径")
    p.add_argument(
        "--trap",
        required=True,
        choices=["slm", "aod", "both"],
        help="选择运行的势阱：slm / aod / both",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="输出目录（覆盖配置中 output.directory，并写入 config_used.yaml）",
    )
    return p


def run(args: argparse.Namespace) -> int:
    try:
        cfg = load_config(args.config, output_dir_override=args.output_dir)
    except ConfigError as exc:
        print(f"[配置错误] {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"[文件缺失] {exc}", file=sys.stderr)
        return 2

    # 决定要运行哪些势阱
    if args.trap == "both":
        traps = [("slm", cfg.slm), ("aod", cfg.aod)]
    elif args.trap == "slm":
        traps = [("slm", cfg.slm)]
    else:
        traps = [("aod", cfg.aod)]

    results = {}
    try:
        for key, trap in traps:
            results[key] = run_single_trap(trap, cfg)
    except BoundTrajectoryError as exc:
        print(f"[束缚条件错误] {exc}", file=sys.stderr)
        return 3
    except ValueError as exc:
        print(f"[运行错误] {exc}", file=sys.stderr)
        return 4

    # 写出
    try:
        image_results = write_outputs(cfg, args.trap, results)
    except Exception as exc:  # pragma: no cover - 输出错误
        print(f"[输出错误] {exc}", file=sys.stderr)
        traceback.print_exc()
        return 5

    # 汇总 Markdown
    if cfg.save_md:
        from .io_utils import write_markdown

        md = build_summary_md(cfg, args.trap, results, image_results)
        write_markdown(os.path.join(cfg.output_dir, "summary.md"), md)

    print(f"[完成] 输出目录：{cfg.output_dir}")
    print(f"[完成] 图像检查：{sum(1 for r in image_results.values() if r['checks'].get('overall_pass'))}"
          f"/{len(image_results)} 张通过")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
