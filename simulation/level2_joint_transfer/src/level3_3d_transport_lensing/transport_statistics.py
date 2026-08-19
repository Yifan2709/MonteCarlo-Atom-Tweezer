"""统计工具：复用 Level 1 Wilson 区间与 Level 2 配对 bootstrap。"""
from __future__ import annotations

from level1_transfer_1d.analysis import wilson_interval
from level2_joint_transfer.statistics import (paired_bootstrap_difference,
                                              paired_counts)


def summarize_retention(retained: list, name: str = "") -> dict:
    """留阱率与 Wilson 95% 区间。"""
    successes = int(sum(bool(r) for r in retained))
    total = len(retained)
    low, high = wilson_interval(successes, total)
    return {"name": name, "retained": successes, "shots": total,
            "retention_fraction": successes / total if total else None,
            "wilson_low": low, "wilson_high": high}


__all__ = ["wilson_interval", "paired_bootstrap_difference", "paired_counts",
           "summarize_retention"]
