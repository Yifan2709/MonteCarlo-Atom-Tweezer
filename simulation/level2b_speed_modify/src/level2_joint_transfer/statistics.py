"""Wilson 区间、配对 bootstrap 与配对计数统计。"""
from __future__ import annotations

import numpy as np

from level1_transfer_1d.analysis import wilson_interval  # 复用 Level 1 实现


def paired_counts(captured_a, captured_b):
    """同一初态下两波形的四种配对结局计数。"""
    a = np.asarray(captured_a, dtype=bool)
    b = np.asarray(captured_b, dtype=bool)
    if a.shape != b.shape:
        raise ValueError("配对统计要求两个波形使用完全相同的初态数组")
    return {
        "both_captured": int(np.sum(a & b)),
        "only_a_captured": int(np.sum(a & ~b)),
        "only_b_captured": int(np.sum(~a & b)),
        "both_failed": int(np.sum(~a & ~b)),
        "shots": int(a.size),
    }


def paired_bootstrap_difference(captured_a, captured_b, num_bootstrap, seed):
    """配对 bootstrap：对同一初态重采样，返回俘获率差的点估计与 95% 区间。"""
    a = np.asarray(captured_a, dtype=bool)
    b = np.asarray(captured_b, dtype=bool)
    if a.shape != b.shape:
        raise ValueError("配对 bootstrap 要求相同初态数组")
    point_estimate = float(a.mean() - b.mean())
    rng = np.random.default_rng(seed)
    n = a.size
    indices = rng.integers(0, n, size=(int(num_bootstrap), n))
    diffs = a[indices].mean(axis=1) - b[indices].mean(axis=1)
    return {
        "difference": point_estimate,
        "ci_low": float(np.percentile(diffs, 2.5)),
        "ci_high": float(np.percentile(diffs, 97.5)),
        "num_bootstrap": int(num_bootstrap),
        "seed": int(seed),
        "ci_excludes_zero": bool(np.percentile(diffs, 2.5) > 0.0 or np.percentile(diffs, 97.5) < 0.0),
    }


def summarize_capture(captured, name=""):
    """俘获计数与 Wilson 95% 区间摘要。"""
    array = np.asarray(captured, dtype=bool)
    count = int(array.sum())
    low, high = wilson_interval(count, array.size)
    return {
        "name": name, "shots": int(array.size), "captured": count,
        "capture_fraction": count / array.size,
        "wilson_low": float(low), "wilson_high": float(high),
    }
