"""论文 Fig. 6d 数字化参考曲线的读取与对比工具（Level 2C 阶段 0/5）。"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

SERIES_LABELS = {
    "ml_400us": "0.4 ms (ML-optimized)",
    "manual_600us": "0.6 ms manual",
    "manual_400us": "0.4 ms manual",
    "manual_200us": "0.2 ms manual",
}


def load_reference_curves(csv_path: Path) -> dict:
    """读取数字化 CSV；返回 {series: (n, survival)}，空值剔除。"""
    path = Path(csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"参考曲线不存在：{path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    n_all = np.array([int(r["n_one_way_transfers"]) for r in rows])
    curves = {}
    for key in SERIES_LABELS:
        column = f"survival_{key}"
        values = np.array([float(r[column]) if r[column] else np.nan for r in rows])
        mask = np.isfinite(values)
        curves[key] = (n_all[mask], values[mask])
    return curves


def compare_to_reference(n_sim, survival_sim, reference: tuple, tolerance_band: float):
    """模拟曲线与参考曲线的逐点对比（在参考有数据的公共 n 上）。

    返回最大绝对偏差、落在 ±tolerance_band 内的比例与逐点偏差数组。
    """
    n_ref, s_ref = reference
    n_sim = np.asarray(n_sim, dtype=float)
    survival_sim = np.asarray(survival_sim, dtype=float)
    common = np.isin(n_sim, n_ref)
    if not common.any():
        raise ValueError("模拟与参考曲线无公共 n 点")
    sim_at_ref = np.interp(n_ref, n_sim, survival_sim)
    deviation = sim_at_ref - s_ref
    within = np.abs(deviation) <= tolerance_band
    return {
        "n": n_ref.tolist(),
        "reference": s_ref.tolist(),
        "simulated": sim_at_ref.tolist(),
        "deviation": deviation.tolist(),
        "max_abs_deviation": float(np.max(np.abs(deviation))),
        "fraction_within_band": float(within.mean()),
        "tolerance_band": float(tolerance_band),
        "matches_reference": bool(within.all()),
    }
