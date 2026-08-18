"""候选波形评分：俘获为主、激发与光滑度为辅，保存全部组成项。"""
from __future__ import annotations

import numpy as np


def score_candidate(evaluation_summary, smoothness, optimization_cfg):
    """由单次评估摘要与光滑度指标构造可解释的加权总分。

    组成项：
      capture_rate            严格 E_f<0 俘获率（首要）
      margin_quantile         俘获裕量低分位数（越大越好）
      excitation_term         1 - 已俘获原子平均归一化末态激发
      smooth_penalty          无量纲加速度/jerk 惩罚
    排序按 (captured_count, total) 字典序，避免小分差淹没俘获差。
    """
    cfg = optimization_cfg
    capture_rate = float(evaluation_summary["capture_rate"])
    captured_excitation = evaluation_summary.get("mean_final_excitation_normalized_captured")
    if captured_excitation is None:
        excitation_term = 0.0
    else:
        excitation_term = float(np.clip(1.0 - captured_excitation, 0.0, 1.0))
    margin_quantile = float(evaluation_summary.get("margin_quantile", 0.0))
    smooth_penalty = smoothness["dimensionless_accel"] / 50.0 + smoothness["dimensionless_jerk"] / 500.0
    total = (cfg.capture_weight * capture_rate
             + cfg.excitation_weight * excitation_term
             + cfg.smoothness_weight * (-smooth_penalty))
    return {
        "capture_rate": capture_rate,
        "captured_count": int(evaluation_summary["captured"]),
        "shots": int(evaluation_summary["shots"]),
        "margin_quantile": margin_quantile,
        "excitation_term": excitation_term,
        "smooth_penalty": float(smooth_penalty),
        "smoothness": smoothness,
        "objective_total": float(total),
        "rank_key": (int(evaluation_summary["captured"]), float(total)),
    }
