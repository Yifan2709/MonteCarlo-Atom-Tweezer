"""Level 4 统计：Wilson 区间、配对 bootstrap、生存曲线拟合与加热斜率。

复用 Level 1/2 的 wilson_interval 与 paired bootstrap；新增经验生存曲线
的 clipped-Boltzmann / exponential 拟合（含 AIC/BIC 与不可辨识保护）。
"""
from __future__ import annotations

import numpy as np

from level1_transfer_1d.analysis import wilson_interval
from level2_joint_transfer.statistics import (paired_bootstrap_difference,
                                              paired_counts)

__all__ = ["wilson_interval", "paired_bootstrap_difference", "paired_counts",
           "summarize_success", "survival_fits", "heating_slope"]


def summarize_success(flags, name="") -> dict:
    """二项成功摘要与 Wilson 95% 区间。"""
    array = np.asarray(flags, dtype=bool)
    count = int(array.sum())
    low, high = wilson_interval(count, array.size)
    return {"name": name, "successes": count, "shots": int(array.size),
            "probability": count / array.size if array.size else None,
            "wilson_low": float(low), "wilson_high": float(high)}


def evaluate_outcomes(checkpoint_matches, final_retained,
                      checkpoint_names=None):
    """checkpoint 匹配序列 → first-failure/吸收成功/recaptured 状态机。

    checkpoint_matches 为按时间排序的 bool 数组列表（每个 checkpoint 一项）；
    final_retained 为最终 E_SLM<0 标志。返回逐 shot 的：
    all_checkpoints_ok、roundtrip_success、final_retained、recaptured、
    first_failure_stage、first_failure_checkpoint。
    """
    from .protocol_segments import CHECKPOINT_STAGES
    n_shots = len(np.asarray(checkpoint_matches[0]))
    n_ck = len(checkpoint_matches)
    names = checkpoint_names or [f"checkpoint_{k}" for k in range(n_ck)]
    all_ok = np.ones(n_shots, dtype=bool)
    first_stage = np.full(n_shots, "none", dtype=object)
    first_ck = np.full(n_shots, None, dtype=object)
    for k, match in enumerate(checkpoint_matches):
        match = np.asarray(match, dtype=bool)
        newly_failed = all_ok & ~match
        idx = np.where(newly_failed)[0]
        if idx.size:
            first_stage[idx] = CHECKPOINT_STAGES.get(names[k], names[k])
            first_ck[idx] = names[k]
        all_ok &= match
    final_retained = np.asarray(final_retained, dtype=bool)
    roundtrip_success = all_ok & final_retained
    recaptured = (~roundtrip_success) & final_retained
    return {
        "all_checkpoints_ok": all_ok,
        "roundtrip_success": roundtrip_success,
        "final_retained": final_retained,
        "recaptured": recaptured,
        "absorbing_failure": ~roundtrip_success,
        "first_failure_stage": first_stage,
        "first_failure_checkpoint": first_ck,
    }


def _fit_survival_model(model, rounds, survival, p0):
    """最小二乘拟合生存模型（模型参数→曲线），失败返回 None。"""
    from scipy.optimize import curve_fit
    try:
        popt, pcov = curve_fit(model, rounds, survival, p0=p0, maxfev=20000)
    except Exception:
        return None
    residual = np.asarray(survival) - model(np.asarray(rounds), *popt)
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((np.asarray(survival) - np.mean(survival)) ** 2))
    if not np.all(np.isfinite(pcov)) or not np.all(np.isfinite(popt)):
        return None
    return {"params": [float(v) for v in popt],
            "param_std": [float(np.sqrt(abs(pcov[i, i])))
                          for i in range(len(popt))],
            "ss_residual": ss_res,
            "r_squared": 1.0 - ss_res / ss_tot if ss_tot > 0 else None,
            "cov_finite": bool(np.all(np.isfinite(pcov)))}


def survival_fits(rounds, survival, identifiable_tol=1e-6):
    """拟合经验生存曲线：exponential、clipped-Boltzmann 与组合模型。

    数据几乎全 0/全 1、参数不可辨识或协方差奇异时跳过并说明。
    """
    rounds = np.asarray(rounds, dtype=float)
    survival = np.asarray(survival, dtype=float)
    result = {"n_points": int(len(rounds)), "skipped": {}, "fits": {},
              "nonparametric": survival.tolist()}
    values = np.unique(survival)
    if len(values) < 2 or np.isclose(values.max() - values.min(), 0.0,
                                     atol=identifiable_tol):
        result["skipped"]["reason"] = (
            "生存曲线几乎恒定（全成功或全失败），参数不可辨识，跳过拟合")
        return result
    span = float(np.max(rounds))
    if span < 2:
        result["skipped"]["reason"] = "非平凡轮次不足，跳过拟合"
        return result

    def exponential(n, p):
        return p ** n

    def clipped_boltzmann(n, tau, f_inf):
        return f_inf + (1.0 - f_inf) * np.exp(-n / max(tau, 1e-9))

    def exp_times_clip(n, p, tau, f_inf):
        return f_inf + (1.0 - f_inf) * (p ** n) * np.exp(-n / max(tau, 1e-9))

    n_param = {"exponential": 1, "clipped_boltzmann": 2, "exponential_times_clipped_boltzmann": 3}
    models = {
        "exponential": (exponential, [0.99]),
        "clipped_boltzmann": (clipped_boltzmann, [span, survival[-1]]),
        "exponential_times_clipped_boltzmann": (exp_times_clip,
                                                [0.99, span, survival[-1]]),
    }
    n = len(rounds)
    for name, (model, p0) in models.items():
        fit = _fit_survival_model(model, rounds, survival, p0)
        if fit is None:
            result["skipped"][name] = "拟合失败或协方差奇异"
            continue
        k = n_param[name]
        ss = fit["ss_residual"] if fit["ss_residual"] > 0 else 1e-300
        aic = n * np.log(ss / n) + 2 * k
        bic = n * np.log(ss / n) + k * np.log(n)
        fit["aic"] = float(aic)
        fit["bic"] = float(bic)
        fit["n_parameters"] = k
        result["fits"][name] = fit
    result["note"] = ("拟合参数为经验曲线参数，不得解释为实验温度或单次量子"
                      "通道 fidelity；小样本尾部 S_{n+1}/S_n 波动不代表瞬时 fidelity")
    return result


def heating_slope(rounds, excitation, residual_rtol=0.35):
    """激发 vs 轮数的线性斜率（仅在残差合理时报告）。"""
    rounds = np.asarray(rounds, dtype=float)
    values = np.asarray(excitation, dtype=float)
    if len(rounds) < 3 or np.allclose(values, values[0]):
        return {"slope_uK_per_round": None, "reason": "点数不足或激发恒定"}
    slope, intercept = np.polyfit(rounds, values, 1)
    pred = slope * rounds + intercept
    ss_res = float(np.sum((values - pred) ** 2))
    ss_tot = float(np.sum((values - np.mean(values)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    ok = r2 is None or (1.0 - r2) < residual_rtol * 3
    return {"slope_uK_per_round": float(slope) if ok else None,
            "intercept_uK": float(intercept),
            "r_squared": None if r2 is None else float(r2),
            "linear_fit_reasonable": bool(ok),
            "note": None if ok else "线性拟合残差偏大，不报告线性加热率"}
