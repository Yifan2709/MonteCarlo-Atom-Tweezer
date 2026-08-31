"""临界速度估计：bootstrap 分段线性插值与 log-logistic 拟合两种方法。

两条约定：
- 俘获率随速度单调不增是对本问题的物理假设，插值在排序后的
  (log v, capture) 折线上进行；
- 曲线未穿越目标水平（全程 ≥ 目标或全程 < 目标）时返回 None 并给出
  方向标记，由 boundary_covered 汇总判定扫描是否覆盖失效边界。
"""
from __future__ import annotations

import numpy as np


def _interp_crossing(speeds, captures, level):
    """在 (log v, capture) 折线上找 capture == level 的穿越速度；未穿越返回 None。"""
    order = np.argsort(speeds)
    v = np.asarray(speeds, dtype=float)[order]
    p = np.asarray(captures, dtype=float)[order]
    if p.size == 0:
        return None, "no_data"
    if np.all(p >= level):
        return None, "above_scanned"
    if np.all(p < level):
        return None, "below_scanned"
    log_v = np.log10(v)
    # 单调化：capture 随速度只允许下降，取前缀最小以抑制数值抖动
    p = np.minimum.accumulate(p)
    for i in range(p.size - 1):
        if p[i] >= level >= p[i + 1] or p[i] > level >= p[i + 1]:
            if p[i] == p[i + 1]:
                return float(v[i + 1]), "ok"
            frac = (p[i] - level) / (p[i] - p[i + 1])
            return float(10.0 ** (log_v[i] + frac * (log_v[i + 1] - log_v[i]))), "ok"
    return None, "no_crossing"


def bootstrap_interpolation_speeds(speeds, captured_flags, level,
                                   num_bootstrap=2000, seed=0):
    """对每个速度点的 shot 重采样后插值，返回临界速度点估计与 95% 区间。

    captured_flags: 长度与 speeds 相同的列表，每项为布尔数组（CRN 下同长度）。
    """
    rng = np.random.default_rng(seed)
    point = _interp_crossing(speeds, [np.mean(f) for f in captured_flags], level)
    replicates = []
    for _ in range(int(num_bootstrap)):
        rates = []
        for flags in captured_flags:
            idx = rng.integers(0, flags.size, size=flags.size)
            rates.append(float(np.mean(np.asarray(flags, dtype=bool)[idx])))
        crossing, status = _interp_crossing(speeds, rates, level)
        if crossing is not None:
            replicates.append(crossing)
    result = {
        "method": "bootstrap_interpolation",
        "level": float(level),
        "point_estimate": point[0],
        "status": point[1],
        "num_bootstrap": int(num_bootstrap),
        "seed": int(seed),
    }
    if replicates:
        result.update({
            "ci_low": float(np.percentile(replicates, 2.5)),
            "ci_high": float(np.percentile(replicates, 97.5)),
            "valid_replicates": len(replicates),
        })
    return result


def _fit_logistic(speeds, captured_flags):
    """MLE 拟合 p(v)=1/(1+exp(−(a+b·log10 v)))；返回 (a, b) 或 None。"""
    from scipy.optimize import minimize
    x = np.log10(np.concatenate([np.full(f.size, s) for s, f in zip(speeds, captured_flags)]))
    y = np.concatenate([np.asarray(f, dtype=float) for f in captured_flags])
    if np.all(y == y[0]):
        return None  # 全 0 或全 1，logistic 不可辨识

    def nll(theta):
        a, b = theta
        z = a + b * x
        # 数值稳定的 log(1+exp(z))
        log_like = np.sum(y * (-np.logaddexp(0.0, -z))
                          + (1.0 - y) * (-np.logaddexp(0.0, z)))
        return -log_like

    best = None
    for b0 in (-4.0, -1.0, 1.0):  # 多起点避免局部极小
        fit = minimize(nll, np.array([0.0, b0]), method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
        if best is None or fit.fun < best.fun:
            best = fit
    return (float(best.x[0]), float(best.x[1])) if best is not None else None


def logistic_fit_speeds(speeds, captured_flags, level,
                        num_bootstrap=500, seed=0):
    """logistic 拟合（MLE）给出临界速度；bootstrap 给出参数不确定度。"""
    theta = _fit_logistic(speeds, captured_flags)
    result = {
        "method": "logistic_fit",
        "level": float(level),
        "seed": int(seed),
    }
    if theta is None:
        result.update({"point_estimate": None, "status": "not_identifiable"})
        return result
    a, b = theta
    if b >= 0:
        result.update({"point_estimate": None, "status": "non_decreasing_fit"})
        return result
    v_level = 10.0 ** ((np.log(level / (1.0 - level)) - a) / b)
    v_min, v_max = float(np.min(speeds)), float(np.max(speeds))
    status = "ok" if v_min <= v_level <= v_max else "extrapolated"
    result.update({"point_estimate": float(v_level), "status": status,
                   "intercept_a": a, "slope_b": b})
    rng = np.random.default_rng(seed)
    replicates = []
    for _ in range(int(num_bootstrap)):
        resampled = [np.asarray(f, dtype=bool)[rng.integers(0, f.size, size=f.size)]
                     for f in captured_flags]
        theta_b = _fit_logistic(speeds, resampled)
        if theta_b is not None and theta_b[1] < 0:
            replicates.append(10.0 ** ((np.log(level / (1.0 - level)) - theta_b[0]) / theta_b[1]))
    if replicates:
        result.update({
            "ci_low": float(np.percentile(replicates, 2.5)),
            "ci_high": float(np.percentile(replicates, 97.5)),
            "valid_replicates": len(replicates),
        })
    return result


def boundary_covered(capture_rates, level=0.95) -> bool:
    """曲线是否真正穿越目标水平：存在 ≥level 的点且存在 <level 的点。"""
    rates = np.asarray(capture_rates, dtype=float)
    return bool(rates.size > 0 and np.any(rates >= level) and np.any(rates < level))
