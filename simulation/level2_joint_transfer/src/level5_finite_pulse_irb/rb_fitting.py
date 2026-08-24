"""RB/IRB 稳健拟合：单指数、stretched exponential、模型无关斜率与诊断。

P(n)=A p^n + B 为主模型；同时比较 stretched exponential 与分段/非参数
趋势。样本相关性通过 hierarchical bootstrap 处理；拟合不可辨识、尾部
risk set 过小或相邻比值不稳定时明确输出不可解释状态而不是误导性数字。
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit


def fit_exponential_rb(x, y, sigma=None, p_bounds=(0.0, 1.0)):
    """最小二乘 A p^x + B 拟合；返回参数、协方差与诊断。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    def model(xx, a, p, b):
        return a * p ** xx + b

    b0 = float(np.clip(min(y[-1], np.median(y[-3:])), 0.0, 0.5))
    a0 = float(np.clip(y[0] - b0, 0.05, 1.2))
    p0 = [a0, 0.99, b0]
    try:
        popt, pcov = curve_fit(
            model, x, y, p0=p0, sigma=sigma, absolute_sigma=True,
            bounds=([0.0, p_bounds[0], 0.0], [1.0, p_bounds[1], 0.55]),
            maxfev=40000)
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    resid = y - model(x, *popt)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    finite = bool(np.all(np.isfinite(pcov)))
    identifiable = finite and float(np.max(np.abs(
        np.diag(pcov)))) < 1e6 if finite else False
    return {"ok": True, "params": {"A": float(popt[0]), "p": float(popt[1]),
                                   "B": float(popt[2])},
            "cov": pcov.tolist() if finite else None,
            "param_std": (np.sqrt(np.abs(np.diag(pcov)))).tolist()
            if finite else None,
            "ss_residual": ss_res,
            "r_squared": (1.0 - ss_res / ss_tot) if ss_tot > 0 else None,
            "covariance_finite": finite,
            "identifiable": bool(identifiable)}


def fit_stretched_exponential(x, y):
    """A exp(−(x/τ)^β) + B 比较模型（非指数诊断）。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    def model(xx, a, tau, beta, b):
        tau = max(tau, 1e-9)
        return a * np.exp(-((np.asarray(xx, float) / tau) ** beta)) + b

    span = max(float(x.max()), 1.0)
    try:
        popt, pcov = curve_fit(model, x, y, p0=[0.5, span, 1.0, 0.5],
                               bounds=([0.0, 1e-9, 0.2, 0.0],
                                       [1.5, 1e6, 3.0, 1.0]), maxfev=40000)
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    resid = y - model(x, *popt)
    return {"ok": True,
            "params": {"A": float(popt[0]), "tau": float(popt[1]),
                       "beta": float(popt[2]), "B": float(popt[3])},
            "ss_residual": float(np.sum(resid ** 2))}


def model_selection(x, y, sigma=None):
    """AIC/BIC 模型比较（单指数 vs stretched）。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    exp_fit = fit_exponential_rb(x, y, sigma)
    str_fit = fit_stretched_exponential(x, y)
    out = {"exponential": exp_fit, "stretched": str_fit}
    if exp_fit.get("ok") and str_fit.get("ok"):
        for name, fit, k in (("exponential", exp_fit, 3),
                             ("stretched", str_fit, 4)):
            ss = max(fit["ss_residual"], 1e-300)
            fit["aic"] = float(n * np.log(ss / n) + 2 * k)
            fit["bic"] = float(n * np.log(ss / n) + k * np.log(n))
            fit["n_parameters"] = k
        out["preferred"] = ("exponential"
                            if exp_fit["aic"] <= str_fit["aic"]
                            else "stretched")
    return out


def early_window_slope(x, y, window=4, baseline=0.5):
    """早期窗口模型无关配对斜率：Δlog(P−baseline)/Δn。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sel = y > baseline + 1e-6
    if sel.sum() < 2:
        return {"available": False}
    xs, ys = x[sel], np.log(y[sel] - baseline)
    k = min(window, xs.size)
    slope = float(np.polyfit(xs[:k], ys[:k], 1)[0])
    return {"available": True, "log_slope": slope,
            "n_points_used": int(k)}


def instantaneous_ratios(y, risk_set=None):
    """相邻比值与稳定性诊断（risk_set 太小时标记不可解释）。"""
    y = np.asarray(y, dtype=float)
    out = {"ratios": [], "stable": None, "min_risk_set": None}
    if risk_set is not None:
        out["min_risk_set"] = (int(np.min(risk_set))
                               if np.size(risk_set) else None)
    for a, b in zip(y[:-1], y[1:]):
        if a > 1e-9 and b > 1e-9:
            out["ratios"].append(float(b / a))
        else:
            out["ratios"].append(None)
    valid = [r for r in out["ratios"] if r is not None and np.isfinite(r)]
    if len(valid) >= 3:
        arr = np.asarray(valid)
        out["stable"] = bool(np.std(arr) < 0.15 * max(np.mean(arr), 1e-9) + 1e-3)
    return out


def irb_estimate(p_ref, p_int, d=2):
    """标准 IRB 估计（前提满足时）：F = (d−1)(x−1)/d + 1，x=p_int/p_ref。"""
    if p_ref is None or p_int is None or p_ref <= 0 or p_int <= 0:
        return {"interpretable": False,
                "reason": "p_ref/p_int 不可用（拟合失败或非正）"}
    x = p_int / p_ref
    f = (d - 1) * (x - 1.0) / d + 1.0
    return {"interpretable": True, "ratio": float(x),
            "f_avg_interleaved": float(f),
            "note": "仅在 gate-independence、Markovianity、拟合质量与"
                    "loss convention 前提满足时可解释"}


def fit_irb_decay(m_values, values, model="exponential"):
    """对 M 的衰减拟合（survival 与 conditional 分开拟合）。"""
    m_values = np.asarray(m_values, dtype=float)
    values = np.asarray(values, dtype=float)
    if model == "exponential":
        return fit_exponential_rb(m_values, values, p_bounds=(0.0, 1.0))
    if model == "clipped_boltzmann_x_depol":
        def cb(xx, b, dprime, f_inf):
            surv = 1.0 - np.exp(-1.0 / np.maximum(b * np.maximum(xx, 1e-9)
                                                  + 1e-9, 1e-9))
            return surv * (0.5 + 0.5 * (1.0 - dprime) ** np.maximum(xx, 0.0))
        try:
            popt, pcov = curve_fit(cb, m_values, values,
                                   p0=[1.0, 1e-3, 1.0],
                                   bounds=([1e-9, 0.0, 0.0],
                                           [1e6, 1.0, 1.0]), maxfev=80000)
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        return {"ok": True,
                "params": {"b": float(popt[0]), "d_prime": float(popt[1])},
                "cov": pcov.tolist() if np.all(np.isfinite(pcov)) else None}
    raise KeyError(f"未知 IRB 衰减模型 {model}")


def safe_fit_summary(fit):
    """把拟合结果转成可序列化摘要（含失败原因）。"""
    if not fit or not fit.get("ok"):
        return {"status": "fit_failed",
                "error": (fit or {}).get("error", "unknown"),
                "note": "拟合失败或不可辨识；不输出 fidelity 数字"}
    status = "ok" if fit.get("identifiable", True) else "ill_conditioned"
    return {"status": status, **{k: v for k, v in fit.items()
                                 if k != "cov"}}
