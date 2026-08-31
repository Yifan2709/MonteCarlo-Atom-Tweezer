"""survival(n) 曲线拟合：论文同款 S_n = p0·p^n·(1 − e^(−1/(bn))) 与参数 bootstrap。

不确定度通过对 shot 重抽样（而非对曲线点加噪）获得：每次重抽 shots 的
逐 shot 存活序列，重算 survival 曲线再拟合——与论文"先 S_n 后 D_n 的
两级 bootstrap"中的第一级对应（D_n 属 IRB 范围，不在 Level 2C）。
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit


def survival_model(n, p0, p, b):
    """S_n = p0 · p^n · (1 − exp(−1/(b·n)))；n≥1。b>0 为每次转移累积的归一化温升参数。"""
    n = np.asarray(n, dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        heat = 1.0 - np.exp(-1.0 / np.maximum(b * n, 1e-300))
    return p0 * np.power(p, n) * heat


def fit_survival_curve(n_values, survival, shot_count=None):
    """拟合 survival(n)；返回参数、协方差与诊断。

    全存活（survival 恒为 1）时模型不可辨识：不做拟合，标记 degenerate，
    改用"零丢失"统计给出 p 的单侧 Wilson 95% 下界（总试验次数 = shots × n_max）。
    """
    n = np.asarray(n_values, dtype=float)
    s = np.asarray(survival, dtype=float)
    mask = np.isfinite(s) & (n >= 1)
    n, s = n[mask], s[mask]
    if s.size < 4:
        raise ValueError("拟合至少需要 4 个有效点")
    if np.all(s >= 1.0 - 1e-15):
        result = {"degenerate_full_survival": True, "p0": 1.0, "p": 1.0, "b": 0.0,
                  "p0_ci": None, "p_ci": None, "b_ci": None, "covariance": None}
        if shot_count is not None:
            from level1_transfer_1d.analysis import wilson_interval
            trials = int(shot_count) * int(n.max())
            low, _ = wilson_interval(trials, trials)
            result["p_lower_wilson95"] = float(low)
            result["total_transfer_trials"] = trials
        return result

    def model(nv, p0, p, b):
        return survival_model(nv, p0, p, b)

    p0_guess = float(s[0])
    p_guess = float(np.exp(np.mean(np.diff(np.log(np.clip(s, 1e-12, None))))))
    b_guess = 0.05
    try:
        popt, pcov = curve_fit(model, n, s, p0=[p0_guess, p_guess, b_guess],
                               bounds=([0.0, 0.0, 1e-9], [1.5, 1.0, 1e6]), maxfev=20000)
    except RuntimeError:
        # 初值退化时退到纯指数起步
        popt, pcov = curve_fit(model, n, s, p0=[1.0, 0.999, 0.1],
                               bounds=([0.0, 0.0, 1e-9], [1.5, 1.0, 1e6]), maxfev=50000)
    perr = np.sqrt(np.diag(pcov))
    return {
        "degenerate_full_survival": False,
        "p0": float(popt[0]), "p": float(popt[1]), "b": float(popt[2]),
        "p0_ci": [float(popt[0] - 1.96 * perr[0]), float(popt[0] + 1.96 * perr[0])],
        "p_ci": [float(popt[1] - 1.96 * perr[1]), float(popt[1] + 1.96 * perr[1])],
        "b_ci": [float(popt[2] - 1.96 * perr[2]), float(popt[2] + 1.96 * perr[2])],
        "covariance": pcov.tolist(),
        "residual_rms": float(np.sqrt(np.mean((s - model(n, *popt)) ** 2))),
    }


def bootstrap_fit(transfers_survived, max_transfers, num_bootstrap, seed):
    """对 shot 存活序列做 bootstrap：重抽 shots → 重算曲线 → 重拟合。

    transfers_survived：逐 shot 存活转移数数组（长度 N）。
    返回 p0/p/b 的 bootstrap 分布分位数区间（95%）。
    """
    ts = np.asarray(transfers_survived, dtype=int)
    n_grid = np.arange(1, max_transfers + 1)
    rng = np.random.default_rng(seed)
    samples = {"p0": [], "p": [], "b": []}
    failures = 0
    for _ in range(int(num_bootstrap)):
        draw = rng.choice(ts, size=ts.size, replace=True)
        curve = np.array([(draw >= k).mean() for k in n_grid])
        try:
            fit = fit_survival_curve(n_grid, curve)
        except (RuntimeError, ValueError):
            failures += 1
            continue
        if fit["degenerate_full_survival"]:
            failures += 1
            continue
        samples["p0"].append(fit["p0"])
        samples["p"].append(fit["p"])
        samples["b"].append(fit["b"])
    out = {"num_bootstrap": int(num_bootstrap), "fit_failures": int(failures),
           "seed": int(seed)}
    for key, values in samples.items():
        if values:
            out[f"{key}_ci95"] = [float(np.percentile(values, 2.5)),
                                  float(np.percentile(values, 97.5))]
            out[f"{key}_median"] = float(np.median(values))
        else:
            out[f"{key}_ci95"] = None
            out[f"{key}_median"] = None
    return out
