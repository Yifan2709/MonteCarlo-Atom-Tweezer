"""多站点并行模拟：站点集合上的 ensemble 统计与站点分辨输出。

站点 = 条件独立的单原子动力学 + 相关经典噪声（common/local 结构）；
不得称为多体量子模拟或包含原子间串扰。输出 site-resolved survival、
channel 指标、array average、站点间/站内方差分解、common-mode 贡献、
worst/best decile 与站点参数相关性。
"""
from __future__ import annotations

import numpy as np

from .spatial_noise import (correlated_site_noise,
                            empirical_cross_site_correlation,
                            rectangular_site_positions,
                            sample_site_parameters,
                            summarize_site_parameters)


def build_site_ensemble(cfg_section: dict, n_sites: int, seed: int):
    """按配置构建站点集合（参数 + 坐标）。"""
    positions = rectangular_site_positions(n_sites)
    params = sample_site_parameters(
        n_sites, seed,
        trap_depth_relative_std=float(cfg_section.get(
            "trap_depth_relative_std", 0.0)),
        rabi_relative_std=float(cfg_section.get("rabi_relative_std", 0.0)),
        static_detuning_std_hz=float(cfg_section.get(
            "static_detuning_std_Hz", 0.0)),
        eta_relative_std=float(cfg_section.get("eta_relative_std", 0.0)),
        positions=positions,
        spatial_correlation=float(cfg_section.get(
            "spatial_noise_correlation", 0.0)))
    return {"positions_m": positions, "params": params,
            "summary": summarize_site_parameters(params)}


def site_resolved_summary(per_site_values: np.ndarray, site_params: dict,
                          param_keys=("depth_scale", "rabi_scale",
                                      "static_detuning_hz", "eta_scale")):
    """站点分辨统计：均值、方差分解、decile 与参数相关性。"""
    values = np.asarray(per_site_values, dtype=float)
    out = {
        "n_sites": int(values.shape[-1]) if values.ndim > 1 else 1,
        "array_average": float(np.mean(values)),
        "between_site_std": float(np.std(np.mean(values, axis=0)))
        if values.ndim > 1 else 0.0,
        "within_site_std": float(np.mean(np.std(values, axis=0)))
        if values.ndim > 1 else float(np.std(values)),
        "worst_decile": float(np.quantile(np.mean(values, axis=0), 0.1))
        if values.ndim > 1 else float(np.quantile(values, 0.1)),
        "best_decile": float(np.quantile(np.mean(values, axis=0), 0.9))
        if values.ndim > 1 else float(np.quantile(values, 0.9)),
    }
    cors = {}
    site_means = np.mean(values, axis=0) if values.ndim > 1 else values
    for key in param_keys:
        arr = np.asarray(site_params[key], dtype=float)
        if np.std(arr) > 1e-12 and np.std(site_means) > 1e-12:
            cors[key] = float(np.corrcoef(arr, site_means)[0, 1])
    out["parameter_correlations"] = cors
    return out


def common_local_decomposition(traces_by_site: np.ndarray):
    """把站点噪声/信号分解为 common 与 local 贡献（经验相关法）。"""
    traces = np.asarray(traces_by_site, dtype=float)
    corr = empirical_cross_site_correlation(traces)
    off = corr[np.triu_indices_from(corr, k=1)]
    common_fraction = float(np.mean(off)) if off.size else 0.0
    return {"mean_off_diagonal_correlation": common_fraction,
            "correlation_matrix": corr.tolist(),
            "common_mode_fraction_estimate": max(0.0, common_fraction),
            "local_fraction_estimate": 1.0 - max(0.0, common_fraction)}


def generate_common_local_noise(n_sites, n_samples, rho, rms, seed):
    """生成 common/local 结构的站点噪声并验证经验相关矩阵。"""
    rng = np.random.default_rng(int(seed))
    traces = correlated_site_noise(rng, n_sites, n_samples, rho, rms)
    validation = common_local_decomposition(traces)
    validation["target_rho"] = float(rho)
    validation["rms"] = float(rms)
    return traces, validation
