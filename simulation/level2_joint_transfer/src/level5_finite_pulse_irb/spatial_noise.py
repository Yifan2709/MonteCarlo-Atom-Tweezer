"""站点间共同/局域空间相关噪声与站点参数不均匀性。

δ_i(t) = sqrt(ρ)·δ_common(t) + sqrt(1-ρ)·δ_local,i(t)，0≤ρ≤1；
站点静态参数（阱深比例、Rabi 比例、静态失谐、η 比例）按配置的
相对标准差抽样。经验相关矩阵由采样直接验证。
"""
from __future__ import annotations

import numpy as np


def correlated_site_noise(rng, n_sites, n_samples, rho, rms):
    """生成 [n_sites, n_samples] 噪声，站点间相关系数 ≈ ρ。"""
    rho = float(np.clip(rho, 0.0, 1.0))
    common = rng.normal(0.0, rms, size=(1, n_samples))
    local = rng.normal(0.0, rms, size=(n_sites, n_samples))
    return np.sqrt(rho) * common + np.sqrt(1.0 - rho) * local


def empirical_cross_site_correlation(traces: np.ndarray) -> np.ndarray:
    """站点间经验相关矩阵（对角置 1）。"""
    traces = np.atleast_2d(np.asarray(traces, dtype=float))
    if traces.shape[0] < 2:
        return np.eye(traces.shape[0])
    corr = np.corrcoef(traces)
    return np.nan_to_num(corr, nan=0.0)


def distance_kernel_correlation(positions_m, length_scale_m):
    """基于站点距离的相关核 exp(-d/ℓ) 矩阵。"""
    pos = np.asarray(positions_m, dtype=float)
    diff = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    return np.exp(-diff / max(length_scale_m, 1e-12))


def sample_site_parameters(n_sites: int, seed: int,
                           trap_depth_relative_std: float = 0.0,
                           rabi_relative_std: float = 0.0,
                           static_detuning_std_hz: float = 0.0,
                           eta_relative_std: float = 0.0,
                           temperature_relative_std: float = 0.0,
                           positions=None, spatial_correlation: float = 0.0,
                           length_scale_m: float = 10e-6):
    """抽样站点参数（返回 dict of arrays，含描述）。

    站点参数只进入量子侧（失谐/η/Rabi 比例；阱深比例按势能线性缩放
    u_SLM/u_AOD 轨迹）。经典轨迹池为冻结重放，不按站点重生成——这是
    frozen-replay 的明确边界，报告与 metrics 中必须注明。
    """
    rng = np.random.default_rng(int(seed))
    n_sites = int(n_sites)
    if positions is None:
        positions = None
    if spatial_correlation > 0 and positions is not None:
        kern = distance_kernel_correlation(positions, length_scale_m)
        kern = spatial_correlation * kern
        np.fill_diagonal(kern, 1.0)
        chol = np.linalg.cholesky(kern + 1e-12 * np.eye(n_sites))
        depth_scale = 1.0 + chol @ rng.normal(0.0, trap_depth_relative_std,
                                              n_sites) \
            if trap_depth_relative_std > 0 else np.ones(n_sites)
    else:
        depth_scale = (1.0 + rng.normal(0.0, trap_depth_relative_std, n_sites)
                       if trap_depth_relative_std > 0 else np.ones(n_sites))
    rabi_scale = (1.0 + rng.normal(0.0, rabi_relative_std, n_sites)
                  if rabi_relative_std > 0 else np.ones(n_sites))
    static_detuning = (rng.normal(0.0, static_detuning_std_hz, n_sites)
                       if static_detuning_std_hz > 0 else np.zeros(n_sites))
    eta_scale = (1.0 + rng.normal(0.0, eta_relative_std, n_sites)
                 if eta_relative_std > 0 else np.ones(n_sites))
    temp_scale = (1.0 + rng.normal(0.0, temperature_relative_std, n_sites)
                  if temperature_relative_std > 0 else np.ones(n_sites))
    return {
        "site_index": np.arange(n_sites),
        "depth_scale": depth_scale,
        "rabi_scale": rabi_scale,
        "static_detuning_hz": static_detuning,
        "eta_scale": eta_scale,
        "temperature_scale": temp_scale,
        "seed": int(seed),
        "description": {
            "trap_depth_relative_std": trap_depth_relative_std,
            "rabi_relative_std": rabi_relative_std,
            "static_detuning_std_hz": static_detuning_std_hz,
            "eta_relative_std": eta_relative_std,
            "spatial_correlation": spatial_correlation,
            "note": "站点参数仅进入量子通道与势能线性缩放；经典轨迹池为"
                    "冻结重放，未按站点重新生成（frozen-replay 边界）"},
    }


def rectangular_site_positions(n_sites: int, pitch_m: float = 9.6e-6,
                               aspect: float | None = None):
    """矩形阵列站点坐标（aspect=n_cols/n_rows 提示，或自动近方形）。"""
    n = int(n_sites)
    n_cols = int(np.ceil(np.sqrt(n * (aspect or 1.0))))
    n_rows = int(np.ceil(n / n_cols))
    xs, ys = np.meshgrid(np.arange(n_cols) * pitch_m,
                         np.arange(n_rows) * pitch_m)
    pos = np.stack([xs.ravel()[:n], ys.ravel()[:n]], axis=1)
    return pos


def summarize_site_parameters(params: dict) -> dict:
    """站点参数分布摘要。"""
    out = {"n_sites": int(len(params["site_index"]))}
    for key in ("depth_scale", "rabi_scale", "eta_scale", "temperature_scale"):
        arr = np.asarray(params[key], dtype=float)
        out[key] = {"mean": float(arr.mean()), "std": float(arr.std()),
                    "min": float(arr.min()), "max": float(arr.max())}
    det = np.asarray(params["static_detuning_hz"], dtype=float)
    out["static_detuning_hz"] = {"mean": float(det.mean()),
                                 "std": float(det.std()),
                                 "min": float(det.min()),
                                 "max": float(det.max())}
    return out
