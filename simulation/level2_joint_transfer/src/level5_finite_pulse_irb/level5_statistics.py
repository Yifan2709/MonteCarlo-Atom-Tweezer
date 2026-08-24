"""统计设计：种子空间隔离、hierarchical bootstrap 与区间报告。

seed 空间至少分离 unit_test / noise_calibration / channel_validation /
reference_rb / interleaved_rb / site_validation / sensitivity 七个池；
校准与 validation observables 不同。hierarchical bootstrap 依次重采样
sequence → site → classical/noise realization，不得把相关样本当独立
Bernoulli。二项计数（逐 shot 0/1 survival/loss）另报 Wilson 95% 区间，
仅作有限样本描述性统计。
"""
from __future__ import annotations

import numpy as np

SEED_POOLS = {
    "unit_test_seed_space": 51000,
    "noise_calibration_pool": 51100,
    "channel_validation_pool": 51200,
    "reference_rb_pool": 51300,
    "interleaved_rb_pool": 51400,
    "site_validation_pool": 51500,
    "sensitivity_pool": 51600,
}


def seed_for(pool: str, offset: int = 0) -> int:
    """按池名取隔离种子（offset 内部派生）。"""
    if pool not in SEED_POOLS:
        raise KeyError(f"未知种子池 {pool}；可用：{list(SEED_POOLS)}")
    return int(SEED_POOLS[pool] + int(offset))


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054):
    """Wilson 95% 二项区间（描述性统计）。"""
    from level1_transfer_1d.analysis import wilson_interval as _wilson
    return _wilson(int(successes), int(n))


def hierarchical_bootstrap(values: np.ndarray, n_seq: int, n_site: int,
                           n_noise: int, n_boot: int = 2000, seed: int = 0,
                           statistic=np.mean):
    """三层 bootstrap：sequence → site → noise realization。

    values 形状 [n_seq, n_site, n_noise]（或可 reshape 成该形状的一维）。
    """
    values = np.asarray(values, dtype=float).reshape(n_seq, n_site, n_noise)
    rng = np.random.default_rng(int(seed))
    stats = np.empty(int(n_boot))
    for b in range(int(n_boot)):
        seq_idx = rng.integers(0, n_seq, size=n_seq)
        vals = values[seq_idx]              # [n_seq, n_site, n_noise]
        site_idx = rng.integers(0, n_site, size=(n_seq, n_site))
        vals = np.take_along_axis(
            vals, site_idx[..., None].repeat(n_noise, axis=-1), axis=1)
        noise_idx = rng.integers(0, n_noise, size=(n_seq, n_site, n_noise))
        vals = np.take_along_axis(vals, noise_idx, axis=-1)
        stats[b] = statistic(vals.ravel())
    return {"mean": float(statistic(values.ravel())),
            "ci_low": float(np.percentile(stats, 2.5)),
            "ci_high": float(np.percentile(stats, 97.5)),
            "n_bootstrap": int(n_boot), "seed": int(seed),
            "levels": ["sequence", "site", "noise_realization"]}


def independent_bernoulli_vs_hierarchical_demo(n_seq=20, n_site=5,
                                               n_noise=8, true_p=0.9,
                                               n_boot=2000, seed=7):
    """对比把相关样本当独立 Bernoulli 与层级 bootstrap 的区间宽度。"""
    rng = np.random.default_rng(seed)
    seq_offsets = rng.normal(0, 0.03, n_seq)
    site_offsets = rng.normal(0, 0.02, (n_seq, n_site))
    p = np.clip(true_p + seq_offsets[:, None, None]
                + site_offsets[:, :, None], 0.0, 1.0)
    draws = rng.binomial(1, np.broadcast_to(p, (n_seq, n_site, n_noise)))
    hier = hierarchical_bootstrap(draws, n_seq, n_site, n_noise,
                                   n_boot=n_boot, seed=seed)
    flat = draws.ravel()
    k = int(flat.sum())
    lo, hi = wilson_interval(k, flat.size)
    return {"hierarchical": hier, "naive_wilson": {"low": float(lo),
                                                   "high": float(hi)},
            "n_samples": int(flat.size),
            "note": "Wilson 把全部样本当独立 Bernoulli 会低估不确定性"}
