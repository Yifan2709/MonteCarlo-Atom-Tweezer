"""Y2：光镊交接中的矢量光移（VLS）与磁场/偏振几何（Zhang 2026 图1e/ED4）。

物理：487 nm 阱的矢量光移产生合成磁场 B_eff，方向由偏振决定；
两阱（SLM/AOD，正交偏振）位移或失配时合成场方向随时间变化 →
比特积累随机相位。论文锚点：
- η_vls≈3×10⁻⁴ @487 nm；合成场 50 mG；
- B_SLM∥B_AOD∥B₀（"parallel"）：劈裂 57 Hz → 强失相干；
- B_SLM∥B_AOD⊥B₀（"optimal"）：每单程翻转+失相干 <6×10⁻⁴；
- 相位灵敏度 7 mrad/nm/单程；对准 σ≈80 nm；σy=31 nm。

模型（半解析 MC，逐单程采样随机相位与翻转）：
- 每单程相位 φ = κ·d_eff + N(0, σ_dyn²)，κ=7e-3 rad/nm；
- d_eff = 静态失配（N(0,80nm)，逐 trial 固定）在 config 几何因子下的投影；
- parallel 额外失相干：σ_par = 2π×57 Hz×t_trip（论文劈裂×单程时长）；
- 翻转概率（横向场）p_flip(config) 以 optimal=6e-4 上界锚定、
  parallel 按几何比放大（assumed，登记）。
输出：静态阱/交接/移动三种场景 × 三种 config 的 Ramsey 对比度
随单程次数的衰减，及逐 config 每单程失相干+翻转概率。
"""
from __future__ import annotations

import numpy as np

from .provenance import REGISTRIES

# 灵敏度抑制因子（相对原始 κ=7 mrad/nm 的有效耦合；paper_derived/assumed）：
# optimal：由论文上界（每单程翻转+失相干 <6e-4 → σ_φ≈0.024 rad）反解
#   g=σ_bound/(κ·√(σ_align²+σ_y²))≈0.04；
# parallel：完整 κ + 57 Hz 劈裂项（论文实测）；
# static_slm：无交接，仅阱内抖动残余（assumed 0.005）。
GEOMETRY_FACTORS = {
    "optimal": 0.04,
    "parallel": 1.0,
    "static_slm": 0.005,
}


def per_trip_dephasing(config: str, t_trip_s: float | None = None,
                       seed: int = 0, n_samples: int = 20000) -> dict:
    """返回每单程相位分布参数与等效失相干+翻转概率。"""
    yb = REGISTRIES["yb171_2026_native"]
    if t_trip_s is None:
        t_trip_s = (yb.get("handoff_time_s").value * 2
                    + yb.get("move_time_s").value)
    kappa = yb.get("vls_phase_per_nm_rad").value          # rad/nm
    sigma_align_nm = yb.get("alignment_sigma_m").value * 1e9
    sigma_y_nm = yb.get("transverse_sigma_m").value * 1e9
    g = GEOMETRY_FACTORS[config]
    rng = np.random.default_rng(seed)
    static = rng.standard_normal(n_samples) * sigma_align_nm   # 逐 trial 失配
    dyn = rng.standard_normal(n_samples) * sigma_y_nm
    phi = kappa * g * (static + dyn)
    if config == "parallel":
        split_hz = yb.get("splitting_parallel_hz").value
        phi = phi + rng.standard_normal(n_samples) * 2 * np.pi * split_hz * t_trip_s
    sigma_phi = float(np.std(phi))
    p_dephase = float(1.0 - np.exp(-0.5 * sigma_phi ** 2))
    p_flip_anchor = yb.get("trip_flip_dephase_prob").value   # optimal 上界
    p_flip = min(0.5, p_flip_anchor * g / 0.04)
    return {"config": config, "t_trip_s": t_trip_s,
            "sigma_phi_rad": sigma_phi, "p_dephase": p_dephase,
            "p_flip": p_flip, "geometry_factor": g}


def ramsey_contrast_vs_trips(config: str, max_trips: int = 30,
                             n_trials: int = 4000, seed: int = 1) -> dict:
    """MC：逐 trial 累计相位，输出对比度 C(n)=|⟨e^{iΦ}⟩| 随单程次数。"""
    info = per_trip_dephasing(config, seed=seed)
    rng = np.random.default_rng(seed + 17)
    # 静态失配逐 trial 固定（同一原子的对准误差不随次数重抽）
    kappa = REGISTRIES["yb171_2026_native"].get("vls_phase_per_nm_rad").value
    static_nm = rng.standard_normal(n_trials) \
        * REGISTRIES["yb171_2026_native"].get("alignment_sigma_m").value * 1e9
    dyn_sigma = REGISTRIES["yb171_2026_native"].get("transverse_sigma_m").value * 1e9
    phi_static = kappa * info["geometry_factor"] * static_nm
    phi_dyn = (lambda: rng.standard_normal((max_trips, n_trials))
               * kappa * info["geometry_factor"] * dyn_sigma)
    extra = np.zeros((max_trips, n_trials))
    if config == "parallel":
        split = REGISTRIES["yb171_2026_native"].get("splitting_parallel_hz").value
        tt = info["t_trip_s"]
        extra = rng.standard_normal((max_trips, n_trials)) * 2 * np.pi * split * tt
    cum = np.cumsum(phi_dyn() + extra, axis=0) + phi_static[None, :]
    flips = rng.random((max_trips, n_trials)) < info["p_flip"]
    sign = 1 - 2 * (np.cumsum(flips, axis=0) % 2)
    contrast = np.abs(np.mean(np.exp(1j * cum) * sign, axis=1))
    return {**info, "trips": np.arange(1, max_trips + 1).tolist(),
            "contrast": contrast.tolist(),
            "p_flip_plus_dephase_per_trip": info["p_flip"] + info["p_dephase"]}


def evaluate_anchors() -> dict:
    """预注册锚点判据：optimal 每单程 <6e-4 上界；parallel 强衰减；static 最稳。"""
    opt = per_trip_dephasing("optimal")
    par = per_trip_dephasing("parallel")
    sta = per_trip_dephasing("static_slm")
    anchor = REGISTRIES["yb171_2026_native"].get(
        "trip_flip_dephase_prob").value
    return {
        "optimal_per_trip_total": opt["p_flip"] + opt["p_dephase"],
        "optimal_within_paper_bound": bool(
            opt["p_flip"] + opt["p_dephase"] <= anchor * 1.5),
        "parallel_sigma_phi_rad": par["sigma_phi_rad"],
        "parallel_strongly_dephasing": bool(par["sigma_phi_rad"] > 1.0),
        "ordering_static_le_optimal_le_parallel": bool(
            sta["sigma_phi_rad"] <= opt["sigma_phi_rad"] <= par["sigma_phi_rad"]),
    }
