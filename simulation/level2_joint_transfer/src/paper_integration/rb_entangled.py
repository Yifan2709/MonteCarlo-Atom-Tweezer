"""R2：两原子纠缠态运输与测量（Bluvstein 2022 图1d 的结构复现）。

协议（R 图1）：两原子各自在独立阱中（同一 AOD 组内相邻音，成对 CZ 制备
Bell 态 |Φ+⟩ 后其中一阱把原子移动 110 μm；另一阱静止）。
与 Cs 的 pick-up/放回不同：原子全程在各自阱内，无交接。

内部态模型（Pauli 轨迹 MC，逐 trial 从经典动力学取噪声实现）：
- 理想制备 |Φ+⟩=（|00⟩+|11⟩）/√2（论文允许的理想制备基准口径）；
- 运输激发：移动阱中原子的剩余振荡能 ΔE（R1 动力学）→ 差分光移相位
  δφ = c_φ·ΔN（c_φ 在训练锚点 D=110 μm, T=300 μs 上标定并冻结；
  登记 calibrated_on_training，验证只看形状：长 T 平台、短 T 下降）；
- 静止阱原子受静态阱噪声（小，T2* 模型）；
- XY16 DD（理想回波）把静态失相干压低 T2*/T2 因子（开关对照）；
- 机械丢失：两原子各自的 R1 动力学存活；
- 读出：宇称振荡（分析旋转角 θ 扫描）与 ⟨ZZ⟩；损失读为 |1⟩（dark）
  与双存活后选两种口径分别报告。

结构判据：
C1 T≥300 μs（b≤0.37 μm/μs < 0.55）保真度平台（相对平台降幅 < 2pp）；
C2 存在 T_fast 使保真度显著低于平台（快速运输下降）；
C3 无 DD 时保真度整体更低（回波对照）。
"""
from __future__ import annotations

import numpy as np

from .provenance import REGISTRIES
from .quantum_mc import BatchedState
from .rb_transport import rb_physics, run_single_transport


def calibrate_c_phi(distance_m: float, time_s: float, seed: int = 7,
                    shots: int = 512, target_raw_fidelity: float
                    | None = None) -> dict:
    """在论文锚点（110 μm/300 μs，raw Bell 94.8%）上标定 c_φ 并冻结。

    raw 口径含：SPAM(≈1%)、单程机械损失、失相干。模型：
    F_raw ≈ F_spam · P_surv_pair · (1+exp(−2σ_φ²))/2 → 解 σ_φ，
    c_φ = σ_φ / ΔN(锚点)。
    """
    rb = REGISTRIES["rb87_2022_native"]
    if target_raw_fidelity is None:
        target_raw_fidelity = rb.get("bell_raw_fidelity").value
    f_spam = 0.99                              # assumed：读出+制备 ≈1%
    res = run_single_transport(distance_m, time_s, shots, seed, 0.05e-6)
    p_pair = res["survival"] ** 2              # 两原子都存活（近似：静止者≈1）
    dn = res["mean_delta_N"]
    # F = f_spam·p_pair·(1+e^{−2σ²})/2 = target
    inner = 2 * target_raw_fidelity / (f_spam * max(p_pair, 1e-6)) - 1
    sigma_phi = float(np.sqrt(-0.5 * np.log(max(inner, 1e-12))))
    c_phi = sigma_phi / max(dn, 1e-9)
    return {"c_phi_rad_per_quanta": c_phi, "sigma_phi_anchor": sigma_phi,
            "delta_N_anchor": dn, "p_survival_anchor": res["survival"],
            "f_spam_assumed": f_spam,
            "provenance": "calibrated_on_training(110um/300us, raw 94.8%)"}


def entangled_transport(distance_m: float, time_s: float, shots: int,
                        seed: int, c_phi: float, dd: bool = True,
                        analysis_angles: tuple[float, ...] = (
                            0.0, np.pi / 8, np.pi / 4, 3 * np.pi / 8,
                            np.pi / 2)) -> dict:
    """一次 Bell 运输 MC：返回各口径保真度与宇称振荡对比度。"""
    rb = REGISTRIES["rb87_2022_native"]
    res = run_single_transport(distance_m, time_s, shots, seed, 0.05e-6)
    dn = res["mean_delta_N"]
    if not np.isfinite(dn) or not res["alive"].any():
        dn = 1e3          # 全灭/异常：占据大失相干，F→0
    sigma_phi = abs(c_phi) * max(dn, 0.0)
    t2s = rb.get("t2_star_s").value
    t2 = rb.get("t2_dd_s").value if dd else t2s
    # 静态失相干残余：DD 把 σ_static 压到 T2*/T2 水平（近似模型）
    sigma_static = max(time_s * 2 / t2, 0.0) if time_s > 0 else 0.0
    sigma_total = float(np.hypot(sigma_phi, sigma_static))
    rng = np.random.default_rng(seed + 31)
    # 逐 trial 相位差（高斯）+ 机械丢失（读为 |1⟩）
    dphi = rng.standard_normal(shots) * sigma_total
    lost_moving = ~res["alive"]
    lost_static = rng.random(shots) < 0.004     # 静止阱保持损失（assumed 小）
    state = BatchedState(2, shots, rng)
    state.h(0)
    state.cx(0, 1)
    # 施加逐 trial 相位：Z1 旋转 −dphi/2、Z2 +dphi/2（等效 Φ+ 失相干）
    fz1 = dphi < 0   # 符号化近似：相位连续分布用两段 Pauli 近似不可行 →
    # 改用确定性相位旋转的等价统计：直接对 |Φ+⟩ 计算（见下）
    # —— 用解析读出替代态矢量（|Φ+⟩ + 已知相位差的关联完全解析）：
    # P(00)=P(11)=(1+cos dphi)/4? 对 Φ+ 在 ZZ 读出：P(同)= (1+cos dphi)/2
    cosphi = np.cos(dphi)
    # ZZ 关联（读出含 FP/FN=0 & 损失读|1⟩：lost→bit=1）
    def zz_outcome(dark_read_as_one: bool):
        b0 = np.where(lost_static, 1, np.zeros(shots, int))
        b1 = np.where(lost_moving, 1, np.zeros(shots, int))
        # 未丢失部分：Z 读出关联采样
        u = rng.random(shots)
        same = u < (1 + cosphi) / 2
        z0 = np.where(lost_static, b0, np.where(same, 0, 0))
        z1 = np.where(lost_moving, b1, np.where(same, 0, 1))
        if dark_read_as_one:
            return z0, z1
        # 后选口径标记
        keep = ~lost_static & ~lost_moving
        return z0, z1, keep

    z0, z1 = zz_outcome(True)
    zz = 1 - 2 * np.mean(z0 ^ z1)
    # 宇称振荡：分析旋转 θ 后 X 基关联 ⟨X(θ)X(−θ)⟩ = cos(dphi − 2θ) 期望
    contrasts = []
    for theta in analysis_angles:
        c = float(np.mean(np.cos(dphi - 2 * theta)))
        contrasts.append(c)
    c_max = float(max(contrasts))
    fid_raw = float(0.5 * (zz + c_max)) * 0.99   # 含 SPAM ≈1%（assumed）
    _, _, keep = zz_outcome(False)
    fid_ps = float(0.5 * ((1 - 2 * np.mean(z0[keep] ^ z1[keep]))
                          + float(np.mean(np.cos(dphi[keep]))))) \
        if keep.any() else float("nan")
    return {
        "distance_m": distance_m, "time_s": time_s, "shots": shots,
        "seed": seed, "dd": dd,
        "mean_delta_N": dn, "sigma_phi_rad": sigma_total,
        "zz_correlation": float(zz), "parity_contrast": c_max,
        "fidelity_raw_loss_as_one": fid_raw,
        "fidelity_postselected_both_alive": fid_ps,
        "both_alive_frac": float(np.mean(keep)),
        "paper_raw_fidelity": REGISTRIES["rb87_2022_native"].get(
            "bell_raw_fidelity").value,
        "paper_flat_speed_m_per_s": REGISTRIES["rb87_2022_native"].get(
            "bell_flat_speed_m_per_s").value,
    }


def scan_entangled(distance_m: float, times_s: list[float], shots: int,
                   seed: int, c_phi: float, dd: bool = True) -> list[dict]:
    return [entangled_transport(distance_m, t, shots, seed + 13 * i, c_phi, dd)
            for i, t in enumerate(times_s)]


def evaluate_structure(scan: list[dict]) -> dict:
    fids = [r["fidelity_raw_loss_as_one"] for r in scan]
    ts = [r["time_s"] for r in scan]
    idx_slow = [i for i, t in enumerate(ts) if t >= 300e-6]
    platform = float(np.mean([fids[i] for i in idx_slow])) if idx_slow else float("nan")
    spread = float(np.std([fids[i] for i in idx_slow])) if idx_slow else float("nan")
    c1 = bool(spread < 0.02)
    c2 = bool(min(fids) < platform - 0.05)
    return {"C1_flat_platform": c1, "platform": platform, "platform_spread": spread,
            "C2_fast_transport_drop": c2, "min_fid": float(min(fids)),
            "C3_note": "DD 对照另行运行（entangled_transport(dd=False)）"}
