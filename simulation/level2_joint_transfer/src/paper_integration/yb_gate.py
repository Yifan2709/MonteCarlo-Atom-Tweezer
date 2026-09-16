"""Y4：Yb 两比特门（AR CZ / TO CZ）与鲁棒性比较（Zhang 2026 图1c/ED6）。

模型（每 trial）：
- 比特位于门区内半径 r（均匀盘，束腰 w_b）；强度误差 ΔI/I = 1−exp(−2r²/w_b²)
  + 指向噪声 N(0, σ_pt)；
- ϵ_TO = ϵ0_TO + c_TO·(ΔI/I)^2；ϵ_AR = ϵ0_AR + c_AR·(ΔI/I)^4
  （指数为论文 ED6c 实测；ϵ0 为论文实测门误差；c 为 assumed，取值使
  代表性边缘 ΔI/I≈10% 时 AR 的强度项贡献 ≈ ϵ0_AR 的 1/3——登记于合同）；
- 误差构成：AR 38% 可擦除检出（Rydberg 衰变→基态/³P₀ 分支）、其余 Pauli；
  TO 75% 为原子丢失、其余 Pauli；
- Rydberg 自发衰变：τ_Ryd=88 μs，门时 T=20.4/Ω_r，衰变概率并入 ϵ0
  （对强度项不重复计）。

输出：ϵ vs ΔI/I 扫描（两种门）、拟合标度指数、平均误差构成、
擦除比例与丢失比例。结构判据：拟合指数 TO≈2、AR≈4；AR 擦除占比 ≈38%；
TO 丢失占比 ≈75%。
"""
from __future__ import annotations

import numpy as np

from .provenance import REGISTRIES

# 强度响应系数（assumed_sensitivity；登记于合同）
C_TO = 0.30   # ΔI/I=10% 时 TO 强度项 ≈ 0.003（≈ ϵ0_TO 的 27%）
C_AR = 300.0  # ΔI/I=10% 时 AR 强度项 ≈ 0.003（≈ ϵ0_AR 的 19%）


def gate_error_prob(gate: str, delta_I_over_I: np.ndarray) -> np.ndarray:
    yb = REGISTRIES["yb171_2026_native"]
    if gate == "AR":
        e0 = yb.get("gate_ar_error").value
        return e0 + C_AR * delta_I_over_I ** 4
    if gate == "TO":
        e0 = yb.get("gate_to_error").value
        return e0 + C_TO * delta_I_over_I ** 2
    raise ValueError(gate)


def run_zone_mc(n_trials: int, seed: int, zone_radius_factor: float = 1.0,
                beam_waist_m: float = 12e-6, pointing_sigma: float = 0.005,
                ) -> dict:
    """门区内均匀分布比特的 CZ 误差 MC，输出两种门的统计与构成。"""
    rng = np.random.default_rng(seed)
    r = np.sqrt(rng.random(n_trials)) * zone_radius_factor * beam_waist_m
    dI = (1 - np.exp(-2 * (r / beam_waist_m) ** 2)
          + rng.standard_normal(n_trials) * pointing_sigma)
    dI = np.clip(dI, 0.0, None)
    out = {"n_trials": n_trials, "mean_delta_I": float(dI.mean())}
    yb = REGISTRIES["yb171_2026_native"]
    for gate in ("AR", "TO"):
        eps = gate_error_prob(gate, dI)
        err = rng.random(n_trials) < eps
        if gate == "AR":
            frac_erasure = yb.get("gate_ar_erasure_frac").value
            frac_loss = 0.10
        else:
            frac_erasure = 0.10
            frac_loss = yb.get("gate_to_loss_frac").value
        u = rng.random(n_trials)
        kind = np.where(~err, "none",
                np.where(u < frac_erasure, "erasure",
                np.where(u < frac_erasure + frac_loss, "loss", "pauli")))
        out[gate] = {
            "mean_error": float(err.mean()),
            "wilson95_halfwidth": _wilson(err.mean(), n_trials),
            "erasure_frac_of_errors": float(
                (kind[kind != "none"] == "erasure").mean()) if err.any() else float("nan"),
            "loss_frac_of_errors": float(
                (kind[kind != "none"] == "loss").mean()) if err.any() else float("nan"),
            "pauli_frac_of_errors": float(
                (kind[kind != "none"] == "pauli").mean()) if err.any() else float("nan"),
        }
    return out


def scan_scaling(gate: str, dI_grid: np.ndarray | None = None,
                 n_trials: int = 200000, seed: int = 3) -> dict:
    """强制 ΔI/I 扫描并拟合 log-log 斜率（应恢复 2 / 4）。

    拟合对象为超出本底 ϵ0 的过剩误差（与论文响应曲线的拟合口径一致），
    仅用过剩 > 5%·ϵ0 的格点。
    """
    if dI_grid is None:
        dI_grid = np.linspace(0.02, 0.20, 10)
    yb = REGISTRIES["yb171_2026_native"]
    eps0 = yb.get("gate_ar_error").value if gate == "AR" \
        else yb.get("gate_to_error").value
    eps_mean = []
    for dI in dI_grid:
        eps = gate_error_prob(gate, np.full(n_trials, dI))
        eps_mean.append(float(eps.mean()))
    eps_arr = np.array(eps_mean) - eps0
    mask = eps_arr > 0.05 * eps0
    if mask.sum() < 3:
        raise ValueError("标度拟合网格不足")
    slope = float(np.polyfit(np.log(dI_grid[mask]), np.log(eps_arr[mask]), 1)[0])
    return {"gate": gate, "dI_grid": dI_grid.tolist(),
            "eps_mean": eps_mean, "fitted_slope": slope}


def rydberg_decay_fraction(gate_time_s: float | None = None) -> float:
    """门期间 Rydberg 自发衰变概率 1−exp(−T/τ_Ryd)。"""
    yb = REGISTRIES["yb171_2026_native"]
    if gate_time_s is None:
        gate_time_s = 20.4 / yb.get("rydberg_omega_hz").value
    tau = yb.get("rydberg_lifetime_s").value
    return float(1 - np.exp(-gate_time_s / tau))


def evaluate_anchors(zone_result: dict, scaling: dict) -> dict:
    yb = REGISTRIES["yb171_2026_native"]
    exp_ar = yb.get("gate_ar_scaling_exp").value
    exp_to = yb.get("gate_to_scaling_exp").value
    fr_ar = yb.get("gate_ar_erasure_frac").value
    fr_to = yb.get("gate_to_loss_frac").value
    return {
        "model_nature": ("given_error_model_sampling：ϵ=ϵ0+c·(ΔI/I)^k 为"
                         "登记的假设有效通道；标度指数拟合是对预置公式的"
                         "自洽检查，非门动力学预测（验证 D5）"),
        "scaling_AR": scaling["AR"]["fitted_slope"],
        "scaling_TO": scaling["TO"]["fitted_slope"],
        "scaling_AR_self_consistent": bool(
            abs(scaling["AR"]["fitted_slope"] - exp_ar) < 0.6),
        "scaling_TO_self_consistent": bool(
            abs(scaling["TO"]["fitted_slope"] - exp_to) < 0.6),
        "scaling_AR_ok": bool(abs(scaling["AR"]["fitted_slope"] - exp_ar) < 0.6),
        "scaling_TO_ok": bool(abs(scaling["TO"]["fitted_slope"] - exp_to) < 0.6),
        "AR_erasure_frac": zone_result["AR"]["erasure_frac_of_errors"],
        "AR_erasure_frac_paper": fr_ar,
        "TO_loss_frac": zone_result["TO"]["loss_frac_of_errors"],
        "TO_loss_frac_paper": fr_to,
        "rydberg_decay_fraction": rydberg_decay_fraction(),
    }


def _wilson(p: float, n: int) -> float:
    z = 1.96
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return float(half)
