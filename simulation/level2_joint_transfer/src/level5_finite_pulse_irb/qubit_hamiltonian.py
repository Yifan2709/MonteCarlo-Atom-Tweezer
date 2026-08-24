"""旋转框架单量子比特哈密顿量（RWA）。

H_i(t) = ħ/2 [Δ_i(t)σ_z + Ω_i(t)(cosφ σ_x + sinφ σ_y)]，
Δ_i(t) = Δ_drive + Δ_LS,i(t) + Δ_B,i(t) + Δ_site,i，
其中轨迹相关差分光移 Δ_LS = (η_SLM·U_SLM + η_AOD·U_AOD)/ħ 沿 Level 4
经典轨迹取值。大载频（9.2 GHz 钟频）不进入步进，只体现在旋转框架失谐中。
不含 Bloch-Siegert 位移、多能级泄漏与自发 Raman 散射（模型边界见报告）。
"""
from __future__ import annotations

import numpy as np

HBAR = 1.054571817e-34
TWO_PI = 2.0 * np.pi

# Pauli 矩阵（convention: σ_z|0>=+|0>, σ_z|1>=-|1>，|0>=|F=3,mF=0>）
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
PAULI_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
IDENTITY = np.eye(2, dtype=complex)

# 六个 PTM 表征输入态
INPUT_STATES = {
    "z_plus": np.array([1.0, 0.0], dtype=complex),
    "z_minus": np.array([0.0, 1.0], dtype=complex),
    "x_plus": np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0),
    "x_minus": np.array([1.0, -1.0], dtype=complex) / np.sqrt(2.0),
    "y_plus": np.array([1.0, 1j], dtype=complex) / np.sqrt(2.0),
    "y_minus": np.array([1.0, -1j], dtype=complex) / np.sqrt(2.0),
}


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """[A, B]。"""
    return a @ b - b @ a


def rotation_unitary(axis: str, angle: float) -> np.ndarray:
    """绕 x/y/z 轴的解析旋转 exp(-iθσ_a/2)。"""
    sigma = {"x": PAULI_X, "y": PAULI_Y, "z": PAULI_Z}[axis]
    return np.cos(angle / 2.0) * IDENTITY - 1j * np.sin(angle / 2.0) * sigma


def equatorial_rotation_unitary(phi: float, angle: float) -> np.ndarray:
    """赤道面旋转 R_φ(θ)=exp[-iθ(cosφσ_x+sinφσ_y)/2]。"""
    nx, ny = np.cos(phi), np.sin(phi)
    generator = nx * PAULI_X + ny * PAULI_Y
    return np.cos(angle / 2.0) * IDENTITY - 1j * np.sin(angle / 2.0) * generator


def hamiltonian_matrix(delta_rad_s, omega_rad_s, phi_rad) -> np.ndarray:
    """给定瞬时控制量返回 2×2 H/ħ（标量输入，测试/验证用）。"""
    return 0.5 * (delta_rad_s * PAULI_Z
                  + omega_rad_s * (np.cos(phi_rad) * PAULI_X
                                   + np.sin(phi_rad) * PAULI_Y))


def light_shift_detuning(u_slm_j, u_aod_j, eta_slm, eta_aod) -> np.ndarray:
    """差分光移角频率 δω=(η_SLM·U_SLM+η_AOD·U_AOD)/ħ（rad/s）。

    u_slm_j/u_aod_j 为沿轨迹的瞬时阱深能量（J，负值），支持数组。
    """
    return (eta_slm * np.asarray(u_slm_j, dtype=float)
            + eta_aod * np.asarray(u_aod_j, dtype=float)) / HBAR


def statevector_to_bloch(psi: np.ndarray) -> np.ndarray:
    """态矢 → Bloch 向量（支持 [..., 2] 批量）。"""
    psi = np.asarray(psi)
    a, b = psi[..., 0], psi[..., 1]
    return np.stack([2.0 * np.real(np.conj(a) * b),
                     2.0 * np.imag(np.conj(a) * b),
                     np.abs(a) ** 2 - np.abs(b) ** 2], axis=-1)


def return_probability_z(psi: np.ndarray) -> np.ndarray:
    """回到 |0> 的概率 |⟨0|ψ⟩|²（支持批量）。"""
    psi = np.asarray(psi)
    return np.abs(psi[..., 0]) ** 2 / np.maximum(
        np.abs(psi[..., 0]) ** 2 + np.abs(psi[..., 1]) ** 2, 1e-300)


def validate_statevector(psi: np.ndarray, norm_tolerance=1e-10) -> dict:
    """态矢数值检查：有限性、归一。"""
    psi = np.asarray(psi)
    norm = np.sqrt(np.sum(np.abs(psi) ** 2, axis=-1))
    finite = bool(np.all(np.isfinite(psi)))
    return {"finite": finite,
            "norm_max_abs_deviation": float(np.max(np.abs(norm - 1.0)))
            if psi.size else 0.0,
            "norm_ok": finite and bool(np.all(np.abs(norm - 1.0)
                                              <= norm_tolerance))}


def analytic_rabi_population(delta_rad_s, omega_rad_s, t_s) -> float:
    """常数 H 下 P_1(t) 的解析 Rabi 振荡（含失谐）。"""
    om_eff = np.sqrt(omega_rad_s ** 2 + delta_rad_s ** 2)
    return (omega_rad_s ** 2 / om_eff ** 2) * np.sin(0.5 * om_eff * t_s) ** 2


def analytic_ramsey_phase(delta_rad_s, t_s, echo=False) -> float:
    """常数失谐下的 Ramsey/echo 累积相位（echo 时符号相消）。"""
    if echo:
        return 0.0
    return delta_rad_s * t_s
