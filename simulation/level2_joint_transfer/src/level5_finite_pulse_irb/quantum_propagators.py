"""批量 SU(2) 解析传播器与 Lindblad 密度矩阵传播器。

核心约定：一批独立实例的态矢存为 [B, 2] 复数数组，所有单步传播在
步内 Hamiltonian 视为常数的假设下使用解析 SU(2) 指数（中点取样）。
自由演化段（Ω=0，纯 σ_z）用相位矩阵向量积一次性累积，脉冲段逐子步推进；
两条路径共用同一 Hamiltonian 组装，不为任何 pulse family 复制传播器。
"""
from __future__ import annotations

import numpy as np

from .qubit_hamiltonian import (HBAR, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                validate_statevector)


def su2_propagator(delta, omega, phi, dt):
    """批量解析 SU(2) 步进矩阵元。

    delta/omega/phi 为 [B] 控制量，dt 标量（秒）。返回 (U00, U01, U10, U11)
    四个 [B] 复数数组，U=exp[-i H dt/ħ]（H=ħ/2[Δσ_z+Ω(cosφσ_x+sinφσ_y)]）。
    """
    delta = np.asarray(delta, dtype=float)
    omega = np.asarray(omega, dtype=float)
    nx = omega * np.cos(phi)
    ny = omega * np.sin(phi)
    n_norm = np.sqrt(nx * nx + ny * ny + delta * delta)
    theta = n_norm * dt
    c = np.cos(0.5 * theta)
    # sin(θ/2)/|n|，|n|→0 时解析极限为 dt/2（零角度安全）
    tiny = n_norm < 1e-300
    safe = np.where(tiny, 1.0, n_norm)
    s = np.where(tiny, 0.5 * dt, np.sin(0.5 * theta) / safe)
    u00 = c - 1j * s * delta
    u11 = c + 1j * s * delta
    u01 = -1j * s * (nx - 1j * ny)
    u10 = -1j * s * (nx + 1j * ny)
    return u00, u01, u10, u11


def apply_su2(psi, u00, u01, u10, u11):
    """把 su2_propagator 的矩阵元作用在 [B, 2] 复数态矢上（就地修改）。"""
    if not np.iscomplexobj(psi):
        raise TypeError("psi 必须是复数数组（实数 dtype 会静默丢弃相位）")
    a = psi[..., 0].copy()
    b = psi[..., 1].copy()
    psi[..., 0] = u00 * a + u01 * b
    psi[..., 1] = u10 * a + u11 * b
    return psi


def apply_z_rotation(psi, angle):
    """绕 z 的批量旋转 exp(-iθσ_z/2)（angle 为 [B] 数组）。"""
    half = 0.5 * np.asarray(angle, dtype=float)
    rot = np.exp(-1j * half)
    psi[..., 0] *= rot
    psi[..., 1] *= np.conj(rot)
    return psi


def propagate_free_segments(psi, delta_per_segment, seg_durations):
    """纯 σ_z 自由演化的批量相位累积。

    delta_per_segment[k] 为第 k 段每个样本点（0.5 μs 网格中点）的 Δ [B, L_k]
    或标量；seg_durations[k] 为段时长。相位 Φ=∫Δdt 用矩形中点权重累积后
    一次性施加 z 旋转（与逐步施加在分段常数假设下严格一致）。
    """
    for delta, duration in zip(delta_per_segment, seg_durations):
        arr = np.asarray(delta, dtype=float)
        if arr.ndim == 0:
            phi_total = float(arr) * duration
        else:
            dt = duration / arr.shape[-1]
            phi_total = arr.sum(axis=-1) * dt
        apply_z_rotation(psi, phi_total)
    return psi


def propagate_pulse_segment(psi, omega_fn, delta_fn, t_start, t_end, dt_max,
                            phi):
    """一个有限时长脉冲的逐子步传播。

    omega_fn(t_mid) 与 delta_fn(t_mid) 返回 [B] 控制量；子步长不超过
    dt_max，段边界精确落在 t_start/t_end。envelope 由 omega_fn 内实现。
    """
    duration = t_end - t_start
    n_sub = max(1, int(np.ceil(duration / dt_max - 1e-12)))
    dt = duration / n_sub
    for k in range(n_sub):
        t_mid = t_start + (k + 0.5) * dt
        delta = np.asarray(delta_fn(t_mid), dtype=float)
        omega = np.asarray(omega_fn(t_mid), dtype=float)
        u = su2_propagator(delta, omega, phi, dt)
        apply_su2(psi, *u)
    return psi


def unitarity_error(u00, u01, u10, u11):
    """批量 2×2 传播子的 unitarity 误差 max|U†U − I|（抽样诊断用）。"""
    u = np.empty(np.shape(u00) + (2, 2), dtype=complex)
    u[..., 0, 0] = u00
    u[..., 0, 1] = u01
    u[..., 1, 0] = u10
    u[..., 1, 1] = u11
    ud = np.conj(np.swapaxes(u, -1, -2))
    prod = ud @ u
    eye = np.eye(2, dtype=complex)
    return float(np.max(np.abs(prod - eye)))


# ---------------------------------------------------------------- Lindblad
def lindblad_step(rho, delta, omega, phi, dt, collapse_ops=()):
    """单步密度矩阵传播（RK4 哈密顿部分 + 解析 Lindblad 一阶合并）。

    rho: [2, 2]；collapse_ops 为 (L, rate) 列表，其中 L 为无量纲算符、
    rate 具有时间倒数单位（每个 rate 必须有配置 provenance）。仅用于
    小规模可选模型（T1 / Markovian 纯退相位 / leakage surrogate）。
    """
    h = 0.5 * (delta * PAULI_Z + omega * (np.cos(phi) * PAULI_X
                                          + np.sin(phi) * PAULI_Y))
    k1 = -1j * (h @ rho - rho @ h)
    k2 = -1j * (h @ (rho + 0.5 * dt * k1) - (rho + 0.5 * dt * k1) @ h)
    k3 = -1j * (h @ (rho + 0.5 * dt * k2) - (rho + 0.5 * dt * k2) @ h)
    k4 = -1j * (h @ (rho + dt * k3) - (rho + dt * k3) @ h)
    rho = rho + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    for op, rate in collapse_ops:
        op_dag = op.conj().T
        rho = rho + dt * rate * (op @ rho @ op_dag
                                 - 0.5 * (op_dag @ op @ rho
                                          + rho @ op_dag @ op))
    return rho


def density_matrix_checks(rho, eigen_tolerance=-1e-10):
    """密度矩阵 trace/Hermiticity/正定性检查。"""
    rho = np.asarray(rho)
    trace_dev = abs(np.trace(rho) - 1.0)
    herm_dev = float(np.max(np.abs(rho - rho.conj().T)))
    eig_min = float(np.min(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))))
    return {"trace_deviation": trace_dev,
            "hermiticity_max_deviation": herm_dev,
            "min_eigenvalue": eig_min,
            "positive": bool(eig_min >= eigen_tolerance),
            "trace_preserved": bool(trace_dev <= 1e-8),
            "hermitian": bool(herm_dev <= 1e-10)}


def depolarizing_ptm(p_depol):
    """已知去极化通道的 3×3 PTM（用于 preflight 参数回收）。"""
    lam = 1.0 - p_depol
    return lam * np.eye(3)


def dephasing_ptm(p_dephase):
    """已知纯退相位通道的 3×3 PTM。"""
    lam = 1.0 - p_dephase
    return np.diag([lam, lam, 1.0])
