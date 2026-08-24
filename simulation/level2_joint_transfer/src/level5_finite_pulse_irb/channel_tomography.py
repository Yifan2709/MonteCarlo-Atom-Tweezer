"""单次通道的过程表征：survival-conditioned PTM、保真度与诊断。

对 {|0⟩,|1⟩,|+x⟩,|−x⟩,+y,−y} 六个输入态沿冻结轨迹传播，重建仿射
Bloch 映射 r→Mr+c，输出平均态保真度、unitarity、相干旋转分量、收缩轴、
survival 与条件/联合指标。非 trace-preserving 的损失映射不被归一化：
survival、条件自旋通道与 survival-weighted 性能并列报告。
"""
from __future__ import annotations

import numpy as np

from .qubit_hamiltonian import INPUT_STATES, statevector_to_bloch

CARDINAL_INPUTS = ("z_plus", "z_minus", "x_plus", "x_minus",
                   "y_plus", "y_minus")
AXIS_VECTORS = {"x": np.array([1.0, 0.0, 0.0]),
                "y": np.array([0.0, 1.0, 0.0]),
                "z": np.array([0.0, 0.0, 1.0])}


def input_basis_vectors():
    """六个输入态的 Bloch 向量。"""
    out = {}
    for name, psi in INPUT_STATES.items():
        out[name] = statevector_to_bloch(psi)
    return out


def affine_fit(bloch_in: np.ndarray, bloch_out: np.ndarray):
    """六个基输入 → 仿射映射 M, c（含非仿射残差诊断）。

    bloch_in/bloch_out: [6, 3]，行序 = CARDINAL_INPUTS，分量序 (x,y,z)。
    M 列序与 Bloch 分量一致：列 0/1/2 对应 x/y/z 输入轴。
    """
    bloch_in = np.asarray(bloch_in, dtype=float)
    bloch_out = np.asarray(bloch_out, dtype=float)
    pairs = ((2, 3), (4, 5), (0, 1))  # (x+,x−), (y+,y−), (z+,z−) 行索引
    m = np.zeros((3, 3))
    for i, (ip, im) in enumerate(pairs):
        m[:, i] = 0.5 * (bloch_out[ip] - bloch_out[im])
    c = bloch_out.mean(axis=0)
    residual = bloch_out - (bloch_in @ m.T + c)
    return m, c, float(np.max(np.abs(residual)))


def average_state_fidelity_from_affine(m: np.ndarray, c: np.ndarray):
    """Haar 平均纯态保真度 F_avg = 1/2 + Tr(M)/6（对任意仿射映射成立）。

    若通道非 unital（‖c‖ 显著非零）或非仿射，需并列报告诊断而不是只报
    该数字。
    """
    return 0.5 + float(np.trace(m)) / 6.0


def unitarity_from_affine(m: np.ndarray, c: np.ndarray):
    """Unitarity（Wall et al.）：unital 块 Tr(M†M)/3；c≠0 时注明包含
    非保范分量（对 TP but non-unital 通道此为广义估计）。"""
    u = float(np.trace(m.T @ m)) / 3.0
    return {"unitarity": u, "translation_norm": float(np.linalg.norm(c)),
            "unital": bool(np.linalg.norm(c) < 1e-6),
            "note": "c≠0（T1 类非 unital 通道）时 Tr(M†M)/3 为广义估计"}


def coherent_rotation_component(m: np.ndarray):
    """M 的最近旋转（极分解正交部分）与收缩奇异值。"""
    u2, s2, vt2 = np.linalg.svd(m)
    det = np.linalg.det(u2 @ vt2)
    rotation = u2 @ np.diag([1.0, 1.0, det]) @ vt2
    angle = float(np.arccos(np.clip(0.5 * (np.trace(rotation) - 1.0),
                                    -1.0, 1.0)))
    axis = None
    if angle > 1e-9:
        raw = np.array([rotation[2, 1] - rotation[1, 2],
                        rotation[0, 2] - rotation[2, 0],
                        rotation[1, 0] - rotation[0, 1]])
        norm = np.linalg.norm(raw)
        if norm > 1e-12:
            axis = raw / norm
    return {"coherent_angle_rad": angle,
            "coherent_axis": None if axis is None else axis.tolist(),
            "singular_values": sorted(s2.tolist(), reverse=True),
            "contraction_anisotropy": float(s2.max() - s2.min())}


def fidelity_applicability(m, c, affine_residual):
    """F_avg 公式适用性与通道口径诊断。"""
    diag = {
        "affine_residual": affine_residual,
        "affine_ok": bool(affine_residual < 1e-6),
        "unital": bool(np.linalg.norm(c) < 1e-6),
        "mixture_of_unitaries_tp": bool(np.linalg.norm(c) < 1e-6
                                        and affine_residual < 1e-6),
        "note": ("条件自旋通道为幺正系综平均（TP+unital）时 "
                 "F_avg=1/2+Tr(M)/6 严格适用；含 Lindblad T1 时 c≠0，"
                 "该值解释为 Haar 平均纯态保真度"),
    }
    return diag


def characterize_from_states(bloch_out_matrix, survival_labels=None):
    """从 [6, 3]（或带 realization 维 [R, 6, 3]）输出 Bloch 向量估计通道。"""
    bloch_out_matrix = np.asarray(bloch_out_matrix, dtype=float)
    if bloch_out_matrix.ndim == 2:
        bloch_out_matrix = bloch_out_matrix[None, ...]
    m_list, c_list, res_list = [], [], []
    for r in range(bloch_out_matrix.shape[0]):
        m, c, res = affine_fit(_reference_inputs(), bloch_out_matrix[r])
        m_list.append(m)
        c_list.append(c)
        res_list.append(res)
    m_mean = np.mean(m_list, axis=0)
    c_mean = np.mean(c_list, axis=0)
    out = {
        "ptm_m": m_mean.tolist(),
        "affine_c": c_mean.tolist(),
        "f_avg": average_state_fidelity_from_affine(m_mean, c_mean),
        "unitarity": unitarity_from_affine(m_mean, c_mean),
        "coherent": coherent_rotation_component(m_mean),
        "applicability": fidelity_applicability(
            m_mean, c_mean, float(np.mean(res_list))),
        "n_realizations": int(bloch_out_matrix.shape[0]),
    }
    if survival_labels is not None:
        surv = float(np.mean(np.asarray(survival_labels, dtype=bool)))
        out["survival"] = surv
        out["f_avg_survival_weighted"] = surv * out["f_avg"]
        out["f_avg_loss_detected_convention"] = surv * out["f_avg"]
        out["conventions"] = {
            "conditional_spin_channel": "以幸存为条件的自旋通道 PTM/F_avg",
            "survival": "经典往返成功概率（Level 4 标签）",
            "survival_weighted": "survival × F_avg_conditional",
            "loss_detected": "把丢失计为返回失败 → 与 survival-weighted "
                             "数值一致，但语义为返回概率口径",
        }
    return out


def _reference_inputs():
    """六输入态的 Bloch 向量矩阵（顺序与 CARDINAL_INPUTS 一致）。"""
    vecs = input_basis_vectors()
    return np.stack([vecs[name] for name in CARDINAL_INPUTS])


def ptm_bootstrap(bloch_by_realization: np.ndarray, survival_labels,
                  n_boot=2000, seed=0):
    """跨 classical/noise realization 的非参数 bootstrap 区间。"""
    rng = np.random.default_rng(seed)
    bloch_by_realization = np.asarray(bloch_by_realization, dtype=float)
    n = bloch_by_realization.shape[0]
    labels = np.asarray(survival_labels, dtype=bool)
    f_vals, survs = [], []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        sub = bloch_by_realization[idx].mean(axis=0)
        m, c, _ = affine_fit(_reference_inputs(), sub)
        f_vals.append(average_state_fidelity_from_affine(m, c))
        survs.append(float(labels[idx].mean()))
    f_vals = np.asarray(f_vals)
    survs = np.asarray(survs)
    return {"f_avg_ci_low": float(np.percentile(f_vals, 2.5)),
            "f_avg_ci_high": float(np.percentile(f_vals, 97.5)),
            "survival_ci_low": float(np.percentile(survs, 2.5)),
            "survival_ci_high": float(np.percentile(survs, 97.5)),
            "f_avg_survival_weighted_ci_low": float(np.percentile(
                f_vals * survs, 2.5)),
            "f_avg_survival_weighted_ci_high": float(np.percentile(
                f_vals * survs, 97.5)),
            "n_bootstrap": int(n_boot), "seed": int(seed)}


def known_channel_checks():
    """已知去极化/退相位通道的 PTM 与 F_avg（preflight 用）。"""
    p = 0.1
    m_depol = (1.0 - p) * np.eye(3)
    f_depol = average_state_fidelity_from_affine(m_depol, np.zeros(3))
    f_depol_exact = 1.0 - p / 2.0
    q = 0.2
    m_deph = np.diag([1.0 - q, 1.0 - q, 1.0])
    f_deph = average_state_fidelity_from_affine(m_deph, np.zeros(3))
    f_deph_exact = 1.0 - q / 3.0
    return {"depolarizing": {"p": p, "f_avg": f_depol,
                             "f_avg_exact": f_depol_exact,
                             "match": abs(f_depol - f_depol_exact) < 1e-12},
            "dephasing": {"p": q, "f_avg": f_deph,
                          "f_avg_exact": f_deph_exact,
                          "match": abs(f_deph - f_deph_exact) < 1e-12}}
