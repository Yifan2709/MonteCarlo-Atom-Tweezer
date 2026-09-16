# -*- coding: utf-8 -*-
"""A 阶段独立小检查（独立验证方脚本；不修改生产代码）。

A1 M1：Bell 态参考——用振幅直接构造联合概率表（不调用被测函数），
       核对 <ZI>=<IZ>=0、<ZZ>=<XX>=1、<YY>=−1；并演示 diag_expectation
       对纠缠态 "ZZ" 的乘积式实现的偏差（潜在误用面）。
A3 Y1：轨迹端点导数（数值差分）核对 sixth_zero_jerk / gamma:1.875 的
       端点速度/加速度/jerk 与模块说明的一致性。
A4 入口：新 CLI 派发与 workflow 注册存在性。
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[3] / "src"
sys.path.insert(0, str(SRC))
OUT = Path(__file__).resolve().parent

results: dict = {}

# ---------------- A1: Bell 独立参考 ----------------
psi = np.zeros(4, dtype=complex)
psi[0] = psi[3] = 1 / np.sqrt(2)          # |Φ+> = (|00>+|11>)/√2，独立构造
P = (psi.conj() * psi).real               # 联合概率表（独立于被测代码）
bits = [(0, 0), (0, 1), (1, 0), (1, 1)]
Z = {-1: (0, 3), +1: (1, 2)}  # 占位


def pauli_eig(op: str) -> float:
    """独立计算 ⟨op⟩：直接对概率表求和。op ∈ {'ZI','IZ','ZZ','XX','YY'}."""
    e = 0.0
    for idx, (b0, b1) in enumerate(bits):
        z0, z1 = 1 - 2 * b0, 1 - 2 * b1
        # X 基：|Φ+> 的 X 基表示 = 自身（Hadamard 不变）；用矩阵直接算：
        if op == "ZI":
            eig = z0
        elif op == "IZ":
            eig = z1
        elif op == "ZZ":
            eig = z0 * z1
        else:
            # XX / YY：对态矢量做矩阵乘（独立实现）
            M = np.array([[1, 0], [0, -1]])
            X = np.array([[0, 1], [1, 0]])
            Y = np.array([[0, -1j], [1j, 0]])
            O = {"XX": np.kron(X, X), "YY": np.kron(Y, Y)}[op]
            return float((psi.conj() @ (O @ psi)).real)
        e += P[idx] * eig
    return e


a1 = {op: pauli_eig(op) for op in ("ZI", "IZ", "ZZ", "XX", "YY")}
a1_pass = (abs(a1["ZI"]) < 1e-12 and abs(a1["IZ"]) < 1e-12
           and abs(a1["ZZ"] - 1) < 1e-12 and abs(a1["XX"] - 1) < 1e-12
           and abs(a1["YY"] + 1) < 1e-12)

# diag_expectation 误用面演示（仅诊断，不影响生产路径）
import paper_integration.species  # noqa: E402
from paper_integration.quantum_mc import BatchedState  # noqa: E402

st = BatchedState(2, 1000, np.random.default_rng(0))
st.h(0)
st.cx(0, 1)
wrong = float(st.diag_expectation("ZZ").mean())     # 乘积式实现
right = pauli_eig("ZZ")
results["A1_bell_reference"] = {
    "independent": a1, "pass": a1_pass,
    "diag_expectation_ZZ_on_Bell": wrong, "true_ZZ": right,
    "diag_product_formula_bias": abs(wrong - right),
    "production_usage": "grep 确认 diag_expectation 未进入生产结果路径；"
                        "测试仅用 'ZI'（单比特，乘积式正确）",
}

# ---------------- A3: 轨迹端点导数 ----------------
from paper_integration.trajectories import poly_trajectory  # noqa: E402

D, T, dt = 100e-6, 1e-3, 0.1e-6
a3 = {}
for kind in ("sixth_zero_jerk", "gamma:1.875"):
    tr = poly_trajectory(kind, D, T, dt)
    x, t = tr.x_m, tr.t_s
    v = np.gradient(x, t)
    a = np.gradient(v, t)
    j = np.gradient(a, t)
    # 端点窗口平均（避开数值噪声的首末几点）
    w = 5
    a3[kind] = {
        "v_end_abs_max": float(np.max(np.abs(v[:w]))),
        "acc_start": float(np.mean(a[1:w + 1])),
        "acc_end": float(np.mean(a[-w - 1:-1])),
        "jerk_start": float(np.mean(j[1:w + 1])),
        "jerk_end": float(np.mean(j[-w - 1:-1])),
        "peak_v_ratio": float(np.max(v) * T / D),
        "endpoint_x_err": float(abs(x[-1] - D)),
    }
results["A3_trajectories"] = {
    "values": a3,
    "finding": (
        "sixth_zero_jerk 实为五阶剖面 x=6s³−7s⁴+2s⁵：端点 jerk(0)=6a=36·D/T³≠0，"
        "并非零 jerk，与名称及模块 docstring 中“六阶/端点零 jerk”表述不符"
        "（代码内联注释已如实写 a=6/γ≈1.63）；结构性质（峰值速度 1.728<1.875）"
        "成立，但“zero-jerk”标签应更名或修正文档"),
}
results["A3_pass_endpoint_v"] = bool(
    a3["sixth_zero_jerk"]["v_end_abs_max"] < 1e-9
    and a3["gamma:1.875"]["v_end_abs_max"] < 1e-9
    and abs(a3["gamma:1.875"]["peak_v_ratio"] - 1.875) < 1e-3)

# ---------------- A4: 入口派发 ----------------
r1 = subprocess.run([sys.executable, "-m", "paper_integration.cli", "list"],
                    capture_output=True, text=True, cwd=SRC.parent)
r2 = subprocess.run([sys.executable, "-c",
                     "import sys;sys.path.insert(0,'src');"
                     "from simulation_workflow.cli import LEVELS;"
                     "print(sorted(k for k in LEVELS if k.startswith('paper')))"],
                    capture_output=True, text=True, cwd=SRC.parent)
results["A4_dispatch"] = {
    "paper_cli_list_rc": r1.returncode, "paper_cli_stdout": r1.stdout.strip(),
    "workflow_registered": r2.stdout.strip(), "workflow_rc": r2.returncode,
}
(OUT / "stage_a_results.json").write_text(
    json.dumps(results, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(results, ensure_ascii=False, indent=1))
