"""批量态矢量 Clifford MC：小系统（n≤10）线路级仿真核心。

设计（服务 R3/Y5/Y6）：
- 形状 (2,)*n + (n_trials,)，一次演化全部 trial；
- Clifford 门 H/S/X/Z/CX/CZ 以轴运算实现，精确；
- 逐 trial 的 Pauli 误差 ⊗_q X^fx Z^fz（Y=i·X·Z，全幅补 i）；
- 丢失比特：门恒等、测量返回 dark(−1) 且不做坍缩——幸存者约化态
  = 对丢失指标求迹，保持正确；
- 破坏性单比特测量（实验口径），含坍缩与非破坏性测量不引入。
"""
from __future__ import annotations

import numpy as np

_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)
_S = np.array([[1, 0], [0, 1j]], dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)


class BatchedState:
    def __init__(self, n: int, n_trials: int, rng: np.random.Generator):
        if n < 1 or n > 10:
            raise ValueError("batched state-vector 仅用于 1≤n≤10")
        self.n = n
        self.trials = n_trials
        self.rng = rng
        self.psi = np.zeros((2,) * n + (n_trials,), dtype=complex)
        self.psi[(0,) * n] = 1.0
        self.lost = np.zeros((n, n_trials), dtype=bool)

    # ---------------- Clifford 门 ----------------
    def _qubit_apply(self, q: int, op: np.ndarray) -> None:
        moved = np.moveaxis(self.psi, q, 0)
        out = op @ moved.reshape(2, -1)
        out = out.reshape(moved.shape)
        if self.lost[q].any():
            out = np.where(self.lost[q], moved, out)
        self.psi = np.moveaxis(out, 0, q)

    def h(self, q: int) -> None:
        self._qubit_apply(q, _H)

    def s(self, q: int) -> None:
        self._qubit_apply(q, _S)

    def x(self, q: int) -> None:
        self._qubit_apply(q, _X)

    def z(self, q: int) -> None:
        self._qubit_apply(q, _Z)

    def cx(self, c: int, t: int) -> None:
        skip = self.lost[c] | self.lost[t]
        moved = np.moveaxis(self.psi, (c, t), (0, 1))
        out = moved.copy()
        out[1, 1], out[1, 0] = moved[1, 0], moved[1, 1]
        if skip.any():
            out = np.where(skip, moved, out)
        self.psi = np.moveaxis(out, (0, 1), (c, t))

    def cz(self, a: int, b: int) -> None:
        skip = self.lost[a] | self.lost[b]
        moved = np.moveaxis(self.psi, (a, b), (0, 1))
        out = moved.copy()
        out[1, 1] = -moved[1, 1]
        if skip.any():
            out = np.where(skip, moved, out)
        self.psi = np.moveaxis(out, (0, 1), (a, b))

    # ---------------- 误差与丢失 ----------------
    def apply_pauli_errors(self, fx: np.ndarray, fz: np.ndarray) -> None:
        """逐 trial 施加 ⊗_q X^fx[q] Z^fz[q]。fx/fz: (n, trials) bool。"""
        for q in range(self.n):
            both = fx[q] & fz[q]
            if fz[q].any():
                moved = np.moveaxis(self.psi, q, 0)
                moved[1] = moved[1] * np.where(fz[q], -1.0, 1.0)
                self.psi = np.moveaxis(moved, 0, q)
            if fx[q].any():
                moved = np.moveaxis(self.psi, q, 0)
                swapped = np.moveaxis(moved[[1, 0]], 0, q)
                self.psi = np.where(fx[q], swapped, self.psi)
            if both.any():
                self.psi = self.psi * np.where(both, 1j, 1.0)

    def mark_lost(self, lost_mask: np.ndarray) -> None:
        self.lost |= lost_mask

    # ---------------- 测量 ----------------
    def measure_z(self, q: int) -> np.ndarray:
        """破坏性 Z 测量。丢失 trial 返回 −1（dark）且不坍缩幸存者状态。"""
        lost = self.lost[q]
        moved = np.moveaxis(np.abs(self.psi) ** 2, q, 0)
        p1 = moved[1].sum(axis=tuple(range(self.n - 1)))
        p1 = np.where(lost, 0.0, p1)
        bits = (self.rng.random(self.trials) < p1).astype(np.int8)
        amp = np.moveaxis(self.psi, q, 0)
        active = ~lost
        keep0 = active & (bits == 0)
        keep1 = active & (bits == 1)
        amp[0] = np.where(keep1, 0, amp[0])
        amp[1] = np.where(keep0, 0, amp[1])
        n0 = np.sqrt((np.abs(amp[0]) ** 2).sum(axis=tuple(range(self.n - 1))))
        n1 = np.sqrt((np.abs(amp[1]) ** 2).sum(axis=tuple(range(self.n - 1))))
        amp[0] = amp[0] / np.maximum(np.where(keep0, n0, 1.0), 1e-300)
        amp[1] = amp[1] / np.maximum(np.where(keep1, n1, 1.0), 1e-300)
        self.psi = np.moveaxis(amp, 0, q)
        return np.where(lost, np.int8(-1), bits)

    def measure_x(self, q: int) -> np.ndarray:
        self.h(q)   # _qubit_apply 已保护丢失位
        return self.measure_z(q)

    # ---------------- 诊断观测 ----------------
    def diag_expectation(self, spec: str) -> np.ndarray:
        """真实联合对角期望 ⟨⊗_q O_q⟩（I/Z 串），逐 trial。

        修正（验证 D1）：旧实现为逐比特边缘期望的乘积，对纠缠态错误
        （Bell 态 "ZZ" 返回 0，真值 1）。现按概率表逐串求和。
        仅支持对角（I/Z）观测；非对角请用测量或专用旋转。
        """
        prob = np.abs(self.psi) ** 2          # (2,)*n + (trials,)
        flat = prob.reshape(2 ** self.n, -1)
        acc = np.zeros(self.trials)
        for idx in range(2 ** self.n):
            bits = [(idx >> (self.n - 1 - q)) & 1 for q in range(self.n)]
            eig = 1.0
            for q, ch in enumerate(spec):
                if ch == "Z":
                    eig *= 1 - 2 * bits[q]
                elif ch != "I":
                    raise ValueError(f"diag_expectation 仅支持 I/Z，收到 {ch!r}")
            acc += eig * flat[idx]
        return acc

    def apply_z_rotation(self, q: int, angles: np.ndarray) -> None:
        """逐 trial 施加 R_z(θ_q)：|0⟩→e^{−iθ/2}|0⟩，|1⟩→e^{+iθ/2}|1⟩。

        angles 形状 (trials,)。丢失位上恒等。
        """
        theta = np.asarray(angles, dtype=float)
        moved = np.moveaxis(self.psi, q, 0)
        ph0 = np.exp(-0.5j * theta)
        ph1 = np.exp(+0.5j * theta)
        out0 = moved[0] * ph0
        out1 = moved[1] * ph1
        moved = np.stack([out0, out1], axis=0)
        if self.lost[q].any():
            orig = np.moveaxis(self.psi, q, 0)
            moved = np.where(self.lost[q], orig, moved)
        self.psi = np.moveaxis(moved, 0, q)

    def probabilities(self) -> np.ndarray:
        return np.abs(self.psi.reshape(2 ** self.n, -1)) ** 2
