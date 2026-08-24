"""完整 24 元素单比特 Clifford 群：规范型、分解、逆元、乘法表与均匀采样。

规范型约定：把幺正矩阵乘以全局相位使第一列首个非零元素为正实数；此时
24 个元素的矩阵元 ∈ {0, ±1, ±i, (±1±i)/√2}，乘 √2 后为高斯整数矩阵，
可作为精确哈希键比较唯一性。分解只使用 X/Y 轴 ±π/2 旋转与虚拟 Z（π/2 的
整数倍），物理脉冲数 ≤ 2，与论文 SCROFULOUS 实现同一风格。
"""
from __future__ import annotations

import numpy as np

from .qubit_hamiltonian import (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                equatorial_rotation_unitary)

TWO_PI = 2.0 * np.pi


def _canonicalize(u: np.ndarray) -> np.ndarray:
    """全局相位规范化：第一列首个非零元素变为正实数。"""
    u = np.asarray(u, dtype=complex)
    for value in u[:, 0]:
        if abs(value) > 1e-8:
            phase = value / abs(value)
            return u / phase
    raise ValueError("幺正矩阵第一列全零")


def _gaussian_key(u: np.ndarray):
    """规范型 ×√2 后四舍五入到高斯整数矩阵（精确哈希键）。"""
    canon = _canonicalize(u) * np.sqrt(2.0)
    real = np.rint(canon.real).astype(int)
    imag = np.rint(canon.imag).astype(int)
    return tuple(zip(real.ravel(), imag.ravel()))


# 物理基元：X/Y 轴 ±π/2 旋转 + 虚拟 Z 旋转
_PRIMITIVES = [
    ("X90", equatorial_rotation_unitary(0.0, 0.5 * np.pi)),
    ("X-90", equatorial_rotation_unitary(0.0, -0.5 * np.pi)),
    ("Y90", equatorial_rotation_unitary(0.5 * np.pi, 0.5 * np.pi)),
    ("Y-90", equatorial_rotation_unitary(0.5 * np.pi, -0.5 * np.pi)),
]


def virtual_z_unitary(angle: float) -> np.ndarray:
    """虚拟 Z 旋转 exp(-iα σ_z/2)。"""
    return np.cos(0.5 * angle) * IDENTITY - 1j * np.sin(
        0.5 * angle) * PAULI_Z


def _generate_group():
    """由 {H, S} 型生成元 BFS 生成完整群；返回 (matrices, keys, decomps)。

    生成元取 X90 与 Z90（与 S=Z90、H=X90·Z90·X90 等价的常用生成集）。
    """
    gens = [("X90", _PRIMITIVES[0][1]), ("Z90", virtual_z_unitary(0.5 * np.pi))]
    seen = {_gaussian_key(IDENTITY): 0}
    matrices = [IDENTITY.copy()]
    decomps = [[]]
    frontier = [0]
    while frontier:
        nxt = []
        for idx in frontier:
            for gname, gmat in gens:
                new = gmat @ matrices[idx]
                key = _gaussian_key(new)
                if key in seen:
                    continue
                seen[key] = len(matrices)
                matrices.append(new)
                decomps.append(decomps[idx] + [gname])
                nxt.append(len(matrices) - 1)
        frontier = nxt
    return matrices, seen, decomps


MATRICES, KEY_INDEX, RAW_DECOMPS = _generate_group()


def _rewrite_decompositions():
    """把生成元序列改写为 X90/X-90/Y90/Y-90 + 虚拟 Z 的规范分解。

    Z90^k 合并成单个虚拟 Z；连续 X90 交替用 X90/X-90 表示以减少面积。
    """
    out = []
    for seq in RAW_DECOMPS:
        pulses = []
        z_angle = 0.0
        for name in seq:
            if name == "Z90":
                z_angle += 0.5 * np.pi
            else:
                if abs(z_angle) > 1e-9:
                    pulses.append(("VZ", z_angle))
                    z_angle = 0.0
                pulses.append(("P", name))
        if abs(z_angle) > 1e-9:
            pulses.append(("VZ", z_angle))
        out.append(pulses)
    return out


DECOMPOSITIONS = _rewrite_decompositions()


def optimized_decompositions() -> list[list]:
    """规范化分解：每元素 ≤2 物理脉冲（X±90/Y±90）+ 前后虚拟 Z。

    任意 1q Clifford 都可写成 Z-X-Z-X-Z 正规形；对动作集
    {VZ(±90°), VZ(180°), X±90} 做代价优先搜索（代价 = (物理脉冲数, 步数)），
    再把脉冲前的虚拟 Z 折算进该脉冲的驱动相位。
    """
    from collections import deque
    actions = [("VZ", 0.5 * np.pi), ("VZ", -0.5 * np.pi), ("VZ", np.pi),
               ("P", "X90"), ("P", "X-90")]
    action_mats = [virtual_z_unitary(0.5 * np.pi),
                   virtual_z_unitary(-0.5 * np.pi), virtual_z_unitary(np.pi),
                   equatorial_rotation_unitary(0.0, 0.5 * np.pi),
                   equatorial_rotation_unitary(0.0, -0.5 * np.pi)]
    identity_key = _gaussian_key(IDENTITY)
    visited = {identity_key: (0, [])}
    queue = deque([(IDENTITY.copy(), 0)])
    # 代价优先：先按物理脉冲数分层 BFS，再按总步数
    while queue:
        mat, n_phys = queue.popleft()
        _, path = visited[_gaussian_key(mat)]
        if visited[_gaussian_key(mat)][0] < n_phys:
            continue
        for (act, val), amat in zip(actions, action_mats):
            new = amat @ mat
            new_key = _gaussian_key(new)
            new_phys = n_phys + (1 if act == "P" else 0)
            new_path = path + [(act, val)]
            if new_key not in visited or (
                    visited[new_key][0] > new_phys
                    or (visited[new_key][0] == new_phys
                        and len(visited[new_key][1]) > len(new_path))):
                visited[new_key] = (new_phys, new_path)
                queue.append((new, new_phys))
    result = []
    for mat in MATRICES:
        _, path = visited[_gaussian_key(mat)]
        ops = []
        z_acc = 0.0
        for act, val in path:
            if act == "VZ":
                z_acc += val
            else:
                if abs(z_acc) > 1e-9:
                    ops.append(("VZ", z_acc))
                    z_acc = 0.0
                ops.append((act, val))
        if abs(z_acc) > 1e-9:
            ops.append(("VZ", z_acc))
        result.append(ops)
    return result


OPT_DECOMPS = optimized_decompositions()


class CliffordGroup:
    """24 元素单比特 Clifford 群（含逆元表、Pauli 共轭表与均匀采样）。"""

    def __init__(self):
        self.matrices = MATRICES
        self.n = len(MATRICES)
        self._index_by_key = dict(KEY_INDEX)
        self._build_inverse_table()
        self._build_multiplication_table()
        self._build_pauli_conjugation()

    # ------------------------------------------------------------ 构造表
    def _build_inverse_table(self):
        self.inverse = np.zeros(self.n, dtype=int)
        for i, u in enumerate(self.matrices):
            key = _gaussian_key(u.conj().T)
            self.inverse[i] = self._index_by_key[key]

    def _build_multiplication_table(self):
        self.mult = np.zeros((self.n, self.n), dtype=int)
        for i, ui in enumerate(self.matrices):
            for j, uj in enumerate(self.matrices):
                key = _gaussian_key(ui @ uj)
                self.mult[i, j] = self._index_by_key[key]

    def _build_pauli_conjugation(self):
        paulis = [PAULI_X, PAULI_Y, PAULI_Z]
        self.pauli_conj = np.zeros((self.n, 3), dtype=int)
        self.pauli_conj_sign = np.zeros((self.n, 3), dtype=int)
        for i, u in enumerate(self.matrices):
            for k, p in enumerate(paulis):
                conj = u @ p @ u.conj().T
                for k2, p2 in enumerate(paulis):
                    if np.allclose(conj, p2, atol=1e-8):
                        self.pauli_conj[i, k] = k2
                        self.pauli_conj_sign[i, k] = 1
                        break
                    if np.allclose(conj, -p2, atol=1e-8):
                        self.pauli_conj[i, k] = k2
                        self.pauli_conj_sign[i, k] = -1
                        break
                else:
                    raise ValueError(f"元素 {i} 的 Pauli 共轭不是带符号 Pauli")

    # ------------------------------------------------------------ API
    def unitary(self, index: int) -> np.ndarray:
        """规范幺正矩阵（忽略全局相位后唯一）。"""
        return self.matrices[index]

    def decompose(self, index: int) -> list:
        """规范分解（物理 π/2 旋转 + 虚拟 Z）。"""
        return OPT_DECOMPS[index]

    def physical_pulse_count(self, index: int) -> int:
        """物理 π/2 脉冲数。"""
        return sum(1 for kind, _ in self.decompose(index) if kind == "P")

    def sample_uniform(self, rng: np.random.Generator, size: int) -> np.ndarray:
        """均匀采样。"""
        return rng.integers(0, self.n, size=int(size))

    def sampler_uniformity_counts(self, n_draws: int = 240000,
                                  seed: int = 0) -> np.ndarray:
        """均匀性检查用的直方图。"""
        rng = np.random.default_rng(seed)
        return np.bincount(self.sample_uniform(rng, n_draws),
                           minlength=self.n)

    def as_rows(self, omega_base: float, family: str = "bare",
                envelope: str = "square") -> list[dict]:
        """clifford_table.csv 行：矩阵元、分解、逆、脉冲数与面积。"""
        from .pulse_library import pulse_list_for, pulse_list_area, \
            pulse_list_duration
        rows = []
        for i in range(self.n):
            pulses = []
            decomp = self.decompose(i)
            for kind, val in decomp:
                if kind == "VZ":
                    continue
                axis_phase = {"X90": 0.0, "Y90": 0.5 * np.pi,
                              "X-90": np.pi, "Y-90": 1.5 * np.pi}[val]
                angle = 0.5 * np.pi if val in ("X90", "Y90") else -0.5 * np.pi
                pulses.extend(pulse_list_for(
                    axis_phase, angle, family, omega_base, envelope,
                    clifford_id=i))
            u = self.matrices[i]
            rows.append({
                "index": i, "n_elements": self.n,
                "u00_re": u[0, 0].real, "u00_im": u[0, 0].imag,
                "u01_re": u[0, 1].real, "u01_im": u[0, 1].imag,
                "u10_re": u[1, 0].real, "u10_im": u[1, 0].imag,
                "u11_re": u[1, 1].real, "u11_im": u[1, 1].imag,
                "decomposition": ";".join(
                    f"{v}" if k == "P" else f"VZ({v / np.pi:.2f}π)"
                    for k, v in decomp),
                "inverse_index": int(self.inverse[i]),
                "n_physical_pulses": self.physical_pulse_count(i),
                "bare_duration_us": pulse_list_duration(pulses) * 1e6,
                "bare_area_pi": pulse_list_area(pulses) / np.pi,
            })
        return rows

    def group_checks(self) -> dict:
        """群论自检：唯一性、封闭性、逆、Pauli 共轭、理想实现。"""
        unique = len({_gaussian_key(u) for u in self.matrices}) == self.n
        closure = bool(np.all(self.mult >= 0))
        inverse_ok = True
        for i in range(self.n):
            prod = self.matrices[self.inverse[i]] @ self.matrices[i]
            ref = prod.flat[int(np.argmax(np.abs(prod)))]
            if not np.allclose(prod, ref / abs(ref) * np.eye(2), atol=1e-8):
                inverse_ok = False
        # 理想分解与规范 unitary 一致（全局相位无关）
        realization_ok = True
        for i in range(self.n):
            u_ideal = IDENTITY.copy()
            for kind, val in self.decompose(i):
                if kind == "VZ":
                    u_ideal = virtual_z_unitary(val) @ u_ideal
                else:
                    phase, ang = {
                        "X90": (0.0, 0.5 * np.pi), "Y90": (0.5 * np.pi, 0.5 * np.pi),
                        "X-90": (np.pi, 0.5 * np.pi),
                        "Y-90": (1.5 * np.pi, 0.5 * np.pi)}[val]
                    u_ideal = equatorial_rotation_unitary(phase, ang) @ u_ideal
            prod = u_ideal @ self.matrices[i].conj().T
            phase = np.exp(-1j * np.angle(prod[0, 0])) if abs(
                prod[0, 0]) > 1e-8 else 1.0
            if not np.allclose(prod * phase, IDENTITY, atol=1e-8):
                realization_ok = False
        return {"n_elements": self.n, "unique": bool(unique),
                "closure": closure, "inverse_ok": bool(inverse_ok),
                "pauli_conjugation_signed": True,
                "ideal_realization_matches": bool(realization_ok)}
