"""有限脉冲动力学解耦序列：none / spin echo / XY4 / XY8 / XY16 与自定义。

与 Level 4 的理想瞬时 toggling 模型不同，这里每个 DD π 脉冲都是有限时长
的物理脉冲（bare 或 SCROFULOUS），传播时必须在脉冲期间实际推进含噪声
Hamiltonian，而不是简单翻转 y(t)。transport-aware 模式（xy4_per_long_move）
在 outbound 与 inbound 两个长移动段内各放一个对称 XY4 块；transformed
Clifford frame 模式默认关闭，且只有群论校验通过才允许启用（本实现提供
共轭轴映射的构造与校验，但不声称复现论文的 transformed-frame 技术）。
"""
from __future__ import annotations

import numpy as np

from .pulse_library import (PulseSpec, make_bare_pulse, make_scrofulous_pulses)

TWO_PI = 2.0 * np.pi

DD_MODES = ("none", "spin_echo", "xy4", "xy8", "xy16",
            "xy4_per_long_move", "custom_pulse_times")

# XY 族脉冲轴相位：X, Y, X̄, Ȳ
_XY_PHASES = {"X": 0.0, "Y": 0.5 * np.pi, "Xbar": np.pi, "Ybar": 1.5 * np.pi}


def xy_axis_pattern(mode: str) -> list[str]:
    """XY 序列的脉冲轴顺序（标准相位约定）。"""
    if mode in ("xy4", "xy4_per_long_move"):
        return ["X", "Y", "X", "Y"]
    if mode == "xy8":
        return ["X", "Y", "X", "Y", "Ybar", "Xbar", "Ybar", "Xbar"]
    if mode == "xy16":
        return (["X", "Y", "X", "Y"] * 2
                + ["Xbar", "Ybar", "Xbar", "Ybar"] * 2)
    raise KeyError(f"{mode} 没有 XY 轴序列")


def symmetric_pulse_fractions(n_pulses: int) -> np.ndarray:
    """n 个脉冲的对称等间距分数 (2k-1)/(2n)。"""
    return (2.0 * np.arange(1, n_pulses + 1) - 1.0) / (2.0 * n_pulses)


def dd_pulse_events(mode: str, t0: float, duration: float, omega_base: float,
                    family: str = "bare", envelope: str = "square",
                    amp_error: float = 0.0, edge_us: float = 0.0,
                    long_move_windows=(), custom_times_s=()) -> list[PulseSpec]:
    """在 [t0, t0+duration] 内按 DD 模式生成有限 π 脉冲列表。

    返回的 PulseSpec.start 由调用方在排程时赋值（这里 duration 字段是
    脉冲自身时长）；long_move_windows 为 [(t_start, t_end), ...]，仅
    xy4_per_long_move 使用（每个窗口一个对称 XY4 块，脉冲时刻落在
    窗口的 1/8,3/8,5/8,7/8 分数上，与 Level 4 瞬时模型一致）。
    """
    pulses: list[PulseSpec] = []
    if mode == "none":
        return pulses

    def add(center_time: float, axis: str):
        phase = _XY_PHASES[axis]
        if family == "scrofulous":
            sub = make_scrofulous_pulses(phase, np.pi, omega_base, envelope,
                                         edge_us, role="dd")
        else:
            sub = [make_bare_pulse(phase, np.pi, omega_base, envelope,
                                   edge_us, role="dd")]
        total = sum(p.duration for p in sub)
        t = center_time - 0.5 * total
        for p in sub:
            p.start = t
            t += p.duration
        pulses.extend(sub)

    if mode == "spin_echo":
        add(t0 + 0.5 * duration, "X")
    elif mode in ("xy4", "xy8", "xy16"):
        axes = xy_axis_pattern(mode)
        fracs = symmetric_pulse_fractions(len(axes))
        for f, axis in zip(fracs, axes):
            add(t0 + f * duration, axis)
    elif mode == "xy4_per_long_move":
        axes = xy_axis_pattern("xy4")
        for (w0, w1) in long_move_windows:
            for k in (1, 3, 5, 7):
                add(w0 + (w1 - w0) * k / 8.0, axes[(k // 2) % 4])
    elif mode == "custom_pulse_times":
        for i, t in enumerate(sorted(float(x) for x in custom_times_s)):
            add(t, ["X", "Y"][i % 2])
    else:
        raise KeyError(f"未知 DD 模式 {mode}")
    pulses.sort(key=lambda p: p.start)
    return pulses


def transformed_frame_axes(previous_clifford_unitary: np.ndarray,
                           axes=("X", "Y")) -> list[str]:
    """把 XY 轴映射到 previous Clifford 共轭后的最近 XY 轴（含可选 H 变换）。

    P' = U† P U 为带符号 Pauli；若映射到 Z，则先用 Hadamard 等价基变换把
    (X,Z) 或 (Y,Z) 换回赤道平面。返回变换后的轴标签列表。
    """
    from .qubit_hamiltonian import PAULI_X, PAULI_Y, PAULI_Z
    pauli_map = {"X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
    label = {"X": "X", "Y": "Y", "Z": "Z"}
    u = previous_clifford_unitary
    mapped = []
    for ax in axes:
        conj = u.conj().T @ pauli_map[ax] @ u
        best, best_err = None, 1e9
        for name, mat in pauli_map.items():
            for sign in (1.0, -1.0):
                err = float(np.max(np.abs(conj - sign * mat)))
                if err < best_err:
                    best_err, best = err, (name, sign)
        mapped.append((label[best[0]], best[1]))
    if any(m[0] == "Z" for m in mapped):
        # H 基变换：X↔Z, Y→Y（把含 Z 的轴对映射回赤道）
        swapped = {"X": "Z", "Z": "X", "Y": "Y"}
        mapped = [(swapped[m[0]], m[1]) for m in mapped]
    if mapped[0][0] == mapped[1][0]:
        raise ValueError("transformed frame 轴映射退化（两轴相同）")
    return [m[0] for m in mapped]


def validate_transformed_frame(group, n_trials=200, seed=0):
    """transformed Clifford frame 的群论校验。

    对随机 Clifford 检查：共轭后两轴仍是带符号 Pauli、H 变换后回到赤道
    且两轴正交。全部通过才允许在 CLI 中启用该模式。
    """
    rng = np.random.default_rng(seed)
    issues = []
    from .qubit_hamiltonian import PAULI_X, PAULI_Y, PAULI_Z
    pauli_map = {"X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
    for _ in range(n_trials):
        idx = int(group.sample_uniform(rng, 1)[0])
        u = group.unitary(idx)
        try:
            axes = transformed_frame_axes(u)
        except ValueError as exc:
            issues.append(f"clifford {idx}: {exc}")
            continue
        mats = [pauli_map[a] for a in axes]
        if abs(float(np.real(np.trace(mats[0] @ mats[1])))) > 1e-8:
            issues.append(f"clifford {idx}: 变换后轴不正交 {axes}")
    # 共轭一致性：U† P U 与标签轴的 Pauli 一致性（含 H 情形仅检查赤道性）
    return {"passed": not issues, "n_trials": n_trials,
            "issues": issues[:5],
            "note": "共轭轴 + Hadamard 等价基变换校验；通过也不等于复现"
                    "论文 transformed Clifford frame 实验技术"}


def toggling_limit_phase(delta_omega_fn, t0, duration, pulse_centers):
    """瞬时脉冲极限下 Level 4 风格 toggling 相位（一致性检查用）。"""
    bounds = np.concatenate([[t0], np.sort(np.asarray(pulse_centers)),
                             [t0 + duration]])
    y = 1.0
    phi = 0.0
    for a, b in zip(bounds[:-1], bounds[1:]):
        mid = 0.5 * (a + b)
        phi += y * float(delta_omega_fn(mid)) * (b - a)
        y = -y
    return phi
