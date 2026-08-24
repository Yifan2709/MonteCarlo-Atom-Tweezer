"""bare 旋转、SCROFULOUS 复合脉冲与脉冲包络。

SCROFULOUS（Cummins & Jones 2003；论文 Methods 2403.12021v4）实现绕赤道
轴 φ 转 θ 的三段对称复合脉冲 (θc)_{φ1} (π)_{φ2} (θc)_{φ1}：

- sinc 约定：sinc(x)=sin(x)/x（未归一化）；
- θc = arcsinc(2cos(θ/2)/π)，根取 [π/2, π] 单调分支（唯一）；
- φ1 = φ3 = φ + δ，cosδ = −π·cosθc/(2θc·sin(θ/2))；
- φ2 = φ1 − arccos(−π/(2θc))。

对脉冲长度（幅度）误差一阶鲁棒；对失谐不鲁棒。identity 与近零角度
按配置阈值退化为空脉冲列表，不构造病态根。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

TWO_PI = 2.0 * np.pi


def sinc_unnormalized(x):
    """sinc(x)=sin(x)/x（论文 arcsinc 使用的归一化约定）。"""
    x = np.asarray(x, dtype=float)
    out = np.where(np.abs(x) < 1e-12, 1.0 - x * x / 6.0, np.sin(x) / np.where(
        np.abs(x) < 1e-12, 1.0, x))
    return out


def arcsinc_upper(value: float) -> float:
    """arcsinc 于 [π/2, π] 单调分支的根：sin(x)/x = value。

    value 必须落在 (0, 2/π] 内；否则抛出明确错误（分支不存在/不唯一）。
    """
    value = float(value)
    if not (0.0 < value <= 2.0 / np.pi + 1e-12):
        raise ValueError(
            f"arcsinc 输入 {value} 超出 [π/2,π] 分支的值域 (0, 2/π]；"
            "该分支无解，拒绝构造 SCROFULOUS")
    if value > 1.0:
        raise ValueError("arcsinc 输入 >1 无实根")
    return brentq(lambda x: float(sinc_unnormalized(x)) - value,
                  0.5 * np.pi, np.pi, xtol=1e-14, rtol=1e-15)


@dataclass
class PulseSpec:
    """一段物理微波脉冲（复合脉冲的子脉冲也用它表示）。

    duration 为物理时长（秒）；nominal_area 为名义旋转角（rad）；实际面积
    = nominal_area × (1+amp_error) × site_rabi_scale（后者在传播时叠加）。
    """

    axis_phase: float          # φ（rad），驱动相位
    nominal_area: float        # θ（rad）
    duration: float            # s
    start: float = 0.0         # s，排程后生效（列表内相对协议 t0）
    family: str = "bare"
    envelope: str = "square"
    edge_us: float = 0.0
    amp_error: float = 0.0     # 静态分数幅度误差（本脉冲）
    composite_id: str = ""
    clifford_id: int | None = None
    role: str = "gate"

    @property
    def area_pi_units(self) -> float:
        """名义面积（π 为单位）。"""
        return self.nominal_area / np.pi

    def as_row(self) -> dict:
        """CSV 行。"""
        return {"family": self.family, "role": self.role,
                "axis_phase_rad": self.axis_phase,
                "axis_label": phase_axis_label(self.axis_phase),
                "nominal_area_rad": self.nominal_area,
                "area_pi_units": self.area_pi_units,
                "duration_us": self.duration * 1e6,
                "envelope": self.envelope, "edge_us": self.edge_us,
                "amp_error": self.amp_error,
                "composite_id": self.composite_id,
                "clifford_id": self.clifford_id}


def phase_axis_label(phi: float) -> str:
    """把驱动相位映射成 X/Y/X̄/Ȳ 风格标签。"""
    key = int(np.round(np.mod(phi, TWO_PI) / (0.5 * np.pi))) % 4
    return {0: "X", 1: "Y", 2: "Xbar", 3: "Ybar"}[key]


AXIS_PHASE = {"X": 0.0, "Y": 0.5 * np.pi, "Xbar": np.pi, "Ybar": 1.5 * np.pi,
              "Z": None}


def envelope_value(envelope: str, frac: float, edge_us: float, t_us: float,
                   duration_us: float):
    """归一化包络（峰值 1）。frac∈[0,1] 为脉冲内归一化时间。

    square：常数 1；cosine_edge：长度 edge_us 的升/降余弦边沿；
    gaussian：σ 由总时长决定的截断高斯（无多能级模型，仅作形状扫描，
    不声称 DRAG 式泄漏抑制）。
    """
    if envelope == "square":
        return 1.0
    if envelope == "cosine_edge":
        edge = float(edge_us)
        if edge <= 0:
            return 1.0
        up = np.clip(t_us / edge, 0.0, 1.0)
        down = np.clip((duration_us - t_us) / edge, 0.0, 1.0)
        return 0.5 * (1.0 - np.cos(np.pi * up)) * 0.5 * (
            1.0 - np.cos(np.pi * down))
    if envelope == "gaussian":
        sigma = max(duration_us / 6.0, 1e-9)
        return float(np.exp(-0.5 * ((frac - 0.5) * duration_us / sigma) ** 2))
    raise KeyError(f"未知包络 {envelope}")


def envelope_area_integral(envelope: str, duration_s: float, edge_us: float,
                           n_samples: int = 20001) -> float:
    """包络 ∫env dt（用于把名义面积换算成峰值 Rabi）。"""
    if envelope == "square":
        return duration_s
    t = np.linspace(0.0, duration_s * 1e6, n_samples)
    fracs = np.linspace(0.0, 1.0, n_samples)
    vals = np.array([envelope_value(envelope, f, edge_us, tt,
                                    duration_s * 1e6)
                     for f, tt in zip(fracs, t)])
    return float(np.trapezoid(vals, t)) * 1e-6


def make_bare_pulse(axis_phase: float, angle: float, omega_base: float,
                    envelope: str = "square", edge_us: float = 0.0,
                    amp_error: float = 0.0, composite_id: str = "",
                    clifford_id=None, role: str = "gate") -> PulseSpec:
    """bare 有限时长旋转脉冲：duration=θ/Ω0；非方包络在传播时按面积补偿峰值。"""
    if angle <= 0:
        raise ValueError("bare 脉冲角度必须为正（identity 用空列表表示）")
    return PulseSpec(axis_phase=axis_phase, nominal_area=angle,
                     duration=float(angle / omega_base), family="bare",
                     envelope=envelope, edge_us=edge_us, amp_error=amp_error,
                     composite_id=composite_id, clifford_id=clifford_id,
                     role=role)


@dataclass(frozen=True)
class ScrofulousAngles:
    """SCROFULOUS 的三段角度/相位。"""

    theta_c: float
    delta_phase: float
    gamma_phase: float
    theta: float
    axis_phase: float

    @property
    def pulse_phases(self):
        """(φ1, φ2, φ1)（时间顺序）。"""
        phi1 = self.axis_phase + self.delta_phase
        return phi1, phi1 + self.gamma_phase, phi1

    @property
    def total_area(self) -> float:
        """总名义脉冲面积（rad）。"""
        return 2.0 * self.theta_c + np.pi

    def describe(self) -> dict:
        """参数描述（含分支与约定标注）。"""
        phi1, phi2, _ = self.pulse_phases
        return {"theta_target_rad": self.theta,
                "theta_c_rad": self.theta_c,
                "theta_c_pi_units": self.theta_c / np.pi,
                "sinc_convention": "sinc(x)=sin(x)/x (unnormalized)",
                "root_branch": "[pi/2, pi] monotone, unique",
                "phi1_rad": phi1, "phi2_rad": phi2,
                "delta_rad": self.delta_phase,
                "gamma_rad": self.gamma_phase,
                "total_area_rad": self.total_area,
                "total_area_pi_units": self.total_area / np.pi}


def scrofulous_angles(theta: float, axis_phase: float = 0.0) -> ScrofulousAngles:
    """求解 SCROFULOUS 角度/相位（θ∈(min_angle, π]，近零见 min_angle_theta）。"""
    theta = float(theta)
    if not (0.0 < theta <= np.pi + 1e-12):
        raise ValueError(f"SCROFULOUS 目标角 θ={theta} 超出 (0, π]")
    theta = min(theta, np.pi)
    theta_c = arcsinc_upper(2.0 * np.cos(0.5 * theta) / np.pi)
    cos_delta = -np.pi * np.cos(theta_c) / (2.0 * theta_c * np.sin(0.5 * theta))
    if abs(cos_delta) > 1.0 + 1e-9:
        raise ValueError(
            f"SCROFULOUS δ 分支不存在：cosδ={cos_delta}；检查 θ 与 sinc 根")
    cos_delta = float(np.clip(cos_delta, -1.0, 1.0))
    delta_phase = float(np.arccos(cos_delta))
    gamma_phase = -float(np.arccos(float(np.clip(-np.pi / (2.0 * theta_c),
                                                 -1.0, 1.0))))
    return ScrofulousAngles(theta_c=theta_c, delta_phase=delta_phase,
                            gamma_phase=gamma_phase, theta=theta,
                            axis_phase=axis_phase)


SCROFULOUS_MIN_THETA = 1e-3  # rad；小于此角度退化为 bare（避免病态相位）


def make_scrofulous_pulses(axis_phase: float, angle: float, omega_base: float,
                           envelope: str = "square", edge_us: float = 0.0,
                           amp_error: float = 0.0, composite_id: str = "",
                           clifford_id=None, role: str = "gate"
                           ) -> list[PulseSpec]:
    """SCROFULOUS 三段脉冲列表（时间顺序）；近零角度退化为 bare。

    amp_error 作用于全部三段（共同的脉冲长度误差模型）。
    """
    angle = abs(float(angle))
    if angle < SCROFULOUS_MIN_THETA:
        if angle < 1e-12:
            return []
        return [make_bare_pulse(axis_phase, angle, omega_base, envelope,
                                edge_us, amp_error, composite_id,
                                clifford_id, role)]
    ang = scrofulous_angles(angle, axis_phase)
    phi1, phi2, _ = ang.pulse_phases
    dur_c = ang.theta_c / omega_base
    dur_pi = np.pi / omega_base
    return [
        PulseSpec(axis_phase=float(phi1), nominal_area=float(ang.theta_c),
                  duration=float(dur_c), family="scrofulous",
                  envelope=envelope, edge_us=edge_us, amp_error=amp_error,
                  composite_id=composite_id, clifford_id=clifford_id,
                  role=role),
        PulseSpec(axis_phase=float(phi2), nominal_area=float(np.pi),
                  duration=float(dur_pi), family="scrofulous",
                  envelope=envelope, edge_us=edge_us, amp_error=amp_error,
                  composite_id=composite_id, clifford_id=clifford_id,
                  role=role),
        PulseSpec(axis_phase=float(phi1), nominal_area=float(ang.theta_c),
                  duration=float(dur_c), family="scrofulous",
                  envelope=envelope, edge_us=edge_us, amp_error=amp_error,
                  composite_id=composite_id, clifford_id=clifford_id,
                  role=role),
    ]


@dataclass
class VirtualZ:
    """虚拟 Z 旋转（瞬时帧更新，零时长、零误差）。"""

    angle: float
    clifford_id: int | None = None
    role: str = "frame"

    def as_row(self) -> dict:
        """CSV 行。"""
        return {"family": "virtual_z", "role": self.role,
                "axis_phase_rad": "", "axis_label": "Z",
                "nominal_area_rad": self.angle,
                "area_pi_units": self.angle / np.pi, "duration_us": 0.0,
                "envelope": "none", "edge_us": 0.0, "amp_error": 0.0,
                "composite_id": "", "clifford_id": self.clifford_id}


def pulse_list_for(axis_phase: float, angle: float, family: str,
                   omega_base: float, envelope: str = "square",
                   edge_us: float = 0.0, amp_error: float = 0.0,
                   composite_id: str = "", clifford_id=None,
                   role: str = "gate") -> list[PulseSpec]:
    """按 family 构造一段旋转的物理脉冲列表（bare/scrofulous）。"""
    if abs(angle) < 1e-12:
        return []
    if family == "bare":
        phase = axis_phase + (np.pi if angle < 0 else 0.0)
        return [make_bare_pulse(phase, abs(angle), omega_base, envelope,
                                edge_us, amp_error, composite_id,
                                clifford_id, role)]
    if family == "scrofulous":
        phase = axis_phase + (np.pi if angle < 0 else 0.0)
        return make_scrofulous_pulses(phase, abs(angle), omega_base, envelope,
                                      edge_us, amp_error, composite_id,
                                      clifford_id, role)
    raise KeyError(f"未知 pulse family {family}")


def pulse_list_duration(pulses: list) -> float:
    """脉冲列表总物理时长（虚拟 Z 不计）。"""
    return float(sum(p.duration for p in pulses
                     if isinstance(p, PulseSpec)))


def pulse_list_area(pulses: list) -> float:
    """脉冲列表总名义面积（rad）。"""
    return float(sum(p.nominal_area for p in pulses
                     if isinstance(p, PulseSpec)))
