"""三类光镊噪声的平均强度（确定性）加热模型。

把噪声的系综平均加热功率作为确定性功率注入轨迹：不含任何随机性，
默认实现为常数，注入能量随时间单调累积，系统随时间延长越发不稳定。

- 光子反冲: P_rec = 2·E_r·Γ_sc —— 常数功率，能量线性增长（扩散型）；
- 强度噪声（参量加热）: dE/dt = Γ_par·ε_osc，
  Γ_par = (π²/4)·ν₀²·S_RIN(2ν₀) —— 功率正比于当前振荡能量，能量指数增长；
- 指向/位置噪声: P_pt = m·ω₀⁴·S_x(ν₀)/8 —— 常数功率，线性增长。
  PSD 均为单边每 Hz 约定；纯音锚点：ε_I = ε₀cos(2ω₀t) 时严格解
  E(t) = E₀·exp(ε₀ω₀t/2)，可用于锁定系数约定。

三个噪声均以函数（callable）形式外放为输入参数，签名统一为
``(x_m, v_m_s, t_s, e_osc_J) -> 加热功率 [W]``；默认为物理常数构造的
常数函数（functools.partial，可 pickle 以支持多进程），调用方可注入
任意自定义函数（如时变平均强度；lambda 仅限单进程路径）。

所有噪声幅度参数均为 assumed_sensitivity_only，不声称来自论文实验数据。
"""
from __future__ import annotations

import math
from functools import partial

import numpy as np
from scipy.constants import hbar, pi, speed_of_light

from level0_static_trap.constants import microkelvin_to_joule
from level0_static_trap.potentials import trap_frequency_analytic

# Cs D2 线常数（NIST）：波长与自然线宽（频率，与失谐同单位代入比值）。
CS_D2_WAVELENGTH_M = 852.347e-9
CS_D2_GAMMA_HZ = 5.234e6


def constant_power(x_m, v_m_s, t_s, e_osc_J, power_w):
    """常数加热功率（光子反冲与指向噪声通道的默认实现）。"""
    return power_w


def parametric_rate_power(x_m, v_m_s, t_s, e_osc_J, rate_per_s):
    """参量加热默认实现：功率 = Γ_par × 当前振荡能量 ε_osc。"""
    return rate_per_s * max(0.0, float(e_osc_J))


def recoil_energy_j(trap_wavelength_m: float, mass_kg: float) -> float:
    """光子反冲能量 E_r = ħ²(2π/λ)²/(2m)。"""
    wave_number = 2.0 * pi / trap_wavelength_m
    return hbar ** 2 * wave_number ** 2 / (2.0 * mass_kg)


def recoil_scattering_rate_cs(trap_wavelength_m: float, depth_j: float,
                              d2_wavelength_m: float = CS_D2_WAVELENGTH_M,
                              d2_gamma_hz: float = CS_D2_GAMMA_HZ) -> float:
    """两能级近似散射率 Γ_sc = (Γ/|Δ|)·(D/ħ)，Δ 相对 Cs D2 线。"""
    detuning_hz = abs(speed_of_light * (1.0 / d2_wavelength_m - 1.0 / trap_wavelength_m))
    return d2_gamma_hz / detuning_hz * depth_j / hbar


def parametric_rate_from_rin(nu0_hz: float, rin_psd_per_hz: float) -> float:
    """参量加热指数增长率 Γ_par = (π²/4)·ν₀²·S_RIN(2ν₀)。

    约定：S_RIN 为单边 PSD（每 Hz）。纯音自检：ε_I = ε₀cos(2ω₀t) 时
    严格解 E(t) = E₀·exp(ε₀ω₀t/2)，可用于锁定耦合系数约定。
    """
    return 0.25 * pi ** 2 * nu0_hz ** 2 * rin_psd_per_hz


def pointing_power_from_psd(mass_kg: float, nu0_hz: float,
                            position_psd_m2_per_hz: float) -> float:
    """指向噪声恒定加热功率 P = m·ω₀⁴·S_x(ν₀)/8，S_x 为单边每 Hz PSD。

    约定经 Green 函数法推导（速度冲击响应核 + 双边/单边谱换算）；
    共振处几十 pm/√Hz 即可在亚毫秒内显著加热，nm/√Hz 对应 K/s 量级。
    """
    omega0 = 2.0 * pi * nu0_hz
    return mass_kg * omega0 ** 4 * position_psd_m2_per_hz / 8.0


class MeanHeating:
    """三个可注入噪声函数的容器；通道为 None 表示关闭。

    channel_powers 返回各通道功率（W）；e_osc_J 由调用方（积分器）提供，
    仅参量通道使用。
    """

    def __init__(self, recoil=None, parametric=None, pointing=None, labels=None):
        self.recoil = recoil
        self.parametric = parametric
        self.pointing = pointing
        self.labels = dict(labels or {})

    def channel_powers(self, x_m: float, v_m_s: float, t_s: float,
                       e_osc_J: float) -> dict:
        """返回 {recoil, parametric, pointing, total} 各通道加热功率（W）。"""
        powers = {
            "recoil": float(self.recoil(x_m, v_m_s, t_s, e_osc_J)) if self.recoil else 0.0,
            "parametric": float(self.parametric(x_m, v_m_s, t_s, e_osc_J)) if self.parametric else 0.0,
            "pointing": float(self.pointing(x_m, v_m_s, t_s, e_osc_J)) if self.pointing else 0.0,
        }
        powers["total"] = powers["recoil"] + powers["parametric"] + powers["pointing"]
        return powers

    def total_power(self, x_m: float, v_m_s: float, t_s: float, e_osc_J: float) -> float:
        """三通道总加热功率（W）。"""
        return self.channel_powers(x_m, v_m_s, t_s, e_osc_J)["total"]

    def describe(self) -> dict:
        """噪声参数描述（全部标注 mean_rate_deterministic / assumed）。"""
        return {
            "model": "mean_rate_deterministic",
            "channels": {
                "recoil": self.labels.get("recoil"),
                "parametric": self.labels.get("parametric"),
                "pointing": self.labels.get("pointing"),
            },
            "provenance": "assumed_sensitivity_only",
            "randomness": False,
        }


def build_mean_heating(recoil_power_w=None, parametric_rate_s=None,
                       pointing_power_w=None) -> MeanHeating:
    """由常数参数构造 MeanHeating；None 通道关闭。标签记录 SI 值与单位换算。"""
    labels = {}
    if recoil_power_w is not None:
        labels["recoil"] = {"kind": "constant_power", "power_w": float(recoil_power_w),
                            "power_uK_per_s": float(recoil_power_w) / microkelvin_to_joule(1.0)}
    if parametric_rate_s is not None:
        labels["parametric"] = {"kind": "rate_times_e_osc", "rate_per_s": float(parametric_rate_s)}
    if pointing_power_w is not None:
        labels["pointing"] = {"kind": "constant_power", "power_w": float(pointing_power_w),
                              "power_uK_per_s": float(pointing_power_w) / microkelvin_to_joule(1.0)}
    return MeanHeating(
        recoil=None if recoil_power_w is None else partial(constant_power, power_w=float(recoil_power_w)),
        parametric=None if parametric_rate_s is None else partial(parametric_rate_power,
                                                                  rate_per_s=float(parametric_rate_s)),
        pointing=None if pointing_power_w is None else partial(constant_power, power_w=float(pointing_power_w)),
        labels=labels)


def heating_from_config(cfg) -> MeanHeating | None:
    """由 Level2Config 的 noise_mean_heating 段构造模型；关闭时返回 None。

    未显式给出速率时按物理常数计算：反冲取标称 AOD 深度的 Γ_sc；
    参量/指向以 SLM 解析阱频为参考频率（可用 reference_nu0_hz 覆盖）。
    """
    noise = cfg.noise_mean_heating
    if not noise.enabled:
        return None
    mass_kg = cfg.atom.mass_kg
    trap_wavelength_m = noise.trap_wavelength_nm * 1e-9

    if noise.recoil_scattering_rate_s is None:
        scattering_rate = recoil_scattering_rate_cs(trap_wavelength_m, cfg.aod.initial_depth_j)
    else:
        scattering_rate = float(noise.recoil_scattering_rate_s)
    recoil_power = 2.0 * recoil_energy_j(trap_wavelength_m, mass_kg) * scattering_rate

    if noise.reference_nu0_hz is None:
        nu0_hz = trap_frequency_analytic(cfg.slm.depth_j, cfg.slm.waist_m, mass_kg)[1]
    else:
        nu0_hz = float(noise.reference_nu0_hz)

    if noise.parametric_rate_s is None:
        parametric_rate = parametric_rate_from_rin(nu0_hz, noise.parametric_rin_psd_per_hz)
    else:
        parametric_rate = float(noise.parametric_rate_s)

    if noise.pointing_heating_rate_uK_per_s is None:
        pointing_power = pointing_power_from_psd(mass_kg, nu0_hz, noise.pointing_position_psd_m2_per_hz)
    else:
        pointing_power = microkelvin_to_joule(float(noise.pointing_heating_rate_uK_per_s))

    heating = build_mean_heating(recoil_power_w=recoil_power,
                                 parametric_rate_s=parametric_rate,
                                 pointing_power_w=pointing_power)
    heating.labels["derived_from"] = {
        "trap_wavelength_nm": float(noise.trap_wavelength_nm),
        "recoil_scattering_rate_s": float(scattering_rate),
        "recoil_energy_kB_nK": float(recoil_energy_j(trap_wavelength_m, mass_kg)
                                     / microkelvin_to_joule(1.0) * 1e3),
        "reference_nu0_hz": float(nu0_hz),
        "parametric_rin_psd_per_hz": float(noise.parametric_rin_psd_per_hz),
        "pointing_position_psd_m2_per_hz": float(noise.pointing_position_psd_m2_per_hz),
    }
    return heating


def _scaled_call(x_m, v_m_s, t_s, e_osc_J, inner, factor):
    """内部通道函数的整体缩放（保持可 pickle 以支持多进程）。"""
    return inner(x_m, v_m_s, t_s, e_osc_J) * factor


def scaled_heating(base: MeanHeating, factor: float) -> MeanHeating:
    """把基础模型三个通道的功率/速率整体放大 factor 倍（敏感性扫描用）。"""
    def scale(fn):
        if fn is None:
            return None
        return partial(_scaled_call, inner=fn, factor=float(factor))

    scaled = MeanHeating(recoil=scale(base.recoil), parametric=scale(base.parametric),
                         pointing=scale(base.pointing), labels=dict(base.labels))
    scaled.labels["scale_factor"] = float(factor)
    return scaled


def velocity_verlet_time_dependent_heating(force_function, mass, x_initial, v_initial,
                                           times, heating: MeanHeating, e_osc_fn):
    """带确定性平均加热的 velocity-Verlet（其余与 Level 1 积分器一致）。

    每步速度更新后施加加热踢：ΔE = P·dt，沿当前速度方向缩放
    ``v ← sign(v)·sqrt(v² + 2ΔE/m)``（v=0 取正方向）；逐通道累计注入能量。
    返回 (time, x, v, a, noise)，noise 为各通道与总的累积能量数组（J）。
    """
    time = np.asarray(times, dtype=float)
    if time.ndim != 1 or time.size < 2 or not np.all(np.isfinite(time)):
        raise ValueError("invalid times")
    steps = np.diff(time)
    if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-11,
                                             atol=max(1e-18, steps[0] * 1e-13)):
        raise ValueError("times must be uniformly increasing")
    if mass <= 0 or not np.isfinite(mass):
        raise ValueError("invalid mass")
    dt = float(steps[0])
    x = np.empty_like(time)
    v = np.empty_like(time)
    a = np.empty_like(time)
    work = {key: np.zeros_like(time) for key in ("recoil_J", "parametric_J",
                                                 "pointing_J", "total_J")}
    x[0], v[0] = x_initial, v_initial
    a[0] = float(force_function(x[0], time[0])) / mass
    for i in range(time.size - 1):
        x[i + 1] = x[i] + v[i] * dt + 0.5 * a[i] * dt ** 2
        a[i + 1] = float(force_function(x[i + 1], time[i + 1])) / mass
        v[i + 1] = v[i] + 0.5 * (a[i] + a[i + 1]) * dt
        e_osc = float(e_osc_fn(x[i + 1], v[i + 1], time[i + 1]))
        powers = heating.channel_powers(float(x[i + 1]), float(v[i + 1]),
                                        float(time[i + 1]), e_osc)
        delta = powers["total"] * dt
        if delta != 0.0:
            speed_squared = max(0.0, v[i + 1] ** 2 + 2.0 * delta / mass)
            v[i + 1] = math.copysign(math.sqrt(speed_squared), v[i + 1])
        work["recoil_J"][i + 1] = work["recoil_J"][i] + powers["recoil"] * dt
        work["parametric_J"][i + 1] = work["parametric_J"][i] + powers["parametric"] * dt
        work["pointing_J"][i + 1] = work["pointing_J"][i] + powers["pointing"] * dt
        work["total_J"][i + 1] = work["total_J"][i] + delta
    return time, x, v, a, work
