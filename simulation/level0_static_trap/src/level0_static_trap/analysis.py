"""数值分析：平衡位置、数值曲率、轨迹频率估计与非线性周期积分。

误差被严格区分为：
- 数值离散化误差（数值曲率频率 vs 解析小振幅频率）
- 有限振幅非简谐频移（小振幅频率 vs 非线性积分频率）
- 轨迹积分误差（轨迹频率 vs 同能量非线性积分频率）
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import integrate, optimize

from .constants import j_to_uk
from .potentials import (
    gaussian_curvature_analytic,
    gaussian_force,
    gaussian_potential,
    trap_frequency_analytic,
)

ForceFunction = Callable[[float | np.ndarray], float | np.ndarray]
PotentialFunction = Callable[[float | np.ndarray], float | np.ndarray]


# ---------------------------------------------------------------------------
# 相对误差（显式处理零参考量）
# ---------------------------------------------------------------------------

def calculate_relative_error(
    value: float, reference: float, zero_reference_action: str = "return_nan"
) -> float:
    """计算 ``|value - reference| / |reference|``。

    当 ``reference == 0`` 时按 ``zero_reference_action`` 处理：
    - ``"return_nan"``：返回 ``nan``（默认）；
    - ``"raise"``：抛出 ValueError；
    - ``"return_signed_inf"``：当 value≠0 时返回 ``inf``，否则返回 0。
    """
    if reference == 0.0:
        if zero_reference_action == "raise":
            raise ValueError("参考量为零，无法计算相对误差")
        if zero_reference_action == "return_signed_inf":
            return math.inf if value != 0.0 else 0.0
        return float("nan")
    return abs(value - reference) / abs(reference)


# ---------------------------------------------------------------------------
# 平衡位置
# ---------------------------------------------------------------------------

@dataclass
class EquilibriumResult:
    """平衡位置搜索结果。"""

    x_eq: float            # 平衡位置 (m)
    u_eq: float            # 平衡位置处势能 (J)
    force_residual: float  # 平衡位置处力的绝对残差 (N)


def find_equilibrium_position(
    potential_fn: PotentialFunction,
    force_fn: ForceFunction,
    center_guess: float,
    waist: float,
    bracket_waist: float = 0.5,
) -> EquilibriumResult:
    """寻找势能极小（力为零）的平衡位置。

    采用有界、无量纲化的策略：以 ``center_guess`` 为中心，在
    ``[center_guess - bracket_waist*w, center_guess + bracket_waist*w]`` 区间内
    对力做有根搜索（``brentq``，要求两侧力符号相反）。对于这里使用的中心对称
    高斯势，势阱中心就是极小，力在两侧符号相反，符号变号法稳定可靠。

    返回平衡位置、势能与力残差。
    """
    if waist <= 0.0:
        raise ValueError(f"waist 必须为正，实际值：{waist:g} m")
    a = center_guess - bracket_waist * waist
    b = center_guess + bracket_waist * waist
    fa = float(np.asarray(force_fn(a)).item())
    fb = float(np.asarray(force_fn(b)).item())
    # brentq 需要符号相反。对中心对称高斯阱 a<xc<b 处处成立：左侧力>0，右侧力<0。
    if fa == 0.0:
        x_eq = a
    elif fb == 0.0:
        x_eq = b
    else:
        if fa * fb > 0.0:
            # 退化情况：用标量极小化势能作为兜底（不依赖默认绝对容差量级）。
            scale = waist
            res = optimize.minimize_scalar(
                lambda s: float(np.asarray(potential_fn(center_guess + s * scale)).item()),
                bracket=(-bracket_waist, 0.0, bracket_waist),
                method="brent",
                options={"xtol": 1.0e-12},
            )
            x_eq = center_guess + res.x * scale
        else:
            x_eq = optimize.brentq(
                lambda xx: float(np.asarray(force_fn(xx)).item()),
                a,
                b,
                xtol=1.0e-15 * waist,
                rtol=4.0 * np.finfo(float).eps,
                maxiter=200,
            )
    u_eq = float(np.asarray(potential_fn(x_eq)).item())
    f_res = float(np.asarray(force_fn(x_eq)).item())
    return EquilibriumResult(x_eq=x_eq, u_eq=u_eq, force_residual=abs(f_res))


# ---------------------------------------------------------------------------
# 数值曲率（差分，h 与 h/2）
# ---------------------------------------------------------------------------

@dataclass
class CurvatureResult:
    """数值曲率结果。"""

    analytic: float       # 解析曲率 (N/m)
    numerical_h: float    # 步长 h 的数值曲率 (N/m)
    numerical_h2: float   # 步长 h/2 的数值曲率 (N/m)
    convergence_diff: float  # |kappa_h - kappa_{h/2}| (N/m)
    step_h: float         # 实际使用的 h (m)


def numerical_curvature(
    potential_fn: PotentialFunction,
    x0: float,
    h: float,
) -> float:
    """中心二阶差分估计 ``d^2U/dx^2`` 在 ``x0`` 处的值 (N/m)。

    ``kappa = (U(x0+h) - 2 U(x0) + U(x0-h)) / h^2``。
    """
    if h <= 0.0:
        raise ValueError(f"差分步长 h 必须为正，实际值：{h:g}")
    up = float(np.asarray(potential_fn(x0 + h)).item())
    mid = float(np.asarray(potential_fn(x0)).item())
    dn = float(np.asarray(potential_fn(x0 - h)).item())
    return (up - 2.0 * mid + dn) / (h * h)


def numerical_curvature_pair(
    potential_fn: PotentialFunction,
    depth: float,
    waist: float,
    center: float,
    step_waist: float,
) -> CurvatureResult:
    """在势阱中心用 ``h`` 与 ``h/2`` 两个步长计算数值曲率，并报告收敛差异。

    ``step_waist`` 是以 ``waist`` 为单位的无量纲步长；实际 h = step_waist * waist。
    """
    if waist <= 0.0:
        raise ValueError(f"waist 必须为正，实际值：{waist:g} m")
    if step_waist <= 0.0:
        raise ValueError(f"step_waist 必须为正，实际值：{step_waist:g}")
    h = step_waist * waist
    k_h = numerical_curvature(potential_fn, center, h)
    k_h2 = numerical_curvature(potential_fn, center, 0.5 * h)
    k_analytic = gaussian_curvature_analytic(depth, waist)
    return CurvatureResult(
        analytic=k_analytic,
        numerical_h=k_h,
        numerical_h2=k_h2,
        convergence_diff=abs(k_h - k_h2),
        step_h=h,
    )


def numerical_trap_frequency(kappa: float, mass: float) -> float:
    """由曲率 ``kappa`` (N/m) 与质量 ``mass`` (kg) 计算普通频率 (Hz)。"""
    if kappa <= 0.0:
        raise ValueError(f"曲率必须为正，实际值：{kappa:g} N/m")
    if mass <= 0.0:
        raise ValueError(f"质量必须为正，实际值：{mass:g} kg")
    return math.sqrt(kappa / mass) / (2.0 * math.pi)


# ---------------------------------------------------------------------------
# 轨迹频率估计（带插值的过零点）
# ---------------------------------------------------------------------------

@dataclass
class TrajectoryFrequencyResult:
    """轨迹频率估计结果。"""

    frequency_hz: float           # 估计频率 (Hz)
    method: str
    num_zero_crossings: int       # 使用的过零点数
    complete_cycles_estimate: int # 估计的完整周期数
    fft_resolution_hz: float      # FFT 原始分辨率 1/T (Hz)
    fft_peak_freq_hz: float       # FFT 最高 bin 频率 (Hz，仅诊断)
    duration_s: float
    cycle_count_sufficient: bool


def estimate_frequency_from_trajectory(
    times: np.ndarray,
    displacement: np.ndarray,
    method: str = "interpolated_zero_crossings",
    minimum_complete_cycles: int = 20,
) -> TrajectoryFrequencyResult:
    """从位移轨迹估计振动频率。

    方法：以位移减去其均值后的过零点为基础。每两个连续同向（上升沿）过零点
    对应一个完整周期；用线性插值精确定位过零时刻，取所有相邻上升沿间距的
    平均作为周期估计，频率 = 1 / <period>。同时输出 FFT 诊断：最高 bin 频率
    与原始分辨率 ``1/T``（不得仅按 FFT bin 声称高精度）。

    ``displacement`` 应当是相对势阱中心的位移（轨迹分析层提供）。
    """
    times = np.asarray(times, dtype=float)
    disp = np.asarray(displacement, dtype=float)
    if times.shape != disp.shape or times.ndim != 1:
        raise ValueError("times 与 displacement 必须是相同形状的一维数组")
    if times.size < 4:
        raise ValueError("轨迹过短，无法估计频率")
    d = disp - np.mean(disp)
    # 上升沿过零点：sign 从 <=0 变 >0
    sign = np.sign(d)
    # 把恰好为零的点当作正
    sign[sign == 0] = 1.0
    crossings = []
    for i in range(len(d) - 1):
        s0 = sign[i]
        s1 = sign[i + 1]
        if s0 <= 0.0 and s1 > 0.0:
            # 线性插值过零时刻
            y0 = d[i]
            y1 = d[i + 1]
            frac = -y0 / (y1 - y0)
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            crossings.append(t_cross)
    duration_s = float(times[-1] - times[0])
    fft_res = float("nan") if duration_s <= 0 else 1.0 / duration_s

    # FFT 诊断（不去除均值，找最高非零 bin）
    n = d.size
    fft_freqs = np.fft.rfftfreq(n, d=(times[1] - times[0]))
    fft_amp = np.abs(np.fft.rfft(d))
    if fft_amp.size > 1:
        peak_idx = int(np.argmax(fft_amp[1:]) + 1)
        fft_peak = float(fft_freqs[peak_idx])
    else:
        fft_peak = 0.0

    if method != "interpolated_zero_crossings":
        raise ValueError(f"未知的频率估计方法：{method!r}")

    if len(crossings) < 2:
        # 过零点不足：退回到 FFT 诊断频率，但标记不足。
        freq = fft_peak if fft_peak > 0 else float("nan")
        return TrajectoryFrequencyResult(
            frequency_hz=freq,
            method=method + "(fallback_fft)",
            num_zero_crossings=len(crossings),
            complete_cycles_estimate=max(len(crossings) - 1, 0),
            fft_resolution_hz=fft_res,
            fft_peak_freq_hz=fft_peak,
            duration_s=duration_s,
            cycle_count_sufficient=False,
        )

    rising = np.array(crossings, dtype=float)
    periods = np.diff(rising)
    period_mean = float(np.mean(periods))
    freq = 1.0 / period_mean if period_mean > 0 else float("nan")
    complete_cycles = len(crossings) - 1
    return TrajectoryFrequencyResult(
        frequency_hz=freq,
        method=method,
        num_zero_crossings=len(crossings),
        complete_cycles_estimate=complete_cycles,
        fft_resolution_hz=fft_res,
        fft_peak_freq_hz=fft_peak,
        duration_s=duration_s,
        cycle_count_sufficient=complete_cycles >= minimum_complete_cycles,
    )


# ---------------------------------------------------------------------------
# 非线性周期（一维能量积分）
# ---------------------------------------------------------------------------

def estimate_nonlinear_frequency_by_quadrature(
    depth: float,
    waist: float,
    center: float,
    mass: float,
    amplitude: float,
) -> float:
    """用一维能量积分计算有限振幅非线性频率 (Hz)。

    对势 ``U(x) = -U0 exp(-2 ((x-xc)/w)^2)``，给定振幅 ``A``（相对中心的位移，
    即转折点 ``x_t = xc + A``），能量 ``E = U(x_t)``。周期

        T = 4 ∫_0^A du / sqrt( 2 (E - U(xc+u)) / m ) .

    用无量纲坐标 ``s = u / w``，被积函数在 ``s → A/w`` 处有可积奇点，因此在
    ``[0, A/w)`` 上用 ``quad``（含 ``points`` 处理奇点，并把上限略向内取一点
    单独积分）。

    谐振子极限（A→0）应回到 ``f0``。
    """
    if depth <= 0.0 or waist <= 0.0 or mass <= 0.0:
        raise ValueError("depth、waist、mass 必须为正")
    if amplitude <= 0.0:
        # 零振幅：返回解析小振幅频率
        return trap_frequency_analytic(depth, waist, mass)

    s_max = amplitude / waist  # 无量纲转折点
    e_turn = float(gaussian_potential(center + amplitude, depth, waist, center))

    def integrand(s: float) -> float:
        # u = s*w；E - U(xc + s w) >= 0 在 [0, s_max]
        u = s * waist
        pot = float(gaussian_potential(center + u, depth, waist, center))
        denom = e_turn - pot
        if denom <= 0.0:
            return 0.0
        return waist / math.sqrt(2.0 * denom / mass)

    # 在 s_max 附近奇点，把积分拆成 [0, s_max - eps] + [s_max - eps, s_max*(1-1e-12)]
    eps = 1.0e-9 * max(s_max, 1.0)
    a1, b1 = 0.0, max(s_max - eps, 0.0)
    integral1, _ = integrate.quad(integrand, a1, b1, limit=400)
    # 末段：变量替换 s = s_max - t^2 消除 1/sqrt 奇点
    a2 = max(s_max - eps, 0.0)
    b2 = s_max * (1.0 - 1.0e-12)
    if b2 > a2:
        # t = sqrt(s_max - s)；ds = -2t dt
        def integrand_t(t: float) -> float:
            s = s_max - t * t
            if s < 0.0:
                return 0.0
            u = s * waist
            pot = float(gaussian_potential(center + u, depth, waist, center))
            denom = e_turn - pot
            if denom <= 0.0:
                return 0.0
            g = waist / math.sqrt(2.0 * denom / mass)
            return g * 2.0 * t  # ds = -2t dt，方向吸收负号

        t_lo = math.sqrt(max(s_max - b2, 0.0))
        t_hi = math.sqrt(max(s_max - a2, 0.0))
        integral2, _ = integrate.quad(integrand_t, t_lo, t_hi, limit=400)
    else:
        integral2 = 0.0

    quarter = integral1 + integral2
    if quarter <= 0.0:
        return float("nan")
    period = 4.0 * quarter
    return 1.0 / period


# ---------------------------------------------------------------------------
# 静态势阱汇总
# ---------------------------------------------------------------------------

def summarize_static_trap(
    name: str,
    depth_J: float,
    waist_m: float,
    center_m: float,
    mass_kg: float,
    times: np.ndarray,
    positions: np.ndarray,
    velocities: np.ndarray,
    accelerations: np.ndarray,
    step_waist: float,
    frequency_estimator: str,
    minimum_complete_cycles: int,
    curv_tol: float,
    harm_tol: float,
    traj_tol: float,
    energy_tol: float,
    sources: dict | None = None,
    waist_assumed: bool = False,
) -> dict:
    """汇总单势阱的全部指标（见模块文档字符串的误差分类）。"""
    sources = sources or {}

    def U(x):
        return gaussian_potential(x, depth_J, waist_m, center_m)

    def F(x):
        return gaussian_force(x, depth_J, waist_m, center_m)

    # 1) 平衡位置
    eq = find_equilibrium_position(U, F, center_m, waist_m)

    # 2) 解析与数值曲率（h 与 h/2）
    curv = numerical_curvature_pair(U, depth_J, waist_m, center_m, step_waist)

    # 3) 四个频率
    f_analytic = trap_frequency_analytic(depth_J, waist_m, mass_kg)
    f_num_curv_h2 = numerical_trap_frequency(curv.numerical_h2, mass_kg)
    # 初始振幅 = |x(0) - center|；用于非线性积分基准
    amplitude = abs(float(positions[0]) - center_m)
    f_nonlinear = estimate_nonlinear_frequency_by_quadrature(
        depth_J, waist_m, center_m, mass_kg, amplitude
    )

    # 4) 能量
    ke = 0.5 * mass_kg * velocities * velocities
    pe = gaussian_potential(positions, depth_J, waist_m, center_m)
    total = ke + pe
    e0 = float(total[0])
    dE = total - e0
    dE_over_U0 = dE / depth_J
    abs_dE_over_U0 = np.abs(dE_over_U0)
    max_abs_err = float(np.max(abs_dE_over_U0))
    rms_err = float(np.sqrt(np.mean(dE_over_U0 * dE_over_U0)))
    pkpk_err = float(np.max(dE_over_U0) - np.min(dE_over_U0))
    # 可选：能量误差线性拟合斜率（仅作诊断；不得把有界振荡误称漂移）
    slope, intercept = np.polyfit(times, dE_over_U0, 1)

    # 5) 轨迹频率
    disp = positions - center_m
    tfreq = estimate_frequency_from_trajectory(
        times, disp, method=frequency_estimator, minimum_complete_cycles=minimum_complete_cycles
    )

    # 6) 误差分类
    discret_err = calculate_relative_error(f_num_curv_h2, f_analytic)
    anharmonic_shift = (f_analytic - f_nonlinear) / f_analytic if f_analytic > 0 else float("nan")
    traj_vs_nonlinear = calculate_relative_error(tfreq.frequency_hz, f_nonlinear)
    traj_vs_small = calculate_relative_error(tfreq.frequency_hz, f_analytic)

    # 7) 验收状态
    is_bound = e0 < 0.0
    validation = {
        "curvature_relative_error_within_tolerance": bool(discret_err < curv_tol),
        "harmonic_frequency_within_tolerance": bool(discret_err < harm_tol),
        "trajectory_frequency_within_tolerance": bool(traj_vs_small < traj_tol),
        "max_energy_error_within_tolerance": bool(max_abs_err < energy_tol),
        "is_bound": bool(is_bound),
        "cycle_count_sufficient": bool(tfreq.cycle_count_sufficient),
    }

    metrics = {
        "trap": name,
        "inputs": {
            "depth_uK": float(j_to_uk(depth_J)),
            "waist_um": waist_m * 1.0e6,
            "center_um": center_m * 1.0e6,
            "mass_u": None,  # 由调用者填充
            "depth_J": depth_J,
            "waist_m": waist_m,
            "center_m": center_m,
            "mass_kg": mass_kg,
            "sources": sources,
            "waist_assumed": bool(waist_assumed),
        },
        "equilibrium": {
            "x_eq_m": eq.x_eq,
            "center_deviation_m": eq.x_eq - center_m,
            "potential_at_eq_J": eq.u_eq,
            "potential_at_eq_uK": float(j_to_uk(eq.u_eq)),
            "force_residual_N": eq.force_residual,
        },
        "curvature": {
            "analytic_N_per_m": curv.analytic,
            "numerical_h_N_per_m": curv.numerical_h,
            "numerical_h2_N_per_m": curv.numerical_h2,
            "convergence_diff_N_per_m": curv.convergence_diff,
            "step_h_m": curv.step_h,
            "step_h_waist": step_waist,
        },
        "frequencies": {
            "analytic_small_amp_Hz": f_analytic,
            "numerical_curvature_Hz": f_num_curv_h2,
            "nonlinear_quadrature_Hz": f_nonlinear,
            "trajectory_estimated_Hz": float(tfreq.frequency_hz),
            "trajectory_method": tfreq.method,
        },
        "fft_diagnostics": {
            "fft_resolution_Hz": tfreq.fft_resolution_hz,
            "fft_peak_freq_Hz": tfreq.fft_peak_freq_hz,
            "complete_cycles": tfreq.complete_cycles_estimate,
            "num_zero_crossings": tfreq.num_zero_crossings,
            "cycle_count_sufficient": bool(tfreq.cycle_count_sufficient),
        },
        "error_classification": {
            "discretization_error_relative": discret_err,
            "anharmonic_shift_relative": anharmonic_shift,
            "trajectory_vs_nonlinear_relative": traj_vs_nonlinear,
            "trajectory_vs_small_amp_relative": traj_vs_small,
        },
        "energy": {
            "initial_total_energy_J": e0,
            "initial_total_energy_uK": float(j_to_uk(e0)),
            "is_bound": bool(is_bound),
            "amplitude_m": amplitude,
            "amplitude_waist": amplitude / waist_m,
            "max_abs_energy_error_over_depth": max_abs_err,
            "rms_energy_error_over_depth": rms_err,
            "peak_to_peak_energy_error_over_depth": pkpk_err,
            "linear_fit_slope_per_s": float(slope),
            "linear_fit_intercept": float(intercept),
        },
        "extremes": {
            "max_displacement_from_center_m": float(np.max(np.abs(disp))),
            "max_velocity_m_per_s": float(np.max(np.abs(velocities))),
            "max_acceleration_m_per_s2": float(np.max(np.abs(accelerations))),
        },
        "validation": validation,
        "tolerances": {
            "curvature_relative_tolerance": curv_tol,
            "harmonic_frequency_relative_tolerance": harm_tol,
            "trajectory_frequency_relative_tolerance": traj_tol,
            "max_energy_error_over_depth": energy_tol,
        },
    }
    return metrics
