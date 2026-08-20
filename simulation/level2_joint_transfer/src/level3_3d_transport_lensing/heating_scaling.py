"""简谐、无透镜、小激发极限的加热标度律检验（任务书 §5A/§5B）。

方法：对移动简谐阱 ξ̈ + ω²ξ = −c̈(t)（ξ 为相对阱心坐标，冷初态 ξ=ξ̇=0），
末端激发 ε_f = ½m(ξ̇² + ω²ξ²) 由解析频域求积
ξ(T) = −∫c̈(t)sin[ω(T−t)]dt/ω、ξ̇(T) = −∫c̈(t)cos[ω(T−t)]dt 给出，
并与显式 velocity-Verlet 时间步进互证。

ε_f ∝ |A(ω)|²，A = ∫c̈e^{iωt}dt 的渐近阶由轨迹在端点的光滑阶 m 决定
（首个非零跳变的导数阶为 m+1）：
  m=2（加速度连续、端点为零：sine、constant-jerk、minimum-jerk）
      Δε ∝ L²/(ω⁴T⁶)，即 ΔN = Δε/ħω ∝ L²/(ω⁵T⁶)；
  m=1（加速度分段常量、三角速度）
      Δε ∝ L²/(ω²T⁴)，即 ΔN ∝ L²/(ω³T⁴)。
响应含 |e^{iωT}−1| 类干涉相位因子，在特定 T 处有节点，因此 T 标度按
上包络拟合；L、ω 标度在固定 T 的密集小窗口内取包络后拟合。

注意：论文 Methods（arXiv:2403.12021v4 行 726 附近）声称 constant-jerk
满足 ΔN ∝ D²/(ω³T⁴)，但该指数属于 m=1（加速度不连续）族；本模块用
const_accel 参照轨迹证明方法能同时复现两类指数，并如实报告四段
constant-jerk（按任务书 §5B 定义，加速度连续）实测为 T⁻⁶ω⁻⁴（Δε）。
"""
from __future__ import annotations

import numpy as np

from .gaussian_3d import KB
from .transport_profiles import get_profile

CONST_ACCEL = "const_accel_reference"  # m=1 参照族：三角速度，a = ±4L/T²
M2_FAMILIES = ("adiabatic_sine", "constant_jerk", "minimum_jerk")
T_GRID_US = np.arange(350.0, 1201.0, 25.0)
T_WINDOW_US = np.arange(750.0, 856.0, 25.0)  # L/ω 标度固定 T 的包络窗口
L_GRID_UM = (270.0, 510.0, 610.0)
OMEGA_DEPTHS_UK = (160.0, 280.0, 920.0)  # ω = ω₀√(D/D₀)
FLOOR_UK = 1e-7  # 包络拟合下限（数值地板之上）
ENV_HALF_WIN = 5


def _cddot(profile, s, distance_m, duration_s):
    """归一化加速度 c̈(t)；const_accel 参照族 a = ±4L/T²。"""
    if profile == CONST_ACCEL:
        return np.where(np.asarray(s) < 0.5, 4.0, -4.0) * distance_m / duration_s ** 2
    return get_profile(profile).q_ddot(s) * distance_m / duration_s ** 2


def excitation_quadrature(profile, duration_us, distance_um, omega, mass_kg, n=20001):
    """解析频域求积的末端激发能量（J）。"""
    T = duration_us * 1e-6
    L = distance_um * 1e-6
    t = np.linspace(0.0, T, n)
    a = _cddot(profile, t / T, L, T)
    w = np.full(n, T / (n - 1))
    w[0] = w[-1] = 0.5 * T / (n - 1)
    phase = omega * (T - t)
    xi = -np.sum(a * np.sin(phase) * w) / omega
    xi_dot = -np.sum(a * np.cos(phase) * w)
    return 0.5 * mass_kg * (xi_dot ** 2 + omega ** 2 * xi ** 2)


def excitation_verlet(profile, duration_us, distance_um, omega, mass_kg, dt_us=0.05):
    """velocity-Verlet 时间步进互证（与主积分器同格式的一维简谐推进）。"""
    T = duration_us * 1e-6
    L = distance_um * 1e-6
    dt = dt_us * 1e-6
    n = int(round(T / dt))
    xi = 0.0
    vi = 0.0

    def acc(t):
        s = min(max(t / T, 0.0), 1.0)
        return -omega ** 2 * xi - float(_cddot(profile, np.array(s), L, T))

    a_i = acc(0.0)
    for i in range(n):
        xi += vi * dt + 0.5 * a_i * dt * dt
        a_new = acc((i + 1) * dt)
        vi += 0.5 * (a_i + a_new) * dt
        a_i = a_new
    return 0.5 * mass_kg * (vi ** 2 + omega ** 2 * xi ** 2)


def _envelope(values):
    """滑动窗口最大值包络（压制干涉相位节点）。"""
    win = ENV_HALF_WIN
    return np.array([values[max(0, i - win):i + win + 1].max()
                     for i in range(len(values))])


def _slope(x, y, floor):
    env = _envelope(y)
    sel = env > floor
    if sel.sum() < 3:
        return float("nan"), 0
    return float(np.polyfit(np.log(x[sel]), np.log(env[sel]), 1)[0]), int(sel.sum())


def _t_slope(profile, omega, mass_kg, distance_um=610.0):
    eps = np.array([excitation_quadrature(profile, T, distance_um, omega, mass_kg)
                    for T in T_GRID_US]) / KB * 1e6
    return _slope(T_GRID_US, eps, FLOOR_UK)


def _l_slope(profile, omega, mass_kg, duration_us=800.0):
    vals = []
    for L in L_GRID_UM:
        eps = np.array([excitation_quadrature(profile, T, L, omega, mass_kg)
                        for T in T_WINDOW_US]) / KB * 1e6
        vals.append(eps.max())
    return float(np.polyfit(np.log(L_GRID_UM), np.log(vals), 1)[0])


def _omega_slope(profile, mass_kg, waist_m):
    """用阱深改变 ω = √(4D/mw²)，固定 T 窗口包络后拟合 Δε 对 ω 的斜率。"""
    omegas, vals = [], []
    for depth_uK in OMEGA_DEPTHS_UK:
        depth_j = depth_uK * 1e-6 * KB
        omega = float(np.sqrt(4.0 * depth_j / (mass_kg * waist_m ** 2)))
        eps = np.array([excitation_quadrature(profile, T, 610.0, omega, mass_kg)
                        for T in T_WINDOW_US]) / KB * 1e6
        omegas.append(omega)
        vals.append(eps.max())
    return float(np.polyfit(np.log(omegas), np.log(vals), 1)[0])


def run_heating_scaling_check(physics, light=False) -> dict:
    """执行 §5A/§5B 加热标度律检验；light=True 用粗网格供测试快速回归。

    通过判据（门控）：求积与步进互证；m=2 族 T 斜率 ≈ −6、L ≈ +2；
    参照族 const_accel T ≈ −4、ω ≈ −2（证明方法能复现两类解析指数）。
    constant_jerk 对任务书 §5B/论文声称的 T⁻⁴ 的一致性单独如实报告，
    不作为门控（见模块 docstring 的差异说明）。
    """
    mass = physics["mass_kg"]
    omega0 = physics["omega_r"]

    # light 模式必须保持 25 μs 步长（干涉周期 2π/ω≈27.8 μs，更粗会混叠），
    # 只缩短范围以减少运行时间。
    t_grid = np.arange(350.0, 876.0, 25.0) if light else T_GRID_US
    l_grid = (270.0, 610.0) if light else L_GRID_UM
    t_window = (750.0, 800.0) if light else T_WINDOW_US
    fam_t = {}
    for fam in (*M2_FAMILIES, CONST_ACCEL):
        eps = np.array([excitation_quadrature(fam, T, 610.0, omega0, mass)
                        for T in t_grid]) / KB * 1e6
        slope, n_pts = _slope(t_grid, eps, FLOOR_UK)
        fam_t[fam] = {"t_slope": slope, "n_fit_points": n_pts}

    cross = {}
    ok_cross = True
    for fam in M2_FAMILIES:
        rel = []
        for T in ((600.0,) if light else (600.0, 900.0)):
            q = excitation_quadrature(fam, T, 610.0, omega0, mass)
            v = excitation_verlet(fam, T, 610.0, omega0, mass)
            rel.append(abs(q - v) / max(q, v))
        cross[fam] = {"max_rel_diff": float(max(rel))}
        ok_cross = ok_cross and max(rel) < 0.02

    l_slopes = {}
    for fam in M2_FAMILIES:
        vals = [max(excitation_quadrature(fam, T, L, omega0, mass)
                    for T in t_window) / KB * 1e6 for L in l_grid]
        l_slopes[fam] = float(np.polyfit(np.log(l_grid), np.log(vals), 1)[0])

    omega_slopes = {}
    if not light:
        for fam in ("adiabatic_sine", CONST_ACCEL):
            omega_slopes[fam] = _omega_slope(fam, mass, physics["waist_m"])

    sine_ok = (abs(fam_t["adiabatic_sine"]["t_slope"] + 6.0) <= 1.0
               and abs(l_slopes["adiabatic_sine"] - 2.0) <= 0.5)
    ref_ok = abs(fam_t[CONST_ACCEL]["t_slope"] + 4.0) <= 1.0
    jerk_ok = (abs(fam_t["constant_jerk"]["t_slope"] + 6.0) <= 1.0
               and abs(l_slopes["constant_jerk"] - 2.0) <= 0.5)
    mjerk_ok = abs(fam_t["minimum_jerk"]["t_slope"] + 6.0) <= 1.0
    if not light:  # ω 斜率仅在完整模式下测量（light 供快速回归）
        sine_ok = sine_ok and abs(omega_slopes["adiabatic_sine"] + 4.0) <= 1.0
        ref_ok = ref_ok and abs(omega_slopes[CONST_ACCEL] + 2.0) <= 1.0

    jerk_claim = {
        "prompt_5B_claimed_delta_n": "L^2/(omega^3 T^4)",
        "measured_delta_eps_exponents": {
            "t_slope": fam_t["constant_jerk"]["t_slope"],
            "l_slope": l_slopes["constant_jerk"],
            "omega_slope": omega_slopes.get("adiabatic_sine"),
            "note": "constant_jerk 的 omega 斜率与 sine 同为 m=2 类，未单独重测"},
        "consistent_with_claim": bool(
            abs(fam_t["constant_jerk"]["t_slope"] + 4.0) <= 1.0),
        "explanation": (
            "四段 constant-jerk（任务书 §5B 定义：加速度连续、端点为零）属于"
            "端点光滑阶 m=2 族，解析与实测均为 ΔN ∝ L²/(ω⁵T⁶)；声称的"
            " T⁻⁴/ω³ 指数属于 m=1（分段恒加速度、三角速度）族，已由 "
            f"{CONST_ACCEL} 参照轨迹实测复现（方法自洽）。差异如实报告，"
            "不以门控强制通过或隐藏。"),
    }

    passed = bool(ok_cross and sine_ok and ref_ok and jerk_ok and mjerk_ok)
    return {
        "crosscheck_quadrature_vs_verlet": cross,
        "families": fam_t,
        "l_slopes": l_slopes,
        "omega_slopes_delta_eps": omega_slopes,
        "expected": {"m2_families": "T^-6 L^+2 (Δε ω^-4, ΔN ω^-5)",
                     CONST_ACCEL: "T^-4 (Δε ω^-2, ΔN ω^-3)"},
        "constant_jerk_vs_prompt_5B": jerk_claim,
        "plot_data": {
            "t_grid_us": t_grid.tolist(),
            "envelope_uK": {fam: (_envelope(
                np.array([excitation_quadrature(fam, T, 610.0, omega0, mass)
                          for T in t_grid]) / KB * 1e6)).tolist()
                for fam in (*M2_FAMILIES, CONST_ACCEL)},
        },
        "passed": passed,
    }
